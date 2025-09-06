
# shor_qiskit_singlefile.py
import math
import random
import argparse
from fractions import Fraction

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator

def gcd(a, b):
    while b:
        a, b = b, a % b
    return abs(a)

def is_coprime(a, N):
    return gcd(a, N) == 1

def int_to_bits(x, n):
    return [(x >> i) & 1 for i in range(n)]

def continued_fraction_phase(bits):
    """Given measured phase bits (msb to lsb), return best Fraction s/r."""
    m = len(bits)
    val = 0
    for b in bits:
        val = (val << 1) | b
    phase = val / (2**m)
    return Fraction(phase).limit_denominator(2**m)

def try_factors_from_order(a, r, N):
    if r is None or r == 0 or r % 2 == 1:
        return None
    x = pow(a, r // 2, N)
    if x == 1 or x == N - 1:
        return None
    p = gcd(x - 1, N)
    q = gcd(x + 1, N)
    if p * q == N and p not in (1, N) and q not in (1, N):
        return tuple(sorted((p, q)))
    return None

def cuccaro_adder(qc: QuantumCircuit, a, b, cin, cout):
    n = len(a)
    qc.cx(a[0], b[0])
    qc.cx(cin, b[0])
    qc.ccx(a[0], cin, b[0])
    for i in range(1, n):
        qc.cx(a[i], b[i])
        qc.cx(b[i-1], b[i])
        qc.ccx(a[i], b[i-1], b[i])
    qc.cx(b[n-1], cout)
    for i in reversed(range(1, n)):
        qc.ccx(a[i], b[i-1], b[i])
        qc.cx(b[i-1], b[i])
        qc.cx(a[i], b[i])
    qc.ccx(a[0], cin, b[0])
    qc.cx(cin, b[0])
    qc.cx(a[0], b[0])

def add_constant_inplace(qc: QuantumCircuit, target, C_bits, anc_a, cin, cout):
    for i, bit in enumerate(C_bits):
        if bit:
            qc.x(anc_a[i])
    cuccaro_adder(qc, anc_a, target, cin, cout)
    for i, bit in enumerate(C_bits):
        if bit:
            qc.x(anc_a[i])

def add_constant_mod_N(qc: QuantumCircuit, target, C, N, anc_a, cin, cout, anc_flag):
    n = len(target)
    C_bits = int_to_bits(C % N, n)
    add_constant_inplace(qc, target, C_bits, anc_a, cin, cout)

    N_comp_bits = int_to_bits((2**n - N) % (2**n), n)
    add_constant_inplace(qc, target, N_comp_bits, anc_a, cin, cout)

    threshold = (2**n - N) % (2**n)
    ctrl_list = [q for i, q in enumerate(target) if ((threshold >> i) & 1)]
    if ctrl_list:
        qc.mcx(ctrl_list, anc_flag)

    if ctrl_list:
        qc.x(anc_flag)
    N_bits = int_to_bits(N % (2**n), n)
    for i, bit in enumerate(N_bits):
        if bit:
            qc.cx(anc_flag, anc_a[i])
    cuccaro_adder(qc, anc_a, target, cin, cout)
    for i, bit in enumerate(N_bits):
        if bit:
            qc.cx(anc_flag, anc_a[i])
    if ctrl_list:
        qc.x(anc_flag)

def multiply_by_const_mod_N(qc: QuantumCircuit, target, a, N, acc, cin, cout, anc_flag):
    n = len(target)
    ctrls = list(target)
    for i in range(n):
        addend = (a * (1 << i)) % N
        if addend == 0:
            continue
        qc.cx(ctrls[i], anc_flag)
        add_constant_mod_N(qc, acc, addend, N, target, cin, cout, anc_flag)
        qc.cx(ctrls[i], anc_flag)
    for i in range(n):
        qc.cx(acc[i], target[i])
        qc.cx(target[i], acc[i])
        qc.cx(acc[i], target[i])

def controlled_modexp_2k(qc: QuantumCircuit, ctrl_qubit, work, a, N, k, acc, cin, cout, anc_flag):
    n = len(work)
    exp = pow(a, 1 << k, N)
    sub = QuantumCircuit(n + n + 3, name=f"mul_{exp}_mod{N}")
    t = list(range(n))
    s = list(range(n, 2*n))
    c_in, c_out, f = 2*n, 2*n+1, 2*n+2
    multiply_by_const_mod_N(sub, [sub.qubits[i] for i in t], exp, N,
                            [sub.qubits[i] for i in s],
                            sub.qubits[c_in], sub.qubits[c_out], sub.qubits[f])
    gate = sub.to_gate(label=f"*{exp} mod {N}")
    cgate = gate.control(1)
    qc.append(cgate, qargs=[ctrl_qubit, *work, *acc, cin, cout, anc_flag])

def iterative_phase_estimation(a, N, m=None, seed=None, verbose=False):
    assert 1 < a < N and is_coprime(a, N)
    n = math.ceil(math.log2(N))
    m = m if m is not None else 2*n

    work = QuantumRegister(n, "w")
    acc = QuantumRegister(n, "s")
    cin = QuantumRegister(1, "cin")
    cout = QuantumRegister(1, "cout")
    flag = QuantumRegister(1, "f")
    probe = QuantumRegister(1, "p")
    cbit = ClassicalRegister(1, "c")
    qc = QuantumCircuit(work, acc, cin, cout, flag, probe, cbit, name="IPEA")

    qc.x(work[0])
    backend = AerSimulator(method="statevector")

    bits = []
    for k in range(m):
        qc.h(probe)
        power_index = m - 1 - k
        controlled_modexp_2k(qc, probe[0], work, a, N, power_index, acc, cin[0], cout[0], flag[0])

        if bits:
            angle = 0.0
            for j, b in enumerate(bits, start=1):
                if b == 1:
                    angle += math.pi / (2**j)
            qc.p(-2*angle, probe[0])

        qc.h(probe)
        qc.measure(probe, cbit)

        expanded = qc.decompose(reps=10)
        res = backend.run(expanded, shots=1).result()
        bit = int(res.get_counts().most_frequent())
        bits.append(bit)
        if verbose:
            print(f"[ipea] bit {k+1}/{m}: {bit}")
        qc.reset(probe)

    return bits

def shor_once(N, m=None, seed=None, verbose=True):
    if N % 2 == 0:
        return (2, N // 2)
    for small in [3,5,7,11,13,17,19,23,29,31,37,41,43]:
        if N % small == 0 and small not in (1, N):
            return (small, N // small)
    rng = random.Random(seed)
    while True:
        a = rng.randrange(2, N-1)
        if is_coprime(a, N):
            break
    if verbose:
        print(f"[info] N={N}, a={a}")
    bits = iterative_phase_estimation(a, N, m=m, seed=seed, verbose=verbose)
    frac = continued_fraction_phase(bits)
    r = frac.denominator
    if verbose:
        print(f"[info] phase bits={bits} -> approx {frac.numerator}/{frac.denominator}, r≈{r}")
    fac = try_factors_from_order(a, r, N)
    return fac

def shor(N, max_trials=10, m=None, seed=None, verbose=True):
    for t in range(1, max_trials+1):
        if verbose:
            print(f"--- Trial {t}/{max_trials} ---")
        fac = shor_once(N, m=m, seed=(None if seed is None else seed+t), verbose=verbose)
        if fac:
            p, q = fac
            if verbose:
                print(f"[success] {p} × {q} = {N}")
            return fac
        if verbose:
            print("[retry] unsuitable order; retrying...\n")
    if verbose:
        print("[fail] no nontrivial factors found.")
    return None

def main():
    parser = argparse.ArgumentParser(description="Single-file Shor (Qiskit 2.1.1)")
    parser.add_argument("--N", type=int, default=16837)
    parser.add_argument("--max-trials", type=int, default=10)
    parser.add_argument("--m", type=int, default=None, help="phase bits (default 2*ceil(log2 N))")
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()

    fac = shor(args.N, max_trials=args.max_trials, m=args.m, seed=args.seed, verbose=True)
    if fac:
        print(f"Result: N={args.N} -> {fac[0]} * {fac[1]}")
    else:
        print("Failed to factor N.")

if __name__ == "__main__":
    main()
