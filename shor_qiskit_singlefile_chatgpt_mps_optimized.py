
# shor_qiskit_singlefile_optimized.py
# Qiskit 2.1.1 / Python 3.11
# - Iterative Phase Estimation (1 counting qubit)
# - Modular arithmetic via ripple-carry add-constant with O(1) ancillas
# - Controlled modular exponentiation by repeated squaring
# - Aer Matrix Product State (MPS) backend to avoid 2^n memory
#
# Note: For general N, at least one n-qubit work register is required for
# reversible modular multiplication. This design uses exactly one such register
# plus two single-qubit ancillas, and reuses them across IPEA bits.

import math
import random
import argparse
from fractions import Fraction

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator

# -----------------------------
# Classical helpers
# -----------------------------

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

# -----------------------------
# Cuccaro ripple-carry adder: b <- a + b, a and b same length (little-endian)
# Uses a single carry ancilla qubit.
# -----------------------------

def cuccaro_adder(qc: QuantumCircuit, a, b, carry):
    n = len(a)
    # Majority
    qc.cx(a[0], b[0])
    qc.cx(carry, b[0])
    qc.ccx(a[0], carry, b[0])
    for i in range(1, n):
        qc.cx(a[i], b[i])
        qc.cx(b[i-1], b[i])
        qc.ccx(a[i], b[i-1], b[i])
    # Unmajority and add
    for i in reversed(range(1, n)):
        qc.ccx(a[i], b[i-1], b[i])
        qc.cx(b[i-1], b[i])
        qc.cx(a[i], b[i])
    qc.ccx(a[0], carry, b[0])
    qc.cx(carry, b[0])
    qc.cx(a[0], b[0])

def add_constant_inplace(qc: QuantumCircuit, target, C_bits, temp, carry):
    """target <- target + C using Cuccaro, with temp prepared as the constant."""
    for i, bit in enumerate(C_bits):
        if bit:
            qc.x(temp[i])
    cuccaro_adder(qc, temp, target, carry)
    for i, bit in enumerate(C_bits):
        if bit:
            qc.x(temp[i])

def add_constant_mod_N(qc: QuantumCircuit, target, C, N, temp, carry, flag):
    """Compute target = (target + C) mod N.
       target,temp: n-qubit little-endian; carry,flag: single ancillas |0>.
    """
    n = len(target)
    C_mod = C % N
    if C_mod == 0:
        return
    # target += C
    add_constant_inplace(qc, target, int_to_bits(C_mod, n), temp, carry)
    # target -= N  via adding (2^n - N)
    N_comp = (1 << n) - N
    add_constant_inplace(qc, target, int_to_bits(N_comp, n), temp, carry)
    # Detect underflow: if target >= 2^n - N then set flag
    threshold = N_comp
    ctrls = [q for i, q in enumerate(target) if ((threshold >> i) & 1)]
    if ctrls:
        qc.mcx(ctrls, flag)
    # If no underflow (flag=1), add N back
    if ctrls:
        qc.x(flag)
    N_bits = int_to_bits(N, n)
    for i, bit in enumerate(N_bits):
        if bit:
            qc.cx(flag, temp[i])
    cuccaro_adder(qc, temp, target, carry)
    for i, bit in enumerate(N_bits):
        if bit:
            qc.cx(flag, temp[i])
    if ctrls:
        qc.x(flag)

def multiply_by_const_mod_N(qc: QuantumCircuit, y, a, N, temp, carry, flag):
    """In-place y <- (a * y) mod N. Uses temp (n-qubit), carry, flag ancillas."""
    n = len(y)
    # Accumulate into temp based on bits of y
    for i in range(n):
        addend = (a * (1 << i)) % N
        if addend == 0:
            continue
        qc.cx(y[i], flag)
        add_constant_mod_N(qc, temp, addend, N, y, carry, flag)
        qc.cx(y[i], flag)
    # XOR swap temp <-> y to move result into y and clear temp
    for i in range(n):
        qc.cx(temp[i], y[i]); qc.cx(y[i], temp[i]); qc.cx(temp[i], y[i])

def controlled_modexp_2k(qc: QuantumCircuit, ctrl, y, a, N, k, temp, carry, flag):
    """Apply |x> -> |(a^(2^k) * x) mod N> controlled by ctrl."""
    n = len(y)
    exp = pow(a, 1 << k, N)
    sub = QuantumCircuit(2*n + 2, name=f"mul_{exp}_mod{N}")
    ya = [sub.qubits[i] for i in range(n)]
    ta = [sub.qubits[n + i] for i in range(n)]
    c = sub.qubits[2*n]
    f = sub.qubits[2*n + 1]
    multiply_by_const_mod_N(sub, ya, exp, N, ta, c, f)
    gate = sub.to_gate(label=f"*{exp} mod {N}")
    cgate = gate.control(1)
    qc.append(cgate, qargs=[ctrl, *y, *temp, carry, flag])

# -----------------------------
# Iterative Phase Estimation (1 probe qubit)
# -----------------------------

def iterative_phase_estimation(a, N, m=None, verbose=False):
    assert 1 < a < N and is_coprime(a, N)
    n = (N - 1).bit_length()
    m = m if m is not None else 2*n

    y = QuantumRegister(n, 'y')      # target
    t = QuantumRegister(n, 't')      # work/temp
    carry = QuantumRegister(1, 'c')  # carry ancilla
    flag = QuantumRegister(1, 'f')   # flag ancilla
    probe = QuantumRegister(1, 'p')  # counting qubit (iterative)
    cbit = ClassicalRegister(1, 'cb')
    qc = QuantumCircuit(y, t, carry, flag, probe, cbit, name="IPEA")

    # Prepare |1> (simple eigenvector surrogate)
    qc.x(y[0])

    backend = AerSimulator(method="matrix_product_state")

    bits = []
    for k in range(m):
        qc.h(probe)
        power_index = m - 1 - k
        controlled_modexp_2k(qc, probe[0], y, a, N, power_index, t, carry[0], flag[0])

        # Feed-forward phase correction based on prior bits
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
        counts = res.get_counts()
        bit = int(max(counts, key=counts.get)) if counts else 0
        bits.append(bit)
        if verbose:
            print(f"[ipea] bit {k+1}/{m}: {bit}")
        qc.reset(probe)

    return bits

# -----------------------------
# Shor driver
# -----------------------------

def shor_once(N, m=None, seed=None, verbose=True):
    if N % 2 == 0:
        return (2, N // 2)
    for small in [3,5,7,11,13,17,19,23,29,31,37,41,43,47]:
        if N % small == 0 and small not in (1, N):
            return (small, N // small)

    rng = random.Random(seed)
    while True:
        a = rng.randrange(2, N-1)
        if is_coprime(a, N):
            break
    if verbose:
        print(f"[info] N={N}, a={a}")

    bits = iterative_phase_estimation(a, N, m=m, verbose=verbose)
    frac = continued_fraction_phase(bits)
    r = frac.denominator
    if verbose:
        print(f"[info] bits={bits} -> approx {frac.numerator}/{frac.denominator}  (r≈{r})")
    fac = try_factors_from_order(a, r, N)
    return fac

def shor(N, max_trials=8, m=None, seed=None, verbose=True):
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
            print("[retry] unsuitable order; retrying...\\n")
    if verbose:
        print("[fail] no nontrivial factors found.")
    return None

def main():
    parser = argparse.ArgumentParser(description="Single-file Shor (Qiskit 2.1.1, optimized, MPS)")
    parser.add_argument("--N", type=int, default=16837, help="Semiprime to factor")
    parser.add_argument("--max-trials", type=int, default=8)
    parser.add_argument("--m", type=int, default=None, help="phase bits (default 2*ceil(log2 N))")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--quiet", action="store_true", help="suppress verbose logs")
    args = parser.parse_args()

    fac = shor(args.N, max_trials=args.max_trials, m=args.m, seed=args.seed, verbose=(not args.quiet))
    if fac:
        print(f"Result: N={args.N} -> {fac[0]} * {fac[1]}")
    else:
        print("Failed to factor N.")

if __name__ == "__main__":
    main()
