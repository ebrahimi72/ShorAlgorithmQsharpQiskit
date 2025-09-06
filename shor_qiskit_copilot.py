from math import gcd
from fractions import Fraction
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import Aer

def inverse_qft(n):
    qc = QuantumCircuit(n)
    for j in range(n // 2):
        qc.swap(j, n - j - 1)
    for j in range(n):
        for k in range(j):
            qc.cp(-np.pi / float(2 ** (j - k)), k, j)
        qc.h(j)
    return qc.to_gate(label="QFT†")

def modular_exp_gate(a, N):
    # Placeholder modular exponentiation gate
    qc = QuantumCircuit(4 + 4)  # 4 work qubits + 4 ancillas
    for i in range(4):
        qc.cx(i, 4 + i)
    return qc.to_gate(label=f"ModExp({a} mod {N})")

def estimate_phase(a, N, n_counting=12):
    n_work = 4
    qc = QuantumCircuit(n_counting + n_work * 2)

    qc.h(range(n_counting))

    for i in range(n_counting):
        power = 2 ** i
        gate = modular_exp_gate(a, N)
        qc.append(gate.control(1), [i] + list(range(n_counting, n_counting + n_work * 2)))

    qc.append(inverse_qft(n_counting), range(n_counting))
    qc.save_statevector()

    backend = Aer.get_backend('aer_simulator')
    compiled = transpile(qc, backend)
    job = backend.run(compiled)
    result = job.result()
    statevector = result.get_statevector()

    probs = np.abs(statevector) ** 2
    top_indices = np.argsort(probs)[-5:]
    avg_phase = np.mean([i / 2**n_counting for i in top_indices])
    return avg_phase

def estimate_order(a, N, n_counting=12):
    phase = estimate_phase(a, N, n_counting)
    frac = Fraction(phase).limit_denominator(N)
    return frac.denominator

def factor_semiprime(N, max_attempts=20):
    if N % 2 == 0:
        return [2, N // 2]

    found = set()
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
        if gcd(a, N) != 1:
            found.add(gcd(a, N))
            continue

        try:
            r = estimate_order(a, N)
        except Exception:
            continue

        if r % 2 != 0:
            continue

        x = pow(a, r // 2, N)
        if x == 1 or x == N - 1:
            continue

        for f in [gcd(x - 1, N), gcd(x + 1, N)]:
            if f > 1 and f < N:
                found.add(f)

        if len(found) == 2:
            break

    return sorted(list(found)) if found else None

def main():
    N = 16837  # Try 15, 21, 33, 253 or 16837
    factors = factor_semiprime(N)
    print(f"Found factors of {N}: {factors if factors else 'None'}")

if __name__ == "__main__":
    main()