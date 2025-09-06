################################################################################
# Aer’s Matrix Product State (MPS) simulator
# it uses the same algorithm but the MPS method in Aer,
# which stores states as tensor networks rather than dense vectors.
# This typically handles 30–40 qubits on a laptop if entanglement isn’t extreme.
################################################################################
import math
import random
import argparse
from fractions import Fraction

from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import QFTGate
from qiskit_aer import AerSimulator


# Controlled modular multiplication by a mod N (toy placeholder)
def c_mult_mod_N(a, N, n):
    # Controlled multiplication by a mod N on n-qubit register (simplified)
    qc = QuantumCircuit(n)
    for i in range(n-1):
        qc.swap(i, i+1)
    gate = qc.to_gate(label=f"c_mult_{a}_mod{N}")
    return gate.control()


def qpe_amodN(a, N, n_count):
    n = math.ceil(math.log2(N))
    qc = QuantumCircuit(n_count + n, n_count)
    # Counting qubits
    qc.h(range(n_count))
    # Target register initialized to |1>
    qc.x(n_count + n - 1)

    # Controlled-U operations
    for i in range(n_count):
        qc.append(c_mult_mod_N(pow(a, 2**i, N), N, n),
                  [i] + list(range(n_count, n_count + n)))

    # Apply inverse QFT using QFTGate().inverse()
    qc.append(QFTGate(n_count).inverse(), range(n_count))

    qc.measure(range(n_count), range(n_count))
    return qc


def shor(N, max_trials=5):
    if N % 2 == 0:
        return 2
    for trial in range(max_trials):
        a = random.randrange(2, N)
        g = math.gcd(a, N)
        if g > 1:
            return g

        n_count = 2 * math.ceil(math.log2(N))
        qc = qpe_amodN(a, N, n_count)
        backend = AerSimulator(method="matrix_product_state")
        qc = transpile(qc, backend)
        job = backend.run(qc, shots=1)
        result = job.result()
        counts = result.get_counts()
        measured = max(counts, key=counts.get) if counts else "0"*n_count
        decimal = int(measured, 2)
        frac = Fraction(decimal, 2**n_count).limit_denominator(N)
        r = frac.denominator
        if r % 2 != 0:
            continue
        guess = math.gcd(a**(r//2) - 1, N)
        if guess not in [1, N]:
            return guess
    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=16837, help="Number to factor")
    parser.add_argument("--max-trials", type=int, default=5, help="Max random trials")
    args = parser.parse_args()

    factor = shor(args.N, args.max_trials)
    if factor is None:
        print(f"Failed to factor {args.N}")
    else:
        print(f"One non-trivial factor of {args.N} is {factor}")
        print(f"The other factor is {args.N // factor}")
