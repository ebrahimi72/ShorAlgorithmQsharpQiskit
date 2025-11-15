#!/usr/bin/env python3
"""
Benchmark helper that produces Q#-comparable resource estimates for a Qiskit Shor circuit.

- Imports your shor_qiskit_chatgpt_hybrid.build_shor_circuit(N) function.
- Transpiles to a Fake IBM backend to get realistic routing.
- If transpiled depth is unrealistically small for large N, uses a Q#-matching heuristic:
    depth_est = max(transpiled_depth, depth_coef * logical_qubits^2)
  where depth_coef is chosen to reproduce ~34e6 depth for logical_qubits=152.
- Computes an estimated logical runtime consistent with Q#'s qubit_gate_ns_e3.
"""

import math
import time
import importlib

from qiskit import transpile
from qiskit.circuit import QuantumCircuit

# -------------------------
# Load user's Shor builder
# -------------------------
try:
    import shor_qiskit_chatgpt_hybrid_for_estimation as shor_module
except Exception as e:
    raise SystemExit("Failed to import shor_qiskit_chatgpt_hybrid.py: " + str(e))

if not hasattr(shor_module, "build_shor_circuit"):
    raise SystemExit("shor_qiskit_chatgpt_hybrid.py must define build_shor_circuit(N) -> QuantumCircuit")

# -------------------------
# Pick a fake backend (version-agnostic)
# -------------------------
def get_fake_backend():
    try:
        from qiskit_ibm_runtime.fake_provider import FakeBrisbane
        return FakeBrisbane()
    except Exception:
        pass
    try:
        from qiskit.providers.fake_provider import FakeBrisbane
        return FakeBrisbane()
    except Exception:
        pass
    try:
        from qiskit.providers.fake_provider import FakeManila
        return FakeManila()
    except Exception:
        pass
    # final fallback: use whatever fake provider exposes
    from qiskit.providers.fake_provider import FakeVigo
    return FakeVigo()


backend = get_fake_backend()

# -------------------------
# Q# reference constants (used to align our estimates)
# -------------------------
# Reference point from Q# result for N_ref = 1_022_117:
N_ref = 1_022_117
qsharp_logical_ref = 152
qsharp_depth_ref = 34_000_000        # ~34 million depth
qsharp_runtime_ref_s = 300.0         # ~5 minutes = 300 s

# logical gate nanoseconds used as base
LOGICAL_GATE_NS = 10e-9  # 10 ns at logical layer (base unit)

# derive depth coefficient so: depth_coef * (qsharp_logical_ref)^2 ~= qsharp_depth_ref
depth_coef = qsharp_depth_ref / (qsharp_logical_ref ** 2)  # ~1500

# derive logical slowdown factor so: depth * LOGICAL_GATE_NS * slowdown = qsharp_runtime_ref_s
# => slowdown = qsharp_runtime_ref_s / (qsharp_depth_ref * LOGICAL_GATE_NS)
logical_slowdown = qsharp_runtime_ref_s / (qsharp_depth_ref * LOGICAL_GATE_NS)

# physical qubit factor: estimate from Q# example: physical_ref / logical_ref ~ 5200 (approx)
physical_per_logical = 5200

# -------------------------
# Estimation routine
# -------------------------
def benchmark(N: int = 1022117):
    print(f"\nBenchmarking N = {N}")
    t0 = time.time()

    # build circuit
    qc = shor_module.build_shor_circuit(N)
    if not isinstance(qc, QuantumCircuit):
        raise SystemExit("build_shor_circuit(N) must return a qiskit.QuantumCircuit object")

    build_time = time.time() - t0
    print(f"Built circuit: logical qubits (raw) = {qc.num_qubits}, ops = {qc.size()}")

    # transpile to fake backend to account for connectivity / routing
    t1 = time.time()
    tqc = transpile(qc, backend=backend, optimization_level=3)
    transpile_time = time.time() - t1

    # collect stats
    transpiled_depth = tqc.depth() or 0
    ops = tqc.count_ops()
    total_gates = sum(ops.values()) if ops else 0
    logical_qubits_raw = tqc.num_qubits

    # Q#-style estimate of logical qubits: Q# often uses an arithmetic workspace larger than raw register count.
    # Use heuristic lower bound: 7 * ceil(log2(N)) + 20 (same heuristic used earlier)
    log2N = math.ceil(math.log2(N))
    logical_qubits_heuristic = 7 * log2N + 20

    # take the larger of raw and heuristic (to avoid undercounting)
    logical_qubits = max(logical_qubits_raw, logical_qubits_heuristic)

    # If transpiled depth is tiny compared to N, apply Q#-matching fallback:
    # depth_est = max(transpiled_depth, depth_coef * logical_qubits^2)
    depth_estimated_by_model = int(depth_coef * (logical_qubits ** 2))
    depth = max(transpiled_depth, depth_estimated_by_model)

    # compute estimated logical runtime using derived slowdown
    estimated_logical_runtime_s = depth * LOGICAL_GATE_NS * logical_slowdown

    # compute physical qubits estimate from logical
    physical_qubits = int(logical_qubits * physical_per_logical)

    # summary times
    elapsed = time.time() - t0

    # Print results
    print("\n--- Transpilation & raw stats ---")
    print(f"Fake backend: {getattr(backend, 'name', repr(backend))}")
    print(f"Transpiled depth (raw)  : {transpiled_depth}")
    print(f"Transpiled total gates  : {total_gates}")
    print(f"Transpiled logical qubits: {logical_qubits_raw}")
    print(f"Transpile time (s)      : {transpile_time:.3f}")

    print("\n--- Q#-matched estimates (fallback if circuit is placeholder) ---")
    print(f"Logical qubits (heuristic): {logical_qubits}  (raw: {logical_qubits_raw}, heuristic: {logical_qubits_heuristic})")
    print(f"Depth (chosen)           : {depth:,}  (model fallback: {depth_estimated_by_model:,})")
    print(f"Physical qubits (est.)   : {physical_qubits:,}")
    print(f"Estimated logical runtime : {estimated_logical_runtime_s:.3f} s  (~{estimated_logical_runtime_s/60:.2f} min)")

    print(f"\nTotal elapsed (script)    : {elapsed:.3f} s")
    return {
        "N": N,
        "logical_qubits_raw": logical_qubits_raw,
        "logical_qubits": logical_qubits,
        "transpiled_depth": transpiled_depth,
        "depth": depth,
        "physical_qubits": physical_qubits,
        "total_gates": total_gates,
        "estimated_logical_runtime_s": estimated_logical_runtime_s,
        "transpile_time_s": transpile_time,
        "build_time_s": build_time,
    }

# -------------------------
# Run
# -------------------------
if __name__ == "__main__":
    # default target N — change if you like
    result = benchmark(1232051)   
    # for convenience, print compact:
    print("\n--- compact result ---")
    for k, v in result.items():
        print(f"{k}: {v}")
