# CROSS-POLLINATION.md — ternary-hamiltonian

> **Conservation Law Connection:** γ + η = C AS Hamiltonian energy conservation

## Role in the Conservation Law

`ternary-hamiltonian` is the conservation law wearing physics clothing. The Hamiltonian
H = T + V (kinetic + potential energy) maps directly:

- **Kinetic energy T** → γ (productive computation in motion)
- **Potential energy V** → η (stored coordination overhead, ready to become waste)
- **Total energy H** → C (the conserved constant)
- **Phase space** → ternary {−1, 0, +1}ⁿ (each dimension is a Trit)

The Hamiltonian formulation proves that γ + η = C is not an empirical observation
but a **Noether symmetry consequence**: the fleet is invariant under time translation,
therefore energy (C) is conserved.

## delta-clt Verification Results

The delta-clt suite empirically confirms what the Hamiltonian proves theoretically:
- δ(n) is the **fluctuation-dissipation** rate — energy radiated into η as the system equilibrates
- The c₂ = 15/8 correction term in delta_expanded() corresponds to a higher-order
  Hamiltonian perturbation (anharmonic correction)
- The colony drift (26% at n=15) indicates the system is NOT in thermal equilibrium —
  it has not reached the temperature where Noether's theorem applies cleanly

**Implication:** The Hamiltonian formulation predicts that fleet "temperature" (average
agent energy) determines when the conservation law becomes exact. Below a critical
temperature, η fluctuations dominate.

## Cross-Repo Connections

### → ternary-types
`ternary-types` defines the discrete phase space axes. Each Trit is one dimension
of the Hamiltonian phase space. The Hamiltonian acts on `Vec<Trit>`.

**Shared:** Trit is the unit of phase space in both.
**Different:** `ternary-types` is static; `ternary-hamiltonian` adds dynamics and conservation.

### → ternary-pid
`ternary-pid` is the control theory application of the Hamiltonian. PID controllers
drive the fleet toward Hamiltonian equilibrium. The P term responds to current γ/η
imbalance, I accumulates past drift, D predicts future trajectory.

**Shared:** Both are dynamical systems over ternary states.
**Different:** Hamiltonian is descriptive (what IS); PID is prescriptive (what to DO).

### → ternary-rhythm
`ternary-rhythm` studies oscillatory patterns in ternary signals. In Hamiltonian
mechanics, these are the **normal modes** of the system. Each rhythm frequency
corresponds to an eigenvalue of the Hamiltonian.

**Shared:** Both study temporal behavior of ternary systems.
**Different:** Hamiltonian gives the equations of motion; rhythm extracts observed frequencies.

## Fleet Position

```
┌───────────────────────────────────────────────────────┐
│  ternary-hamiltonian — THE PHYSICS LAYER               │
│                                                        │
│  H = T(γ) + V(η) = C     (Noether → time symmetry)    │
│                                                        │
│  Phase space: {−1, 0, +1}ⁿ  (from ternary-types)      │
│  Dynamics:    ċ = ∂H/∂p, ṗ = −∂H/∂c                   │
│  Control:     ternary-pid drives toward H-equilibrium  │
│  Modes:       ternary-rhythm extracts eigenfrequencies │
│                                                        │
│  The conservation law is not imposed — it EMERGES      │
│  from time-translation symmetry.                       │
└───────────────────────────────────────────────────────┘
```

