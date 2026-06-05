# ternary-hamiltonian

**Hamiltonian mechanics on ternary phase space: symplectic integration that preserves energy forever.**

Hamilton's equations describe conservative dynamical systems where energy is exactly conserved. This crate implements symplectic integrators (symplectic Euler, Störmer-Verlet) on ternary phase space grids where positions and momenta live in {-1, 0, +1}.

The key property: **symplectic integrators preserve the symplectic 2-form ω = dp ∧ dq exactly**. This means energy oscillates but never drifts — unlike Runge-Kutta methods where energy slowly leaks away.

---

## The Mathematics

Hamilton's equations:
```
dq/dt =  ∂H/∂p    (position update)
dp/dt = -∂H/∂q    (momentum update)
```

Where H(q,p) = T(p) + V(q) is the Hamiltonian (kinetic + potential energy).

**Symplectic Euler** (1st order):
```
p(n+1) = p(n) - Δt · ∂V/∂q(q(n))
q(n+1) = q(n) + Δt · ∂T/∂p(p(n+1))
```

**Störmer-Verlet** (2nd order):
```
p(n+½) = p(n) - ½Δt · ∂V/∂q(q(n))
q(n+1) = q(n) + Δt · ∂T/∂p(p(n+½))
p(n+1) = p(n+½) - ½Δt · ∂V/∂q(q(n+1))
```

The Störmer-Verlet method is time-reversible and symplectic — properties that ensure long-term stability.

---

## Architecture

- **`Hamiltonian`** — H(q,p) = kinetic(p) + potential(q)
- **`PhaseSpace`** — Positions q and momenta p on ternary grid
- **`SymplecticIntegrator`** — euler_step() and verlet_step()
- **`EnergyConservation`** — Track total energy drift over N steps
- **`PoissonBracket`** — Compute {f, g} = Σ (∂f/∂qᵢ · ∂g/∂pᵢ - ∂f/∂pᵢ · ∂g/∂qᵢ)
- **`LiouvilleTheorem`** — Verify phase space volume preservation

---

## Quick Start

```rust
use ternary_hamiltonian::{Hamiltonian, PhaseSpace, SymplecticIntegrator, EnergyConservation};

let h = Hamiltonian::harmonic_oscillator(1.0); // k=1 spring
let mut state = PhaseSpace::new(vec![1.0], vec![0.0]); // q=1, p=0

let mut tracker = EnergyConservation::new();
for _ in 0..10000 {
    SymplecticIntegrator::verlet_step(&h, &mut state, 0.01);
    tracker.record(h.total_energy(&state));
}

println!("Energy drift: {:.6}", tracker.drift());
// Expected: < 0.001 (symplectic preserves energy!)
```

---

## Connection to Noether's Theorem

Hamiltonian systems with symmetries have conserved quantities (Noether's theorem). Our companion crate **ternary-noether** numerically verifies this: it defines discrete symmetries on ternary grids and checks that the corresponding quantities (energy, momentum, angular momentum) are conserved by the symplectic integrator.

---

## Ecosystem

- **ternary-noether** — Noether's theorem: symmetry → conservation
- **ternary-electromagnetism** — Yee lattice EM (also symplectic!)
- **ternary-symplectic** — Symplectic geometry primitives
- **ternary-thermostat** — Practical thermostat using PID control

## License

MIT
