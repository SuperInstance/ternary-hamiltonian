# PLUG_AND_PLAY — Hamiltonian

> Hamiltonian mechanics on ternary phase space

## 🚀 Quick Start

Add to your `Cargo.toml`:

```toml
[dependencies]
ternary-hamiltonian = { git = "https://github.com/SuperInstance/ternary-hamiltonian" }
```

Use in your code:

```rust
use ternary_hamiltonian::{
    EnergyConservation, Hamiltonian, PhaseSpace, SymplecticIntegrator, TernaryCoupling,
};

fn main() {
    // One degree of freedom: q = +1, p = 0.
    let mut phase = PhaseSpace::new(vec![1], vec![0]);
    let h = Hamiltonian::new(1.0, 1.0);          // kinetic T, potential V
    let coupling = TernaryCoupling::harmonic();  // α = β = 1
    let mut tracker = EnergyConservation::new(h.energy_at(&phase));

    for step in 0..10 {
        phase = SymplecticIntegrator::symplectic_euler(&phase, &coupling);
        tracker.record(h.energy_at(&phase));
        println!("step {step}: q = {}, p = {}", phase.positions[0], phase.momenta[0]);
    }

    println!("max energy drift = {}", tracker.drift());
}
```

> The snippet above is kept in sync with the code by `tests/readme_example.rs`
> (the equivalent of this example is compiled and run by `cargo test`).

## API surface

The crate exposes these public types (no `ternary-types` dependency — phase-space
coordinates are plain `i8` values in `{-1, 0, +1}`):

| Type | Role |
|------|------|
| `PhaseSpace` | A point `(q, p)` with each component in `{-1, 0, +1}`. |
| `Hamiltonian` | Energy `H = T + V`; `energy_at(&phase)` evaluates `½ Σᵢ(T·pᵢ² + V·qᵢ²)`. |
| `TernaryCoupling` | Coupling constants `α`, `β ∈ Z₃` for the discrete flow. |
| `SymplecticIntegrator` | `symplectic_euler` / `stormer_verlet` steps (Z₃ modular arithmetic). |
| `EnergyConservation` | Records energy samples and reports drift `maxₙ |Eₙ − E₀|`. |
| `PoissonBracket` | Discrete Poisson bracket of observables. |
| `LiouvilleTheorem` | Distinct-cell volume counting / conservation check. |

## 🔗 Integration

This crate is part of the [SuperInstance ternary fleet](https://github.com/SuperInstance).
It operates on plain `i8` ternary values directly; it does **not** depend on a
separate `ternary-types` crate.

## 📄 License

MIT
