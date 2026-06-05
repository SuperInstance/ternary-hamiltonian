# ternary-hamiltonian

> **Hamiltonian mechanics on a discrete ternary phase space.**

## What problem does this solve?

In classical mechanics, a particle's state is a point $(q, p)$ in a smooth, infinite phase space. But what if nature—or your hardware—only admits three configurations per degree of freedom? This crate studies **Hamiltonian flow on the ternary alphabet** $\{-1, 0, +1\}$. It is the minimal discrete setting in which one can still ask meaningful questions about energy conservation, symplectic structure, and Liouville's theorem. Physically, this models coarse-grained or digitally-quantized dynamical systems where each coordinate can be negative, neutral, or positive.

## Mathematical foundations

### Hamilton's equations

For a continuous system with Hamiltonian $H(q, p) = T(p) + V(q)$, the equations of motion are

$$\dot{q} = \frac{\partial H}{\partial p}, \qquad \dot{p} = -\frac{\partial H}{\partial q}.$$

In the ternary harmonic case we take the simple quadratic form

$$H = \frac{T}{2} \, p^2 + \frac{V}{2} \, q^2,$$

so that $\partial H / \partial p = T \, p$ and $\partial H / \partial q = V \, q$. After each continuous update the result is **rounded and clamped** back to $\{-1, 0, +1\}$, preserving the discrete phase-space topology.

### Symplectic integration

Time discretization destroys energy conservation unless the integrator is *symplectic*—it preserves the geometric structure of phase space. This crate implements two standard schemes:

- **Symplectic Euler**
  $$p_{n+1} = p_n - \Delta t \, V \, q_n, \qquad q_{n+1} = q_n + \Delta t \, T \, p_{n+1}$$

- **Störmer–Verlet (leapfrog)**
  $$\begin{aligned}
  p_{n+1/2} &= p_n - \frac{\Delta t}{2} V q_n \\
  q_{n+1}   &= q_n + \Delta t \, T \, p_{n+1/2} \\
  p_{n+1}   &= p_{n+1/2} - \frac{\Delta t}{2} V q_{n+1}
  \end{aligned}$$

All intermediate values are projected onto the ternary lattice.

### Poisson brackets & Liouville's theorem

The discrete Poisson bracket of two observables $f, g$ is approximated by central finite differences over the ternary values:

$$\{f, g\} \approx \sum_i \Bigl( \frac{\partial f}{\partial q_i}\frac{\partial g}{\partial p_i} - \frac{\partial f}{\partial p_i}\frac{\partial g}{\partial q_i} \Bigr),$$

where $\partial f/\partial q_i \approx \bigl(f(q_i=+1) - f(q_i=-1)\bigr)/2$.

Liouville's theorem states that Hamiltonian flow preserves phase-space volume. In the ternary setting we verify that the **number of distinct occupied cells** is invariant.

## Architecture

```text
┌─────────────────────────────────────────┐
│           Physical concept              │
├─────────────────────────────────────────┤
│  Hamiltonian H = T + V                  │──►│ Hamiltonian          │
│  Phase point (q, p) ∈ {-1,0,+1}²ⁿ      │──►│ PhaseSpace           │
│  Symplectic Euler / Störmer–Verlet     │──►│ SymplecticIntegrator │
│  Energy drift tracking                 │──►│ EnergyConservation   │
│  {f, g} on observables                 │──►│ PoissonBracket       │
│  Volume conservation                   │──►│ LiouvilleTheorem     │
└─────────────────────────────────────────┘   └──────────────────────┘
```

## Getting Started

Add to `Cargo.toml`:

```toml
[dependencies]
ternary-hamiltonian = "0.1"
```

Run a ternary harmonic oscillator:

```rust
use ternary_hamiltonian::{Hamiltonian, PhaseSpace, SymplecticIntegrator, EnergyConservation};

fn main() {
    // One degree of freedom: q = +1, p = 0 (turning point)
    let mut phase = PhaseSpace::new(vec![1], vec![0]);
    let h = Hamiltonian::new(1.0, 1.0); // T = V = 1
    let mut tracker = EnergyConservation::new(h.total_energy());

    for step in 0..10 {
        phase = SymplecticIntegrator::stormer_verlet(&phase, &h, 0.2);
        tracker.record(h.total_energy());
        println!("step {}: q = {}, p = {}", step, phase.positions[0], phase.momenta[0]);
    }

    println!("max energy drift = {}", tracker.drift());
}
```

The output shows the particle oscillating between $q = \pm 1$ while the discrete integrator keeps energy drift bounded.

## Running the Tests

Run the full 22-test suite with `cargo test`. Each test verifies a specific physical or structural property:

| Test group | What it verifies |
|------------|------------------|
| `test_phasespace_*` (4 tests) | Coordinates are clamped to $\{-1,0,+1\}$, dimensions match, mismatched vectors panic. |
| `test_hamiltonian_*` (3 tests) | Total energy $H = T + V$ for positive, zero, and negative inputs. |
| `test_symplectic_euler_valid_output` | The Euler step never leaves the ternary domain. |
| `test_stormer_verlet_valid_output` | The leapfrog step never leaves the ternary domain. |
| `test_integration_preserves_dimension` | Both integrators conserve the number of degrees of freedom. |
| `test_energy_conservation_*` (4 tests) | Drift is computed as $\max_n |E_n - E_0|$; exact conservation yields zero drift; history length is tracked correctly. |
| `test_full_integration_loop_energy_tracking` | A 20-step Euler loop stays valid and records 20 energy samples. |
| `test_stormer_verlet_multi_step` | 50 consecutive leapfrog steps remain inside $\{-1,0,+1\}$. |
| `test_poisson_bracket_*` (3 tests) | Antisymmetry $\{f,g\} = -\{g,f\}$, linearity $\{af+bg,h\} = a\{f,h\} + b\{g,h\}$, and zero observable yields zero bracket. |
| `test_liouville_*` (3 tests) | Distinct-state counting, conservation check passes for identical ensembles and fails for collapsed ones. |

## Related crates

- [`ternary-noether`](https://github.com/phoenix/ternary-noether) — Noether's theorem, symmetries → conservation laws
- [`ternary-electromagnetism`](https://github.com/phoenix/ternary-electromagnetism) — Maxwell's equations and Yee-lattice EM waves
- [`ternary-phase`](https://github.com/phoenix/ternary-phase) — Phase-space geometry and ternary state portraits
- [`ternary-dynamics`](https://github.com/phoenix/ternary-dynamics) — General discrete dynamical systems on {-1, 0, +1}

## License

MIT
