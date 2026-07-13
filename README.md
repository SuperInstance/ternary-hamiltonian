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

so that $\partial H / \partial p = T \, p$ and $\partial H / \partial q = V \, q$.

> **Note on discreteness (Z₃, not clamping).** The ternary alphabet $\{-1, 0, +1\}$ is *not* treated as a subset of $\mathbb{R}$ to be rounded/clamped to after each continuous step. Such clamping is a many-to-one map: it collapses distinct phase-space points onto the same ternary value and would destroy phase-space volume (Liouville's theorem would fail). Instead, ternary values are identified with $\mathbb{Z}_3 = \{0, 1, 2\}$ via $-1 \mapsto 0,\; 0 \mapsto 1,\; +1 \mapsto 2$, and **every update is performed modulo 3**. An update $x \leftarrow x + c \pmod 3$ is a cyclic rotation — a bijection on the finite state space — so the composition of updates is again a bijection and phase-space volume is preserved exactly by construction. No rounding or clamping occurs during the dynamics.

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

All updates are performed in $\mathbb{Z}_3$ (modular arithmetic); because each line is a $\mathbb{Z}_3$ rotation, every step is a bijection and hence exactly volume-preserving (symplectic).

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
use ternary_hamiltonian::{
    Hamiltonian, PhaseSpace, SymplecticIntegrator, TernaryCoupling, EnergyConservation,
};

fn main() {
    // One degree of freedom: q = +1, p = 0 (a turning point)
    let mut phase = PhaseSpace::new(vec![1], vec![0]);
    let h = Hamiltonian::new(1.0, 1.0);          // T = V = 1
    let coupling = TernaryCoupling::harmonic();  // α = β = 1
    let mut tracker = EnergyConservation::new(h.energy_at(&phase));

    for step in 0..10 {
        phase = SymplecticIntegrator::symplectic_euler(&phase, &coupling);
        tracker.record(h.energy_at(&phase));
        println!("step {}: q = {}, p = {}", step, phase.positions[0], phase.momenta[0]);
    }

    println!("max energy drift = {}", tracker.drift());
}
```

For $\alpha = \beta = 1$ the $\mathbb{Z}_3$ symplectic-Euler flow is periodic with period 2: $(q,p) = (+1,0) \to (0,+1) \to (+1,0)$. The state therefore oscillates forever inside $\{-1,0,+1\}$ and, because it returns to the initial point, the recorded energy drift over a full period is exactly zero. (Note: this exact periodicity is the real conservation statement here — see *Physics checks* below. The Euclidean energy $H = \tfrac{1}{2}(T p^2 + V q^2)$ is **not** a step-wise invariant of the $\mathbb{Z}_3$ rotation in general; what is conserved is phase-space volume, and the flow is exactly periodic.)

> The snippet above is kept in sync with the code by `tests/readme_example.rs`, which `cargo test` compiles and runs.

## Running the Tests

Run the full 30-test unit suite (plus 5 doc-tests) with `cargo test`. Each test verifies a specific physical or structural property:

| Test group | What it verifies |
|------------|------------------|
| `test_z3_*` (4 tests) | $\mathbb{Z}_3$ encode/decode round-trips; `add`/`sub`/`mul` are correct mod 3; `add` is the inverse of `sub`. |
| `test_phasespace_*` (4 tests) | Out-of-range inputs are clamped to $\{-1,0,+1\}$, dimensions match, mismatched vectors panic. |
| `test_hamiltonian_*` (3 tests) | Total energy $H = T + V$ for positive, zero, and negative inputs; `energy_at` matches $H = \tfrac{1}{2}\sum_i(T p_i^2 + V q_i^2)$. |
| `test_symplectic_euler_valid_output` | The Euler step never leaves the ternary domain. |
| `test_stormer_verlet_valid_output` | The leapfrog step never leaves the ternary domain. |
| `test_integration_preserves_dimension` | Both integrators conserve the number of degrees of freedom. |
| `test_stormer_verlet_multi_step` | 50 consecutive leapfrog steps remain inside $\{-1,0,+1\}$. |
| `test_phase_space_volume_preservation_*` (3 tests) | **Liouville's theorem:** evolving the *entire* phase space (all $3^2=9$ states for 1 DOF, all $3^4=81$ for 2 DOF) through many steps never reduces the count of distinct occupied cells — the $\mathbb{Z}_3$ map is a permutation. |
| `test_z3_verlet_is_periodic` | The Verlet map is invertible: every state returns to itself within 100 steps (a true permutation is periodic). |
| `test_energy_conservation_*` (4 tests) | Drift is computed as $\max_n |E_n - E_0|$; exact conservation yields zero drift; history length is tracked correctly. |
| `test_full_integration_loop_energy_tracking` | A 20-step Euler loop stays valid and records 20 energy samples. |
| `test_poisson_bracket_*` (3 tests) | Antisymmetry $\{f,g\} = -\{g,f\}$, linearity $\{af+bg,h\} = a\{f,h\} + b\{g,h\}$, and zero observable yields zero bracket. |
| `test_liouville_*` (3 tests) | Distinct-state counting, conservation check passes for identical ensembles and fails for collapsed ones. |

### Physics checks

The headline physics claim of this crate is **not** that the Euclidean energy $H = \tfrac{1}{2}(T p^2 + V q^2)$ is conserved step-by-step (it is not, in general, under a $\mathbb{Z}_3$ rotation). The real claims are:

1. **Volume preservation (Liouville).** Each integrator step is a $\mathbb{Z}_3$ translation, i.e. a cyclic permutation of the finite state space. The composition of permutations is a permutation, so the number of distinct occupied cells is exactly invariant. Verified exhaustively over the full 9- and 81-state spaces (`test_phase_space_volume_preservation_*`). Equivalently, for 1 DOF the symplectic-Euler update matrix is $\begin{pmatrix}1-\alpha\beta & \beta \\ -\alpha & 1\end{pmatrix}$ with determinant $(1-\alpha\beta) - (-\alpha\beta) = 1 \pmod 3$, so the linear map is volume-preserving.
2. **Exact periodicity.** A permutation of a finite set is invertible, so every trajectory is periodic. The full state space is recycled without loss or duplication.

The energy $H$ is exposed as a *diagnostic* over the (periodic) trajectory: it is bounded (the alphabet is finite) and returns to its initial value over a full period.

## Related crates

- [`ternary-noether`](https://github.com/SuperInstance/ternary-noether) — Noether's theorem, symmetries → conservation laws
- [`ternary-electromagnetism`](https://github.com/SuperInstance/ternary-electromagnetism) — Maxwell's equations and Yee-lattice EM waves
- [`ternary-phase`](https://github.com/SuperInstance/ternary-phase) — Phase-space geometry and ternary state portraits
- [`ternary-dynamics`](https://github.com/SuperInstance/ternary-dynamics) — General discrete dynamical systems on {-1, 0, +1}

## License

MIT
