//! # ternary-hamiltonian
//!
//! Hamiltonian mechanics on ternary phase space — all positions and momenta
//! live in the discrete set {-1, 0, +1}, treated as elements of Z₃.
//!
//! ## Fix: Z₃ Symplectic Dynamics
//!
//! **Original bug:** The integrators used clamping (round-then-truncate via
//! `clamp_ternary_f64`) to project continuous Hamiltonian flow back onto
//! {-1, 0, +1}. Clamping is a many-to-one map — it collapses distinct phase
//! space points onto the same ternary value, destroying the symplectic 2-form
//! and violating Liouville's theorem. The core claim of a symplectic integrator
//! on ternary phase space was invalid.
//!
//! **Fix:** Ternary values are now elements of Z₃ = {0, 1, 2} (via the
//! bijection −1↦0, 0↦1, +1↦2). All dynamics use modular arithmetic in Z₃.
//! Each integration step is a Z₃ translation (addition mod 3), which is a
//! cyclic permutation — a bijection on the finite state space. The composition
//! of permutations is a permutation, so phase space volume is preserved exactly.
//! No clamping occurs during dynamics; the symplectic structure is intact.
//!
//! **Key insight:** ternary {-1, 0, +1} is a rotation group (Z₃), not a subset
//! of ℝ to be truncated to. Treating it as Z₃ makes every dynamical step a
//! rotation, which is automatically volume-preserving.

// ─── Z₃ Arithmetic ────────────────────────────────────────────────────────

/// Modular arithmetic over Z₃ = {0, 1, 2}.
///
/// Ternary phase space values {-1, 0, +1} are mapped to Z₃ via −1↦0, 0↦1, +1↦2.
/// All operations are performed mod 3 so the result always lies in Z₃.
pub mod z3 {
    /// Map a ternary value {−1, 0, +1} to Z₃ = {0, 1, 2}.
    #[inline]
    pub fn encode(v: i8) -> u8 {
        ((v as i32 + 4) % 3) as u8
    }

    /// Map a Z₃ value {0, 1, 2} back to ternary {−1, 0, +1}.
    #[inline]
    pub fn decode(z: u8) -> i8 {
        (z as i8) - 1
    }

    /// Add in Z₃.
    #[inline]
    pub fn add(a: u8, b: u8) -> u8 {
        ((a as u32 + b as u32) % 3) as u8
    }

    /// Subtract in Z₃: a − b mod 3.
    #[inline]
    pub fn sub(a: u8, b: u8) -> u8 {
        ((a as i32 + 3 - b as i32) % 3) as u8
    }

    /// Multiply in Z₃.
    #[inline]
    pub fn mul(a: u8, b: u8) -> u8 {
        ((a as u32 * b as u32) % 3) as u8
    }
}

// ─── TernaryCoupling ──────────────────────────────────────────────────────

/// Coupling constants for Z₃ Hamiltonian dynamics.
///
/// Each constant is an element of Z₃ = {0, 1, 2}:
/// - `alpha` (force coupling): how strongly position q drives momentum p
/// - `beta` (velocity coupling): how strongly momentum p drives position q
///
/// The discrete Hamiltonian flow is:
///   p ← p − α·q   (mod 3)
///   q ← q + β·p   (mod 3)
///
/// # Example
/// ```
/// use ternary_hamiltonian::TernaryCoupling;
/// let c = TernaryCoupling::new(1, 1);  // isotropic harmonic
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TernaryCoupling {
    pub alpha: u8,
    pub beta: u8,
}

impl TernaryCoupling {
    /// Create coupling constants. Values are reduced mod 3.
    pub fn new(alpha: u8, beta: u8) -> Self {
        Self {
            alpha: alpha % 3,
            beta: beta % 3,
        }
    }

    /// Isotropic harmonic coupling: α = β = 1.
    pub fn harmonic() -> Self {
        Self { alpha: 1, beta: 1 }
    }
}

// ─── PhaseSpace ───────────────────────────────────────────────────────────

/// A point in ternary phase space. Both positions and momenta are elements of
/// {-1, 0, +1}. Out-of-range values supplied to the constructor are rounded
/// to the nearest ternary value (input validation only; dynamics never clamp).
///
/// # Example
/// ```
/// use ternary_hamiltonian::PhaseSpace;
/// let ps = PhaseSpace::new(vec![1, 0, -1], vec![0, 1, -1]);
/// assert_eq!(ps.dimension(), 3);
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct PhaseSpace {
    pub positions: Vec<i8>,
    pub momenta: Vec<i8>,
}

/// Clamp a value to the nearest element of {-1, 0, +1} (input validation only).
fn clamp_ternary(v: i8) -> i8 {
    v.signum()
}

impl PhaseSpace {
    /// Construct a PhaseSpace, clamping all values to {-1, 0, +1}.
    ///
    /// # Panics
    /// Panics if `positions` and `momenta` have different lengths.
    pub fn new(positions: Vec<i8>, momenta: Vec<i8>) -> Self {
        assert_eq!(
            positions.len(),
            momenta.len(),
            "positions and momenta must have the same length"
        );
        Self {
            positions: positions.into_iter().map(clamp_ternary).collect(),
            momenta: momenta.into_iter().map(clamp_ternary).collect(),
        }
    }

    /// The number of degrees of freedom.
    pub fn dimension(&self) -> usize {
        self.positions.len()
    }

    /// Verify that every stored value is a valid ternary value.
    pub fn is_valid(&self) -> bool {
        self.positions
            .iter()
            .chain(self.momenta.iter())
            .all(|&v| v == -1 || v == 0 || v == 1)
    }

    /// Encode to Z₃ representation: positions and momenta as Vec<u8> in {0,1,2}.
    fn encode_z3(&self) -> (Vec<u8>, Vec<u8>) {
        let q = self.positions.iter().map(|&v| z3::encode(v)).collect();
        let p = self.momenta.iter().map(|&v| z3::encode(v)).collect();
        (q, p)
    }

    /// Decode from Z₃ representation.
    fn decode_z3(q: Vec<u8>, p: Vec<u8>) -> Self {
        Self {
            positions: q.into_iter().map(z3::decode).collect(),
            momenta: p.into_iter().map(z3::decode).collect(),
        }
    }
}

// ─── Hamiltonian ──────────────────────────────────────────────────────────

/// Represents the Hamiltonian (total mechanical energy) as kinetic + potential.
///
/// # Example
/// ```
/// use ternary_hamiltonian::Hamiltonian;
/// let h = Hamiltonian::new(1.5, -0.5);
/// assert_eq!(h.total_energy(), 1.0);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Hamiltonian {
    pub kinetic: f64,
    pub potential: f64,
}

impl Hamiltonian {
    /// Create a new Hamiltonian with the given kinetic and potential energies.
    pub fn new(kinetic: f64, potential: f64) -> Self {
        Self { kinetic, potential }
    }

    /// Return the total energy H = T + V.
    pub fn total_energy(&self) -> f64 {
        self.kinetic + self.potential
    }

    /// Compute the Hamiltonian energy at a specific phase space point.
    ///
    /// H(q, p) = Σᵢ (T·pᵢ² + V·qᵢ²) / 2
    ///
    /// Since qᵢ, pᵢ ∈ {−1, 0, +1}, their squares are in {0, 1}.
    pub fn energy_at(&self, phase: &PhaseSpace) -> f64 {
        let mut h = 0.0;
        for i in 0..phase.dimension() {
            let q2 = (phase.positions[i] as f64).powi(2);
            let p2 = (phase.momenta[i] as f64).powi(2);
            h += self.kinetic * p2 + self.potential * q2;
        }
        h / 2.0
    }
}

// ─── SymplecticIntegrator ─────────────────────────────────────────────────

/// Symplectic integration over Z₃ ternary phase space.
///
/// All updates use modular arithmetic in Z₃, ensuring each step is a
/// permutation of the finite state space. Phase space volume is preserved
/// by construction (Liouville's theorem holds exactly).
pub struct SymplecticIntegrator;

impl SymplecticIntegrator {
    /// Symplectic Euler step over Z₃.
    ///
    ///   p_{n+1} = p_n − α·q_n      (mod 3)
    ///   q_{n+1} = q_n + β·p_{n+1}  (mod 3)
    ///
    /// Each line is a Z₃ translation (cyclic permutation), so the map is
    /// bijective and preserves phase space volume.
    pub fn symplectic_euler(phase: &PhaseSpace, coupling: &TernaryCoupling) -> PhaseSpace {
        let (q, mut p) = phase.encode_z3();
        let n = phase.dimension();

        // Update momentum: p -= α·q  (mod 3)
        for i in 0..n {
            p[i] = z3::sub(p[i], z3::mul(coupling.alpha, q[i]));
        }

        // Update position: q += β·p_new  (mod 3, using updated momentum)
        let mut q_new = q;
        for i in 0..n {
            q_new[i] = z3::add(q_new[i], z3::mul(coupling.beta, p[i]));
        }

        PhaseSpace::decode_z3(q_new, p)
    }

    /// Störmer–Verlet (leapfrog) step over Z₃.
    ///
    ///   p_{½}   = p_n − α·q_n        (mod 3)
    ///   q_{n+1} = q_n + β·p_{½}      (mod 3)
    ///   p_{n+1} = p_{½} − α·q_{n+1}  (mod 3)
    ///
    /// This is second-order accurate (in the discrete sense) and exactly
    /// symplectic: the composition of three Z₃ permutations is a permutation.
    pub fn stormer_verlet(phase: &PhaseSpace, coupling: &TernaryCoupling) -> PhaseSpace {
        let (q, mut p) = phase.encode_z3();
        let n = phase.dimension();

        // Half-step momentum: p -= α·q  (mod 3)
        for i in 0..n {
            p[i] = z3::sub(p[i], z3::mul(coupling.alpha, q[i]));
        }

        // Full-step position: q += β·p_half  (mod 3)
        let mut q_new = q;
        for i in 0..n {
            q_new[i] = z3::add(q_new[i], z3::mul(coupling.beta, p[i]));
        }

        // Complete momentum step: p -= α·q_new  (mod 3)
        for i in 0..n {
            p[i] = z3::sub(p[i], z3::mul(coupling.alpha, q_new[i]));
        }

        PhaseSpace::decode_z3(q_new, p)
    }
}

// ─── EnergyConservation ──────────────────────────────────────────────────

/// Track total energy values over time and measure drift from the initial value.
///
/// # Example
/// ```
/// use ternary_hamiltonian::EnergyConservation;
/// let mut ec = EnergyConservation::new(1.0);
/// ec.record(1.1);
/// ec.record(0.9);
/// assert!((ec.drift() - 0.1).abs() < 1e-12);
/// ```
#[derive(Debug, Clone)]
pub struct EnergyConservation {
    pub initial_energy: f64,
    pub history: Vec<f64>,
}

impl EnergyConservation {
    /// Create a new tracker with the given initial energy.
    pub fn new(initial: f64) -> Self {
        Self {
            initial_energy: initial,
            history: Vec::new(),
        }
    }

    /// Record a new energy measurement.
    pub fn record(&mut self, energy: f64) {
        self.history.push(energy);
    }

    /// Return the maximum absolute deviation from the initial energy across
    /// all recorded samples. Returns `0.0` if no samples have been recorded.
    pub fn drift(&self) -> f64 {
        self.history
            .iter()
            .map(|&e| (e - self.initial_energy).abs())
            .fold(0.0_f64, f64::max)
    }
}

// ─── PoissonBracket ──────────────────────────────────────────────────────

/// Discrete Poisson bracket for observables defined on ternary phase space.
///
/// The Poisson bracket is approximated by a finite-difference sum over the
/// degrees of freedom:
///
///   {f, g} = Σᵢ (∂f/∂qᵢ · ∂g/∂pᵢ − ∂f/∂pᵢ · ∂g/∂qᵢ)
///
/// Partial derivatives are estimated by the central finite difference over the
/// ternary alphabet:
///
///   ∂f/∂qᵢ ≈ (f(qᵢ=+1) − f(qᵢ=−1)) / 2
///
/// The observable slices `f` and `g` must each have length `2 * n` where `n =
/// phase.dimension()`. The layout is:
///
///   index i       → value of observable at position component +1
///   index i + n   → value of observable at momentum component +1
///
/// Under the odd-function assumption, ∂f/∂qᵢ = f[i] and ∂f/∂pᵢ = f[i+n].
pub struct PoissonBracket;

impl PoissonBracket {
    /// Compute the discrete Poisson bracket {f, g}.
    ///
    /// # Panics
    /// Panics if `f` or `g` do not have exactly `2 * phase.dimension()` elements.
    pub fn compute(f: &[f64], g: &[f64], phase: &PhaseSpace) -> f64 {
        let n = phase.dimension();
        assert_eq!(f.len(), 2 * n, "f must have length 2 * dimension");
        assert_eq!(g.len(), 2 * n, "g must have length 2 * dimension");

        let mut bracket = 0.0;
        for i in 0..n {
            let df_dq = f[i];
            let df_dp = f[i + n];
            let dg_dq = g[i];
            let dg_dp = g[i + n];
            bracket += df_dq * dg_dp - df_dp * dg_dq;
        }
        bracket
    }
}

// ─── LiouvilleTheorem ────────────────────────────────────────────────────

/// Verify discrete Liouville's theorem: the number of distinct occupied cells
/// in ternary phase space is preserved by Hamiltonian flow.
///
/// # Example
/// ```
/// use ternary_hamiltonian::{LiouvilleTheorem, PhaseSpace};
/// let states = vec![
///     PhaseSpace::new(vec![1], vec![0]),
///     PhaseSpace::new(vec![-1], vec![1]),
///     PhaseSpace::new(vec![1], vec![0]),  // duplicate
/// ];
/// assert_eq!(LiouvilleTheorem::volume(&states), 2.0);
/// ```
pub struct LiouvilleTheorem;

impl LiouvilleTheorem {
    /// Count the number of distinct phase space cells occupied by the given
    /// collection of states.
    pub fn volume(states: &[PhaseSpace]) -> f64 {
        use std::collections::HashSet;
        let unique: HashSet<PhaseSpace> = states.iter().cloned().collect();
        unique.len() as f64
    }

    /// Check that the number of distinct occupied cells is preserved between
    /// an initial and a final ensemble of states.
    pub fn check_conservation(initial: &[PhaseSpace], final_states: &[PhaseSpace]) -> bool {
        let v_init = Self::volume(initial);
        let v_final = Self::volume(final_states);
        (v_init - v_final).abs() < 0.5
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Z₃ Arithmetic ────────────────────────────────────────────────

    #[test]
    fn test_z3_encode_decode_roundtrip() {
        for v in [-1_i8, 0, 1] {
            assert_eq!(z3::decode(z3::encode(v)), v);
        }
    }

    #[test]
    fn test_z3_add() {
        // 0+0=0, 0+1=1, 0+2=2, 1+1=2, 1+2=0, 2+2=1
        assert_eq!(z3::add(0, 0), 0);
        assert_eq!(z3::add(1, 2), 0);
        assert_eq!(z3::add(2, 2), 1);
    }

    #[test]
    fn test_z3_sub_inverse_of_add() {
        for a in 0..3u8 {
            for b in 0..3u8 {
                assert_eq!(z3::add(z3::sub(a, b), b), a);
            }
        }
    }

    #[test]
    fn test_z3_mul() {
        assert_eq!(z3::mul(0, 2), 0);
        assert_eq!(z3::mul(1, 2), 2);
        assert_eq!(z3::mul(2, 2), 1); // 4 mod 3 = 1
    }

    // ── PhaseSpace ────────────────────────────────────────────────────

    #[test]
    fn test_phasespace_clamps_out_of_range() {
        let ps = PhaseSpace::new(vec![5, -3, 0, 2], vec![-10, 1, 0, 4]);
        assert!(ps.is_valid());
        assert_eq!(ps.positions, vec![1, -1, 0, 1]);
        assert_eq!(ps.momenta, vec![-1, 1, 0, 1]);
    }

    #[test]
    fn test_phasespace_valid_values_unchanged() {
        let ps = PhaseSpace::new(vec![1, 0, -1], vec![-1, 0, 1]);
        assert_eq!(ps.positions, vec![1, 0, -1]);
        assert_eq!(ps.momenta, vec![-1, 0, 1]);
        assert!(ps.is_valid());
    }

    #[test]
    fn test_phasespace_dimension() {
        let ps = PhaseSpace::new(vec![1, 0, -1, 1, 0], vec![0, 1, -1, 0, 1]);
        assert_eq!(ps.dimension(), 5);
    }

    #[test]
    #[should_panic]
    fn test_phasespace_panics_on_mismatch() {
        let _ = PhaseSpace::new(vec![1, 0], vec![1]);
    }

    // ── Hamiltonian ────────────────────────────────────────────────────

    #[test]
    fn test_hamiltonian_total_energy() {
        let h = Hamiltonian::new(2.5, -1.5);
        assert!((h.total_energy() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_hamiltonian_zero_energy() {
        let h = Hamiltonian::new(0.0, 0.0);
        assert_eq!(h.total_energy(), 0.0);
    }

    #[test]
    fn test_hamiltonian_energy_at() {
        // H = T*p² + V*q²) / 2, with T=2, V=3, q=1, p=-1
        // H = (2*1 + 3*1) / 2 = 2.5
        let h = Hamiltonian::new(2.0, 3.0);
        let ps = PhaseSpace::new(vec![1], vec![-1]);
        assert!((h.energy_at(&ps) - 2.5).abs() < 1e-12);
    }

    // ── SymplecticIntegrator ──────────────────────────────────────────

    #[test]
    fn test_symplectic_euler_valid_output() {
        let phase = PhaseSpace::new(vec![1, -1, 0], vec![0, 1, -1]);
        let coupling = TernaryCoupling::harmonic();
        let result = SymplecticIntegrator::symplectic_euler(&phase, &coupling);
        assert!(result.is_valid());
        assert_eq!(result.dimension(), 3);
    }

    #[test]
    fn test_stormer_verlet_valid_output() {
        let phase = PhaseSpace::new(vec![1, 0, -1], vec![-1, 1, 0]);
        let coupling = TernaryCoupling::new(1, 2);
        let result = SymplecticIntegrator::stormer_verlet(&phase, &coupling);
        assert!(result.is_valid());
        assert_eq!(result.dimension(), 3);
    }

    #[test]
    fn test_integration_preserves_dimension() {
        let phase = PhaseSpace::new(vec![1, 0, -1, 1, 0], vec![0, -1, 1, 0, 1]);
        let coupling = TernaryCoupling::harmonic();
        let after_euler = SymplecticIntegrator::symplectic_euler(&phase, &coupling);
        let after_verlet = SymplecticIntegrator::stormer_verlet(&phase, &coupling);
        assert_eq!(after_euler.dimension(), phase.dimension());
        assert_eq!(after_verlet.dimension(), phase.dimension());
    }

    #[test]
    fn test_stormer_verlet_multi_step() {
        let mut phase = PhaseSpace::new(vec![1, 0, -1, 1], vec![-1, 1, 0, 0]);
        let coupling = TernaryCoupling::new(1, 2);

        for _ in 0..50 {
            phase = SymplecticIntegrator::stormer_verlet(&phase, &coupling);
            assert!(phase.is_valid(), "Phase space left ternary domain");
        }
        assert_eq!(phase.dimension(), 4);
    }

    // ── Phase Space Volume Preservation (Liouville's Theorem) ──────────

    /// The critical test: evolve the ENTIRE phase space (all 9 states for
    /// 1 degree of freedom) through many steps. If any two states collide
    /// (map to the same point), volume decreases. If the map is a true
    /// permutation, volume stays at 9.0 forever.
    #[test]
    fn test_phase_space_volume_preservation_euler() {
        let coupling = TernaryCoupling::harmonic();
        let values = [-1_i8, 0, 1];

        // Enumerate all 3² = 9 states for 1 DOF
        let mut ensemble: Vec<PhaseSpace> = Vec::new();
        for &q in &values {
            for &p in &values {
                ensemble.push(PhaseSpace::new(vec![q], vec![p]));
            }
        }
        assert_eq!(LiouvilleTheorem::volume(&ensemble), 9.0);

        // Evolve through many steps
        for step in 0..30 {
            ensemble = ensemble
                .iter()
                .map(|ps| SymplecticIntegrator::symplectic_euler(ps, &coupling))
                .collect();

            let vol = LiouvilleTheorem::volume(&ensemble);
            assert_eq!(
                vol,
                9.0,
                "Volume lost at Euler step {}: got {} (expected 9.0)",
                step + 1,
                vol
            );
        }
    }

    #[test]
    fn test_phase_space_volume_preservation_verlet() {
        let coupling = TernaryCoupling::new(1, 2);
        let values = [-1_i8, 0, 1];

        let mut ensemble: Vec<PhaseSpace> = Vec::new();
        for &q in &values {
            for &p in &values {
                ensemble.push(PhaseSpace::new(vec![q], vec![p]));
            }
        }
        assert_eq!(LiouvilleTheorem::volume(&ensemble), 9.0);

        for step in 0..30 {
            ensemble = ensemble
                .iter()
                .map(|ps| SymplecticIntegrator::stormer_verlet(ps, &coupling))
                .collect();

            let vol = LiouvilleTheorem::volume(&ensemble);
            assert_eq!(
                vol,
                9.0,
                "Volume lost at Verlet step {}: got {} (expected 9.0)",
                step + 1,
                vol
            );
        }
    }

    /// Multi-DOF volume preservation: all 3^4 = 81 states for 2 DOFs.
    #[test]
    fn test_phase_space_volume_preservation_2dof() {
        let coupling = TernaryCoupling::harmonic();
        let values = [-1_i8, 0, 1];

        let mut ensemble: Vec<PhaseSpace> = Vec::new();
        for &q0 in &values {
            for &q1 in &values {
                for &p0 in &values {
                    for &p1 in &values {
                        ensemble.push(PhaseSpace::new(vec![q0, q1], vec![p0, p1]));
                    }
                }
            }
        }
        assert_eq!(LiouvilleTheorem::volume(&ensemble), 81.0);

        for step in 0..20 {
            ensemble = ensemble
                .iter()
                .map(|ps| SymplecticIntegrator::stormer_verlet(ps, &coupling))
                .collect();

            let vol = LiouvilleTheorem::volume(&ensemble);
            assert_eq!(
                vol,
                81.0,
                "Volume lost at 2-DOF Verlet step {}: got {}",
                step + 1,
                vol
            );
        }
    }

    /// The Z₃ map is a permutation, so it must be invertible. Verify that
    /// enough steps bring every state back to itself (periodicity).
    #[test]
    fn test_z3_verlet_is_periodic() {
        let coupling = TernaryCoupling::harmonic();
        let start = PhaseSpace::new(vec![1], vec![0]);

        let mut state = start.clone();
        for _ in 0..100 {
            state = SymplecticIntegrator::stormer_verlet(&state, &coupling);
            if state == start {
                return; // Found the period — map is invertible
            }
        }
        panic!("State did not return to initial value in 100 steps; map may not be a permutation");
    }

    // ── EnergyConservation ────────────────────────────────────────────

    #[test]
    fn test_energy_conservation_drift() {
        let mut ec = EnergyConservation::new(1.0);
        ec.record(1.1);
        ec.record(0.8);
        ec.record(1.05);
        assert!((ec.drift() - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_energy_conservation_zero_drift_no_records() {
        let ec = EnergyConservation::new(5.0);
        assert_eq!(ec.drift(), 0.0);
    }

    #[test]
    fn test_energy_conservation_exact() {
        let mut ec = EnergyConservation::new(2.0);
        for _ in 0..10 {
            ec.record(2.0);
        }
        assert!(ec.drift() < 1e-12);
    }

    #[test]
    fn test_energy_conservation_history_length() {
        let mut ec = EnergyConservation::new(0.0);
        for i in 0..100 {
            ec.record(i as f64 * 0.01);
        }
        assert_eq!(ec.history.len(), 100);
        assert!((ec.drift() - 0.99).abs() < 1e-10);
    }

    // ── Full integration loop with energy tracking ────────────────────

    #[test]
    fn test_full_integration_loop_energy_tracking() {
        let mut phase = PhaseSpace::new(vec![1, -1, 0], vec![0, 1, -1]);
        let h = Hamiltonian::new(1.0, 1.0);
        let coupling = TernaryCoupling::harmonic();
        let initial_energy = h.energy_at(&phase);
        let mut ec = EnergyConservation::new(initial_energy);

        for _ in 0..20 {
            phase = SymplecticIntegrator::symplectic_euler(&phase, &coupling);
            ec.record(h.energy_at(&phase));
            assert!(phase.is_valid());
        }

        // Z₃ dynamics is periodic, so energy should return to initial value
        // at some point. We just verify the tracker works.
        assert_eq!(ec.history.len(), 20);
    }

    // ── PoissonBracket ────────────────────────────────────────────────

    #[test]
    fn test_poisson_bracket_antisymmetry() {
        let phase = PhaseSpace::new(vec![1, -1], vec![0, 1]);
        let f = vec![1.0, 2.0, 0.5, -1.0];
        let g = vec![-0.5, 1.5, 2.0, 0.0];
        let fg = PoissonBracket::compute(&f, &g, &phase);
        let gf = PoissonBracket::compute(&g, &f, &phase);
        assert!(
            (fg + gf).abs() < 1e-12,
            "antisymmetry violated: fg={fg}, gf={gf}"
        );
    }

    #[test]
    fn test_poisson_bracket_linearity() {
        let phase = PhaseSpace::new(vec![1, 0], vec![-1, 1]);
        let f = vec![1.0, -1.0, 0.5, 2.0];
        let g = vec![0.5, 1.5, -1.0, 0.0];
        let h = vec![2.0, 0.0, 1.0, -0.5];
        let a = 3.0_f64;
        let b = -2.0_f64;

        let afpbg: Vec<f64> = f
            .iter()
            .zip(g.iter())
            .map(|(fi, gi)| a * fi + b * gi)
            .collect();

        let lhs = PoissonBracket::compute(&afpbg, &h, &phase);
        let rhs = a * PoissonBracket::compute(&f, &h, &phase)
            + b * PoissonBracket::compute(&g, &h, &phase);
        assert!(
            (lhs - rhs).abs() < 1e-12,
            "Linearity violated: lhs={lhs}, rhs={rhs}"
        );
    }

    #[test]
    fn test_poisson_bracket_zero_observable() {
        let phase = PhaseSpace::new(vec![1, 0, -1], vec![0, 1, -1]);
        let f = vec![0.0; 6];
        let g = vec![1.0, 2.0, 3.0, 0.5, -0.5, 1.5];
        let result = PoissonBracket::compute(&f, &g, &phase);
        assert!(result.abs() < 1e-12);
    }

    // ── LiouvilleTheorem ──────────────────────────────────────────────

    #[test]
    fn test_liouville_volume() {
        let states = vec![
            PhaseSpace::new(vec![1], vec![0]),
            PhaseSpace::new(vec![-1], vec![1]),
            PhaseSpace::new(vec![1], vec![0]), // duplicate
            PhaseSpace::new(vec![0], vec![0]),
        ];
        assert!((LiouvilleTheorem::volume(&states) - 3.0).abs() < 0.5);
    }

    #[test]
    fn test_liouville_conservation_identical() {
        let initial = vec![
            PhaseSpace::new(vec![1, 0], vec![0, 1]),
            PhaseSpace::new(vec![-1, 1], vec![1, -1]),
        ];
        let final_states = initial.clone();
        assert!(LiouvilleTheorem::check_conservation(
            &initial,
            &final_states
        ));
    }

    #[test]
    fn test_liouville_conservation_fails() {
        let initial = vec![
            PhaseSpace::new(vec![1], vec![0]),
            PhaseSpace::new(vec![-1], vec![1]),
            PhaseSpace::new(vec![0], vec![-1]),
        ];
        let final_states = vec![
            PhaseSpace::new(vec![1], vec![0]),
            PhaseSpace::new(vec![1], vec![0]),
        ];
        assert!(!LiouvilleTheorem::check_conservation(
            &initial,
            &final_states
        ));
    }
}
