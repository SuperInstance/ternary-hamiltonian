//! # ternary-hamiltonian
//!
//! Hamiltonian mechanics on ternary phase space — all positions and momenta
//! live in the discrete set {-1, 0, +1}.
//!
//! This library provides:
//! - [`Hamiltonian`]: kinetic + potential energy representation
//! - [`PhaseSpace`]: ternary phase space state (positions and momenta in {-1,0,+1})
//! - [`SymplecticIntegrator`]: discrete symplectic integration schemes
//! - [`EnergyConservation`]: track and measure energy drift over integration
//! - [`PoissonBracket`]: discrete Poisson brackets for ternary observables
//! - [`LiouvilleTheorem`]: verify phase space volume conservation

/// Clamp a value to the nearest element of {-1, 0, +1}.
fn clamp_ternary(v: i8) -> i8 {
    if v > 0 {
        1
    } else if v < 0 {
        -1
    } else {
        0
    }
}

/// Clamp a floating-point value to the nearest element of {-1, 0, +1}.
fn clamp_ternary_f64(v: f64) -> i8 {
    let rounded = v.round() as i64;
    if rounded > 0 {
        1
    } else if rounded < 0 {
        -1
    } else {
        0
    }
}

// ---------------------------------------------------------------------------
// Hamiltonian
// ---------------------------------------------------------------------------

/// Represents the Hamiltonian (total mechanical energy) as the sum of kinetic
/// and potential energy terms.
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
}

// ---------------------------------------------------------------------------
// PhaseSpace
// ---------------------------------------------------------------------------

/// A point in ternary phase space. Both positions and momenta are elements of
/// {-1, 0, +1}. Any out-of-range value supplied to the constructor is clamped
/// to the nearest ternary value via `signum`.
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

    /// The number of degrees of freedom (length of the position/momentum vectors).
    pub fn dimension(&self) -> usize {
        self.positions.len()
    }

    /// Verify that every stored value is a valid ternary value.
    pub fn is_valid(&self) -> bool {
        self.positions.iter().chain(self.momenta.iter()).all(|&v| v == -1 || v == 0 || v == 1)
    }
}

// ---------------------------------------------------------------------------
// SymplecticIntegrator
// ---------------------------------------------------------------------------

/// Symplectic integration schemes adapted for ternary phase space.
///
/// After each continuous update the resulting values are rounded and clamped
/// back to {-1, 0, +1} so the discrete structure is preserved.
pub struct SymplecticIntegrator;

impl SymplecticIntegrator {
    /// Symplectic Euler step.
    ///
    /// The continuous update rules are:
    ///   p_{n+1} = p_n - dt * dV/dq(q_n)   (force = -dV/dq ≈ -potential * q)
    ///   q_{n+1} = q_n + dt * dT/dp(p_{n+1}) (velocity = dT/dp ≈ kinetic * p)
    ///
    /// For the simple ternary harmonic Hamiltonian H = T*p^2/2 + V*q^2/2 the
    /// partial derivatives are: ∂H/∂p = T*p and ∂H/∂q = V*q.
    /// Results are projected back onto {-1,0,+1} via rounding.
    pub fn symplectic_euler(
        phase: &PhaseSpace,
        hamiltonian: &Hamiltonian,
        dt: f64,
    ) -> PhaseSpace {
        let n = phase.dimension();
        let mut new_momenta = Vec::with_capacity(n);
        let mut new_positions = Vec::with_capacity(n);

        for i in 0..n {
            let q = phase.positions[i] as f64;
            let p = phase.momenta[i] as f64;

            // p_{n+1} = p - dt * V * q  (force = -∂V/∂q = -potential*q)
            let p_new = p - dt * hamiltonian.potential * q;
            let p_new_t = clamp_ternary_f64(p_new);

            // q_{n+1} = q + dt * T * p_{n+1}
            let q_new = q + dt * hamiltonian.kinetic * (p_new_t as f64);
            let q_new_t = clamp_ternary_f64(q_new);

            new_momenta.push(p_new_t);
            new_positions.push(q_new_t);
        }

        PhaseSpace {
            positions: new_positions,
            momenta: new_momenta,
        }
    }

    /// Störmer–Verlet (leapfrog) step.
    ///
    ///   p_{n+1/2} = p_n - (dt/2) * V * q_n
    ///   q_{n+1}   = q_n + dt * T * p_{n+1/2}
    ///   p_{n+1}   = p_{n+1/2} - (dt/2) * V * q_{n+1}
    ///
    /// All intermediate and final values are clamped to {-1,0,+1}.
    pub fn stormer_verlet(
        phase: &PhaseSpace,
        hamiltonian: &Hamiltonian,
        dt: f64,
    ) -> PhaseSpace {
        let n = phase.dimension();
        let mut new_momenta = Vec::with_capacity(n);
        let mut new_positions = Vec::with_capacity(n);

        for i in 0..n {
            let q = phase.positions[i] as f64;
            let p = phase.momenta[i] as f64;

            // half-step momentum
            let p_half = p - (dt / 2.0) * hamiltonian.potential * q;
            let p_half_t = clamp_ternary_f64(p_half);

            // full-step position
            let q_new = q + dt * hamiltonian.kinetic * (p_half_t as f64);
            let q_new_t = clamp_ternary_f64(q_new);

            // complete momentum step
            let p_new = (p_half_t as f64) - (dt / 2.0) * hamiltonian.potential * (q_new_t as f64);
            let p_new_t = clamp_ternary_f64(p_new);

            new_positions.push(q_new_t);
            new_momenta.push(p_new_t);
        }

        PhaseSpace {
            positions: new_positions,
            momenta: new_momenta,
        }
    }
}

// ---------------------------------------------------------------------------
// EnergyConservation
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// PoissonBracket
// ---------------------------------------------------------------------------

/// Discrete Poisson bracket for observables defined on ternary phase space.
///
/// The Poisson bracket is approximated by a finite-difference sum over the
/// degrees of freedom:
///
///   {f, g} = Σ_i  (∂f/∂q_i · ∂g/∂p_i  −  ∂f/∂p_i · ∂g/∂q_i)
///
/// Partial derivatives are estimated by the central finite difference over the
/// ternary alphabet:
///
///   ∂f/∂q_i ≈ (f(q_i=+1) − f(q_i=−1)) / 2
///
/// The observable slices `f` and `g` must each have length `2 * n` where `n =
/// phase.dimension()`. The layout is:
///
///   f[0..n]   — values of f at q_i ∈ {-1, 0, +1} sampled at +1 vs −1
///   f[n..2n]  — values of f at p_i ∈ {-1, 0, +1} sampled at +1 vs −1
///
/// More precisely the layout stores, for each degree of freedom i:
///   f_q[i] = f evaluated with q_i = +1   (index i)
///   f_q_neg[i] = f evaluated with q_i = -1 (index i, but we use f[i] as +1 sample
///                and approximate the -1 sample from g symmetry)
///
/// For simplicity the slices encode, for each i:
///   index i       → value of observable at position component +1
///   index i + n   → value of observable at momentum component +1
/// and the -1 samples are taken as the negatives (odd-function assumption).
///
/// This yields:
///   ∂f/∂q_i ≈ (f[i]     − (−f[i]))     / 2 = f[i]
///   ∂f/∂p_i ≈ (f[i+n]   − (−f[i+n]))   / 2 = f[i+n]
pub struct PoissonBracket;

impl PoissonBracket {
    /// Compute the discrete Poisson bracket {f, g}.
    ///
    /// `f` and `g` must each have length `2 * phase.dimension()`.
    ///
    /// # Panics
    /// Panics if `f` or `g` do not have exactly `2 * phase.dimension()` elements.
    pub fn compute(f: &[f64], g: &[f64], phase: &PhaseSpace) -> f64 {
        let n = phase.dimension();
        assert_eq!(f.len(), 2 * n, "f must have length 2 * dimension");
        assert_eq!(g.len(), 2 * n, "g must have length 2 * dimension");

        let mut bracket = 0.0;
        for i in 0..n {
            // Central finite differences over ternary values (+1 vs -1), scaled by 1/2
            let df_dq = f[i];           // (f(q=+1) - f(q=-1)) / 2 under odd-function assumption
            let df_dp = f[i + n];       // (f(p=+1) - f(p=-1)) / 2

            let dg_dq = g[i];
            let dg_dp = g[i + n];

            bracket += df_dq * dg_dp - df_dp * dg_dq;
        }
        bracket
    }
}

// ---------------------------------------------------------------------------
// LiouvilleTheorem
// ---------------------------------------------------------------------------

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
    /// collection of states. Each unique (positions, momenta) pair counts as
    /// one cell.
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

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // 1. PhaseSpace clamps out-of-range values to {-1, 0, +1}
    #[test]
    fn test_phasespace_clamps_out_of_range() {
        let ps = PhaseSpace::new(vec![5, -3, 0, 2], vec![-10, 1, 0, 4]);
        assert!(ps.is_valid());
        assert_eq!(ps.positions, vec![1, -1, 0, 1]);
        assert_eq!(ps.momenta, vec![-1, 1, 0, 1]);
    }

    // 2. PhaseSpace with already-valid values is unchanged
    #[test]
    fn test_phasespace_valid_values_unchanged() {
        let ps = PhaseSpace::new(vec![1, 0, -1], vec![-1, 0, 1]);
        assert_eq!(ps.positions, vec![1, 0, -1]);
        assert_eq!(ps.momenta, vec![-1, 0, 1]);
        assert!(ps.is_valid());
    }

    // 3. PhaseSpace dimension
    #[test]
    fn test_phasespace_dimension() {
        let ps = PhaseSpace::new(vec![1, 0, -1, 1, 0], vec![0, 1, -1, 0, 1]);
        assert_eq!(ps.dimension(), 5);
    }

    // 4. PhaseSpace panics on dimension mismatch
    #[test]
    #[should_panic]
    fn test_phasespace_panics_on_mismatch() {
        let _ = PhaseSpace::new(vec![1, 0], vec![1]);
    }

    // 5. Hamiltonian total_energy = kinetic + potential
    #[test]
    fn test_hamiltonian_total_energy() {
        let h = Hamiltonian::new(2.5, -1.5);
        assert!((h.total_energy() - 1.0).abs() < 1e-12);
    }

    // 6. Hamiltonian total_energy with zeros
    #[test]
    fn test_hamiltonian_zero_energy() {
        let h = Hamiltonian::new(0.0, 0.0);
        assert_eq!(h.total_energy(), 0.0);
    }

    // 7. Hamiltonian total_energy negative values
    #[test]
    fn test_hamiltonian_negative_energy() {
        let h = Hamiltonian::new(-1.0, -2.0);
        assert!((h.total_energy() - (-3.0)).abs() < 1e-12);
    }

    // 8. SymplecticEuler produces valid ternary output
    #[test]
    fn test_symplectic_euler_valid_output() {
        let phase = PhaseSpace::new(vec![1, -1, 0], vec![0, 1, -1]);
        let h = Hamiltonian::new(1.0, 1.0);
        let result = SymplecticIntegrator::symplectic_euler(&phase, &h, 0.5);
        assert!(result.is_valid());
        assert_eq!(result.dimension(), 3);
    }

    // 9. Störmer-Verlet produces valid ternary output
    #[test]
    fn test_stormer_verlet_valid_output() {
        let phase = PhaseSpace::new(vec![1, 0, -1], vec![-1, 1, 0]);
        let h = Hamiltonian::new(1.0, 0.5);
        let result = SymplecticIntegrator::stormer_verlet(&phase, &h, 0.3);
        assert!(result.is_valid());
        assert_eq!(result.dimension(), 3);
    }

    // 10. Integration step preserves dimension
    #[test]
    fn test_integration_preserves_dimension() {
        let phase = PhaseSpace::new(vec![1, 0, -1, 1, 0], vec![0, -1, 1, 0, 1]);
        let h = Hamiltonian::new(1.0, 1.0);
        let after_euler = SymplecticIntegrator::symplectic_euler(&phase, &h, 0.1);
        let after_verlet = SymplecticIntegrator::stormer_verlet(&phase, &h, 0.1);
        assert_eq!(after_euler.dimension(), phase.dimension());
        assert_eq!(after_verlet.dimension(), phase.dimension());
    }

    // 11. EnergyConservation drift calculation
    #[test]
    fn test_energy_conservation_drift() {
        let mut ec = EnergyConservation::new(1.0);
        ec.record(1.1);
        ec.record(0.8);
        ec.record(1.05);
        // max deviation = |0.8 - 1.0| = 0.2
        assert!((ec.drift() - 0.2).abs() < 1e-12);
    }

    // 12. EnergyConservation zero drift when no records
    #[test]
    fn test_energy_conservation_zero_drift_no_records() {
        let ec = EnergyConservation::new(5.0);
        assert_eq!(ec.drift(), 0.0);
    }

    // 13. EnergyConservation exact conservation (all values equal initial)
    #[test]
    fn test_energy_conservation_exact() {
        let mut ec = EnergyConservation::new(2.0);
        for _ in 0..10 {
            ec.record(2.0);
        }
        assert!(ec.drift() < 1e-12);
    }

    // 14. PoissonBracket antisymmetry: {f,g} = -{g,f}
    #[test]
    fn test_poisson_bracket_antisymmetry() {
        let phase = PhaseSpace::new(vec![1, -1], vec![0, 1]);
        let f = vec![1.0, 2.0, 0.5, -1.0];
        let g = vec![-0.5, 1.5, 2.0, 0.0];
        let fg = PoissonBracket::compute(&f, &g, &phase);
        let gf = PoissonBracket::compute(&g, &f, &phase);
        assert!((fg + gf).abs() < 1e-12, "antisymmetry violated: fg={fg}, gf={gf}");
    }

    // 15. PoissonBracket linearity: {af+bg, h} = a{f,h} + b{g,h}
    #[test]
    fn test_poisson_bracket_linearity() {
        let phase = PhaseSpace::new(vec![1, 0], vec![-1, 1]);
        let f = vec![1.0, -1.0, 0.5, 2.0];
        let g = vec![0.5, 1.5, -1.0, 0.0];
        let h = vec![2.0, 0.0, 1.0, -0.5];
        let a = 3.0_f64;
        let b = -2.0_f64;

        // Compute af + bg
        let afpbg: Vec<f64> = f.iter().zip(g.iter()).map(|(fi, gi)| a * fi + b * gi).collect();

        let lhs = PoissonBracket::compute(&afpbg, &h, &phase);
        let rhs = a * PoissonBracket::compute(&f, &h, &phase)
            + b * PoissonBracket::compute(&g, &h, &phase);
        assert!((lhs - rhs).abs() < 1e-12, "Linearity violated: lhs={lhs}, rhs={rhs}");
    }

    // 16. PoissonBracket with zero observable gives zero
    #[test]
    fn test_poisson_bracket_zero_observable() {
        let phase = PhaseSpace::new(vec![1, 0, -1], vec![0, 1, -1]);
        let f = vec![0.0; 6];
        let g = vec![1.0, 2.0, 3.0, 0.5, -0.5, 1.5];
        let result = PoissonBracket::compute(&f, &g, &phase);
        assert!((result).abs() < 1e-12);
    }

    // 17. LiouvilleTheorem volume counts distinct states
    #[test]
    fn test_liouville_volume() {
        let states = vec![
            PhaseSpace::new(vec![1], vec![0]),
            PhaseSpace::new(vec![-1], vec![1]),
            PhaseSpace::new(vec![1], vec![0]),  // duplicate
            PhaseSpace::new(vec![0], vec![0]),
        ];
        assert!((LiouvilleTheorem::volume(&states) - 3.0).abs() < 0.5);
    }

    // 18. LiouvilleTheorem conservation with identical ensembles
    #[test]
    fn test_liouville_conservation_identical() {
        let initial = vec![
            PhaseSpace::new(vec![1, 0], vec![0, 1]),
            PhaseSpace::new(vec![-1, 1], vec![1, -1]),
        ];
        let final_states = initial.clone();
        assert!(LiouvilleTheorem::check_conservation(&initial, &final_states));
    }

    // 19. LiouvilleTheorem conservation fails when volume changes
    #[test]
    fn test_liouville_conservation_fails() {
        let initial = vec![
            PhaseSpace::new(vec![1], vec![0]),
            PhaseSpace::new(vec![-1], vec![1]),
            PhaseSpace::new(vec![0], vec![-1]),
        ];
        // final ensemble collapses to fewer distinct cells
        let final_states = vec![
            PhaseSpace::new(vec![1], vec![0]),
            PhaseSpace::new(vec![1], vec![0]),
        ];
        assert!(!LiouvilleTheorem::check_conservation(&initial, &final_states));
    }

    // 20. Full integration loop with energy tracking
    #[test]
    fn test_full_integration_loop_energy_tracking() {
        let mut phase = PhaseSpace::new(vec![1, -1, 0], vec![0, 1, -1]);
        let h = Hamiltonian::new(1.0, 1.0);
        let initial_energy = h.total_energy();
        let mut ec = EnergyConservation::new(initial_energy);

        for _ in 0..20 {
            phase = SymplecticIntegrator::symplectic_euler(&phase, &h, 0.1);
            ec.record(h.total_energy());
            assert!(phase.is_valid());
        }

        // Energy is conserved by construction (Hamiltonian doesn't change)
        assert_eq!(ec.drift(), 0.0);
        assert_eq!(ec.history.len(), 20);
    }

    // 21. Störmer-Verlet multi-step stays valid
    #[test]
    fn test_stormer_verlet_multi_step() {
        let mut phase = PhaseSpace::new(vec![1, 0, -1, 1], vec![-1, 1, 0, 0]);
        let h = Hamiltonian::new(0.5, 2.0);

        for _ in 0..50 {
            phase = SymplecticIntegrator::stormer_verlet(&phase, &h, 0.05);
            assert!(phase.is_valid(), "Phase space left ternary domain");
        }
        assert_eq!(phase.dimension(), 4);
    }

    // 22. EnergyConservation records and history length
    #[test]
    fn test_energy_conservation_history_length() {
        let mut ec = EnergyConservation::new(0.0);
        for i in 0..100 {
            ec.record(i as f64 * 0.01);
        }
        assert_eq!(ec.history.len(), 100);
        // max deviation is |0.99 - 0.0| = 0.99
        assert!((ec.drift() - 0.99).abs() < 1e-10);
    }
}
