//! Mirrors the "Getting Started" example in README.md.
//!
//! The README is not compiled as a doctest, so without this integration test
//! the documented API could silently drift from the real one (which is exactly
//! what happened: the README previously called a non-existent
//! `stormer_verlet(&phase, &h, 0.2)` signature). If this test fails to compile
//! or panics, the README example must be updated to match.

use ternary_hamiltonian::{
    EnergyConservation, Hamiltonian, PhaseSpace, SymplecticIntegrator, TernaryCoupling,
};

#[test]
fn readme_example_runs() {
    // One degree of freedom: q = +1, p = 0 (a turning point).
    let mut phase = PhaseSpace::new(vec![1], vec![0]);
    let h = Hamiltonian::new(1.0, 1.0); // T = V = 1
    let coupling = TernaryCoupling::harmonic(); // α = β = 1
    let mut tracker = EnergyConservation::new(h.energy_at(&phase));

    for _step in 0..10 {
        phase = SymplecticIntegrator::symplectic_euler(&phase, &coupling);
        tracker.record(h.energy_at(&phase));
        // The Z3 flow must never leave the ternary alphabet {-1, 0, +1}.
        assert!(phase.is_valid());
    }

    // One energy sample is recorded per step.
    assert_eq!(tracker.history.len(), 10);

    // The harmonic flow on Z3 is periodic, so the trajectory must return to the
    // initial state (q = +1, p = 0). For α = β = 1 the symplectic-Euler period
    // is 2: (+1,0) -> (0,+1) -> (+1,0). After 10 (even) steps we are back home.
    assert_eq!(phase.positions, vec![1]);
    assert_eq!(phase.momenta, vec![0]);
    // Identical state => identical energy => zero drift over a full period.
    assert!(tracker.drift() < 1e-12);
}
