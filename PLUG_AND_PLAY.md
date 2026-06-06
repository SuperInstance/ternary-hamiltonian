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
use ternary_hamiltonian::TernaryHamiltonian;

let mut h = TernaryHamiltonian::new(3);
h.set_potential(|x| x * x / 2.0);
h.simulate(100);
```

## 🔗 Integration

This crate is part of the [SuperInstance ternary fleet](https://github.com/SuperInstance). It uses the canonical `Ternary` type from `ternary-types` for cross-crate compatibility.

## 📄 License

MIT
