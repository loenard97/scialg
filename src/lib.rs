//! This crate is a collection of common algorithms used in data science written in pure Rust.
//! It includes modules for
//!  - Evaluation of constants
//!  - Filter functions
//!  - Fourier analysis
//!  - Special function evaluations
//!  - Integration
//!  - Inter- (and extra)polation
//!  - Linear algebra
//!  - Dataset modelling (aka fitting)
//!  - Differential equations
//!  - Optimization
//!  - Root finding
//!  - Sorting
//!  - Statistics
//!  - Vector math

pub mod consts;
pub mod filter;
pub mod fourier;
pub mod function;
pub mod integration;
pub mod interpolation;
pub mod linalg;
pub mod model;
pub mod ode;
pub mod optimize;
pub mod root;
pub mod sort;
pub mod statistic;
pub mod vector;
