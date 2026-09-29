//! Shared, platform-independent analysis code for clinker.
//!
//! This crate deliberately accepts GenBank bytes rather than filesystem paths,
//! so both the native CLI and the browser/WASM frontend can use it.

mod model;
mod parse_genbank;

pub use model::{Cluster, Gene, Locus};
pub use parse_genbank::{ParseError, parse_genbank};
