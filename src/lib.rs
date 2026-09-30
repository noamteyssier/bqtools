#![allow(clippy::module_inception)]
// The crate is a CLI; its pub API only exists so `bqdoc` can reach `cli`.
#![allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::must_use_candidate,
    clippy::iter_without_into_iter,
    clippy::len_without_is_empty
)]

pub mod cli;
pub mod commands;
pub mod types;

#[cfg(test)]
mod testutils;
