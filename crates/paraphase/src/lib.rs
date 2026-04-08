#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_possible_wrap)]

// Core functionality
pub mod assembly;
pub mod phaser;

// Paraphase-specific utilities
pub mod detail;
pub mod io;

// BAM utilities
pub mod depth;
pub mod realign;

// Configuration
pub mod config;
