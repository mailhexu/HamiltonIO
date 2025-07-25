# HamiltonIO Product Overview

HamiltonIO is a Python library for reading and writing Hamiltonian matrices from density functional theory (DFT) codes. It provides a unified interface for handling tight-binding-like Hamiltonians with localized basis sets.

## Core Purpose
- IO operations for Hamiltonian files from various DFT codes (SIESTA, Wannier90, ABACUS, etc.)
- Unified API for working with electronic structure data
- Support for different spin configurations (non-polarized, collinear, non-collinear/SOC)
- k-space and real-space Hamiltonian transformations

## Target Users
- Computational physicists and chemists
- Electronic structure researchers
- Materials science developers working with DFT calculations

## Key Features
- Abstract Hamiltonian base class with consistent interface
- Support for orthogonal and non-orthogonal basis sets
- Spin-orbit coupling (SOC) handling
- k-point sampling and eigenvalue calculations
- Integration with scientific Python ecosystem (NumPy, SciPy, ASE)

## Development Status
Currently in Alpha stage (version 0.2.5) with active development focused on expanding DFT code support and improving the API.