# FermionCommute

## Overview

FermionCommute is a small C++ program that builds on the symbolic_operators library to compute commutators and expectation-value expansions of fermionic (and bosonic) operator expressions using Wick's theorem. 
The program supports model-specific definitions (Hamiltonians, contraction templates, operator bases and symmetries) via a `DefinitionsBase`-style interface and ships with an example implementation for the half-filled extended Hubbard model. 

## Features

- Constructs and manipulates symbolic operator Terms (creation/annihilation operators, coefficients, sums, Kronecker deltas, etc.) via the `mrock::symbolic_operators` library. 
- Computes commutators and nested commutators (e.g. $[H, A]$ and $[B, [H, A]]$). 
- Applies Wick's theorem to expand higher-order operator products into Wick operators using user-provided templates, and simplifies results with symmetry classes (spin, inversion, phase, ...). 
- Serializes computed Wick-term collections to disk (binary Boost archives by default). 

## Requirements

- A C++20-capable compiler (the library and examples target C++20). 
- Boost (iostream and serialization) for saving/loading results; the program uses Boost binary archives for output by default.
- The symbolic_operators library (must be installed to a standard CMake location or to `~/usr/local/` for plug-and-play)
- (Recommended for plug-and-play) CMake 3.30 or newer for building and running tests. 

## Building

If the plug-and-play requirements are fulfilled, you can simply type
```
make build/FermionCommute
```

Running
```
make
```
without any arguments will build the program and then run it for all implemented models.

## Running

The executable expects at least two command-line arguments: an execution type and a model identifier. The documented syntax is:
```
./build/FermionCommute <XP/std> <model>
```
and the program prints a usage message when arguments are missing. 

- Execution types:
  - `XP` - run with the "XP" basis (uses `XP_basis()` from the model definition).  
  - `std` - run with the standard basis (`STD_basis()`). 
The meaning of `XP` is laid out in the context of the iEoM, but in short, the XP basis contains only Hermitian and anti-Hermitian operators, while no such restrictions exist for the `STD` basis.

- Model identifiers shipped with the program include: `hubbard`, `continuum`, `hubbard_dispersions`, `lattice_cut`. The model selector returns a concrete `DefinitionsBase` implementation for the chosen model. 

Example:
```
./build/main XP hubbard
```
This runs the code generation for the (half-filled) extended Hubbard model using the XP basis (see https://doi.org/10.1103/PhysRevB.109.205153). 

## What the program does internally (high-level workflow)

1. It loads the chosen model via `get_model(model_type)`, which instantiates a class like `Hubbard` that supplies Hamiltonian terms, templates and bases.   
2. The program queries the model for its Hamiltonian (`hamiltonian()`), Wick operator templates (`templates()`), the operator basis (`XP_basis()` or `STD_basis()`), and symmetry objects (`symmetries()`).  
3. For each pair of basis vectors the program computes the commutators $[H, \text{basis[i]}]$ and $[\text{basis}_j^\dagger, [H, \text{basis}_i]]$.  
4. The raw results are cleaned (normal ordering, removal of trivial sums/multiplicities) and then transformed into a collection of `WickTerm` objects (representing expectation values) by `wicks_theorem(...)`, using the model-supplied Wick operator templates. 
5. The resulting `WickTermCollector` is post-processed: specific coefficient manipulations or channel rearrangements (model-dependent) are applied, and symmetry simplifications (e.g., spin symmetry, inversion symmetry, phase symmetry) are applied via `clean_wicks(...)`. 
6. The results are serialized to disk as binary boost archives in `../commutators/<model_subfolder>/` by default. The program constructs filenames such as `wick_M_j_i.bin` and `wick_N_j_i.bin`.

## Overview of the files

- `FermionCommute.cpp` - main driver that orchestrates model loading, commutator evaluation, Wick expansion, cleaning and serialization.  
- `Definitions/<model>.hpp` / `Definitions/<model>.cpp` - implementation of various models
  - a Hamiltonian built from `Term` and `Operator` objects;   
  - a set of Wick operator templates used to identify number, SC (superconducting), eta, and CDW operator types;   
  - XP and STD operator bases used to generate matrices \(M\) and \(N\);   
  - model-specific symmetries (inversion, phase symmetry for specified operator types) and a `get_subfolder()` function for output layout.   
- `symbolic_operators` library headers - core types such as `Term`, `Operator`, `Coefficient`, `Momentum`, `WickOperator`, `WickTerm`, `WickTermCollector`, `WickOperatorTemplate`, and symmetry classes live in namespace `mrock::symbolic_operators`.

## Extending / Adding a new model

To add a new model, implement a class derived from the project's `DefinitionsBase` interface and provide at least the following methods, analogous to the included `Hubbard` example:

- `std::vector<Term> hamiltonian() const` - return the Hamiltonian as a vector of `Term` objects.
- `std::vector<WickOperatorTemplate> templates() const` - return templates that define how ordinary operators map to Wick operator types for your problem. 
- `std::vector<std::vector<Term>> XP_basis() const` and `STD_basis() const` - define your operator bases used for computing matrix elements. 
- `std::vector<std::unique_ptr<WickSymmetry>> symmetries() const` - return symmetry objects to be applied during `clean_wicks(...)`.
- `std::string get_subfolder() const` - optional: subfolder name for serialization output.

See the existing files for concrete examples of all these functions.

## Output and serialization

- The program serializes `WickTermCollector` objects into binary files using Boost serialization; filenames include an `M` or `N` marker plus indices derived from the basis loop.
- By default, files are written under `../commutators/<model_subfolder>/` where `model_subfolder` is returned by the model's `get_subfolder()` method (for the Hubbard example, this is `"hubbard/"`).

## Debugging and tests

- Running with `test` (internally supported by the example main) triggers smaller test workflows and prints intermediate results for verification. 
- The source contains helper functions and examples (e.g., a `boson_test()` function) that illustrate bosonic operator handling and commutator behavior.

## Related publications

The following implementations of models have been used in publications:

- `Hubbard`: 
	- Collective excitations in competing phases in two and three dimensions, J. Althüser & G. S. Uhrig, 
 		Physical Review B 109, 205153 (2024), https://doi.org/10.1103/PhysRevB.109.205153
- `Continuum`: 
	- Collective modes in superconductors including Coulomb repulsion, J. Althüser & G. S. Uhrig, 
 		SciPost Physics 19, 067 (2025) https://doi.org/10.21468/SciPostPhys.19.3.067
- `LatticeCut`: 
	- Enhanced Superconductivity in Proximity to Peaks in Densities of States, J. Althüser, I. M. Eremin & G. S. Uhrig,
 		(preprint) arXiv:2512.11451 (2025) https://doi.org/10.48550/arXiv.2512.11451
 	- Secondary Collective Excitations in Intermediate to Strong-Coupling Superconductors, J. Althüser & G. S. Uhrig, 
 		(preprint) arXiv:2605.20059 (2026) https://doi.org/10.48550/arXiv.2605.20059
