# Stokes Constant Calculator

Numerical computation of the Stokes constant in the He-Cu scattering model with a Morse potential

## Overview

This repository contains code for a high-precision numerical study of the Stokes constant in a Hamiltonian system modeling the scattering of helium (He) atoms off a copper (Cu) surface.
This work supports evidence that the Stokes constant is non-zero, validating chaotic behavior in the scattering model by calculating the splitting of invariant manifolds. 
This project uses arbitrary-precision arithmetic to ensure reliability in exponentially small computations.

## Key Features

Arbitrary Precision Calculations: Handles exponentially small terms with arbitrary precision using heyoka and mpmath.
Custom Polynomial Classes: Implements TrigPol, ExpPol, and CustomComplex for specialized handling of trigonometric, exponential polynomials, and complex numbers.
Extrapolation Techniques: Minimizes error terms through extrapolation to achieve greater numerical accuracy in Stokes constant estimation.

## Installation

To use this repository, clone it and install the required dependencies:

```bash
# Clone the repository
git clone https://github.com/AnnaDaSa/StokesConstant_calculator.git

# Navigate into the repository
cd StokesConstant_calculator

# Install dependencies
pip install -r requirements.txt
```

## Usage

Here’s a quick example of how to use this code:

```bash
# First:
# Change the precision and the parameters as wanted in constants.py

# Then the Stokes constants can be calculated by
python calculation_f1.py
```

## Acknowledgments

This repository is based on Anna Danot’s master’s thesis on the numerical study of the He-Cu scattering model, supervised by Prof. Pau Martín and Prof. Mercè Ollé at the Universitat Politècnica de Catalunya.
