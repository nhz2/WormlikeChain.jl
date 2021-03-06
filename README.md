# WormlikeChain WIP
Started by Nathan Zimmerberg on Mar 6, 2021

**Authors** Nathan Zimmerberg

Latest Revision: Mar 6, 2021


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://nhz2.github.io/WormlikeChain.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nhz2.github.io/WormlikeChain.jl/dev)
[![Build Status](https://github.com/nhz2/WormlikeChain.jl/workflows/CI/badge.svg)](https://github.com/nhz2/WormlikeChain.jl/actions)
[![Coverage](https://codecov.io/gh/nhz2/WormlikeChain.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/nhz2/WormlikeChain.jl)


A project to simulate a twisting discrete wormlike chain model, specifically plectoneme formation.

## Simulation rough design specifications and scope

1. dWLC gas phase (implicit solvent) where non bonded interactions are relatively uncommon.
2. Less than 10,000 particles.
3. 1000 particle system should be able to step in <10 micro seconds of computer time.
    1. The full system step must happen "on chip".
    2. Each simulation is compiled before being run.
    3. Written in Julia for both CPU and GPU.
5. Double precision positions.
6. Point masses, not rigid bodies or rigid constraints.
7. Fine tunable Langevin dynamics and optional monte carlo steps.
10. Deterministic simulation when run on the same hardware.
11. Additive force field.
    1. Bonds.
    2. Finite range non bonded with exclusions.
    3. External.

## Secondary Goals
1. 100 particle system step in <1 micro seconds of computer time.
2. Fancy periodic boundary conditions.
3. Good user interface for setting up simulations.
4. Good interface for modifying simulation behavior/ adding new features. 
5. Adding support for constraints/rigid bodies.

## Rough Todo list
1. Figure out what the system step function should look like to get the required performance for some toy systems.
    1. Ideal gas
        1. Get a stateless RNG function to get gaussian random numbers.
        2. Store system positions and velocities.
    2. SHO
        1. Calculate and store forces and potentials.
    2. SHO with monte carlo step
        1. Use energy and temp to accept or reject moves.
    3. Single chain with harmonic length bonds.
    3. Single chain with harmonic length and angle bonds.
    1. Two chains with harmonic length and angle bonds, and shifted 6-12 potential
        1. Neighbor list.

2. Implement simple twistable wormlike chain model for DNA without mismatch. And compare to other DNA models and experiments.


## Data Structures

????

## Functions

????