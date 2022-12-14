# qrot
Demo of tennis racket theorem.

## Overview
Ground-up library-free standard `C++` program that simulates free rigid body motion, using quaternions to apply rotations. The inertia matrix is recalculated in lab frame and Cholesky factorized each time step. The demo script `qrot-demo.sh` compiles the code, runs three simulations spinning a set of particles around three different (perturbed) axes but with same angular velocity, and then produces the plots stored in the `readme-figures` folder.

## Dzhanibekov effect
The stability of the different spins is illustrated by the following figure.

![intermediate axis theorem](/readme-figures/qrot-plot-theta.png) 

The orientation of the angular velocity as seen in the body frame flips periodically. This is only true for one of the simulations (axis 2). The other two simulations show stable rotation (axis 1, axis 3).

![energy conservation](/readme-figures/qrot-plot-K.png)

The kinetic energy is conserved in all simulations. The unstable simulation spins around the so-called intermediate axis (axis 2, middle energy). The angular momentum is also conserved (see additional figures in folder).

The angular velocity vector has an interesting evolution for the unstable simulation.

![energy conservation](/readme-figures/qrot-plot-omega-axis2.png)

## Browser animation export
The `--x3d-file` option (with related configurations such as `--x3d-template`) allows saving snapshots to a file that can be animated in a web browser. Here is an [example](https://raw.githack.com/olofer/qrot/main/x3dom-generated/qrot-anim.html) showing the unstable (axis 2) spin graphically.

## Things to improve
Among other things, the integrator should be replaced with a discrete variational technique, or similar. The text file import is implemented and tested but not demonstrated..

