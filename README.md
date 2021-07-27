# CholUp
This repository contains a basic implementation 
for the method described in ["Factor Once: Reusing Cholesky Factorizations on Sub-Meshes"](https://people.inf.ethz.ch/herholzp/SigAsia2018/files/HerholzSigAsia2018.pdf),
presented at SIGGRAPH Asia 2018.

## Installation
The code is depends only on [Eigen](https://gitlab.com/libeigen/eigen). Best performance will be achieved when linking to a fast parallel BLAS implementation like Intel MKL. The code can also use Eigen to mimic BLAS. Compile with the option `-DUSE_EIGEN_FOR_BLAS` to select this option.

## Acknowledgements
Large parts of the code a heavily inspired by [Cholmod](https://github.com/DrTimothyAldenDavis/SuiteSparse) authored by Tim Davis and William Hager.
