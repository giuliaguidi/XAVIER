# Xavier
Xavier: High-Performance X-Drop Adaptive Banded Pairwise Alignment 

## Introduction

## Usage

### Compilation

LOGAN requires CUDA 10 and C++14. To compile LOGAN simply type:
```
make demo_v100
```
LOGAN has been optimized to run on the NVIDIA Tesla V100 (GB), but can run on any NVIDA GPU.
To compile to use other GPUs simply type:
```
make demo
```
This command disables Tesla V100 GPU optmimization flags. 

### Demo

LOGAN generates an executable called **demo**.
To check everything has been compiled properly type:
```
./demo inputs_demo/example.txt 17 21 1
```
This command executes LOGAN on our example dataset with a k-mer length of 17, an X-drop value of 21 using a single GPU.
If everything executes correctly you can start using LOGAN with any input, any X-drop, and any number of GPUs.

The command line inputs are:
```
./demo [input] [k-mer-length] [X-drop] [#GPUS]
```
The input format for this demo is:
```
[seqV] [posV] [seqH] [posH] [strand]
```
**Each line of the input contains a pair of sequences to align**: the query sequence (seqV), the starting position of the seed on the query sequence (posV), the target sequence (seqH), the starting position of the seed on the target sequence (posH), and the relative strand ("c" if on opposite strands, "n" otherwise). Tab separated.

## Performance Analysis

TBD

## Copyright Notice

TBD

## Acknowledgments

Funding provided in part by DOE ASCR through the [Exascale Computing Project](https://www.exascaleproject.org/) and computing provided by [NERSC](https://www.nersc.gov/).

