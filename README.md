# MPI Adaptive Numerical Integration
This project implements an **Adaptive Simpson’s** Rule using MPI to evaluate parallel efficiency through static and dynamic load balancing strategies.

## Compilation
The program is compiled using the Open MPI wrapper for C, utilizing high-level optimization flags for numerical performance.
```Bash
mpicc -O3 -fopenmp -o integration integration.c -lm
```

### Flags:
- `-O3`: Enables aggressive code optimization.
- `-fopenmp`: Supports OpenMP pragmas if hybrid threading is utilized.
- `-lm`: Links the standard math library.

## Execution Examples
The implementation supports two modes: **Static (1)** and **Dynamic (2)**.

### Mode 1: 
Static PartitioningDivide the range $[0, 1]$ into equal segments across all available processes.
```Bash
# Run with 4 processes in Static Mode
mpiexec --oversubscribe -np 4 ./integration 1 1 1e-6
```

### Mode 2: 
Dynamic Load BalancingUses a Master-Worker (Rank 0 as Master) architecture to handle irregular workloads.
```Bash
# Run with 8 processes in Dynamic Mode
mpiexec --oversubscribe -np 8 ./integration 1 2 1e-8
```
**Note:**
- The `--oversubscribe` flag is used to handle environment-specific process limits.
- `mpiexec` can be replaced with `mpirun`

## Machine Specifications
The benchmarks were conducted on the following hardware environment:

|Component | Specification|
| ------------- | ------------- |
|CPU Model | Intel(R) Core(TM) i5-8250U @ 1.60GHz |
|Architecture | x86_64 (4 Cores / 8 Threads) |
|OS Environment | WSL Ubuntu 24.04 |
|Byte Order | Little Endian|

## Software Dependencies
- **MPI Implementation**: Open MPI 5.0.9
- **Compiler**: mpicc (Language: C)
- **Mathematical Constraints**: Requires `<math.h>` for sine and oscillatory functions.