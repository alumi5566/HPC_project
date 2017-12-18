# HPC_project
basic HPC implement (Cache reuse/ Blocking/ MPI )

this is the implement of high performance computing, the trick included:
cache reuse/ register reuse/ blocking/ MPI multi processing
using for different problem included:
matrix multiplication/ GEPP/ SeivingPrimes

environment:
Tardis server provided by UCR
gcc-4.7.2 (mvapich2-1.9a2/gnu-4.6.2 for project_3)
https://kaimingweb.wordpress.com/2016/10/15/tardis-tutorial/

project1:
cache reuse ad reister reuse to accelerate matrix multiplication
notice that the oder we use the loop (ijk/ kij/...) effect the cache we can use

project2:
implement the GEPP decomposite and compare with LAPACK library
use the block version GEPP decomposite and compare again

project3:
parallel sieve of eratosthenes for finding all prime numbers within 10^10

