module purge
module load gcc-4.7.2

g++ -o $1 $1.c -I/opt/lapack/include /opt/lapack/lib/liblapacke.a /opt/lapack/lib/liblapack.a /opt/lapack/lib/librefblas.a -lgfortran -lm
