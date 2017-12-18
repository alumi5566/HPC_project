all:
	g++ -o callMypack callMypack.c -I/opt/lapack/include /opt/lapack/lib/liblapacke.a /opt/lapack/lib/liblapack.a /opt/lapack/lib/librefblas.a -lgfortran -lm
	g++ -o callpack callpack.c -I/opt/lapack/include /opt/lapack/lib/liblapacke.a /opt/lapack/lib/liblapack.a /opt/lapack/lib/librefblas.a -lgfortran -lm -O3
clean:
	rm -f callMypack
	rm -f callpack

