#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
//#include "MyMPI.h"
#define MIN(a,b)  ((a)<(b)?(a):(b))


#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1 )
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW((id),p,n))
#define BLOCK_OWNER(index,p,n) (( ((p)*(index)+1)-1 ) / (n) )

int main (int argc, char *argv[]){
   long long int i, n;
   int id, p;
   long long int low_value, high_value;
   long long int first,index,prime;
   long long int size, proc0_size;
   long long int count;      //local sum
   long long int global_count;       //global sum
   double elapsed_time;
   double max_elapsed_time;

   MPI_Init (&argc, &argv);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();
   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);

    if (argc != 2) {
      if (!id) printf ("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
      return 0;
}

   n = atoll(argv[1]);
   low_value = 2 + BLOCK_LOW(id,p,n-1);
   high_value = 2 + BLOCK_HIGH(id,p,n-1);
   size = BLOCK_SIZE(id,p,n-1);
   printf("process %d, [%lld:%lld]:%lld\n",id,low_value,high_value,size);
   proc0_size = (n-1)/p;
   if ((2 + proc0_size) < (int) sqrt((double) n)) {
      if (!id) printf ("Too many processes\n");
      MPI_Finalize();
      return 0;
   }

   char *marked = (char *)malloc(size*sizeof(char));
   if (marked == NULL) {
      printf ("Cannot allocate enough memory\n");
      MPI_Finalize();
      return 0;
   }

   for (i = 0; i < size; i++)
	marked[i] = 0;
   if (!id)
	index = 0;
   prime = 2;
   do {
      if (prime * prime > low_value)
         first = prime * prime - low_value;
      else {
         if (!(low_value % prime)) first = 0;
         else first = prime - (low_value % prime);
      }
      printf("\tfirst = %lld(low:%lld)\n",first,low_value);
      for (i = first; i < size; i += prime) marked[i] = 1;
      if (!id) {
         while (marked[++index]);
         prime = index + 2;
	 printf("\tbroadcast prime:%lld\n",prime);
      }
      MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
   } while (prime * prime <= n);

   count = 0;
   for (i = 0; i < size; i++)
      if (!marked[i]) count++;
   MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,
      0, MPI_COMM_WORLD);
   elapsed_time += MPI_Wtime();
   if (!id) {
      printf ("%lld primes are less than or equal to %lld\n",
         global_count, n);
      printf ("Total elapsed time: %10.6f\n", elapsed_time);
   }
   free(marked);
   MPI_Finalize ();
   return 0;
}



