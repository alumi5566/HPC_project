#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stdio.h>
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
   long long int global_count=0;       //global sum
   int node;
   double elapsed_time;
   double max_elapsed_time;
   
   long long int low_value_0, high_value_0;
   MPI_Init (&argc, &argv);
   MPI_Barrier(MPI_COMM_WORLD);
   //elapsed_time = -MPI_Wtime();
   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);

    if (argc != 3) {
      if (!id) printf ("Command line: %s <m> <node>\n", argv[0]);
        MPI_Finalize();
      return 0;
}

   n = atoll(argv[1]);
   node = atoi(argv[2]);
   low_value = 2 + BLOCK_LOW(id,p,n-1);
   high_value = 2 + BLOCK_HIGH(id,p,n-1);
   low_value_0 = 2 + BLOCK_LOW(0,p,n-1);
   high_value_0 = 2 + BLOCK_HIGH(0,p,n-1);

   size = BLOCK_SIZE(id,p,n-1);
   //printf("process %d, [%lld:%lld]:%lld\n",id,low_value,high_value,size);
   proc0_size = (n-1)/p;
   if ((2 + proc0_size) < (int) sqrt((double) n)) {
      if (!id) printf ("Too many processes\n");
      MPI_Finalize();
      return 0;
   }

   //char *marked = (char *)malloc(size*sizeof(char)); //part1
   char *marked = (char *)malloc((size/2)*sizeof(char));
   char *marked_prime = (char *)malloc((sqrt(n))*sizeof(char));
   if (marked == NULL) {
      printf ("Cannot allocate enough memory\n");
      MPI_Finalize();
      return 0;
   }


   //for (i = 0; i < size; i++) //part1
   for (i = 0; i < (size/2); i++)
	marked[i] = 0;
   //marked[0]=1; // part1: value 1 is in mark[0]
#if 0
   if (!id)
	index = 0;
#endif
   for(i=0; i< (sqrt(n)); i++)
	marked_prime[i] = 0;
   if(1)
	index = 0;
   //prime = 2; //part1
   prime = 3;
   long long int tmp;
   elapsed_time = -MPI_Wtime();
   do {
      if (prime * prime > low_value)
         first = prime * prime - low_value;
      else {
         if (!(low_value % prime)) first = 0;
         else first = prime - (low_value % prime);
      }
     //printf("\tfirst = %lld(low:%lld)\n",first,low_value);
     //printf("\tfirst = %d\n",first);
     //for (i = first; i < size; i += prime) //part1
     //    marked[i] = 1; 
     for (i = first; i < size; i += prime){
	    //if(id == 1)
		//printf("\tmark:%d\n",i+low_value);
	    tmp = i+low_value;
	    if((tmp)%2==0){
	    	//printf("\tits an even!!\n");
            }else{
                marked[(tmp-low_value)/2]=1;
		//if(id == 7)
		//	printf("\tmark:%lld\n",tmp);
	    }
     }
#if 0 
     if(prime * prime > low_value_0)
	first = prime * prime - low_value;
     else{
         if (!(low_value_0 % prime)) first = 0;
         else first = prime - (low_value_0 % prime);
     }
#endif
     long long int first_0 = prime*prime - low_value_0;
     for(i = first_0; i < (sqrt(n)); i += prime){
            tmp = i+low_value_0;
            if((tmp)%2==0){
                //printf("\tits an even!!\n");
                }else{
                    marked_prime[(tmp-low_value_0)/2]=1;
                }
     }
#if 0
      if (!id) {
         while (marked[++index]);
         //prime = index + 2; //part1
	 prime = (index*2) + 3;
	 printf("\tbroadcast prime:%lld\n",prime);
      }
      MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
      if(1){
	 while(marked_prime[++index]);
	 prime = (index*2) + 3;
      }
      //printf("\tp = %lld\n",prime);
   } while (prime * prime <= n);

   count = 0;
   //for (i = 0; i < size; i++) //part1
   for (i = 0; i < size; i++){
      //if(!id) printf("[%d]:%d ",i+low_value,marked[i]);
      //if(i%10 == 0) printf("\n");
      tmp = i+low_value;
          if((tmp%2) == 0){
      }else{
          if(!marked[(tmp-low_value)/2]){
		//if(id == 7)
		//	printf("\tcnt mark:%lld\n",tmp);
		count++;
	  }
	  //if(id == 7)
	  //	printf("\tcnt mark:%lld\n",tmp);
      }
   }
   //printf("cnt = %lld\n",count);
   MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,
      0, MPI_COMM_WORLD);
   elapsed_time += MPI_Wtime();
   if (!id) {
      //printf ("%lld primes are less than or equal to %lld\n",global_count, n);
      //printf ("Total elapsed time: %10.6f\n", elapsed_time);
      if(node == 1)
        printf("Part2:\n");
      printf("The total number of prime: %lld, total time: %10.6f, total node %d\n",global_count+1,elapsed_time,node);
   }
   free(marked);
   free(marked_prime);
   MPI_Finalize ();
   return 0;
}



