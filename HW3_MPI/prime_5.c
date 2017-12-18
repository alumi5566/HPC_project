#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stdio.h>
//#include "MyMPI.h"
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MAX4(x, y, s, t) MAX( MAX(x,y), MAX(s,t) )

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1 )
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW((id),p,n))
#define BLOCK_OWNER(index,p,n) (( ((p)*(index)+1)-1 ) / (n) )



int main (int argc, char *argv[]){
    long long int i, n;
    int id, p;
    long long int low_value, high_value;
    long long int first,index,prime;
    long long int first2,prime2;
    long long int first3,prime3;
    long long int first4,prime4;
    
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
        index = 3;
    //prime = 2; //part1
    prime = 3;
    prime2 = 5;
    prime3 = 7;
    prime4 = 11;
    register long long int tmp;
    int Flag_prime = 0;
    register long long int r = prime;
    register long long int r2 = prime2;
    register long long int r3 = prime3;
    register long long int r4 = prime4;
elapsed_time = -MPI_Wtime();
    do {
        
	if (r * r > low_value)
            first = r * r - low_value;
        else {
            if (!(low_value % r)) first = 0;
            else first = r - (low_value % r);
        }

if(!Flag_prime){
        if (r2 * r2 > low_value)
            first2 = r2 * r2 - low_value;
        else {
            if (!(low_value % r2)) first2 = 0;
            else first2 = r2 - (low_value % r2);
        }
        if (r3 * r3 > low_value)
            first3 = r3 * r3 - low_value;
        else {
            if (!(low_value % r3)) first3 = 0;
            else first3 = r3 - (low_value % r3);
        }
        
        if (r4 * r4 > low_value)
            first4 = r4 * r4 - low_value;
        else {
            if (!(low_value % r4)) first4 = 0;
            else first4 = r4 - (low_value % r4);
        }
}

	//int Flag_prime = 0;
	long long int j = 0;
	long long int piv = size/100;
	i = size;
	register long long int B = 14000;
	//for(i = piv+size%100; i <= size; i += piv ){
if(!Flag_prime){
	for(i = 0; i< size; i+=B){
	//while(i<=size){
            for (j = first; j < MIN(size,i+B); j += r){
		tmp = j+low_value;
            	    if((tmp)%2==0){
		    }else{
			marked[(tmp-low_value)/2]=1;
		    }
	    }
	    first = j;
            for (j = first2; j < MIN(size,i+B); j += r2){
                tmp = j+low_value;
                    if((tmp)%2==0){
                    }else{
                        marked[(tmp-low_value)/2]=1;
                    }
            }
	    first2 = j;
            for (j = first3; j < MIN(size,i+B); j += r3){
                tmp = j+low_value;
                    if((tmp)%2==0){
                    }else{
                        marked[(tmp-low_value)/2]=1;
                    }
            }
	    first3 = j;
            for (j = first4; j < MIN(size,i+B); j += r4){
                tmp = j+low_value;
                    if((tmp)%2==0){
                    }else{
                        marked[(tmp-low_value)/2]=1;
                    }
            }
	    first4 = j;
	
	//if(size%2 == 1)	i+=(size/2+1);
	//else	i+=size/2;
	}// 
}else{
            for (j = first; j < size; j += r){
                tmp = j+low_value;
                    if((tmp)%2==0){
                    }else{
                        marked[(tmp-low_value)/2]=1;
                    }
            }

}

        long long int first_0 = r*r - low_value_0;
        for(i = first_0; i < (sqrt(n)); i += r){
            tmp = i+low_value_0;
            if((tmp)%2==0){
                //printf("\tits an even!!\n");
            }else{
                marked_prime[(tmp-low_value_0)/2]=1;
            }
        }
if(!Flag_prime){
	first_0 = r2*r2 - low_value_0;
        for(i = first_0; i < (sqrt(n)); i += r2){
            tmp = i+low_value_0;
            if((tmp)%2==0){
                //printf("\tits an even!!\n");
            }else{
                marked_prime[(tmp-low_value_0)/2]=1;
            }
        }
	first_0 = r3*r3 - low_value_0;
        for(i = first_0; i < (sqrt(n)); i += r3){
            tmp = i+low_value_0;
            if((tmp)%2==0){
                //printf("\tits an even!!\n");
            }else{
                marked_prime[(tmp-low_value_0)/2]=1;
            }
        }
	first_0 = r4*r4 - low_value_0;
        for(i = first_0; i < (sqrt(n)); i += r4){
            tmp = i+low_value_0;
            if((tmp)%2==0){
                //printf("\tits an even!!\n");
            }else{
                marked_prime[(tmp-low_value_0)/2]=1;
            }
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
	//int Flag_prime = 0;
            while(marked_prime[++index] && index<sqrt(n));
	    r = (index*2) + 3;
if(r>=1000000)
	Flag_prime = 0;
if(!Flag_prime){
            while(marked_prime[++index] && index<sqrt(n));
            r2 = (index*2) + 3;
	    if(r2*r2>n)
                r2 = r;
            while(marked_prime[++index] && index<sqrt(n) );
            r3 = (index*2) + 3;
	    if(r3*r3>n)
                r3 = r;
            while(marked_prime[++index] && index<sqrt(n));
	    r4 = (index*2) + 3;
	    tmp = r4;
	    if(r4*r4>n)
		r4 = r;
	    //printf("p1 = %lld, p2 = %lld, p3 = %lld, p4 = %lld\n", prime, prime2, prime3, prime4);
}
    //} while (prime * prime <= n);
    } while ( r * r <= n);
    
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
        }
    }
    //printf("\tid: %d cnt = %lld\n",id, count);
    MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,
                0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    if (!id) {
        //printf ("%lld primes are less than or equal to %lld\n",global_count, n);
        //printf ("Total elapsed time: %10.6f\n", elapsed_time);
      if(node == 1)
        printf("Part3:\n");
      printf("The total number of prime: %lld, total time: %10.6f, total node %d\n",global_count+1,elapsed_time,node);

    }
    free(marked);
    free(marked_prime);
    MPI_Finalize ();
    return 0;
}


