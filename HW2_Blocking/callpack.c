

/* Calling DGELS using row-major order */
#include <stdlib.h>
#include <stdio.h>
//#include </home/cyu059/HW2/include/lapacke.h>
#include "sys/time.h"
#include <iostream>
#include "lapacke.h"
#include "blas.h"

#define Matrix_size 1600
#define LDA_Size 3
#define LDB_Size 3

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define swap(x,y) { x = x + y; y = x - y; x = x - y; }

#define debugFlag	0

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}


int mySwap(double *a, int n){
    //only use for n*n matrix
    int i,j;
    for(i=0;i<n;i++){
	for(j=0;j<n;j++){
	    if(i>j){
		swap(a[i*n+j],a[j*n+i]);
	    }else{
		//do nothing
	    }
	}
    }
    return 0;

}

int mydgetrf(double *a, int n, int *IPIV){
    int i,j;
    int t,s,k;
    int maxInd;
    double max;
    int z;
#if debugFlag
    print_matrix((char*)"in mydgetrf",n,n,a,n);
#endif
    for(i=0;i<n-1;i++){
	//PIVOT
	maxInd = i;
	max = abs(a[i*n+i]);
	for(t=i+1;t<n;t++){
	    if(abs(a[t*n+i])>max){
		maxInd = t;
		max = abs(a[t*n+i]);
	    }
	}
	//printf("need swap, maxInd=%d\n",maxInd);
	if(max == 0){
	    printf("LU goes wrong");
	}else{
	    if(maxInd != i){
		//save PIVOT info
		//swap(IPIV[i],IPIV[maxInd]);
		IPIV[i]=IPIV[maxInd];
		//swap rows
		for(k=0;k<n;k++)
		    swap(a[i*n+k],a[maxInd*n+k]);
	    }
	}//else
	//print_matrix((char*)"after swap",n,n,a,n);	
#if debugFlag
        //printf("Pivot:");
        //for(z=0;z<n;z++)
        //    printf("[%d]%d",z,IPIV[z]);
#endif	
	for(j=i+1;j<n;j++){
	    a[j*n+i] = a[j*n+i]/a[i*n+i];
		for(k=i+1;k<n;k++)
		    a[j*n+k]-=a[j*n+i]*a[i*n+k];
	}
	//print_matrix((char*)"after factorization",n,n,a,n);
        //printf("\n\n");
    }
    return 0;
}

int myBlockdgetrf(double *a, int n, int *IPIV){

    int ib,end;
    int b=2;
    int i,j,t,k;
    int maxInd;
    double max; 
    for(ib=0;ib<n;ib+=b){
	end = MIN(ib+b-1,n);
#if debugFlag
	printf("handle block: ib=%d,end=%d\n",ib,end);
#endif
	for(i=ib;i<=end;i++){
	    //find pivot row
	    maxInd = i;
            max = abs(a[i*n+i]);
	    for(t=i+1;t<n;t++){
		if(abs(a[t*n+i])>max){
            	    maxInd = t;
            	    max = abs(a[t*n+i]);
                }//need swap
	    }
	    if(max == 0){
            printf("LU goes wrong");
            }else{
                if(maxInd != i){
                    //save PIVOT info
                    //swap(IPIV[i],IPIV[maxInd]);
                    IPIV[i]=IPIV[maxInd];
                    //swap rows
                    for(k=0;k<n;k++)
                    	swap(a[i*n+k],a[maxInd*n+k]);
		}
     	    }//else
#if debugFlag
	    print_matrix((char*)"[Block] after swap(1)(2)",n,n,a,n);  
#endif
	    //register double t = a[i*n+i];
	    for(j=i+1;j<n;j++){
            	a[j*n+i] = a[j*n+i]/a[i*n+i];
		//register double s = a[j*n+i];
                for(k=i+1;k<=end;k++){
                    //for(k=i+1;k<n;k++){ was gonna do the whole row
		    a[j*n+k]-=a[j*n+i]*a[i*n+k];
		}
            }
#if debugFlag
	    print_matrix((char*)"[Block] after(3)(4)",n,n,a,n);
#endif
	    ///system("pause");
	    //getchar();
	}//for(i=ib;i<end;i++){
	for(i=ib+1;i<=end;i++){
	    for(j=end+1;j<n;j++){
		a[i*n+j]-=a[i*n+ib]*a[ib*n+j];
	    }
	}
#if debugFlag
	print_matrix((char*)"[Block] after(7)(8)",n,n,a,n);
#endif
        /*for(i=end+1;i<n;i++){
	    for(j=end+1;j<n;j++){
		register double t = a[i*n+j];
		for(k=ib;k<=end;k++)
		    t -= a[i*n+k]*a[k*n+j];
	    	a[i*n+j] = t;
	    }
	}*/

	for(i=end+1;i<n;i+=2){
	    for(j=end+1;j<n;j+=2){
		register double c00 = a[i*n+j];		register double c01 = a[i*n+(j+1)];
		register double c10 = a[(i+1)*n+j];	register double c11 = a[(i+1)*n+(j+1)];
		
		register double a00 = a[i*n+ib];	register double a01 = a[i*n+ib+1];
		register double a10 = a[(i+1)*n+ib];	register double a11 = a[(i+1)*n+(ib+1)];
		register double b00 = a[ib*n+j];	register double b01 = a[ib*n+(j+1)];
		register double b10 = a[(ib+1)*n+j];	register double b11 = a[(ib+1)*n+(j+1)];

		c00 -= a00*b00 + a01*b10;	c01 -= a00*b01 + a01*b11;
		c10 -= a10*b00 + a11*b10;	c11 -= a10*b01 + a11*b11;

		a[i*n+j] = c00;	a[i*n+(j+1)] = c01;
		a[(i+1)*n+j] = c10;	a[(i+1)*n+(j+1)] = c11;

	    }

	}
#if debugFlag
	print_matrix((char*)"[Block] after(9)(10)(11)",n,n,a,n);
#endif
    //getchar();
    }//for(ib=0;ib<n;ib++){


}


//AX = b
/*int mydtrsm(double *a, double *b, int n, bool upDown){
    //upDown;
    int i,j;
    double *tmp = (double *)malloc(n*sizeof(double));
    if(upDown){
	tmp[0] = b[0];
	for(i=0;i<n;i++){
	    //tmp[i] = b[i];
	    for(j=0;j<i;j++){
		tmp[i] = b[i];
		tmp[i]-=a[i*n+j]*tmp[j];
	    }
	    tmp[i]/=a[i*n+i];
	}
	
	//b<-tmp
	for(i=0;i<n;i++)
	    b[i]=tmp[i];
    }else{
	for(i=n-1;i>=0;i--){
	    tmp[i] = b[i];
	    for(j=i+1;j<n;j++){
		tmp[i]-=a[i*n+j]*tmp[j];

	    }
	    tmp[i]/=a[i*n+i];
	}
	
	for(i=0;i<n;i++)
	    b[i]=tmp[i];
    
    }
    //free(tmp);

}*/

int mydtrsm2(double *a, double *b, int n){
	double *y = (double *)malloc(n*sizeof(double));
	double *x = (double *)malloc(n*sizeof(double));
	int i,j;
	for(i=0;i<n;i++){
		//y[i] = b[p[i]];
		y[i] = b[i];
		for(j=0;j<i;j++){
			y[i]-=a[i*n+j]*y[j];
		}
		y[i]/=a[i*n+i];
	}

	for(i=n-1;i>=0;i--){
		x[i]=y[i];
		for(j=i+1;j<n;j++){
			x[i]-=a[i*n+j]*x[j];
		}
		x[i]/=a[i*n+i];
	}
	for(i=0;i<n;i++)
            b[i]=x[i];
}

int mydtrsm(double *a, double *b, int n, int upDown){
    int i,j;
    //double *y = (double *)malloc(n*sizeof(double));
    if(upDown == 1){
	for(i=0;i<n;i++){
	    //y[i] = b[i];
	    for(j=0;j<i;j++){
		b[i]-=a[i*n+j]*b[j];
	    }
	}	
    }else{//upDown = 0
	for(i=n-1;i>=0;i--){
            for(j=i+1;j<n;j++){
            	b[i]-=a[i*n+j]*b[j];
            }
            b[i]/=a[i*n+i];
	}
    }
    //for(i=0;i<n;i++)
	//b[i]=y[i];

    //free(y);
    return 0;
}
//extern void print_matrix( char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda );

int Blockdtrsm(double *a, double *b, int n, int upDown){
    int i,j;

    if(upDown == 1){
        for(i=0;i<n;i++){
            register double t = b[i];
	    //register double t;
	    for(j=0;j<i;j++){
                t-=a[i*n+j]*b[j];
            }
	    b[i]=t;
        }
    }else{
        for(i=n-1;i>=0;i--){
            //register double t;
	    register double t = b[i];
	    for(j=i+1;j<n;j++){
                //t-=a[i*n+j]*b[j];
		t-=a[i*n+j]*b[j];
            }
            //t/=a[i*n+i];
	    //b[i] = t;
	    t/=a[i*n+i];
	    b[i]=t;
        }
    }

    return 0;


}

int correctVerify(double *b1, double *b2, int n){
    int error = -1;
    int i;
    for(i=0;i<n;i++){
        error = MAX(abs(b1[i] - b2[i]),error);
    }
    if(error<1e-3)
        printf("\tcorrectness verify error= %d\n",error);
    else
        printf("\toops!!!!!!!! error= %d\n",error);

}

int correctVerifyReverse(double *b1, double *b2, int n){
    int error = -1;
    int i,j;

    int Maxi,Maxj;
    for(i=0;i<n;i++){
	for(j=0;j<n;j++)
	    if(i>=j){
		error = MAX(abs(abs(b1[i*n+j])-abs(b2[j*n+i])),error);
	    	if(error>1e-3){
			//printf("error found:[%d][%d]\n",i,j);
			//printf("b1=%f,b2=%f\n",b1[i*n+j],b2[j*n+i]);
		}
	    }
    }

    if(error<1e-3)
        printf("\tcorrectness verify error= %d\n",error);
    else
        printf("\toops!!!!!!!! error= %d\n",error);

}


int main (int argc, char * argv[])
{
    //double a[16] = {2,5,3,11,5,6,-4,7,-9,-4,2,4,3,2,7,-8};
    //double b[4] = {151,103,16,-32};
    int k, n = Matrix_size;

    if(argc>1){
        n = strtol(argv[1], NULL, 10);
    }else{
	n = Matrix_size;
    }

    double *a = (double *)malloc(n*n*sizeof(double));
    double *b = (double *)malloc(n*sizeof(double));
	for (k = 0;k<n*n;k++)
		a[k] = (abs(rand()*k))%10;
	for(k=0;k<n;k++)
		b[k] = (abs(rand())*k)%100;
    
    //double a[9] = {1,2,1,2,1,1,3,1,1};
    //double b[3] = {1,2,3};
    double *Mya = (double *)malloc(n*n*sizeof(double));
    double *Myb = (double *)malloc(n*sizeof(double));
	for (k = 0;k<n*n;k++)
                Mya[k] = a[k];
        for(k=0;k<n;k++)
                Myb[k] = b[k];
    mySwap(Mya, n);

    double *Blocka = (double *)malloc(n*n*sizeof(double));
    double *Blockb = (double *)malloc(n*sizeof(double));
        for (k = 0;k<n*n;k++)
                Blocka[k] = a[k];
        for(k=0;k<n;k++)
                Blockb[k] = b[k];
    mySwap(Blocka, n);

    // matrix reverse
    //double a[9] = {1,2,3,2,5,7,3,5,3};
    //double b[3] = {140,230,220}; 
    int i,j;
    char    TRANS = 'N';
    int     INFO;
    int     LDA = n;//Matrix_size;
    int     LDB = n;//Matrix_size;
    //int     n = Matrix_size;
    int     NRHS;
    //int     IPIV[Matrix_size];
    //int	    MyIPIV[Matrix_size];
    //int     BlockIPIV[Matrix_size];
    int     IPIV[n];
    int     MyIPIV[n];
    int     BlockIPIV[n];

struct timeval stop, start;

    for(i=0;i<n;i++)
	MyIPIV[i]=i;
    for(i=0;i<n;i++)
        BlockIPIV[i]=i;

//==========================
//==========LACK============
//==========================
#if debugFlag
    print_matrix((char*)"Entry Matrix A", n, n, a, LDA );
#endif
    //LAPACK_dgetrf( lapack_int* m, lapack_int* n, double* a, lapack_int* lda, lapack_int* ipiv, lapack_int *info );
    //info = LAPACKE_dgels(LAPACK_RO`W_MAJOR,'N',m,n,nrhs,*a,lda,*b,ldb);
    //INFO = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,n,n, a ,LDA, IPIV); 
gettimeofday(&start, NULL);
    LAPACK_dgetrf(&n,&n,a,&LDA,IPIV,&INFO);
gettimeofday(&stop, NULL);
float time_dgemmLAPCK_1 = (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f;//(stop.tv_sec - start.tv_sec) * 1000.0f;
    // A should be {3, 5, 5, 2/3, 5/3, 5, 1/3, 1/5, 1}
#if debugFlag
    if(INFO){
        //cout << "an error occured : "<< INFO << endl << endl;
	printf("LAPACK_dgetrf go wrong\n");
    }else{
	print_matrix((char *)"after dgetrf Matrix A", n, n, a, LDA );
        printf("Pivot:");
        for(i=0;i<n;i++)
            printf("[%d]%d",i,IPIV[i]);
        printf("\n");
    }
#endif

//substitution
    for(int i = 0; i < n; i++)
    {
        double tmp = b[IPIV[i]-1];
	b[IPIV[i]-1] = b[i];
	b[i] = tmp;
    }
#if debugFlag
print_matrix((char *)"Matrix B after pivot swap", 1, n, b, LDA );
#endif
    //print_matrix((char *)"Matrix B after swap", 1, n, b, LDA );
    char     SIDE = 'L';
    char     UPLO = 'L';
    char     DIAG = 'U';
    int      M = 1;
    double   aaa    = 1.0;
gettimeofday(&start, NULL);
    dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&n,&M,&aaa,a, &n, b, &n);
    //print_matrix((char *)"Matrix B after first substitution", 1, n, b, LDA );
    UPLO = 'U';
    DIAG = 'N';
    // backward Ux = y
    dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&n,&M,&aaa, a, &n, b, &n);
gettimeofday(&stop, NULL);
float time_dgemmLAPCK_2 = (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f;//(stop.tv_sec - start.tv_sec) * 1000.0f;
#if debugFlag
    print_matrix((char *)"Matrix B after 2nd substitution", 1, n, b, LDA );
#endif

//==========================
//============My============
//==========================
gettimeofday(&start, NULL);
    mydgetrf(&Mya[0],n,&MyIPIV[0]);
gettimeofday(&stop, NULL);
float time_dgemmMY_1 = (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f;//(stop.tv_sec - start.tv_sec) * 1000.0f;
#if debugFlag
    print_matrix((char *)"after mydgetrf Matrix A", n, n, Mya, LDA );
    printf("Pivot:");
        for(i=0;i<n;i++)
            printf("[%d]%d",i,MyIPIV[i]);
        printf("\n");//end of LUU
#endif

//mySwap(Mya, n);
///correctVerify(Mya,a,n*n);
//mySwap(Mya, n);
printf("\tcorrectness Mya & a(reverse)\n");
correctVerifyReverse(Mya,a,n);
//Correctess(&a[0],&Mya[0],n*n);
//CorrectessPIV(&IPIV[0],&MyIPIV[0],n);

#if debugFlag
print_matrix((char *)"myB before pivot swap", 1, n, Myb, LDA );
#endif
    for(i = 0; i < n; i++){
	double tmp = Myb[MyIPIV[i]];
        Myb[MyIPIV[i]] = Myb[i];
        Myb[i] = tmp;
//swap(Myb[i],Myb[MyIPIV[i]]);// swap row by row
    }
#if debugFlag
print_matrix((char *)"myB after pivot swap", 1, n, Myb, LDA );
#endif

    //mydtrsm2(Mya, Myb, n);
    //print_matrix((char *)"My B after 2nd substitution", 1, n, Myb, LDA );
gettimeofday(&start, NULL);
    mydtrsm(Mya, Myb, n,1);
    //print_matrix((char *)"My B after 1st substitution", 1, n, Myb, LDA );
    mydtrsm(Mya, Myb, n,0);
gettimeofday(&stop, NULL);
float time_dgemmMY_2 = (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f;//(stop.tv_sec - start.tv_sec) * 1000.0f;
#if debugFlag
    print_matrix((char *)"My B after 2nd substitution", 1, n, Myb, LDA );
#endif

//correctness (LACKE-My)
    printf("\tcorrectness Myb & b\n");
    correctVerify(Myb,b,n);


//==========================
//==========BLOCk===========
//==========================
//Block algorithm
//correctness (LACKE-Block)

gettimeofday(&start, NULL);
    myBlockdgetrf(Blocka , n, BlockIPIV);
gettimeofday(&stop, NULL);
float time_dgemmB_1 = (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f;//(stop.tv_sec - start.tv_sec) * 1000.0f;
    printf("\tcorrectness Blocka & Mya\n");
    correctVerify(Blocka,Mya,n*n);

    for(int i = 0; i < n; i++)
    {
        double tmp = Blockb[BlockIPIV[i]];
        Blockb[BlockIPIV[i]] = Blockb[i];
        Blockb[i] = tmp;
    }
#if debugFlag
    print_matrix((char *)"BlockB after pivot swap", 1, n, Blockb, LDA );
#endif
gettimeofday(&start, NULL);
    Blockdtrsm(Blocka, Blockb, n,1);
    Blockdtrsm(Blocka, Blockb, n,0);
gettimeofday(&stop, NULL);
float time_dgemmB_2 = (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f;//(stop.tv_sec - start.tv_sec) * 1000.0f;
#if debugFlag
    print_matrix((char *)"BlockB after 2nd substitution", 1, n, Blockb, LDA );
#endif
    printf("\tcorrectness Blockb & b\n");
    correctVerify(Myb,Blockb,n);


//performance
printf("\tN = %d\n",n);
printf("LAPCKE:\n");
printf("\ttook %f\tmilli-sec\n",time_dgemmLAPCK_1);
printf("\ttook %f\tmilli-sec\n",time_dgemmLAPCK_2);
double Gflops;
Gflops = 4*n*((n/6)*n/1000000000.000);
Gflops -= 3*((n/6)*n/1000000000.000);
Gflops -= (n/6)/1000000000.000;
Gflops /= (time_dgemmLAPCK_1/1000);
printf("\tGflops = %f\n", Gflops);


printf("my:\n");
printf("\ttook %f\tmilli-sec\n",time_dgemmMY_1);
printf("\ttook %f\tmilli-sec\n",time_dgemmMY_2);
//double Gflops
Gflops = 0;
Gflops = 4*n*((n/6)*n/1000000000.000);
Gflops -= 3*((n/6)*n/1000000000.000);
Gflops -= (n/6)/1000000000.000;
Gflops /= (time_dgemmMY_1/1000);
printf("\tGflops = %f\n", Gflops);

printf("Blocked:\n");
printf("\ttook %f\tmilli-sec\n",time_dgemmB_1);
printf("\ttook %f\tmilli-sec\n",time_dgemmB_2);
Gflops = 0;
Gflops = 4*n*((n/6)*n/1000000000.000);
Gflops -= 3*((n/6)*n/1000000000.000);
Gflops -= (n/6)/1000000000.000;
Gflops /= (time_dgemmB_1/1000);
printf("\tGflops = %f\n", Gflops);

    //free(Mya);
    //free(Myb);
    //free(a);
    //free(b);
    return 0;
}

