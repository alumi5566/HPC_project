

/* Calling DGELS using row-major order */
#include <stdlib.h>
#include <stdio.h>
//#include </home/cyu059/HW2/include/lapacke.h>
#include <iostream>
#include "lapacke.h"
#include "blas.h"

#define Matrix_size 3000
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
	    for(j=i+1;j<n;j++){
            	a[j*n+i] = a[j*n+i]/a[i*n+i];
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
        for(i=end+1;i<n;i++){
	    for(j=end+1;j<n;j++){
		for(k=ib;k<=end;k++)
		    a[i*n+j] -= a[i*n+k]*a[k*n+j];
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
    int     LDA = Matrix_size;
    int     LDB = Matrix_size;
    //int     n = Matrix_size;
    int     NRHS;
    int     IPIV[Matrix_size];
    int	    MyIPIV[Matrix_size];
    int     BlockIPIV[Matrix_size];

    for(i=0;i<n;i++)
	MyIPIV[i]=i;
    for(i=0;i<n;i++)
        BlockIPIV[i]=i;
#if debugFlag
    print_matrix((char*)"Entry Matrix A", n, n, a, LDA );
#endif
    //LAPACK_dgetrf( lapack_int* m, lapack_int* n, double* a, lapack_int* lda, lapack_int* ipiv, lapack_int *info );
    //info = LAPACKE_dgels(LAPACK_RO`W_MAJOR,'N',m,n,nrhs,*a,lda,*b,ldb);
    //INFO = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,n,n, a ,LDA, IPIV); 
    LAPACK_dgetrf(&n,&n,a,&LDA,IPIV,&INFO);
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

    mydgetrf(&Mya[0],n,&MyIPIV[0]);
#if debugFlag
    print_matrix((char *)"after mydgetrf Matrix A", n, n, Mya, LDA );
    printf("Pivot:");
        for(i=0;i<n;i++)
            printf("[%d]%d",i,MyIPIV[i]);
        printf("\n");//end of LUU
#endif

//mySwap(Mya, n);
//correctVerify(Mya,a,n*n);
//mySwap(Mya, n);
correctVerifyReverse(Mya,a,n);
//Correctess(&a[0],&Mya[0],n*n);
//CorrectessPIV(&IPIV[0],&MyIPIV[0],n);

#if debugFlag
print_matrix((char *)"myB before pivot swap", 1, n, Myb, LDA );
#endif
    for(i = 0; i < n; i++)
    {
        double tmp = Myb[MyIPIV[i]];
        Myb[MyIPIV[i]] = Myb[i];
        Myb[i] = tmp;
        //swap(Myb[i],Myb[MyIPIV[i]]);// swap row by row
    }
#if debugFlag
print_matrix((char *)"myB after pivot swap", 1, n, Myb, LDA );
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
    dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&n,&M,&aaa,a, &n, b, &n);
    //print_matrix((char *)"Matrix B after first substitution", 1, n, b, LDA );
    UPLO = 'U';
    DIAG = 'N';
    // backward Ux = y
    dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&n,&M,&aaa, a, &n, b, &n);
    print_matrix((char *)"Matrix B after 2nd substitution", 1, n, b, LDA );


//My substitution
    //mydtrsm2(Mya, Myb, n);
    //print_matrix((char *)"My B after 2nd substitution", 1, n, Myb, LDA );
    mydtrsm(Mya, Myb, n,1);
    //print_matrix((char *)"My B after 1st substitution", 1, n, Myb, LDA );
    mydtrsm(Mya, Myb, n,0);
    print_matrix((char *)"My B after 2nd substitution", 1, n, Myb, LDA );
    


//
//correctness
    printf("Myb & b\n");
    correctVerify(Myb,b,n);

    myBlockdgetrf(Blocka , n, BlockIPIV);
    printf("Blocka & Mya\n");
    correctVerify(Blocka,Mya,n*n);

    //free(Mya);
    //free(Myb);
    //free(a);
    //free(b);
    return 0;
}

