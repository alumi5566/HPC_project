
#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"

#define sizeofMatrix 2048//516// 2048
#define sizeofBlock 10
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))


int ijkB_reg_reuse(double *a, double *b, double *c, int n, int B){
//suppose input n = multiple of 36 (best block size from last problem)
int i,j,k;
int i1,j1,k1;
for (i = 0; i < n; i+=B){
    for (j = 0; j < n; j+=B){
        for (k = 0; k < n; k+=B){
		
		for (i1=i; i1<i+B; i1+=2){
	            for (j1=j; j1<j+B; j1+=2) {
        	    register int t = i1*n+j1; register int tt = t+n;
        	    register double c00 = c[t]; register double c01 = c[t+1];  register double c10 = c[tt]; register double c11 = c[tt+1];
        		for (k1=k; k1<k+B; k1+=2){
            		    register int ta = i1*n+k1; register int tta = ta+n; register int tb = k1*n+j1; register int ttb = tb+n;
            		    register double a00 = a[ta]; register double a01 = a[ta+1]; register double a10 = a[tta]; register double a11 = a[tta+1];
            		    register double b00 = b[tb]; register double b01 = b[tb+1]; register double b10 = b[ttb]; register double b11 = b[ttb+1];
            		    c00 += a00*b00 + a01*b10;
            		    c01 += a00*b01 + a01*b11;
            		    c10 += a10*b00 + a11*b10;
            		    c11 += a10*b01 + a11*b11;
        		}
        	    c[t] = c00;
        	    c[t+1] = c01;
        	    c[tt] = c10;
        	    c[tt+1] = c11;
        	    }
    		}// for i

	}//k
    }//j
}//i
}

int kjiB(double *a, double *b, double *c, int n, int B){
int i,j,k;
int i1,j1,k1;
for (k = 0; k < n; k+=B){
    for (j = 0; j < n; j+=B){
        for (i = 0; i < n; i+=B){    
            for (k1 = k; k1 < min(k+B,n); k1++){
     		for (j1 = j; j1 < min(j+B,n); j1++){       
	        		register double r=b[k1*n+j1];
                    for (i1= i; i1 < min(i+B,n); i1++) {
                        c[i1*n+j1] += r * a[i1*n + k1];
                    }
                }
            }
        }
    }
}
}

int jkiB(double *a, double *b, double *c, int n, int B){
int i,j,k;
int i1,j1,k1;
for (j = 0; j < n; j+=B){
    for (k = 0; k < n; k+=B){
        for (i = 0; i < n; i+=B){
	    for (j1 = j; j1 < min(j+B,n); j1++){    
                for (k1 = k; k1 < min(k+B,n); k1++){
                    register double r=b[k1*n+j1];
                    for (i1= i; i1 < min(i+B,n); i1++) {
                        c[i1*n+j1] += r * a[i1*n + k1];
                    }
                }
            }
        }
    }
}
}

int ikjB(double *a, double *b, double *c, int n, int B){
int i,j,k;
int i1,j1,k1;
for (i = 0; i < n; i+=B){
    for (k = 0; k < n; k+=B){
        for (j = 0; j < n; j+=B){
            for (i1= i; i1 < min(i+B,n); i1++) {
                for (k1 = k; k1 < min(k+B,n); k1++){
                    register double r=a[i1*n+k1];
                    for (j1 = j; j1 < min(j+B,n); j1++){
                        c[i1*n+j1] += r * b[k1*n + j1];
                    }
                }
            }
        }
    }
}
}


int kijB(double *a, double *b, double *c, int n, int B){
int i,j,k;
int i1,j1,k1;
for (k = 0; k < n; k+=B){
    for (i = 0; i < n; i+=B){
        for (j = 0; j < n; j+=B){
            for (k1 = k; k1 < min(k+B,n); k1++){
                for (i1= i; i1 < min(i+B,n); i1++) {
                    register double r=a[i1*n+k1];
                    for (j1 = j; j1 < min(j+B,n); j1++){
                        c[i1*n+j1] += r * b[k1*n + j1];
                    }
                }
            }
        }
    }
}
}

int jikB(double *a, double *b, double *c, int n, int B){
int i,j,k;
int i1,j1,k1;
for (j = 0; j < n; j+=B){
    for (i = 0; i < n; i+=B){
        for (k = 0; k < n; k+=B){
            for (j1 = j; j1 < min(j+B,n); j1++){
                for (i1= i; i1 < min(i+B,n); i1++) {
                    register double r=c[i1*n+j1];
                    for (k1 = k; k1 < min(k+B,n); k1++){
                        r += a[i1*n + k1]*b[k1*n + j1];
                    }
                    c[i1*n+j1]=r;
                }
            }
        }
    }
}
}

int ijkB(double *a, double *b, double *c, int n, int B){
int i,j,k;
int i1,j1,k1;
for (i = 0; i < n; i+=B){
    for (j = 0; j < n; j+=B){
        for (k = 0; k < n; k+=B){
            for (i1 = i; i1 < min(i+B,n); i1++){
                for (j1= j; j1 < min(j+B,n); j1++) {
                    register double r=c[i1*n+j1];
                    for (k1 = k; k1 < min(k+B,n); k1++){
                        r += a[i1*n + k1]*b[k1*n + j1];
                    }
                    c[i1*n+j1]=r;
                }
            }
	}//k
    }//j
}//i
}

int kji(double *a, double *b, double *c, int n){
    int i,j,k;
    for (k=0; k<n; k++)
        for (j=0; j<n; j++) {
        register double r = b[k*n+j];
        for (i=0; i<n; i++)
            c[i*n+j] += r * a[i*n+k];
        }
    return 0;
}

int jki(double *a, double *b, double *c, int n){
int i,j,k;
    for (j=0; j<n; j++)
        for (k=0; k<n; k++) {
        register double r = b[k*n+j];
        for (i=0; i<n; i++) 
            c[i*n+j] += r * a[i*n+k];
        }
    return 0;
}


int ikj(double *a, double *b, double *c, int n){
int i,j,k;
    for (i=0; i<n; i++)
        for (k=0; k<n; k++) {
        register double r = a[i*n+k];
        for (j=0; j<n; j++) 
            c[i*n+j] += r * b[k*n+j];
        }
    return 0;
}

int kij(double *a, double *b, double *c, int n){
int i,j,k;
    for (k=0; k<n; k++)
        for (i=0; i<n; i++) {
        register double r = a[i*n+k];
        for (j=0; j<n; j++)
            c[i*n+j] += r * b[k*n+j];
        }
    return 0;
}

int jik(double *a, double *b, double *c, int n){
int i,j,k;
   for (j=0; j<n; j++)
        for (i=0; i<n; i++) {
        register double r = c[i*n+j];
        for (k=0; k<n; k++)
            r+= a[i*n+k] * b[k*n+j];
        c[i*n+j]= r;
        }
    return 0;
}

int ijk(double *a, double *b, double *c, int n){
int i,j,k; 
   for (i=0; i<n; i++)
        for (j=0; j<n; j++) {
        register double r = c[i*n+j];
        for (k=0; k<n; k++)
            r+= a[i*n+k] * b[k*n+j];
        c[i*n+j]= r;
        }
    return 0;
}

int correctVerify(double *a, double *b, int n){
//Verify correctness
    int error = -1;
    int total = n*n;
    int i; 
    for(i=0;i<total;i++){
	error = MAX(abs(a[i] - b[i]),error);
    }
    //printf("\terror = %d\n",error);
    //printf("\tdiffTime =:%f\n",time_dgemm0-time_dgemm1);
    if(error<1e-3)
	printf("\tcorrectness verify error= %d\n",error);
    else
	printf("\toops!!!!!!!! error= $d\n",error);
}

int time_and_Gflops(float timeDiff,int n){
    //float time_dgemm0 = (stop.tv_sec - start.tv_sec) * 1000.0f;
    //printf("\ttook %f\tmilli-sec\n",n, (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f);
    //float time_dgemm = (stop-start)*1000.0f;
    printf("\ttook %f\tmilli-sec\n", timeDiff);
    double Gflops;
    if(n>100){
	Gflops = 2*(n*n/1000000000.000)/(timeDiff/1000);
	Gflops = Gflops*n;
    }else{
	Gflops = (2*(n*n*n))/(timeDiff/1000)/1000000000.000;
    }
    printf("\tGflops = %f\n", Gflops);
}
                                                                
int main(int argc, char* argv[]){

//do something
struct timeval stop, start;
double diff_t;
int i,j,k,n;
int m;
int total; // n*n
int B = sizeofBlock;

n = sizeofMatrix; // initial for 2048

//if(n%36!=0)
//	n = n-(n%36);
//printf("new n = %d",n);

double *a = (double *)malloc(n * n * sizeof(double));
double *b = (double *)malloc(n * n * sizeof(double));
double *c = (double *)malloc(n * n * sizeof(double));
double *c2 = (double *)malloc(n * n * sizeof(double));
//double *c3 = (double *)malloc(n * n * sizeof(double));

//for (m=0;m<1;m++){

    //n = n/(m+1); //2048, 1024, 512, 256, 128, 64
    total = n*n;
    for(i=0;i<total;i++){     //fill the array with random numbers
	a[i]=rand()%100;
	b[i]=rand()%100;
	c[i]=0;
	c2[i]=0;
    }

printf("=====Block_1:=====n=%d\n",n);
gettimeofday(&start, NULL);
    for (i=0; i<n; i++)
        for (j=0; j<n; j++) {
        register double r = c[i*n+j];
        for (k=0; k<n; k++)
            r+= a[i*n+k] * b[k*n+j];
        c[i*n+j]= r;
        }

gettimeofday(&stop, NULL);
float time_dgemm0 = (stop.tv_sec - start.tv_sec) * 1000.0f;
printf("\ttook %f\tmilli-sec\n",n, (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f);
double Gflops;
if(n>100){
    Gflops = 2*(n*n/1000000000.000)/(time_dgemm0/1000);
    Gflops = Gflops*n;
}else{
    Gflops = (2*(n*n*n))/(time_dgemm0/1000)/1000000000.000;
}
printf("\tGflops = %f\n", Gflops);


printf("=====ijk:=====n=%d\n",n);
double *c_6algo = (double *)malloc(n * n * sizeof(double));
for(i=0;i<total;i++)
    c_6algo[i]=0;
gettimeofday(&start, NULL);
ijk(&a[0],&b[0],&c_6algo[0],n);
gettimeofday(&stop, NULL);
correctVerify(&c[0],&c_6algo[0],n);
float timeDiff = (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f;
time_and_Gflops(timeDiff, n);
free(c_6algo);

B = 64;
printf("=====Block_2:=====n=%d B=%d\n",n,B);
printf("=====ijk:=====n=%d\n",n);
gettimeofday(&start, NULL);
//ijkB(&a[0],&b[0],&c2[0],n,B);
ijkB_reg_reuse(&a[0],&b[0],&c2[0],n,B);
gettimeofday(&stop, NULL);
timeDiff = (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f;
time_and_Gflops(timeDiff, n);
correctVerify(&c[0],&c2[0],n);
free(c2);

//free(c2);
free(c);
free(b);
free(a);

return 0;

}
