
#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"

#define sizeofMatrix 2048
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define swap(x,y) { x = x + y; y = x - y; x = x - y; }

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

int main(int argc, char* argv[]){

//do something
struct timeval stop, start;
double diff_t;
int i,j,k,n;
int m;
int total; // n*n

n = sizeofMatrix; // initial for 2048

//double *a = (double *)malloc(n * n * sizeof(double));
//double *b = (double *)malloc(n * n * sizeof(double));
//double *c = (double *)malloc(n * n * sizeof(double));
//double *c2 = (double *)malloc(n * n * sizeof(double));
//double *c3 = (double *)malloc(n * n * sizeof(double));

for (m=0;m<6;m++){

//n/=2  2048,1024,512,256,128,64
    

double *a = (double *)malloc(n * n * sizeof(double));
double *b = (double *)malloc(n * n * sizeof(double));

    //n = n/(m+1); //2048, 1024, 512, 256, 128, 64
    total = n*n;
    for(i=0;i<n*n;i++){     //fill the array with random numbers
	a[i]=rand()%100;
	b[i]=rand()%100;
	//c[i]=rand();
    }

gettimeofday(&start, NULL);
double *c = (double *)malloc(n * n * sizeof(double));
for(i=0;i<n*n;i++)
	c[i]=0;

printf("=====dgemm0:=====n=%d\n",n);
    for (i=0; i<n; i++){
	for (j=0; j<n; j++){
    	    for (k=0; k<n; k++){ 
        	c[i*n+j] += a[i*n+k] * b[k*n+j];
            }
    	}
    }

gettimeofday(&stop, NULL);
float time_dgemm0 = (stop.tv_sec - start.tv_sec) * 1000.0f;
printf("\ttook %f\tmilli-sec\n",n, (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f);
//float Gflops = ((2*n*n*n)/(time_dgemm0/1000) )/ 1000000000;
//float Gflops = (n/(time_dgemm0/1000))*(n/(time_dgemm0/1000))*(n/(time_dgemm0/1000))*(2/(time_dgemm0/1000));
// Gflops = (2/(time_dgemm0/1000))*((n/2)/time_dgemm0)*((n/2)/time_dgemm0)*((n/2)/time_dgemm0);
double Gflops;
if(n>100){
    Gflops = 2*(n*n/1000000000.000)/(time_dgemm0/1000);
    Gflops = Gflops*n;
}else{
    Gflops = (2*(n*n*n))/(time_dgemm0/1000)/1000000000.000;
}
printf("\tGflops = %f\n", Gflops);

double *c2 = (double *)malloc(n * n * sizeof(double));
for(i=0;i<total;i++)
        c2[i]=0;
gettimeofday(&start, NULL);
printf("=====dgemm1:=====n=%d\n",n);
    for (i=0; i<n; i++)  
	for (j=0; j<n; j++) {
	register double r =
	c2[i*n+j];
	for (k=0; k<n; k++) 
	    r += a[i*n+k] * b[k*n+j];
	c2[i*n+j] = r;
	}

gettimeofday(&stop, NULL);
float time_dgemm1 = (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f;
printf("\ttook %f\tmilli-sec\n",n, (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f);
//Gflops = ((2*n*n*n)/(time_dgemm1/1000) )/ 1000000000;
//Gflops = ( ( (n/(time_dgemm1/1000))^3 )*(2/(time_dgemm1/1000)) )/1000000000;
//Gflops = (n/(time_dgemm1/1000))*(n/(time_dgemm1/1000))*(n/(time_dgemm1/1000))*(2/(time_dgemm1/1000));
//Gflops = (2/(time_dgemm1/1000))*((n/2)/time_dgemm1)*((n/2)/time_dgemm1)*((n/2)/time_dgemm1);
if(n>100){
    Gflops = 2*(n*n/1000000000.000)/(time_dgemm1/1000);
    Gflops = Gflops*n;
}else{
    Gflops = (2*(n*n*n)/(time_dgemm1/1000))/1000000000.000;
}
printf("\tGflops = %f\n", Gflops);

correctVerify(&c[0],&c2[0],n);

double *c3 = (double *)malloc(n * n * sizeof(double));
for(i=0;i<n*n;i++)
	c3[i]=0;

gettimeofday(&start, NULL);
printf("=====dgemm2:=====n=%d\n",n);
    for (i=0; i<n; i+=2){
	for (j=0; j<n; j+=2) {
        register int t = i*n+j; register int tt = t+n;
	register double c00 = c3[t]; register double c01 = c3[t+1];  register double c10 = c3[tt]; register double c11 = c3[tt+1];
	for (k=0; k<n; k+=2){
	    // 2 by 2 mini matrix multiplication using registers
	    register int ta = i*n+k; register int tta = ta+n; register int tb = k*n+j; register int ttb = tb+n;
	    register double a00 = a[ta]; register double a01 = a[ta+1]; register double a10 = a[tta]; register double a11 = a[tta+1];
     	    register double b00 = b[tb]; register double b01 = b[tb+1]; register double b10 = b[ttb]; register double b11 = b[ttb+1];
      	    c00 += a00*b00 + a01*b10;
       	    c01 += a00*b01 + a01*b11;
       	    c10 += a10*b00 + a11*b10;
      	    c11 += a10*b01 + a11*b11;
	}//k loop
	c3[t] = c00;
        c3[t+1] = c01;
        c3[tt] = c10;
        c3[tt+1] = c11;
	}// j loop
    }// i loop
gettimeofday(&stop, NULL);
float time_dgemm2 = (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f;
printf("\ttook %f\tmilli-sec\n",n, (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f);
//Gflops = ( (16*(n/2)*(n/2)*(n/2) )/(time_dgemm2/1000) )/ 1000000000;
//Gflops = ( ( ((n/2)/(time_dgemm2/1000))^3 )*(16/(time_dgemm2/1000)) )/1000000000;
//Gflops = ((n/2)/(time_dgemm2/1000))*((n/2)/(time_dgemm2/1000))*((n/2)/(time_dgemm2/1000))*(16/(time_dgemm2/1000));
//Gflops = (16/(time_dgemm2/1000))*((n/2)/time_dgemm2)*((n/2)/time_dgemm2)*((n/2)/time_dgemm2);
//Gflops = 16*((n/2)*(n/2)/1000000000)*(n/2)/(time_dgemm2/1000);
if(n>100){
    Gflops = 16*((n/2)*(n/2)/1000000000.000)/(time_dgemm2/1000);
    Gflops = Gflops*(n/2);
}else{
    //Gflops = 16*((n/2)*(n/2)*(n/2)/1000000000.000)/(time_dgemm2/1000);
    Gflops = (16*((n/2)*(n/2)*(n/2))/(time_dgemm2/1000))/1000000000.000;
}
printf("\tGflops = %f\n", Gflops);

correctVerify(&c[0],&c3[0],n);

n = n/2;

free(c3);
free(c2);
free(c);
free(b);
free(a);
}//end of m loop

//free(c3);
//free(c2);
//free(c);
//free(b);
//free(a);

return 0;

}
