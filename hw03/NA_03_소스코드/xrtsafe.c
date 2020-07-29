
/* Driver for routine rtsafe */
#include <math.h>
#include <stdio.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define N 100
#define NBMAX 20
#define X1 1.0
#define X2 10.0

//fn 01
float func01(float x){
    float ans = (10.0*expf(-1*x)*sin(2.0*M_PI*x) - 2.0);
    return ans;
}

float dfunc01(float x){
    float ans = (10*-expf(-1*x)*(sin(2*M_PI*x)-2*M_PI*cos(2*M_PI*x)));
    return ans;
}

static float fx01(float x){
    return func01(x);
}

static void function01(float x, float* fn, float* df){
    *fn = func01(x);
    *df = dfunc01(x);
}

float func02(float x){
    float ans = ((x - expf(-1*x))*(x - expf(-1*x)));
    //float ans = x*x - 2*x*exp(-1*x) + exp(-2*x);
    return ans;
}

float dfunc02(float x){
    float ans = (2*(x - expf(-1*x))*(expf(-1*x) + 1));
    return ans;
}

static float fx02(float x){
    return func02(x);
}

static void function02(float x, float* fn, float* df){
    *fn = func02(x);
    *df = dfunc02(x);
}

float func03(float x){
    float ans = (cos(x+sqrt(2))+x*(x/2 + sqrt(2)));
    return ans;
}

float dfunc03(float x){
    float ans = (-sin(x+sqrt(2))+x+sqrt(2));
    return ans;
}

static float fx03(float x){
    return func03(x);
}

static void function03(float x, float* fn, float* df){
    *fn = func03(x);
    *df = dfunc03(x);
}

float funcMyOwn(float x){
    float ans = (expf(2*x)+x-5);
    return ans;
}

float dfuncMyOwn(float x){
    float ans = (2*expf(2*x)+1);
    return ans;
}

static float fxMyOwn(float x){
    return funcMyOwn(x);
}

static void functionMyOwn(float x, float* fn, float* df){
    *fn = funcMyOwn(x);
    *df = dfuncMyOwn(x);
}

static float fx(float x)
{
	return bessj0(x);
}

static void funcd(float x,float *fn, float *df)
{
	*fn=bessj0(x);
	*df = -bessj1(x);
}

int main(void)
{
	int i,nb=NBMAX;
	float xacc,root,*xb1,*xb2;

    float b0,b1;
    b0 = 0.0;
    b1 = 2.0;
    
    
	xb1=vector(1,NBMAX);
	xb2=vector(1,NBMAX);
	zbrak(fx,X1,X2,N,xb1,xb2,&nb);
    //zbrak(fx,b0,b1,N,xb1,xb2,&nb);
    
    printf("by Newton w/ Bracketing Method\n");
	printf("\nRoots of bessj0 on [1, 10]:\n\n");
    //printf("\nRoots of 10e^(-x)sin(2PIx)-2 on [0.1, 1]:\n\n");
    //printf("\nRoots of x^2-2xe^(-x)+e^(-2x) on [0, 1]:\n\n");
    //printf("\nRoots of cos(x+sqrt(2))+x(x/2+sqrt(2)) on [-2, -1]:\n\n");
    //printf("\nRoots of e^(2x)+x-5 on [0, 2]:\n");
	printf("%21s %15s\n","x","f(x)\n");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtsafe(funcd,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %14.6f\n\n",i,root,fx(root));
	}
	free_vector(xb2,1,NBMAX);
	free_vector(xb1,1,NBMAX);
	return 0;
}
#undef NRANSI
