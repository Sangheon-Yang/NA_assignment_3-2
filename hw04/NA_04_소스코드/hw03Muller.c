#include "nr.h"
#include "nrutil.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXIT 100
#define N 100
#define NBMAX 20
#define X1 1.0
#define X2 10.0


float func01(float x){
    return (10.0*exp(-x)*sin(2.0*M_PI*x) - 2.0);
}

float dfunc01(float x){
    return (10*-exp(-x)*(sin(2*M_PI*x)-2*M_PI*cos(2*M_PI*x)));
}

static float fx01(x){
    return func01(x);
}

static void function01(float x, float* fn, float* df){
    *fn = func01(x);
    *df = dfunc01(x);
}

float func02(float x){
    return ((x-exp(-x))*(x-exp(-x)));
}

float dfunc02(float x){
    return (2*(x-exp(-x))*(exp(-1)+1));
}

static float fx02(x){
    return func02(x);
}

static void function02(float x, float* fn, float* df){
    *fn = func02(x);
    *df = dfunc02(x);
}

float func03(float x){
    return (cos(x+sqrt(2))+x*(x/2 + sqrt(2)));
}

float dfunc03(float x){
    return (-sin(x+sqrt(2))+x+sqrt(2));
}

static float fx03(x){
    return func03(x);
}

static void function03(float x, float* fn, float* df){
    *fn = func03(x);
    *df = dfunc03(x);
}

float funcMyOwn(float x){
    return (exp(2*x)+x-5);
}

float dfuncMyOwn(float x){
    return (2*exp(2*x)+1);
}

static float fxMyOwn(x){
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

float muller(float (*func)(float), float x0, float x1, float x2, float bound0, float bound1, float xacc){
    void nrerror(char error_text[]);
    int j;
    int sgn;
    float p0, p1, p2, a, b, c, rts, f;
    
    p0 = x0;
    p1 = x1;
    p2 = x2;
    
    for(j=1; j<=MAXIT; j++){
        c = (*func)(p2);
        b = ((p0-p2)*(p0-p2)*((*func)(p1)-(*func)(p2))-(p1-p2)*(p1-p2)*((*func)(p0)-(*func)(p2)))/((p0-p2)*(p1-p2)*(p0-p1));
        a = ((p1-p2)*((*func)(p0)-(*func)(p2))-(p0-p2)*((*func)(p1)-(*func)(p2)))/((p0-p2)*(p1-p2)*(p0-p1));
        
        if (b < 0.0)
            sgn = -1;
        else
            sgn = 1;
        rts = p2 - ((2*c)/(b+sgn*sqrtf(b*b-4*a*c)));
        f = (*func)(rts);
        if (rts < bound0 || rts > bound1)
            nrerror("Jumped out of brackets in muller");
        if(fabs(rts-p2) < xacc || f == 0.0){
            printf("Total Iteration : %d\n",j);
            return rts;
        }
        else{
            p0 = p1;
            p1 = p2;
            p2 = rts;
        }
    }
    nrerror("Maximum number of iterations exceeded in muller");
    return 0.0;
}



int main(void)
{
    int i,nb=NBMAX;
    float xacc,root,*xb1,*xb2;
    float tmp;
    
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(fx,X1,X2,N,xb1,xb2,&nb);
    printf("by Muller Method\n");
    printf("\nRoots of bessj0: on [1, 10]\n");
    printf("%21s %15s\n","x","f(x)\n");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
        tmp = xb2[i]-xb1[i];
        root=muller(fx,xb1[i],xb1[i]+tmp/10.0, xb1[i]+tmp/2.0, xb1[i], xb2[i], xacc);
        printf("root %3d %14.6f %14.6f\n\n",i,root,fx(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    return 0;
}
