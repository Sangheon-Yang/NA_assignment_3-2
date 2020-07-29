//
//  hw04.c
//  NA_hw04
//
//  Created by 양상헌 on 24/09/2019.
//  Copyright © 2019 양상헌. All rights reserved.
//

#include "nr.h"
#include "nrutil.h"
#include <math.h>
#include <stdio.h>

#define N 1
#define NBMAX 20
#define X1 0.0
#define X2 400.0
#define MAXIT 1000

float func01(float x){
    float ans = (expf(-0.005*x)*cosf(0.05*sqrtf(2000-0.01*x*x))-0.01);
    return ans;
}
float dfunc01(float x){
    float ans = (expf(-0.005*x)*((0.0005*x*sinf(0.05*sqrtf(2000-0.01*x*x))/sqrtf(2000-0.01*x*x) )-0.005*cosf(0.05*sqrtf(2000-0.01*x*x))));
    return ans;
}
static float fx01(float x){
    return func01(x);
}
static void funcSet01(float x, float* fn, float* df){
    *fn = func01(x);
    *df = dfunc01(x);
}

static float funcEx32(float x){
    float ans = (x/powf((x*x+0.81),3.0/2.0) - 0.0885*M_PI);
    return ans;
}

static float funcEx36(float x){
    float ans = (0.99403 - 1.2 + 1.671*x/10000.0 + 9.7215*x*x/100000000.0 - 9.5838*x*x*x/100000000000.0 + 1.9520*x*x*x*x/100000000000000.0);
    return ans;
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


int main(){
    
    int i,nb=NBMAX;
    float xacc,root,*xb1,*xb2;
    float tmp;
    printf("---------------------------------------------------------------------------------------------------------\n");
    //Bisection Method 0.0001% r.e.
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(funcEx32,0,100,100,xb1,xb2,&nb);
    printf("--- < by Bisection Method > ---\n");
    printf("\nRoots of Ex8.32 on [0, 100] with r.e. 0.0001%% \n");
    printf("%21s %15s\n","x","f(x)");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
        root=rtbis(funcEx32,xb1[i],xb2[i],xacc);
        printf("root %3d %14.9f %14.9f\n\n",i,root,funcEx32(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    
    printf("---------------------------------------------------------------------------------------------------------\n");
    //Bisection Method 0.0001% r.e.
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(funcEx36,0,2000,100,xb1,xb2,&nb);
    printf("--- < by Bisection Method > ---\n");
    printf("\nRoots of Ex8.36 on [0,2000] with r.e. 0.0001%% \n");
    printf("%21s %15s\n","x","f(x)");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
        root=rtbis(funcEx36,xb1[i],xb2[i],xacc);
        printf("root %3d %14.9f %14.9f\n\n",i,root,funcEx36(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    
    printf("---------------------------------------------------------------------------------------------------------\n\n");
    //Bisection Method 0.0001% r.e.
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(fx01,X1,X2,N,xb1,xb2,&nb);
    printf("--- < by Bisection Method > ---\n");
    printf("\nRoots of e^(-0.005*x)*cos(0.05*sqrt(2000-0.01*x^2))-0.01 on [0, 400] with r.e. 0.0001%% \n");
    printf("%21s %15s\n","x","f(x)");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
        root=rtbis(fx01,xb1[i],xb2[i],xacc);
        printf("root %3d %14.9f %14.9f\n\n",i,root,fx01(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    
    //Bisection Method 0.000001% r.e.
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(fx01,X1,X2,N,xb1,xb2,&nb);
    //printf("by Bisection Method\n");
    printf("\nRoots of e^(-0.005*x)*cos(0.05*sqrt(2000-0.01*x^2))-0.01 on [0, 400] with r.e. 0.000001%% \n");
    printf("%21s %15s\n","x","f(x)");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-8)*(xb1[i]+xb2[i])/2.0;
        root=rtbis(fx01,xb1[i],xb2[i],xacc);
        printf("root %3d %14.9f %14.9f\n\n",i,root,fx01(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    
    printf("---------------------------------------------------------------------------------------------------------\n\n");
    
    //Linear Interpolation Method 0.0001% r.e.
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(fx01,X1,X2,N,xb1,xb2,&nb);
    printf("--- < by Linear Interpolation Method(False Position) > ---\n");
    printf("\nRoots of e^(-0.005*x)*cos(0.05*sqrt(2000-0.01*x^2))-0.01 on [0, 400] with r.e. 0.0001%% \n");
    printf("%21s %15s\n","x","f(x)");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
        root=rtflsp(fx01,xb1[i],xb2[i],xacc);
        printf("root %3d %14.9f %14.9f\n\n",i,root,fx01(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    
    //Linear Interpolation Method 0.000001% r.e.
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(fx01,X1,X2,N,xb1,xb2,&nb);
   // printf("by Linear Interpolation Method(False Position)\n");
    printf("\nRoots of e^(-0.005*x)*cos(0.05*sqrt(2000-0.01*x^2))-0.01 on [0, 400] with r.e. 0.000001%% \n");
    printf("%21s %15s\n","x","f(x)");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-8)*(xb1[i]+xb2[i])/2.0;
        root=rtflsp(fx01,xb1[i],xb2[i],xacc);
        printf("root %3d %14.9f %14.9f\n\n",i,root,fx01(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    
    printf("---------------------------------------------------------------------------------------------------------\n\n");
    
    //Newton-Raphson Method 0.0001% r.e.
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(fx01,X1,X2,N,xb1,xb2,&nb);
    printf("--- < by Newton-Raphson Method > ---\n");
    printf("\nRoots of e^(-0.005*x)*cos(0.05*sqrt(2000-0.01*x^2))-0.01 on [0, 400] with r.e. 0.0001%% \n");
    printf("%21s %15s\n","x","f(x)");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
        root=rtnewt(funcSet01,xb1[i],xb2[i],xacc);
        printf("root %3d %14.9f %14.9f\n\n",i,root,fx01(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    
    //Newton-Raphson Method 0.000001% r.e.
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(fx01,X1,X2,N,xb1,xb2,&nb);
    //printf("by Newton-Raphson Method\n");
    printf("\nRoots of e^(-0.005*x)*cos(0.05*sqrt(2000-0.01*x^2))-0.01 on [0, 400] with r.e. 0.000001%% \n");
    printf("%21s %15s\n","x","f(x)");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-8)*(xb1[i]+xb2[i])/2.0;
        root=rtnewt(funcSet01,xb1[i],xb2[i],xacc);
        printf("root %3d %14.9f %14.9f\n\n",i,root,fx01(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    
    printf("---------------------------------------------------------------------------------------------------------\n\n");
    
    //Secant Method 0.0001%
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(fx01,X1,X2,N,xb1,xb2,&nb);
    printf("--- < by Secant Method > ---\n");
    printf("\nRoots of e^(-0.005*x)*cos(0.05*sqrt(2000-0.01*x^2))-0.01 on [0, 400] with r.e. 0.0001%% \n");
    printf("%21s %15s\n","x","f(x)");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
        root=rtsec(fx01,xb1[i],xb2[i],xacc);
        printf("root %3d %14.9f %14.9f\n\n",i,root,fx01(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    
    //Secant Method 0.000001%
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(fx01,X1,X2,N,xb1,xb2,&nb);
    //printf("by Secant Method\n");
    printf("\nRoots of e^(-0.005*x)*cos(0.05*sqrt(2000-0.01*x^2))-0.01 on [0, 400] with r.e. 0.000001%% \n");
    printf("%21s %15s\n","x","f(x)");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-8)*(xb1[i]+xb2[i])/2.0;
        root=rtsec(fx01,xb1[i],xb2[i],xacc);
        printf("root %3d %14.9f %14.9f\n\n",i,root,fx01(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    
    printf("---------------------------------------------------------------------------------------------------------\n\n");
    
    //Newton-Raphson with Bracketin Method 0.0001%
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(fx01,X1,X2,N,xb1,xb2,&nb);
    printf("--- < by Newton w/ Bracketing Method > ---\n");
    printf("\nRoots of e^(-0.005*x)*cos(0.05*sqrt(2000-0.01*x^2))-0.01 on [0, 400] with r.e. 0.0001%% \n");
    printf("%21s %15s\n","x","f(x)");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
        root=rtsafe(funcSet01,xb1[i],xb2[i],xacc);
        printf("root %3d %14.9f %14.9f\n\n",i,root,fx01(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    
    //Newton-Raphson with Bracketin Method 0.000001%
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(fx01,X1,X2,N,xb1,xb2,&nb);
    //printf("by Newton w/ Bracketing Method\n");
    printf("\nRoots of e^(-0.005*x)*cos(0.05*sqrt(2000-0.01*x^2))-0.01 on [0, 400] with r.e. 0.000001%% \n");
    printf("%21s %15s\n","x","f(x)");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-8)*(xb1[i]+xb2[i])/2.0;
        root=rtsafe(funcSet01,xb1[i],xb2[i],xacc);
        printf("root %3d %14.9f %14.9f\n\n",i,root,fx01(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    
    printf("---------------------------------------------------------------------------------------------------------\n\n");
    
    //Muller Method 0.0001%
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(fx01,X1,X2,N,xb1,xb2,&nb);
    printf("--- < by Muller Method > ---\n");
    printf("\nRoots of e^(-0.005*x)*cos(0.05*sqrt(2000-0.01*x^2))-0.01 on [0, 400] with r.e. 0.0001%% \n");
    printf("%21s %15s\n","x","f(x)");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
        tmp = xb2[i]-xb1[i];
        root=muller(fx01,100.0,200.0, 300.0, xb1[i], xb2[i], xacc);
        //root=muller(fx01,xb1[i],xb1[i]+tmp/10.0, xb1[i]+tmp/5.0, xb1[i], xb2[i], xacc);
        printf("root %3d %14.9f %14.9f\n\n",i,root,fx01(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    
    //Muller Method 0.000001%
    xb1=vector(1,NBMAX);
    xb2=vector(1,NBMAX);
    zbrak(fx01,X1,X2,N,xb1,xb2,&nb);
    //printf("by Muller Method\n");
    printf("\nRoots of e^(-0.005*x)*cos(0.05*sqrt(2000-0.01*x^2))-0.01 on [0, 400] with r.e. 0.000001%% \n");
    printf("%21s %15s\n","x","f(x)");
    for (i=1;i<=nb;i++) {
        xacc=(1.0e-8)*(xb1[i]+xb2[i])/2.0;
        tmp = xb2[i]-xb1[i];
        root=muller(fx01,100.0,200.0, 300.0, xb1[i], xb2[i], xacc);
        //root=muller(fx01,xb1[i],xb1[i]+tmp/10.0, xb1[i]+tmp/5.0, xb1[i], xb2[i], xacc);
        printf("root %3d %14.9f %14.9f\n\n",i,root,fx01(root));
    }
    free_vector(xb2,1,NBMAX);
    free_vector(xb1,1,NBMAX);
    
    printf("---------------------------------------------------------------------------------------------------------\n\n");
    
    return 0;
}
