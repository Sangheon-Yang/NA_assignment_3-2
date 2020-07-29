//
//  main.c
//  NA_NR_practice
//
//  Created by 양상헌 on 09/09/2019.
//  Copyright © 2019 양상헌. All rights reserved.
//  2015004693

#include <stdio.h>
#include <math.h>

#define CONVF(i) ((float)(i))
void macharFloat(int *ibeta, int *it, int *irnd, int *ngrd, int *machep, int *negep,
                 int *iexp, int *minexp, int *maxexp, float *eps, float *epsneg,
                 float *xmin, float *xmax)
{
    int i,itemp,iz,j,k,mx,nxres;
    float a,b,beta,betah,betain,one,t,temp,temp1,tempa,two,y,z,zero;
    
    one=CONVF(1);
    two=one+one;
    zero=one-one;
    a=one;
    do {
        a += a;
        temp=a+one;
        temp1=temp-a;
    } while (temp1-one == zero);
    b=one;
    do {
        b += b;
        temp=a+b;
        itemp=(int)(temp-a);
    } while (itemp == 0);
    *ibeta=itemp;
    beta=CONVF(*ibeta);
    *it=0;
    b=one;
    do {
        ++(*it);
        b *= beta;
        temp=b+one;
        temp1=temp-b;
    } while (temp1-one == zero);
    *irnd=0;
    betah=beta/two;
    temp=a+betah;
    if (temp-a != zero) *irnd=1;
    tempa=a+beta;
    temp=tempa+betah;
    if (*irnd == 0 && temp-tempa != zero) *irnd=2;
    *negep=(*it)+3;
    betain=one/beta;
    a=one;
    for (i=1;i<=(*negep);i++) a *= betain;
    b=a;
    for (;;) {
        temp=one-a;
        if (temp-one != zero) break;
        a *= beta;
        --(*negep);
    }
    *negep = -(*negep);
    *epsneg=a;
    *machep = -(*it)-3;
    a=b;
    for (;;) {
        temp=one+a;
        if (temp-one != zero) break;
        a *= beta;
        ++(*machep);
    }
    *eps=a;
    *ngrd=0;
    temp=one+(*eps);
    if (*irnd == 0 && temp*one-one != zero) *ngrd=1;
    i=0;
    k=1;
    z=betain;
    t=one+(*eps);
    nxres=0;
    for (;;) {
        y=z;
        z=y*y;
        a=z*one;
        temp=z*t;
        if (a+a == zero || fabs(z) >= y) break;
        temp1=temp*betain;
        if (temp1*beta == z) break;
        ++i;
        k += k;
    }
    if (*ibeta != 10) {
        *iexp=i+1;
        mx=k+k;
    } else {
        *iexp=2;
        iz=(*ibeta);
        while (k >= iz) {
            iz *= *ibeta;
            ++(*iexp);
        }
        mx=iz+iz-1;
    }
    for (;;) {
        *xmin=y;
        y *= betain;
        a=y*one;
        temp=y*t;
        if (a+a != zero && fabs(y) < *xmin) {
            ++k;
            temp1=temp*betain;
            if (temp1*beta == y && temp != y) {
                nxres=3;
                *xmin=y;
                break;
            }
        }
        else break;
    }
    *minexp = -k;
    if (mx <= k+k-3 && *ibeta != 10) {
        mx += mx;
        ++(*iexp);
    }
    *maxexp=mx+(*minexp);
    *irnd += nxres;
    if (*irnd >= 2) *maxexp -= 2;
    i=(*maxexp)+(*minexp);
    if (*ibeta == 2 && !i) --(*maxexp);
    if (i > 20) --(*maxexp);
    if (a != y) *maxexp -= 2;
    *xmax=one-(*epsneg);
    if ((*xmax)*one != *xmax) *xmax=one-beta*(*epsneg);
    *xmax /= (*xmin*beta*beta*beta);
    i=(*maxexp)+(*minexp)+3;
    for (j=1;j<=i;j++) {
        if (*ibeta == 2) *xmax += *xmax;
        else *xmax *= beta;
    }
}
#undef CONVF


#define CONVD(i) ((double)(i))
void macharDouble(int *ibeta, int *it, int *irnd, int *ngrd, int *machep, int *negep,
                  int *iexp, int *minexp, int *maxexp, double *eps, double *epsneg,
                  double *xmin, double *xmax)
{
    int i,itemp,iz,j,k,mx,nxres;
    double a,b,beta,betah,betain,one,t,temp,temp1,tempa,two,y,z,zero;
    
    one=CONVD(1);
    two=one+one;
    zero=one-one;
    a=one;
    do {
        a += a;
        temp=a+one;
        temp1=temp-a;
    } while (temp1-one == zero);
    b=one;
    do {
        b += b;
        temp=a+b;
        itemp=(int)(temp-a);
    } while (itemp == 0);
    *ibeta=itemp;
    beta=CONVD(*ibeta);
    *it=0;
    b=one;
    do {
        ++(*it);
        b *= beta;
        temp=b+one;
        temp1=temp-b;
    } while (temp1-one == zero);
    *irnd=0;
    betah=beta/two;
    temp=a+betah;
    if (temp-a != zero) *irnd=1;
    tempa=a+beta;
    temp=tempa+betah;
    if (*irnd == 0 && temp-tempa != zero) *irnd=2;
    *negep=(*it)+3;
    betain=one/beta;
    a=one;
    for (i=1;i<=(*negep);i++) a *= betain;
    b=a;
    for (;;) {
        temp=one-a;
        if (temp-one != zero) break;
        a *= beta;
        --(*negep);
    }
    *negep = -(*negep);
    *epsneg=a;
    *machep = -(*it)-3;
    a=b;
    for (;;) {
        temp=one+a;
        if (temp-one != zero) break;
        a *= beta;
        ++(*machep);
    }
    *eps=a;
    *ngrd=0;
    temp=one+(*eps);
    if (*irnd == 0 && temp*one-one != zero) *ngrd=1;
    i=0;
    k=1;
    z=betain;
    t=one+(*eps);
    nxres=0;
    for (;;) {
        y=z;
        z=y*y;
        a=z*one;
        temp=z*t;
        if (a+a == zero || fabs(z) >= y) break;
        temp1=temp*betain;
        if (temp1*beta == z) break;
        ++i;
        k += k;
    }
    if (*ibeta != 10) {
        *iexp=i+1;
        mx=k+k;
    } else {
        *iexp=2;
        iz=(*ibeta);
        while (k >= iz) {
            iz *= *ibeta;
            ++(*iexp);
        }
        mx=iz+iz-1;
    }
    for (;;) {
        *xmin=y;
        y *= betain;
        a=y*one;
        temp=y*t;
        if (a+a != zero && fabs(y) < *xmin) {
            ++k;
            temp1=temp*betain;
            if (temp1*beta == y && temp != y) {
                nxres=3;
                *xmin=y;
                break;
            }
        }
        else break;
    }
    *minexp = -k;
    if (mx <= k+k-3 && *ibeta != 10) {
        mx += mx;
        ++(*iexp);
    }
    *maxexp=mx+(*minexp);
    *irnd += nxres;
    if (*irnd >= 2) *maxexp -= 2;
    i=(*maxexp)+(*minexp);
    if (*ibeta == 2 && !i) --(*maxexp);
    if (i > 20) --(*maxexp);
    if (a != y) *maxexp -= 2;
    *xmax=one-(*epsneg);
    if ((*xmax)*one != *xmax) *xmax=one-beta*(*epsneg);
    *xmax /= (*xmin*beta*beta*beta);
    i=(*maxexp)+(*minexp)+3;
    for (j=1;j<=i;j++) {
        if (*ibeta == 2) *xmax += *xmax;
        else *xmax *= beta;
    }
}
#undef CONVD

// self_made get_eps()
/**/
int get_eps_float(void){
    int result = 0;
    float one = 1.0;
    float temp = 1.0;
    //printf("%f",i);
    //float one = 1.0;
    do{
        result = result+1;
        temp = temp/2;
        //printf("%d\n",result);
    }while(one + temp != one);
    
    return result;
}

int get_eps_double(void){
    int result = 0;
    double one = 1.0;
    double temp = 1.0;
    //printf("%f",i);
    //float one = 1.0;
    do{
        result = result+1;
        temp = temp/2;
        //printf("%d\n",result);
    }while(one + temp != one);
    
    return result;
}

// in machar()
/*
 *it=0;
 b=one;
 do {
 ++(*it);
 b *= beta;
 temp=b+one;
 temp1=temp-b;
 } while (temp1-one == zero);
 */

int main(void)
{
    int ibeta,iexp,irnd,it,machep,maxexp,minexp,negep,ngrd;
    float eps,epsneg,xmax,xmin;
    int ibetaD,iexpD,irndD,itD,machepD,maxexpD,minexpD,negepD,ngrdD;
    double epsD,epsnegD,xmaxD,xminD;
    
    // result of machar() for float
    macharFloat(&ibeta,&it,&irnd,&ngrd,&machep,&negep,&iexp,&minexp,&maxexp,
                &eps,&epsneg,&xmin,&xmax);
    printf("//---for Float---//\n");
    printf("ibeta = %d\n",ibeta);
    printf("it = %d\n",it);
    printf("irnd = %d\n",irnd);
    printf("ngrd = %d\n",ngrd);
    printf("machep = %d\n",machep);
    printf("negep = %d\n",negep);
    printf("iexp = %d\n",iexp);
    printf("minexp = %d\n",minexp);
    printf("maxexp = %d\n",maxexp);
    printf("eps = %12.6g\n",eps);
    printf("epsneg = %12.6g\n",epsneg);
    printf("xmin = %12.6g\n",xmin);
    printf("xmax = %12.6g\n",xmax);
    //printf("//---------------//\n\n\n");
    
    //result of machar() for double
    macharDouble(&ibetaD,&itD,&irndD,&ngrdD,&machepD,&negepD,&iexpD,&minexpD,&maxexpD,
                 &epsD,&epsnegD,&xminD,&xmaxD);
    printf("//---for Double---//\n");
    printf("ibeta = %d\n",ibetaD);
    printf("it = %d\n",itD);
    printf("irnd = %d\n",irndD);
    printf("ngrd = %d\n",ngrdD);
    printf("machep = %d\n",machepD);
    printf("negep = %d\n",negepD);
    printf("iexp = %d\n",iexpD);
    printf("minexp = %d\n",minexpD);
    printf("maxexp = %d\n",maxexpD);
    printf("eps = %12.6g\n",epsD);
    printf("epsneg = %12.6g\n",epsnegD);
    printf("xmin = %12.6g\n",xminD);
    printf("xmax = %12.6g\n",xmaxD);
    //printf("//---------------//\n");
    
    //result of self-made get_eps()
    printf("//---------------'N' by 'get_eps()'--------------//\n");
    printf("'N' epsilon for float using get_eps_float(): %d\n", get_eps_float());
    printf("'N' epsilon for float using get_eps_Double(): %d\n", get_eps_double());
    printf("//---------------------------------------------------//\n");
    return 0;
}


