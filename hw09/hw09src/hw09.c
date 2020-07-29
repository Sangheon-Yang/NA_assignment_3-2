//
//  main.c
//  NA_hw09
//
//  Created by 양상헌 on 22/11/2019.
//  Copyright © 2019 양상헌. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"
#define NP 77
#define MP 20
#define MAXSTR 80
int main(void)
{
    int i,j,k,N,f;
    float *S_Vec123, *S_Vec456, *X, *Y, *Xprime, *Yprime, **S_Mat, **S_Mat_, **S_RVmat;
    //char dummy[MAXSTR];
    FILE *fp[3];
    fp[0]= fopen("fitdata1.dat","r");
    fp[1]= fopen("fitdata2.dat","r");
    fp[2]= fopen("fitdata3.dat","r");
    
    X = vector(1, NP);
    Y = vector(1, NP);
    Xprime = vector(1, NP);
    Yprime = vector(1, NP);
    S_Mat = matrix(1,3,1,3);
    S_Mat_ = matrix(1,3,1,3);
    S_Vec123 = vector(1,3);
    S_Vec456 = vector(1,3);
    S_RVmat = matrix(1,3,1,2);
    
    N = NP;
    for(f=0;f<3;f++){
        printf("\n>-------------------------FITDATA (%d)--------------------------<\n",f+1);
        if (fp[f] == NULL)
            nrerror("Data file matrx1.dat not found\n");
        
        //Saving Data
        for (k=1;k<=N;k++){
            fscanf(fp[f],"%f",&X[k]);
            fscanf(fp[f],"%f",&Y[k]);
            fscanf(fp[f],"%f",&Xprime[k]);
            fscanf(fp[f],"%f",&Yprime[k]);
        }
        
        //print data
        printf("\nTotal number of data: %d\n", N);
        printf(">List of data below:\n        X            Y            X'           Y' \n");
        for (k=1;k<=N;k++){
            printf("%12.6f %12.6f %12.6f %12.6f\n", X[k],Y[k],Xprime[k],Yprime[k]);
        }
        
        // initialize S_Matrix & S_Vec123, S_Vec456
        for (i=1;i<=3;i++){
            S_Vec123[i] = S_Vec456[i] = 0.0;
            for (j=1;j<=3;j++)
                S_Mat[i][j] = 0.0;
        }
        
        // initialize Sigma Matrices
        for (k=1;k<=N;k++){
            // S_Matrix
            S_Mat[1][1] += X[k]*X[k];//x*x
            S_Mat[1][2] += X[k]*Y[k];//x*y
            S_Mat[1][3] += X[k];//x
            S_Mat[2][2] += Y[k]*Y[k];//y*y
            S_Mat[2][3] += Y[k];//y
            S_Mat[3][3] += 1.0;//1.0
            // S_Vector for a1, a2, a3
            S_Vec123[1] += Xprime[k]*X[k];//x'*x
            S_Vec123[2] += Xprime[k]*Y[k];//x'*y
            S_Vec123[3] += Xprime[k];//x'
            // S_Vector for a4, a5, a6
            S_Vec456[1] += Yprime[k]*X[k];//y'*x
            S_Vec456[2] += Yprime[k]*Y[k];//y'*y
            S_Vec456[3] += Yprime[k];//y'
        }
        
        // S_Mat is a Symmetric Matrix
        S_Mat[2][1] = S_Mat[1][2];
        S_Mat[3][1] = S_Mat[1][3];
        S_Mat[3][2] = S_Mat[2][3];
        
        // initialize 3x2 Sigma Matrix
        S_RVmat[1][1] = S_Vec123[1];
        S_RVmat[1][2] = S_Vec456[1];
        S_RVmat[2][1] = S_Vec123[2];
        S_RVmat[2][2] = S_Vec456[2];
        S_RVmat[3][1] = S_Vec123[3];
        S_RVmat[3][2] = S_Vec456[3];
        
        //printing Formula
        printf("\n <Formula>: \n");
        for(i=1;i<=3;i++){
            printf(" |%15.6f %15.6f %15.6f|", S_Mat[i][1],S_Mat[i][2],S_Mat[i][3]);
            if(i==2)
                printf("*");
            else
                printf(" ");
            printf("|a%d a%d|" ,i,i+3 );
            if(i==2)
                printf(" = ");
            else
                printf("   ");
            printf("|%15.6f %15.6f|", S_Vec123[i], S_Vec456[i]);
            printf("\n");
        }
        
        // get Solution by using inverse matrix of S_Mat
        gaussj(S_Mat,3,S_RVmat,2);
        
        // print Solution
        printf("\n solution: \n");
        printf(" for X' : %12.6f * X + %12.6f * Y + %12.6f\n", S_RVmat[1][1], S_RVmat[2][1], S_RVmat[3][1]);
        printf(" for Y' : %12.6f * X + %12.6f * Y + %12.6f\n", S_RVmat[1][2], S_RVmat[2][2], S_RVmat[3][2]);
        printf("\n > a1 = %12.6f\n",S_RVmat[1][1]);
        printf(" > a2 = %12.6f\n",S_RVmat[2][1]);
        printf(" > a3 = %12.6f\n",S_RVmat[3][1]);
        printf(" > a4 = %12.6f\n",S_RVmat[1][2]);
        printf(" > a5 = %12.6f\n",S_RVmat[2][2]);
        printf(" > a6 = %12.6f\n\n",S_RVmat[3][2]);
        printf(">------------------------------------------------------------------<\n\n");
        fclose(fp[f]);
    }
    free_matrix(S_RVmat,1,3,1,2);
    free_vector(S_Vec456, 1, NP);
    free_vector(S_Vec123, 1, NP);
    free_matrix(S_Mat_,1,NP,1,NP);
    free_matrix(S_Mat,1,NP,1,NP);
    free_vector(Yprime, 1, NP);
    free_vector(Xprime, 1, NP);
    free_vector(Y, 1, NP);
    free_vector(X, 1, NP);
    return 0;
}
#undef NRANSI
