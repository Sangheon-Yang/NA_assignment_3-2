#include <stdio.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"
#include <time.h>
#define NP 11

int main(void)
{
    int i,j,k,kk,l,ll,nrot;
    float *d,*r,**v,**e;
    float a[NP][NP] = {0.0, };
    long idum = -523;
    clock_t T;
    
    // random한 number 로 11x11 symmetric matrix 생성
    // (Gaussian distribution) N(0,1)
    for(i= 0; i < NP; i++){
        for(j=i;j<NP;j++){
            T = clock();
            idum = (long)T + idum;
            a[i][j] = a[j][i]= gasdev(&idum);
        }
    }

    //생성된 Matrix 출력
    printf("\nRANDOM NUMBER MATRIX : \n");
    for(i = 0;i<NP;i++){
        for(j = 0;j<NP;j++){
            printf("%10.6f", a[i][j]);
        }
        printf("\n");
    }

    d=vector(1,NP);
    r=vector(1,NP);
    v=matrix(1,NP,1,NP);
    e=convert_matrix(&a[0][0],1,NP,1,NP);
    
    //jacobi transformation 실행
    jacobi(e,NP,d,v,&nrot);
    printf("\n number of JACOBI rotations: %3d\n\n",nrot);
    
    //eigenvalue를 descending order로 sort
    eigsrt(d,v,NP);
    
    printf("eigenvalues in descending order: \n");
    for (j=1;j<=NP;j++) {
        printf("%12.6f\n",d[j]);
    }
    
    printf("\neigenvectors:\n");
    for (j=1;j<=NP;j++) {
        printf("%9s %3d \n","number",j);
        for (k=1;k<=NP;k++) {
            printf("%12.6f\n",v[k][j]);
        }
        printf("\n");
    }
    
    printf("eigenvector test\n");
    for (j=1;j<=NP;j++) {
        for (l=1;l<=NP;l++) {
            r[l]=0.0;
            for (k=1;k<=NP;k++) {
                if (k > l) {
                    kk=l;
                    ll=k;
                } else {
                    kk=k;
                    ll=l;
                }
                r[l] += (e[ll][kk]*v[k][j]);
            }
        }
        printf("vector number %3d\n",j);
        printf("%11s %14s %10s %14s\n",
            "vector","mtrx*vec.","ratio", "eigenValue");
        for (l=1;l<=NP;l++)
            printf("%12.6f %12.6f %12.6f %12.6f\n",
                v[l][j],r[l],r[l]/v[l][j],d[j]);
    }
    
    free_convert_matrix(e,1,NP,1,NP);
    free_matrix(v,1,NP,1,NP);
    free_vector(r,1,NP);
    free_vector(d,1,NP);
    return 0;
}
#undef NRANSI
