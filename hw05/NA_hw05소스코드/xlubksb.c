
/* Driver for routine lubksb */

#include <stdio.h>
#include <stdlib.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define NP 20
#define MAXSTR 80

int main(void)
{
	int i,j,k,l,m,n,*indx;
	float determ,p,*x,**a,**b,**c,**xu,**xl,*bb,*bbb;
    long idum=(-13);
	//char dummy[MAXSTR];
	FILE *fp[3];
    indx=ivector(1,NP);
    x=vector(1,NP);
    a=matrix(1,NP,1,NP);
    b=matrix(1,NP,1,NP);
    c=matrix(1,NP,1,NP);
    xu=matrix(1, NP, 1, NP);
    xl=matrix(1, NP, 1, NP);
    bb=vector(1,NP);
    bbb=vector(1,NP);
	
    if ((fp[0] = fopen("lineq1.dat","r")) == NULL)
        nrerror("Data file matrx1.dat not found\n");
    if ((fp[1] = fopen("lineq2.dat","r")) == NULL)
        nrerror("Data file matrx1.dat not found\n");
    if ((fp[2] = fopen("lineq3.dat","r")) == NULL)
        nrerror("Data file matrx1.dat not found\n");
    for(i=0;i<3;i++){
        printf("\n*************<lineq%d.dat> ***********************************\n",i+1);
	//while (!feof(fp)) {
		//fgets(dummy,MAXSTR,fp);
		//fgets(dummy,MAXSTR,fp);
		fscanf(fp[i],"%d ",&n);
        fscanf(fp[i],"%d ",&n);
        m=1;
		//fgets(dummy,MAXSTR,fp);
		for (k=1;k<=n;k++)
			for (l=1;l<=n;l++) fscanf(fp[i],"%f ",&a[k][l]);
		//fgets(dummy,MAXSTR,fp);
		for (l=1;l<=m;l++)
			for (k=1;k<=n;k++) fscanf(fp[i],"%f ",&b[k][l]);
		/* Save matrix a for later testing */
		for (l=1;l<=n;l++)
			for (k=1;k<=n;k++) c[k][l]=a[k][l];
        
        
        printf("\nOriginal matrix a : \n");
        for (k=1;k<=n;k++) {
            for (l=1;l<=n;l++) printf("%12.6f",a[k][l]);
            printf("\n");
        }
        
        printf("\nOriginal matrix b : \n");
        for (k=1;k<=n;k++) {
            bb[k]=bbb[k] =b[k][1];
            printf("%12.6f",b[k][1]);
            printf("\n");
        }
        
        
		/* Do LU decomposition */
		ludcmp(c,n,indx,&p);
		/* Solve equations for each right-hand vector */
        
        /* Compose separately the lower and upper matrices */
        for (k=1;k<=n;k++) {
            for (l=1;l<=n;l++) {
                if (l > k) {
                    xu[k][l]=c[k][l];
                    xl[k][l]=0.0;
                } else if (l < k) {
                    xu[k][l]=0.0;
                    xl[k][l]=c[k][l];
                } else {
                    xu[k][l]=c[k][l];
                    xl[k][l]=1.0;
                }
            }
        }
        
        printf("\nlower matrix of the decomposition:\n");
        for (k=1;k<=n;k++) {
            for (l=1;l<=n;l++) printf("%12.6f",xl[k][l]);
            printf("\n");
        }
        printf("\nupper matrix of the decomposition:\n");
        for (k=1;k<=n;k++) {
            for (l=1;l<=n;l++) printf("%12.6f",xu[k][l]);
            printf("\n");
        }
        determ = 1.0;
        for(k=1;k<=n;k++){
            determ *= xu[k][k];
        }
        printf("\nDeterminant of Matrix a : %12.6f\n",determ);
        
		for (k=1;k<=m;k++) {
            
			for (l=1;l<=n;l++) x[l]=b[l][k];
			lubksb(c,n,indx,x);
            
            printf("\n<Solution vector X> : \n");
            for(k=1;k<=n;k++){
                printf("%12.6f\n",x[k]);
            }
            printf("\n");
            
			/* Test results with original matrix */
			printf("right-hand side vector:\n");
			for (l=1;l<=n;l++)
				printf("%12.6f",b[l][1]);
            
			printf("\n\n%s%s\n","result of matrix applied",
				" to sol'n vector");
			for (l=1;l<=n;l++) {
				b[l][k]=0.0;
				for (j=1;j<=n;j++)
					b[l][k] += (a[l][j]*x[j]);
			}
			for (l=1;l<=n;l++)
				printf("%12.6f",b[l][k]);
            
            
            for (l=1;l<=n;l++) x[l] *= (1.0+0.2*ran3(&idum));
            printf("\n\nSolution vector with noise added:\n");
            for (l=1;l<=n;l++) printf("%12.6f",x[l]);
            printf("\n");
            mprove(a,c,n,indx,bb,x);
            printf("\nSolution vector recovered by mprove:\n");
            for (l=1;l<=n;l++) printf("%12.6f",x[l]);
            printf("\n");
            
            printf("\n\n%s%s\n","result of matrix applied",
                " to sol'n vector with mprove");
            for (l=1;l<=n;l++) {
                bbb[l]=0.0;
                for (j=1;j<=n;j++)
                    bbb[l] += (a[l][j]*x[j]);
            }
            for (l=1;l<=n;l++)
                printf("%12.6f",b[l][k]);
            printf("\n*******************************************************\n\n");
		}
		//printf("press RETURN for next problem:\n");
		//(void) getchar();
	//}
        fclose(fp[i]);
        
    }
    free_vector(bb,1,NP);
    free_vector(bbb,1,NP);
    free_matrix(c,1,NP,1,NP);
    free_matrix(b,1,NP,1,NP);
    free_matrix(a,1,NP,1,NP);
    free_vector(x,1,NP);
    free_ivector(indx,1,NP);
    free_matrix(xu,1,NP,1,NP);
    free_matrix(xl,1,NP,1,NP);
	return 0;
}
#undef NRANSI
