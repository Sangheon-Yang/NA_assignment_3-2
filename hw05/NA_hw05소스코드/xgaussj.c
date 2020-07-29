
/* Driver for routine gaussj */

#include <stdio.h>
#include <stdlib.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define NP 20
#define MP 20
#define MAXSTR 80

int main(void)
{
	int i,j,k,l,m,n;
	float **a,**ai,**u,**b,**x,**t;
	//char dummy[MAXSTR];
	FILE *fp[3];

	a=matrix(1,NP,1,NP);
	ai=matrix(1,NP,1,NP);
	u=matrix(1,NP,1,NP);
	b=matrix(1,NP,1,MP);
	x=matrix(1,NP,1,MP);
	t=matrix(1,NP,1,MP);
    if ((fp[0] = fopen("lineq1.dat","r")) == NULL)
        nrerror("Data file matrx1.dat not found\n");
    if ((fp[1] = fopen("lineq2.dat","r")) == NULL)
        nrerror("Data file matrx1.dat not found\n");
    if ((fp[2] = fopen("lineq3.dat","r")) == NULL)
        nrerror("Data file matrx1.dat not found\n");
    
    for(i=1;i<3;i++){
	printf("\n*************<lineq%d.dat> ***********************************\n",i+1);
	//while (!feof(fp)) {
		//fgets(dummy,MAXSTR,fp);
		//fgets(dummy,MAXSTR,fp);
		fscanf(fp[i],"%d",&n);
    fscanf(fp[i],"%d",&n);
    
		//fgets(dummy,MAXSTR,fp);
        m = 1;
    
		for (k=1;k<=n;k++)
			for (l=1;l<=n;l++) fscanf(fp[i],"%f ",&a[k][l]);
		//fgets(dummy,MAXSTR,fp);
		for (l=1;l<=m;l++)
			for (k=1;k<=n;k++) fscanf(fp[i],"%f ",&b[k][l]);
		/* save matrices for later testing of results */
		for (l=1;l<=n;l++) {
			for (k=1;k<=n;k++) ai[k][l]=a[k][l];
			for (k=1;k<=m;k++) x[l][k]=b[l][k];
		}
    //----print original A & B---//
    printf("\nOriginal matrix a : \n");
    for (k=1;k<=n;k++) {
        for (l=1;l<=n;l++) printf("%12.6f",a[k][l]);
        printf("\n");
    }
    
    printf("\nOriginal matrix b : \n");
    for (k=1;k<=n;k++) {
        printf("%12.6f",b[k][1]);
        printf("\n");
    }
		/* invert matrix */
		gaussj(ai,n,x,m);
		printf("\nInverse of matrix a : \n");
		for (k=1;k<=n;k++) {
			for (l=1;l<=n;l++) printf("%12.6f",ai[k][l]);
			printf("\n");
		}
		/* check inverse */
		printf("\na times a-inverse:\n");
		for (k=1;k<=n;k++) {
			for (l=1;l<=n;l++) {
				u[k][l]=0.0;
				for (j=1;j<=n;j++)
					u[k][l] += (a[k][j]*ai[j][l]);
			}
			for (l=1;l<=n;l++) printf("%12.6f",u[k][l]);
			printf("\n");
		}
    
        //--------edited--------//
    printf("\n<Solution vector X> : \n");
    for(k=1;k<=n;k++){
        printf("%12.6f\n",x[k][1]);
    }
        //----------------------//
    
		/* check vector solutions */
		printf("\nCheck the following for equality:\n");
		printf("%21s %14s\n","original","matrix*sol'n");
		for (l=1;l<=m;l++) {
			printf("vector %2d: \n",l);
			for (k=1;k<=n;k++) {
				t[k][l]=0.0;
				for (j=1;j<=n;j++)
					t[k][l] += (a[k][j]*x[j][l]);
				printf("%8s %12.6f %12.6f\n"," ",
					b[k][l],t[k][l]);
			}
		}
		printf("*************************************************************\n");
	//	printf("press RETURN for next problem:\n");
	//	(void) getchar();
	//}
	fclose(fp[i]);
    }
	free_matrix(t,1,NP,1,MP);
	free_matrix(x,1,NP,1,MP);
	free_matrix(b,1,NP,1,MP);
	free_matrix(u,1,NP,1,NP);
	free_matrix(ai,1,NP,1,NP);
	free_matrix(a,1,NP,1,NP);
	return 0;
}
#undef NRANSI
