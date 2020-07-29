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
    int i,j,k,N,f,iteration, NotUpdated;
    float *X, *Y, *Xprime, *Yprime;
    float **Grad, **A_Mat, **Ai_Mat, *tempPD, *deltaA;
    float **H_Mat, **H_Mat_Save;
    float rambda, ErrCur, ErrNew, initialErr, tmpXsqr, tmpYsqr;
    float tempX, tempY, tempXmX, tempYmY;
    float tempPD31_1, tempPD31_2, tempPD32_1, tempPD32_2;
    
    FILE *fp[8];
    fp[0]= fopen("01.MatchingPoints_Bad_without_Noise.txt","r");
    fp[1]= fopen("02.MatchingPoints_Bad_with_GaussianNoise_SD=1.txt","r");
    fp[2]= fopen("03.MatchingPoints_Bad_with_GaussianNoise_SD=10.txt","r");
    fp[3]= fopen("04.MatchingPoints_Bad_with_GaussianNoise_SD=30.txt","r");
    fp[4]= fopen("05.MatchingPoints_Good_without_Noise.txt","r");
    fp[5]= fopen("06.MatchingPoints_Good_with_GaussianNoise_SD=1.txt","r");
    fp[6]= fopen("07.MatchingPoints_Good_with_GaussianNoise_SD=10.txt","r");
    fp[7]= fopen("08.MatchingPoints_Good_with_GaussianNoise_SD=30.txt","r");
    
    for(f=0;f<8;f++){
        if(f == 0)
            printf("\n>------- 01. Matching Points of Bad LEFT RIGHT without Gaussian Noise  --------------------<\n");
        else if(f == 1)
            printf("\n>------- 02. Matching Points of Bad LEFT RIGHT with Gaussian Noise SD = 1 -----------------<\n");
        else if(f == 2)
            printf("\n>------- 03. Matching Points of Bad LEFT RIGHT with Gaussian Noise SD = 10 ----------------<\n");
        else if(f == 3)
            printf("\n>------- 04. Matching Points of Bad LEFT RIGHT with Gaussian Noise SD = 30 ----------------<\n");
        else if(f == 4)
            printf("\n>------- 05. Matching Points of Good LEFT RIGHT without Gaussian Noise  -------------------<\n");
        else if(f == 5)
            printf("\n>------- 06. Matching Points of Good LEFT RIGHT with Gaussian Noise SD = 1 ----------------<\n");
        else if(f == 6)
            printf("\n>------- 07. Matching Points of Good LEFT RIGHT with Gaussian Noise SD = 10 ---------------<\n");
        else //if(f == 7)
            printf("\n>------- 08. Matching Points of Good LEFT RIGHT with Gaussian Noise SD = 30 ---------------<\n");
        
        if (fp[f] == NULL)
            nrerror("Data file is not found\n");
        
        fscanf(fp[f],"%d",&N);
        printf("Number of Matching Points: %d\n",N);
        
        X = vector(1, N);
        Y = vector(1, N);
        Xprime = vector(1, N);
        Yprime = vector(1, N);
        Grad = matrix(1,8,1,1);
        A_Mat = matrix(1,3,1,3);
        Ai_Mat =matrix(1,3,1,3);
        tempPD = vector(1,8);
        deltaA = vector(1,8);
        H_Mat = matrix(1,8,1,8);
        H_Mat_Save = matrix(1,8,1,8);
        
        //Saving Data set
        for (k=1;k<=N;k++){
            fscanf(fp[f],"%f",&X[k]);
            fscanf(fp[f],"%f",&Y[k]);
            fscanf(fp[f],"%f",&Xprime[k]);
            fscanf(fp[f],"%f",&Yprime[k]);
        }
        
        //print data set
        printf("\nTotal number of data: %d\n", N);
        printf(">List of data below:\n        X            Y            X'           Y' \n");
        for (k=1;k<=N;k++){
            printf("%12.6f %12.6f %12.6f %12.6f\n", X[k],Y[k],Xprime[k],Yprime[k]);
        }
        printf("\n");
        
        // [1] - inititial guess of Parameter Set A
        for(i=1;i<=3;i++)
            for(j=1;j<=3;j++)
                A_Mat[i][j] = 1.0;
       
        // [2] - calculating initial error
        initialErr = 0.0;
        for(i=1;i<=N ;i++){
            tmpXsqr = (Xprime[i] - ((A_Mat[1][1]*X[i]+A_Mat[1][2]*Y[i]+A_Mat[1][3])/(A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0)));
            tmpYsqr = (Yprime[i] - ((A_Mat[2][1]*X[i]+A_Mat[2][2]*Y[i]+A_Mat[2][3])/(A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0)));
            initialErr += (tmpXsqr*tmpXsqr+tmpYsqr*tmpYsqr);
        }
        ErrNew = initialErr;
        
        // [3] - initialize rambda value to 0.001
        rambda = 0.001;
        
        //iteration flague for counting
        iteration = 0;
        NotUpdated = 0;
        
        // [4] - ---------Start of Iteration------------//
        while( rambda > 0 && NotUpdated < 10){
            iteration++;
            ErrCur = 0.0;

            // [4]-(1)- initialize Hessian and Gradient to 0.0;
            for(k=1;k<=8;k++){
                Grad[k][1] = 0.0;
                for(j=1;j<=8;j++){
                    H_Mat[k][j] = 0.0;
                }
            }
            
            // [4]-(2) calculate Gradient and Hessian with current A;
            for(i=1;i<=N;i++){
                // for error calculating
                tempX = ((A_Mat[1][1]*X[i]+A_Mat[1][2]*Y[i]+A_Mat[1][3])/(A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0));
                tempY = ((A_Mat[2][1]*X[i]+A_Mat[2][2]*Y[i]+A_Mat[2][3])/(A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0));
                tempXmX = (Xprime[i] - tempX);
                tempYmY = (Yprime[i] - tempY);
                
                // for dy/dak
                tempPD[1] = X[i]/(A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0);
                tempPD[2] = Y[i]/(A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0);
                tempPD[3] = 1.0/(A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0);
                tempPD[4] = tempPD[1]; //X[i]/(A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0);
                tempPD[5] = tempPD[2]; //Y[i]/(A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0);
                tempPD[6] = tempPD[3]; //1.0/(A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0);
                tempPD31_1 = -1.0*tempPD[1]*tempX;
                //(-1.0*X[i]*(A_Mat[1][1]*X[i]+A_Mat[1][2]*Y[i]+A_Mat[1][3]))/((A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0)*(A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0));
                tempPD31_2 = -1.0*tempPD[1]*tempY;
                //(-1.0*X[i]*(A_Mat[2][1]*X[i]+A_Mat[2][2]*Y[i]+A_Mat[2][3]))/((A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0)*(A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0));
                tempPD[7] = tempPD31_1+tempPD31_2;
                // = -1.0*tempPD[1]
                tempPD32_1 = -1.0*tempPD[2]*tempX;
                //(-1.0*Y[i]*(A_Mat[1][1]*X[i]+A_Mat[1][2]*Y[i]+A_Mat[1][3]))/((A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0)*(A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0));
                tempPD32_2 = -1.0*tempPD[2]*tempY;
                //(-1.0*Y[i]*(A_Mat[2][1]*X[i]+A_Mat[2][2]*Y[i]+A_Mat[2][3]))/((A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0)*(A_Mat[3][1]*X[i]+A_Mat[3][2]*Y[i]+1.0));
                tempPD[8] = tempPD32_1+tempPD32_2;
                
                //Gradient;
                Grad[1][1] += (-2.0*tempXmX*tempPD[1]);
                Grad[2][1] += (-2.0*tempXmX*tempPD[2]);
                Grad[3][1] += (-2.0*tempXmX*tempPD[3]);
                Grad[4][1] += (-2.0*tempYmY*tempPD[4]);
                Grad[5][1] += (-2.0*tempYmY*tempPD[5]);
                Grad[6][1] += (-2.0*tempYmY*tempPD[6]);
                Grad[7][1] += (-2.0*(tempXmX*tempPD31_1+tempYmY*tempPD31_2));
                Grad[8][1] += (-2.0*(tempXmX*tempPD32_1+tempYmY*tempPD32_2));
            
                //calculate Hessian in current Iteration step
                for(j=1;j<=8;j++){
                    for(k=1;k<=8;k++){
                        H_Mat[j][k] += (2.0*tempPD[j]*tempPD[k]);
                    }
                }
                //calculating current error
                ErrCur += (tempXmX*tempXmX + tempYmY*tempYmY);
            }
            
            //make -gradient for Matrix operation
            for(i=1;i<=8;i++)
                Grad[i][1] *= -1.0;
            
            //H' = H + rambda*I
            for(i=1;i<=8;i++)
                H_Mat[i][i] += rambda;
            
            // Grad -> new delta A
            // H_Mat -> H_Mat^-1
            gaussj(H_Mat,8,Grad,1);
            for(i=1;i<=8;i++){
                deltaA[i]=Grad[i][1];
            }
            printf("%7dth iteration, [RAMBDA]: %20.15f,  [Err]: %12.6f  ",iteration,rambda,ErrCur);
            
            //ERROR 비교
            Ai_Mat[1][1] = A_Mat[1][1]+deltaA[1];
            Ai_Mat[1][2] = A_Mat[1][2]+deltaA[2];
            Ai_Mat[1][3] = A_Mat[1][3]+deltaA[3];
            Ai_Mat[2][1] = A_Mat[2][1]+deltaA[4];
            Ai_Mat[2][2] = A_Mat[2][2]+deltaA[5];
            Ai_Mat[2][3] = A_Mat[2][3]+deltaA[6];
            Ai_Mat[3][1] = A_Mat[3][1]+deltaA[7];
            Ai_Mat[3][2] = A_Mat[3][2]+deltaA[8];
            ErrNew = 0.0;
            for(i=1;i<=N;i++){
                tempX = ((Ai_Mat[1][1]*X[i]+Ai_Mat[1][2]*Y[i]+Ai_Mat[1][3])/(Ai_Mat[3][1]*X[i]+Ai_Mat[3][2]*Y[i]+1.0));
                tempY = ((Ai_Mat[2][1]*X[i]+Ai_Mat[2][2]*Y[i]+Ai_Mat[2][3])/(Ai_Mat[3][1]*X[i]+Ai_Mat[3][2]*Y[i]+1.0));
                tempXmX = (Xprime[i] - tempX);
                tempYmY = (Yprime[i] - tempY);
                ErrNew += (tempXmX*tempXmX + tempYmY*tempYmY);
            }
            
            if(ErrNew < ErrCur){
                rambda *= 0.1;
                //update A_MAt[][];
                A_Mat[1][1] += deltaA[1];
                A_Mat[1][2] += deltaA[2];
                A_Mat[1][3] += deltaA[3];
                A_Mat[2][1] += deltaA[4];
                A_Mat[2][2] += deltaA[5];
                A_Mat[2][3] += deltaA[6];
                A_Mat[3][1] += deltaA[7];
                A_Mat[3][2] += deltaA[8];
                printf(",  A is updated\n");
                NotUpdated = 0;
            }else{
                rambda *= 10.0;
                NotUpdated++;
                printf("\n");
            }
        
        //--while close
        }
        // print Solution
        printf("\n solution: ");
        printf("\n > a11 = %12.6f",A_Mat[1][1]);
        printf("\n > a12 = %12.6f",A_Mat[1][2]);
        printf("\n > a13 = %12.6f",A_Mat[1][3]);
        printf("\n > a21 = %12.6f",A_Mat[2][1]);
        printf("\n > a22 = %12.6f",A_Mat[2][2]);
        printf("\n > a23 = %12.6f",A_Mat[2][3]);
        printf("\n > a31 = %12.6f",A_Mat[3][1]);
        printf("\n > a32 = %12.6f\n\n",A_Mat[3][2]);
        printf(">-------------------------------------------------------------------------------<\n\n");
        fclose(fp[f]);
        
        free_matrix(H_Mat_Save, 1, 8, 1, 8);
        free_matrix(H_Mat, 1, 8, 1, 8);
        free_vector(deltaA, 1, 8);
        free_vector(tempPD, 1, 8);
        free_matrix(Ai_Mat, 1, 3, 1, 3);
        free_matrix(A_Mat, 1, 3, 1, 3);
        free_matrix(Grad, 1, 8, 1, 1);
        
        free_vector(Yprime, 1, NP);
        free_vector(Xprime, 1, NP);
        free_vector(Y, 1, NP);
        free_vector(X, 1, NP);
    }
    return 0;
}
#undef NRANSI
