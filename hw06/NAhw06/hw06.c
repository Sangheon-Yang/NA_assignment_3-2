//
//  hw06.c
//  NA_HW06
//
//  Created by 양상헌 on 31/10/2019.
//  Copyright © 2019 양상헌. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "nr.h"
int main(){
    int i , j, k, inBound;
    int sample[4] = { 100, 1000, 10000, 100000 };
    float sampleDiv[4] = {1.0,1.0,2.0,20.0};
    int distributionU[4][100] = {0, };
    int distributionG[4][102] = {0, };
    int temp = 0;
    float temp1 = 0.0;
    float tempF = 0.0;
    
    
    clock_t T;
    long Te = -500;
    float rN, erN;
    
    //uniform distribution
    printf("\n ----< Uniform Distribution >---- \n\n");
    for(j = 0; j < 4; j++){
        printf("> Sample : %d \n", sample[j]);
        for(i = 0; i< sample[j];i++){
            T = clock();
            Te = (long)T + Te;
            rN = ran1(&Te);
            
            temp = (int)(rN*100);
            distributionU[j][temp]++;

        }
        printf("\n> distribution histogram : sample = %d , [-3,2] in 100 intervals \n", sample[j]);
        
        for(i = 0;i<100;i++){
            printf(" %4d", distributionU[j][i]);
            if(i%20 == 19){
                printf("\n");
            }
        }
        printf("\n Graph: \n");
        for(i = 0;i<100;i++){
            tempF = (float)distributionU[j][i]/(float)sampleDiv[j];
            tempF = tempF+0.5;//반올림
            printf("%5.2f~%5.2f : ", -3.0+0.05*i, -3.0+0.05*(i+1));
            for(k = 0; k<(int)tempF;k++){
                printf("*");
            }
            printf("\n");
        }
        printf("\n------------------------------------------------------------\n");
    }
    
    tempF = 0.0;
    Te = -500;
    //Gaussian distribution
    printf("\n ----< Gaussian Distribution >---- \n\n");
    for(j = 0; j < 4; j++){
        inBound = 0;
        printf("> Sample : %d \n", sample[j]);
        for(i = 0; i < sample[j] ;i++){
            T = clock();
            Te = (long)T + Te;
            rN = gasdev(&Te);
            
            erN = rN*1.5-0.5;
            
            temp1 = (erN+3.0)*20;
            
            if(temp1 < 0){
                distributionG[j][100]++;//out of bound negative
            }
            else if(temp1 >= 100){
                distributionG[j][101]++;//out of bound positive
            }
            else{
                distributionG[j][(int)temp1]++;
                inBound++;
            }
        }
        printf("\n> distribution histogram : sample = %d , [-3,2] in 100 intervals \n", sample[j]);
        for(i = 0;i<100;i++){
            printf(" %4d", distributionG[j][i]);
            if(i%20 == 19){
                printf("\n");
            }
        }
        printf("in bound [-3,2]: %d\n",inBound);
        printf("out of negative bound: %d && out of positive bound: %d\n",distributionG[j][100], distributionG[j][101]);
        printf("\n\n Graph: \n");
        for(i = 0;i<100;i++){
            tempF = (float)distributionG[j][i]/(float)sampleDiv[j];
            tempF = tempF+0.5;//반올림
            printf("%5.2f~%5.2f : ", -3.0+0.05*i, -3.0+0.05*(i+1));
            for(k = 0; k<(int)tempF;k++){
                printf("*");
            }
            printf("\n");
        }
        printf("\n--------------------------------------------------------------\n");
    }
    
    return 0;
}


