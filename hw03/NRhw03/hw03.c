#include "nr.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double func01(double x){
    return (10.0*exp(-x)*sin(2.0*M_PI*x) - 2.0);
}

double func02(double x){
    return ((x-exp(-x))*(x-exp(-x)));
}

double func03(double x){
    return (cos(x+sqrt(2))+x*(x/2 + sqrt(2)));
}

double funcMyOwn(double x){
    return 0.0;
}

int main(int agrc, char** argv){
    printf("start\n");
    
	return 0;
}
