#include "wgmsolver.h"

/*
 * sparams.c
 *
 *  Created on: Nov 25, 2011
 *      Author: alberto
 */


double complex *sparams(){

}

/**
 * simple matrix multiply for matrices m1xn1xNfreqs by m2xn2xNfreqs,
 * defined as m1xn1 times m2xn2 for every freq
 * n1 must be equal to m2.
 */
double complex *mmatrix( int Nfreqs, int m1, int n1,double complex matrix1[m1][n1][Nfreqs], int m2, int n2,double complex matrix2[m2][n2][Nfreqs]){
	double complex retMatrix[][]=(double complex(*)[m1][n2][Nfreqs])malloc(m1*n2*Nfreqs*sizeof(double complex));
	for(int i=0;i<m1;i++){
		double complex partsum[Nfreqs]=0;
		for(int j=0;j<n2;j++){
			for(int k=0;k<Nfreqs;k++){
				retMatrix[i][j][k]=matrix1[i][j][k]*matrix2[i][j][k];
			}
		}
	}
}
