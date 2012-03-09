#include "wgmsolver.h"

/*
 * sparams.c
 *
 *  Created on: Nov 25, 2011
 *      Author: alberto
 */


double complex *sparams(){
return;
}

/**
 * simple matrix multiply for matrices m1xn1xNfreqs by m2xn2xNfreqs,
 * defined as m1xn1 times m2xn2 for every freq
 * n1 must be equal to m2.
 *
 * retMatrix MUST BE FREED AFTER RETURNING
 */
double complex *mmatrix( int Nfreqs, int m1, int n1,double complex matrix1[m1][n1][Nfreqs], int m2, int n2,double complex matrix2[m2][n2][Nfreqs]){
	double complex (*retMatrix)[m1][n2][Nfreqs]=(double complex(*)[m1][n2][Nfreqs])malloc(m1*n2*Nfreqs*sizeof(double complex));
	for(int i=0;i<m1;i++){
		double complex partsum[Nfreqs];
		memset(&partsum,0,Nfreqs);
		for(int j=0;j<n2;j++){
			for(int p=0;p<n1;p++){
				for(int k=0;k<Nfreqs;k++){
					partsum[k]+=matrix1[i][p][k]*matrix2[p][j][k];
				}
			}
			for(int k=0;k<Nfreqs;k++){
				(*retMatrix)[i][j][k]=partsum[k];
			}
		}
	}
	return (double complex *)retMatrix;
}

/**
 * Method to multiply the gammas by a dense matrix of size nxn
 * The gammas are not placed in a diagonal matrix, but in a unidimensional array
 * (it is equivalent, as the other elemnts are zero)
 * This method works in O(n^2)
 *
 * Params:
 * 	- Nfreqs Number of frequencies
 * 	- m number of rows
 * 	- n number of columns
 * 	- matrix mxnxNfreqs dense matrix
 * 	- gammas mxNfreqs matrix
 *
 * Returns a pointer to a mxnxNfreqs double complex array.
 */
double complex *mdiagonaltimes(int Nfreqs, int m, int n, double complex matrix[m][n][Nfreqs], double complex gammas[m][Nfreqs]){
	double complex (*retMatrix)[m][n][Nfreqs]=(double complex (*)[m][n][Nfreqs])malloc(m*n*Nfreqs*sizeof(double complex));
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			for(int k=0;k<Nfreqs;k++){
				(*retMatrix)[i][j][k]=gammas[j][k]*matrix[i][j][k];
			}
		}
	}
	return (double complex *)retMatrix;
}

/**
 * Method to multiply a dense matrix of size nxn times the gammas
 * The gammas are not placed in a diagonal matrix, but in a unidimensional array
 * (it is equivalent, as the other elemnts are zero)
 * This method works in O(n^2)
 *
 * Params:
 * 	- Nfreqs Number of frequencies
 * 	- m number of rows
 * 	- n number of columns
 * 	- matrix mxnxNfreqs dense matrix
 * 	- gammas mxNfreqs matrix
 *
 * Returns a pointer to a mxnxNfreqs double complex array.
 */
double complex *mtimesdiagonal(int Nfreqs, int m, int n, double complex matrix[m][n][Nfreqs], double complex gammas[n][Nfreqs]){
	double complex (*retMatrix)[m][n][Nfreqs]=(double complex (*)[m][n][Nfreqs])malloc(m*n*Nfreqs*sizeof(double complex));
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			for(int k=0;k<Nfreqs;k++){
				(*retMatrix)[i][j][k]=gammas[i][k]*matrix[i][j][k];
			}
		}
	}
	return (double complex *)retMatrix;
}


/**
 * Method to multiply gammas by a dense matrix of size nxn times other set of gammas gammas
 * The gammas are not placed in a diagonal matrix, but in a unidimensional array
 * (it is equivalent, as the other elemnts are zero)
 * This method works in O(n^2)
 *
 * Params:
 * 	- Nfreqs Number of frequencies
 * 	- m number of rows
 * 	- n number of columns
 * 	- matrix mxnxNfreqs dense matrix
 * 	- gammas1 mxNfreqs matrix
 * 	- gammas2 nxNfreqs matrix
 *
 * Returns a pointer to a mxnxNfreqs double complex array.
 */
double complex *mdiagonaltimesdiagonal(int Nfreqs, int m, int n, double complex matrix[m][n][Nfreqs],
									   double complex gammas1[m][Nfreqs], double complex gammas2[n][Nfreqs]){
	double complex (*retMatrix)[m][n][Nfreqs]=(double complex (*)[m][n][Nfreqs])malloc(m*n*Nfreqs*sizeof(double complex));
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			for(int k=0;k<Nfreqs;k++){
				(*retMatrix)[i][j][k]=gammas1[j][k]*gammas2[i][k]*matrix[i][j][k];
			}
		}
	}
	return (double complex *)retMatrix;
}

/*
 * LDU decomposition of a rectangular matrix
 */
double complex **LDUdecomposition(double complex **LDU,int m, int Nfreqs, double complex matrix[m][m][Nfreqs]){
	  LDU= (double complex **)malloc(3*sizeof(double complex));

	  //Memmory needed to store a triangular matrix n*(n+1)/2
	  double complex (*L)[m][m][Nfreqs]= (double complex (*)[m][m][Nfreqs])malloc(m*m*Nfreqs*sizeof(double complex));
	  memset(&L,0,m*m*Nfreqs*sizeof(double complex));
	  *LDU=(double complex *)L;
	  double complex (*D)[m][Nfreqs]= (double complex (*)[m][Nfreqs])malloc(m*Nfreqs*sizeof(double complex));
	  *(LDU+1)=(double complex *)D;
	  double complex (*U)[m][m][Nfreqs]= (double complex (*)[m][m][Nfreqs])malloc(m*m*Nfreqs*sizeof(double complex));
	  memset(&U,0,m*m*Nfreqs*sizeof(double complex));
	  *(LDU+1)=(double complex *)U;

//	  //Fill the diagonals of the triangular matrices
//	  for(int i=0;i<m;i++){
//		  for(int f=0;f<Nfreqs;f++){
//		    //(*L)[i][i][f]=1;
//		    //(*U)[i][i][f]=1;
//		    //(*D)[i][f]=matrix[i][i][f];
//		  }
//	  }

	  for(int i=0; i<m;i++){
		  for(int j=0;j<m;j++){
			  double complex partSum[m][Nfreqs];
			  memset(&partSum,0,m*m*Nfreqs*sizeof(double complex));
			  int p=0;
			  //partial sum (undo row times column)
			  while(p<((j<i)?j:i)){
				  for(int f=0;f<Nfreqs;f++){
					  partSum[j][f]+=(*L)[i][p][f]*(*U)[p][j][f];
				  }
				  p++;
			  }
			  //substract partial sum to get last element
			  if (j>i){ //<--this case upper matrix
				  for(int f=0; f<Nfreqs;f++){
					  (*U)[i][j][f]=matrix[i][j][f]-partSum[j][f];
				  }
			  }else if(i==j){ //<-- Lower AND diagonal is filled
				  for(int f=0; f<Nfreqs;f++){
					  (*D)[i][f]=matrix[i][j][f]-partSum[j][f];
					  (*L)[i][j][f]=1;
				  }
			  }else{ //<-- Lower AND normalized
				  for(int f=0; f<Nfreqs;f++){
					  (*L)[i][j][f]=(matrix[i][j][f]-partSum[j][f])/(*D)[i][f];
				  }
			  }
		  }
	  }
	  return LDU;
}

double complex* invDiagonal(int m, int Nfreqs, double complex matrix[m][m][Nfreqs]){


	return;
}
