#include "wgmsolver.h"

/*
 * sparams.c
 *
 *  Created on: Nov 25, 2011
 *      Author: alberto
 */

/**
 * simple matrix multiply for matrices m1xn1xNfreqs by m2xn2xNfreqs,
 * defined as m1xn1 times m2xn2 for every freq
 * n1 must be equal to m2.
 *
 * retMatrix MUST BE FREED AFTER RETURNING
 */
double complex *mmatrix( int m1, int n1, int Nfreqs,double complex matrix1[m1][n1][Nfreqs], int n2, int m2, double complex matrix2[m2][n2][Nfreqs]){
	double complex (*retMatrix)[m1][n2][Nfreqs]=(double complex(*)[m1][n2][Nfreqs])malloc(m1*n2*Nfreqs*sizeof(double complex));
	for(int i=0;i<m1;i++){
		double complex partsum[Nfreqs];
		memset(&partsum,0,Nfreqs*sizeof(double complex));
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
	return (double complex (*) [m1][n2][Nfreqs])retMatrix;
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
	memset(&retMatrix,0,Nfreqs*m*n*sizeof(double complex));
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
double complex **LUdecomposition(int m, int Nfreqs, double complex matrix[m][m][Nfreqs]){
	  //Lines needing editing for LDU instead of LU are marked with <------------------

	  //change to 2 to 3 for LDU
	  double complex **LDU= (double complex **)malloc(2*sizeof(double complex)); //<------------------------------

	  //Memmory needed to store a triangular matrix n*(n+1)/2
	  double complex (*L)[m][m][Nfreqs]= (double complex (*)[m][m][Nfreqs])malloc(m*m*Nfreqs*sizeof(double complex));
	  memset(&L,0,m*m*Nfreqs*sizeof(double complex));
	  *LDU=(double complex *)L;
	  //double complex (*D)[m][Nfreqs]= (double complex (*)[m][Nfreqs])malloc(m*Nfreqs*sizeof(double complex));
	  //*(LDU+1)=(double complex *)D;
	  double complex (*U)[m][m][Nfreqs]= (double complex (*)[m][m][Nfreqs])malloc(m*m*Nfreqs*sizeof(double complex));
	  memset(&U,0,m*m*Nfreqs*sizeof(double complex));
	  //change to LDU+2 for LDU
	  *(LDU+1)=(double complex *)U; //<-------------------------------------------------

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
			  if (j>=i){ //<--this case upper matrix
				  for(int f=0; f<Nfreqs;f++){
					  (*U)[i][j][f]=matrix[i][j][f]-partSum[j][f];
				  }
				  /*
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
			  */
			  }else{
				  for(int f=0; f<Nfreqs;f++){
					  (*L)[i][j][f]=(matrix[i][j][f]-partSum[j][f]);
				  }
			  }
		  }
	  }
	  return LDU;
}


/*
 * Inverse of an square matrix decomposed using a LU decomposition.
 * [L][U][A]^-1=[I]
 *
 */
double complex* LUinverse(int m, int Nfreqs, double complex ** LDU){
	// [L][Z]=[I]
	  double complex (*Z)[m][m][Nfreqs]= (double complex (*)[m][m][Nfreqs])malloc(m*m*Nfreqs*sizeof(double complex));
	  memset(&Z,0,m*m*Nfreqs*sizeof(double complex));

	  double complex (*L)[m][m][Nfreqs]=(double complex (*)[m][m][Nfreqs]) *(LDU);
	  //double complex (*D)[m][Nfreqs]=(double complex (*)[m][Nfreqs]) *(LDU+1);
	  double complex (*U)[m][m][Nfreqs]=(double complex (*)[m][m][Nfreqs]) *(LDU+2);

	  //column of the identity matrix
	  double complex Icol[m];
	  //partial sum
	  double complex partSum[m][Nfreqs];

	  for(int i=0;i<m;i++){
		  memset(&Icol, 0, m*sizeof(double complex));
		  Icol[i]=1;
		  for(int j=0;j<m;j++){
			  memset(&partSum, 0, m*Nfreqs*sizeof(double complex));
			  for(int k=0;k<j;k++){
				  for(int f=0;f<Nfreqs;f++){
					  partSum[j][f]+=(*L)[j][k][f]*(*Z)[k][i][f];
				  }
			  }
			  for(int f=0;f<Nfreqs;f++){
				  (*Z)[i][j][f]=(Icol[j]-partSum[j][f])/(*L)[j][j][f];
			  }
		  }
	  }

	  double complex (*X)[m][m][Nfreqs]= (double complex (*)[m][m][Nfreqs])malloc(m*m*Nfreqs*sizeof(double complex));
	  memset(&X,0,m*m*Nfreqs*sizeof(double complex));
	  //[U][X]=[Z]
	  for(int i=0;i<m;i++){
		  for(int j=m-1;j<0;j--){
			  memset(&partSum, 0, m*Nfreqs*sizeof(double complex));
			  for(int k=m-1;k>j;k--){
				  for(int f=0;f<Nfreqs;f++){
					  partSum[j][f]+=(*U)[j][k][f]*(*X)[k][i][f];
				  }
			  }
			  for(int f=0;f<Nfreqs;f++){
				  (*X)[i][j][f]=((*Z)[i][j][f]-partSum[j][f])/(*U)[j][j][f];
			  }
		  }
	  }

	return (double complex *)X;
}



double complex *cascade( int Nfreqs, int m1,double complex matrix1[m1][m1][Nfreqs], int m2,double complex matrix2[m2][m2][Nfreqs], int N){
	if(N>m1){
		fprintf(stderr,"Too many connections between matrices\n (Matrix 1 has size %d, but the connections are %d)\n",m1,N);
		exit(EXIT_FAILURE);
	}
	if(N>m2){
		fprintf(stderr,"Too many connections between matrices\n (Matrix 2 has size %d, but the connections are %d)\n",m2,N);
		exit(EXIT_FAILURE);
	}
	if(N==m1){
		fprintf(stderr,"Too many connections between matrices\n (Matrix 1 has size %d, same as the connections)\n",m1);
		exit(EXIT_FAILURE);
	}
	if(N==m2){
		fprintf(stderr,"Too many connections between matrices\n (Matrix 2 has size %d, same as the connections)\n",m2);
		exit(EXIT_FAILURE);
	}
	if(N<0){
		fprintf(stderr,"The number of interconnections must be possitive, not %d",N);
		exit(EXIT_FAILURE);
	}
	double complex SL11[m1-N][m1-N][Nfreqs];
	double complex SL12[m1-N][N][Nfreqs];
	double complex SL21[N][m1-N][Nfreqs];
	double complex SL22[N][N][Nfreqs];

	double complex SR11[N][N][Nfreqs];
	double complex SR12[N][m2-N][Nfreqs];
	double complex SR21[m2-N][N][Nfreqs];
	double complex SR22[m2-N][m2-N][Nfreqs];

	//SL11
	for(int i=0;i<m1-N;i++){
		for(int j=0;j<m1-N;j++){
			for(int f=0;f<Nfreqs;f++){
				SL11[i][j][f]=matrix1[i][j][f];
			}
		}
	}
	//SL12
	for(int i=0;i<m1-N;i++){
		for(int j=0;j<N;j++){
			for(int f=0;f<Nfreqs;f++){
				SL12[i][j][f]=matrix1[i][m1-N+j][f];
			}
		}
	}
	//SL21
	for(int i=0;i<N;i++){
		for(int j=0;j<m1-N;j++){
			for(int f=0;f<Nfreqs;f++){
				SL21[i][j][f]=matrix1[m1-N+i][j][f];
			}
		}
	}
	//SL22 && SR11
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			for(int f=0;f<Nfreqs;f++){
				SL22[i][j][f]=matrix1[m1-N+i][m1-N+j][f];
				SR11[i][j][f]=matrix2[i][j][f];
			}
		}
	}
	//SR12
	for(int i=0;i<N;i++){
		for(int j=0;j<m2-N;j++){
			for(int f=0;f<Nfreqs;f++){
				SR12[i][j][f]=matrix1[m2-N+i][j][f];
			}
		}
	}
	//SR21
	for(int i=0;i<m2-N;i++){
		for(int j=0;j<N;j++){
			for(int f=0;f<Nfreqs;f++){
				SR21[i][j][f]=matrix1[i][m2-N+j][f];
			}
		}
	}
	//SR22
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			for(int f=0;f<Nfreqs;f++){
				SR22[i][j][f]=matrix1[m2-N+i][m2-N+j][f];
			}
		}
	}

	double complex (*auxMatrix1)[N][N][Nfreqs]=(double complex (*)[N][N][Nfreqs])mmatrix(Nfreqs, N, N, SR11,N,N,SL22);
	memset(&auxMatrix1,0,N*N*Nfreqs*sizeof(double complex));
//	//double complex auxMatrix2[N][N][Nfreqs];
//	//memcpy(&auxMatrix2,auxMatrix1,N*N*Nfreqs*sizeof(double complex));
//
//	for(int i=0;i<N;i++){
//		for(int j=0;j<N;j++){
//			if(i==j){
//				for(int f=0;f<Nfreqs;f++){
//					auxMatrix2[i][j][f]=1-(*auxMatrix1)[i][j][f];
//				}
//			}else{
//				for(int f=0;f<Nfreqs;f++){
//					auxMatrix2[i][j][f]=-(*auxMatrix1)[i][j][f];
//				}
//			}
//		}
//	}
//	free(auxMatrix1);
	for(int i=0;i<N;i++){
		for(int f=0;f<Nfreqs;f++){
			(*auxMatrix1)[i][i][f]+=1;
		}
	}
	complex double **LU = LUdecomposition(N,Nfreqs,(double complex (*)[N][N][Nfreqs])auxMatrix1);
	free(auxMatrix1);
	complex double (*W)[N][N][Nfreqs] = (double complex (*)[N][N][Nfreqs])LUinverse(N,Nfreqs, LU);
	free((*LU)+1);
	free((*LU));
	free(LU);

	//double complex (*ST11)[N][N][Nfreqs]=(double complex (*)[N][N][Nfreqs])mmatrix(Nfreqs, N, N, SR11,N,N,SL22);

	//ST11 calculation
	//auxMatrix3 stored to compute ST21
	double complex (*auxMatrix3)[N][m1-N][Nfreqs]=(double complex (*)[N][m1-N][Nfreqs])mmatrix(Nfreqs, N, N, W,N,m1-N,SL21);
	double complex (*auxMatrix4)[N][m1-N][Nfreqs]=(double complex (*)[N][m1-N][Nfreqs])mmatrix(Nfreqs, N, N, SR11,N,m1-N,auxMatrix3);
	//free(auxMatrix3);
    double complex (*auxMatrix5)[m1-N][m1-N][Nfreqs]=(double complex (*)[m1-N][m1-N][Nfreqs])mmatrix(Nfreqs, m1-N,N,SL12,N,m1-N,auxMatrix4);
    free(auxMatrix4);

    double complex (*ST)[m1+m2-2*N][m1+m2-2*N][Nfreqs]=(double complex (*)[m1+m2-2*N][m1+m2-2*N][Nfreqs])malloc((m1+m2-2*N)*(m1+m2-2*N)*Nfreqs*sizeof(double complex));

    for(int i=0; i<m1-N;i++){
    	for(int j=0;j<m1-N;j++){
    		for(int f=0;f<Nfreqs;f++){
    			(*ST)[i][j][f]=(*auxMatrix5)[i][j][f]+SL11[i][j][f];
    		}
    	}
    }
    free(auxMatrix5);

    //ST12 calculation
    //auxMatrix6 will be reused in ST22 to save 1 matrix multiplication
    double complex (*auxMatrix6)[N][N][Nfreqs]=(double complex (*)[N][N][Nfreqs])mmatrix(Nfreqs, N, N, W,N,N,SL22);
    free(W);
    double complex (*auxMatrix7)[N][N][Nfreqs]=(double complex (*)[N][N][Nfreqs])mmatrix(Nfreqs, N,N,SR11,N,N,auxMatrix6);
    //plus Identity matrix
    for(int i=0;i<N;i++){
    	for(int f=0;f<Nfreqs;f++){
    		(*auxMatrix7)[i][i][f]+=1;
    	}
    }
    double complex (*auxMatrix8)[m1-N][N][Nfreqs]=(double complex (*)[m1-N][N][Nfreqs])mmatrix(Nfreqs, m1-N,N,SL12,N,N,auxMatrix7);
    free(auxMatrix7);
    double complex (*auxMatrix9)[m1-N][m2-N][Nfreqs]=(double complex (*)[m1-N][m2-N][Nfreqs])mmatrix(Nfreqs, m1-N,N,auxMatrix4,N,m2-N,SL21);
    free(auxMatrix8);

    for(int i=0; i<m1-N;i++){
    	for(int j=0;j<m1-N;j++){
    		for(int f=0;f<Nfreqs;f++){
    			(*ST)[i][N+j][f]=(*auxMatrix9)[i][j][f];
    		}
    	}
    }
    free(auxMatrix9);

    //ST21
    double complex (*auxMatrix10)[m2-N][m1-N][Nfreqs]=(double complex (*)[m2-N][m1-N][Nfreqs])mmatrix(Nfreqs,m2-N,N,SR21,N,m1-N,auxMatrix3);
    free(auxMatrix3);

    for(int i=0; i<m1-N;i++){
    	for(int j=0;j<m1-N;j++){
    		for(int f=0;f<Nfreqs;f++){
    			(*ST)[N+i][j][f]=(*auxMatrix10)[i][j][f];
    		}
    	}
    }
    free(auxMatrix10);

    //ST22
    double complex (*auxMatrix11)[N][m2-N][Nfreqs]=(double complex (*)[N][m2-N][Nfreqs])mmatrix(Nfreqs,N,N,auxMatrix6,N,m2-N,SR12);
    free(auxMatrix6);
    double complex (*auxMatrix12)[m2-N][m2-N][Nfreqs]=(double complex (*)[m2-N][m2-N][Nfreqs])mmatrix(Nfreqs,m2-N,N,SR21,N,m2-N,auxMatrix11);
    free(auxMatrix11);

    for(int i=0; i<m1-N;i++){
    	for(int j=0;j<m1-N;j++){
    		for(int f=0;f<Nfreqs;f++){
    			(*ST)[N+i][N+j][f]=(*auxMatrix12)[i][j][f]+SR22[i][j][f];
    		}
    	}
    }
    free(auxMatrix12);

    //double complex (*L)[m][m][Nfreqs]= (double complex (*)[m][m][Nfreqs])malloc(m*m*Nfreqs*sizeof(double complex));
    //memset(&L,0,m*m*Nfreqs*sizeof(double complex));

	return (double complex *)ST;
}

double complex * SfromCoupling(int Nfreqs, int N, double complex X[N][N][Nfreqs]){
	//X transpose
	double complex Xt[N][N][Nfreqs];
    for(int i=0; i<N;i++){
    	for(int j=0;j<N;j++){
    		for(int f=0;f<Nfreqs;f++){
    			Xt[i][j][f]=X[j][i][f];
    		}
    	}
    }
    //X * Xt +I or X * Xt -I
    double complex (*XXtP)[N][N][Nfreqs]=(double complex (*)[N][N][Nfreqs])mmatrix(Nfreqs, N,N,X,N,N,Xt);
    double complex (*XXtM)[N][N][Nfreqs]=(double complex (*)[N][N][Nfreqs])malloc(N*N*Nfreqs*sizeof(double complex));
    memcpy(XXtM,XXtP,N*N*Nfreqs*sizeof(double complex));
    for(int i=0; i<N;i++){
		for(int f=0;f<Nfreqs;f++){
			(*XXtP)[i][i][f]+=1;
			(*XXtM)[i][i][f]-=1;
		}
 	}
    //(X * Xt +I)^-1
    double complex** XXtPLU=LUdecomposition(N,Nfreqs,XXtP);
    double complex (*XXtPInverse)[N][N][Nfreqs]=(double complex (*)[N][N][Nfreqs])LUinverse(N,Nfreqs,XXtPLU);
    free(XXtP);
    free(*XXtPLU);
    free(*XXtPLU+1);
    free(XXtPLU);
    //S11=[X * Xt +I]^-1 * [X * Xt-I]
    double complex (*S11)[N][N][Nfreqs]=(double complex (*) [N][N][Nfreqs]) mmatrix(Nfreqs,N,N,XXtPInverse, N,N,XXtM);
    //S12=2*[X * Xt +I]^-1*X
    //the 2 is missing, incorporated at S
    double complex (*S12)[N][N][Nfreqs]=(double complex (*) [N][N][Nfreqs]) mmatrix(Nfreqs,N,N,XXtPInverse, N,N,X);
    free(XXtPInverse);
    free(XXtM);
    //S22=I-Xt*S12
    //the change of sign will be added at S (end of method)
    double complex (*S22)[N][N][Nfreqs]=(double complex (*) [N][N][Nfreqs]) mmatrix(Nfreqs,N,N,Xt,N,N,S12);
    //S11-I
    double complex (*ImS11)[N][N][Nfreqs]=(double complex (*) [N][N][Nfreqs])malloc(N*N*Nfreqs*sizeof(double complex));
    memcpy(ImS11,S11,N*N*Nfreqs*sizeof(double complex));
    for(int i=0;i<N;i++){
    	for(int f=0;f<Nfreqs;f++){
    		(*ImS11)[i][i][f]-=1;
    		(*S22)[i][i][f]-=1;
    	}
    }
    //S21=Xt[I-S11]=Xt*-ImS11=-Xt*ImS11, the minus at S
    double complex (*S21)[N][N][Nfreqs]=(double complex (*) [N][N][Nfreqs]) mmatrix(Nfreqs,N,N,Xt,N,N,ImS11);
    free(ImS11);

    //S is filled, the missing 2 and -s, are added.
    double complex (*S)[2*N][2*N][Nfreqs]=(double complex (*) [2*N][2*N][Nfreqs])malloc(4*N*N*Nfreqs*sizeof(double complex));
    for(int i=0;i<2*N;i++){
    	for(int j=0;j<2*N;j++){
			if(i<N && j<N){
				for(int f=0;f<Nfreqs;f++)
					(*S)[i][j][f]=(*S11)[i][j][f];
			}else if(i<N && j>=N){
				for(int f=0;f<Nfreqs;f++)
					(*S)[i][j][f]=2*(*S12)[i][j-N][f];
			}else if(i>=N && j<N){
				for(int f=0;f<Nfreqs;f++)
					(*S)[i][j][f]=-(*S21)[i-N][j][f];
			}else{
				for(int f=0;f<Nfreqs;f++)
					(*S)[i][j][f]=-(*S22)[i][j][f];
			}
    	}
    }
    free(S11);free(S12);free(S21);free(S22);
    return (double complex *)S;
}
