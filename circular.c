/*
 * circular.c
 *
 *  Created on: 03/12/2011
 *      Author:
 */

#include "wgmsolver.h"


//array zeros Bessel (TMnm) & zeros BesselD (TEnm) sorted
double pq[190]={1.8412,2.4048,3.0542,3.8317,3.8317,4.2012,5.1356,5.3175,5.3314,5.5201,
				6.3802,6.4156,6.7061,7.0156,7.0156,7.5013,7.5883,8.0152,8.4172,8.5363,
				8.5778,8.6537,8.7715,9.2824,9.6474,9.761,9.9361,9.9695,10.173,10.173,
				10.52,11.065,11.086,11.346,11.62,11.706,11.735,11.792,12.225,12.339,
				12.682,12.932,13.015,13.17,13.324,13.324,13.354,13.589,13.987,14.116,
				14.373,14.586,14.796,14.821,14.864,14.931,15.268,15.7,15.964,16.038,
				16.223,16.348,16.471,16.471,16.529,17.004,17.241,17.313,17.616,17.774,
				17.789,17.96,18.016,18.071,18.288,18.637,18.98,19.196,19.409,19.513,
				19.555,19.616,19.616,19.942,20.321,20.576,20.807,20.827,20.972,21.117,
				21.164,21.212,21.229,21.642,21.932,22.218,22.401,22.583,22.672,22.76,
				22.76,22.945,23.268,23.586,23.804,24.019,24.145,24.234,24.27,24.311,
				24.352,24.587,24.935,25.184,25.43,25.59,25.748,25.826,25.904,25.904,
				26.267,26.545,26.82,27.01,27.199,27.31,27.421,27.457,27.493,27.584,
				27.889,28.191,28.41,28.627,28.768,28.908,28.978,29.047,29.047,29.546,
				29.791,30.034,30.203,30.371,30.47,30.569,30.602,30.635,30.885,31.155,
				31.423,31.618,31.812,31.939,32.065,32.127,32.19,32.19,32.796,33.015,
				33.233,33.385,33.537,33.627,33.717,34.154,34.397,34.637,34.813,34.989,
				35.104,35.219,36.026,36.224,36.422,36.561,36.699,37.4,37.62,37.839,
				38,38.16,39.24,39.422,39.603,40.629,40.83,41.031,42.444,43.844};

//array index 1 (n) sorted
int index1[190]={1,0,2,1,0,3,2,4,1,0,3,5,2,1,0,6,4,3,2,1,
				 7,0,5,4,8,3,6,2,1,0,5,4,7,3,2,1,6,0,8,5,
				 4,7,3,2,1,0,9,6,5,8,4,3,2,7,1,0,6,5,4,8,
				 3,2,1,0,7,6,9,5,4,8,3,2,1,0,7,6,5,4,3,2,
				 8,1,0,7,6,5,9,4,3,2,1,0,8,7,6,5,4,3,2,1,
				 0,8,7,6,5,4,3,9,2,1,0,8,7,6,5,4,3,2,1,0,
				 8,7,6,5,4,3,2,1,0,9,8,7,6,5,4,3,2,1,0,8,
				 7,6,5,4,3,2,1,0,9,8,7,6,5,4,3,2,1,0,8,7,
				 6,5,4,3,2,9,8,7,6,5,4,3,8,7,6,5,4,9,8,7,
				 6,5,8,7,6,9,8,7,8,9};

//array index 2 (m) sorted
int index2[190]={1,1,1,1,1,1,1,1,2,2,1,1,2,2,2,1,1,2,2,3,
				 1,3,1,2,1,2,1,3,3,3,2,2,1,3,3,4,2,4,1,2,
				 3,2,3,4,4,4,1,2,3,2,3,4,4,2,5,5,3,3,4,2,
				 4,5,5,5,3,3,2,4,4,3,5,5,6,6,3,4,4,5,5,6,
				 3,6,6,4,4,5,3,5,6,6,7,7,4,4,5,5,6,6,7,7,
				 7,4,5,5,6,6,7,4,7,8,8,5,5,6,6,7,7,8,8,8,
				 5,6,6,7,7,8,8,9,9,5,6,6,7,7,8,8,9,9,9,6,
				 7,7,8,8,9,9,10,10,6,7,7,8,8,9,9,10,10,10,
				 7,8,8,9,9,10,10,7,8,8,9,9,10,10,8,9,9,10,
				 10,8,9,9,10,10,9,10,10,9,10,10,10,10};

//array type of mode (TE 0 or TM 1) sorted
int type[190]={0,1,0,1,0,0,1,0,0,1,1,0,0,1,0,0,1,0,1,0,
			   0,1,1,0,0,1,1,0,1,0,0,1,1,0,1,0,0,1,1,1,
			   0,0,1,0,1,0,1,1,0,0,1,0,1,1,0,1,0,1,0,1,
			   1,0,1,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,1,0,
			   1,1,0,0,1,0,1,1,0,1,0,1,0,1,0,1,0,1,0,1,
			   0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0,1,0,1,0,
			   1,0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0,1,0,1,
			   0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0,1,0,1,0,
			   1,0,1,0,1,1,0,1,0,1,0,1,1,0,1,0,1,1,0,1,
			   0,1,1,0,1,1,0,1,1,1};




/*
Function to search for the first N modes in a circular waveguide.

Parameters
- radius of the guide
- N number of modes searched

Returns
- Modes struct
*/

section *circularNum(section *sect, int N){

	//copy the dimensions
	double radius= sect->width;

	//number of modes found
	int ix=0;

	/*
	 * ptr array of modes
	 * */
	mode *ptr= (mode *)malloc(N*sizeof(mode));

	while(ix<N){

		mode currentMode;
		currentMode.cutFrequency=1.0/((sqrt(e*u*er*ur))*2.0*pi)*pq[ix]/radius;
		currentMode.firstcoor=index1[ix];
		currentMode.secondcoor=index2[ix];
		currentMode.type=type[ix];

		*(ptr+ix)=currentMode;

		ix++;
	}

	sect->modes=ptr;
	sect->Nmodes=N;

	//free the pointers to pointers
	free(ptr);
	return sect;
}




/*
Function to search all the modes below the cut-off frequency.

Parameters
- guide section
- cut-off frequency

Returns
- Modes struct
*/

section *circularFreq(section *sect, double cutFreq){

	double radius= sect->width;

	//number of modes found
	int ix=0;

	/*
	* ptr array of modes
	* */
	mode *ptr= (mode *)malloc(190*sizeof(mode));

	//initial frequency
	double freq_M=1.0/((sqrt(e*u*er*ur))*2.0*pi)*pq[ix]/radius;

	while(freq_M<cutFreq){

		mode currentMode;
		currentMode.cutFrequency=freq_M;
		currentMode.firstcoor=index1[ix];
		currentMode.secondcoor=index2[ix];
		currentMode.type=type[ix];

		*(ptr+ix)=currentMode;

		ix++;

		freq_M=1.0/((sqrt(e*u*er*ur))*2.0*pi)*pq[ix]/radius;

	}

	sect->modes=ptr;
	sect->Nmodes=ix;

	//free the pointers to pointers
	free(ptr);
	return sect;
	}




/*
Function to compute the junction between two circular waveguides

Parameters
- section 1, 2

Returns
- X matrix
*/

complex double *circularcircular(section *sect1,section *sect2){

	//Radius
	double R1= sect1->width;
	double R2= sect2->width;

	//Polarization x (0) or y (1)
	int pol1=0;
	int pol2=0;


	complex double X[sect1->Nmodes][sect2->Nmodes];

	for(int i=0;i<sect1->Nmodes;i++){
		for(int j=0;j<sect2->Nmodes;j++){

			//double kx1=(*sect1->modes+i)->firstcoor*pi/sect1->width;
			//double kx2=(*sect2->modes+j)->firstcoor*pi/sect2->width;
			//double ky1=(*sect1->modes+i)->secondcoor*pi/sect1->height;
			//double ky2=(*sect2->modes+j)->secondcoor*pi/sect2->height;

			//types (TE or TM)
			int type1=(sect1->modes+i)->type;
			int type2=(sect2->modes+j)->type;

			//indexes
			int n1=(sect1->modes+i)->firstcoor;
			int m1=(sect1->modes+i)->secondcoor;
			int n2=(sect2->modes+j)->firstcoor;
			int m2=(sect2->modes+j)->secondcoor;

			//k
			int k1=pq[i]/R1;
			int k2=pq[j]/R2;

			//Cross matrix

			if(n1==n2){ //same 1 coordinate

				if(type1==type2){ //same type

					if(pol1==pol2){ //same polarization

						if(type1==0){ //TE-TE
							X[i][j]=pi/(e*er)*2*R1*k1*jn(n1,k1*R1)*Djn(n1,k2*R1)/(k1^2-k2^2);

						}else{ //TM-TM
							X[i][j]=-pi/(e*er)*2*R1*k2*Djn(n1,k1*R1)*jn(n1,k2*R1)/(k1^2-k2^2);
						}
					}else{ //different polarization
						X[i][j]=0;
					}

				}else{ //different type

					if(type1==0){ //TE-TM

						if(pol1!=pol2){ //different polarization
							X[i][j]=-pi/(e*er)*2*n1*jn(n1,k1*R1)*jn(n1,k2*R1)/(k1*k2);

							if(pol1==1){ //y-x
								X[i][j]=-X[i][j];
							}

						}else{ //same polarization
							X[i][j]=0;
						}

					}else{ //TM-TE
						X[i][j]=0;
					}
				}

			}else{
				X[i][j]=0; //different first coordinate

			}
		}
	}

	return X;

}




/*
 * Function that calculates the first derivative of BesselJN
 *
 * Parameters
 * n order
 * r argument
 */
double Djn (int n, double r){
	return (jn(n-1,r)-jn(n+1,r))/2;
}
