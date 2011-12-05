/*
 * circular.c
 *
 *  Created on: 01/11/2011
 *      Author:
 */

#include "wgmsolver.h"



//array zeros Bessel (TMnm) (n=index, m=index+1)
double pnm[10][10]={{2.40483,5.52008,8.65373,11.79153,14.93092,18.07106,21.21164,24.35247,27.49348,30.63461},
					{3.83171,7.01559,10.17347,13.32369,16.47063,19.61586,22.76008,25.90367,29.04683,32.18968},
					{5.13562,8.41724,11.61984,14.79595,17.95982,21.11700,24.27011,27.42057,30.56920,33.71652},
					{6.38016,9.76102,13.01520,16.22347,19.40942,22.58273,25.74817,28.90835,32.06485,35.21867},
					{7.58834,11.06471,14.37254,17.61597,20.82693,24.01902,27.19909,30.37101,33.53714,36.69900},
					{8.77148,12.33860,15.70017,18.98013,22.21780,25.43034,28.62662,31.81172,34.98878,38.15987},
					{9.93611,13.58929,17.00382,20.32079,23.58608,26.82015,30.03372,33.23304,36.42202,39.60324},
					{11.08637,14.82127,18.28758,21.64154,24.93493,28.19119,31.42279,34.63709,37.83872,41.03077},
					{12.22509,16.03777,19.55454,22.94517,26.26681,29.54566,32.79580,36.02562,39.24045,42.44389},
					{13.35430,17.2412, 20.807, 24.2339, 27.5837, 30.8854, 34.1544, 37.4001, 40.6286, 43.8438}};


//array zeros BesselD (TEnm) (n=index, m=index+1)
double  qnm[10][10]={{3.83171,7.01559,10.17347,13.32370,16.47063,19.61586,22.76008,25.90367,29.04683,32.18968},
            		{1.84118,5.33144,8.53632,11.70600,14.86359,18.01553,21.16437,24.31133,27.45705,30.60192},
            		{3.05424,6.70613,9.96947,13.17037,16.34752,19.51291,22.67158,25.82604,28.97767,32.12733},
            		{4.20119,8.01524,11.34592,14.58585,17.78875,20.97248,24.14490,27.31006,30.47027,33.62695},
            		{5.31755,9.28240,12.68191,15.96411,19.19603,22.40103,25.58976,28.76784,31.93854,35.10392},
            		{6.41562,10.51986,13.98719,17.31284,20.57551,23.80358,27.01031,30.20285,33.38544,36.56078},
            		{7.50127,11.73494,15.26818,18.63744,21.93172,25.18393,28.40978,31.61788,34.81339,37.99964},
            		{8.57784,12.93239,16.52937,19.94185,23.26805,26.54503,29.79075,33.01518,36.22438,39.42227},
            		{9.64742,14.11552,17.77401,21.22906,24.58720,27.88927,31.15533,34.39663,37.62008,40.83018},
            		{}};

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
double radius= sect->radius;

//number of modes found
int Nmodes=0;

/*
* ptr array of modes
* */
mode *ptr= (mode *)malloc(N*sizeof(mode));

int n=0;
int m=1;
while(Nmodes<N){

		double q = qnm[n][m-1];
		double p = pnm[n][m-1];

		mode currentModeTE;
		currentModeTE.cutFrequency=cutoff(radius,0,n,m);
		currentModeTE.firstcoor=n;
		currentModeTE.secondcoor=m;
		currentModeTE.type=0;

		mode currentModeTM;
		currentModeTM.cutFrequency=cutoff(radius,0,n,m);
		currentModeTM.firstcoor=n;
		currentModeTM.secondcoor=m;
		currentModeTM.type=0;

		*(ptr+Nmodes)=currentModeTE;
		*(ptr+Nmodes+1)=currentModeTM;

		Nmodes+=2;
		n++;
		m++;

}

//sort cutFrequencies
qsort (ptr, N, sizeof(mode), compareF);

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

	double radius= sect->radius;

	//number of modes found
	int Nmodes=0;

	/*
	* ptr array of modes
	* */
	mode *ptr= (mode *)malloc(100*sizeof(mode));

	int n=0;
	int m=1;
	for(int n=0; n<=9; ++n){

		int Nmodes0=Nmodes-1;

		while(Nmodes0<Nmodes){

			Nmodes0=Nmodes;

			double q = qnm[n][m-1];
			double p = pnm[n][m-1];

			mode currentModeTE;
			currentModeTE.cutFrequency=cutoff(radius,0,n,m);
			currentModeTE.firstcoor=n;
			currentModeTE.secondcoor=m;
			currentModeTE.type=0;

			mode currentModeTM;
			currentModeTM.cutFrequency=cutoff(radius,0,n,m);
			currentModeTM.firstcoor=n;
			currentModeTM.secondcoor=m;
			currentModeTM.type=0;

			if(currentModeTE.cutFrequency<cutFreq){
				*(ptr+Nmodes)=currentModeTE;
				Nmodes++;
			}
			if(currentModeTM.cutFrequency<cutFreq){
				*(ptr+Nmodes)=currentModeTM;
				Nmodes++;
			}

			m++;

		}

	}

	//sort cutFrequencies
	qsort (ptr, 100, sizeof(mode), compareF);

	sect->modes=ptr;
	sect->Nmodes=Nmodes;

	//free the pointers to pointers
	free(ptr);
	return sect;
	}




double circularcircular(section *sect1,section *sect2){

	//Radius
	double R1= sect1->radius;
	double R2= sect1->radius;

	//Polarization x (0) or y (1)
	int pol1;
	int pol2;




	double X[sect1->Nmodes][sect2->Nmodes];

	for(int i=0;i<sect1->Nmodes;i++){
		for(int j=0;j<sect2->Nmodes;j++){

			//double kx1=(*sect1->modes+i)->firstcoor*pi/sect1->width;
			//double kx2=(*sect2->modes+j)->firstcoor*pi/sect2->width;
			//double ky1=(*sect1->modes+i)->secondcoor*pi/sect1->height;
			//double ky2=(*sect2->modes+j)->secondcoor*pi/sect2->height;

			//types (TE or TM)
			int type1=(*sect1->modes+i)->type;
			int type2=(*sect2->modes+i)->type;

			//indexes
			int n1=(*sect1->modes+i)->firstcoor;
			int m1=(*sect1->modes+i)->secondcoor;
			int n2=(*sect2->modes+i)->firstcoor;
			int m2=(*sect2->modes+i)->secondcoor;

			//k
			int k1;
			int k2;

			//Cross matrix

			if(n1==n2){ //same 1 coordinate

				if(type1==type2){ //same type

					if(pol1==pol2){ //same polarization

						if(type1==0){ //TE-TE
							k1=qnm[n1][m1-1]/R1;
							k2=qnm[n2][m2-1]/R2;
							X(i,j)=pi/(e*er)*2*R1*k1*jn(n1,k1*R1)*Djn(n1,k2*R1)/(k1^2-k2^2);

						}else{ //TM-TM
							k1=pnm[n1][m1-1]/R1;
							k2=pnm[n2][m2-1]/R2;
							X(i,j)=-pi/(e*er)*2*R1*k2*Djn(n1,k1*R1)*jn(n1,k2*R1)/(k1^2-k2^2);
						}

					}else{ //different polarization
						X(i,j)=0;
					}

				}else{ //different type

					if(type1==0){ //TE-TM

						if(pol1!=pol2){ //different polarization
							k1=qnm[n1][m1-1]/R1;
							k2=pnm[n2][m2-1]/R2;
							X(i,j)=-pi/(e*er)*2*n1*jn(n1,k1*R1)*jn(n1,k2*R1)/(k1*k2);

							if(pol1==1){ //y-x
								X(i,j)=-X(i,j);
							}

						}else{ //same polarization
							X(i,j)=0;
						}

					}else{ //TM-TE
						X(i,j)=0;
					}
				}

			}else{
				X(i,j)=0; //different first coordinate

			}
		}
	}

}






/*
Calculate the cutoff frequency of TEnm or TMnm
  Parameters
- radius of the guide
- type (0 TE, 1 TM)
- n,m
 */

double cutoff(double radius, int type, int n, int m)
{
if(type!=0){
	return 1.0/((sqrt(e*u*er*ur))*2.0*pi)*pnm[n][m-1]/radius; //TMnm
}
return 1.0/((sqrt(e*u*er*ur))*2.0*pi)*qnm[n][m-1]/radius; //TEnm
}


/*
 * Function that compares two cutoff frequencies
 */
//double compareF (const void *f1, const void *f2)

double compareF (const mode *m1, const mode *m2)
{
  return ( *(mode*)m1.cutFrequency - *(mode*)m2.cutFrequency );
}


/*
 * Function that calculates the first derivative of BesselJN
 */

double Djn (int n, double r)
{
	return (jn(n-1,r)-jn(n+1,r))/2;
}
