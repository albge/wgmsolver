#include "wgmsolver.h"
/*
 * Touchstone.c
 *
 *  Created on: Nov 27, 2011
 *      Author: alberto
 */

/**
 * Function writeTouchstone
 * Writes a touchstone file with a full matrix, Real and Imaginary part, for Nports, Nfreqs, in the file specified in filename.
 */
char writeTouchstone(char* filename, double complex* Sparams,double *freqValues, int Nports, int Nfreqs){
	double complex (*data)[Nports][Nports][Nfreqs]=(double complex (*) [Nports][Nports][Nfreqs])Sparams;

	FILE *fp;

	fp=fopen(filename, "w");
	if (fp == NULL) {
	  fprintf(stderr, "Can't open output file %s!\n",
	          filename);
	  exit(1);
	}

	//first line, ! n-port S parameters, f frequency points.
	fprintf(fp,"! %d-port S parameters, %d frequency points.\n", Nports, Nfreqs);

	fprintf(fp, "[Version] 2.0\n");
	fprintf(fp, "# GHz S RI R 50\n");
	fprintf(fp,"[Number of Ports] %d \n", Nports);
	fprintf(fp,"[Number of Frequencies] %d \n",Nfreqs);
	fprintf(fp,"[Matrix Format] Full\n");

	char frequency[12];
	for(int f=0;f<Nfreqs;f++){
		fprintf(fp,"%.10e ",*(freqValues+f));
		for(int i=0;i<Nports;i++){
			fprintf(fp,"! %d row",i);
			char notfirst=FALSE;
			int jlim=Nports/4;
			int j=0;
			for(int j4=0;j<jlim;j4++){
				if(notfirst)
					memset(frequency, 0, sizeof(frequency));
				//precission 8
				fprintf(fp,"%s %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g \n",frequency,
						creal(*data[i][j][f]) , cimag(*data[i][j][f]),creal(*data[i][j+1][f]) , cimag(*data[i][j+1][f]),
						creal(*data[i][j+2][f]) , cimag(*data[i][j+2][f]),creal(*data[i][j+3][f]) , cimag(*data[i][j+3][f]));
				//fprintf(fp,"%s", buffer);
				//memset(buffer, 0, sizeof(buffer));
				j+=4;
			}
			jlim=Nports%4;
			switch(jlim){
			case 3:
				fprintf(fp,"%s %.8g %.8g %.8g %.8g %.8g %.8g \n",frequency,
				creal(*data[i][j][f]) , cimag(*data[i][j][f]),creal(*data[i][j+1][f]) , cimag(*data[i][j+1][f]),
				creal(*data[i][j+2][f]) , cimag(*data[i][j+2][f]));
				//fprintf(fp,"%s", buffer);
				break;
			case 2:
				fprintf(fp,"%s %.8g %.8g %.8g %.8g \n",frequency,
				creal(*data[i][j][f]) , cimag(*data[i][j][f]),creal(*data[i][j+1][f]) , cimag(*data[i][j+1][f]));
				//fprintf(fp,"%s", buffer);
				break;
			case 1:
				fprintf(fp,"%s %.8g %.8g \n",frequency,
				creal(*data[i][j][f]) , cimag(*data[i][j][f]));
				//fprintf(fp,"%s", buffer);
				break;
			default: break;
			}
			memset(frequency, 0, sizeof(frequency));
		}
	}

	fclose(fp);

	return TRUE;
}
