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
	double complex (*data)[Nports][Nports][Nfreqs]=Sparams;
	//max line length
	char buffer[1024];
	//max 9999999999 frequencies
	char freqs[10];
	//max 9999999999 ports
	char ports[10];
	FILE *fp;

	fp=fopen("filename", "w");
	if (fp == NULL) {
	  fprintf(stderr, "Can't open output file %s!\n",
	          filename);
	  exit(1);
	}

	sprintf(ports, "%d",Nports);
	sprintf(freqs, "%d",Nfreqs);

	//first line, ! n-port S parameters, f frequency points.
	strcpy(buffer,"! ");
	strcat(buffer,ports);
	strcat(buffer,"-port S parameters, " );
	strcat(buffer, freqs);
	strcat(buffer, " frequency points.\n");

	fprintf(fp,"%s", buffer);
	memset(buffer, 0, sizeof(buffer));

	fprintf(fp, "[Version] 2.0\n");
	fprintf(fp, "# GHz S RI R 50\n");

	strcpy(buffer,"[Number of Ports]");
	strcat(buffer, ports);
	strcat(buffer, "\n");

	fprintf(fp,"%s", buffer);

	memset(buffer, 0, sizeof(buffer));

	strcpy(buffer,"[Number of Frequencies] ");
	strcat(buffer, freqs);
	strcat(buffer, "\n");

	fprintf(fp,"%s", buffer);

	memset(buffer, 0, sizeof(buffer));

	fprintf(fp,"[Matrix Format] Full\n");

	char frequency[12];
	for(int f=0;f<Nfreqs;f++){
		sprintf(frequency,"%.10e ",*(freqValues+f));
		for(int i=0;i<Nports;i++){
			sprintf(buffer,"! %d row",i);
			fprintf(fp,"%s", buffer);
			memset(buffer, 0, sizeof(buffer));

			char notfirst=FALSE;
			int jlim=Nports/4;
			int j=0;
			for(int j4=0;j<jlim;j4++){
				if(notfirst)
				memset(frequency, 0, sizeof(frequency));
				//precission 8
				sprintf(buffer,"%s %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g \n",frequency,
						creal(*data[i][j][f]) , cimag(*data[i][j][f]),creal(*data[i][j+1][f]) , cimag(*data[i][j+1][f]),
						creal(*data[i][j+2][f]) , cimag(*data[i][j+2][f]),creal(*data[i][j+3][f]) , cimag(*data[i][j+3][f]));
				fprintf(fp,"%s", buffer);
				memset(buffer, 0, sizeof(buffer));
				j+=4;
			}
			jlim=Nports%4;
			switch(jlim){
			case 3:
				sprintf(buffer,"%s %.8g %.8g %.8g %.8g %.8g %.8g \n",frequency,
				creal(*data[i][j][f]) , cimag(*data[i][j][f]),creal(*data[i][j+1][f]) , cimag(*data[i][j+1][f]),
				creal(*data[i][j+2][f]) , cimag(*data[i][j+2][f]));
				fprintf(fp,"%s", buffer);
				break;
			case 2:
				sprintf(buffer,"%s %.8g %.8g %.8g %.8g \n",frequency,
				creal(*data[i][j][f]) , cimag(*data[i][j][f]),creal(*data[i][j+1][f]) , cimag(*data[i][j+1][f]));
				fprintf(fp,"%s", buffer);
				break;
			case 1:
				sprintf(buffer,"%s %.8g %.8g \n",frequency,
				creal(*data[i][j][f]) , cimag(*data[i][j][f]));
				fprintf(fp,"%s", buffer);
				break;
			default: break;
			}
			memset(frequency, 0, sizeof(frequency));
		}
	}

	fclose(fp);

	return TRUE;
}
