#include "wgmsolver.h"
/*
 *
 * readfile.c
 *
 *  Created on: Jun 29, 2012
 *      Author: alberto
 */

#define COMMENTSYMBOL '#'
#define RESERVEDWORDS "Units ","Frequencies ", "FStep ", "MaxFreq "
#define NWORDS         4

#define     NUNITS          8
#define 	UNITS		    "m","dm","cm","mm","yd","ft","in","th"
#define 	UNITS_VALUE    1,0.1,0.01,0.001,0.9144,0.3048,0.0254,0.0000254



void configRead(char* filename)
{
        FILE*fp ;

     	fp=fopen(filename, "r");
     	if (fp == NULL) {
     	  fprintf(stderr, "Can't open output file %s!\n",
     	          filename);
     	  exit(1);
     	}

        char *line = NULL;
        size_t len = 0;
        ssize_t read;

	char * ready;
        while ((read = getline(&line, &len, fp)) != -1) {
            ready = removeComments(line);
            //printf("Retrieved line of length %zu :\n", read);
            if (ready!=NULL){
	    	printf("%s \n", ready);
	    }
	    free(ready);
        }
}


char* removeComments(char* line){
	//Remove comments
	char comment= COMMENTSYMBOL;
	size_t commentStart = strcspn(line, &comment);

	if(commentStart==0){
		return NULL;
	}

	char * decommented = (char *)malloc(commentStart*sizeof(char));
	strncpy(decommented, line, commentStart);

	return decommented;
}

void printReservedWords(FILE * stream){
	char *words[] = {RESERVEDWORDS};
	fprintf(stream, "Reserved words are: \n");
	for(int i=0; i<NWORDS; i++){
		fprintf(stream, "%s", *(words+i));
	}
	fprintf(stream,"\n");
	return;
}

int charCount(char *string, char character){
	int found=0;
	for(int i=0; string[i]!='\0';i++){
		if(string[i]==character){
			found++;
		}
	}
	return found;
}

int loadUnits(char *line){
	char * unit = "Units";
	char **ptr=NULL;

	char * config = strdup(line);
	char * tokenized=strtok_r(config," []",ptr);

	if(tokenized==NULL){
	return -1;
	}
	if(strcmp(tokenized, unit)!=0){
	return 1;
	}
	tokenized = strtok_r(NULL," []",ptr);

	double unit_values [NUNITS] = {UNITS_VALUE};
	char *units[]={UNITS};
    extern parameters Params;
	for(int i=0; i<NUNITS; i++){
		if(strcasecmp(*units[i],tokenized)==0){
			Params.units=unit_values[i];
		}
	}
	return 0;
}
