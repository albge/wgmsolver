/*
 ============================================================================
 Name        : wgmsolver.c
 Author      : Alberto González
 Version     :
 Copyright   : 
 Description : Main
 ============================================================================
 */

#include "wgmsolver.h"

parameters Params;

topology Topolog;

int main( int argc,char *argv[]){
  if(argc<2){
	printf( "No arguments\n");
  	usage();
  	return EXIT_FAILURE;
  }
  int option;
  char *fileInputName=NULL;
  char *fileOutputName=NULL;
  /* process command line options */
  while( (option = getopt(argc, argv, "hvf:o: ") ) != -1 )
  {
    switch( option )
    {   
      case 'h' : /* print usage and exit */
		usage();
		return EXIT_SUCCESS;

      case 'v' :
		puts( VERSION );
		return EXIT_SUCCESS;
      case 'f':
    	  fileInputName=optarg;
    	  printf ("option f with value '%s'\n", optarg);
    	  configRead(optarg);
    	  break;
      case 'o':
    	  fileOutputName=optarg;
    	  printf ("option o with value '%s'\n", optarg);
    	  break;
      default: /* print usage and exit */
		usage();
		return EXIT_FAILURE;
    } /* end of switch( option ) */
  } /* while( (option = getopt(argc, argv, "i:o:hv") ) != -1 ) */
}
