/*
 ============================================================================
 Name        : wgmsolver.c
 Author      : 
 Version     :
 Copyright   : 
 Description : Main
 ============================================================================
 */

#include "wgmsolver.h"

int main( int argc,char *argv[])
{
  if(argc<2){
	printf( "No arguments\n");
  	usage();
  	return EXIT_FAILURE;
  }
  
  int option;
/* process command line options */
  while( (option = getopt(argc, argv, "i:o:hv") ) != -1 )
  {
    switch( option )
    {   
      case 'h' : /* print usage and exit */
		usage();
		return EXIT_SUCCESS;

      case 'v' : /* print nec2c version */
		puts( VERSION );
		return EXIT_SUCCESS;

      default: /* print usage and exit */
		usage();
		return EXIT_FAILURE;
    } /* end of switch( option ) */
  } /* while( (option = getopt(argc, argv, "i:o:hv") ) != -1 ) */
}
