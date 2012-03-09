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

int main( int argc,char *argv[]){
  if(argc<2){
	printf( "No arguments\n");
  	usage();
  	return EXIT_FAILURE;
  }
	point p={1,1};
	section s={"rectangular",0.48,0.24,2,{2,2}};
  int option;
/* process command line options */
  while( (option = getopt(argc, argv, "i:o:hv") ) != -1 )
  {
    switch( option )
    {   
      case 'h' : /* print usage and exit */
		usage();
		return EXIT_SUCCESS;

      case 'v' :
		//puts( VERSION );
/*
		typedef struct {
			double cutFrequency;
			int firstcoor;
			int secondcoor;
			int type;
		} mode;
*/
		s= *rectangularFreq(&s, 1000000000);
		printf("hello %d \n",(s.Nmodes));
		return EXIT_SUCCESS;
      default: /* print usage and exit */
		usage();
		return EXIT_FAILURE;
    } /* end of switch( option ) */
  } /* while( (option = getopt(argc, argv, "i:o:hv") ) != -1 ) */
}
