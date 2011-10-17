#include "wgmsolver.h"

int main( int argc,char *argv[])
{
  if(argc<2){
	printf( "No arguments\n");
  	usage();
  	exit(-1);
  }
  
  int option;
/* process command line options */
  while( (option = getopt(argc, argv, "i:o:hv") ) != -1 )
  {
    switch( option )
    {   
      case 'h' : /* print usage and exit */
		usage();
		exit(0);

      case 'v' : /* print nec2c version */
		puts( version );
		exit(0);

      default: /* print usage and exit */
		usage();
		exit(-1);
    } /* end of switch( option ) */
  } /* while( (option = getopt(argc, argv, "i:o:hv") ) != -1 ) */
  
  
  
}
