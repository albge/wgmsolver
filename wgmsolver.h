#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


#define		TRUE		1
#define		FALSE		0

#define 	pi		3.141592653589793//238462643383279502884197169399375
#define 	c		299792458


#define         u 		0.0000004*pi
#define         e               1.0/(u*c*c)
#define         er		1.0
#define         ur		1.0


#define		version 	"wgmsolver v0"

/*struct for each guide, to have the differents modes sorted according cutoff frequency
The type is 0 for TE, 1 for TM.
*/
typedef struct {
	float cutFrequency[100];
	int firstcoor[100];
	int secondcoor[100];
	int type[100];
} Modes;


// wgmsolver.c

int main(int argc, char *argv[]);

// misc.c

void 	usage(void);

//rectangular.c
Modes rectangularmodes (double width, double height, int N);
double cutoff(double width, double height, int m, int n);

