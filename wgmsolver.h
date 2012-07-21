#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
//#include <ctype.h>
#include <unistd.h>
#include <getopt.h>
//#include <features.h>

#define		TRUE		1
#define		FALSE		0

#define 	pi		3.141592653589793//238462643383279502884197169399375
#define 	c		299792458


#define         u 		0.0000004*pi
#define         e       1.0/(u*c*c)
#define         er		1.0
#define         ur		1.0

#define			I		_Complex_I

#define		VERSION 	"wgmsolver v0"

typedef struct {
	double x;
	double y;
} point;

typedef struct {
	double cutFrequency;
	int firstcoor;
	int secondcoor;
	int type;
} mode;

typedef struct {
	char* type;
	double width; //radius for circular
	double height; //angle for circular
	double length;
	point center;
	int Nmodes;
	mode *modes;
} section;

typedef struct{
<<<<<<< HEAD
	section* source;
	int connexions;
	section* connected;
} node;

=======
	int Nsections;	
	section *sections;
	
} topology;

typedef struct{
	double height;
	double width;
	char* type; //DS, Eplane, Hplane, empty --> ignore optimizations
	section rightSection;
	section leftSection;
} interserction;
>>>>>>> be0a697cae05dc58a67eb1ab904671ce7f2d2314
// wgmsolver.c

int main(int argc, char *argv[]);
int getopt(int argc, char * const argv[], const char *optstring);

// misc.c

void usage(void);

//rectangular.c
section *rectangularNum(section *sect, int N);
section *rectangularFreq(section *sect, double fc);
double complex *rectangularrectangular(section *sect1, section *sect2, double *freqs,double complex *Sparams);
//Modes rectangularmodes (double width, double height, int N);

//circular.c
double Djn (int n, double r);
