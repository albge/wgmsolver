#include "wgmsolver.h"

/*
Function to search for the first N modes in a rectangular waveguide.
It only returns cutoff frequencies and m,n indixes, do not specify if the
mode is TE or TM
Parameters
 - width wide side of the guide
 - height narrow side (NOTE: can be bigger than width)
 - N number of modes searched

Returns
 - Modes struct
*/

Modes rectangularmodes (double width, double height, int N)
{
int found =0;

int m_max=0;
int n_max=0;

Modes guide;
guide.cutFrequency[0]=0;
guide.firstcoor[0]=0;
guide.secondcoor[0]=0;
guide.type[0]=0;


while (found < N)
{
	int n;
	double possiblemodes [4][n_max+1];
	for (n=0;n<=n_max+1;n++)
	{
		int m=0;
		while (cutoff(width,height,m,n)<=guide.cutFrequency[found])
	 		m++;
		if (m!=0 && n!=0)
		{
		possiblemodes[1][n]=cutoff(width,height,m,n);
		possiblemodes[2][n]=m;
		possiblemodes[3][n]=n;
		possiblemodes[4][n]=2;
		}
		else
		{
		possiblemodes[1][n]=cutoff(width,height,m,n);
		possiblemodes[2][n]=m;
		possiblemodes[3][n]=n;
		possiblemodes[4][n]=1;
		}
	}/*for (p=0:n_max:p++)*/
	
	
	//The next loop works if there is no multiple modes with the same cutoff. 
	//A sort is needed
    int i;
    int next=0;
    for(i=0;i<n_max+1;i++)	
	{
		if (possiblemodes[1][i]<possiblemodes[1][next])
		{
			next=i;
			if (i>n_max)
				n_max++;
		}
	}
	if (possiblemodes[4][next] == 2){
		//degenerated modes, TE and TM are added
		guide.cutFrequency[found]=possiblemodes[1][next];
		guide.firstcoor[found]=possiblemodes[2][next];
		guide.secondcoor[found]=possiblemodes[3][next];
		guide.type[found]=0;
		
		guide.cutFrequency[found+1]=possiblemodes[1][next];
		guide.firstcoor[found+1]=possiblemodes[2][next];
		guide.secondcoor[found+1]=possiblemodes[3][next];
		guide.type[found+1]=1;
		
		found+=2;
	}else{
		//only TE
		guide.cutFrequency[found]=possiblemodes[1][next];
		guide.firstcoor[found]=possiblemodes[2][next];
		guide.secondcoor[found]=possiblemodes[3][next];
		guide.type[found]=0;
		found++;
	}
}/*while (found < N)*/

}
        
double cutoff(double width, double height, int m, int n)
{
	return 1.0/((sqrt(e*u*er*ur))*2.0*pi)*sqrt((m*pi)/width*(m*pi)/width+(n*pi)/height*(n*pi)/height);
}
