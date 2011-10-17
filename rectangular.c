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
int m=0
int n=0;
int found =0;

int m_max=0;
int n_max=0;

Modes guide;
guide.cutFrequency[0]=0;
guide.firstcoor[0]=0;
guide.secondcoo[0]=0;
guide.type[0]=0;


while (found < N)
{
	int p;
	double possiblemodes [4][n_max+1];
	for (p=0,p<=n_max+1:p++)
	{
		int k=0;
		while (cutoff(width,height,m,n)<=guide.cutFrequency[found])
	 		k++;
		if (k!=0 && p!=0)
		{
		possiblemodes[1][p]=cutoff(width,height,m,n);
		possiblemodes[2][p]=k;
		possiblemodes[3][p]=p;
		possiblemodes[4][p]=2;
		}
		else
		{
		possiblemodes[1][p]=cutoff(width,height,m,n);
		possiblemodes[2][p]=k;
		possiblemodes[3][p]=p;
		possiblemodes[4][p]=1;
		}
	}/*for (p=0:n_max:p++)*/


	
	guide.cutFrequency[0]=0;
	guide.firstcoor[0]=0;
	guide.secondcoo[0]=0;
	guide.type[0]=0;
}/*while (found < N)*/

}
    
            %all the modes get sortes according their cutoff frequencies
            m=sortrows(thisIteration')';
    
            %If we have multiple modes with the same cutoff, we select them all.
            for i=1:length(m(1,:))
                if m(1,1) == m(1,i)
                    modes = [modes m(:,i)];
                           
                    found =found+1;
                end
            end
        end
        %So we skip the 0 at the beginning.
        m=modes(:,2:N+1);
        end
        
double cutoff(double width, double height, int m, int n)
{
	1/((sqrt(e*u*er*ur))*2*pi)*sqrt(((m*pi)/width)^2+((n*pi)/height)^2);
}
