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

section *rectangularNum(section *sect, int N){
	//copy the dimensions
	double width= sect->width;
	double height= sect->height;

	//calculate the ratio between dimensions
	double alpha = width/height;

	/* 	max m (TEmn notation), using the elipse formula to find a maximum
	 *	because the cutFreq formula can be seen as an elipse formula.
	 *	Understanding this, the number of modes is related to the area covered
	 *	by the curve with cutoff frequency constant.
	 *	Out of that area/modes relatio, the max m is found.*/
	int max_m =(int) (1+alpha+sqrt((1+alpha)*(1+alpha)+32*N*pi*alpha))/(4*pi*(1+alpha));
	//max n
	//double max_n = max_m/alpha;

	double cutFreq = cutoff(width, height, max_m, 0);

	//double max_m= floor(2*width*sqrt(e*er*u*ur*cutFreq));
	/*Once the max fc has been found, the maximum n is inmediate*/
	int max_n=(int)floor(2*height*sqrt(e*er*u*ur*cutFreq));//floor(max_m*height/width);
	//array of indexes for the num of elements in each row (max_n rows)
	int nElems[max_n];
	memset(nElems, 0, max_n*sizeof(int)); //init to zero
	int Nmodes=0; //zero modes found so far

	//sect.modes=(mode **)malloc(max_n*sizeof(*mode));
	/*
	 * A pointer to pointers to the modes is allocated
	 * a+1 represents the row with n=1, a+7 n=7
	 * *(a+1)+2 points at the mode with n=1, m=2, ...&al
	 * */

	mode **a= (mode **)malloc(max_n*sizeof(mode*));
	for(int n=0; n<=max_n; ++n){
		//the maximum m FOR THIS ROW! (local maximum)
		int m_max=(int)width*sqrt((4*e*u*er*ur*cutFreq*cutFreq-(n*n/(height*height))));
		nElems[n]=m_max;
		Nmodes +=m_max;
		//the TM modes are taken into account
		if (n!=0){
			Nmodes += (m_max-1>0)? m_max-1:0;
			nElems[n]+=(m_max-1>0)? m_max-1:0;
		}
		//the row is allocated in memory
		*(a+n)=(mode *)malloc(m_max*sizeof(mode));
		//the row is filled with the modes
		for(int m=0;m<=m_max;m++){
			mode currentMode;
			//TE
			currentMode.cutFrequency=cutoff(width, height, m, n);
			currentMode.firstcoor=m;
			currentMode.secondcoor=n;
			currentMode.type=0;
			*(*(a+n)+m)=currentMode;
			//and TM if m and n are not zero
			if(m!=0 && n!=0){
				Nmodes++;
				currentMode.type=1;
				*(*(a+n)+m+1)=currentMode;
				m++;
			}
		}
	}
	//Now the modes are sorted, using the fact that each row is ordered.
	sect->modes=(mode*)malloc(N*sizeof(mode));
	int el=0;
	int currentEl[max_n];
	memset(currentEl,0,max_n*sizeof(int));
	//the search stops when the required number of modes is found
	//each time, the lowest unselected mode in each row is selected,
	//then the lowest among them is selected, and then updated so next time the next mode will be taken.
	while(el< N){
		int min=0;
		for(int n=0;n<=max_n;n++){
			if(currentEl[min]==(nElems[min]-1)) continue;
			min=((*(a+min)+currentEl[min])->cutFrequency < (*(a+n)+currentEl[n])->cutFrequency)?min:n;
		}
		currentEl[min]++;
		*((sect->modes)+el)=*(*(a+min)+currentEl[min]);
		el++;
	}
	sect->Nmodes=N;
	//free the rows
	for(int n=0;n<=max_n;n++){
		free(*(a+n));
	}
	//free the pointers to pointers
	free(a);
	return sect;
}

section *rectangularFreq(section *sect, double cutFreq){
	double width= sect->width;
	double height= sect->height;

	//double max_m= floor(2*width*sqrt(e*er*u*ur*cutFreq));
	//Maximum n given the cutoff frequency.
	int max_n=(int)floor(2*height*sqrt(e*er*u*ur*cutFreq));//floor(max_m*height/width);
	//array of indexes for the num of elements in each row (max_n rows)
	int nElems[max_n];
	memset(nElems, 0, max_n*sizeof(int)); //init to zero
	sect->Nmodes=0; //zero modes found so far

	//sect.modes=(mode **)malloc(max_n*sizeof(*mode));
	mode **a= (mode **)malloc(max_n*sizeof(mode*));
	for(int n=0; n<=max_n; ++n){
		//the maximum m FOR THIS ROW! (local maximum)
		int m_max=(int)width*sqrt((4*e*u*er*ur*cutFreq*cutFreq-(n*n/(height*height))));
		nElems[n]=m_max;
		sect->Nmodes +=m_max;
		//the TM modes are taken into account
		if (n!=0){
			sect->Nmodes += (m_max-1>0)? m_max-1:0;
			nElems[n]+=(m_max-1>0)? m_max-1:0;
		}
		//the row is allocated in memory
		*(a+n)=(mode *)malloc(m_max*sizeof(mode));
		//the row is filed with the modes
		for(int m=0;m<=m_max;m++){
			mode currentMode;
			//TE
			currentMode.cutFrequency=cutoff(width, height, m, n);
			currentMode.firstcoor=m;
			currentMode.secondcoor=n;
			currentMode.type=0;
			*(*(a+n)+m)=currentMode;
			//and TM if m and n are not zero
			if(m!=0 && n!=0){
				sect->Nmodes++;
				currentMode.type=1;
				*(*(a+n)+m+1)=currentMode;
				m++;
			}
		}
	}
	//Now the modes are sorted, using the fact that each row is ordered.
	sect->modes=(mode*)malloc(sect->Nmodes*sizeof(mode));
	int el=0;
	int currentEl[max_n];
	memset(currentEl,0,max_n*sizeof(int));
	//All the modes are sorted
	//each time, the lowest unselected mode in each row is selected,
	//then the lowest among them is selected, and then updated so next time the next mode will be taken.
	while(el< (sect->Nmodes)){
		int min=0;
		for(int n=0;n<=max_n;n++){
			if(currentEl[min]==(nElems[min]-1)) continue;
			min=((*(a+min)+currentEl[min])->cutFrequency < (*(a+n)+currentEl[n])->cutFrequency)?min:n;
		}
		currentEl[min]++;
		*((sect->modes)+el)=*(*(a+min)+currentEl[min]);
		el++;
	}
	//free the rows
	for(int n=0;n<=max_n;n++){
		free(*(a+n));
	}
	//free the pointers to pointers
	free(a);
	return sect;
}

double rectangularrectangular(section *sect1,section *sect2){
	int Eplane = (sect1->height != sect2->height)? TRUE:FALSE;
	int Hplane = (sect1->width != sect2->width)? TRUE:FALSE;
//integral_x0^x1 sin(k (x-x2)) cos(d (x-x3)) dx
//(-d*sin(d*(x0-x3))*sin(k*(x0-x2))-k*cos(d*(x0-x3))*cos(k*(x0-x2))+d*sin(d*(x1-x3))*sin(k*(x1-x2))+k*cos(d*(x1-x3))*cos(k*(x1-x2)))/(2*(d-k)*(d+k));
//integral_x0^x1 cos(k (x-x2)) cos(d (x-x3)) dx
//-((d-k)*(sin(x0*(d+k)-d*x3-k*x2)-sin(x1*(d+k)-d*x3-k*x2))+(d+k)*sin(d*(x0-x3)+k*(x2-x0))-(d+k)*sin(d*(x1-x3)+k*(x2-x1)))/(2*(d-k)*(d+k));
	double X[sect1->Nmodes][sect2->Nmodes];

	for(int i=0;i<sect1->Nmodes;i++){
		for(int j=0;j<sect2->Nmodes;j++){

			double kx1=(*sect1->modes+i)->firstcoor*pi/sect1->width;
			double kx2=(*sect2->modes+j)->firstcoor*pi/sect2->width;
			double ky1=(*sect1->modes+i)->secondcoor*pi/sect1->height;
			double ky2=(*sect2->modes+j)->secondcoor*pi/sect2->height;
			//integral_x0^x1 cos(kx1 (x-x2)) cos(kx2 (x-x3)) dx
			double i1=-((kx2-kx1)*(sin(x0*(kx2+kx1)-kx2*x3-kx1*x2)-sin(x1*(kx2+kx1)-kx2*x3-kx1*x2))+(kx2+kx1)*sin(kx2*(x0-x3)+kx1*(x2-x0))-(kx2+kx1)*sin(kx2*(x1-x3)+kx1*(x2-x1)))/(2*(kx2-kx1)*(kx2+kx1));
			//integral_y0^y1 sin(ky1(y-y2)) cos(ky2(y-y3)) dy
			double i2=(-ky2*sin(ky2*(y0-y3))*sin(ky1*(y0-y2))-ky1*cos(ky2*(y0-y3))*cos(ky1*(y0-y2))+ky2*sin(ky2*(y1-y3))*sin(ky1*(y1-y2))+ky1*cos(ky2*(y1-y3))*cos(ky1*(y1-y2)))/(2*(ky2-ky1)*(ky2+ky1));
			//integral_x0^x1 sin(kx1(x-x2)) cos(kx2(x-x3)) dx
			double i3=(-kx2*sin(kx2*(x0-x3))*sin(kx1*(x0-x2))-kx1*cos(kx2*(x0-x3))*cos(kx1*(x0-x2))+kx2*sin(kx2*(x1-x3))*sin(kx1*(x1-x2))+kx1*cos(kx2*(x1-x3))*cos(kx1*(x1-x2)))/(2*(kx2-kx1)*(kx2+kx1));
			//integral_y0^y1 cos(ky1 (y-y2)) cos(ky2 (y-y3)) dy
			double i4=-((ky2-ky1)*(sin(y0*(ky2+ky1)-ky2*y3-ky1*y2)-sin(y1*(ky2+ky1)-ky2*y3-ky1*y2))+(ky2+ky1)*sin(ky2*(y0-y3)+ky1*(y2-y0))-(ky2+ky1)*sin(ky2*(y1-y3)+ky1*(y2-y1)))/(2*(ky2-ky1)*(ky2+ky1));

			if (Eplane && Hplane){

			}else if (Eplane){

			}else{

			}

		}
	}
	return 0;
}


double cutoff(double width, double height, int m, int n)
{
	return 1.0/((sqrt(e*u*er*ur))*2.0)*sqrt((m)/width*(m)/width+(n)/height*(n)/height);
}
