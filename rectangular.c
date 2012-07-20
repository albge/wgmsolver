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
	double max_m = (1+alpha+sqrt((1+alpha)*(1+alpha)+32*N*pi*alpha))/(4*pi*(1+alpha));
	//max n
	//double max_n = max_m/alpha;

	double cutFreq = 1.0/((sqrt(e*u*er*ur))*2.0)*max_m/width;

	/*FROM HERE ON, ALMOST IDENTICAL TO rectangularFreq(section *sect, int N)*/


	//double max_m= floor(2*width*sqrt(e*er*u*ur*cutFreq));
	/*Once the max fc has been found, the maximum n is inmediate*/
	int max_n=(int)floor(2*height*sqrt(e*er*u*ur*cutFreq));//floor(max_m*height/width);
	//array of indexes for the num of elements in each row (max_n rows)
	int nElems[max_n+1];
	memset(&nElems, 0, (max_n+1)*sizeof(int)); //init to zero

	int Nmodes=-1; //zero modes, so TE00 is skipped
	//sect->Nmodes=0; //zero modes found so far
	//sect.modes=(mode **)malloc(max_n*sizeof(*mode));
	mode **a= (mode **)malloc((max_n+1)*sizeof(mode *));
	if(a == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
	}

	for(int n=0; n<=max_n; ++n){
		//the maximum m FOR THIS ROW! (local maximum)
		int m_max=(int)floor(width*sqrt(4*e*u*er*ur*cutFreq*cutFreq-((n*n)/(height*height))));
		nElems[n]=m_max+1;
		Nmodes+=m_max+1;
		//the TM modes are taken into account
		if (n!=0){
			nElems[n]+=m_max;
			Nmodes+=m_max;
		}
		//the row is allocated in memory
		//int j=sizeof(struct mode *);
		//int k=sizeof(mode);
		(*(a+n))=(mode *)malloc(nElems[n]*sizeof(mode));
		if((a+n) == NULL){
			fprintf(stderr, "out of memory\n");
			exit(EXIT_FAILURE);
		}
		//the row is filed with the modes
		int sub_m=0;
		for(int m=0;m<=m_max;m++){
			if(n==0 && m==0){ sub_m++; continue;}
			mode currentMode;
			//TE
			currentMode.cutFrequency=1.0/((sqrt(e*u*er*ur))*2.0)*sqrt((m)/width*(m)/width+(n)/height*(n)/height);
			currentMode.firstcoor=m;
			currentMode.secondcoor=n;
			currentMode.type=0;
			*(*(a+n)+sub_m)=currentMode;
			sub_m++;
			//and TM if m and n are not zero
			if(m!=0 && n!=0){
				//sect->Nmodes++;
				currentMode.type=1;
				*(*(a+n)+sub_m)=currentMode;
				sub_m++;
			}
		}
	}

	//Now the modes are sorted, using the fact that each row is ordered.
	sect->modes=(mode*)malloc(N*sizeof(mode));
	sect->Nmodes=N;
	int el=0;
	int currentEl[max_n+1];
	memset(&currentEl,0,(max_n+1)*sizeof(int));
	currentEl[0]=1;//TE00 skip
	//All the modes are sorted
	//each time, the lowest unselected mode in each row is selected,
	//then the lowest among them is selected, and then updated so next time the next mode will be taken.
	while(el<= N){
		int min=0;
		for(int n=0;n<=max_n;n++){
			//if(currentEl[min]==(nElems[min])) continue;
			min=((*(a+min)+currentEl[min])->cutFrequency <= (*(a+n)+currentEl[n])->cutFrequency)?min:n;
		}
		*((sect->modes)+el)=*(*(a+min)+currentEl[min]);
		currentEl[min]+=1;
		el++;
	}
	//free the rows
	for(int n=0;n<=max_n+1;n++){
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
	int max_n=(int)floor(2*height*sqrt(e*er*u*ur)*cutFreq);//floor(max_m*height/width);

	//array of indexes for the num of elements in each row (max_n rows)
	int nElems[max_n+1];
	memset(&nElems, 0, (max_n+1)*sizeof(int)); //init to zero

	int Nmodes=-1; //zero modes, so TE00 is skipped
	//sect->Nmodes=0; //zero modes found so far
	//sect.modes=(mode **)malloc(max_n*sizeof(*mode));
	mode **a= (mode **)malloc((max_n+1)*sizeof(mode *));
	if(a == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
	}

	for(int n=0; n<=max_n; ++n){

		//the maximum m FOR THIS ROW! (local maximum)
		int m_max=(int)floor(width*sqrt(4*e*u*er*ur*cutFreq*cutFreq-((n*n)/(height*height))));
		nElems[n]=m_max+1;
		Nmodes+=m_max+1;
		//the TM modes are taken into account
		if (n!=0){
			nElems[n]+=m_max;
			Nmodes+=m_max;
		}

		//the row is allocated in memory
		//int j=sizeof(struct mode *);
		//int k=sizeof(mode);
		(*(a+n))=(mode *)malloc(nElems[n]*sizeof(mode));
		if((a+n) == NULL){
			fprintf(stderr, "out of memory\n");
			exit(EXIT_FAILURE);
		}
		//the row is filed with the modes

		int sub_m=0;
		for(int m=0;m<=m_max;m++){
			if(n==0 && m==0){ sub_m++; continue;}
			mode currentMode;
			//TE
			currentMode.cutFrequency=1.0/((sqrt(e*u*er*ur))*2.0)*sqrt((m)/width*(m)/width+(n)/height*(n)/height);
			currentMode.firstcoor=m;
			currentMode.secondcoor=n;
			currentMode.type=0;
			*(*(a+n)+sub_m)=currentMode;
			sub_m++;
			//and TM if m and n are not zero
			if(m!=0 && n!=0){
				//sect->Nmodes++;
				currentMode.type=1;
				*(*(a+n)+sub_m)=currentMode;
				sub_m++;
			}
		}
	}

	//Now the modes are sorted, using the fact that each row is ordered.
	sect->modes=(mode*)malloc(Nmodes*sizeof(mode));
	sect->Nmodes=Nmodes;
	int el=0;
	int currentEl[max_n+1];
	memset(&currentEl,0,(max_n+1)*sizeof(int));
	currentEl[0]=1;//TE00 skip
	//All the modes are sorted
	//each time, the lowest unselected mode in each row is selected,
	//then the lowest among them is selected, and then updated so next time the next mode will be taken.
	while(el<= Nmodes){
		int min=0;
		for(int n=0;n<=max_n;n++){
			//if(currentEl[min]==(nElems[min])) continue;
			min=((*(a+min)+currentEl[min])->cutFrequency <= (*(a+n)+currentEl[n])->cutFrequency)?min:n;
		}
		*((sect->modes)+el)=*(*(a+min)+currentEl[min]);
		currentEl[min]+=1;
		el++;
	}
	//free the rows
	for(int n=0;n<=max_n+1;n++){
		free(*(a+n));
	}
	//free the pointers to pointers
	free(a);
	return sect;
}


double complex *rectangularrectangular(section *sect1,section *sect2, double *freqs,double complex *Sparams){
	//num of freq dots
	int elems=0;
	while (*(freqs+elems)!=0)
		elems++;

	double a1 = sect1->height;
	double a2 = sect2->height;
	double b1 = sect1->width;
	double b2 = sect2->width;

	int Eplane = (a1 != a2)? TRUE:FALSE;
	int Hplane = (b1 != b2)? TRUE:FALSE;

	//Maybe &Sparams? npi...
	double complex (*X)[sect1->Nmodes][sect2->Nmodes][elems] = Sparams;

	//for each mode in the guide 1
	for(int i=0;i<sect1->Nmodes;i++){
		//epsilon1&2, for the TE with n=0 or m=0
		int epsi1=((sect1->modes+i)->firstcoor==0)?2:1;
		int epsi2=((sect1->modes+i)->secondcoor==0)?2:1;

		double kx1=(sect1->modes+i)->firstcoor*pi/sect1->width;
		double ky1=(sect1->modes+i)->secondcoor*pi/sect1->height;

		double complex Z1[elems];

		//impedance is filled for each frequency
		for(int f=0;f<elems;f++){
			Z1[f]=((sect1->modes+i)->type==0)?
					(2*pi*f*u*ur)/sqrt((double complex)(4*pi*pi*f*f*u*e*ur*er)-(kx1*kx1)-(ky1*ky1))
					:sqrt((double complex)(4*pi*pi*f*f*u*e*ur*er)-(kx1*kx1)-(ky1*ky1))/(2*pi*f*e*er);
		}

		//idem for the second guide.
		for(int j=0;j<sect2->Nmodes;j++){


			int epsj1=((sect2->modes+j)->firstcoor==0)?2:1;
			int epsj2=((sect2->modes+j)->secondcoor==0)?2:1;

			double kx2=(sect2->modes+j)->firstcoor*pi/sect2->width;
			double ky2=(sect2->modes+j)->secondcoor*pi/sect2->height;

			double kc1=sqrt(kx1*kx1+ky1*ky1);
			double kc2=sqrt(kx2*kx2+ky2*ky2);


			double complex Z2[elems];

			for(int f=0;f<elems;f++){
				Z2[f]=((sect2->modes+j)->type==0)?
						(2*pi*f*u*ur)/sqrt((double complex)(4*pi*pi*f*f*u*e*ur*er)-(kx2*kx2)-(ky2*ky2))
						:sqrt((double complex)(4*pi*pi*f*f*u*e*ur*er)-(kx2*kx2)-(ky2*ky2))/(2*pi*f*e*er);
			}

			//integral_x0^x1 cos(kx1 (x-x2)) cos(kx2 (x-x3)) dx
			//double i1=-((kx2-kx1)*(sin(x0*(kx2+kx1)-kx2*x3-kx1*x2)-sin(x1*(kx2+kx1)-kx2*x3-kx1*x2))+(kx2+kx1)*sin(kx2*(x0-x3)+kx1*(x2-x0))-(kx2+kx1)*sin(kx2*(x1-x3)+kx1*(x2-x1)))/(2*(kx2-kx1)*(kx2+kx1));
			//integral_y0^y1 sin(ky1(y-y2)) cos(ky2(y-y3)) dy
			//double i2=(-ky2*sin(ky2*(y0-y3))*sin(ky1*(y0-y2))-ky1*cos(ky2*(y0-y3))*cos(ky1*(y0-y2))+ky2*sin(ky2*(y1-y3))*sin(ky1*(y1-y2))+ky1*cos(ky2*(y1-y3))*cos(ky1*(y1-y2)))/(2*(ky2-ky1)*(ky2+ky1));
			//integral_x0^x1 sin(kx1(x-x2)) cos(kx2(x-x3)) dx
			//double i3=(-kx2*sin(kx2*(x0-x3))*sin(kx1*(x0-x2))-kx1*cos(kx2*(x0-x3))*cos(kx1*(x0-x2))+kx2*sin(kx2*(x1-x3))*sin(kx1*(x1-x2))+kx1*cos(kx2*(x1-x3))*cos(kx1*(x1-x2)))/(2*(kx2-kx1)*(kx2+kx1));
			//integral_y0^y1 cos(ky1 (y-y2)) cos(ky2 (y-y3)) dy
			//double i4=-((ky2-ky1)*(sin(y0*(ky2+ky1)-ky2*y3-ky1*y2)-sin(y1*(ky2+ky1)-ky2*y3-ky1*y2))+(ky2+ky1)*sin(ky2*(y0-y3)+ky1*(y2-y0))-(ky2+ky1)*sin(ky2*(y1-y3)+ky1*(y2-y1)))/(2*(ky2-ky1)*(ky2+ky1));

			double x0,x1,y0,y1;
			x0= sect1->center.x - sect2->center.x - sect2->width/2;
			x1= x0 + sect2->width;
			y0= sect1->center.y - sect2->center.y - sect2->height/2;
			y1= y0 + sect2->height;

			//integral_x0^x1 cos(kx1 x) cos(kx2 (x-x0)) dx =
			double i1 = (kx2*cos(kx1*x1)*sin(kx2*(x0-x1))+kx1*sin(kx1*x1)*cos(kx2*(x0-x1))-kx1*sin(kx1*x0))/((kx1+kx2)*(kx1-kx2));
			//integral_y0^y1 sin(ky1 y) sin(ky2 (y-y0)) dy =
			double i2 = (ky1*cos(ky1*y1)*sin(ky2*(y0-y1))+ky2*sin(ky1*y1)*cos(ky2*(y0-y1))-ky2*sin(ky1*y0))/((ky1+ky2)*(ky1-ky2));
			//integral_x0^x1 sin(kx1 x) sin(kx2 (x-x0)) dx =
			double i3 = (kx1*cos(kx1*x1)*sin(kx2*(x0-x1))+kx2*sin(kx1*x1)*cos(kx2*(x0-x1))-kx2*sin(kx1*x0))/((kx1+kx2)*(kx1-kx2));
			//integral_y0^y1 cos(ky1 y) cos(ky2 (y-y0)) dy =
			double i4 = (ky2*cos(ky1*y1)*sin(ky2*(y0-y1))+ky1*sin(ky1*y1)*cos(ky2*(y0-y1))-ky1*sin(ky1*y0))/((ky1+ky2)*(ky1-ky2));

			double complex subX[sect1->Nmodes][sect2->Nmodes];
			//double step
			if (Eplane && Hplane){
				//TE with TE
				if((sect1->modes+i)->type==0 && (sect2->modes+i)->type==0){
					subX[i][j]=4*((ky1*ky2)*i1*i2+(kx1*kx2)*i3*i4)
							/(kc1*kc2*sqrt(a1*a2*b1*b2*epsi1*epsi2*epsj1*epsj2));
				//TM with TE
				}else if((sect1->modes+i)->type==1 && (sect2->modes+i)->type==0){
					subX[i][j]=4*(-(kx1*ky2)*i1*i2+(ky1*kx2)*i3*i4)
							/(kc1*kc2*sqrt(a1*a2*b1*b2*epsi1*epsi2));
				//TE with TM
				}else if((sect1->modes+i)->type==0 && (sect2->modes+i)->type==1){
					subX[i][j]=4*(-(ky1*kx2)*i1*i2+(kx1*ky2)*i3*i4)
							/(kc1*kc2*sqrt(a1*a2*b1*b2*epsj1*epsj2));
				//TM with TM
				}else{
					subX[i][j]=4*((kx1*kx2)*i1*i2+(ky1*ky2)*i3*i4)
							/(kc1*kc2*sqrt(a1*a2*b1*b2));
				}

				for(int f=0;f<elems;f++){
					(*X)[i][j][f]=subX[i][j]*Z1[f]/Z2[f];
				}

			}else if (Eplane){

			}else{
				//Hplane
			}

		}
	}
	return (double complex*)X;
}
/**
 * this method returns and array with the propagation constants
 * for all the modes in a section.
 *
 * Params:
 * - *section where the modes belong to
 * - *freqs, frequencies where to calculate the propagation constants
 * - freqSpots number of frequency spots
 */
double complex *recPropConstants(section *section,double *freqs, int Nfreqs){
	double w2[Nfreqs];
	for(int f=0;f<Nfreqs;f++){
		w2[f]=2*pi*((*freqs)*(*freqs));
	}
	int Nmodes=section->Nmodes;
	double complex (*gammas)[Nmodes][Nfreqs]=(double complex (*)[Nmodes][Nfreqs])malloc(Nmodes*Nfreqs*sizeof(double complex));
	if(gammas == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
	}

	for(int p=0;p<Nmodes;p++){
		int m,n;
		m=(section->modes+p)->firstcoor;
		n=(section->modes+p)->secondcoor;
		double a, b;
		a=section->width;
		b=section->height;
		for(int f=0;f<Nfreqs;f++){
			(*gammas)[p][f]=I*sqrt((w2[f]*u*e)-(m*pi/a)*(m*pi/a)-(n*pi/b)*(n*pi/b));
		}
	}
	return gammas;
}
/**
 * this method returns and array with the propagation constant for 1 mode
 *
 * Params:
 * - *section where the modes belong to
 * - mode, the mode of interest
 * - *freqs, frequencies where to calculate the propagation constants
 * - freqSpots number of frequency spots
 */
double complex *recModePropConstant(section *section, mode *mode,double *freqs, int Nfreqs ){
	double w2[Nfreqs];
	double complex (*gammas)[Nfreqs]=(double complex (*)[Nfreqs])malloc(Nfreqs*sizeof(complex));
	if(gammas == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
	}

	int m,n;
	m=mode->firstcoor;
	n=mode->secondcoor;
	double a, b;
	a=section->width;
	b=section->height;
	for(int f=0;f<Nfreqs;f++){
		w2[f]=2*pi*(*freqs*(*freqs));
		(*gammas)[f]=I*sqrt((w2[f]*u*e)-(m*pi/a)*(m*pi/a)-(n*pi/b)*(n*pi/b));
	}
	return gammas;
}
