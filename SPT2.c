#include "stdio.h"
#include "math.h"
#include "stdlib.h"

#define OHS 0
#define HS 1 

double  rho_0, rho_1, rho_2; // rho_0 is the density of the matrix, rho_1 and rho_2 are the densities of fluid.
double  sigma_0, sigma_1, sigma_2; // sigma_0 is the diameter of the matrix, sigma_1 and sigma_2 are the diameters os fluid.
double  xi_m1, xi_m2, xi_m3;
double  xi_f1, xi_f2, xi_f3;
double  x1, x2;

double P1(double sigma, int flag);
double a_alpha(double sigma, int flag);
double b_alpha(double sigma, int flag);
double c_alpha(double sigma, int flag);
double e_alpha(double sigma, int flag);
double Pressure_SPT2a(int flag);
double Pressure_SPT2(int flag);
double miu_SPT2a (double sigma, double rho, int flag);
double miu_SPT2b (double sigma, double rho, int flag);
double miu_SPT2b1 (double sigma, double rho, int flag);
double miu_SPT2c (double sigma, double rho, int flag);
double miu_SPT2d (double sigma, double rho, int flag);
double miu_SPT2 (double sigma, double rho, int flag);


int main (int argc, char **argv){
FILE * fp;
int		flag; // flag = 0 is OHS matrix, flag = 1 is HS matrix.
double  rho; // rho = rho_1 + rho_2.
double  miu1_GCMC, miu2_GCMC;

/*
if((fp=fopen(argv[1],"wt"))==NULL){    
	printf("cannot open %s\n",argv[1]);    
	exit(1);    
}
*/

printf("0 is OHS, 1 is HS:\n");
scanf("%d",&flag);
printf("the density of matrix: \n");
scanf("%lf",&rho_0);
printf("the total density of fluid : \n");
scanf("%lf",&rho);
printf("the concentration of species 1: \n");
scanf("%lf",&x1);
printf("the DIAMETER of matrix:\n");
scanf("%lf",&sigma_0);
printf("the DIAMETER of compoment 1:\n");
scanf("%lf",&sigma_1);
printf("the DIAMETER of compoment 2:\n");
scanf("%lf",&sigma_2);


	rho_1 = rho * x1;
	rho_2 = rho - rho_1;
	x2  = 1-x1;

	xi_f1 = M_PI*(x1 * sigma_1 + x2 * sigma_2)*rho/6.0;
	xi_f2 = M_PI*(x1 * sigma_1 * sigma_1 + x2 * sigma_2 * sigma_2)*rho/6.0;
	xi_f3 = M_PI*(x1 * sigma_1 * sigma_1 * sigma_1 + x2 * sigma_2 * sigma_2 * sigma_2)*rho/6.0;

	xi_m1 = M_PI*rho_0*sigma_0/6.0;
	xi_m2 = M_PI*rho_0*sigma_0*sigma_0/6.0;
	xi_m3 = M_PI*rho_0*sigma_0*sigma_0*sigma_0/6.0;
/*
    printf("component 1: spt2a: %lf ( %lf ) \t spt2b: %lf (%lf) \t spt2c:  %lf (%lf) \t  spt2d:  %lf (%lf) \t  spt2: %lf (%lf)\t \n", miu_SPT2a(sigma_1, rho_1, flag), miu_SPT2a(sigma_1, rho_1, flag)-miu1_GCMC,\
			miu_SPT2b(sigma_1, rho_1, flag),  miu_SPT2b(sigma_1, rho_1, flag)-miu1_GCMC,\
			miu_SPT2c(sigma_1, rho_1, flag),  miu_SPT2c(sigma_1, rho_1, flag)-miu1_GCMC,\
			miu_SPT2d(sigma_1, rho_1, flag),  miu_SPT2d(sigma_1, rho_1, flag)-miu1_GCMC,\
			miu_SPT2(sigma_1, rho_1, flag),   miu_SPT2(sigma_1, rho_1, flag)-miu1_GCMC);

    printf("component 2: spt2a: %lf ( %lf ) \t spt2b: %lf (%lf) \t spt2c:  %lf (%lf) \t  spt2d:  %lf (%lf) \t  spt2: %lf (%lf)\t \n", miu_SPT2a(sigma_2, rho_2, flag), miu_SPT2a(sigma_2, rho_2, flag)-miu2_GCMC,\
			miu_SPT2b(sigma_2, rho_2, flag),  miu_SPT2b(sigma_2, rho_2, flag)-miu2_GCMC,\
			miu_SPT2c(sigma_2, rho_2, flag),  miu_SPT2c(sigma_2, rho_2, flag)-miu2_GCMC,\
			miu_SPT2d(sigma_2, rho_2, flag),  miu_SPT2d(sigma_2, rho_2, flag)-miu2_GCMC,\
			miu_SPT2(sigma_2, rho_2, flag),   miu_SPT2(sigma_2, rho_2, flag)-miu2_GCMC);
*/
	    printf("component 1:\n	spt2a:	%lf	\n	spt2b:	%lf	\n	spt2b1:	%lf	\n	spt2c:	%lf	\n	spt2d:	%lf	\n	spt2:	%lf	\n",
			miu_SPT2a(sigma_1, rho_1, flag),\
			miu_SPT2b(sigma_1, rho_1, flag),\
			miu_SPT2b1(sigma_1, rho_1, flag),\
			miu_SPT2c(sigma_1, rho_1, flag),\
			miu_SPT2d(sigma_1, rho_1, flag),\
			miu_SPT2(sigma_1, rho_1, flag));

       printf("component 2:\n	spt2a:	%lf	\n	spt2b:	%lf	\n	spt2b1:	%lf	\n	spt2c:	%lf	\n	spt2d:	%lf	\n	spt2:	%lf	\n", 
			miu_SPT2a(sigma_2, rho_2, flag),\
			miu_SPT2b(sigma_2, rho_2, flag),\
			miu_SPT2b1(sigma_2, rho_2, flag),\
			miu_SPT2c(sigma_2, rho_2, flag),\
			miu_SPT2d(sigma_2, rho_2, flag),\
			miu_SPT2(sigma_2, rho_2, flag));

	return 0;
}


// The probability of inserting a scaled particle into an empty matrix
double P1(double sigma, int flag) {
	double P0;
	double t0;
	P0 =rho_0 * (1.0+xi_m3+xi_m3*xi_m3)/pow(1.0-xi_m3,3.0);
	if(flag == HS) {
		t0 = -log(1.0-xi_m3)+3.0*xi_m2/(1.0-xi_m3)*sigma+(3.0*xi_m1/(1.0-xi_m3)+9.0/2.0*xi_m2*xi_m2/pow(1.0-xi_m3,2.0))*sigma*sigma+M_PI*sigma*sigma*sigma*P0/6.0;
	}else if(flag == OHS){
		t0 = xi_m3*(1 + sigma / sigma_0)*(1 + sigma / sigma_0)*(1 + sigma / sigma_0);
		}
	return(exp(-1.0*t0));	
}

double a_alpha(double sigma, int flag){
	double x;
	double R;
	double D1, D2;
	R = 0.5 * sigma;
	if(flag == OHS){
		D1 = -3.0 * xi_m3 / (sigma_0 / 2);
		D2 = -6.0 * xi_m3 / (sigma_0 / 2) / (sigma_0 / 2);
		} else if(flag == HS) {
			D1 = -6.0 * xi_m2 / (1.0 - xi_m3);
			D2 = -24.0 * xi_m1 / (1.0 - xi_m3) - 36.0 * xi_m2 * xi_m2 / (1.0 - xi_m3) / (1.0 - xi_m3);
			}
	x = 6.0 * xi_f2 * R / xi_f3	+ 12 * xi_f1 * R * R / xi_f3 - (1.0 + 6.0 * xi_f2 * R / xi_f3) * R * D1 + 0.5 * ((R * D1 * R * D1) - R * R * D2);
	return x;
}

double b_alpha(double sigma, int flag){
	double x;
	double R;
	double D1, D2;
	R = 0.5 * sigma;
	if(flag == OHS){
		D1 = -3.0 * xi_m3 / (sigma_0 / 2);
		D2 = -6.0 * xi_m3 / (sigma_0 / 2) / (sigma_0 / 2);
		} else if(flag == HS) {
			D1 = -6.0 * xi_m2 / (1.0 - xi_m3);
			D2 = -24.0 * xi_m1 / (1.0 - xi_m3) - 36.0 * xi_m2 * xi_m2 / (1.0 - xi_m3) / (1.0 - xi_m3);
			}
	x = 0.5 * (R * D1 - 6.0 * xi_f2 * R / xi_f3) * 	(R * D1 - 6.0 * xi_f2 * R / xi_f3);
	return x;
}

double c_alpha(double sigma, int flag){
	double x;
	double R;
	double D1, D2;
	R = 0.5 * sigma;
	if(flag == OHS){
		D1 = -3.0 * xi_m3 / (sigma_0 / 2);
		D2 = -6.0 * xi_m3 / (sigma_0 / 2) / (sigma_0 / 2);
		} else if(flag == HS) {
			D1 = -6.0 * xi_m2 / (1.0 - xi_m3);
			D2 = -24.0 * xi_m1 / (1.0 - xi_m3) - 36.0 * xi_m2 * xi_m2 / (1.0 - xi_m3) / (1.0 - xi_m3);
			}
	x = -6.0 * R *R * D1 * xi_f2 / xi_f3  + 0.5 * R * R *D1 * D1;
	return x; 
}

double e_alpha(double sigma, int flag){
	double x;
	double R;
	double D1, D2;
	R = 0.5 * sigma;
	if(flag == OHS){
		D1 = -3.0 * xi_m3 / (sigma_0 / 2);
		D2 = -6.0 * xi_m3 / (sigma_0 / 2) / (sigma_0 / 2);
		} else if(flag == HS) {
			D1 = -6.0 * xi_m2 / (1.0 - xi_m3);
			D2 = -24.0 * xi_m1 / (1.0 - xi_m3) - 36.0 * xi_m2 * xi_m2 / (1.0 - xi_m3) / (1.0 - xi_m3);
			}
	x = R * D1 + 0.5 * R * R * D2;
	return x;	
}

double Pressure_SPT2a(int flag){
    double x;
    double phi_0;
    double A, B;
    if(flag == OHS){
		phi_0 = P1(0.0, OHS);
		A = x1 * a_alpha(sigma_1,OHS) + x2 * a_alpha(sigma_2, OHS);
		B = x1 * b_alpha(sigma_1,OHS) + x2 * b_alpha(sigma_2, OHS);
		} else if(flag == HS){
		phi_0 = P1(0.0, HS);
		A = x1 * a_alpha(sigma_1,HS) + x2 * a_alpha(sigma_2, HS);
		B = x1 * b_alpha(sigma_1,HS) + x2 * b_alpha(sigma_2, HS);
		}
   x = 1.0 / (1.0 - xi_f3 / phi_0) + 0.5 * A * (xi_f3 / phi_0)/(1.0 - xi_f3 / phi_0)/(1.0 - xi_f3 / phi_0) + 2.0 * B * 	(xi_f3 / phi_0) * (xi_f3 / phi_0) /(1.0 - xi_f3 / phi_0)/(1.0 - xi_f3 / phi_0)/	(1.0 - xi_f3 / phi_0)/3.0;
   x = (rho_1 + rho_2) * x;
   return x;	
}

double Pressure_SPT2(int flag){
	double phi_0, phi;
	double A, B;
	double y;
	double t1, t2, t3, t4;

	    if(flag == OHS){
		phi_0 = P1(0.0, OHS);
		A = x1 * a_alpha(sigma_1,OHS) + x2 * a_alpha(sigma_2, OHS);
		B = x1 * b_alpha(sigma_1,OHS) + x2 * b_alpha(sigma_2, OHS);
		} else if(flag == HS){
		phi_0 = P1(0.0, HS);
		A = x1 * a_alpha(sigma_1,HS) + x2 * a_alpha(sigma_2, HS);
		B = x1 * b_alpha(sigma_1,HS) + x2 * b_alpha(sigma_2, HS);
		} 
	y = M_PI * sigma_1 * sigma_1 * sigma_1 * rho_1 / 6.0 / P1(sigma_1, flag) +  M_PI * sigma_2 * sigma_2 * sigma_2 * rho_2 / 6.0 / P1(sigma_2, flag) ;
	phi = xi_f3 / y;

	t1 = -1.0 * phi * log(1.0 - xi_f3 / phi) / xi_f3;
	t2 = ((1.0 + A) * phi / (phi - phi_0)) * (phi * log(1.0 - xi_f3 / phi) / xi_f3 - phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3);
	t3 = ((A + 2.0 * B) * phi / (phi - phi_0)) * (1.0 / (1.0 - xi_f3 / phi_0) +  phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3 - (phi / (phi - phi_0)) * (phi * log(1.0 - xi_f3 / phi) / xi_f3 - phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3));
	t4 = (2.0 * B * phi / (phi - phi_0)) * ((xi_f3/phi_0)/2.0/(1.0 - xi_f3 / phi_0)/(1.0 - xi_f3 / phi_0) - 1.0 / (1.0 - xi_f3 / phi_0) -  (phi / (phi - phi_0))*(1.0 / (1.0 - xi_f3 / phi_0) + phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3) + (phi / (phi - phi_0))*(phi / (phi - phi_0))*(phi * log(1.0 - xi_f3 / phi) / xi_f3 - phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3));
	return (t1+t2+t3+t4);

}


double miu_SPT2a (double sigma, double rho, int flag){
	double x;
    double phi_0;
    double C, E;
	double D1, D2;
    double R = 0.5 * sigma;
    if(flag == OHS){
		phi_0 = P1(0.0, OHS);
		C = x1 * c_alpha(sigma_1,OHS) + x2 * c_alpha(sigma_2, OHS);
		E = x1 * e_alpha(sigma_1,OHS) + x2 * e_alpha(sigma_2, OHS);
		D1 = -3.0 * xi_m3 / (sigma_0 / 2);
		D2 = -6.0 * xi_m3 / (sigma_0 / 2) / (sigma_0 / 2);
		} else if(flag == HS){
		phi_0 = P1(0.0, HS);
		C = x1 * c_alpha(sigma_1,HS) + x2 * c_alpha(sigma_2, HS);
		E = x1 * e_alpha(sigma_1,HS) + x2 * e_alpha(sigma_2, HS);
		D1 = -6.0 * xi_m2 / (1.0 - xi_m3);
		D2 = -24.0 * xi_m1 / (1.0 - xi_m3) - 36.0 * xi_m2 * xi_m2 / (1.0 - xi_m3) / (1.0 - xi_m3);
		} 
	
	x = log(rho) - log(P1(sigma, flag)) - log(1.0 - xi_f3/phi_0) + a_alpha(sigma, flag) * (xi_f3 / phi_0)/(1.0 - xi_f3 / phi_0) + b_alpha(sigma, flag)* (xi_f3 / phi_0) * (xi_f3 / phi_0) /(1.0 - xi_f3 / phi_0)/(1.0 - xi_f3 / phi_0)  \
	 + 4.0 * M_PI * R * R * R * Pressure_SPT2a(flag) / 3.0 / phi_0 \
	 + (4.0 * M_PI * R * R * R * (rho_1 + rho_2) * (C - E) / 3.0 / xi_f3 - c_alpha(sigma,flag) + e_alpha(sigma,flag)) * (xi_f3 / phi_0)/(1.0 - xi_f3 / phi_0)/2.0 \
	 + (4.0 * M_PI * R * R * R * (rho_1 + rho_2) * C / 3.0 / xi_f3 - c_alpha(sigma,flag)) * 2.0 * (xi_f3 / phi_0) * (xi_f3 / phi_0) /(1.0 - xi_f3 / phi_0)/(1.0 - xi_f3 / phi_0)/3.0 \
	 + (sigma * sigma * sigma / (xi_f3 / phi_0) - sigma * sigma / (xi_f2 /phi_0)) * ((xi_f3 / phi_0)/(1.0 - xi_f3 / phi_0)/2.0 + (xi_f3 / phi_0) * (xi_f3 / phi_0) /(1.0 - xi_f3 / phi_0)/(1.0 - xi_f3 / phi_0)/3.0) * 3.0 * D1 * (xi_f2 / phi_0) * (xi_f2 / phi_0) / 2.0 / (xi_f3 / phi_0);
	 return x;	   
	}

double miu_SPT2b (double sigma, double rho, int flag){
	double x;
	double phi_0, phi;
	double y;
	double R = 0.5 * sigma;
    if(flag == OHS){
		phi_0 = P1(0.0, OHS);
		} else if(flag == HS){
		phi_0 = P1(0.0, HS);
		} 
	y = M_PI * sigma_1 * sigma_1 * sigma_1 * rho_1 / 6.0 / P1(sigma_1, flag) +  M_PI * sigma_2 * sigma_2 * sigma_2 * rho_2 / 6.0 / P1(sigma_2, flag) ;
	phi = xi_f3 / y;
	x = miu_SPT2a(sigma,rho,flag) - log(1 - xi_f3/phi) + log(1- xi_f3/phi_0) + (4.0 * M_PI * R * R * R * (rho_1 + rho_2) / 3.0 / xi_f3 - 1.0) \
	* (-1.0 * phi * log(1.0 - xi_f3 / phi) / xi_f3 + phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3 ) \
	- (phi * log(1.0 - xi_f3 / phi) / xi_f3 + 1.0) * (phi / P1(sigma,flag) - 1.0) * 4.0 * M_PI * R * R * R * (rho_1 + rho_2) / 3.0 / xi_f3;
	   

	return x;
}

double miu_SPT2b1 (double sigma, double rho, int flag){
	double x;
	double phi_0, phi;
	double y;
	double R = 0.5 * sigma;
    if(flag == OHS){
		phi_0 = P1(0.0, OHS);
		} else if(flag == HS){
		phi_0 = P1(0.0, HS);
		} 
	y = M_PI * sigma_1 * sigma_1 * sigma_1 * rho_1 / 6.0 / P1(sigma_1, flag) +  M_PI * sigma_2 * sigma_2 * sigma_2 * rho_2 / 6.0 / P1(sigma_2, flag) ;
	phi = xi_f3 / y;
	x = miu_SPT2a(sigma,rho,flag) + (phi/phi_0 - 1.0) * (phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3 + 1.0) \
		- (4.0 * M_PI * R * R * R * (rho_1 + rho_2) / 3.0 / xi_f3) *  (phi/phi_0 - 1.0) * (phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3 + 1.0/(1.0 - xi_f3 / phi_0)) \
		+ (4.0 * M_PI * R * R * R * (rho_1 + rho_2) / 3.0 / xi_f3) * (phi/phi_0)*(1.0-phi/P1(sigma,flag))*(phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3 + 1.0);
	
	return x;
}


double miu_SPT2c (double sigma, double rho, int flag){
	double x;
	double phi_0, phi;
	double y;
	double R = 0.5 * sigma;
	double A, D1;

    if(flag == OHS){
		phi_0 = P1(0.0, OHS);
		A = x1 * a_alpha(sigma_1,OHS) + x2 * a_alpha(sigma_2, OHS);
		D1 = -3.0 * xi_m3 / (sigma_0 / 2);
		} else if(flag == HS){
		phi_0 = P1(0.0, HS);
		A = x1 * a_alpha(sigma_1,HS) + x2 * a_alpha(sigma_2, HS);
		D1 = -6.0 * xi_m2 / (1.0 - xi_m3);
		} 
	y = M_PI * sigma_1 * sigma_1 * sigma_1 * rho_1 / 6.0 / P1(sigma_1, flag) +  M_PI * sigma_2 * sigma_2 * sigma_2 * rho_2 / 6.0 / P1(sigma_2, flag) ;
	phi = xi_f3 / y;
	
	x = miu_SPT2b(sigma,rho,flag) + (1 + a_alpha(sigma,flag) + 3.0 * xi_f1 * xi_f2 * (sigma  / xi_f1 + sigma * sigma  / xi_f2 - 2.0 * sigma * sigma * sigma / xi_f3) / xi_f3 + 3 * D1 * xi_f2 * xi_f2 * (sigma * sigma * sigma  / xi_f3 - sigma * sigma / xi_f2) / xi_f3 / 2 )	\
	* (1.0 + phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3 + (phi/(phi-phi_0)) * ((1 - phi / xi_f3) * log(1.0 - xi_f3 / phi) - (1.0 - phi_0 / xi_f3) * log(1.0 - xi_f3 / phi_0))) \
	+ (4.0 * M_PI * R * R * R * (rho_1 + rho_2) / 3.0 / xi_f3) * ((1.0 + A) * (-1.0 * phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3  - 1.0 / (1.0 - xi_f3 / phi_0) + phi * (phi * log(1.0 - xi_f3 / phi) / xi_f3 - phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3) / (phi - phi_0))) \
	+ (1.0 + A) * (4.0 * M_PI * R * R * R * (rho_1 + rho_2)/ 3.0 / xi_f3) * (phi /  P1(sigma_1, flag) -1.0) * (phi / (phi - phi_0)) * (1.0 + phi * log(1.0 - xi_f3 / phi) / xi_f3 + (phi_0 / (phi - phi_0)) * ((1.0 - phi/xi_f3)*log(1.0 - xi_f3 / phi) - (1.0 - phi_0 / xi_f3) * log(1.0 - xi_f3 / phi_0)));
	
	return x;
}

double miu_SPT2d (double sigma, double rho, int flag){
	double x;
	double phi_0, phi;
	double y;
	double R = 0.5 * sigma;
	double A, B, D1;
	double t1, t2, t3, t4;
	
    if(flag == OHS){
		phi_0 = P1(0.0, OHS);
		A = x1 * a_alpha(sigma_1,OHS) + x2 * a_alpha(sigma_2, OHS);
		B = x1 * b_alpha(sigma_1,OHS) + x2 * b_alpha(sigma_2, OHS);
		D1 = -3.0 * xi_m3 / (sigma_0 / 2);
		} else if(flag == HS){
		phi_0 = P1(0.0, HS);
		A = x1 * a_alpha(sigma_1,HS) + x2 * a_alpha(sigma_2, HS);
		B = x1 * b_alpha(sigma_1,HS) + x2 * b_alpha(sigma_2, HS);
		D1 = -6.0 * xi_m2 / (1.0 - xi_m3);
		} 
	y = M_PI * sigma_1 * sigma_1 * sigma_1 * rho_1 / 6.0 / P1(sigma_1, flag) +  M_PI * sigma_2 * sigma_2 * sigma_2 * rho_2 / 6.0 / P1(sigma_2, flag) ;
	phi = xi_f3 / y;	
	
	t1 = a_alpha(sigma,flag) + 2.0 * b_alpha(sigma,flag) + (3.0 * xi_f1 * xi_f2 / xi_f3 ) * (sigma / xi_f1 + sigma * sigma / xi_f2 - 2.0 * sigma * sigma * sigma / xi_f3)\
	 + (9.0 * xi_f2 * xi_f2 / 2.0 / xi_f3) * (D1 - 4.0 * xi_f2 / xi_f3) * (sigma * sigma * sigma / xi_f3 - sigma * sigma / xi_f2);
	 
    t2 = (1.0 + phi/(phi - phi_0)) * (phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3 + 1.0) + xi_f3 / phi_0 / 2.0 / (1.0 - xi_f3 / phi_0) + (phi/(phi - phi_0)) * (phi/(phi - phi_0)) \
    * ((1.0 - phi / xi_f3) * log(1.0 - xi_f3 / phi) - (1.0 - phi_0 / xi_f3) * log(1.0 - xi_f3 / phi_0));
    
    t3 = (A + 2.0 * B) * (4.0 * M_PI * R * R * R * (rho_1 + rho_2) / 3.0 / xi_f3) * (phi /  P1(sigma, flag) - 1.0) * (phi/(phi - phi_0)) * ((phi_0/(phi - phi_0)) * (phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3  + 1.0) \
    + (phi/(phi - phi_0)) * (phi * log(1.0 - xi_f3 / phi) / xi_f3 + 1.0) + (2.0 * phi_0 * phi / (phi - phi_0) / (phi - phi_0)) * ((1.0 - phi / xi_f3) * log(1.0 - xi_f3 / phi) - (1.0 - phi_0 / xi_f3) * log(1.0 - xi_f3 / phi_0)));
    
    t4 = (4.0 * M_PI * R * R * R * (rho_1 + rho_2) / 3.0 / xi_f3) * (A + 2.0 * B) * ((1.0 + phi/(phi - phi_0)) * (phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3 + 1.0 / (1.0 - xi_f3 / phi_0)) \
    - (xi_f3 / phi_0 / 2.0 / (1.0 - xi_f3 / phi_0) / (1.0 - xi_f3 / phi_0))  - (phi/(phi - phi_0)) * (phi/(phi - phi_0)) * ((phi / xi_f3) * log(1.0 - xi_f3 / phi) - (phi_0 / xi_f3) * log(1.0 - xi_f3 / phi_0)));
    
    x = miu_SPT2c(sigma,rho,flag) - t1 * t2 + t3 + t4;

    return x;
}

double miu_SPT2 (double sigma, double rho, int flag) {
	double x;
	double phi_0, phi;
	double y;
	double R = 0.5 * sigma;
	double A, B, D1;
	double t1, t2, t3, t4;
	
    if(flag == OHS){
		phi_0 = P1(0.0, OHS);
		A = x1 * a_alpha(sigma_1,OHS) + x2 * a_alpha(sigma_2, OHS);
		B = x1 * b_alpha(sigma_1,OHS) + x2 * b_alpha(sigma_2, OHS);
		D1 = -3.0 * xi_m3 / (sigma_0 / 2);
		} else if(flag == HS){
		phi_0 = P1(0.0, HS);
		A = x1 * a_alpha(sigma_1,HS) + x2 * a_alpha(sigma_2, HS);
		B = x1 * b_alpha(sigma_1,HS) + x2 * b_alpha(sigma_2, HS);
		D1 = -6.0 * xi_m2 / (1.0 - xi_m3);
		} 
	y = M_PI * sigma_1 * sigma_1 * sigma_1 * rho_1 / 6.0 / P1(sigma_1, flag) +  M_PI * sigma_2 * sigma_2 * sigma_2 * rho_2 / 6.0 / P1(sigma_2, flag) ;
	phi = xi_f3 / y;
	
	t1  = 2.0 * b_alpha(sigma,flag) + (3.0 * xi_f2 * xi_f2 / xi_f3) * (D1 - 6.0 * xi_f2 / xi_f3) * (sigma * sigma * sigma / xi_f3 - sigma * sigma / xi_f2);
	
	t2 = (1.0 + phi/(phi - phi_0) + (phi/(phi - phi_0)) * (phi/(phi - phi_0))) * ((phi_0 / xi_f3) * log(1.0 - xi_f3 / phi_0) + 1.0) \
	+ (1 + phi/(phi - phi_0)) * (xi_f3 / phi_0 / 2.0 / (1.0 - xi_f3 / phi_0)) \
	- (xi_f3 / phi_0 / (1.0 - xi_f3 / phi_0)) * (xi_f3 / phi_0 / (1.0 - xi_f3 / phi_0)) / 6.0 \
	+ ((phi/(phi - phi_0)) * (phi/(phi - phi_0)) * (phi/(phi - phi_0))) * ((1.0 - phi / xi_f3) * log(1.0 - xi_f3 / phi) - (1.0 - phi_0 / xi_f3) * log(1.0 - xi_f3 / phi_0));
	
	t3 = 2.0 * B * (4.0 * M_PI * R * R * R * (rho_1 + rho_2) / 3.0 / xi_f3) * (phi / P1(sigma,flag) - 1.0) * (phi/(phi - phi_0)) \
	* (  (phi_0/(phi - phi_0)) * ( 1.0 + 2.0 * phi / (phi - phi_0) ) * (  phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3 + 1.0)  +  (phi/(phi - phi_0)) * (phi/(phi - phi_0)) * ( phi * log(1.0 - xi_f3 / phi) / xi_f3 + 1.0 ) \
	+ ( phi_0/(phi - phi_0) ) * ( xi_f3 / phi_0 / 2.0 / ( 1.0 - xi_f3/phi_0) ) + ( 3.0 * phi_0 * phi * phi / (phi - phi_0) / (phi - phi_0) / (phi - phi_0) ) * ( (1.0 - phi / xi_f3) * log(1.0 - xi_f3 / phi) - (1.0 - phi_0 / xi_f3) * log(1.0 - xi_f3 / phi_0) ));
	
	t4 = ( 4.0 * M_PI * R * R * R * (rho_1 + rho_2) / 3.0 / xi_f3 ) \
	* 2.0 * B * ( -1.0 * (1.0 + (phi / (phi - phi_0)) + (phi / (phi - phi_0)) * (phi / (phi - phi_0)) ) \
	* (phi_0 * log(1.0 - xi_f3 / phi_0) / xi_f3 + 1.0 / ( 1.0 - xi_f3 / phi_0) ) \
	+ (1.0 + phi / (phi - phi_0)) * (xi_f3 / phi_0) / 2.0 / (1.0 - xi_f3 / phi_0) / (1.0 - xi_f3 / phi_0) \
	- (xi_f3/phi_0)*(xi_f3/phi_0)/3.0/(1.0 - xi_f3 / phi_0)/(1.0 - xi_f3 / phi_0)/(1.0 - xi_f3 / phi_0)\
	+ ((phi/(phi - phi_0)) * (phi/(phi - phi_0)) * (phi/(phi - phi_0))) * ((phi / xi_f3) * log(1.0 - xi_f3 / phi) - (phi_0 / xi_f3) * log(1.0 - xi_f3 / phi_0) )) ;
	
	x = miu_SPT2d(sigma,rho,flag) + t1 * t2 + t3 + t4;
	
	return x;
	}

