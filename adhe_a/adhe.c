/*
 * This code calculates the Hubble constant based on the luminosity distance-redshift relation of a light source. 
 * A perturbation with amplitude A and width L exists between the observer and the light source. 
 * The structure formation is calculated using the adhesion model.
 * Multiple calculations are performed with different values of A.
 * The quantities for each value of A are recorded in the file output_aa.dat.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* fundamental constant */
#define C_SI 2.99792458e8                         /* light speed (m/s) */
#define C (C_SI/3.085677581e22)                   /* light speed (Mpc/s)　*/
#define PI 3.14159265359                          /*　 circle ratio */
#define G_SI 6.6740831e-11                        /* gravitational constant (m^3 kg^{-1} s^{-2}) */
#define G (G_SI/pow(3.085677581e22,3)/C/C)        /* gravitational constant (Mpc/kg) */

/* Hubble constant of background */
#define H 0.67                         /* dimensionless Hubble constant */
#define H0 (100*H/(C_SI*1.0e-3))       /* Hubble constant of background universe (Mpc^{-1}) */

/* present value */
#define RHO0 (3.0*H0*H0/(8.0*PI*G))     /* present energy density (kg Mpc^{-3}) */
#define T0 (2.0/(3.0*H0))               /* present cosmic time (Mpc) */

/* equality time */
#define Z_EQ 3411.0                     /* redshift of equality time */
#define A_EQ (1.0/(1.0 + Z_EQ))         /* scale factor of equality time */

/* parameters of perturbation */
#define L 10.0           /* width of perturbation (Mpc) */

/* functions */

/* scale factor */
double a(double t) {
	return 1.5*H0*pow(T0,1.0/3.0)*pow(t,2.0/3.0);
}

/* growth factor */
double b(double t) { 
    return a(t)/A_EQ;
}

/* da/dt */
double a_dot(double t) { 
    return H0*pow(T0,1.0/3.0)*pow(t,-1.0/3.0);
}

/* db/dt */
double b_dot(double t) { 
    return a_dot(t)/A_EQ;
}

/* d^2 a/dt^2 */
double a_ddot(double t) { 
    return -(1.0/3.0)*H0*pow(T0,1.0/3.0)*pow(t,-4.0/3.0);
}

/* d^2 b/dt^2 */
double b_ddot(double t) { 
    return a_ddot(t)/A_EQ;
}

/* velocity potential of fluid (Mpc^2) */
double Phi(double q, double A) {
	return -A*L*L*exp(-q*q/(2*L*L));
}

/* d \Phi/dq */
double dPhi(double q, double A) {
	return A*q*exp(-q*q/(2*L*L));
}

/* d^2 \Phi/dq^2 */
double ddPhi(double q, double A) {
	return A*(1-q*q/(L*L))*exp(-q*q/(2*L*L)); 
}

/* gravitational potential */
double phi(double t, double q, double A) {
	return 1.5*a_dot(t)*a_dot(t)*b(t)*(Phi(q,A) - 0.5*b(t)*dPhi(q,A)*dPhi(q,A));
}

/* partial derivative of phi with respect ot x */
double phi_x(double t, double q, double A) {
	return 1.5*a_dot(t)*a_dot(t)*b(t)*dPhi(q,A);
}

/* phi at t = 0 */
double phi_in(double q, double A) { 
    return 1.5*1.5*H0*H0*H0*T0/A_EQ*Phi(q,A); 
}

/* velocity of fluid */
double v(double t, double q, double A) {
	return -b_dot(t)*dPhi(q,A);
}

/* partial derivative of v with respect to t */
double vt(double t, double q, double A) {
	return -b_ddot(t)*dPhi(q,A) - b_dot(t)*ddPhi(q,A)*v(t,q,A)/(1 - b(t)*ddPhi(q,A));
}

/* partial derivative of v with respect to x */
double vx(double t, double q, double A) {
	return -b_dot(t)*ddPhi(q,A)/(1 - b(t)*ddPhi(q,A));
}

/* energy density of fluid (kg Mpc^{-3})*/
double rho(double t, double q, double A) {
	return RHO0*pow(a(t),-3.0)/(1 - b(t)*ddPhi(q,A));
}

/* u_{tx} */
double utx(double t, double q, double A) {
	return -a(t)*a(t)*vx(t,q,A)*v(t,q,A) - a_dot(t)*a(t)*v(t,q,A);
}

/* u_{xt} */
double uxt(double t, double q, double A) {
	return a(t)*a_dot(t)*v(t,q,A) + a(t)*a(t)*vt(t,q,A) + phi_x(t,q,A);
}

/* u_{xx} */
double uxx(double t, double q, double A) {
	return a(t)*a(t)*vx(t,q,A) + a(t)*a_dot(t)*(1 - 3*phi(t,q,A));
}

/* k^{t} */
double kt(double t, double q, double z, double A) {
	return (1+z)*(1 - a(t)*v(t,q,A) - phi(t,q,A) + 0.5*a(t)*a(t)*v(t,q,A)*v(t,q,A));
}

/* k^{x} */
double kx(double t, double q, double z, double A) {
	return -(1+z)*(1 - a(t)*v(t,q,A) + phi(t,q,A) + 0.5*a(t)*a(t)*v(t,q,A)*v(t,q,A))/a(t);
}

/* dt/dq */
double fA(double t, double q, double z, double A) {
	return kt(t,q,z,A)*(1 - b(t)*ddPhi(q,A))/(kx(t,q,z,A) + b_dot(t)*dPhi(q,A)*kt(t,q,z,A));
}

/* dz/dq (Mpc^{-1]}) */
double fB(double t, double q, double z, double A) {
	return  -(kt(t,q,z,A)*kx(t,q,z,A)*utx(t,q,A) + kx(t,q,z,A)*kt(t,q,z,A)*uxt(t,q,A) + kx(t,q,z,A)*kx(t,q,z,A)*uxx(t,q,A))*(1 - b(t)*ddPhi(q,A))/(kx(t,q,z,A) + b_dot(t)*dPhi(q,A)*kt(t,q,z,A));
}

/* d D_A/dq */
double fC(double t, double q, double z, double E_A, double A) {
	return  E_A*(1 - b(t)*ddPhi(q, A))/(kx(t,q,z,A) + b_dot(t)*dPhi(q,A)*kt(t,q,z,A));
}

/* d^2 D_A/dq^2 (Mpc^{-1}) */
double fD(double t, double q, double z, double D_A, double A) {
	return -4*PI*G*rho(t,q,A)*D_A*(1+z)*(1+z)*(1 - b(t)*ddPhi(q,A))/(kx(t,q,z,A) + b_dot(t)*dPhi(q,A)*kt(t,q,z,A));
}

/* smaller boundary value of collapsed region (Mpc) */
double q1(double t, double A) {
	return -L*sqrt(2*log(A*b(t)));
}

/* larger boundary value of collapsed region (Mpc) */
double q2(double t, double A) {
	return -q1(t,A);
}


int main(void) {

	/* for A's loop */
	int i_a;
	int N_a = 300;

	/* open output file */
	FILE *fp;

    fp = fopen("output_aa.dat","w");
    if (fp == NULL) {

        perror("Error opening file\n");
        exit(1);
    }

	for (i_a = 0; i_a < N_a; i_a++) {

		/* amplitude of perturbation */
		double A = pow(10,(3.0*i_a/N_a - 5));

		/* spatial resolution */
  		int N_q = 100000;
		double dq = 0.1;

    	/* initial value of q */
    	int j_i = 0;
    	double q0 = -100;   /* (Mpc) */
    	double q = q0 + j_i*dq;

    	/* label of q */
    	int j;

		/* for calculating precise variables of junction and source */
    	int jc1;
    	int N_JC1 = 1000000; 
    	double dq_jc1 = dq*1.0e-3;

    	int jc2;
    	int N_JC2 = 1000000; 
    	double dq_jc2 = dq_jc1*1.0e-3;

		/* redshift of source */
    	double z_s = 0.1;

    	/* initial conditions of light propagation */
    	double t = T0;        /* (Mpc) */
    	double z = 0;
    	double D_A = 0;
    	double E_A = -1;      /* dD_A/dq */


    	/* for fourth-order Runge-Kutta method */
    	double kA1,kA2,kA3,kA4;
    	double kB1,kB2,kB3,kB4;
    	double kC1,kC2,kC3,kC4;
    	double kD1,kD2,kD3,kD4;
    	double dz,dt,dD_A,dE_A;

    	/* light propagation from the present time */
    	for (j = j_i; j < N_q; j++) {
        	if (z > z_s) {
            	break; 
        	}

        	if (A*b(t) > 1) { 
            	if (q > q1(t,A)) {
                	break;
            	}    
        	}

			kA1 = fA(t,q,z,A);
			kB1 = fB(t,q,z,A);
			kC1 = fC(t,q,z,E_A,A);
			kD1 = fD(t,q,z,D_A,A);

			kA2 = fA(t + dq*kA1/2.0, q + dq/2.0, z + dq*kB1/2.0, A);
			kB2 = fB(t + dq*kA1/2.0, q + dq/2.0, z + dq*kB1/2.0, A);
			kC2 = fC(t + dq*kA1/2.0, q + dq/2.0, z + dq*kB1/2.0, E_A + dq*kD1/2.0, A);
			kD2 = fD(t + dq*kA1/2.0, q + dq/2.0, z + dq*kB1/2.0, D_A + dq*kC1/2.0, A);
 
			kA3 = fA(t + dq*kA2/2.0, q + dq/2.0, z + dq*kB2/2.0, A);
			kB3 = fB(t + dq*kA2/2.0, q + dq/2.0, z + dq*kB2/2.0, A);
			kC3 = fC(t + dq*kA2/2.0, q + dq/2.0, z + dq*kB2/2.0, E_A + dq*kD2/2.0, A);
			kD3 = fD(t + dq*kA2/2.0, q + dq/2.0, z + dq*kB2/2.0, D_A + dq*kC2/2.0, A);

			kA4 = fA(t + dq*kA3, q + dq, z + dq*kB3, A);
			kB4 = fB(t + dq*kA3, q + dq, z + dq*kB3, A);
			kC4 = fC(t + dq*kA3, q + dq, z + dq*kB3, E_A + dq*kD3, A);
			kD4 = fD(t + dq*kA3, q + dq, z + dq*kB3, D_A + dq*kC3, A);

			dt = dq*(kA1 + 2*kA2 + 2*kA3 + kA4)/6.0;
			dz = dq*(kB1 + 2*kB2 + 2*kB3 + kB4)/6.0;
			dD_A = dq*(kC1 + 2*kC2 + 2*kC3 + kC4)/6.0;
			dE_A = dq*(kD1 + 2*kD2 + 2*kD3 + kD4)/6.0;

			t = t + dt;
			z = z + dz;
			D_A = D_A + dD_A;
			E_A = E_A + dE_A;
			q = q + dq;
		}

		if (z < z_s) {

			/* go back one step */
  			t = t - dt;
  			z = z - dz;
  			D_A = D_A - dD_A;
  			E_A = E_A - dE_A;
  			q = q - dq;

			/* reduce step size */
			for (jc1 = 0; j < N_JC1; jc1++) {
				if (z > z_s) {
					break;
				}

				if (A*b(t) > 1) {
    				if (q > q1(t,A)) {
						break;
					}
				}

				kA1 = fA(t,q,z,A);
				kB1 = fB(t,q,z,A);
				kC1 = fC(t,q,z,E_A,A);
				kD1 = fD(t,q,z,D_A,A);

				kA2 = fA(t + dq_jc1*kA1/2.0, q + dq_jc1/2.0, z + dq_jc1*kB1/2.0, A);
				kB2 = fB(t + dq_jc1*kA1/2.0, q + dq_jc1/2.0, z + dq_jc1*kB1/2.0, A);
				kC2 = fC(t + dq_jc1*kA1/2.0, q + dq_jc1/2.0, z + dq_jc1*kB1/2.0, E_A + dq_jc1*kD1/2.0, A);
				kD2 = fD(t + dq_jc1*kA1/2.0, q + dq_jc1/2.0, z + dq_jc1*kB1/2.0, D_A + dq_jc1*kC1/2.0, A);
  
				kA3 = fA(t + dq_jc1*kA2/2.0, q + dq_jc1/2.0, z + dq_jc1*kB2/2.0, A);
				kB3 = fB(t + dq_jc1*kA2/2.0, q + dq_jc1/2.0, z + dq_jc1*kB2/2.0, A);
				kC3 = fC(t + dq_jc1*kA2/2.0, q + dq_jc1/2.0, z + dq_jc1*kB2/2.0, E_A + dq_jc1*kD2/2.0, A);
				kD3 = fD(t + dq_jc1*kA2/2.0, q + dq_jc1/2.0, z + dq_jc1*kB2/2.0, D_A + dq_jc1*kC2/2.0, A);
  
				kA4 = fA(t + dq_jc1*kA3, q + dq_jc1, z + dq_jc1*kB3, A);
				kB4 = fB(t + dq_jc1*kA3, q + dq_jc1, z + dq_jc1*kB3, A);
				kC4 = fC(t + dq_jc1*kA3, q + dq_jc1, z + dq_jc1*kB3, E_A + dq_jc1*kD3, A);
				kD4 = fD(t + dq_jc1*kA3, q + dq_jc1, z + dq_jc1*kB3, D_A + dq_jc1*kC3, A);

				dt = dq_jc1*(kA1 + 2*kA2 + 2*kA3 + kA4)/6.0;
				dz = dq_jc1*(kB1 + 2*kB2 + 2*kB3 + kB4)/6.0;
				dD_A = dq_jc1*(kC1 + 2*kC2 + 2*kC3 + kC4)/6.0;
				dE_A = dq_jc1*(kD1 + 2*kD2 + 2*kD3 + kD4)/6.0;
  
				t = t + dt;
				z = z + dz;
				D_A = D_A + dD_A;
				E_A = E_A + dE_A;
				q = q + dq_jc1;
			}
  
			/* go back one step */
			t = t - dt;
			z = z - dz;
			D_A = D_A - dD_A;
			E_A = E_A - dE_A;
			q = q - dq_jc1;

			/* reduce step size */
			for (jc2 = 0; j < N_JC2; jc2++) {
				if (z > z_s) {
					break;
				}

				if (A*b(t) > 1) {
    				if (q > q1(t,A)) {
						break;
					}
				}

				kA1 = fA(t,q,z,A);
   				kB1 = fB(t,q,z,A);
   				kC1 = fC(t,q,z,E_A,A);
   				kD1 = fD(t,q,z,D_A,A);

				kA2 = fA(t + dq_jc2*kA1/2.0, q + dq_jc2/2.0, z + dq_jc2*kB1/2.0, A);
				kB2 = fB(t + dq_jc2*kA1/2.0, q + dq_jc2/2.0, z + dq_jc2*kB1/2.0, A);
				kC2 = fC(t + dq_jc2*kA1/2.0, q + dq_jc2/2.0, z + dq_jc2*kB1/2.0, E_A + dq_jc2*kD1/2.0, A);
				kD2 = fD(t + dq_jc2*kA1/2.0, q + dq_jc2/2.0, z + dq_jc2*kB1/2.0, D_A + dq_jc2*kC1/2.0, A);

				kA3 = fA(t + dq_jc2*kA2/2.0, q + dq_jc2/2.0, z + dq_jc2*kB2/2.0, A);
				kB3 = fB(t + dq_jc2*kA2/2.0, q + dq_jc2/2.0, z + dq_jc2*kB2/2.0, A);
				kC3 = fC(t + dq_jc2*kA2/2.0, q + dq_jc2/2.0, z + dq_jc2*kB2/2.0, E_A + dq_jc2*kD2/2.0, A);
				kD3 = fD(t + dq_jc2*kA2/2.0, q + dq_jc2/2.0, z + dq_jc2*kB2/2.0, D_A + dq_jc2*kC2/2.0, A);
 
				kA4 = fA(t + dq_jc2*kA3, q + dq_jc2, z + dq_jc2*kB3, A);
				kB4 = fB(t + dq_jc2*kA3, q + dq_jc2, z + dq_jc2*kB3, A);
				kC4 = fC(t + dq_jc2*kA3, q + dq_jc2, z + dq_jc2*kB3, E_A + dq_jc2*kD3, A);
				kD4 = fD(t + dq_jc2*kA3, q + dq_jc2, z + dq_jc2*kB3, D_A + dq_jc2*kC3, A);
 
				dt = dq_jc2*(kA1 + 2*kA2 + 2*kA3 + kA4)/6.0;
				dz = dq_jc2*(kB1 + 2*kB2 + 2*kB3 + kB4)/6.0;
				dD_A = dq_jc2*(kC1 + 2*kC2 + 2*kC3 + kC4)/6.0;
				dE_A = dq_jc2*(kD1 + 2*kD2 + 2*kD3 + kD4)/6.0;

				t = t + dt;
				z = z + dz;
				D_A = D_A + dD_A;
				E_A = E_A + dE_A;
				q = q + dq_jc2;
			}

			/* junction condition at sheet */
 			E_A += -4*PI*G*D_A*kt(t,q,z,A)*kt(t,q,z,A)/kx(t,q,z,A)*(q2(t,A)-q1(t,A))*RHO0/a(t)/a(t)/a(t);
 			z += -a(t)*a(t)*kx(t,q,z,A)*(v(t,q2(t,A),A)-v(t,q1(t,A),A));
			q = -q; 

			/* light propagation from sheet */
			for (; j < N_q; j++) {
				if (z > z_s) {
					break;
				}
   
				kA1 = fA(t,q,z,A);
				kB1 = fB(t,q,z,A);
				kC1 = fC(t,q,z,E_A,A);
				kD1 = fD(t,q,z,D_A,A);
  
				kA2 = fA(t + dq*kA1/2.0, q + dq/2.0, z + dq*kB1/2.0, A);
				kB2 = fB(t + dq*kA1/2.0, q + dq/2.0, z + dq*kB1/2.0, A);
				kC2 = fC(t + dq*kA1/2.0, q + dq/2.0, z + dq*kB1/2.0, E_A + dq*kD1/2.0, A);
				kD2 = fD(t + dq*kA1/2.0, q + dq/2.0, z + dq*kB1/2.0, D_A + dq*kC1/2.0, A);
 
				kA3 = fA(t + dq*kA2/2.0, q + dq/2.0, z + dq*kB2/2.0, A);
				kB3 = fB(t + dq*kA2/2.0, q + dq/2.0, z + dq*kB2/2.0, A);
				kC3 = fC(t + dq*kA2/2.0, q + dq/2.0, z + dq*kB2/2.0, E_A + dq*kD2/2.0, A);
				kD3 = fD(t + dq*kA2/2.0, q + dq/2.0, z + dq*kB2/2.0, D_A + dq*kC2/2.0, A);
  
				kA4 = fA(t + dq*kA3, q + dq, z + dq*kB3, A);
				kB4 = fB(t + dq*kA3, q + dq, z + dq*kB3, A);
				kC4 = fC(t + dq*kA3, q + dq, z + dq*kB3, E_A + dq*kD3, A);
				kD4 = fD(t + dq*kA3, q + dq, z + dq*kB3, D_A + dq*kC3, A);

				dt = dq*(kA1 + 2*kA2 + 2*kA3 + kA4)/6.0;
				dz = dq*(kB1 + 2*kB2 + 2*kB3 + kB4)/6.0;
				dD_A = dq*(kC1 + 2*kC2 + 2*kC3 + kC4)/6.0;
				dE_A = dq*(kD1 + 2*kD2 + 2*kD3 + kD4)/6.0;

				t = t + dt;
				z = z + dz;
				D_A = D_A + dD_A;
				E_A = E_A + dE_A;
				q = q + dq;
			}
		}

		/* Go back one step */
		t = t - dt;
		z = z - dz;
		D_A = D_A - dD_A;
		E_A = E_A - dE_A;
		q = q - dq;

		/* reduce step size */
		for (jc1 = 0; j < N_JC1; jc1++) {
  			if (z > z_s) {
				break;
			}

			kA1 = fA(t,q,z,A);
			kB1 = fB(t,q,z,A);
			kC1 = fC(t,q,z,E_A,A);
			kD1 = fD(t,q,z,D_A,A);

			kA2 = fA(t + dq_jc1*kA1/2.0, q + dq_jc1/2.0, z + dq_jc1*kB1/2.0, A);
			kB2 = fB(t + dq_jc1*kA1/2.0, q + dq_jc1/2.0, z + dq_jc1*kB1/2.0, A);
			kC2 = fC(t + dq_jc1*kA1/2.0, q + dq_jc1/2.0, z + dq_jc1*kB1/2.0, E_A + dq_jc1*kD1/2.0, A);
			kD2 = fD(t + dq_jc1*kA1/2.0, q + dq_jc1/2.0, z + dq_jc1*kB1/2.0, D_A + dq_jc1*kC1/2.0, A);

			kA3 = fA(t + dq_jc1*kA2/2.0, q + dq_jc1/2.0, z + dq_jc1*kB2/2.0, A);
			kB3 = fB(t + dq_jc1*kA2/2.0, q + dq_jc1/2.0, z + dq_jc1*kB2/2.0, A);
			kC3 = fC(t + dq_jc1*kA2/2.0, q + dq_jc1/2.0, z + dq_jc1*kB2/2.0, E_A + dq_jc1*kD2/2.0, A);
			kD3 = fD(t + dq_jc1*kA2/2.0, q + dq_jc1/2.0, z + dq_jc1*kB2/2.0, D_A + dq_jc1*kC2/2.0, A);
   
			kA4 = fA(t + dq_jc1*kA3, q + dq_jc1, z + dq_jc1*kB3, A);
			kB4 = fB(t + dq_jc1*kA3, q + dq_jc1, z + dq_jc1*kB3, A);
			kC4 = fC(t + dq_jc1*kA3, q + dq_jc1, z + dq_jc1*kB3, E_A + dq_jc1*kD3, A);
			kD4 = fD(t + dq_jc1*kA3, q + dq_jc1, z + dq_jc1*kB3, D_A + dq_jc1*kC3, A);

			dt = dq_jc1*(kA1 + 2*kA2 + 2*kA3 + kA4)/6.0;
			dz = dq_jc1*(kB1 + 2*kB2 + 2*kB3 + kB4)/6.0;
			dD_A = dq_jc1*(kC1 + 2*kC2 + 2*kC3 + kC4)/6.0;
			dE_A = dq_jc1*(kD1 + 2*kD2 + 2*kD3 + kD4)/6.0;

			t = t + dt;
			z = z + dz;
			D_A = D_A + dD_A;
			E_A = E_A + dE_A;
			q = q + dq_jc1;
		}

		/* go back one step */
		t = t - dt;
		z = z - dz;
		D_A = D_A - dD_A;
		E_A = E_A - dE_A;
		q = q - dq_jc1;

		/* reduce step size */
		for (jc2 = 0; j < N_JC2; jc2++) {
			if (z > z_s) {
				break;
			}

			kA1 = fA(t,q,z,A);
			kB1 = fB(t,q,z,A);
			kC1 = fC(t,q,z,E_A,A);
			kD1 = fD(t,q,z,D_A,A);

			kA2 = fA(t + dq_jc2*kA1/2.0, q + dq_jc2/2.0, z + dq_jc2*kB1/2.0, A);
			kB2 = fB(t + dq_jc2*kA1/2.0, q + dq_jc2/2.0, z + dq_jc2*kB1/2.0, A);
			kC2 = fC(t + dq_jc2*kA1/2.0, q + dq_jc2/2.0, z + dq_jc2*kB1/2.0, E_A + dq_jc2*kD1/2.0, A);
			kD2 = fD(t + dq_jc2*kA1/2.0, q + dq_jc2/2.0, z + dq_jc2*kB1/2.0,  D_A + dq_jc2*kC1/2.0, A);

			kA3 = fA(t + dq_jc2*kA2/2.0, q + dq_jc2/2.0, z + dq_jc2*kB2/2.0, A);
			kB3 = fB(t + dq_jc2*kA2/2.0, q + dq_jc2/2.0, z + dq_jc2*kB2/2.0, A);
			kC3 = fC(t + dq_jc2*kA2/2.0, q + dq_jc2/2.0, z + dq_jc2*kB2/2.0, E_A + dq_jc2*kD2/2.0, A);
			kD3 = fD(t + dq_jc2*kA2/2.0, q + dq_jc2/2.0, z + dq_jc2*kB2/2.0,  D_A + dq_jc2*kC2/2.0, A);

			kA4 = fA(t + dq_jc2*kA3, q + dq_jc2, z + dq_jc2*kB3, A);
			kB4 = fB(t + dq_jc2*kA3, q + dq_jc2, z + dq_jc2*kB3, A);
			kC4 = fC(t + dq_jc2*kA3, q + dq_jc2, z + dq_jc2*kB3, E_A + dq_jc2*kD3, A);
			kD4 = fD(t + dq_jc2*kA3, q + dq_jc2, z + dq_jc2*kB3, D_A + dq_jc2*kC3, A);

			dt = dq_jc2*(kA1 + 2*kA2 + 2*kA3 + kA4)/6.0;
			dz = dq_jc2*(kB1 + 2*kB2 + 2*kB3 + kB4)/6.0;
			dD_A = dq_jc2*(kC1 + 2*kC2 + 2*kC3 + kC4)/6.0;
			dE_A = dq_jc2*(kD1 + 2*kD2 + 2*kD3 + kD4)/6.0;

			t = t + dt;
			z = z + dz;
			D_A = D_A + dD_A;
			E_A = E_A + dE_A;
			q = q + dq_jc2;
		}

		/* find the Hubble constant */
		double D_L;				/* luminosity distance */
		double H0_I;			/* Hubble constant estimated from luminosity distance (Mpc^{-1}) */
		double delta_H0;		/* relative difference between H0 and H0_I */

		D_L = (1+z)*(1+z)*D_A;
		H0_I = (2*(1+z)/D_L)*(1.0 - 1.0/sqrt(1+z));
		delta_H0 = sqrt(pow((H0_I - H0)/H0, 2.0));

		fprintf(fp, "%e %e %e\n", A, delta_H0, sqrt(phi_in(0,A)*phi_in(0,A)));

		printf( "\rProgress : %03d / %d" , i_a+1, N_a );
		fflush(stdout);
	}

  	fclose(fp);

	printf("\nCalculation complete.\n");

    return 0;
}

