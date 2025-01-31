/*
 * This code calculates the Hubble constant based on the luminosity distance-redshift relation of a light source. 
 * A perturbation with amplitude A and width L exists between the observer and the light source. 
 * The structure formation is calculated using the freezing model.
 * The physical quantities at each step from the observer to the light source are recorded in the file output_fo.dat.
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
#define A 1.0e-3         /* amplitude of perturbation */
#define L 10.0           /* width of perturbation (Mpc) */

/* functions */

/* scale factor */
double a(double T) {
    return 1.5*H0*pow(T0,1.0/3.0)*pow(T,2.0/3.0);
}

/* growth factor */
double b(double T) {
    return a(T)/A_EQ;
}

/* da/dT */
double a_dot(double T) {
    return H0*pow(T0,1.0/3.0)*pow(T,-1.0/3.0);
}

/* db/dT */
double b_dot(double T) {
    return a_dot(T)/A_EQ;
}

/* velocity potential of fluid (Mpc^2) */
double Phi(double q) {
    return -A*L*L*exp(-q*q/(2*L*L));
}

/* d \Phi/dq */
double dPhi(double q) {
    return A*q*exp(-q*q/(2*L*L));
}

/* d^2 \Phi/dq^2 */
double ddPhi(double q) {
    return A*(1-q*q/(L*L))*exp(-q*q/(2*L*L)); 
}

/* gravitational potential */
double phi(double T,double q) {
    return 1.5*a_dot(T)*a_dot(T)*b(T)*(Phi(q) - 0.5*b(T)*dPhi(q)*dPhi(q));
}

/* phi at T = 0 */
double phi_in(double q) { 
    return 1.5*1.5*H0*H0*H0*T0/A_EQ*Phi(q); 
}

/* energy density of fluid (kg Mpc^{-3})*/
double rho(double T, double q) {
    return RHO0*pow(a(T),-3.0)*(1 + 2*phi(T,q))/(1 - b(T)*ddPhi(q));
}

/* expansion of fluid */
double Theta(double T, double q) {
    return 3*a_dot(T)/a(T) - b_dot(T)*ddPhi(q)/(1 - b(T)*ddPhi(q));
}

/* shear scalar of fluid */
double Sigma(double T, double q) {
    return (1.0/3.0)*b_dot(T)*ddPhi(q)/(1 - b(T)*ddPhi(q));
}

/* metric sqrt{gamma_{qq}} */
double gqqSR(double T, double q) {
    return a(T)*(1 - (5.0/3.0)*phi(T,q))*(1 - b(T)*ddPhi(q));
}

/* dT/dq */
double fA(double T, double q) {
    return -gqqSR(T,q);
}

/* dz/dq (Mpc^{-1]}) */
double fB(double T, double q, double z) {
    return (-2*Sigma(T,q) + (1.0/3.0)*Theta(T,q))*(1+z)*gqqSR(T,q);
}

/* d D_A/dq */
double fC(double T, double q, double z, double E_A) {
    return -E_A*gqqSR(T,q)/(1+z);
}

/* d^2 D_A/dq^2 (Mpc^{-1}) */
double fD(double T, double q, double z, double D_A) {
    return 4*PI*G*rho(T,q)*D_A*(1+z)*gqqSR(T,q);
}

/* freezing time. This function is defined only for input values where ddPhi(q) > 0.　*/
double T_FRZ(double q) {
    return pow((3.0/4.0)*A_EQ/ddPhi(q),1.5)*T0;
}

/* energy density of fluid in frozen region */
double rho_FRZ(double q) {
    return rho(T_FRZ(q),q);
}

/* shear scalar of fluid in frozen region */
double Sigma_FRZ(double q) {
    return Sigma(T_FRZ(q),q);
}

/* metric sqrt{gamma_{qq}} in frozen region */
double gqqSR_FRZ(double T, double q) {
    return gqqSR(T_FRZ(q),q)*exp(-4*Sigma_FRZ(q)*(T - T_FRZ(q)));
}

/* dT/dq in frozen region */
double fA_FRZ(double T, double q) {
    return -gqqSR_FRZ(T,q);
}

/* dz/dq in frozen region (Mpc^{-1]}) */
double fB_FRZ(double T, double q, double z) {
    return -2*Sigma_FRZ(q)*(1+z)*gqqSR_FRZ(T,q);
}

/* d D_A/dq in frozen region */
double fC_FRZ(double T, double q, double z, double E_A) {
    return - E_A*gqqSR_FRZ(T,q)/(1+z);
}

/* d^2 D_A/dq^2 in frozen region (Mpc^{-1}) */
double fD_FRZ(double T, double q, double z, double D_A) {
    return 4*PI*G*rho_FRZ(q)*D_A*(1+z)*gqqSR_FRZ(T,q);
}

/* collapsing time of adhesion model */
double t_col(double q) {
    return pow(A_EQ/A*exp(q*q/2.0/L/L),1.5)*T0;
}

/* x in Zel'dovich approximation */
double x(double T, double q) {
    return q - b(T)*dPhi(q);
}


int main(void){

    /* spatial resolution */
    int N_q = 100000;
    double dq = 0.1;    /* (Mpc) */

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
    double T = (1 + phi(T0,q))*T0;          /* solve from t = T0 */
    double z = 0;
    double D_A = 0;
    double E_A = -1;                        /* dD_A/dq */

    /* for fourth-order Runge-Kutta method */
    double kA1,kA2,kA3,kA4;
    double kB1,kB2,kB3,kB4;
    double kC1,kC2,kC3,kC4;
    double kD1,kD2,kD3,kD4;
    double dz,dT,dD_A,dE_A;

    /* open output file */
    FILE *fp;

    fp = fopen("output_fo.dat","w");
    if (fp == NULL) {

        perror("Error opening file\n");
        exit(1);
    }

    fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
            q, z, (1+z)*(1+z)*D_A, (2/H0/a(T))*(1.0 - sqrt(a(T))), 
            ((1+z)*(1+z)*D_A-(2/H0/a(T))*(1.0 - sqrt(a(T))))/((2/H0/a(T))*(1.0 - sqrt(a(T)))), 
            (z-(1.0/a(T) - 1.0))/(1.0/a(T) - 1.0), 
            (D_A-(2*a(T)/H0)*(1.0 - sqrt(a(T))))/((2*a(T)/H0)*(1.0 - sqrt(a(T)))), 
            T, ddPhi(q), dPhi(q), T_FRZ(q), t_col(q), x(T,q));

    /* light propagation from the present time */
    for (j = j_i; j < N_q; j++) {
        if (z > z_s) {
            break;
        }

        if (ddPhi(q) > 0 && T >= T_FRZ(q)) {

            kA1 = fA_FRZ(T,q);
            kB1 = fB_FRZ(T,q,z);
            kC1 = fC_FRZ(T,q,z,E_A);
            kD1 = fD_FRZ(T,q,z,D_A);

            kA2 = fA_FRZ(T + dq*kA1/2.0, q + dq/2.0);
            kB2 = fB_FRZ(T + dq*kA1/2.0, q + dq/2.0, z + dq*kB1/2.0);
            kC2 = fC_FRZ(T + dq*kA1/2.0, q + dq/2.0, z + dq*kB1/2.0, E_A + dq*kD1/2.0);
            kD2 = fD_FRZ(T + dq*kA1/2.0, q + dq/2.0, z + dq*kB1/2.0, D_A + dq*kC1/2.0);

            kA3 = fA_FRZ(T + dq*kA2/2.0, q + dq/2.0);
            kB3 = fB_FRZ(T + dq*kA2/2.0, q + dq/2.0, z + dq*kB2/2.0);
            kC3 = fC_FRZ(T + dq*kA2/2.0, q + dq/2.0, z + dq*kB2/2.0, E_A + dq*kD2/2.0);
            kD3 = fD_FRZ(T + dq*kA2/2.0, q + dq/2.0, z + dq*kB2/2.0, D_A + dq*kC2/2.0);

            kA4 = fA_FRZ(T + dq*kA3, q + dq);
            kB4 = fB_FRZ(T + dq*kA3, q + dq, z + dq*kB3);
            kC4 = fC_FRZ(T + dq*kA3, q + dq, z + dq*kB3, E_A + dq*kD3);
            kD4 = fD_FRZ(T + dq*kA3, q + dq, z + dq*kB3, D_A + dq*kC3);

            dT = dq*(kA1 + 2*kA2 + 2*kA3 + kA4)/6.0;
            dz = dq*(kB1 + 2*kB2 + 2*kB3 + kB4)/6.0;
            dD_A = dq*(kC1 + 2*kC2 + 2*kC3 + kC4)/6.0;
            dE_A = dq*(kD1 + 2*kD2 + 2*kD3 + kD4)/6.0;

            T = T + dT;
            z = z + dz;
            D_A = D_A + dD_A;
            E_A = E_A + dE_A;
            q = q + dq;

            fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
                    q, z, (1+z)*(1+z)*D_A, (2/H0/a(T))*(1.0 - sqrt(a(T))), 
                    ((1+z)*(1+z)*D_A-(2/H0/a(T))*(1.0 - sqrt(a(T))))/((2/H0/a(T))*(1.0 - sqrt(a(T)))), 
                    (z-(1.0/a(T) - 1.0))/(1.0/a(T) - 1.0), 
                    (D_A-(2*a(T)/H0)*(1.0 - sqrt(a(T))))/((2*a(T)/H0)*(1.0 - sqrt(a(T)))), 
                    T, ddPhi(q), dPhi(q), T_FRZ(q), t_col(q), x(T,q));
        }   

        else{

            kA1 = fA(T,q);
            kB1 = fB(T,q,z);
            kC1 = fC(T,q,z,E_A);
            kD1 = fD(T,q,z,D_A);

            kA2 = fA(T + dq*kA1/2.0, q + dq/2.0);
            kB2 = fB(T + dq*kA1/2.0, q + dq/2.0, z + dq*kB1/2.0);
            kC2 = fC(T + dq*kA1/2.0, q + dq/2.0, z + dq*kB1/2.0, E_A + dq*kD1/2.0);
            kD2 = fD(T + dq*kA1/2.0, q + dq/2.0, z + dq*kB1/2.0,  D_A + dq*kC1/2.0);

            kA3 = fA(T + dq*kA2/2.0, q + dq/2.0);
            kB3 = fB(T + dq*kA2/2.0, q + dq/2.0, z + dq*kB2/2.0);
            kC3 = fC(T + dq*kA2/2.0, q + dq/2.0, z + dq*kB2/2.0, E_A + dq*kD2/2.0);
            kD3 = fD(T + dq*kA2/2.0, q + dq/2.0, z + dq*kB2/2.0,  D_A + dq*kC2/2.0);
   
            kA4 = fA(T + dq*kA3, q + dq);
            kB4 = fB(T + dq*kA3, q + dq, z + dq*kB3);
            kC4 = fC(T + dq*kA3, q + dq, z + dq*kB3, E_A + dq*kD3);
            kD4 = fD(T + dq*kA3, q + dq, z + dq*kB3, D_A + dq*kC3);

            dT = dq*(kA1 + 2*kA2 + 2*kA3 + kA4)/6.0;
            dz = dq*(kB1 + 2*kB2 + 2*kB3 + kB4)/6.0;
            dD_A = dq*(kC1 + 2*kC2 + 2*kC3 + kC4)/6.0;
            dE_A = dq*(kD1 + 2*kD2 + 2*kD3 + kD4)/6.0;

            T = T + dT;
            z = z + dz;
            D_A = D_A + dD_A;
            E_A = E_A + dE_A;
            q = q + dq;

            fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
                    q, z, (1+z)*(1+z)*D_A, (2/H0/a(T))*(1.0 - sqrt(a(T))), 
                    ((1+z)*(1+z)*D_A-(2/H0/a(T))*(1.0 - sqrt(a(T))))/((2/H0/a(T))*(1.0 - sqrt(a(T)))), 
                    (z-(1.0/a(T) - 1.0))/(1.0/a(T) - 1.0), 
                    (D_A-(2*a(T)/H0)*(1.0 - sqrt(a(T))))/((2*a(T)/H0)*(1.0 - sqrt(a(T)))), 
                    T, ddPhi(q), dPhi(q), T_FRZ(q), t_col(q), x(T,q));
        }
    }

    /* Go back one step */
    T = T - dT;
    z = z - dz;
    D_A = D_A - dD_A;
    E_A = E_A - dE_A;
    q = q - dq;

	/* reduce step size */
    for (jc1 = 0; j < N_JC1; jc1++) {
        if (z > z_s) {
            break;
        }

        if (ddPhi(q) > 0 && T >= T_FRZ(q)) {

            kA1 = fA_FRZ(T,q);
            kB1 = fB_FRZ(T,q,z);
            kC1 = fC_FRZ(T,q,z,E_A);
            kD1 = fD_FRZ(T,q,z,D_A);

            kA2 = fA_FRZ(T + dq_jc1*kA1/2.0, q + dq_jc1/2.0);
            kB2 = fB_FRZ(T + dq_jc1*kA1/2.0, q + dq_jc1/2.0, z + dq_jc1*kB1/2.0);
            kC2 = fC_FRZ(T + dq_jc1*kA1/2.0, q + dq_jc1/2.0, z + dq_jc1*kB1/2.0, E_A + dq_jc1*kD1/2.0);
            kD2 = fD_FRZ(T + dq_jc1*kA1/2.0, q + dq_jc1/2.0, z + dq_jc1*kB1/2.0, D_A + dq_jc1*kC1/2.0);

            kA3 = fA_FRZ(T + dq_jc1*kA2/2.0, q + dq_jc1/2.0);
            kB3 = fB_FRZ(T + dq_jc1*kA2/2.0, q + dq_jc1/2.0, z + dq_jc1*kB2/2.0);
            kC3 = fC_FRZ(T + dq_jc1*kA2/2.0, q + dq_jc1/2.0, z + dq_jc1*kB2/2.0, E_A + dq_jc1*kD2/2.0);
            kD3 = fD_FRZ(T + dq_jc1*kA2/2.0, q + dq_jc1/2.0, z + dq_jc1*kB2/2.0, D_A + dq_jc1*kC2/2.0);

            kA4 = fA_FRZ(T + dq_jc1*kA3, q + dq_jc1);
            kB4 = fB_FRZ(T + dq_jc1*kA3, q + dq_jc1, z + dq_jc1*kB3);
            kC4 = fC_FRZ(T + dq_jc1*kA3, q + dq_jc1, z + dq_jc1*kB3, E_A + dq_jc1*kD3);
            kD4 = fD_FRZ(T + dq_jc1*kA3, q + dq_jc1, z + dq_jc1*kB3, D_A + dq_jc1*kC3);

            dT = dq_jc1*(kA1 + 2*kA2 + 2*kA3 + kA4)/6.0;
            dz = dq_jc1*(kB1 + 2*kB2 + 2*kB3 + kB4)/6.0;
            dD_A = dq_jc1*(kC1 + 2*kC2 + 2*kC3 + kC4)/6.0;
            dE_A = dq_jc1*(kD1 + 2*kD2 + 2*kD3 + kD4)/6.0;

            T = T + dT;
            z = z + dz;
            D_A = D_A + dD_A;
            E_A = E_A + dE_A;
            q = q + dq_jc1;

            fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
                    q, z, (1+z)*(1+z)*D_A, (2/H0/a(T))*(1.0 - sqrt(a(T))), 
                    ((1+z)*(1+z)*D_A-(2/H0/a(T))*(1.0 - sqrt(a(T))))/((2/H0/a(T))*(1.0 - sqrt(a(T)))), 
                    (z-(1.0/a(T) - 1.0))/(1.0/a(T) - 1.0), 
                    (D_A-(2*a(T)/H0)*(1.0 - sqrt(a(T))))/((2*a(T)/H0)*(1.0 - sqrt(a(T)))), 
                    T, ddPhi(q), dPhi(q), T_FRZ(q), t_col(q), x(T,q));
        }

        else{

            kA1 = fA(T,q);
            kB1 = fB(T,q,z);
            kC1 = fC(T,q,z,E_A);
            kD1 = fD(T,q,z,D_A);

            kA2 = fA(T + dq_jc1*kA1/2.0, q + dq_jc1/2.0);
            kB2 = fB(T + dq_jc1*kA1/2.0, q + dq_jc1/2.0, z + dq_jc1*kB1/2.0);
            kC2 = fC(T + dq_jc1*kA1/2.0, q + dq_jc1/2.0, z + dq_jc1*kB1/2.0, E_A + dq_jc1*kD1/2.0);
            kD2 = fD(T + dq_jc1*kA1/2.0, q + dq_jc1/2.0, z + dq_jc1*kB1/2.0, D_A + dq_jc1*kC1/2.0);

            kA3 = fA(T + dq_jc1*kA2/2.0, q + dq_jc1/2.0);
            kB3 = fB(T + dq_jc1*kA2/2.0, q + dq_jc1/2.0, z + dq_jc1*kB2/2.0);
            kC3 = fC(T + dq_jc1*kA2/2.0, q + dq_jc1/2.0, z + dq_jc1*kB2/2.0, E_A + dq_jc1*kD2/2.0);
            kD3 = fD(T + dq_jc1*kA2/2.0, q + dq_jc1/2.0, z + dq_jc1*kB2/2.0, D_A + dq_jc1*kC2/2.0);

            kA4 = fA(T + dq_jc1*kA3, q + dq_jc1);
            kB4 = fB(T + dq_jc1*kA3, q + dq_jc1, z + dq_jc1*kB3);
            kC4 = fC(T + dq_jc1*kA3, q + dq_jc1, z + dq_jc1*kB3, E_A + dq_jc1*kD3);
            kD4 = fD(T + dq_jc1*kA3, q + dq_jc1, z + dq_jc1*kB3, D_A + dq_jc1*kC3);

            dT = dq_jc1*(kA1 + 2*kA2 + 2*kA3 + kA4)/6.0;
            dz = dq_jc1*(kB1 + 2*kB2 + 2*kB3 + kB4)/6.0;
            dD_A = dq_jc1*(kC1 + 2*kC2 + 2*kC3 + kC4)/6.0;
            dE_A = dq_jc1*(kD1 + 2*kD2 + 2*kD3 + kD4)/6.0;

            T = T + dT;
            z = z + dz;
            D_A = D_A + dD_A;
            E_A = E_A + dE_A;
            q = q + dq_jc1;

            fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
                    q, z, (1+z)*(1+z)*D_A, (2/H0/a(T))*(1.0 - sqrt(a(T))), 
                    ((1+z)*(1+z)*D_A-(2/H0/a(T))*(1.0 - sqrt(a(T))))/((2/H0/a(T))*(1.0 - sqrt(a(T)))), 
                    (z-(1.0/a(T) - 1.0))/(1.0/a(T) - 1.0), 
                    (D_A-(2*a(T)/H0)*(1.0 - sqrt(a(T))))/((2*a(T)/H0)*(1.0 - sqrt(a(T)))), 
                    T, ddPhi(q), dPhi(q), T_FRZ(q), t_col(q), x(T,q));
        }
    }

    /* go back one step */
    T = T - dT;
    z = z - dz;
    D_A = D_A - dD_A;
    E_A = E_A - dE_A;
    q = q - dq_jc1;

    /* reduce step size */
    for (jc2 = 0; j < N_JC2; jc2++) {
        if (z > z_s) {
            break;
        }
  
        if (ddPhi(q) > 0 && T >= T_FRZ(q)) {

            kA1 = fA_FRZ(T,q);
            kB1 = fB_FRZ(T,q,z);
            kC1 = fC_FRZ(T,q,z,E_A);
            kD1 = fD_FRZ(T,q,z,D_A);

            kA2 = fA_FRZ(T + dq_jc2*kA1/2.0, q + dq_jc2/2.0);
            kB2 = fB_FRZ(T + dq_jc2*kA1/2.0, q + dq_jc2/2.0, z + dq_jc2*kB1/2.0);
            kC2 = fC_FRZ(T + dq_jc2*kA1/2.0, q + dq_jc2/2.0, z + dq_jc2*kB1/2.0, E_A + dq_jc2*kD1/2.0);
            kD2 = fD_FRZ(T + dq_jc2*kA1/2.0, q + dq_jc2/2.0, z + dq_jc2*kB1/2.0, D_A + dq_jc2*kC1/2.0);

            kA3 = fA_FRZ(T + dq_jc2*kA2/2.0, q + dq_jc2/2.0);
            kB3 = fB_FRZ(T + dq_jc2*kA2/2.0, q + dq_jc2/2.0, z + dq_jc2*kB2/2.0);
            kC3 = fC_FRZ(T + dq_jc2*kA2/2.0, q + dq_jc2/2.0, z + dq_jc2*kB2/2.0, E_A + dq_jc2*kD2/2.0);
            kD3 = fD_FRZ(T + dq_jc2*kA2/2.0, q + dq_jc2/2.0, z + dq_jc2*kB2/2.0, D_A + dq_jc2*kC2/2.0);

            kA4 = fA_FRZ(T + dq_jc2*kA3, q + dq_jc2);
            kB4 = fB_FRZ(T + dq_jc2*kA3, q + dq_jc2, z + dq_jc2*kB3);
            kC4 = fC_FRZ(T + dq_jc2*kA3, q + dq_jc2, z + dq_jc2*kB3, E_A + dq_jc2*kD3);
            kD4 = fD_FRZ(T + dq_jc2*kA3, q + dq_jc2, z + dq_jc2*kB3, D_A + dq_jc2*kC3);

            dT = dq_jc2*(kA1 + 2*kA2 + 2*kA3 + kA4)/6.0;
            dz = dq_jc2*(kB1 + 2*kB2 + 2*kB3 + kB4)/6.0;
            dD_A = dq_jc2*(kC1 + 2*kC2 + 2*kC3 + kC4)/6.0;
            dE_A = dq_jc2*(kD1 + 2*kD2 + 2*kD3 + kD4)/6.0;

            T = T + dT;
            z = z + dz;
            D_A = D_A + dD_A;
            E_A = E_A + dE_A;
            q = q + dq_jc2;

            fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
                    q, z, (1+z)*(1+z)*D_A, (2/H0/a(T))*(1.0 - sqrt(a(T))), 
                    ((1+z)*(1+z)*D_A-(2/H0/a(T))*(1.0 - sqrt(a(T))))/((2/H0/a(T))*(1.0 - sqrt(a(T)))), 
                    (z-(1.0/a(T) - 1.0))/(1.0/a(T) - 1.0), 
                    (D_A-(2*a(T)/H0)*(1.0 - sqrt(a(T))))/((2*a(T)/H0)*(1.0 - sqrt(a(T)))), 
                    T, ddPhi(q), dPhi(q), T_FRZ(q), t_col(q), x(T,q));
        }

        else{

            kA1 = fA(T,q);
            kB1 = fB(T,q,z);
            kC1 = fC(T,q,z,E_A);
            kD1 = fD(T,q,z,D_A);

            kA2 = fA(T + dq_jc2*kA1/2.0, q + dq_jc2/2.0);
            kB2 = fB(T + dq_jc2*kA1/2.0, q + dq_jc2/2.0, z + dq_jc2*kB1/2.0);
            kC2 = fC(T + dq_jc2*kA1/2.0, q + dq_jc2/2.0, z + dq_jc2*kB1/2.0, E_A + dq_jc2*kD1/2.0);
            kD2 = fD(T + dq_jc2*kA1/2.0, q + dq_jc2/2.0, z + dq_jc2*kB1/2.0, D_A + dq_jc2*kC1/2.0);

            kA3 = fA(T + dq_jc2*kA2/2.0, q + dq_jc2/2.0);
            kB3 = fB(T + dq_jc2*kA2/2.0, q + dq_jc2/2.0, z + dq_jc2*kB2/2.0);
            kC3 = fC(T + dq_jc2*kA2/2.0, q + dq_jc2/2.0, z + dq_jc2*kB2/2.0, E_A + dq_jc2*kD2/2.0);
            kD3 = fD(T + dq_jc2*kA2/2.0, q + dq_jc2/2.0, z + dq_jc2*kB2/2.0, D_A + dq_jc2*kC2/2.0);

            kA4 = fA(T + dq_jc2*kA3, q + dq_jc2);
            kB4 = fB(T + dq_jc2*kA3, q + dq_jc2, z + dq_jc2*kB3);
            kC4 = fC(T + dq_jc2*kA3, q + dq_jc2, z + dq_jc2*kB3, E_A + dq_jc2*kD3);
            kD4 = fD(T + dq_jc2*kA3, q + dq_jc2, z + dq_jc2*kB3, D_A + dq_jc2*kC3);

            dT = dq_jc2*(kA1 + 2*kA2 + 2*kA3 + kA4)/6.0;
            dz = dq_jc2*(kB1 + 2*kB2 + 2*kB3 + kB4)/6.0;
            dD_A = dq_jc2*(kC1 + 2*kC2 + 2*kC3 + kC4)/6.0;
            dE_A = dq_jc2*(kD1 + 2*kD2 + 2*kD3 + kD4)/6.0;

            T = T + dT;
            z = z + dz;
            D_A = D_A + dD_A;
            E_A = E_A + dE_A;
            q = q + dq_jc2;

            fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
                    q, z, (1+z)*(1+z)*D_A, (2/H0/a(T))*(1.0 - sqrt(a(T))), 
                    ((1+z)*(1+z)*D_A-(2/H0/a(T))*(1.0 - sqrt(a(T))))/((2/H0/a(T))*(1.0 - sqrt(a(T)))), 
                    (z-(1.0/a(T) - 1.0))/(1.0/a(T) - 1.0), 
                    (D_A-(2*a(T)/H0)*(1.0 - sqrt(a(T))))/((2*a(T)/H0)*(1.0 - sqrt(a(T)))), 
                    T, ddPhi(q), dPhi(q), T_FRZ(q), t_col(q), x(T,q));
        }
    }

    fclose(fp);

	/* find the Hubble constant */
    double D_L;				/* luminosity distance */
    double H0_I;			/* Hubble constant estimated from luminosity distance (Mpc^{-1}) */

    D_L = (1+z)*(1+z)*D_A;
    H0_I = (2*(1+z)/D_L)*(1.0 - 1.0/sqrt(1+z));

    printf("[PERTURBATION]\n");
	printf("  Amplitude A  	    = %e\n", A);
	printf("  phi(t=0,q=0) 	    = %e,\n", phi_in(0));
	printf("  Width L           = %e\n",L);
    printf("\n");

    printf("[FINAL REDSHIFT]\n");
	printf("  z     = %e\n",z);
    printf("\n");

    printf("[COORDINATES]\n");
    printf("  FROM OBSERVER q     = %e\n", q0);
	printf("  TO SOURCE     q     = %e\n", q);
    printf("\n");

    printf("[HUBBLE CONSTANT]\n");
	printf("  H0_I               	= %f km/s/Mpc\n",H0_I*(C_SI*1.0e-3));
    printf("  (H0_I - H0_B)/H0_B 	= %e\n",(H0_I-H0)/H0);

    return 0;
}

