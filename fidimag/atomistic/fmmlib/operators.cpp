#include "operators.h"
#include "math.h"
void P2M_2(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = (1.0/2.0)*q;
M[0] += -Mtmp0;
M[1] += -Mtmp1;
M[2] += -q*z;
M[3] += Mtmp2*pow(x, 2);
M[4] += Mtmp0*y;
M[5] += Mtmp0*z;
M[6] += Mtmp2*pow(y, 2);
M[7] += Mtmp1*z;
M[8] += Mtmp2*pow(z, 2);

}
void M2M_2(double x, double y, double z, double * M, double * Ms) {
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += x*M[0] + M[3];
#pragma omp atomic
Ms[4] += x*M[1] + y*M[0] + M[4];
#pragma omp atomic
Ms[5] += x*M[2] + z*M[0] + M[5];
#pragma omp atomic
Ms[6] += y*M[1] + M[6];
#pragma omp atomic
Ms[7] += y*M[2] + z*M[1] + M[7];
#pragma omp atomic
Ms[8] += z*M[2] + M[8];

}

void M2L_2(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[9];
double Dtmp0 = pow(R, -3);
double Dtmp1 = 1.0*Dtmp0;
double Dtmp2 = 3.0/pow(R, 2);
double Dtmp3 = pow(R, -5);
double Dtmp4 = 3.0*Dtmp3*x;
D[0] = -Dtmp1*x;
D[1] = -Dtmp1*y;
D[2] = -Dtmp1*z;
D[3] = Dtmp0*(Dtmp2*pow(x, 2) - 1.0);
D[4] = Dtmp4*y;
D[5] = Dtmp4*z;
D[6] = Dtmp0*(Dtmp2*pow(y, 2) - 1.0);
D[7] = 3.0*Dtmp3*y*z;
D[8] = -D[3] - D[6];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8];
#pragma omp atomic
L[1] += D[3]*M[0] + D[4]*M[1] + D[5]*M[2];
#pragma omp atomic
L[2] += D[4]*M[0] + D[6]*M[1] + D[7]*M[2];
#pragma omp atomic
L[3] += D[5]*M[0] + D[7]*M[1] + D[8]*M[2];

}

void L2L_2(double x, double y, double z, double * L, double * Ls) {
#pragma omp atomic
Ls[0] += x*L[1] + y*L[2] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += L[1];
#pragma omp atomic
Ls[2] += L[2];
#pragma omp atomic
Ls[3] += L[3];

}

void L2P_2(double x, double y, double z, double * L, double * F) {
#pragma omp atomic
F[0] += -L[1];
#pragma omp atomic
F[1] += -L[2];
#pragma omp atomic
F[2] += -L[3];

}

void M2P_2(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(R, -3);
double Ftmp1 = pow(R, -2);
double Ftmp2 = 3.0*Ftmp1;
double Ftmp3 = y*M[4];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = Ftmp2*x;
double Ftmp6 = y*M[1];
double Ftmp7 = Ftmp4*M[2];
double Ftmp8 = pow(x, 2);
double Ftmp9 = Ftmp2*M[0];
double Ftmp10 = y*M[7];
double Ftmp11 = 15.0/pow(R, 4);
double Ftmp12 = Ftmp11*x*z;
double Ftmp13 = Ftmp11*Ftmp8;
double Ftmp14 = z*M[5];
double Ftmp15 = Ftmp1*x;
double Ftmp16 = pow(y, 2);
double Ftmp17 = 15.0*Ftmp1;
double Ftmp18 = -Ftmp16*Ftmp17;
double Ftmp19 = (Ftmp18 + 3.0)*M[6];
double Ftmp20 = pow(z, 2);
double Ftmp21 = -Ftmp17*Ftmp20;
double Ftmp22 = (Ftmp21 + 3.0)*M[8];
double Ftmp23 = -Ftmp17*Ftmp8;
double Ftmp24 = x*y;
double Ftmp25 = Ftmp11*Ftmp16;
double Ftmp26 = Ftmp1*y;
double Ftmp27 = (Ftmp23 + 3.0)*M[3];
double Ftmp28 = Ftmp11*Ftmp20;
double Ftmp29 = Ftmp1*z;
#pragma omp atomic
F[0] += Ftmp0*(Ftmp10*Ftmp12 + Ftmp13*Ftmp14 + Ftmp13*Ftmp3 - Ftmp15*Ftmp19 - Ftmp15*Ftmp22 - Ftmp15*(Ftmp23 + 9.0)*M[3] - Ftmp2*Ftmp3 - Ftmp4*M[5] - Ftmp5*Ftmp6 - Ftmp7*x - Ftmp8*Ftmp9 + 1.0*M[0]);
#pragma omp atomic
F[1] += Ftmp0*(Ftmp11*Ftmp14*Ftmp24 - Ftmp16*Ftmp2*M[1] - Ftmp22*Ftmp26 - Ftmp24*Ftmp9 + Ftmp25*x*M[4] + Ftmp25*z*M[7] - Ftmp26*Ftmp27 - Ftmp26*(Ftmp18 + 9.0)*M[6] - Ftmp4*M[7] - Ftmp5*M[4] - Ftmp7*y + 1.0*M[1]);
#pragma omp atomic
F[2] += Ftmp0*(-Ftmp10*Ftmp2 + Ftmp10*Ftmp28 + Ftmp12*Ftmp3 - Ftmp19*Ftmp29 - Ftmp2*Ftmp20*M[2] - Ftmp27*Ftmp29 + Ftmp28*x*M[5] - Ftmp29*(Ftmp21 + 9.0)*M[8] - Ftmp4*Ftmp6 - Ftmp4*x*M[0] - Ftmp5*M[5] + 1.0*M[2]);

}

void P2P(double x, double y, double z, double * S, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(R, -3);
double Ftmp1 = 3.0/pow(R, 2);
double Ftmp2 = Ftmp1*S[1];
double Ftmp3 = x*y;
double Ftmp4 = Ftmp1*S[2];
double Ftmp5 = Ftmp4*z;
double Ftmp6 = Ftmp1*S[0];
#pragma omp atomic
F[0] += Ftmp0*(-Ftmp2*Ftmp3 - Ftmp5*x - Ftmp6*pow(x, 2) + 1.0*S[0]);
#pragma omp atomic
F[1] += Ftmp0*(-Ftmp2*pow(y, 2) - Ftmp3*Ftmp6 - Ftmp5*y + 1.0*S[1]);
#pragma omp atomic
F[2] += Ftmp0*(-Ftmp2*y*z - Ftmp4*pow(z, 2) - Ftmp6*x*z + 1.0*S[2]);

}

void P2M_3(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = q*z;
double Mtmp3 = pow(x, 2);
double Mtmp4 = (1.0/2.0)*q;
double Mtmp5 = Mtmp0*y;
double Mtmp6 = pow(y, 2);
double Mtmp7 = pow(z, 2);
double Mtmp8 = (1.0/6.0)*q;
double Mtmp9 = (1.0/2.0)*Mtmp3;
double Mtmp10 = (1.0/2.0)*Mtmp0;
M[0] += -Mtmp0;
M[1] += -Mtmp1;
M[2] += -Mtmp2;
M[3] += Mtmp3*Mtmp4;
M[4] += Mtmp5;
M[5] += Mtmp0*z;
M[6] += Mtmp4*Mtmp6;
M[7] += Mtmp1*z;
M[8] += Mtmp4*Mtmp7;
M[9] += -Mtmp8*pow(x, 3);
M[10] += -Mtmp1*Mtmp9;
M[11] += -Mtmp2*Mtmp9;
M[12] += -Mtmp10*Mtmp6;
M[13] += -Mtmp5*z;
M[14] += -Mtmp10*Mtmp7;
M[15] += -Mtmp8*pow(y, 3);
M[16] += -1.0/2.0*Mtmp2*Mtmp6;
M[17] += -1.0/2.0*Mtmp1*Mtmp7;
M[18] += -Mtmp8*pow(z, 3);

}
void M2M_3(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = x*M[1];
double Mstmp2 = y*M[0];
double Mstmp3 = x*M[2];
double Mstmp4 = y*M[1];
double Mstmp5 = y*M[2];
double Mstmp6 = (1.0/2.0)*pow(x, 2);
double Mstmp7 = pow(y, 2);
double Mstmp8 = (1.0/2.0)*M[0];
double Mstmp9 = pow(z, 2);
double Mstmp10 = (1.0/2.0)*Mstmp7;
double Mstmp11 = (1.0/2.0)*Mstmp9;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += Mstmp0 + M[3];
#pragma omp atomic
Ms[4] += Mstmp1 + Mstmp2 + M[4];
#pragma omp atomic
Ms[5] += Mstmp3 + z*M[0] + M[5];
#pragma omp atomic
Ms[6] += Mstmp4 + M[6];
#pragma omp atomic
Ms[7] += Mstmp5 + z*M[1] + M[7];
#pragma omp atomic
Ms[8] += z*M[2] + M[8];
#pragma omp atomic
Ms[9] += Mstmp6*M[0] + x*M[3] + M[9];
#pragma omp atomic
Ms[10] += Mstmp0*y + Mstmp6*M[1] + x*M[4] + y*M[3] + M[10];
#pragma omp atomic
Ms[11] += Mstmp0*z + Mstmp6*M[2] + x*M[5] + z*M[3] + M[11];
#pragma omp atomic
Ms[12] += Mstmp1*y + Mstmp7*Mstmp8 + x*M[6] + y*M[4] + M[12];
#pragma omp atomic
Ms[13] += Mstmp1*z + Mstmp2*z + Mstmp3*y + x*M[7] + y*M[5] + z*M[4] + M[13];
#pragma omp atomic
Ms[14] += Mstmp3*z + Mstmp8*Mstmp9 + x*M[8] + z*M[5] + M[14];
#pragma omp atomic
Ms[15] += Mstmp10*M[1] + y*M[6] + M[15];
#pragma omp atomic
Ms[16] += Mstmp10*M[2] + Mstmp4*z + y*M[7] + z*M[6] + M[16];
#pragma omp atomic
Ms[17] += Mstmp11*M[1] + Mstmp5*z + y*M[8] + z*M[7] + M[17];
#pragma omp atomic
Ms[18] += Mstmp11*M[2] + z*M[8] + M[18];

}

void M2L_3(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[19];
double Dtmp0 = pow(R, -3);
double Dtmp1 = 1.0*Dtmp0;
double Dtmp2 = pow(x, 2);
double Dtmp3 = pow(R, -2);
double Dtmp4 = 3.0*Dtmp3;
double Dtmp5 = pow(R, -5);
double Dtmp6 = Dtmp5*x;
double Dtmp7 = 3.0*Dtmp6;
double Dtmp8 = pow(y, 2);
double Dtmp9 = Dtmp5*y;
double Dtmp10 = 15.0*Dtmp3;
double Dtmp11 = -Dtmp10*Dtmp2;
double Dtmp12 = Dtmp5*(Dtmp11 + 3.0);
double Dtmp13 = -Dtmp10*Dtmp8;
double Dtmp14 = Dtmp13 + 3.0;
D[0] = -Dtmp1*x;
D[1] = -Dtmp1*y;
D[2] = -Dtmp1*z;
D[3] = Dtmp0*(Dtmp2*Dtmp4 - 1.0);
D[4] = Dtmp7*y;
D[5] = Dtmp7*z;
D[6] = Dtmp0*(Dtmp4*Dtmp8 - 1.0);
D[7] = 3.0*Dtmp9*z;
D[8] = -D[3] - D[6];
D[9] = Dtmp6*(Dtmp11 + 9.0);
D[10] = Dtmp12*y;
D[11] = Dtmp12*z;
D[12] = 1.0*Dtmp14*Dtmp6;
D[13] = -15.0*x*y*z/pow(R, 7);
D[14] = -D[9] - D[12];
D[15] = Dtmp9*(Dtmp13 + 9.0);
D[16] = Dtmp14*Dtmp5*z;
D[17] = -D[10] - D[15];
D[18] = -D[11] - D[16];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8] + D[9]*M[9] + D[10]*M[10] + D[11]*M[11] + D[12]*M[12] + D[13]*M[13] + D[14]*M[14] + D[15]*M[15] + D[16]*M[16] + D[17]*M[17] + D[18]*M[18];
#pragma omp atomic
L[1] += D[3]*M[0] + D[4]*M[1] + D[5]*M[2] + D[9]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[14]*M[8];
#pragma omp atomic
L[2] += D[4]*M[0] + D[6]*M[1] + D[7]*M[2] + D[10]*M[3] + D[12]*M[4] + D[13]*M[5] + D[15]*M[6] + D[16]*M[7] + D[17]*M[8];
#pragma omp atomic
L[3] += D[5]*M[0] + D[7]*M[1] + D[8]*M[2] + D[11]*M[3] + D[13]*M[4] + D[14]*M[5] + D[16]*M[6] + D[17]*M[7] + D[18]*M[8];
#pragma omp atomic
L[4] += D[9]*M[0] + D[10]*M[1] + D[11]*M[2];
#pragma omp atomic
L[5] += D[10]*M[0] + D[12]*M[1] + D[13]*M[2];
#pragma omp atomic
L[6] += D[11]*M[0] + D[13]*M[1] + D[14]*M[2];
#pragma omp atomic
L[7] += D[12]*M[0] + D[15]*M[1] + D[16]*M[2];
#pragma omp atomic
L[8] += D[13]*M[0] + D[16]*M[1] + D[17]*M[2];
#pragma omp atomic
L[9] += D[14]*M[0] + D[17]*M[1] + D[18]*M[2];

}

void L2L_3(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp2*y + (1.0/2.0)*pow(x, 2)*L[4] + x*L[1] + (1.0/2.0)*pow(y, 2)*L[7] + y*L[2] + (1.0/2.0)*pow(z, 2)*L[9] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp2 + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += x*L[6] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += L[4];
#pragma omp atomic
Ls[5] += L[5];
#pragma omp atomic
Ls[6] += L[6];
#pragma omp atomic
Ls[7] += L[7];
#pragma omp atomic
Ls[8] += L[8];
#pragma omp atomic
Ls[9] += L[9];

}

void L2P_3(double x, double y, double z, double * L, double * F) {
#pragma omp atomic
F[0] += -x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -x*L[6] - y*L[8] - z*L[9] - L[3];

}

void M2P_3(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(R, -3);
double Ftmp1 = pow(R, -2);
double Ftmp2 = 3.0*Ftmp1;
double Ftmp3 = y*M[4];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = pow(R, -4);
double Ftmp6 = Ftmp5*y;
double Ftmp7 = 15.0*z;
double Ftmp8 = Ftmp7*M[13];
double Ftmp9 = Ftmp2*x;
double Ftmp10 = y*M[1];
double Ftmp11 = Ftmp4*M[2];
double Ftmp12 = pow(x, 2);
double Ftmp13 = Ftmp2*M[0];
double Ftmp14 = y*M[7];
double Ftmp15 = Ftmp5*x;
double Ftmp16 = Ftmp15*Ftmp7;
double Ftmp17 = Ftmp12*Ftmp5;
double Ftmp18 = Ftmp7*M[5];
double Ftmp19 = 105.0*M[13]/pow(R, 6);
double Ftmp20 = Ftmp19*z;
double Ftmp21 = pow(y, 2);
double Ftmp22 = 15.0*Ftmp1;
double Ftmp23 = -Ftmp21*Ftmp22;
double Ftmp24 = Ftmp1*(Ftmp23 + 3.0);
double Ftmp25 = pow(z, 2);
double Ftmp26 = -Ftmp22*Ftmp25;
double Ftmp27 = Ftmp1*(Ftmp26 + 3.0);
double Ftmp28 = -Ftmp12*Ftmp22;
double Ftmp29 = Ftmp1*(Ftmp28 + 9.0);
double Ftmp30 = Ftmp24*M[6];
double Ftmp31 = Ftmp27*M[8];
double Ftmp32 = 105.0*Ftmp1;
double Ftmp33 = -Ftmp12*Ftmp32;
double Ftmp34 = Ftmp33 + 45.0;
double Ftmp35 = Ftmp6*x;
double Ftmp36 = Ftmp34*Ftmp35;
double Ftmp37 = -Ftmp21*Ftmp32;
double Ftmp38 = Ftmp37 + 45.0;
double Ftmp39 = Ftmp38*M[15];
double Ftmp40 = Ftmp37 + 15.0;
double Ftmp41 = Ftmp40*M[16];
double Ftmp42 = Ftmp15*z;
double Ftmp43 = Ftmp34*Ftmp42;
double Ftmp44 = -Ftmp25*Ftmp32;
double Ftmp45 = Ftmp44 + 45.0;
double Ftmp46 = Ftmp45*M[18];
double Ftmp47 = Ftmp44 + 15.0;
double Ftmp48 = 1.0*M[17];
double Ftmp49 = Ftmp47*Ftmp48;
double Ftmp50 = 1.0*Ftmp17;
double Ftmp51 = Ftmp40*M[12];
double Ftmp52 = Ftmp47*M[14];
double Ftmp53 = x*y;
double Ftmp54 = Ftmp21*Ftmp5;
double Ftmp55 = 15.0*x;
double Ftmp56 = Ftmp1*(Ftmp28 + 3.0);
double Ftmp57 = Ftmp1*(Ftmp23 + 9.0);
double Ftmp58 = Ftmp56*M[3];
double Ftmp59 = Ftmp33 + 15.0;
double Ftmp60 = Ftmp59*M[11];
double Ftmp61 = Ftmp6*z;
double Ftmp62 = Ftmp59*M[10];
double Ftmp63 = 1.0*Ftmp35;
double Ftmp64 = Ftmp25*Ftmp5;
double Ftmp65 = Ftmp1*(Ftmp26 + 9.0);
double Ftmp66 = 1.0*Ftmp42;
#pragma omp atomic
F[0] += Ftmp0*(-Ftmp10*Ftmp9 - Ftmp11*x - Ftmp12*Ftmp13 - Ftmp12*Ftmp20*y + Ftmp14*Ftmp16 + Ftmp17*Ftmp18 + 15.0*Ftmp17*Ftmp3 + Ftmp17*(Ftmp33 + 75.0)*M[9] - Ftmp2*Ftmp3 - Ftmp24*M[12] - Ftmp27*M[14] - Ftmp29*x*M[3] - Ftmp29*M[9] - Ftmp30*x - Ftmp31*x + Ftmp35*Ftmp39 + Ftmp35*Ftmp49 + Ftmp36*M[10] - Ftmp4*M[5] + Ftmp41*Ftmp42 + Ftmp42*Ftmp46 + Ftmp43*M[11] + Ftmp50*Ftmp51 + Ftmp50*Ftmp52 + Ftmp6*Ftmp8 + 1.0*M[0]);
#pragma omp atomic
F[1] += Ftmp0*(-Ftmp11*y - Ftmp13*Ftmp53 + Ftmp15*Ftmp8 + Ftmp18*Ftmp35 - Ftmp2*Ftmp21*M[1] - Ftmp20*Ftmp21*x - Ftmp27*M[17] - Ftmp31*y + Ftmp36*M[9] + Ftmp38*Ftmp61*M[16] + Ftmp38*Ftmp63*M[12] - Ftmp4*M[7] + Ftmp46*Ftmp61 + Ftmp49*Ftmp54 + Ftmp52*Ftmp63 + Ftmp54*Ftmp55*M[4] + Ftmp54*Ftmp62 + Ftmp54*Ftmp7*M[7] + Ftmp54*(Ftmp37 + 75.0)*M[15] - Ftmp56*M[10] - Ftmp57*y*M[6] - Ftmp57*M[15] - Ftmp58*y + Ftmp60*Ftmp61 - Ftmp9*M[4] + 1.0*M[1]);
#pragma omp atomic
F[2] += Ftmp0*(-Ftmp10*Ftmp4 - Ftmp14*Ftmp2 + 15.0*Ftmp14*Ftmp64 + Ftmp16*Ftmp3 - Ftmp19*Ftmp25*Ftmp53 - Ftmp2*Ftmp25*M[2] - Ftmp24*M[16] - Ftmp30*z + 15.0*Ftmp35*M[13] + Ftmp39*Ftmp61 - Ftmp4*x*M[0] + Ftmp41*Ftmp64 + Ftmp43*M[9] + Ftmp45*Ftmp48*Ftmp61 + Ftmp45*Ftmp66*M[14] + Ftmp51*Ftmp66 + Ftmp55*Ftmp64*M[5] - Ftmp56*M[11] - Ftmp58*z + Ftmp60*Ftmp64 + Ftmp61*Ftmp62 + Ftmp64*(Ftmp44 + 75.0)*M[18] - Ftmp65*z*M[8] - Ftmp65*M[18] - Ftmp9*M[5] + 1.0*M[2]);

}

void P2M_4(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = q*z;
double Mtmp3 = pow(x, 2);
double Mtmp4 = (1.0/2.0)*q;
double Mtmp5 = Mtmp0*y;
double Mtmp6 = Mtmp0*z;
double Mtmp7 = pow(y, 2);
double Mtmp8 = Mtmp1*z;
double Mtmp9 = pow(z, 2);
double Mtmp10 = pow(x, 3);
double Mtmp11 = (1.0/6.0)*q;
double Mtmp12 = (1.0/2.0)*Mtmp3;
double Mtmp13 = (1.0/2.0)*Mtmp0;
double Mtmp14 = pow(y, 3);
double Mtmp15 = (1.0/2.0)*Mtmp7;
double Mtmp16 = (1.0/2.0)*Mtmp9;
double Mtmp17 = pow(z, 3);
double Mtmp18 = (1.0/24.0)*q;
double Mtmp19 = (1.0/6.0)*Mtmp10;
double Mtmp20 = (1.0/4.0)*Mtmp3*q;
double Mtmp21 = (1.0/6.0)*Mtmp0;
M[0] += -Mtmp0;
M[1] += -Mtmp1;
M[2] += -Mtmp2;
M[3] += Mtmp3*Mtmp4;
M[4] += Mtmp5;
M[5] += Mtmp6;
M[6] += Mtmp4*Mtmp7;
M[7] += Mtmp8;
M[8] += Mtmp4*Mtmp9;
M[9] += -Mtmp10*Mtmp11;
M[10] += -Mtmp1*Mtmp12;
M[11] += -Mtmp12*Mtmp2;
M[12] += -Mtmp13*Mtmp7;
M[13] += -Mtmp5*z;
M[14] += -Mtmp13*Mtmp9;
M[15] += -Mtmp11*Mtmp14;
M[16] += -Mtmp15*Mtmp2;
M[17] += -Mtmp1*Mtmp16;
M[18] += -Mtmp11*Mtmp17;
M[19] += Mtmp18*pow(x, 4);
M[20] += Mtmp1*Mtmp19;
M[21] += Mtmp19*Mtmp2;
M[22] += Mtmp20*Mtmp7;
M[23] += Mtmp12*Mtmp8;
M[24] += Mtmp20*Mtmp9;
M[25] += Mtmp14*Mtmp21;
M[26] += Mtmp15*Mtmp6;
M[27] += Mtmp16*Mtmp5;
M[28] += Mtmp17*Mtmp21;
M[29] += Mtmp18*pow(y, 4);
M[30] += (1.0/6.0)*Mtmp14*Mtmp2;
M[31] += (1.0/4.0)*Mtmp7*Mtmp9*q;
M[32] += (1.0/6.0)*Mtmp1*Mtmp17;
M[33] += Mtmp18*pow(z, 4);

}
void M2M_4(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = x*M[1];
double Mstmp2 = y*M[0];
double Mstmp3 = x*M[2];
double Mstmp4 = z*M[0];
double Mstmp5 = y*M[1];
double Mstmp6 = y*M[2];
double Mstmp7 = z*M[1];
double Mstmp8 = z*M[2];
double Mstmp9 = x*M[3];
double Mstmp10 = (1.0/2.0)*pow(x, 2);
double Mstmp11 = x*M[4];
double Mstmp12 = y*M[3];
double Mstmp13 = Mstmp0*y;
double Mstmp14 = x*M[5];
double Mstmp15 = x*M[6];
double Mstmp16 = y*M[4];
double Mstmp17 = Mstmp1*y;
double Mstmp18 = pow(y, 2);
double Mstmp19 = (1.0/2.0)*M[0];
double Mstmp20 = x*M[7];
double Mstmp21 = y*M[5];
double Mstmp22 = Mstmp3*y;
double Mstmp23 = x*M[8];
double Mstmp24 = pow(z, 2);
double Mstmp25 = y*M[6];
double Mstmp26 = (1.0/2.0)*Mstmp18;
double Mstmp27 = y*M[7];
double Mstmp28 = y*M[8];
double Mstmp29 = (1.0/2.0)*Mstmp24;
double Mstmp30 = (1.0/6.0)*pow(x, 3);
double Mstmp31 = pow(y, 3);
double Mstmp32 = (1.0/6.0)*M[0];
double Mstmp33 = pow(z, 3);
double Mstmp34 = (1.0/6.0)*Mstmp31;
double Mstmp35 = (1.0/6.0)*Mstmp33;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += Mstmp0 + M[3];
#pragma omp atomic
Ms[4] += Mstmp1 + Mstmp2 + M[4];
#pragma omp atomic
Ms[5] += Mstmp3 + Mstmp4 + M[5];
#pragma omp atomic
Ms[6] += Mstmp5 + M[6];
#pragma omp atomic
Ms[7] += Mstmp6 + Mstmp7 + M[7];
#pragma omp atomic
Ms[8] += Mstmp8 + M[8];
#pragma omp atomic
Ms[9] += Mstmp10*M[0] + Mstmp9 + M[9];
#pragma omp atomic
Ms[10] += Mstmp10*M[1] + Mstmp11 + Mstmp12 + Mstmp13 + M[10];
#pragma omp atomic
Ms[11] += Mstmp0*z + Mstmp10*M[2] + Mstmp14 + z*M[3] + M[11];
#pragma omp atomic
Ms[12] += Mstmp15 + Mstmp16 + Mstmp17 + Mstmp18*Mstmp19 + M[12];
#pragma omp atomic
Ms[13] += Mstmp1*z + Mstmp2*z + Mstmp20 + Mstmp21 + Mstmp22 + z*M[4] + M[13];
#pragma omp atomic
Ms[14] += Mstmp19*Mstmp24 + Mstmp23 + Mstmp3*z + z*M[5] + M[14];
#pragma omp atomic
Ms[15] += Mstmp25 + Mstmp26*M[1] + M[15];
#pragma omp atomic
Ms[16] += Mstmp26*M[2] + Mstmp27 + Mstmp5*z + z*M[6] + M[16];
#pragma omp atomic
Ms[17] += Mstmp28 + Mstmp29*M[1] + Mstmp6*z + z*M[7] + M[17];
#pragma omp atomic
Ms[18] += Mstmp29*M[2] + z*M[8] + M[18];
#pragma omp atomic
Ms[19] += Mstmp10*M[3] + Mstmp30*M[0] + x*M[9] + M[19];
#pragma omp atomic
Ms[20] += Mstmp10*Mstmp2 + Mstmp10*M[4] + Mstmp30*M[1] + Mstmp9*y + x*M[10] + y*M[9] + M[20];
#pragma omp atomic
Ms[21] += Mstmp10*Mstmp4 + Mstmp10*M[5] + Mstmp30*M[2] + Mstmp9*z + x*M[11] + z*M[9] + M[21];
#pragma omp atomic
Ms[22] += Mstmp0*Mstmp26 + Mstmp10*Mstmp5 + Mstmp10*M[6] + Mstmp11*y + Mstmp26*M[3] + x*M[12] + y*M[10] + M[22];
#pragma omp atomic
Ms[23] += Mstmp10*Mstmp6 + Mstmp10*Mstmp7 + Mstmp10*M[7] + Mstmp11*z + Mstmp12*z + Mstmp13*z + Mstmp14*y + x*M[13] + y*M[11] + z*M[10] + M[23];
#pragma omp atomic
Ms[24] += Mstmp0*Mstmp29 + Mstmp10*Mstmp8 + Mstmp10*M[8] + Mstmp14*z + Mstmp29*M[3] + x*M[14] + z*M[11] + M[24];
#pragma omp atomic
Ms[25] += Mstmp1*Mstmp26 + Mstmp15*y + Mstmp26*M[4] + Mstmp31*Mstmp32 + x*M[15] + y*M[12] + M[25];
#pragma omp atomic
Ms[26] += Mstmp15*z + Mstmp16*z + Mstmp17*z + Mstmp20*y + Mstmp26*Mstmp3 + Mstmp26*Mstmp4 + Mstmp26*M[5] + x*M[16] + y*M[13] + z*M[12] + M[26];
#pragma omp atomic
Ms[27] += Mstmp1*Mstmp29 + Mstmp2*Mstmp29 + Mstmp20*z + Mstmp21*z + Mstmp22*z + Mstmp23*y + Mstmp29*M[4] + x*M[17] + y*M[14] + z*M[13] + M[27];
#pragma omp atomic
Ms[28] += Mstmp23*z + Mstmp29*Mstmp3 + Mstmp29*M[5] + Mstmp32*Mstmp33 + x*M[18] + z*M[14] + M[28];
#pragma omp atomic
Ms[29] += Mstmp26*M[6] + Mstmp34*M[1] + y*M[15] + M[29];
#pragma omp atomic
Ms[30] += Mstmp25*z + Mstmp26*Mstmp7 + Mstmp26*M[7] + Mstmp34*M[2] + y*M[16] + z*M[15] + M[30];
#pragma omp atomic
Ms[31] += Mstmp26*Mstmp8 + Mstmp26*M[8] + Mstmp27*z + Mstmp29*Mstmp5 + Mstmp29*M[6] + y*M[17] + z*M[16] + M[31];
#pragma omp atomic
Ms[32] += Mstmp28*z + Mstmp29*Mstmp6 + Mstmp29*M[7] + Mstmp35*M[1] + y*M[18] + z*M[17] + M[32];
#pragma omp atomic
Ms[33] += Mstmp29*M[8] + Mstmp35*M[2] + z*M[18] + M[33];

}

void M2L_4(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[34];
double Dtmp0 = pow(R, -3);
double Dtmp1 = 1.0*Dtmp0;
double Dtmp2 = pow(x, 2);
double Dtmp3 = pow(R, -2);
double Dtmp4 = 3.0*Dtmp3;
double Dtmp5 = pow(R, -5);
double Dtmp6 = Dtmp5*x;
double Dtmp7 = 3.0*Dtmp6;
double Dtmp8 = pow(y, 2);
double Dtmp9 = Dtmp5*y;
double Dtmp10 = 15.0*Dtmp3;
double Dtmp11 = -Dtmp10*Dtmp2;
double Dtmp12 = Dtmp11 + 3.0;
double Dtmp13 = Dtmp12*Dtmp5;
double Dtmp14 = -Dtmp10*Dtmp8;
double Dtmp15 = Dtmp14 + 3.0;
double Dtmp16 = pow(R, -7);
double Dtmp17 = Dtmp16*y;
double Dtmp18 = Dtmp17*z;
double Dtmp19 = 105.0/pow(R, 4);
double Dtmp20 = Dtmp2*Dtmp3;
double Dtmp21 = -105.0*Dtmp20;
double Dtmp22 = x*(Dtmp21 + 45.0);
double Dtmp23 = Dtmp16*z;
double Dtmp24 = Dtmp3*Dtmp8;
double Dtmp25 = -105.0*Dtmp24;
double Dtmp26 = Dtmp25 + 45.0;
double Dtmp27 = 1.0*x;
D[0] = -Dtmp1*x;
D[1] = -Dtmp1*y;
D[2] = -Dtmp1*z;
D[3] = Dtmp0*(Dtmp2*Dtmp4 - 1.0);
D[4] = Dtmp7*y;
D[5] = Dtmp7*z;
D[6] = Dtmp0*(Dtmp4*Dtmp8 - 1.0);
D[7] = 3.0*Dtmp9*z;
D[8] = -D[3] - D[6];
D[9] = Dtmp6*(Dtmp11 + 9.0);
D[10] = Dtmp13*y;
D[11] = Dtmp13*z;
D[12] = 1.0*Dtmp15*Dtmp6;
D[13] = -15.0*Dtmp18*x;
D[14] = -D[9] - D[12];
D[15] = Dtmp9*(Dtmp14 + 9.0);
D[16] = Dtmp15*Dtmp5*z;
D[17] = -D[10] - D[15];
D[18] = -D[11] - D[16];
D[19] = Dtmp5*(Dtmp19*pow(x, 4) - 90.0*Dtmp20 + 9.0);
D[20] = -Dtmp17*Dtmp22;
D[21] = -Dtmp22*Dtmp23;
D[22] = Dtmp5*(Dtmp12 + Dtmp14 + Dtmp19*Dtmp2*Dtmp8);
D[23] = -Dtmp18*(Dtmp21 + 15.0);
D[24] = -D[19] - D[22];
D[25] = -Dtmp17*Dtmp26*Dtmp27;
D[26] = -Dtmp23*Dtmp27*(Dtmp25 + 15.0);
D[27] = -D[20] - D[25];
D[28] = -D[21] - D[26];
D[29] = Dtmp5*(Dtmp19*pow(y, 4) - 90.0*Dtmp24 + 9.0);
D[30] = -Dtmp18*Dtmp26;
D[31] = -D[22] - D[29];
D[32] = -D[23] - D[30];
D[33] = -D[24] - D[31];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8] + D[9]*M[9] + D[10]*M[10] + D[11]*M[11] + D[12]*M[12] + D[13]*M[13] + D[14]*M[14] + D[15]*M[15] + D[16]*M[16] + D[17]*M[17] + D[18]*M[18] + D[19]*M[19] + D[20]*M[20] + D[21]*M[21] + D[22]*M[22] + D[23]*M[23] + D[24]*M[24] + D[25]*M[25] + D[26]*M[26] + D[27]*M[27] + D[28]*M[28] + D[29]*M[29] + D[30]*M[30] + D[31]*M[31] + D[32]*M[32] + D[33]*M[33];
#pragma omp atomic
L[1] += D[3]*M[0] + D[4]*M[1] + D[5]*M[2] + D[9]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[14]*M[8] + D[19]*M[9] + D[20]*M[10] + D[21]*M[11] + D[22]*M[12] + D[23]*M[13] + D[24]*M[14] + D[25]*M[15] + D[26]*M[16] + D[27]*M[17] + D[28]*M[18];
#pragma omp atomic
L[2] += D[4]*M[0] + D[6]*M[1] + D[7]*M[2] + D[10]*M[3] + D[12]*M[4] + D[13]*M[5] + D[15]*M[6] + D[16]*M[7] + D[17]*M[8] + D[20]*M[9] + D[22]*M[10] + D[23]*M[11] + D[25]*M[12] + D[26]*M[13] + D[27]*M[14] + D[29]*M[15] + D[30]*M[16] + D[31]*M[17] + D[32]*M[18];
#pragma omp atomic
L[3] += D[5]*M[0] + D[7]*M[1] + D[8]*M[2] + D[11]*M[3] + D[13]*M[4] + D[14]*M[5] + D[16]*M[6] + D[17]*M[7] + D[18]*M[8] + D[21]*M[9] + D[23]*M[10] + D[24]*M[11] + D[26]*M[12] + D[27]*M[13] + D[28]*M[14] + D[30]*M[15] + D[31]*M[16] + D[32]*M[17] + D[33]*M[18];
#pragma omp atomic
L[4] += D[9]*M[0] + D[10]*M[1] + D[11]*M[2] + D[19]*M[3] + D[20]*M[4] + D[21]*M[5] + D[22]*M[6] + D[23]*M[7] + D[24]*M[8];
#pragma omp atomic
L[5] += D[10]*M[0] + D[12]*M[1] + D[13]*M[2] + D[20]*M[3] + D[22]*M[4] + D[23]*M[5] + D[25]*M[6] + D[26]*M[7] + D[27]*M[8];
#pragma omp atomic
L[6] += D[11]*M[0] + D[13]*M[1] + D[14]*M[2] + D[21]*M[3] + D[23]*M[4] + D[24]*M[5] + D[26]*M[6] + D[27]*M[7] + D[28]*M[8];
#pragma omp atomic
L[7] += D[12]*M[0] + D[15]*M[1] + D[16]*M[2] + D[22]*M[3] + D[25]*M[4] + D[26]*M[5] + D[29]*M[6] + D[30]*M[7] + D[31]*M[8];
#pragma omp atomic
L[8] += D[13]*M[0] + D[16]*M[1] + D[17]*M[2] + D[23]*M[3] + D[26]*M[4] + D[27]*M[5] + D[30]*M[6] + D[31]*M[7] + D[32]*M[8];
#pragma omp atomic
L[9] += D[14]*M[0] + D[17]*M[1] + D[18]*M[2] + D[24]*M[3] + D[27]*M[4] + D[28]*M[5] + D[31]*M[6] + D[32]*M[7] + D[33]*M[8];
#pragma omp atomic
L[10] += D[19]*M[0] + D[20]*M[1] + D[21]*M[2];
#pragma omp atomic
L[11] += D[20]*M[0] + D[22]*M[1] + D[23]*M[2];
#pragma omp atomic
L[12] += D[21]*M[0] + D[23]*M[1] + D[24]*M[2];
#pragma omp atomic
L[13] += D[22]*M[0] + D[25]*M[1] + D[26]*M[2];
#pragma omp atomic
L[14] += D[23]*M[0] + D[26]*M[1] + D[27]*M[2];
#pragma omp atomic
L[15] += D[24]*M[0] + D[27]*M[1] + D[28]*M[2];
#pragma omp atomic
L[16] += D[25]*M[0] + D[29]*M[1] + D[30]*M[2];
#pragma omp atomic
L[17] += D[26]*M[0] + D[30]*M[1] + D[31]*M[2];
#pragma omp atomic
L[18] += D[27]*M[0] + D[31]*M[1] + D[32]*M[2];
#pragma omp atomic
L[19] += D[28]*M[0] + D[32]*M[1] + D[33]*M[2];

}

void L2L_4(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[14];
double Lstmp4 = Lstmp3*y;
double Lstmp5 = (1.0/2.0)*pow(x, 2);
double Lstmp6 = (1.0/2.0)*pow(y, 2);
double Lstmp7 = (1.0/2.0)*pow(z, 2);
double Lstmp8 = x*L[13];
double Lstmp9 = x*L[15];
double Lstmp10 = y*L[11];
double Lstmp11 = z*L[12];
double Lstmp12 = y*L[18];
double Lstmp13 = z*L[17];
double Lstmp14 = y*L[13];
double Lstmp15 = y*L[14];
double Lstmp16 = z*L[15];
double Lstmp17 = z*L[18];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp10*Lstmp5 + Lstmp11*Lstmp5 + Lstmp12*Lstmp7 + Lstmp13*Lstmp6 + Lstmp2*y + Lstmp4*x + Lstmp5*L[4] + Lstmp6*Lstmp8 + Lstmp6*L[7] + Lstmp7*Lstmp9 + Lstmp7*L[9] + (1.0/6.0)*pow(x, 3)*L[10] + x*L[1] + (1.0/6.0)*pow(y, 3)*L[16] + y*L[2] + (1.0/6.0)*pow(z, 3)*L[19] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + Lstmp10*x + Lstmp11*x + Lstmp4 + Lstmp5*L[10] + Lstmp6*L[13] + Lstmp7*L[15] + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp13*y + Lstmp14*x + Lstmp2 + Lstmp3*x + Lstmp5*L[11] + Lstmp6*L[16] + Lstmp7*L[18] + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += Lstmp15*x + Lstmp16*x + Lstmp17*y + Lstmp5*L[12] + Lstmp6*L[17] + Lstmp7*L[19] + x*L[6] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += Lstmp10 + Lstmp11 + x*L[10] + L[4];
#pragma omp atomic
Ls[5] += Lstmp14 + Lstmp3 + x*L[11] + L[5];
#pragma omp atomic
Ls[6] += Lstmp15 + Lstmp16 + x*L[12] + L[6];
#pragma omp atomic
Ls[7] += Lstmp13 + Lstmp8 + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += Lstmp17 + x*L[14] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += Lstmp12 + Lstmp9 + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += L[10];
#pragma omp atomic
Ls[11] += L[11];
#pragma omp atomic
Ls[12] += L[12];
#pragma omp atomic
Ls[13] += L[13];
#pragma omp atomic
Ls[14] += L[14];
#pragma omp atomic
Ls[15] += L[15];
#pragma omp atomic
Ls[16] += L[16];
#pragma omp atomic
Ls[17] += L[17];
#pragma omp atomic
Ls[18] += L[18];
#pragma omp atomic
Ls[19] += L[19];

}

void L2P_4(double x, double y, double z, double * L, double * F) {
double Ftmp0 = x*y;
double Ftmp1 = x*z;
double Ftmp2 = y*z;
double Ftmp3 = (1.0/2.0)*pow(x, 2);
double Ftmp4 = (1.0/2.0)*pow(y, 2);
double Ftmp5 = (1.0/2.0)*pow(z, 2);
#pragma omp atomic
F[0] += -Ftmp0*L[11] - Ftmp1*L[12] - Ftmp2*L[14] - Ftmp3*L[10] - Ftmp4*L[13] - Ftmp5*L[15] - x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -Ftmp0*L[13] - Ftmp1*L[14] - Ftmp2*L[17] - Ftmp3*L[11] - Ftmp4*L[16] - Ftmp5*L[18] - x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -Ftmp0*L[14] - Ftmp1*L[15] - Ftmp2*L[18] - Ftmp3*L[12] - Ftmp4*L[17] - Ftmp5*L[19] - x*L[6] - y*L[8] - z*L[9] - L[3];

}

void M2P_4(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(R, -3);
double Ftmp1 = pow(R, -2);
double Ftmp2 = 3.0*Ftmp1;
double Ftmp3 = y*M[4];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = pow(R, -4);
double Ftmp6 = Ftmp5*y;
double Ftmp7 = 15.0*z;
double Ftmp8 = Ftmp7*M[13];
double Ftmp9 = Ftmp2*x;
double Ftmp10 = Ftmp9*y;
double Ftmp11 = Ftmp4*M[2];
double Ftmp12 = pow(x, 2);
double Ftmp13 = Ftmp1*Ftmp12;
double Ftmp14 = y*M[7];
double Ftmp15 = Ftmp5*x;
double Ftmp16 = Ftmp15*Ftmp7;
double Ftmp17 = Ftmp12*Ftmp5;
double Ftmp18 = Ftmp7*M[5];
double Ftmp19 = pow(R, -6);
double Ftmp20 = Ftmp12*Ftmp19;
double Ftmp21 = Ftmp20*z;
double Ftmp22 = y*M[13];
double Ftmp23 = 105.0*Ftmp22;
double Ftmp24 = pow(y, 2);
double Ftmp25 = 15.0*Ftmp1;
double Ftmp26 = -Ftmp24*Ftmp25;
double Ftmp27 = Ftmp1*(Ftmp26 + 3.0);
double Ftmp28 = pow(z, 2);
double Ftmp29 = -Ftmp25*Ftmp28;
double Ftmp30 = Ftmp1*(Ftmp29 + 3.0);
double Ftmp31 = -Ftmp12*Ftmp25;
double Ftmp32 = Ftmp1*(Ftmp31 + 9.0);
double Ftmp33 = -105.0*Ftmp13;
double Ftmp34 = Ftmp33 + 45.0;
double Ftmp35 = Ftmp34*Ftmp5;
double Ftmp36 = y*M[20];
double Ftmp37 = z*M[21];
double Ftmp38 = Ftmp1*Ftmp28;
double Ftmp39 = 3.0*M[27];
double Ftmp40 = Ftmp39*(5.0 - 35.0*Ftmp38);
double Ftmp41 = 105.0*Ftmp1;
double Ftmp42 = -Ftmp24*Ftmp41;
double Ftmp43 = Ftmp42 + 45.0;
double Ftmp44 = Ftmp43*Ftmp6;
double Ftmp45 = 1.0*M[25];
double Ftmp46 = Ftmp5*z;
double Ftmp47 = 1.0*Ftmp46;
double Ftmp48 = Ftmp42 + 15.0;
double Ftmp49 = Ftmp48*M[26];
double Ftmp50 = -Ftmp28*Ftmp41;
double Ftmp51 = Ftmp50 + 45.0;
double Ftmp52 = Ftmp47*Ftmp51;
double Ftmp53 = Ftmp27*M[6];
double Ftmp54 = Ftmp30*M[8];
double Ftmp55 = Ftmp15*Ftmp34;
double Ftmp56 = Ftmp55*y;
double Ftmp57 = Ftmp15*Ftmp43;
double Ftmp58 = Ftmp57*y;
double Ftmp59 = Ftmp48*M[16];
double Ftmp60 = Ftmp15*z;
double Ftmp61 = Ftmp55*z;
double Ftmp62 = Ftmp51*M[18];
double Ftmp63 = Ftmp50 + 15.0;
double Ftmp64 = 1.0*M[17];
double Ftmp65 = Ftmp63*Ftmp64;
double Ftmp66 = Ftmp15*y;
double Ftmp67 = 1.0*Ftmp17;
double Ftmp68 = Ftmp48*M[12];
double Ftmp69 = Ftmp63*M[14];
double Ftmp70 = -945.0*Ftmp13;
double Ftmp71 = x*(Ftmp70 + 315.0);
double Ftmp72 = Ftmp19*y*z;
double Ftmp73 = Ftmp1*Ftmp24;
double Ftmp74 = -945.0*Ftmp73;
double Ftmp75 = Ftmp74 + 315.0;
double Ftmp76 = Ftmp75*M[30];
double Ftmp77 = Ftmp72*x;
double Ftmp78 = -945.0*Ftmp38;
double Ftmp79 = Ftmp78 + 315.0;
double Ftmp80 = 1.0*M[32];
double Ftmp81 = Ftmp79*Ftmp80;
double Ftmp82 = Ftmp20*y;
double Ftmp83 = 315.0*Ftmp1;
double Ftmp84 = -Ftmp28*Ftmp83;
double Ftmp85 = Ftmp39*(Ftmp84 + 35.0);
double Ftmp86 = Ftmp45*Ftmp75;
double Ftmp87 = Ftmp20*(Ftmp70 + 525.0);
double Ftmp88 = 1.0*Ftmp21;
double Ftmp89 = (Ftmp74 + 105.0)*M[26];
double Ftmp90 = Ftmp79*M[28];
double Ftmp91 = 945.0*Ftmp5;
double Ftmp92 = Ftmp91*pow(y, 4);
double Ftmp93 = 630.0*Ftmp1;
double Ftmp94 = (-Ftmp24*Ftmp93 + Ftmp92 + 45.0)*M[29];
double Ftmp95 = Ftmp91*pow(z, 4);
double Ftmp96 = (-Ftmp28*Ftmp93 + Ftmp95 + 45.0)*M[33];
double Ftmp97 = Ftmp91*pow(x, 4);
double Ftmp98 = Ftmp24*Ftmp91;
double Ftmp99 = Ftmp28*Ftmp98;
double Ftmp100 = -Ftmp24*Ftmp83;
double Ftmp101 = Ftmp12*Ftmp98;
double Ftmp102 = Ftmp12*Ftmp28*Ftmp91;
double Ftmp103 = 15.0*Ftmp24;
double Ftmp104 = Ftmp19*Ftmp24;
double Ftmp105 = Ftmp104*z;
double Ftmp106 = Ftmp1*(Ftmp31 + 3.0);
double Ftmp107 = Ftmp1*(Ftmp26 + 9.0);
double Ftmp108 = Ftmp33 + 15.0;
double Ftmp109 = Ftmp108*M[23];
double Ftmp110 = Ftmp43*M[30];
double Ftmp111 = Ftmp106*M[3];
double Ftmp112 = Ftmp108*M[11];
double Ftmp113 = Ftmp6*z;
double Ftmp114 = Ftmp44*z;
double Ftmp115 = Ftmp24*Ftmp5;
double Ftmp116 = Ftmp108*M[10];
double Ftmp117 = 1.0*Ftmp15;
double Ftmp118 = Ftmp19*Ftmp71;
double Ftmp119 = 1.0*Ftmp77;
double Ftmp120 = Ftmp104*x;
double Ftmp121 = Ftmp74 + 525.0;
double Ftmp122 = (Ftmp70 + 105.0)*M[23];
double Ftmp123 = (-Ftmp12*Ftmp93 + Ftmp97 + 45.0)*M[19];
double Ftmp124 = -315.0*Ftmp13;
double Ftmp125 = 15.0*Ftmp15;
double Ftmp126 = Ftmp19*Ftmp28;
double Ftmp127 = Ftmp126*x;
double Ftmp128 = Ftmp1*(Ftmp29 + 9.0);
double Ftmp129 = Ftmp28*Ftmp5;
double Ftmp130 = 1.0*Ftmp60;
double Ftmp131 = 1.0*Ftmp127;
double Ftmp132 = Ftmp78 + 525.0;
double Ftmp133 = Ftmp126*y;
#pragma omp atomic
F[0] += Ftmp0*(-Ftmp10*M[1] - Ftmp11*x - 3.0*Ftmp13*M[0] + Ftmp14*Ftmp16 + Ftmp15*Ftmp94 + Ftmp15*Ftmp96 + Ftmp15*(Ftmp100 + Ftmp101 + Ftmp34)*M[22] + Ftmp15*(Ftmp102 + Ftmp34 + Ftmp84)*M[24] + Ftmp15*(-1050.0*Ftmp13 + Ftmp97 + 225.0)*M[19] + Ftmp15*(Ftmp48 + Ftmp50 + Ftmp99)*M[31] + Ftmp17*Ftmp18 + 15.0*Ftmp17*Ftmp3 + Ftmp17*(Ftmp33 + 75.0)*M[9] - Ftmp2*Ftmp3 - Ftmp21*Ftmp23 - Ftmp27*M[12] - Ftmp30*M[14] - Ftmp32*x*M[3] - Ftmp32*M[9] + Ftmp35*Ftmp36 + Ftmp35*Ftmp37 - Ftmp36*Ftmp87 - Ftmp37*Ftmp87 - Ftmp4*M[5] + Ftmp40*Ftmp6 + Ftmp44*Ftmp45 + Ftmp47*Ftmp49 + Ftmp52*M[28] - Ftmp53*x - Ftmp54*x + Ftmp56*M[10] + Ftmp58*M[15] + Ftmp59*Ftmp60 + Ftmp6*Ftmp8 + Ftmp60*Ftmp62 + Ftmp61*M[11] + Ftmp65*Ftmp66 + Ftmp67*Ftmp68 + Ftmp67*Ftmp69 - Ftmp71*Ftmp72*M[23] - Ftmp76*Ftmp77 - Ftmp77*Ftmp81 - Ftmp82*Ftmp85 - Ftmp82*Ftmp86 - Ftmp88*Ftmp89 - Ftmp88*Ftmp90 + 1.0*M[0]);
#pragma omp atomic
F[1] += Ftmp0*(-Ftmp10*M[0] + Ftmp103*Ftmp15*M[4] + Ftmp103*Ftmp46*M[7] - Ftmp104*Ftmp71*M[20] - Ftmp105*Ftmp121*M[30] - Ftmp105*Ftmp122 - Ftmp105*Ftmp81 - 105.0*Ftmp105*x*M[13] - Ftmp106*M[10] - Ftmp107*y*M[6] - Ftmp107*M[15] + Ftmp109*Ftmp46 - Ftmp11*y + Ftmp110*Ftmp46 - Ftmp111*y + Ftmp112*Ftmp113 + Ftmp113*Ftmp62 + Ftmp114*M[16] + Ftmp115*Ftmp116 + Ftmp115*Ftmp65 + Ftmp115*(Ftmp42 + 75.0)*M[15] + Ftmp117*Ftmp69*y - Ftmp118*Ftmp37*y - Ftmp119*Ftmp75*M[26] - Ftmp119*Ftmp90 - Ftmp120*Ftmp121*Ftmp45 - Ftmp120*Ftmp85 + Ftmp123*Ftmp6 + Ftmp15*Ftmp40 + Ftmp15*Ftmp8 + Ftmp18*Ftmp66 - Ftmp30*M[17] - Ftmp4*M[7] + Ftmp45*Ftmp57 + Ftmp52*M[32] - Ftmp54*y + Ftmp55*M[20] + Ftmp56*M[9] + 1.0*Ftmp58*M[12] + Ftmp6*Ftmp96 + Ftmp6*(Ftmp101 + Ftmp124 + Ftmp43)*M[22] + Ftmp6*(Ftmp102 + Ftmp33 + Ftmp63)*M[24] + Ftmp6*(Ftmp43 + Ftmp84 + Ftmp99)*M[31] + Ftmp6*(-1050.0*Ftmp73 + Ftmp92 + 225.0)*M[29] - 3.0*Ftmp73*M[1] - Ftmp9*M[4] + 1.0*M[1]);
#pragma omp atomic
F[2] += Ftmp0*(-Ftmp106*M[11] + Ftmp109*Ftmp6 + Ftmp110*Ftmp6 - Ftmp111*z + Ftmp112*Ftmp129 + Ftmp113*Ftmp116 + Ftmp113*Ftmp51*Ftmp64 + Ftmp114*M[15] + Ftmp117*Ftmp49 + Ftmp117*Ftmp51*M[28] - Ftmp118*Ftmp36*z - Ftmp122*Ftmp133 + Ftmp123*Ftmp46 + Ftmp125*Ftmp22 + Ftmp125*Ftmp28*M[5] - Ftmp126*Ftmp71*M[21] - Ftmp127*Ftmp23 - Ftmp128*z*M[8] - Ftmp128*M[18] + Ftmp129*Ftmp59 + Ftmp129*(Ftmp50 + 75.0)*M[18] + Ftmp130*Ftmp51*M[14] + Ftmp130*Ftmp68 - Ftmp131*Ftmp132*M[28] - Ftmp131*Ftmp89 - Ftmp132*Ftmp133*Ftmp80 - Ftmp133*Ftmp76 - Ftmp14*Ftmp2 + Ftmp16*Ftmp3 - Ftmp27*M[16] + 15.0*Ftmp28*Ftmp6*M[7] - 3.0*Ftmp38*M[2] - Ftmp39*Ftmp77*(Ftmp84 + 105.0) - Ftmp4*x*M[0] - Ftmp4*y*M[1] + Ftmp46*Ftmp94 + Ftmp46*(Ftmp100 + Ftmp51 + Ftmp99)*M[31] + Ftmp46*(Ftmp101 + Ftmp33 + Ftmp48)*M[22] + Ftmp46*(Ftmp102 + Ftmp124 + Ftmp51)*M[24] + Ftmp46*(-1050.0*Ftmp38 + Ftmp95 + 225.0)*M[33] + Ftmp51*Ftmp6*Ftmp80 - Ftmp53*z + Ftmp55*M[21] + Ftmp61*M[9] - Ftmp77*Ftmp86 - Ftmp9*M[5] + 1.0*M[2]);

}

void P2M_5(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = q*z;
double Mtmp3 = pow(x, 2);
double Mtmp4 = (1.0/2.0)*q;
double Mtmp5 = Mtmp0*y;
double Mtmp6 = Mtmp0*z;
double Mtmp7 = pow(y, 2);
double Mtmp8 = Mtmp1*z;
double Mtmp9 = pow(z, 2);
double Mtmp10 = pow(x, 3);
double Mtmp11 = (1.0/6.0)*q;
double Mtmp12 = (1.0/2.0)*Mtmp3;
double Mtmp13 = (1.0/2.0)*Mtmp0;
double Mtmp14 = pow(y, 3);
double Mtmp15 = (1.0/2.0)*Mtmp7;
double Mtmp16 = (1.0/2.0)*Mtmp9;
double Mtmp17 = pow(z, 3);
double Mtmp18 = pow(x, 4);
double Mtmp19 = (1.0/24.0)*q;
double Mtmp20 = (1.0/6.0)*Mtmp10;
double Mtmp21 = Mtmp7*q;
double Mtmp22 = (1.0/4.0)*Mtmp3;
double Mtmp23 = Mtmp9*q;
double Mtmp24 = (1.0/6.0)*Mtmp0;
double Mtmp25 = pow(y, 4);
double Mtmp26 = (1.0/6.0)*Mtmp14;
double Mtmp27 = (1.0/4.0)*Mtmp9;
double Mtmp28 = (1.0/6.0)*Mtmp17;
double Mtmp29 = pow(z, 4);
double Mtmp30 = (1.0/120.0)*q;
double Mtmp31 = (1.0/24.0)*Mtmp18;
double Mtmp32 = (1.0/12.0)*Mtmp10;
double Mtmp33 = (1.0/12.0)*Mtmp14;
double Mtmp34 = Mtmp3*q;
double Mtmp35 = (1.0/12.0)*Mtmp17;
double Mtmp36 = (1.0/24.0)*Mtmp0;
M[0] += -Mtmp0;
M[1] += -Mtmp1;
M[2] += -Mtmp2;
M[3] += Mtmp3*Mtmp4;
M[4] += Mtmp5;
M[5] += Mtmp6;
M[6] += Mtmp4*Mtmp7;
M[7] += Mtmp8;
M[8] += Mtmp4*Mtmp9;
M[9] += -Mtmp10*Mtmp11;
M[10] += -Mtmp1*Mtmp12;
M[11] += -Mtmp12*Mtmp2;
M[12] += -Mtmp13*Mtmp7;
M[13] += -Mtmp5*z;
M[14] += -Mtmp13*Mtmp9;
M[15] += -Mtmp11*Mtmp14;
M[16] += -Mtmp15*Mtmp2;
M[17] += -Mtmp1*Mtmp16;
M[18] += -Mtmp11*Mtmp17;
M[19] += Mtmp18*Mtmp19;
M[20] += Mtmp1*Mtmp20;
M[21] += Mtmp2*Mtmp20;
M[22] += Mtmp21*Mtmp22;
M[23] += Mtmp12*Mtmp8;
M[24] += Mtmp22*Mtmp23;
M[25] += Mtmp14*Mtmp24;
M[26] += Mtmp15*Mtmp6;
M[27] += Mtmp16*Mtmp5;
M[28] += Mtmp17*Mtmp24;
M[29] += Mtmp19*Mtmp25;
M[30] += Mtmp2*Mtmp26;
M[31] += Mtmp21*Mtmp27;
M[32] += Mtmp1*Mtmp28;
M[33] += Mtmp19*Mtmp29;
M[34] += -Mtmp30*pow(x, 5);
M[35] += -Mtmp1*Mtmp31;
M[36] += -Mtmp2*Mtmp31;
M[37] += -Mtmp21*Mtmp32;
M[38] += -Mtmp20*Mtmp8;
M[39] += -Mtmp23*Mtmp32;
M[40] += -Mtmp33*Mtmp34;
M[41] += -Mtmp2*Mtmp22*Mtmp7;
M[42] += -Mtmp1*Mtmp22*Mtmp9;
M[43] += -Mtmp34*Mtmp35;
M[44] += -Mtmp25*Mtmp36;
M[45] += -Mtmp26*Mtmp6;
M[46] += -Mtmp0*Mtmp27*Mtmp7;
M[47] += -Mtmp28*Mtmp5;
M[48] += -Mtmp29*Mtmp36;
M[49] += -Mtmp30*pow(y, 5);
M[50] += -1.0/24.0*Mtmp2*Mtmp25;
M[51] += -Mtmp23*Mtmp33;
M[52] += -Mtmp21*Mtmp35;
M[53] += -1.0/24.0*Mtmp1*Mtmp29;
M[54] += -Mtmp30*pow(z, 5);

}
void M2M_5(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = x*M[1];
double Mstmp2 = y*M[0];
double Mstmp3 = x*M[2];
double Mstmp4 = z*M[0];
double Mstmp5 = y*M[1];
double Mstmp6 = y*M[2];
double Mstmp7 = z*M[1];
double Mstmp8 = z*M[2];
double Mstmp9 = x*M[3];
double Mstmp10 = pow(x, 2);
double Mstmp11 = (1.0/2.0)*Mstmp10;
double Mstmp12 = x*M[4];
double Mstmp13 = y*M[3];
double Mstmp14 = Mstmp0*y;
double Mstmp15 = x*M[5];
double Mstmp16 = z*M[3];
double Mstmp17 = Mstmp0*z;
double Mstmp18 = x*M[6];
double Mstmp19 = y*M[4];
double Mstmp20 = Mstmp1*y;
double Mstmp21 = pow(y, 2);
double Mstmp22 = (1.0/2.0)*M[0];
double Mstmp23 = x*M[7];
double Mstmp24 = y*M[5];
double Mstmp25 = z*M[4];
double Mstmp26 = Mstmp3*y;
double Mstmp27 = Mstmp1*z;
double Mstmp28 = Mstmp2*z;
double Mstmp29 = x*M[8];
double Mstmp30 = z*M[5];
double Mstmp31 = Mstmp3*z;
double Mstmp32 = pow(z, 2);
double Mstmp33 = y*M[6];
double Mstmp34 = (1.0/2.0)*Mstmp21;
double Mstmp35 = y*M[7];
double Mstmp36 = z*M[6];
double Mstmp37 = Mstmp5*z;
double Mstmp38 = y*M[8];
double Mstmp39 = z*M[7];
double Mstmp40 = Mstmp6*z;
double Mstmp41 = (1.0/2.0)*Mstmp32;
double Mstmp42 = z*M[8];
double Mstmp43 = x*M[9];
double Mstmp44 = (1.0/6.0)*pow(x, 3);
double Mstmp45 = x*M[10];
double Mstmp46 = y*M[9];
double Mstmp47 = Mstmp9*y;
double Mstmp48 = x*M[11];
double Mstmp49 = x*M[12];
double Mstmp50 = y*M[10];
double Mstmp51 = Mstmp12*y;
double Mstmp52 = x*M[13];
double Mstmp53 = y*M[11];
double Mstmp54 = Mstmp15*y;
double Mstmp55 = x*M[14];
double Mstmp56 = x*M[15];
double Mstmp57 = y*M[12];
double Mstmp58 = Mstmp18*y;
double Mstmp59 = pow(y, 3);
double Mstmp60 = (1.0/6.0)*M[0];
double Mstmp61 = x*M[16];
double Mstmp62 = y*M[13];
double Mstmp63 = Mstmp23*y;
double Mstmp64 = x*M[17];
double Mstmp65 = y*M[14];
double Mstmp66 = Mstmp29*y;
double Mstmp67 = x*M[18];
double Mstmp68 = pow(z, 3);
double Mstmp69 = y*M[15];
double Mstmp70 = (1.0/6.0)*Mstmp59;
double Mstmp71 = y*M[16];
double Mstmp72 = y*M[17];
double Mstmp73 = y*M[18];
double Mstmp74 = (1.0/6.0)*Mstmp68;
double Mstmp75 = (1.0/24.0)*pow(x, 4);
double Mstmp76 = (1.0/4.0)*Mstmp10;
double Mstmp77 = Mstmp76*M[0];
double Mstmp78 = Mstmp21*Mstmp76;
double Mstmp79 = Mstmp32*Mstmp76;
double Mstmp80 = pow(y, 4);
double Mstmp81 = (1.0/24.0)*M[0];
double Mstmp82 = (1.0/4.0)*Mstmp21*Mstmp32;
double Mstmp83 = pow(z, 4);
double Mstmp84 = (1.0/24.0)*Mstmp80;
double Mstmp85 = (1.0/24.0)*Mstmp83;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += Mstmp0 + M[3];
#pragma omp atomic
Ms[4] += Mstmp1 + Mstmp2 + M[4];
#pragma omp atomic
Ms[5] += Mstmp3 + Mstmp4 + M[5];
#pragma omp atomic
Ms[6] += Mstmp5 + M[6];
#pragma omp atomic
Ms[7] += Mstmp6 + Mstmp7 + M[7];
#pragma omp atomic
Ms[8] += Mstmp8 + M[8];
#pragma omp atomic
Ms[9] += Mstmp11*M[0] + Mstmp9 + M[9];
#pragma omp atomic
Ms[10] += Mstmp11*M[1] + Mstmp12 + Mstmp13 + Mstmp14 + M[10];
#pragma omp atomic
Ms[11] += Mstmp11*M[2] + Mstmp15 + Mstmp16 + Mstmp17 + M[11];
#pragma omp atomic
Ms[12] += Mstmp18 + Mstmp19 + Mstmp20 + Mstmp21*Mstmp22 + M[12];
#pragma omp atomic
Ms[13] += Mstmp23 + Mstmp24 + Mstmp25 + Mstmp26 + Mstmp27 + Mstmp28 + M[13];
#pragma omp atomic
Ms[14] += Mstmp22*Mstmp32 + Mstmp29 + Mstmp30 + Mstmp31 + M[14];
#pragma omp atomic
Ms[15] += Mstmp33 + Mstmp34*M[1] + M[15];
#pragma omp atomic
Ms[16] += Mstmp34*M[2] + Mstmp35 + Mstmp36 + Mstmp37 + M[16];
#pragma omp atomic
Ms[17] += Mstmp38 + Mstmp39 + Mstmp40 + Mstmp41*M[1] + M[17];
#pragma omp atomic
Ms[18] += Mstmp41*M[2] + Mstmp42 + M[18];
#pragma omp atomic
Ms[19] += Mstmp11*M[3] + Mstmp43 + Mstmp44*M[0] + M[19];
#pragma omp atomic
Ms[20] += Mstmp11*Mstmp2 + Mstmp11*M[4] + Mstmp44*M[1] + Mstmp45 + Mstmp46 + Mstmp47 + M[20];
#pragma omp atomic
Ms[21] += Mstmp11*Mstmp4 + Mstmp11*M[5] + Mstmp44*M[2] + Mstmp48 + Mstmp9*z + z*M[9] + M[21];
#pragma omp atomic
Ms[22] += Mstmp0*Mstmp34 + Mstmp11*Mstmp5 + Mstmp11*M[6] + Mstmp34*M[3] + Mstmp49 + Mstmp50 + Mstmp51 + M[22];
#pragma omp atomic
Ms[23] += Mstmp11*Mstmp6 + Mstmp11*Mstmp7 + Mstmp11*M[7] + Mstmp12*z + Mstmp13*z + Mstmp14*z + Mstmp52 + Mstmp53 + Mstmp54 + z*M[10] + M[23];
#pragma omp atomic
Ms[24] += Mstmp0*Mstmp41 + Mstmp11*Mstmp8 + Mstmp11*M[8] + Mstmp15*z + Mstmp41*M[3] + Mstmp55 + z*M[11] + M[24];
#pragma omp atomic
Ms[25] += Mstmp1*Mstmp34 + Mstmp34*M[4] + Mstmp56 + Mstmp57 + Mstmp58 + Mstmp59*Mstmp60 + M[25];
#pragma omp atomic
Ms[26] += Mstmp18*z + Mstmp19*z + Mstmp20*z + Mstmp3*Mstmp34 + Mstmp34*Mstmp4 + Mstmp34*M[5] + Mstmp61 + Mstmp62 + Mstmp63 + z*M[12] + M[26];
#pragma omp atomic
Ms[27] += Mstmp1*Mstmp41 + Mstmp2*Mstmp41 + Mstmp23*z + Mstmp24*z + Mstmp26*z + Mstmp41*M[4] + Mstmp64 + Mstmp65 + Mstmp66 + z*M[13] + M[27];
#pragma omp atomic
Ms[28] += Mstmp29*z + Mstmp3*Mstmp41 + Mstmp41*M[5] + Mstmp60*Mstmp68 + Mstmp67 + z*M[14] + M[28];
#pragma omp atomic
Ms[29] += Mstmp34*M[6] + Mstmp69 + Mstmp70*M[1] + M[29];
#pragma omp atomic
Ms[30] += Mstmp33*z + Mstmp34*Mstmp7 + Mstmp34*M[7] + Mstmp70*M[2] + Mstmp71 + z*M[15] + M[30];
#pragma omp atomic
Ms[31] += Mstmp34*Mstmp8 + Mstmp34*M[8] + Mstmp35*z + Mstmp41*Mstmp5 + Mstmp41*M[6] + Mstmp72 + z*M[16] + M[31];
#pragma omp atomic
Ms[32] += Mstmp38*z + Mstmp41*Mstmp6 + Mstmp41*M[7] + Mstmp73 + Mstmp74*M[1] + z*M[17] + M[32];
#pragma omp atomic
Ms[33] += Mstmp41*M[8] + Mstmp74*M[2] + z*M[18] + M[33];
#pragma omp atomic
Ms[34] += Mstmp11*M[9] + Mstmp44*M[3] + Mstmp75*M[0] + x*M[19] + M[34];
#pragma omp atomic
Ms[35] += Mstmp11*Mstmp13 + Mstmp11*M[10] + Mstmp2*Mstmp44 + Mstmp43*y + Mstmp44*M[4] + Mstmp75*M[1] + x*M[20] + y*M[19] + M[35];
#pragma omp atomic
Ms[36] += Mstmp11*Mstmp16 + Mstmp11*M[11] + Mstmp4*Mstmp44 + Mstmp43*z + Mstmp44*M[5] + Mstmp75*M[2] + x*M[21] + z*M[19] + M[36];
#pragma omp atomic
Ms[37] += Mstmp11*Mstmp19 + Mstmp11*M[12] + Mstmp21*Mstmp77 + Mstmp34*Mstmp9 + Mstmp34*M[9] + Mstmp44*Mstmp5 + Mstmp44*M[6] + Mstmp45*y + x*M[22] + y*M[20] + M[37];
#pragma omp atomic
Ms[38] += Mstmp11*Mstmp24 + Mstmp11*Mstmp25 + Mstmp11*Mstmp28 + Mstmp11*M[13] + Mstmp44*Mstmp6 + Mstmp44*Mstmp7 + Mstmp44*M[7] + Mstmp45*z + Mstmp46*z + Mstmp47*z + Mstmp48*y + x*M[23] + y*M[21] + z*M[20] + M[38];
#pragma omp atomic
Ms[39] += Mstmp11*Mstmp30 + Mstmp11*M[14] + Mstmp32*Mstmp77 + Mstmp41*Mstmp9 + Mstmp41*M[9] + Mstmp44*Mstmp8 + Mstmp44*M[8] + Mstmp48*z + x*M[24] + z*M[21] + M[39];
#pragma omp atomic
Ms[40] += Mstmp0*Mstmp70 + Mstmp11*Mstmp33 + Mstmp11*M[15] + Mstmp12*Mstmp34 + Mstmp34*M[10] + Mstmp49*y + Mstmp70*M[3] + Mstmp78*M[1] + x*M[25] + y*M[22] + M[40];
#pragma omp atomic
Ms[41] += Mstmp11*Mstmp35 + Mstmp11*Mstmp36 + Mstmp11*Mstmp37 + Mstmp11*M[16] + Mstmp15*Mstmp34 + Mstmp16*Mstmp34 + Mstmp17*Mstmp34 + Mstmp34*M[11] + Mstmp49*z + Mstmp50*z + Mstmp51*z + Mstmp52*y + Mstmp78*M[2] + x*M[26] + y*M[23] + z*M[22] + M[41];
#pragma omp atomic
Ms[42] += Mstmp11*Mstmp38 + Mstmp11*Mstmp39 + Mstmp11*Mstmp40 + Mstmp11*M[17] + Mstmp12*Mstmp41 + Mstmp13*Mstmp41 + Mstmp14*Mstmp41 + Mstmp41*M[10] + Mstmp52*z + Mstmp53*z + Mstmp54*z + Mstmp55*y + Mstmp79*M[1] + x*M[27] + y*M[24] + z*M[23] + M[42];
#pragma omp atomic
Ms[43] += Mstmp0*Mstmp74 + Mstmp11*Mstmp42 + Mstmp11*M[18] + Mstmp15*Mstmp41 + Mstmp41*M[11] + Mstmp55*z + Mstmp74*M[3] + Mstmp79*M[2] + x*M[28] + z*M[24] + M[43];
#pragma omp atomic
Ms[44] += Mstmp1*Mstmp70 + Mstmp18*Mstmp34 + Mstmp34*M[12] + Mstmp56*y + Mstmp70*M[4] + Mstmp80*Mstmp81 + x*M[29] + y*M[25] + M[44];
#pragma omp atomic
Ms[45] += Mstmp23*Mstmp34 + Mstmp25*Mstmp34 + Mstmp27*Mstmp34 + Mstmp3*Mstmp70 + Mstmp34*M[13] + Mstmp4*Mstmp70 + Mstmp56*z + Mstmp57*z + Mstmp58*z + Mstmp61*y + Mstmp70*M[5] + x*M[30] + y*M[26] + z*M[25] + M[45];
#pragma omp atomic
Ms[46] += Mstmp18*Mstmp41 + Mstmp19*Mstmp41 + Mstmp20*Mstmp41 + Mstmp29*Mstmp34 + Mstmp30*Mstmp34 + Mstmp31*Mstmp34 + Mstmp34*M[14] + Mstmp41*M[12] + Mstmp61*z + Mstmp62*z + Mstmp63*z + Mstmp64*y + Mstmp82*M[0] + x*M[31] + y*M[27] + z*M[26] + M[46];
#pragma omp atomic
Ms[47] += Mstmp1*Mstmp74 + Mstmp2*Mstmp74 + Mstmp23*Mstmp41 + Mstmp24*Mstmp41 + Mstmp26*Mstmp41 + Mstmp41*M[13] + Mstmp64*z + Mstmp65*z + Mstmp66*z + Mstmp67*y + Mstmp74*M[4] + x*M[32] + y*M[28] + z*M[27] + M[47];
#pragma omp atomic
Ms[48] += Mstmp29*Mstmp41 + Mstmp3*Mstmp74 + Mstmp41*M[14] + Mstmp67*z + Mstmp74*M[5] + Mstmp81*Mstmp83 + x*M[33] + z*M[28] + M[48];
#pragma omp atomic
Ms[49] += Mstmp34*M[15] + Mstmp70*M[6] + Mstmp84*M[1] + y*M[29] + M[49];
#pragma omp atomic
Ms[50] += Mstmp34*Mstmp36 + Mstmp34*M[16] + Mstmp69*z + Mstmp7*Mstmp70 + Mstmp70*M[7] + Mstmp84*M[2] + y*M[30] + z*M[29] + M[50];
#pragma omp atomic
Ms[51] += Mstmp33*Mstmp41 + Mstmp34*Mstmp39 + Mstmp34*M[17] + Mstmp41*M[15] + Mstmp70*Mstmp8 + Mstmp70*M[8] + Mstmp71*z + Mstmp82*M[1] + y*M[31] + z*M[30] + M[51];
#pragma omp atomic
Ms[52] += Mstmp34*Mstmp42 + Mstmp34*M[18] + Mstmp35*Mstmp41 + Mstmp41*M[16] + Mstmp5*Mstmp74 + Mstmp72*z + Mstmp74*M[6] + Mstmp82*M[2] + y*M[32] + z*M[31] + M[52];
#pragma omp atomic
Ms[53] += Mstmp38*Mstmp41 + Mstmp41*M[17] + Mstmp6*Mstmp74 + Mstmp73*z + Mstmp74*M[7] + Mstmp85*M[1] + y*M[33] + z*M[32] + M[53];
#pragma omp atomic
Ms[54] += Mstmp41*M[18] + Mstmp74*M[8] + Mstmp85*M[2] + z*M[33] + M[54];

}

void M2L_5(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[55];
double Dtmp0 = pow(R, -3);
double Dtmp1 = 1.0*Dtmp0;
double Dtmp2 = pow(x, 2);
double Dtmp3 = pow(R, -2);
double Dtmp4 = 3.0*Dtmp3;
double Dtmp5 = pow(R, -5);
double Dtmp6 = Dtmp5*x;
double Dtmp7 = 3.0*Dtmp6;
double Dtmp8 = pow(y, 2);
double Dtmp9 = Dtmp5*y;
double Dtmp10 = 15.0*Dtmp3;
double Dtmp11 = -Dtmp10*Dtmp2;
double Dtmp12 = Dtmp11 + 3.0;
double Dtmp13 = Dtmp12*Dtmp5;
double Dtmp14 = -Dtmp10*Dtmp8;
double Dtmp15 = Dtmp14 + 3.0;
double Dtmp16 = pow(R, -7);
double Dtmp17 = Dtmp16*x;
double Dtmp18 = y*z;
double Dtmp19 = pow(x, 4);
double Dtmp20 = pow(R, -4);
double Dtmp21 = 105.0*Dtmp20;
double Dtmp22 = Dtmp2*Dtmp3;
double Dtmp23 = -105.0*Dtmp22;
double Dtmp24 = Dtmp23 + 45.0;
double Dtmp25 = Dtmp17*Dtmp24;
double Dtmp26 = Dtmp2*Dtmp8;
double Dtmp27 = Dtmp23 + 15.0;
double Dtmp28 = Dtmp16*y;
double Dtmp29 = Dtmp28*z;
double Dtmp30 = Dtmp3*Dtmp8;
double Dtmp31 = -105.0*Dtmp30;
double Dtmp32 = Dtmp31 + 45.0;
double Dtmp33 = 1.0*Dtmp17;
double Dtmp34 = pow(y, 4);
double Dtmp35 = 945.0*Dtmp20;
double Dtmp36 = Dtmp19*Dtmp35;
double Dtmp37 = Dtmp16*(-630.0*Dtmp22 + Dtmp36 + 45.0);
double Dtmp38 = Dtmp26*Dtmp35;
double Dtmp39 = Dtmp18*x/pow(R, 9);
double Dtmp40 = Dtmp16*z;
double Dtmp41 = Dtmp34*Dtmp35;
double Dtmp42 = -630.0*Dtmp30 + Dtmp41 + 45.0;
D[0] = -Dtmp1*x;
D[1] = -Dtmp1*y;
D[2] = -Dtmp1*z;
D[3] = Dtmp0*(Dtmp2*Dtmp4 - 1.0);
D[4] = Dtmp7*y;
D[5] = Dtmp7*z;
D[6] = Dtmp0*(Dtmp4*Dtmp8 - 1.0);
D[7] = 3.0*Dtmp9*z;
D[8] = -D[3] - D[6];
D[9] = Dtmp6*(Dtmp11 + 9.0);
D[10] = Dtmp13*y;
D[11] = Dtmp13*z;
D[12] = 1.0*Dtmp15*Dtmp6;
D[13] = -15.0*Dtmp17*Dtmp18;
D[14] = -D[9] - D[12];
D[15] = Dtmp9*(Dtmp14 + 9.0);
D[16] = Dtmp15*Dtmp5*z;
D[17] = -D[10] - D[15];
D[18] = -D[11] - D[16];
D[19] = Dtmp5*(Dtmp19*Dtmp21 - 90.0*Dtmp22 + 9.0);
D[20] = -Dtmp25*y;
D[21] = -Dtmp25*z;
D[22] = Dtmp5*(Dtmp12 + Dtmp14 + Dtmp21*Dtmp26);
D[23] = -Dtmp27*Dtmp29;
D[24] = -D[19] - D[22];
D[25] = -Dtmp32*Dtmp33*y;
D[26] = -Dtmp33*z*(Dtmp31 + 15.0);
D[27] = -D[20] - D[25];
D[28] = -D[21] - D[26];
D[29] = Dtmp5*(Dtmp21*Dtmp34 - 90.0*Dtmp30 + 9.0);
D[30] = -Dtmp29*Dtmp32;
D[31] = -D[22] - D[29];
D[32] = -D[23] - D[30];
D[33] = -D[24] - D[31];
D[34] = -Dtmp17*(-1050.0*Dtmp22 + Dtmp36 + 225.0);
D[35] = -Dtmp37*y;
D[36] = -Dtmp37*z;
D[37] = -Dtmp17*(Dtmp24 - 315.0*Dtmp30 + Dtmp38);
D[38] = Dtmp39*(315.0 - 945.0*Dtmp22);
D[39] = -D[34] - D[37];
D[40] = -Dtmp28*(-315.0*Dtmp22 + Dtmp32 + Dtmp38);
D[41] = -Dtmp40*(Dtmp27 + Dtmp31 + Dtmp38);
D[42] = -D[35] - D[40];
D[43] = -D[36] - D[41];
D[44] = -Dtmp33*Dtmp42;
D[45] = 1.0*Dtmp39*(315.0 - 945.0*Dtmp30);
D[46] = -D[37] - D[44];
D[47] = -D[38] - D[45];
D[48] = -D[39] - D[46];
D[49] = -Dtmp28*(-1050.0*Dtmp30 + Dtmp41 + 225.0);
D[50] = -Dtmp40*Dtmp42;
D[51] = -D[40] - D[49];
D[52] = -D[41] - D[50];
D[53] = -D[42] - D[51];
D[54] = -D[43] - D[52];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8] + D[9]*M[9] + D[10]*M[10] + D[11]*M[11] + D[12]*M[12] + D[13]*M[13] + D[14]*M[14] + D[15]*M[15] + D[16]*M[16] + D[17]*M[17] + D[18]*M[18] + D[19]*M[19] + D[20]*M[20] + D[21]*M[21] + D[22]*M[22] + D[23]*M[23] + D[24]*M[24] + D[25]*M[25] + D[26]*M[26] + D[27]*M[27] + D[28]*M[28] + D[29]*M[29] + D[30]*M[30] + D[31]*M[31] + D[32]*M[32] + D[33]*M[33] + D[34]*M[34] + D[35]*M[35] + D[36]*M[36] + D[37]*M[37] + D[38]*M[38] + D[39]*M[39] + D[40]*M[40] + D[41]*M[41] + D[42]*M[42] + D[43]*M[43] + D[44]*M[44] + D[45]*M[45] + D[46]*M[46] + D[47]*M[47] + D[48]*M[48] + D[49]*M[49] + D[50]*M[50] + D[51]*M[51] + D[52]*M[52] + D[53]*M[53] + D[54]*M[54];
#pragma omp atomic
L[1] += D[3]*M[0] + D[4]*M[1] + D[5]*M[2] + D[9]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[14]*M[8] + D[19]*M[9] + D[20]*M[10] + D[21]*M[11] + D[22]*M[12] + D[23]*M[13] + D[24]*M[14] + D[25]*M[15] + D[26]*M[16] + D[27]*M[17] + D[28]*M[18] + D[34]*M[19] + D[35]*M[20] + D[36]*M[21] + D[37]*M[22] + D[38]*M[23] + D[39]*M[24] + D[40]*M[25] + D[41]*M[26] + D[42]*M[27] + D[43]*M[28] + D[44]*M[29] + D[45]*M[30] + D[46]*M[31] + D[47]*M[32] + D[48]*M[33];
#pragma omp atomic
L[2] += D[4]*M[0] + D[6]*M[1] + D[7]*M[2] + D[10]*M[3] + D[12]*M[4] + D[13]*M[5] + D[15]*M[6] + D[16]*M[7] + D[17]*M[8] + D[20]*M[9] + D[22]*M[10] + D[23]*M[11] + D[25]*M[12] + D[26]*M[13] + D[27]*M[14] + D[29]*M[15] + D[30]*M[16] + D[31]*M[17] + D[32]*M[18] + D[35]*M[19] + D[37]*M[20] + D[38]*M[21] + D[40]*M[22] + D[41]*M[23] + D[42]*M[24] + D[44]*M[25] + D[45]*M[26] + D[46]*M[27] + D[47]*M[28] + D[49]*M[29] + D[50]*M[30] + D[51]*M[31] + D[52]*M[32] + D[53]*M[33];
#pragma omp atomic
L[3] += D[5]*M[0] + D[7]*M[1] + D[8]*M[2] + D[11]*M[3] + D[13]*M[4] + D[14]*M[5] + D[16]*M[6] + D[17]*M[7] + D[18]*M[8] + D[21]*M[9] + D[23]*M[10] + D[24]*M[11] + D[26]*M[12] + D[27]*M[13] + D[28]*M[14] + D[30]*M[15] + D[31]*M[16] + D[32]*M[17] + D[33]*M[18] + D[36]*M[19] + D[38]*M[20] + D[39]*M[21] + D[41]*M[22] + D[42]*M[23] + D[43]*M[24] + D[45]*M[25] + D[46]*M[26] + D[47]*M[27] + D[48]*M[28] + D[50]*M[29] + D[51]*M[30] + D[52]*M[31] + D[53]*M[32] + D[54]*M[33];
#pragma omp atomic
L[4] += D[9]*M[0] + D[10]*M[1] + D[11]*M[2] + D[19]*M[3] + D[20]*M[4] + D[21]*M[5] + D[22]*M[6] + D[23]*M[7] + D[24]*M[8] + D[34]*M[9] + D[35]*M[10] + D[36]*M[11] + D[37]*M[12] + D[38]*M[13] + D[39]*M[14] + D[40]*M[15] + D[41]*M[16] + D[42]*M[17] + D[43]*M[18];
#pragma omp atomic
L[5] += D[10]*M[0] + D[12]*M[1] + D[13]*M[2] + D[20]*M[3] + D[22]*M[4] + D[23]*M[5] + D[25]*M[6] + D[26]*M[7] + D[27]*M[8] + D[35]*M[9] + D[37]*M[10] + D[38]*M[11] + D[40]*M[12] + D[41]*M[13] + D[42]*M[14] + D[44]*M[15] + D[45]*M[16] + D[46]*M[17] + D[47]*M[18];
#pragma omp atomic
L[6] += D[11]*M[0] + D[13]*M[1] + D[14]*M[2] + D[21]*M[3] + D[23]*M[4] + D[24]*M[5] + D[26]*M[6] + D[27]*M[7] + D[28]*M[8] + D[36]*M[9] + D[38]*M[10] + D[39]*M[11] + D[41]*M[12] + D[42]*M[13] + D[43]*M[14] + D[45]*M[15] + D[46]*M[16] + D[47]*M[17] + D[48]*M[18];
#pragma omp atomic
L[7] += D[12]*M[0] + D[15]*M[1] + D[16]*M[2] + D[22]*M[3] + D[25]*M[4] + D[26]*M[5] + D[29]*M[6] + D[30]*M[7] + D[31]*M[8] + D[37]*M[9] + D[40]*M[10] + D[41]*M[11] + D[44]*M[12] + D[45]*M[13] + D[46]*M[14] + D[49]*M[15] + D[50]*M[16] + D[51]*M[17] + D[52]*M[18];
#pragma omp atomic
L[8] += D[13]*M[0] + D[16]*M[1] + D[17]*M[2] + D[23]*M[3] + D[26]*M[4] + D[27]*M[5] + D[30]*M[6] + D[31]*M[7] + D[32]*M[8] + D[38]*M[9] + D[41]*M[10] + D[42]*M[11] + D[45]*M[12] + D[46]*M[13] + D[47]*M[14] + D[50]*M[15] + D[51]*M[16] + D[52]*M[17] + D[53]*M[18];
#pragma omp atomic
L[9] += D[14]*M[0] + D[17]*M[1] + D[18]*M[2] + D[24]*M[3] + D[27]*M[4] + D[28]*M[5] + D[31]*M[6] + D[32]*M[7] + D[33]*M[8] + D[39]*M[9] + D[42]*M[10] + D[43]*M[11] + D[46]*M[12] + D[47]*M[13] + D[48]*M[14] + D[51]*M[15] + D[52]*M[16] + D[53]*M[17] + D[54]*M[18];
#pragma omp atomic
L[10] += D[19]*M[0] + D[20]*M[1] + D[21]*M[2] + D[34]*M[3] + D[35]*M[4] + D[36]*M[5] + D[37]*M[6] + D[38]*M[7] + D[39]*M[8];
#pragma omp atomic
L[11] += D[20]*M[0] + D[22]*M[1] + D[23]*M[2] + D[35]*M[3] + D[37]*M[4] + D[38]*M[5] + D[40]*M[6] + D[41]*M[7] + D[42]*M[8];
#pragma omp atomic
L[12] += D[21]*M[0] + D[23]*M[1] + D[24]*M[2] + D[36]*M[3] + D[38]*M[4] + D[39]*M[5] + D[41]*M[6] + D[42]*M[7] + D[43]*M[8];
#pragma omp atomic
L[13] += D[22]*M[0] + D[25]*M[1] + D[26]*M[2] + D[37]*M[3] + D[40]*M[4] + D[41]*M[5] + D[44]*M[6] + D[45]*M[7] + D[46]*M[8];
#pragma omp atomic
L[14] += D[23]*M[0] + D[26]*M[1] + D[27]*M[2] + D[38]*M[3] + D[41]*M[4] + D[42]*M[5] + D[45]*M[6] + D[46]*M[7] + D[47]*M[8];
#pragma omp atomic
L[15] += D[24]*M[0] + D[27]*M[1] + D[28]*M[2] + D[39]*M[3] + D[42]*M[4] + D[43]*M[5] + D[46]*M[6] + D[47]*M[7] + D[48]*M[8];
#pragma omp atomic
L[16] += D[25]*M[0] + D[29]*M[1] + D[30]*M[2] + D[40]*M[3] + D[44]*M[4] + D[45]*M[5] + D[49]*M[6] + D[50]*M[7] + D[51]*M[8];
#pragma omp atomic
L[17] += D[26]*M[0] + D[30]*M[1] + D[31]*M[2] + D[41]*M[3] + D[45]*M[4] + D[46]*M[5] + D[50]*M[6] + D[51]*M[7] + D[52]*M[8];
#pragma omp atomic
L[18] += D[27]*M[0] + D[31]*M[1] + D[32]*M[2] + D[42]*M[3] + D[46]*M[4] + D[47]*M[5] + D[51]*M[6] + D[52]*M[7] + D[53]*M[8];
#pragma omp atomic
L[19] += D[28]*M[0] + D[32]*M[1] + D[33]*M[2] + D[43]*M[3] + D[47]*M[4] + D[48]*M[5] + D[52]*M[6] + D[53]*M[7] + D[54]*M[8];
#pragma omp atomic
L[20] += D[34]*M[0] + D[35]*M[1] + D[36]*M[2];
#pragma omp atomic
L[21] += D[35]*M[0] + D[37]*M[1] + D[38]*M[2];
#pragma omp atomic
L[22] += D[36]*M[0] + D[38]*M[1] + D[39]*M[2];
#pragma omp atomic
L[23] += D[37]*M[0] + D[40]*M[1] + D[41]*M[2];
#pragma omp atomic
L[24] += D[38]*M[0] + D[41]*M[1] + D[42]*M[2];
#pragma omp atomic
L[25] += D[39]*M[0] + D[42]*M[1] + D[43]*M[2];
#pragma omp atomic
L[26] += D[40]*M[0] + D[44]*M[1] + D[45]*M[2];
#pragma omp atomic
L[27] += D[41]*M[0] + D[45]*M[1] + D[46]*M[2];
#pragma omp atomic
L[28] += D[42]*M[0] + D[46]*M[1] + D[47]*M[2];
#pragma omp atomic
L[29] += D[43]*M[0] + D[47]*M[1] + D[48]*M[2];
#pragma omp atomic
L[30] += D[44]*M[0] + D[49]*M[1] + D[50]*M[2];
#pragma omp atomic
L[31] += D[45]*M[0] + D[50]*M[1] + D[51]*M[2];
#pragma omp atomic
L[32] += D[46]*M[0] + D[51]*M[1] + D[52]*M[2];
#pragma omp atomic
L[33] += D[47]*M[0] + D[52]*M[1] + D[53]*M[2];
#pragma omp atomic
L[34] += D[48]*M[0] + D[53]*M[1] + D[54]*M[2];

}

void L2L_5(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[14];
double Lstmp4 = Lstmp3*y;
double Lstmp5 = pow(x, 2);
double Lstmp6 = (1.0/2.0)*Lstmp5;
double Lstmp7 = (1.0/6.0)*pow(x, 3);
double Lstmp8 = pow(y, 2);
double Lstmp9 = (1.0/2.0)*Lstmp8;
double Lstmp10 = (1.0/6.0)*pow(y, 3);
double Lstmp11 = pow(z, 2);
double Lstmp12 = (1.0/2.0)*Lstmp11;
double Lstmp13 = (1.0/6.0)*pow(z, 3);
double Lstmp14 = x*L[13];
double Lstmp15 = x*L[26];
double Lstmp16 = x*L[15];
double Lstmp17 = x*L[29];
double Lstmp18 = y*L[11];
double Lstmp19 = z*L[12];
double Lstmp20 = y*L[21];
double Lstmp21 = z*L[22];
double Lstmp22 = y*L[18];
double Lstmp23 = y*L[33];
double Lstmp24 = z*L[17];
double Lstmp25 = z*L[31];
double Lstmp26 = y*L[28];
double Lstmp27 = Lstmp26*x;
double Lstmp28 = z*L[27];
double Lstmp29 = Lstmp28*x;
double Lstmp30 = z*L[24];
double Lstmp31 = Lstmp30*y;
double Lstmp32 = (1.0/4.0)*Lstmp5;
double Lstmp33 = x*L[23];
double Lstmp34 = x*L[25];
double Lstmp35 = y*L[13];
double Lstmp36 = Lstmp28*y;
double Lstmp37 = x*L[28];
double Lstmp38 = y*L[23];
double Lstmp39 = y*L[32];
double Lstmp40 = y*L[14];
double Lstmp41 = z*L[15];
double Lstmp42 = z*L[18];
double Lstmp43 = z*L[28];
double Lstmp44 = Lstmp43*y;
double Lstmp45 = x*L[27];
double Lstmp46 = y*L[24];
double Lstmp47 = z*L[25];
double Lstmp48 = z*L[32];
double Lstmp49 = y*L[26];
double Lstmp50 = y*L[27];
double Lstmp51 = z*L[29];
double Lstmp52 = z*L[33];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp10*Lstmp15 + Lstmp10*Lstmp25 + Lstmp10*L[16] + Lstmp11*Lstmp32*L[25] + (1.0/4.0)*Lstmp11*Lstmp8*L[32] + Lstmp12*Lstmp16 + Lstmp12*Lstmp22 + Lstmp12*Lstmp27 + Lstmp12*L[9] + Lstmp13*Lstmp17 + Lstmp13*Lstmp23 + Lstmp13*L[19] + Lstmp14*Lstmp9 + Lstmp18*Lstmp6 + Lstmp19*Lstmp6 + Lstmp2*y + Lstmp20*Lstmp7 + Lstmp21*Lstmp7 + Lstmp24*Lstmp9 + Lstmp29*Lstmp9 + Lstmp31*Lstmp6 + Lstmp32*Lstmp8*L[23] + Lstmp4*x + Lstmp6*L[4] + Lstmp7*L[10] + Lstmp9*L[7] + (1.0/24.0)*pow(x, 4)*L[20] + x*L[1] + (1.0/24.0)*pow(y, 4)*L[30] + y*L[2] + (1.0/24.0)*pow(z, 4)*L[34] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + Lstmp10*L[26] + Lstmp12*Lstmp26 + Lstmp12*Lstmp34 + Lstmp12*L[15] + Lstmp13*L[29] + Lstmp18*x + Lstmp19*x + Lstmp20*Lstmp6 + Lstmp21*Lstmp6 + Lstmp28*Lstmp9 + Lstmp31*x + Lstmp33*Lstmp9 + Lstmp4 + Lstmp6*L[10] + Lstmp7*L[20] + Lstmp9*L[13] + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp10*L[30] + Lstmp12*Lstmp37 + Lstmp12*Lstmp39 + Lstmp12*L[18] + Lstmp13*L[33] + Lstmp15*Lstmp9 + Lstmp2 + Lstmp24*y + Lstmp25*Lstmp9 + Lstmp3*x + Lstmp30*Lstmp6 + Lstmp35*x + Lstmp36*x + Lstmp38*Lstmp6 + Lstmp6*L[11] + Lstmp7*L[21] + Lstmp9*L[16] + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += Lstmp10*L[31] + Lstmp12*Lstmp17 + Lstmp12*Lstmp23 + Lstmp12*L[19] + Lstmp13*L[34] + Lstmp40*x + Lstmp41*x + Lstmp42*y + Lstmp44*x + Lstmp45*Lstmp9 + Lstmp46*Lstmp6 + Lstmp47*Lstmp6 + Lstmp48*Lstmp9 + Lstmp6*L[12] + Lstmp7*L[22] + Lstmp9*L[17] + x*L[6] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += Lstmp12*L[25] + Lstmp18 + Lstmp19 + Lstmp20*x + Lstmp21*x + Lstmp31 + Lstmp6*L[20] + Lstmp9*L[23] + x*L[10] + L[4];
#pragma omp atomic
Ls[5] += Lstmp12*L[28] + Lstmp3 + Lstmp30*x + Lstmp35 + Lstmp36 + Lstmp38*x + Lstmp6*L[21] + Lstmp9*L[26] + x*L[11] + L[5];
#pragma omp atomic
Ls[6] += Lstmp12*L[29] + Lstmp40 + Lstmp41 + Lstmp44 + Lstmp46*x + Lstmp47*x + Lstmp6*L[22] + Lstmp9*L[27] + x*L[12] + L[6];
#pragma omp atomic
Ls[7] += Lstmp12*L[32] + Lstmp14 + Lstmp24 + Lstmp25*y + Lstmp29 + Lstmp49*x + Lstmp6*L[23] + Lstmp9*L[30] + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += Lstmp12*L[33] + Lstmp42 + Lstmp43*x + Lstmp48*y + Lstmp50*x + Lstmp6*L[24] + Lstmp9*L[31] + x*L[14] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += Lstmp12*L[34] + Lstmp16 + Lstmp22 + Lstmp27 + Lstmp51*x + Lstmp52*y + Lstmp6*L[25] + Lstmp9*L[32] + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += Lstmp20 + Lstmp21 + x*L[20] + L[10];
#pragma omp atomic
Ls[11] += Lstmp30 + Lstmp38 + x*L[21] + L[11];
#pragma omp atomic
Ls[12] += Lstmp46 + Lstmp47 + x*L[22] + L[12];
#pragma omp atomic
Ls[13] += Lstmp28 + Lstmp33 + Lstmp49 + L[13];
#pragma omp atomic
Ls[14] += Lstmp43 + Lstmp50 + x*L[24] + L[14];
#pragma omp atomic
Ls[15] += Lstmp26 + Lstmp34 + Lstmp51 + L[15];
#pragma omp atomic
Ls[16] += Lstmp15 + Lstmp25 + y*L[30] + L[16];
#pragma omp atomic
Ls[17] += Lstmp45 + Lstmp48 + y*L[31] + L[17];
#pragma omp atomic
Ls[18] += Lstmp37 + Lstmp39 + Lstmp52 + L[18];
#pragma omp atomic
Ls[19] += Lstmp17 + Lstmp23 + z*L[34] + L[19];
#pragma omp atomic
Ls[20] += L[20];
#pragma omp atomic
Ls[21] += L[21];
#pragma omp atomic
Ls[22] += L[22];
#pragma omp atomic
Ls[23] += L[23];
#pragma omp atomic
Ls[24] += L[24];
#pragma omp atomic
Ls[25] += L[25];
#pragma omp atomic
Ls[26] += L[26];
#pragma omp atomic
Ls[27] += L[27];
#pragma omp atomic
Ls[28] += L[28];
#pragma omp atomic
Ls[29] += L[29];
#pragma omp atomic
Ls[30] += L[30];
#pragma omp atomic
Ls[31] += L[31];
#pragma omp atomic
Ls[32] += L[32];
#pragma omp atomic
Ls[33] += L[33];
#pragma omp atomic
Ls[34] += L[34];

}

void L2P_5(double x, double y, double z, double * L, double * F) {
double Ftmp0 = x*y;
double Ftmp1 = x*z;
double Ftmp2 = y*z;
double Ftmp3 = Ftmp0*z;
double Ftmp4 = (1.0/2.0)*pow(x, 2);
double Ftmp5 = (1.0/6.0)*pow(x, 3);
double Ftmp6 = (1.0/2.0)*pow(y, 2);
double Ftmp7 = (1.0/6.0)*pow(y, 3);
double Ftmp8 = (1.0/2.0)*pow(z, 2);
double Ftmp9 = (1.0/6.0)*pow(z, 3);
double Ftmp10 = Ftmp6*x;
double Ftmp11 = Ftmp8*x;
double Ftmp12 = Ftmp4*y;
double Ftmp13 = Ftmp4*z;
double Ftmp14 = Ftmp8*y;
double Ftmp15 = Ftmp6*z;
#pragma omp atomic
F[0] += -Ftmp0*L[11] - Ftmp1*L[12] - Ftmp10*L[23] - Ftmp11*L[25] - Ftmp12*L[21] - Ftmp13*L[22] - Ftmp14*L[28] - Ftmp15*L[27] - Ftmp2*L[14] - Ftmp3*L[24] - Ftmp4*L[10] - Ftmp5*L[20] - Ftmp6*L[13] - Ftmp7*L[26] - Ftmp8*L[15] - Ftmp9*L[29] - x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -Ftmp0*L[13] - Ftmp1*L[14] - Ftmp10*L[26] - Ftmp11*L[28] - Ftmp12*L[23] - Ftmp13*L[24] - Ftmp14*L[32] - Ftmp15*L[31] - Ftmp2*L[17] - Ftmp3*L[27] - Ftmp4*L[11] - Ftmp5*L[21] - Ftmp6*L[16] - Ftmp7*L[30] - Ftmp8*L[18] - Ftmp9*L[33] - x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -Ftmp0*L[14] - Ftmp1*L[15] - Ftmp10*L[27] - Ftmp11*L[29] - Ftmp12*L[24] - Ftmp13*L[25] - Ftmp14*L[33] - Ftmp15*L[32] - Ftmp2*L[18] - Ftmp3*L[28] - Ftmp4*L[12] - Ftmp5*L[22] - Ftmp6*L[17] - Ftmp7*L[31] - Ftmp8*L[19] - Ftmp9*L[34] - x*L[6] - y*L[8] - z*L[9] - L[3];

}

void M2P_5(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(R, -3);
double Ftmp1 = pow(R, -2);
double Ftmp2 = 3.0*Ftmp1;
double Ftmp3 = y*M[4];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = pow(R, -4);
double Ftmp6 = Ftmp5*z;
double Ftmp7 = 15.0*Ftmp6;
double Ftmp8 = y*M[13];
double Ftmp9 = Ftmp2*x;
double Ftmp10 = Ftmp9*y;
double Ftmp11 = Ftmp4*M[2];
double Ftmp12 = pow(x, 2);
double Ftmp13 = Ftmp1*Ftmp12;
double Ftmp14 = y*M[7];
double Ftmp15 = Ftmp7*x;
double Ftmp16 = Ftmp12*Ftmp5;
double Ftmp17 = 15.0*Ftmp16;
double Ftmp18 = pow(R, -6);
double Ftmp19 = Ftmp12*Ftmp18;
double Ftmp20 = 105.0*Ftmp8;
double Ftmp21 = pow(y, 2);
double Ftmp22 = 15.0*Ftmp1;
double Ftmp23 = -Ftmp21*Ftmp22;
double Ftmp24 = Ftmp1*(Ftmp23 + 3.0);
double Ftmp25 = pow(z, 2);
double Ftmp26 = -Ftmp22*Ftmp25;
double Ftmp27 = Ftmp1*(Ftmp26 + 3.0);
double Ftmp28 = -15.0*Ftmp13;
double Ftmp29 = Ftmp1*(Ftmp28 + 9.0);
double Ftmp30 = -105.0*Ftmp13;
double Ftmp31 = Ftmp30 + 45.0;
double Ftmp32 = Ftmp31*Ftmp5;
double Ftmp33 = y*M[20];
double Ftmp34 = Ftmp32*M[21];
double Ftmp35 = Ftmp5*y;
double Ftmp36 = Ftmp1*Ftmp25;
double Ftmp37 = 3.0*M[27];
double Ftmp38 = Ftmp37*(5.0 - 35.0*Ftmp36);
double Ftmp39 = 105.0*Ftmp1;
double Ftmp40 = -Ftmp21*Ftmp39;
double Ftmp41 = Ftmp40 + 45.0;
double Ftmp42 = Ftmp35*Ftmp41;
double Ftmp43 = 1.0*M[25];
double Ftmp44 = 1.0*Ftmp6;
double Ftmp45 = Ftmp40 + 15.0;
double Ftmp46 = Ftmp45*M[26];
double Ftmp47 = -Ftmp25*Ftmp39;
double Ftmp48 = Ftmp47 + 45.0;
double Ftmp49 = Ftmp44*Ftmp48;
double Ftmp50 = Ftmp24*M[6];
double Ftmp51 = Ftmp27*M[8];
double Ftmp52 = Ftmp32*x;
double Ftmp53 = Ftmp52*y;
double Ftmp54 = Ftmp42*x;
double Ftmp55 = Ftmp45*M[16];
double Ftmp56 = Ftmp6*x;
double Ftmp57 = Ftmp52*z;
double Ftmp58 = Ftmp48*M[18];
double Ftmp59 = y*z;
double Ftmp60 = Ftmp18*Ftmp59;
double Ftmp61 = 315.0*Ftmp1;
double Ftmp62 = -Ftmp25*Ftmp61;
double Ftmp63 = Ftmp62 + 105.0;
double Ftmp64 = 3.0*M[47];
double Ftmp65 = Ftmp63*Ftmp64;
double Ftmp66 = -945.0*Ftmp13;
double Ftmp67 = Ftmp66 + 315.0;
double Ftmp68 = Ftmp67*M[38];
double Ftmp69 = Ftmp1*Ftmp21;
double Ftmp70 = -945.0*Ftmp69;
double Ftmp71 = Ftmp70 + 315.0;
double Ftmp72 = Ftmp71*M[45];
double Ftmp73 = 1.0*Ftmp60;
double Ftmp74 = Ftmp47 + 15.0;
double Ftmp75 = 1.0*Ftmp74*M[17];
double Ftmp76 = 1.0*Ftmp16;
double Ftmp77 = Ftmp45*M[12];
double Ftmp78 = Ftmp74*M[14];
double Ftmp79 = Ftmp18*x;
double Ftmp80 = Ftmp79*y;
double Ftmp81 = Ftmp80*z;
double Ftmp82 = Ftmp71*M[30];
double Ftmp83 = 1.0*Ftmp80;
double Ftmp84 = -945.0*Ftmp36;
double Ftmp85 = Ftmp84 + 315.0;
double Ftmp86 = Ftmp85*z*M[32];
double Ftmp87 = Ftmp19*y;
double Ftmp88 = Ftmp37*(Ftmp62 + 35.0);
double Ftmp89 = Ftmp43*Ftmp71;
double Ftmp90 = Ftmp66 + 525.0;
double Ftmp91 = Ftmp19*Ftmp90;
double Ftmp92 = 1.0*Ftmp19;
double Ftmp93 = Ftmp92*z;
double Ftmp94 = Ftmp70 + 105.0;
double Ftmp95 = Ftmp94*M[26];
double Ftmp96 = Ftmp85*M[28];
double Ftmp97 = z*M[21];
double Ftmp98 = -10395.0*Ftmp13;
double Ftmp99 = pow(R, -8);
double Ftmp100 = Ftmp99*M[38];
double Ftmp101 = Ftmp12*Ftmp59;
double Ftmp102 = Ftmp101*Ftmp99;
double Ftmp103 = -3465.0*Ftmp36;
double Ftmp104 = Ftmp64*(Ftmp103 + 945.0);
double Ftmp105 = -10395.0*Ftmp69;
double Ftmp106 = 1.0*M[45];
double Ftmp107 = Ftmp106*(Ftmp105 + 2835.0);
double Ftmp108 = pow(y, 4);
double Ftmp109 = 945.0*Ftmp5;
double Ftmp110 = Ftmp108*Ftmp109;
double Ftmp111 = 630.0*Ftmp1;
double Ftmp112 = Ftmp5*(Ftmp110 - Ftmp111*Ftmp21 + 45.0);
double Ftmp113 = pow(z, 4);
double Ftmp114 = Ftmp109*Ftmp113;
double Ftmp115 = Ftmp5*(-Ftmp111*Ftmp25 + Ftmp114 + 45.0);
double Ftmp116 = pow(x, 4);
double Ftmp117 = Ftmp109*Ftmp116;
double Ftmp118 = Ftmp5*(Ftmp117 - 1050.0*Ftmp13 + 225.0);
double Ftmp119 = Ftmp112*M[29];
double Ftmp120 = Ftmp115*M[33];
double Ftmp121 = 10395.0*Ftmp5;
double Ftmp122 = Ftmp113*Ftmp121;
double Ftmp123 = Ftmp122 - 5670.0*Ftmp36 + 315.0;
double Ftmp124 = Ftmp123*M[53];
double Ftmp125 = Ftmp116*Ftmp121;
double Ftmp126 = Ftmp125 - 9450.0*Ftmp13 + 1575.0;
double Ftmp127 = Ftmp126*Ftmp80;
double Ftmp128 = Ftmp108*Ftmp121;
double Ftmp129 = Ftmp128 - 9450.0*Ftmp69 + 1575.0;
double Ftmp130 = Ftmp129*M[49];
double Ftmp131 = Ftmp128 - 5670.0*Ftmp69 + 315.0;
double Ftmp132 = Ftmp131*M[50];
double Ftmp133 = Ftmp79*z;
double Ftmp134 = Ftmp126*Ftmp133;
double Ftmp135 = Ftmp122 - 9450.0*Ftmp36 + 1575.0;
double Ftmp136 = Ftmp135*M[54];
double Ftmp137 = Ftmp131*M[44];
double Ftmp138 = Ftmp123*M[48];
double Ftmp139 = Ftmp109*Ftmp21;
double Ftmp140 = Ftmp139*Ftmp25;
double Ftmp141 = Ftmp5*(Ftmp140 + Ftmp45 + Ftmp47);
double Ftmp142 = -Ftmp21*Ftmp61;
double Ftmp143 = Ftmp12*Ftmp139;
double Ftmp144 = Ftmp5*(Ftmp142 + Ftmp143 + Ftmp31);
double Ftmp145 = Ftmp12*Ftmp25;
double Ftmp146 = Ftmp109*Ftmp145;
double Ftmp147 = Ftmp5*(Ftmp146 + Ftmp31 + Ftmp62);
double Ftmp148 = Ftmp121*Ftmp145;
double Ftmp149 = -2835.0*Ftmp36;
double Ftmp150 = Ftmp148 + Ftmp149;
double Ftmp151 = Ftmp80*(Ftmp150 + Ftmp67);
double Ftmp152 = Ftmp121*Ftmp21;
double Ftmp153 = Ftmp152*Ftmp25;
double Ftmp154 = Ftmp149 + Ftmp153;
double Ftmp155 = Ftmp80*(Ftmp154 + Ftmp71);
double Ftmp156 = Ftmp12*Ftmp152;
double Ftmp157 = -2835.0*Ftmp69;
double Ftmp158 = Ftmp156 + Ftmp157;
double Ftmp159 = -2835.0*Ftmp13;
double Ftmp160 = Ftmp159 + 945.0;
double Ftmp161 = Ftmp80*(Ftmp158 + Ftmp160);
double Ftmp162 = Ftmp133*(Ftmp158 + Ftmp67);
double Ftmp163 = Ftmp133*(Ftmp153 + Ftmp157 + Ftmp85);
double Ftmp164 = Ftmp133*(Ftmp150 + Ftmp160);
double Ftmp165 = 4725.0*Ftmp1;
double Ftmp166 = -Ftmp165*Ftmp21;
double Ftmp167 = -Ftmp165*Ftmp25;
double Ftmp168 = x*M[13];
double Ftmp169 = Ftmp21*Ftmp5;
double Ftmp170 = 15.0*x;
double Ftmp171 = Ftmp18*Ftmp21;
double Ftmp172 = Ftmp171*z;
double Ftmp173 = Ftmp1*(Ftmp28 + 3.0);
double Ftmp174 = Ftmp1*(Ftmp23 + 9.0);
double Ftmp175 = Ftmp30 + 15.0;
double Ftmp176 = Ftmp175*M[23];
double Ftmp177 = Ftmp41*M[30];
double Ftmp178 = Ftmp5*x;
double Ftmp179 = Ftmp173*M[3];
double Ftmp180 = Ftmp175*M[11];
double Ftmp181 = Ftmp6*y;
double Ftmp182 = Ftmp181*Ftmp41;
double Ftmp183 = Ftmp175*M[10];
double Ftmp184 = 1.0*Ftmp133;
double Ftmp185 = 1.0*Ftmp35;
double Ftmp186 = Ftmp83*z;
double Ftmp187 = Ftmp171*x;
double Ftmp188 = Ftmp67*x;
double Ftmp189 = Ftmp70 + 525.0;
double Ftmp190 = Ftmp66 + 105.0;
double Ftmp191 = Ftmp190*M[23];
double Ftmp192 = 1.0*Ftmp171;
double Ftmp193 = Ftmp100*x*(Ftmp98 + 2835.0);
double Ftmp194 = Ftmp21*z;
double Ftmp195 = Ftmp99*x;
double Ftmp196 = Ftmp194*Ftmp195;
double Ftmp197 = Ftmp5*(-Ftmp111*Ftmp12 + Ftmp117 + 45.0);
double Ftmp198 = Ftmp5*(Ftmp110 - 1050.0*Ftmp69 + 225.0);
double Ftmp199 = Ftmp197*M[19];
double Ftmp200 = Ftmp125 - 5670.0*Ftmp13 + 315.0;
double Ftmp201 = Ftmp200*M[36];
double Ftmp202 = Ftmp200*M[35];
double Ftmp203 = Ftmp5*(Ftmp146 + Ftmp30 + Ftmp74);
double Ftmp204 = -315.0*Ftmp13;
double Ftmp205 = Ftmp5*(Ftmp143 + Ftmp204 + Ftmp41);
double Ftmp206 = Ftmp5*(Ftmp140 + Ftmp41 + Ftmp62);
double Ftmp207 = 1.0*M[46];
double Ftmp208 = Ftmp60*(Ftmp156 + Ftmp159 + Ftmp71);
double Ftmp209 = Ftmp60*(Ftmp148 + Ftmp159 + Ftmp85);
double Ftmp210 = Ftmp60*(Ftmp154 + Ftmp157 + 945.0);
double Ftmp211 = -4725.0*Ftmp13;
double Ftmp212 = Ftmp25*Ftmp5;
double Ftmp213 = Ftmp18*Ftmp25;
double Ftmp214 = Ftmp213*x;
double Ftmp215 = Ftmp1*(Ftmp26 + 9.0);
double Ftmp216 = 1.0*Ftmp178;
double Ftmp217 = 1.0*Ftmp214;
double Ftmp218 = Ftmp84 + 525.0;
double Ftmp219 = Ftmp213*y;
double Ftmp220 = Ftmp25*y;
double Ftmp221 = Ftmp195*Ftmp220;
double Ftmp222 = Ftmp5*(Ftmp114 - 1050.0*Ftmp36 + 225.0);
double Ftmp223 = Ftmp5*(Ftmp143 + Ftmp30 + Ftmp45);
double Ftmp224 = Ftmp5*(Ftmp146 + Ftmp204 + Ftmp48);
double Ftmp225 = Ftmp5*(Ftmp140 + Ftmp142 + Ftmp48);
#pragma omp atomic
F[0] += Ftmp0*(-Ftmp10*M[1] + Ftmp100*Ftmp101*(Ftmp98 + 4725.0) + Ftmp102*Ftmp104 + Ftmp102*Ftmp107 - Ftmp11*x + Ftmp112*M[44] + Ftmp115*M[48] + Ftmp118*x*M[19] + Ftmp118*M[34] + Ftmp119*x + Ftmp120*x - Ftmp124*Ftmp83 - Ftmp127*M[35] - 3.0*Ftmp13*M[0] - Ftmp130*Ftmp80 - Ftmp132*Ftmp133 - Ftmp133*Ftmp136 - Ftmp134*M[36] - Ftmp137*Ftmp92 - Ftmp138*Ftmp92 + Ftmp14*Ftmp15 + Ftmp141*x*M[31] + Ftmp141*M[46] + Ftmp144*x*M[22] + Ftmp144*M[37] + Ftmp147*x*M[24] + Ftmp147*M[39] - Ftmp151*M[42] - Ftmp155*M[51] + Ftmp16*(Ftmp30 + 75.0)*M[9] - Ftmp161*M[40] - Ftmp162*M[41] - Ftmp163*M[52] - Ftmp164*M[43] + Ftmp17*Ftmp3 + Ftmp17*z*M[5] - Ftmp19*Ftmp20*z - Ftmp19*(Ftmp125 - 13230.0*Ftmp13 + 3675.0)*M[34] - Ftmp19*(Ftmp148 + Ftmp167 + Ftmp90)*M[39] - Ftmp19*(Ftmp156 + Ftmp166 + Ftmp90)*M[37] - Ftmp2*Ftmp3 - Ftmp24*M[12] - Ftmp27*M[14] - Ftmp29*x*M[3] - Ftmp29*M[9] + Ftmp32*Ftmp33 - Ftmp33*Ftmp91 + Ftmp34*z + Ftmp35*Ftmp38 + Ftmp35*Ftmp75*x - Ftmp4*M[5] + Ftmp42*Ftmp43 + Ftmp44*Ftmp46 + Ftmp49*M[28] - Ftmp50*x - Ftmp51*x + Ftmp53*M[10] + Ftmp54*M[15] + Ftmp55*Ftmp56 + Ftmp56*Ftmp58 + Ftmp57*M[11] - Ftmp60*Ftmp65 - Ftmp60*Ftmp68 - Ftmp67*Ftmp81*M[23] + Ftmp7*Ftmp8 - Ftmp72*Ftmp73 + Ftmp76*Ftmp77 + Ftmp76*Ftmp78 - Ftmp81*Ftmp82 - Ftmp83*Ftmp86 - Ftmp87*Ftmp88 - Ftmp87*Ftmp89 - Ftmp91*Ftmp97 - Ftmp92*(Ftmp153 + Ftmp84 + Ftmp94)*M[46] - Ftmp93*Ftmp95 - Ftmp93*Ftmp96 + 1.0*M[0]);
#pragma omp atomic
F[1] += Ftmp0*(-Ftmp10*M[0] + Ftmp104*Ftmp196 + Ftmp106*Ftmp196*(Ftmp105 + 4725.0) - Ftmp11*y + Ftmp115*M[53] + Ftmp120*y - Ftmp124*Ftmp192 - Ftmp127*M[34] - Ftmp129*Ftmp60*M[50] - Ftmp129*Ftmp83*M[44] - Ftmp133*Ftmp65 - Ftmp133*Ftmp68 - Ftmp136*Ftmp60 - Ftmp138*Ftmp83 + Ftmp15*y*M[5] - Ftmp151*M[39] - Ftmp155*Ftmp207 - Ftmp161*M[37] - 105.0*Ftmp168*Ftmp172 + Ftmp168*Ftmp7 + Ftmp169*Ftmp170*M[4] + Ftmp169*Ftmp183 + Ftmp169*Ftmp75 + Ftmp169*(Ftmp40 + 75.0)*M[15] - Ftmp171*Ftmp188*M[20] - Ftmp171*Ftmp202 - Ftmp171*(Ftmp128 - 13230.0*Ftmp69 + 3675.0)*M[49] - Ftmp171*(Ftmp148 + Ftmp190 + Ftmp84)*M[42] - Ftmp171*(Ftmp153 + Ftmp167 + Ftmp189)*M[51] - Ftmp171*(Ftmp156 + Ftmp189 + Ftmp211)*M[40] - Ftmp172*Ftmp189*M[30] - Ftmp172*Ftmp191 - Ftmp173*M[10] - Ftmp174*y*M[6] - Ftmp174*M[15] + Ftmp176*Ftmp6 + Ftmp177*Ftmp6 + Ftmp178*Ftmp38 + Ftmp178*Ftmp41*Ftmp43 - Ftmp179*y + Ftmp180*Ftmp181 + Ftmp181*Ftmp58 + Ftmp182*M[16] - Ftmp184*Ftmp72 + Ftmp185*Ftmp78*x - Ftmp186*Ftmp71*M[26] - Ftmp186*Ftmp96 - Ftmp187*Ftmp189*Ftmp43 - Ftmp187*Ftmp88 - Ftmp192*Ftmp86 + Ftmp193*Ftmp194 + Ftmp197*M[35] + Ftmp198*y*M[29] + Ftmp198*M[49] + Ftmp199*y - Ftmp201*Ftmp60 + Ftmp203*y*M[24] + Ftmp203*M[42] + Ftmp205*y*M[22] + Ftmp205*M[40] + Ftmp206*y*M[31] + Ftmp206*M[51] - Ftmp208*M[41] - Ftmp209*M[43] + Ftmp21*Ftmp7*M[7] - Ftmp210*M[52] - Ftmp27*M[17] - Ftmp4*M[7] + Ftmp49*M[32] - Ftmp51*y + Ftmp52*M[20] + Ftmp53*M[9] + 1.0*Ftmp54*M[12] - Ftmp67*Ftmp80*Ftmp97 - 3.0*Ftmp69*M[1] - Ftmp9*M[4] + 1.0*M[1]);
#pragma omp atomic
F[2] += Ftmp0*(Ftmp107*Ftmp221 + Ftmp112*M[50] + Ftmp119*z - Ftmp130*Ftmp60 - Ftmp132*Ftmp213 - Ftmp133*Ftmp33*Ftmp67 - Ftmp134*M[34] - Ftmp135*Ftmp184*M[48] - Ftmp135*Ftmp73*M[53] - Ftmp137*Ftmp184 - Ftmp14*Ftmp2 + 15.0*Ftmp14*Ftmp212 + Ftmp15*Ftmp3 - Ftmp162*M[37] - Ftmp163*Ftmp207 - Ftmp164*M[39] + 15.0*Ftmp168*Ftmp35 + Ftmp170*Ftmp212*M[5] - Ftmp173*M[11] + Ftmp176*Ftmp35 + Ftmp177*Ftmp35 - Ftmp179*z + Ftmp180*Ftmp212 + Ftmp181*Ftmp183 + Ftmp182*M[15] + Ftmp185*Ftmp48*M[32] - Ftmp188*Ftmp213*M[21] - Ftmp191*Ftmp219 + Ftmp193*Ftmp220 + Ftmp197*M[36] + Ftmp199*z - Ftmp20*Ftmp214 - Ftmp201*Ftmp213 - Ftmp202*Ftmp60 - Ftmp208*M[40] - Ftmp209*M[42] - Ftmp210*M[51] + Ftmp212*Ftmp55 + Ftmp212*(Ftmp47 + 75.0)*M[18] - Ftmp213*(Ftmp122 - 13230.0*Ftmp36 + 3675.0)*M[54] - Ftmp213*(Ftmp148 + Ftmp211 + Ftmp218)*M[43] - Ftmp213*(Ftmp153 + Ftmp166 + Ftmp218)*M[52] - Ftmp213*(Ftmp156 + Ftmp66 + Ftmp94)*M[41] - Ftmp215*z*M[8] - Ftmp215*M[18] + Ftmp216*Ftmp46 + Ftmp216*Ftmp48*M[28] - Ftmp217*Ftmp218*M[28] - Ftmp217*Ftmp95 - 1.0*Ftmp218*Ftmp219*M[32] - Ftmp219*Ftmp82 + Ftmp221*Ftmp64*(Ftmp103 + 1575.0) + Ftmp222*z*M[33] + Ftmp222*M[54] + Ftmp223*z*M[22] + Ftmp223*M[41] + Ftmp224*z*M[24] + Ftmp224*M[43] + Ftmp225*z*M[31] + Ftmp225*M[52] - Ftmp24*M[16] + Ftmp34*x - 3.0*Ftmp36*M[2] - Ftmp37*Ftmp63*Ftmp81 - Ftmp4*x*M[0] - Ftmp4*y*M[1] + Ftmp44*Ftmp77*x + Ftmp49*x*M[14] + Ftmp49*y*M[17] - Ftmp50*z + Ftmp57*M[9] - Ftmp65*Ftmp80 - Ftmp68*Ftmp80 - Ftmp72*Ftmp83 - Ftmp81*Ftmp89 - Ftmp9*M[5] + 1.0*M[2]);

}

void P2M_6(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = q*z;
double Mtmp3 = pow(x, 2);
double Mtmp4 = (1.0/2.0)*q;
double Mtmp5 = Mtmp0*y;
double Mtmp6 = Mtmp0*z;
double Mtmp7 = pow(y, 2);
double Mtmp8 = Mtmp1*z;
double Mtmp9 = pow(z, 2);
double Mtmp10 = pow(x, 3);
double Mtmp11 = (1.0/6.0)*q;
double Mtmp12 = (1.0/2.0)*Mtmp3;
double Mtmp13 = (1.0/2.0)*Mtmp0;
double Mtmp14 = pow(y, 3);
double Mtmp15 = (1.0/2.0)*Mtmp7;
double Mtmp16 = (1.0/2.0)*Mtmp9;
double Mtmp17 = pow(z, 3);
double Mtmp18 = pow(x, 4);
double Mtmp19 = (1.0/24.0)*q;
double Mtmp20 = (1.0/6.0)*Mtmp10;
double Mtmp21 = Mtmp7*q;
double Mtmp22 = (1.0/4.0)*Mtmp3;
double Mtmp23 = Mtmp9*q;
double Mtmp24 = (1.0/6.0)*Mtmp0;
double Mtmp25 = pow(y, 4);
double Mtmp26 = (1.0/6.0)*Mtmp14;
double Mtmp27 = (1.0/4.0)*Mtmp9;
double Mtmp28 = (1.0/6.0)*Mtmp17;
double Mtmp29 = pow(z, 4);
double Mtmp30 = pow(x, 5);
double Mtmp31 = (1.0/120.0)*q;
double Mtmp32 = (1.0/24.0)*Mtmp18;
double Mtmp33 = (1.0/12.0)*Mtmp10;
double Mtmp34 = (1.0/12.0)*Mtmp14;
double Mtmp35 = Mtmp3*q;
double Mtmp36 = Mtmp2*Mtmp7;
double Mtmp37 = Mtmp1*Mtmp9;
double Mtmp38 = (1.0/12.0)*Mtmp17;
double Mtmp39 = (1.0/24.0)*Mtmp0;
double Mtmp40 = Mtmp0*Mtmp7;
double Mtmp41 = pow(y, 5);
double Mtmp42 = (1.0/24.0)*Mtmp25;
double Mtmp43 = (1.0/24.0)*Mtmp29;
double Mtmp44 = pow(z, 5);
double Mtmp45 = (1.0/720.0)*q;
double Mtmp46 = (1.0/120.0)*Mtmp30;
double Mtmp47 = (1.0/48.0)*Mtmp18;
double Mtmp48 = (1.0/36.0)*Mtmp10*q;
double Mtmp49 = (1.0/48.0)*Mtmp35;
double Mtmp50 = (1.0/120.0)*Mtmp0;
M[0] += -Mtmp0;
M[1] += -Mtmp1;
M[2] += -Mtmp2;
M[3] += Mtmp3*Mtmp4;
M[4] += Mtmp5;
M[5] += Mtmp6;
M[6] += Mtmp4*Mtmp7;
M[7] += Mtmp8;
M[8] += Mtmp4*Mtmp9;
M[9] += -Mtmp10*Mtmp11;
M[10] += -Mtmp1*Mtmp12;
M[11] += -Mtmp12*Mtmp2;
M[12] += -Mtmp13*Mtmp7;
M[13] += -Mtmp5*z;
M[14] += -Mtmp13*Mtmp9;
M[15] += -Mtmp11*Mtmp14;
M[16] += -Mtmp15*Mtmp2;
M[17] += -Mtmp1*Mtmp16;
M[18] += -Mtmp11*Mtmp17;
M[19] += Mtmp18*Mtmp19;
M[20] += Mtmp1*Mtmp20;
M[21] += Mtmp2*Mtmp20;
M[22] += Mtmp21*Mtmp22;
M[23] += Mtmp12*Mtmp8;
M[24] += Mtmp22*Mtmp23;
M[25] += Mtmp14*Mtmp24;
M[26] += Mtmp15*Mtmp6;
M[27] += Mtmp16*Mtmp5;
M[28] += Mtmp17*Mtmp24;
M[29] += Mtmp19*Mtmp25;
M[30] += Mtmp2*Mtmp26;
M[31] += Mtmp21*Mtmp27;
M[32] += Mtmp1*Mtmp28;
M[33] += Mtmp19*Mtmp29;
M[34] += -Mtmp30*Mtmp31;
M[35] += -Mtmp1*Mtmp32;
M[36] += -Mtmp2*Mtmp32;
M[37] += -Mtmp21*Mtmp33;
M[38] += -Mtmp20*Mtmp8;
M[39] += -Mtmp23*Mtmp33;
M[40] += -Mtmp34*Mtmp35;
M[41] += -Mtmp22*Mtmp36;
M[42] += -Mtmp22*Mtmp37;
M[43] += -Mtmp35*Mtmp38;
M[44] += -Mtmp25*Mtmp39;
M[45] += -Mtmp26*Mtmp6;
M[46] += -Mtmp27*Mtmp40;
M[47] += -Mtmp28*Mtmp5;
M[48] += -Mtmp29*Mtmp39;
M[49] += -Mtmp31*Mtmp41;
M[50] += -Mtmp2*Mtmp42;
M[51] += -Mtmp23*Mtmp34;
M[52] += -Mtmp21*Mtmp38;
M[53] += -Mtmp1*Mtmp43;
M[54] += -Mtmp31*Mtmp44;
M[55] += Mtmp45*pow(x, 6);
M[56] += Mtmp1*Mtmp46;
M[57] += Mtmp2*Mtmp46;
M[58] += Mtmp21*Mtmp47;
M[59] += Mtmp32*Mtmp8;
M[60] += Mtmp23*Mtmp47;
M[61] += Mtmp14*Mtmp48;
M[62] += Mtmp33*Mtmp36;
M[63] += Mtmp33*Mtmp37;
M[64] += Mtmp17*Mtmp48;
M[65] += Mtmp25*Mtmp49;
M[66] += Mtmp2*Mtmp3*Mtmp34;
M[67] += (1.0/8.0)*Mtmp21*Mtmp3*Mtmp9;
M[68] += Mtmp1*Mtmp3*Mtmp38;
M[69] += Mtmp29*Mtmp49;
M[70] += Mtmp41*Mtmp50;
M[71] += Mtmp42*Mtmp6;
M[72] += Mtmp0*Mtmp34*Mtmp9;
M[73] += Mtmp38*Mtmp40;
M[74] += Mtmp43*Mtmp5;
M[75] += Mtmp44*Mtmp50;
M[76] += Mtmp45*pow(y, 6);
M[77] += (1.0/120.0)*Mtmp2*Mtmp41;
M[78] += (1.0/48.0)*Mtmp23*Mtmp25;
M[79] += (1.0/36.0)*Mtmp14*Mtmp17*q;
M[80] += (1.0/48.0)*Mtmp21*Mtmp29;
M[81] += (1.0/120.0)*Mtmp1*Mtmp44;
M[82] += Mtmp45*pow(z, 6);

}
void M2M_6(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = x*M[1];
double Mstmp2 = y*M[0];
double Mstmp3 = x*M[2];
double Mstmp4 = z*M[0];
double Mstmp5 = y*M[1];
double Mstmp6 = y*M[2];
double Mstmp7 = z*M[1];
double Mstmp8 = z*M[2];
double Mstmp9 = x*M[3];
double Mstmp10 = pow(x, 2);
double Mstmp11 = (1.0/2.0)*Mstmp10;
double Mstmp12 = x*M[4];
double Mstmp13 = y*M[3];
double Mstmp14 = Mstmp0*y;
double Mstmp15 = x*M[5];
double Mstmp16 = z*M[3];
double Mstmp17 = Mstmp0*z;
double Mstmp18 = x*M[6];
double Mstmp19 = y*M[4];
double Mstmp20 = Mstmp1*y;
double Mstmp21 = pow(y, 2);
double Mstmp22 = (1.0/2.0)*M[0];
double Mstmp23 = x*M[7];
double Mstmp24 = y*M[5];
double Mstmp25 = z*M[4];
double Mstmp26 = Mstmp3*y;
double Mstmp27 = Mstmp1*z;
double Mstmp28 = Mstmp2*z;
double Mstmp29 = x*M[8];
double Mstmp30 = z*M[5];
double Mstmp31 = Mstmp3*z;
double Mstmp32 = pow(z, 2);
double Mstmp33 = y*M[6];
double Mstmp34 = (1.0/2.0)*Mstmp21;
double Mstmp35 = y*M[7];
double Mstmp36 = z*M[6];
double Mstmp37 = Mstmp5*z;
double Mstmp38 = y*M[8];
double Mstmp39 = z*M[7];
double Mstmp40 = Mstmp6*z;
double Mstmp41 = (1.0/2.0)*Mstmp32;
double Mstmp42 = z*M[8];
double Mstmp43 = x*M[9];
double Mstmp44 = pow(x, 3);
double Mstmp45 = (1.0/6.0)*Mstmp44;
double Mstmp46 = x*M[10];
double Mstmp47 = y*M[9];
double Mstmp48 = Mstmp9*y;
double Mstmp49 = x*M[11];
double Mstmp50 = z*M[9];
double Mstmp51 = Mstmp9*z;
double Mstmp52 = x*M[12];
double Mstmp53 = y*M[10];
double Mstmp54 = Mstmp12*y;
double Mstmp55 = x*M[13];
double Mstmp56 = y*M[11];
double Mstmp57 = z*M[10];
double Mstmp58 = Mstmp15*y;
double Mstmp59 = Mstmp12*z;
double Mstmp60 = Mstmp13*z;
double Mstmp61 = x*M[14];
double Mstmp62 = z*M[11];
double Mstmp63 = Mstmp15*z;
double Mstmp64 = x*M[15];
double Mstmp65 = y*M[12];
double Mstmp66 = Mstmp18*y;
double Mstmp67 = pow(y, 3);
double Mstmp68 = (1.0/6.0)*M[0];
double Mstmp69 = x*M[16];
double Mstmp70 = y*M[13];
double Mstmp71 = z*M[12];
double Mstmp72 = Mstmp23*y;
double Mstmp73 = Mstmp18*z;
double Mstmp74 = Mstmp19*z;
double Mstmp75 = x*M[17];
double Mstmp76 = y*M[14];
double Mstmp77 = z*M[13];
double Mstmp78 = Mstmp29*y;
double Mstmp79 = Mstmp23*z;
double Mstmp80 = Mstmp24*z;
double Mstmp81 = x*M[18];
double Mstmp82 = z*M[14];
double Mstmp83 = Mstmp29*z;
double Mstmp84 = pow(z, 3);
double Mstmp85 = y*M[15];
double Mstmp86 = (1.0/6.0)*Mstmp67;
double Mstmp87 = y*M[16];
double Mstmp88 = z*M[15];
double Mstmp89 = Mstmp33*z;
double Mstmp90 = y*M[17];
double Mstmp91 = z*M[16];
double Mstmp92 = Mstmp35*z;
double Mstmp93 = y*M[18];
double Mstmp94 = z*M[17];
double Mstmp95 = Mstmp38*z;
double Mstmp96 = (1.0/6.0)*Mstmp84;
double Mstmp97 = z*M[18];
double Mstmp98 = x*M[19];
double Mstmp99 = (1.0/24.0)*pow(x, 4);
double Mstmp100 = x*M[20];
double Mstmp101 = y*M[19];
double Mstmp102 = Mstmp43*y;
double Mstmp103 = x*M[21];
double Mstmp104 = x*M[22];
double Mstmp105 = y*M[20];
double Mstmp106 = Mstmp46*y;
double Mstmp107 = (1.0/4.0)*Mstmp10;
double Mstmp108 = Mstmp21*M[0];
double Mstmp109 = x*M[23];
double Mstmp110 = y*M[21];
double Mstmp111 = Mstmp49*y;
double Mstmp112 = x*M[24];
double Mstmp113 = Mstmp107*Mstmp32;
double Mstmp114 = x*M[25];
double Mstmp115 = y*M[22];
double Mstmp116 = Mstmp52*y;
double Mstmp117 = Mstmp107*Mstmp21;
double Mstmp118 = x*M[26];
double Mstmp119 = y*M[23];
double Mstmp120 = Mstmp55*y;
double Mstmp121 = x*M[27];
double Mstmp122 = y*M[24];
double Mstmp123 = Mstmp61*y;
double Mstmp124 = x*M[28];
double Mstmp125 = x*M[29];
double Mstmp126 = y*M[25];
double Mstmp127 = Mstmp64*y;
double Mstmp128 = pow(y, 4);
double Mstmp129 = (1.0/24.0)*M[0];
double Mstmp130 = x*M[30];
double Mstmp131 = y*M[26];
double Mstmp132 = Mstmp69*y;
double Mstmp133 = x*M[31];
double Mstmp134 = y*M[27];
double Mstmp135 = Mstmp75*y;
double Mstmp136 = (1.0/4.0)*Mstmp32;
double Mstmp137 = x*M[32];
double Mstmp138 = y*M[28];
double Mstmp139 = Mstmp81*y;
double Mstmp140 = x*M[33];
double Mstmp141 = pow(z, 4);
double Mstmp142 = y*M[29];
double Mstmp143 = (1.0/24.0)*Mstmp128;
double Mstmp144 = y*M[30];
double Mstmp145 = y*M[31];
double Mstmp146 = Mstmp136*Mstmp21;
double Mstmp147 = y*M[32];
double Mstmp148 = y*M[33];
double Mstmp149 = (1.0/24.0)*Mstmp141;
double Mstmp150 = (1.0/120.0)*pow(x, 5);
double Mstmp151 = (1.0/12.0)*Mstmp44;
double Mstmp152 = Mstmp151*Mstmp32;
double Mstmp153 = (1.0/12.0)*Mstmp10;
double Mstmp154 = Mstmp153*M[0];
double Mstmp155 = Mstmp151*Mstmp21;
double Mstmp156 = Mstmp153*Mstmp67;
double Mstmp157 = Mstmp153*Mstmp84;
double Mstmp158 = pow(y, 5);
double Mstmp159 = (1.0/120.0)*M[0];
double Mstmp160 = (1.0/12.0)*Mstmp32*Mstmp67;
double Mstmp161 = (1.0/12.0)*Mstmp84;
double Mstmp162 = pow(z, 5);
double Mstmp163 = (1.0/120.0)*Mstmp158;
double Mstmp164 = Mstmp161*Mstmp21;
double Mstmp165 = (1.0/120.0)*Mstmp162;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += Mstmp0 + M[3];
#pragma omp atomic
Ms[4] += Mstmp1 + Mstmp2 + M[4];
#pragma omp atomic
Ms[5] += Mstmp3 + Mstmp4 + M[5];
#pragma omp atomic
Ms[6] += Mstmp5 + M[6];
#pragma omp atomic
Ms[7] += Mstmp6 + Mstmp7 + M[7];
#pragma omp atomic
Ms[8] += Mstmp8 + M[8];
#pragma omp atomic
Ms[9] += Mstmp11*M[0] + Mstmp9 + M[9];
#pragma omp atomic
Ms[10] += Mstmp11*M[1] + Mstmp12 + Mstmp13 + Mstmp14 + M[10];
#pragma omp atomic
Ms[11] += Mstmp11*M[2] + Mstmp15 + Mstmp16 + Mstmp17 + M[11];
#pragma omp atomic
Ms[12] += Mstmp18 + Mstmp19 + Mstmp20 + Mstmp21*Mstmp22 + M[12];
#pragma omp atomic
Ms[13] += Mstmp23 + Mstmp24 + Mstmp25 + Mstmp26 + Mstmp27 + Mstmp28 + M[13];
#pragma omp atomic
Ms[14] += Mstmp22*Mstmp32 + Mstmp29 + Mstmp30 + Mstmp31 + M[14];
#pragma omp atomic
Ms[15] += Mstmp33 + Mstmp34*M[1] + M[15];
#pragma omp atomic
Ms[16] += Mstmp34*M[2] + Mstmp35 + Mstmp36 + Mstmp37 + M[16];
#pragma omp atomic
Ms[17] += Mstmp38 + Mstmp39 + Mstmp40 + Mstmp41*M[1] + M[17];
#pragma omp atomic
Ms[18] += Mstmp41*M[2] + Mstmp42 + M[18];
#pragma omp atomic
Ms[19] += Mstmp11*M[3] + Mstmp43 + Mstmp45*M[0] + M[19];
#pragma omp atomic
Ms[20] += Mstmp11*Mstmp2 + Mstmp11*M[4] + Mstmp45*M[1] + Mstmp46 + Mstmp47 + Mstmp48 + M[20];
#pragma omp atomic
Ms[21] += Mstmp11*Mstmp4 + Mstmp11*M[5] + Mstmp45*M[2] + Mstmp49 + Mstmp50 + Mstmp51 + M[21];
#pragma omp atomic
Ms[22] += Mstmp0*Mstmp34 + Mstmp11*Mstmp5 + Mstmp11*M[6] + Mstmp34*M[3] + Mstmp52 + Mstmp53 + Mstmp54 + M[22];
#pragma omp atomic
Ms[23] += Mstmp11*Mstmp6 + Mstmp11*Mstmp7 + Mstmp11*M[7] + Mstmp14*z + Mstmp55 + Mstmp56 + Mstmp57 + Mstmp58 + Mstmp59 + Mstmp60 + M[23];
#pragma omp atomic
Ms[24] += Mstmp0*Mstmp41 + Mstmp11*Mstmp8 + Mstmp11*M[8] + Mstmp41*M[3] + Mstmp61 + Mstmp62 + Mstmp63 + M[24];
#pragma omp atomic
Ms[25] += Mstmp1*Mstmp34 + Mstmp34*M[4] + Mstmp64 + Mstmp65 + Mstmp66 + Mstmp67*Mstmp68 + M[25];
#pragma omp atomic
Ms[26] += Mstmp20*z + Mstmp3*Mstmp34 + Mstmp34*Mstmp4 + Mstmp34*M[5] + Mstmp69 + Mstmp70 + Mstmp71 + Mstmp72 + Mstmp73 + Mstmp74 + M[26];
#pragma omp atomic
Ms[27] += Mstmp1*Mstmp41 + Mstmp2*Mstmp41 + Mstmp26*z + Mstmp41*M[4] + Mstmp75 + Mstmp76 + Mstmp77 + Mstmp78 + Mstmp79 + Mstmp80 + M[27];
#pragma omp atomic
Ms[28] += Mstmp3*Mstmp41 + Mstmp41*M[5] + Mstmp68*Mstmp84 + Mstmp81 + Mstmp82 + Mstmp83 + M[28];
#pragma omp atomic
Ms[29] += Mstmp34*M[6] + Mstmp85 + Mstmp86*M[1] + M[29];
#pragma omp atomic
Ms[30] += Mstmp34*Mstmp7 + Mstmp34*M[7] + Mstmp86*M[2] + Mstmp87 + Mstmp88 + Mstmp89 + M[30];
#pragma omp atomic
Ms[31] += Mstmp34*Mstmp8 + Mstmp34*M[8] + Mstmp41*Mstmp5 + Mstmp41*M[6] + Mstmp90 + Mstmp91 + Mstmp92 + M[31];
#pragma omp atomic
Ms[32] += Mstmp41*Mstmp6 + Mstmp41*M[7] + Mstmp93 + Mstmp94 + Mstmp95 + Mstmp96*M[1] + M[32];
#pragma omp atomic
Ms[33] += Mstmp41*M[8] + Mstmp96*M[2] + Mstmp97 + M[33];
#pragma omp atomic
Ms[34] += Mstmp11*M[9] + Mstmp45*M[3] + Mstmp98 + Mstmp99*M[0] + M[34];
#pragma omp atomic
Ms[35] += Mstmp100 + Mstmp101 + Mstmp102 + Mstmp11*Mstmp13 + Mstmp11*M[10] + Mstmp2*Mstmp45 + Mstmp45*M[4] + Mstmp99*M[1] + M[35];
#pragma omp atomic
Ms[36] += Mstmp103 + Mstmp11*Mstmp16 + Mstmp11*M[11] + Mstmp4*Mstmp45 + Mstmp43*z + Mstmp45*M[5] + Mstmp99*M[2] + z*M[19] + M[36];
#pragma omp atomic
Ms[37] += Mstmp104 + Mstmp105 + Mstmp106 + Mstmp107*Mstmp108 + Mstmp11*Mstmp19 + Mstmp11*M[12] + Mstmp34*Mstmp9 + Mstmp34*M[9] + Mstmp45*Mstmp5 + Mstmp45*M[6] + M[37];
#pragma omp atomic
Ms[38] += Mstmp109 + Mstmp11*Mstmp24 + Mstmp11*Mstmp25 + Mstmp11*Mstmp28 + Mstmp11*M[13] + Mstmp110 + Mstmp111 + Mstmp45*Mstmp6 + Mstmp45*Mstmp7 + Mstmp45*M[7] + Mstmp46*z + Mstmp47*z + Mstmp48*z + z*M[20] + M[38];
#pragma omp atomic
Ms[39] += Mstmp11*Mstmp30 + Mstmp11*M[14] + Mstmp112 + Mstmp113*M[0] + Mstmp41*Mstmp9 + Mstmp41*M[9] + Mstmp45*Mstmp8 + Mstmp45*M[8] + Mstmp49*z + z*M[21] + M[39];
#pragma omp atomic
Ms[40] += Mstmp0*Mstmp86 + Mstmp11*Mstmp33 + Mstmp11*M[15] + Mstmp114 + Mstmp115 + Mstmp116 + Mstmp117*M[1] + Mstmp12*Mstmp34 + Mstmp34*M[10] + Mstmp86*M[3] + M[40];
#pragma omp atomic
Ms[41] += Mstmp11*Mstmp35 + Mstmp11*Mstmp36 + Mstmp11*Mstmp37 + Mstmp11*M[16] + Mstmp117*M[2] + Mstmp118 + Mstmp119 + Mstmp120 + Mstmp15*Mstmp34 + Mstmp16*Mstmp34 + Mstmp17*Mstmp34 + Mstmp34*M[11] + Mstmp52*z + Mstmp53*z + Mstmp54*z + z*M[22] + M[41];
#pragma omp atomic
Ms[42] += Mstmp11*Mstmp38 + Mstmp11*Mstmp39 + Mstmp11*Mstmp40 + Mstmp11*M[17] + Mstmp113*M[1] + Mstmp12*Mstmp41 + Mstmp121 + Mstmp122 + Mstmp123 + Mstmp13*Mstmp41 + Mstmp14*Mstmp41 + Mstmp41*M[10] + Mstmp55*z + Mstmp56*z + Mstmp58*z + z*M[23] + M[42];
#pragma omp atomic
Ms[43] += Mstmp0*Mstmp96 + Mstmp11*Mstmp42 + Mstmp11*M[18] + Mstmp113*M[2] + Mstmp124 + Mstmp15*Mstmp41 + Mstmp41*M[11] + Mstmp61*z + Mstmp96*M[3] + z*M[24] + M[43];
#pragma omp atomic
Ms[44] += Mstmp1*Mstmp86 + Mstmp125 + Mstmp126 + Mstmp127 + Mstmp128*Mstmp129 + Mstmp18*Mstmp34 + Mstmp34*M[12] + Mstmp86*M[4] + M[44];
#pragma omp atomic
Ms[45] += Mstmp130 + Mstmp131 + Mstmp132 + Mstmp23*Mstmp34 + Mstmp25*Mstmp34 + Mstmp27*Mstmp34 + Mstmp3*Mstmp86 + Mstmp34*M[13] + Mstmp4*Mstmp86 + Mstmp64*z + Mstmp65*z + Mstmp66*z + Mstmp86*M[5] + z*M[25] + M[45];
#pragma omp atomic
Ms[46] += Mstmp108*Mstmp136 + Mstmp133 + Mstmp134 + Mstmp135 + Mstmp18*Mstmp41 + Mstmp19*Mstmp41 + Mstmp20*Mstmp41 + Mstmp29*Mstmp34 + Mstmp30*Mstmp34 + Mstmp31*Mstmp34 + Mstmp34*M[14] + Mstmp41*M[12] + Mstmp69*z + Mstmp70*z + Mstmp72*z + z*M[26] + M[46];
#pragma omp atomic
Ms[47] += Mstmp1*Mstmp96 + Mstmp137 + Mstmp138 + Mstmp139 + Mstmp2*Mstmp96 + Mstmp23*Mstmp41 + Mstmp24*Mstmp41 + Mstmp26*Mstmp41 + Mstmp41*M[13] + Mstmp75*z + Mstmp76*z + Mstmp78*z + Mstmp96*M[4] + z*M[27] + M[47];
#pragma omp atomic
Ms[48] += Mstmp129*Mstmp141 + Mstmp140 + Mstmp29*Mstmp41 + Mstmp3*Mstmp96 + Mstmp41*M[14] + Mstmp81*z + Mstmp96*M[5] + z*M[28] + M[48];
#pragma omp atomic
Ms[49] += Mstmp142 + Mstmp143*M[1] + Mstmp34*M[15] + Mstmp86*M[6] + M[49];
#pragma omp atomic
Ms[50] += Mstmp143*M[2] + Mstmp144 + Mstmp34*Mstmp36 + Mstmp34*M[16] + Mstmp7*Mstmp86 + Mstmp85*z + Mstmp86*M[7] + z*M[29] + M[50];
#pragma omp atomic
Ms[51] += Mstmp145 + Mstmp146*M[1] + Mstmp33*Mstmp41 + Mstmp34*Mstmp39 + Mstmp34*M[17] + Mstmp41*M[15] + Mstmp8*Mstmp86 + Mstmp86*M[8] + Mstmp87*z + z*M[30] + M[51];
#pragma omp atomic
Ms[52] += Mstmp146*M[2] + Mstmp147 + Mstmp34*Mstmp42 + Mstmp34*M[18] + Mstmp35*Mstmp41 + Mstmp41*M[16] + Mstmp5*Mstmp96 + Mstmp90*z + Mstmp96*M[6] + z*M[31] + M[52];
#pragma omp atomic
Ms[53] += Mstmp148 + Mstmp149*M[1] + Mstmp38*Mstmp41 + Mstmp41*M[17] + Mstmp6*Mstmp96 + Mstmp93*z + Mstmp96*M[7] + z*M[32] + M[53];
#pragma omp atomic
Ms[54] += Mstmp149*M[2] + Mstmp41*M[18] + Mstmp96*M[8] + z*M[33] + M[54];
#pragma omp atomic
Ms[55] += Mstmp11*M[19] + Mstmp150*M[0] + Mstmp45*M[9] + Mstmp99*M[3] + x*M[34] + M[55];
#pragma omp atomic
Ms[56] += Mstmp11*Mstmp47 + Mstmp11*M[20] + Mstmp13*Mstmp45 + Mstmp150*M[1] + Mstmp2*Mstmp99 + Mstmp45*M[10] + Mstmp98*y + Mstmp99*M[4] + x*M[35] + y*M[34] + M[56];
#pragma omp atomic
Ms[57] += Mstmp11*Mstmp50 + Mstmp11*M[21] + Mstmp150*M[2] + Mstmp16*Mstmp45 + Mstmp4*Mstmp99 + Mstmp45*M[11] + Mstmp98*z + Mstmp99*M[5] + x*M[36] + z*M[34] + M[57];
#pragma omp atomic
Ms[58] += Mstmp100*y + Mstmp108*Mstmp151 + Mstmp11*Mstmp53 + Mstmp11*M[22] + Mstmp117*M[3] + Mstmp19*Mstmp45 + Mstmp34*Mstmp43 + Mstmp34*M[19] + Mstmp45*M[12] + Mstmp5*Mstmp99 + Mstmp99*M[6] + x*M[37] + y*M[35] + M[58];
#pragma omp atomic
Ms[59] += Mstmp100*z + Mstmp101*z + Mstmp102*z + Mstmp103*y + Mstmp11*Mstmp56 + Mstmp11*Mstmp57 + Mstmp11*Mstmp60 + Mstmp11*M[23] + Mstmp24*Mstmp45 + Mstmp25*Mstmp45 + Mstmp28*Mstmp45 + Mstmp45*M[13] + Mstmp6*Mstmp99 + Mstmp7*Mstmp99 + Mstmp99*M[7] + x*M[38] + y*M[36] + z*M[35] + M[59];
#pragma omp atomic
Ms[60] += Mstmp103*z + Mstmp11*Mstmp62 + Mstmp11*M[24] + Mstmp113*M[3] + Mstmp152*M[0] + Mstmp30*Mstmp45 + Mstmp41*Mstmp43 + Mstmp41*M[19] + Mstmp45*M[14] + Mstmp8*Mstmp99 + Mstmp99*M[8] + x*M[39] + z*M[36] + M[60];
#pragma omp atomic
Ms[61] += Mstmp104*y + Mstmp11*Mstmp65 + Mstmp11*M[25] + Mstmp117*M[4] + Mstmp154*Mstmp67 + Mstmp155*M[1] + Mstmp33*Mstmp45 + Mstmp34*Mstmp46 + Mstmp34*M[20] + Mstmp45*M[15] + Mstmp86*Mstmp9 + Mstmp86*M[9] + x*M[40] + y*M[37] + M[61];
#pragma omp atomic
Ms[62] += Mstmp104*z + Mstmp105*z + Mstmp106*z + Mstmp109*y + Mstmp11*Mstmp70 + Mstmp11*Mstmp71 + Mstmp11*Mstmp74 + Mstmp11*M[26] + Mstmp117*Mstmp4 + Mstmp117*M[5] + Mstmp155*M[2] + Mstmp34*Mstmp49 + Mstmp34*Mstmp50 + Mstmp34*Mstmp51 + Mstmp34*M[21] + Mstmp35*Mstmp45 + Mstmp36*Mstmp45 + Mstmp37*Mstmp45 + Mstmp45*M[16] + x*M[41] + y*M[38] + z*M[37] + M[62];
#pragma omp atomic
Ms[63] += Mstmp109*z + Mstmp11*Mstmp76 + Mstmp11*Mstmp77 + Mstmp11*Mstmp80 + Mstmp11*M[27] + Mstmp110*z + Mstmp111*z + Mstmp112*y + Mstmp113*Mstmp2 + Mstmp113*M[4] + Mstmp152*M[1] + Mstmp38*Mstmp45 + Mstmp39*Mstmp45 + Mstmp40*Mstmp45 + Mstmp41*Mstmp46 + Mstmp41*Mstmp47 + Mstmp41*Mstmp48 + Mstmp41*M[20] + Mstmp45*M[17] + x*M[42] + y*M[39] + z*M[38] + M[63];
#pragma omp atomic
Ms[64] += Mstmp11*Mstmp82 + Mstmp11*M[28] + Mstmp112*z + Mstmp113*M[5] + Mstmp152*M[2] + Mstmp154*Mstmp84 + Mstmp41*Mstmp49 + Mstmp41*M[21] + Mstmp42*Mstmp45 + Mstmp45*M[18] + Mstmp9*Mstmp96 + Mstmp96*M[9] + x*M[43] + z*M[39] + M[64];
#pragma omp atomic
Ms[65] += Mstmp0*Mstmp143 + Mstmp11*Mstmp85 + Mstmp11*M[29] + Mstmp114*y + Mstmp117*M[6] + Mstmp12*Mstmp86 + Mstmp143*M[3] + Mstmp156*M[1] + Mstmp34*Mstmp52 + Mstmp34*M[22] + Mstmp86*M[10] + x*M[44] + y*M[40] + M[65];
#pragma omp atomic
Ms[66] += Mstmp11*Mstmp87 + Mstmp11*Mstmp88 + Mstmp11*Mstmp89 + Mstmp11*M[30] + Mstmp114*z + Mstmp115*z + Mstmp116*z + Mstmp117*Mstmp7 + Mstmp117*M[7] + Mstmp118*y + Mstmp15*Mstmp86 + Mstmp156*M[2] + Mstmp16*Mstmp86 + Mstmp17*Mstmp86 + Mstmp34*Mstmp55 + Mstmp34*Mstmp57 + Mstmp34*Mstmp59 + Mstmp34*M[23] + Mstmp86*M[11] + x*M[45] + y*M[41] + z*M[40] + M[66];
#pragma omp atomic
Ms[67] += Mstmp0*Mstmp146 + Mstmp11*Mstmp90 + Mstmp11*Mstmp91 + Mstmp11*Mstmp92 + Mstmp11*M[31] + Mstmp113*Mstmp5 + Mstmp113*M[6] + Mstmp117*Mstmp8 + Mstmp117*M[8] + Mstmp118*z + Mstmp119*z + Mstmp120*z + Mstmp121*y + Mstmp146*M[3] + Mstmp34*Mstmp61 + Mstmp34*Mstmp62 + Mstmp34*Mstmp63 + Mstmp34*M[24] + Mstmp41*Mstmp52 + Mstmp41*Mstmp53 + Mstmp41*Mstmp54 + Mstmp41*M[22] + x*M[46] + y*M[42] + z*M[41] + M[67];
#pragma omp atomic
Ms[68] += Mstmp11*Mstmp93 + Mstmp11*Mstmp94 + Mstmp11*Mstmp95 + Mstmp11*M[32] + Mstmp113*Mstmp6 + Mstmp113*M[7] + Mstmp12*Mstmp96 + Mstmp121*z + Mstmp122*z + Mstmp123*z + Mstmp124*y + Mstmp13*Mstmp96 + Mstmp14*Mstmp96 + Mstmp157*M[1] + Mstmp41*Mstmp55 + Mstmp41*Mstmp56 + Mstmp41*Mstmp58 + Mstmp41*M[23] + Mstmp96*M[10] + x*M[47] + y*M[43] + z*M[42] + M[68];
#pragma omp atomic
Ms[69] += Mstmp0*Mstmp149 + Mstmp11*Mstmp97 + Mstmp11*M[33] + Mstmp113*M[8] + Mstmp124*z + Mstmp149*M[3] + Mstmp15*Mstmp96 + Mstmp157*M[2] + Mstmp41*Mstmp61 + Mstmp41*M[24] + Mstmp96*M[11] + x*M[48] + z*M[43] + M[69];
#pragma omp atomic
Ms[70] += Mstmp1*Mstmp143 + Mstmp125*y + Mstmp143*M[4] + Mstmp158*Mstmp159 + Mstmp18*Mstmp86 + Mstmp34*Mstmp64 + Mstmp34*M[25] + Mstmp86*M[12] + x*M[49] + y*M[44] + M[70];
#pragma omp atomic
Ms[71] += Mstmp125*z + Mstmp126*z + Mstmp127*z + Mstmp130*y + Mstmp143*Mstmp3 + Mstmp143*Mstmp4 + Mstmp143*M[5] + Mstmp23*Mstmp86 + Mstmp25*Mstmp86 + Mstmp27*Mstmp86 + Mstmp34*Mstmp69 + Mstmp34*Mstmp71 + Mstmp34*Mstmp73 + Mstmp34*M[26] + Mstmp86*M[13] + x*M[50] + y*M[45] + z*M[44] + M[71];
#pragma omp atomic
Ms[72] += Mstmp1*Mstmp146 + Mstmp130*z + Mstmp131*z + Mstmp132*z + Mstmp133*y + Mstmp146*M[4] + Mstmp160*M[0] + Mstmp29*Mstmp86 + Mstmp30*Mstmp86 + Mstmp31*Mstmp86 + Mstmp34*Mstmp75 + Mstmp34*Mstmp77 + Mstmp34*Mstmp79 + Mstmp34*M[27] + Mstmp41*Mstmp64 + Mstmp41*Mstmp65 + Mstmp41*Mstmp66 + Mstmp41*M[25] + Mstmp86*M[14] + x*M[51] + y*M[46] + z*M[45] + M[72];
#pragma omp atomic
Ms[73] += Mstmp108*Mstmp161 + Mstmp133*z + Mstmp134*z + Mstmp135*z + Mstmp137*y + Mstmp146*Mstmp3 + Mstmp146*M[5] + Mstmp18*Mstmp96 + Mstmp19*Mstmp96 + Mstmp20*Mstmp96 + Mstmp34*Mstmp81 + Mstmp34*Mstmp82 + Mstmp34*Mstmp83 + Mstmp34*M[28] + Mstmp41*Mstmp69 + Mstmp41*Mstmp70 + Mstmp41*Mstmp72 + Mstmp41*M[26] + Mstmp96*M[12] + x*M[52] + y*M[47] + z*M[46] + M[73];
#pragma omp atomic
Ms[74] += Mstmp1*Mstmp149 + Mstmp137*z + Mstmp138*z + Mstmp139*z + Mstmp140*y + Mstmp149*Mstmp2 + Mstmp149*M[4] + Mstmp23*Mstmp96 + Mstmp24*Mstmp96 + Mstmp26*Mstmp96 + Mstmp41*Mstmp75 + Mstmp41*Mstmp76 + Mstmp41*Mstmp78 + Mstmp41*M[27] + Mstmp96*M[13] + x*M[53] + y*M[48] + z*M[47] + M[74];
#pragma omp atomic
Ms[75] += Mstmp140*z + Mstmp149*Mstmp3 + Mstmp149*M[5] + Mstmp159*Mstmp162 + Mstmp29*Mstmp96 + Mstmp41*Mstmp81 + Mstmp41*M[28] + Mstmp96*M[14] + x*M[54] + z*M[48] + M[75];
#pragma omp atomic
Ms[76] += Mstmp143*M[6] + Mstmp163*M[1] + Mstmp34*M[29] + Mstmp86*M[15] + y*M[49] + M[76];
#pragma omp atomic
Ms[77] += Mstmp142*z + Mstmp143*Mstmp7 + Mstmp143*M[7] + Mstmp163*M[2] + Mstmp34*Mstmp88 + Mstmp34*M[30] + Mstmp36*Mstmp86 + Mstmp86*M[16] + y*M[50] + z*M[49] + M[77];
#pragma omp atomic
Ms[78] += Mstmp143*Mstmp8 + Mstmp143*M[8] + Mstmp144*z + Mstmp146*M[6] + Mstmp160*M[1] + Mstmp34*Mstmp91 + Mstmp34*M[31] + Mstmp39*Mstmp86 + Mstmp41*Mstmp85 + Mstmp41*M[29] + Mstmp86*M[17] + y*M[51] + z*M[50] + M[78];
#pragma omp atomic
Ms[79] += Mstmp145*z + Mstmp146*M[7] + Mstmp160*M[2] + Mstmp164*M[1] + Mstmp33*Mstmp96 + Mstmp34*Mstmp94 + Mstmp34*M[32] + Mstmp41*Mstmp87 + Mstmp41*M[30] + Mstmp42*Mstmp86 + Mstmp86*M[18] + Mstmp96*M[15] + y*M[52] + z*M[51] + M[79];
#pragma omp atomic
Ms[80] += Mstmp146*M[8] + Mstmp147*z + Mstmp149*Mstmp5 + Mstmp149*M[6] + Mstmp164*M[2] + Mstmp34*Mstmp97 + Mstmp34*M[33] + Mstmp35*Mstmp96 + Mstmp41*Mstmp90 + Mstmp41*M[31] + Mstmp96*M[16] + y*M[53] + z*M[52] + M[80];
#pragma omp atomic
Ms[81] += Mstmp148*z + Mstmp149*Mstmp6 + Mstmp149*M[7] + Mstmp165*M[1] + Mstmp38*Mstmp96 + Mstmp41*Mstmp93 + Mstmp41*M[32] + Mstmp96*M[17] + y*M[54] + z*M[53] + M[81];
#pragma omp atomic
Ms[82] += Mstmp149*M[8] + Mstmp165*M[2] + Mstmp41*M[33] + Mstmp96*M[18] + z*M[54] + M[82];

}

void M2L_6(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[83];
double Dtmp0 = pow(R, -3);
double Dtmp1 = 1.0*Dtmp0;
double Dtmp2 = pow(x, 2);
double Dtmp3 = pow(R, -2);
double Dtmp4 = 3.0*Dtmp3;
double Dtmp5 = pow(R, -5);
double Dtmp6 = Dtmp5*x;
double Dtmp7 = 3.0*Dtmp6;
double Dtmp8 = pow(y, 2);
double Dtmp9 = Dtmp5*y;
double Dtmp10 = 15.0*Dtmp3;
double Dtmp11 = -Dtmp10*Dtmp2;
double Dtmp12 = Dtmp11 + 3.0;
double Dtmp13 = Dtmp12*Dtmp5;
double Dtmp14 = -Dtmp10*Dtmp8;
double Dtmp15 = Dtmp14 + 3.0;
double Dtmp16 = pow(R, -7);
double Dtmp17 = Dtmp16*x;
double Dtmp18 = pow(x, 4);
double Dtmp19 = pow(R, -4);
double Dtmp20 = 105.0*Dtmp19;
double Dtmp21 = Dtmp2*Dtmp3;
double Dtmp22 = -105.0*Dtmp21;
double Dtmp23 = Dtmp22 + 45.0;
double Dtmp24 = Dtmp17*Dtmp23;
double Dtmp25 = Dtmp2*Dtmp8;
double Dtmp26 = Dtmp22 + 15.0;
double Dtmp27 = Dtmp16*y;
double Dtmp28 = Dtmp27*z;
double Dtmp29 = Dtmp3*Dtmp8;
double Dtmp30 = -105.0*Dtmp29;
double Dtmp31 = Dtmp30 + 45.0;
double Dtmp32 = 1.0*Dtmp17;
double Dtmp33 = pow(y, 4);
double Dtmp34 = 945.0*Dtmp19;
double Dtmp35 = Dtmp18*Dtmp34;
double Dtmp36 = 630.0*Dtmp21;
double Dtmp37 = Dtmp16*(Dtmp35 - Dtmp36 + 45.0);
double Dtmp38 = 315.0*Dtmp29;
double Dtmp39 = Dtmp25*Dtmp34;
double Dtmp40 = 315.0 - 945.0*Dtmp21;
double Dtmp41 = pow(R, -9);
double Dtmp42 = Dtmp41*y;
double Dtmp43 = Dtmp42*z;
double Dtmp44 = Dtmp43*x;
double Dtmp45 = 315.0*Dtmp21;
double Dtmp46 = Dtmp16*z;
double Dtmp47 = Dtmp33*Dtmp34;
double Dtmp48 = 630.0*Dtmp29;
double Dtmp49 = Dtmp47 - Dtmp48 + 45.0;
double Dtmp50 = 315.0 - 945.0*Dtmp29;
double Dtmp51 = 10395.0/pow(R, 6);
double Dtmp52 = Dtmp18*Dtmp19;
double Dtmp53 = 10395.0*Dtmp52;
double Dtmp54 = x*(-9450.0*Dtmp21 + Dtmp53 + 1575.0);
double Dtmp55 = Dtmp41*z;
double Dtmp56 = Dtmp19*Dtmp25;
double Dtmp57 = -5670.0*Dtmp56 - 45.0;
double Dtmp58 = -2835.0*Dtmp21;
double Dtmp59 = 10395.0*Dtmp56;
double Dtmp60 = -2835.0*Dtmp29 + Dtmp59;
double Dtmp61 = Dtmp42*x;
double Dtmp62 = Dtmp55*x;
double Dtmp63 = Dtmp19*Dtmp33;
double Dtmp64 = 10395.0*Dtmp63;
double Dtmp65 = -9450.0*Dtmp29 + Dtmp64 + 1575.0;
D[0] = -Dtmp1*x;
D[1] = -Dtmp1*y;
D[2] = -Dtmp1*z;
D[3] = Dtmp0*(Dtmp2*Dtmp4 - 1.0);
D[4] = Dtmp7*y;
D[5] = Dtmp7*z;
D[6] = Dtmp0*(Dtmp4*Dtmp8 - 1.0);
D[7] = 3.0*Dtmp9*z;
D[8] = -D[3] - D[6];
D[9] = Dtmp6*(Dtmp11 + 9.0);
D[10] = Dtmp13*y;
D[11] = Dtmp13*z;
D[12] = 1.0*Dtmp15*Dtmp6;
D[13] = -15.0*Dtmp17*y*z;
D[14] = -D[9] - D[12];
D[15] = Dtmp9*(Dtmp14 + 9.0);
D[16] = Dtmp15*Dtmp5*z;
D[17] = -D[10] - D[15];
D[18] = -D[11] - D[16];
D[19] = Dtmp5*(Dtmp18*Dtmp20 - 90.0*Dtmp21 + 9.0);
D[20] = -Dtmp24*y;
D[21] = -Dtmp24*z;
D[22] = Dtmp5*(Dtmp12 + Dtmp14 + Dtmp20*Dtmp25);
D[23] = -Dtmp26*Dtmp28;
D[24] = -D[19] - D[22];
D[25] = -Dtmp31*Dtmp32*y;
D[26] = -Dtmp32*z*(Dtmp30 + 15.0);
D[27] = -D[20] - D[25];
D[28] = -D[21] - D[26];
D[29] = Dtmp5*(Dtmp20*Dtmp33 - 90.0*Dtmp29 + 9.0);
D[30] = -Dtmp28*Dtmp31;
D[31] = -D[22] - D[29];
D[32] = -D[23] - D[30];
D[33] = -D[24] - D[31];
D[34] = -Dtmp17*(-1050.0*Dtmp21 + Dtmp35 + 225.0);
D[35] = -Dtmp37*y;
D[36] = -Dtmp37*z;
D[37] = -Dtmp17*(Dtmp23 - Dtmp38 + Dtmp39);
D[38] = Dtmp40*Dtmp44;
D[39] = -D[34] - D[37];
D[40] = -Dtmp27*(Dtmp31 + Dtmp39 - Dtmp45);
D[41] = -Dtmp46*(Dtmp26 + Dtmp30 + Dtmp39);
D[42] = -D[35] - D[40];
D[43] = -D[36] - D[41];
D[44] = -Dtmp32*Dtmp49;
D[45] = 1.0*Dtmp44*Dtmp50;
D[46] = -D[37] - D[44];
D[47] = -D[38] - D[45];
D[48] = -D[39] - D[46];
D[49] = -Dtmp27*(-1050.0*Dtmp29 + Dtmp47 + 225.0);
D[50] = -Dtmp46*Dtmp49;
D[51] = -D[40] - D[49];
D[52] = -D[41] - D[50];
D[53] = -D[42] - D[51];
D[54] = -D[43] - D[52];
D[55] = Dtmp16*(4725.0*Dtmp21 + Dtmp51*pow(x, 6) - 14175.0*Dtmp52 - 225.0);
D[56] = Dtmp42*Dtmp54;
D[57] = Dtmp54*Dtmp55;
D[58] = Dtmp16*(Dtmp18*Dtmp51*Dtmp8 - Dtmp35 + Dtmp36 + Dtmp38 + Dtmp57);
D[59] = Dtmp43*(-5670.0*Dtmp21 + Dtmp53 + 315.0);
D[60] = -D[55] - D[58];
D[61] = Dtmp61*(Dtmp58 + Dtmp60 + 945.0);
D[62] = Dtmp62*(Dtmp40 + Dtmp60);
D[63] = -D[56] - D[61];
D[64] = -D[57] - D[62];
D[65] = Dtmp16*(Dtmp2*Dtmp33*Dtmp51 + Dtmp45 - Dtmp47 + Dtmp48 + Dtmp57);
D[66] = Dtmp43*(Dtmp50 + Dtmp58 + Dtmp59);
D[67] = -D[58] - D[65];
D[68] = -D[59] - D[66];
D[69] = -D[60] - D[67];
D[70] = 1.0*Dtmp61*Dtmp65;
D[71] = 1.0*Dtmp62*(-5670.0*Dtmp29 + Dtmp64 + 315.0);
D[72] = -D[61] - D[70];
D[73] = -D[62] - D[71];
D[74] = -D[63] - D[72];
D[75] = -D[64] - D[73];
D[76] = Dtmp16*(4725.0*Dtmp29 + Dtmp51*pow(y, 6) - 14175.0*Dtmp63 - 225.0);
D[77] = Dtmp43*Dtmp65;
D[78] = -D[65] - D[76];
D[79] = -D[66] - D[77];
D[80] = -D[67] - D[78];
D[81] = -D[68] - D[79];
D[82] = -D[69] - D[80];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8] + D[9]*M[9] + D[10]*M[10] + D[11]*M[11] + D[12]*M[12] + D[13]*M[13] + D[14]*M[14] + D[15]*M[15] + D[16]*M[16] + D[17]*M[17] + D[18]*M[18] + D[19]*M[19] + D[20]*M[20] + D[21]*M[21] + D[22]*M[22] + D[23]*M[23] + D[24]*M[24] + D[25]*M[25] + D[26]*M[26] + D[27]*M[27] + D[28]*M[28] + D[29]*M[29] + D[30]*M[30] + D[31]*M[31] + D[32]*M[32] + D[33]*M[33] + D[34]*M[34] + D[35]*M[35] + D[36]*M[36] + D[37]*M[37] + D[38]*M[38] + D[39]*M[39] + D[40]*M[40] + D[41]*M[41] + D[42]*M[42] + D[43]*M[43] + D[44]*M[44] + D[45]*M[45] + D[46]*M[46] + D[47]*M[47] + D[48]*M[48] + D[49]*M[49] + D[50]*M[50] + D[51]*M[51] + D[52]*M[52] + D[53]*M[53] + D[54]*M[54] + D[55]*M[55] + D[56]*M[56] + D[57]*M[57] + D[58]*M[58] + D[59]*M[59] + D[60]*M[60] + D[61]*M[61] + D[62]*M[62] + D[63]*M[63] + D[64]*M[64] + D[65]*M[65] + D[66]*M[66] + D[67]*M[67] + D[68]*M[68] + D[69]*M[69] + D[70]*M[70] + D[71]*M[71] + D[72]*M[72] + D[73]*M[73] + D[74]*M[74] + D[75]*M[75] + D[76]*M[76] + D[77]*M[77] + D[78]*M[78] + D[79]*M[79] + D[80]*M[80] + D[81]*M[81] + D[82]*M[82];
#pragma omp atomic
L[1] += D[3]*M[0] + D[4]*M[1] + D[5]*M[2] + D[9]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[14]*M[8] + D[19]*M[9] + D[20]*M[10] + D[21]*M[11] + D[22]*M[12] + D[23]*M[13] + D[24]*M[14] + D[25]*M[15] + D[26]*M[16] + D[27]*M[17] + D[28]*M[18] + D[34]*M[19] + D[35]*M[20] + D[36]*M[21] + D[37]*M[22] + D[38]*M[23] + D[39]*M[24] + D[40]*M[25] + D[41]*M[26] + D[42]*M[27] + D[43]*M[28] + D[44]*M[29] + D[45]*M[30] + D[46]*M[31] + D[47]*M[32] + D[48]*M[33] + D[55]*M[34] + D[56]*M[35] + D[57]*M[36] + D[58]*M[37] + D[59]*M[38] + D[60]*M[39] + D[61]*M[40] + D[62]*M[41] + D[63]*M[42] + D[64]*M[43] + D[65]*M[44] + D[66]*M[45] + D[67]*M[46] + D[68]*M[47] + D[69]*M[48] + D[70]*M[49] + D[71]*M[50] + D[72]*M[51] + D[73]*M[52] + D[74]*M[53] + D[75]*M[54];
#pragma omp atomic
L[2] += D[4]*M[0] + D[6]*M[1] + D[7]*M[2] + D[10]*M[3] + D[12]*M[4] + D[13]*M[5] + D[15]*M[6] + D[16]*M[7] + D[17]*M[8] + D[20]*M[9] + D[22]*M[10] + D[23]*M[11] + D[25]*M[12] + D[26]*M[13] + D[27]*M[14] + D[29]*M[15] + D[30]*M[16] + D[31]*M[17] + D[32]*M[18] + D[35]*M[19] + D[37]*M[20] + D[38]*M[21] + D[40]*M[22] + D[41]*M[23] + D[42]*M[24] + D[44]*M[25] + D[45]*M[26] + D[46]*M[27] + D[47]*M[28] + D[49]*M[29] + D[50]*M[30] + D[51]*M[31] + D[52]*M[32] + D[53]*M[33] + D[56]*M[34] + D[58]*M[35] + D[59]*M[36] + D[61]*M[37] + D[62]*M[38] + D[63]*M[39] + D[65]*M[40] + D[66]*M[41] + D[67]*M[42] + D[68]*M[43] + D[70]*M[44] + D[71]*M[45] + D[72]*M[46] + D[73]*M[47] + D[74]*M[48] + D[76]*M[49] + D[77]*M[50] + D[78]*M[51] + D[79]*M[52] + D[80]*M[53] + D[81]*M[54];
#pragma omp atomic
L[3] += D[5]*M[0] + D[7]*M[1] + D[8]*M[2] + D[11]*M[3] + D[13]*M[4] + D[14]*M[5] + D[16]*M[6] + D[17]*M[7] + D[18]*M[8] + D[21]*M[9] + D[23]*M[10] + D[24]*M[11] + D[26]*M[12] + D[27]*M[13] + D[28]*M[14] + D[30]*M[15] + D[31]*M[16] + D[32]*M[17] + D[33]*M[18] + D[36]*M[19] + D[38]*M[20] + D[39]*M[21] + D[41]*M[22] + D[42]*M[23] + D[43]*M[24] + D[45]*M[25] + D[46]*M[26] + D[47]*M[27] + D[48]*M[28] + D[50]*M[29] + D[51]*M[30] + D[52]*M[31] + D[53]*M[32] + D[54]*M[33] + D[57]*M[34] + D[59]*M[35] + D[60]*M[36] + D[62]*M[37] + D[63]*M[38] + D[64]*M[39] + D[66]*M[40] + D[67]*M[41] + D[68]*M[42] + D[69]*M[43] + D[71]*M[44] + D[72]*M[45] + D[73]*M[46] + D[74]*M[47] + D[75]*M[48] + D[77]*M[49] + D[78]*M[50] + D[79]*M[51] + D[80]*M[52] + D[81]*M[53] + D[82]*M[54];
#pragma omp atomic
L[4] += D[9]*M[0] + D[10]*M[1] + D[11]*M[2] + D[19]*M[3] + D[20]*M[4] + D[21]*M[5] + D[22]*M[6] + D[23]*M[7] + D[24]*M[8] + D[34]*M[9] + D[35]*M[10] + D[36]*M[11] + D[37]*M[12] + D[38]*M[13] + D[39]*M[14] + D[40]*M[15] + D[41]*M[16] + D[42]*M[17] + D[43]*M[18] + D[55]*M[19] + D[56]*M[20] + D[57]*M[21] + D[58]*M[22] + D[59]*M[23] + D[60]*M[24] + D[61]*M[25] + D[62]*M[26] + D[63]*M[27] + D[64]*M[28] + D[65]*M[29] + D[66]*M[30] + D[67]*M[31] + D[68]*M[32] + D[69]*M[33];
#pragma omp atomic
L[5] += D[10]*M[0] + D[12]*M[1] + D[13]*M[2] + D[20]*M[3] + D[22]*M[4] + D[23]*M[5] + D[25]*M[6] + D[26]*M[7] + D[27]*M[8] + D[35]*M[9] + D[37]*M[10] + D[38]*M[11] + D[40]*M[12] + D[41]*M[13] + D[42]*M[14] + D[44]*M[15] + D[45]*M[16] + D[46]*M[17] + D[47]*M[18] + D[56]*M[19] + D[58]*M[20] + D[59]*M[21] + D[61]*M[22] + D[62]*M[23] + D[63]*M[24] + D[65]*M[25] + D[66]*M[26] + D[67]*M[27] + D[68]*M[28] + D[70]*M[29] + D[71]*M[30] + D[72]*M[31] + D[73]*M[32] + D[74]*M[33];
#pragma omp atomic
L[6] += D[11]*M[0] + D[13]*M[1] + D[14]*M[2] + D[21]*M[3] + D[23]*M[4] + D[24]*M[5] + D[26]*M[6] + D[27]*M[7] + D[28]*M[8] + D[36]*M[9] + D[38]*M[10] + D[39]*M[11] + D[41]*M[12] + D[42]*M[13] + D[43]*M[14] + D[45]*M[15] + D[46]*M[16] + D[47]*M[17] + D[48]*M[18] + D[57]*M[19] + D[59]*M[20] + D[60]*M[21] + D[62]*M[22] + D[63]*M[23] + D[64]*M[24] + D[66]*M[25] + D[67]*M[26] + D[68]*M[27] + D[69]*M[28] + D[71]*M[29] + D[72]*M[30] + D[73]*M[31] + D[74]*M[32] + D[75]*M[33];
#pragma omp atomic
L[7] += D[12]*M[0] + D[15]*M[1] + D[16]*M[2] + D[22]*M[3] + D[25]*M[4] + D[26]*M[5] + D[29]*M[6] + D[30]*M[7] + D[31]*M[8] + D[37]*M[9] + D[40]*M[10] + D[41]*M[11] + D[44]*M[12] + D[45]*M[13] + D[46]*M[14] + D[49]*M[15] + D[50]*M[16] + D[51]*M[17] + D[52]*M[18] + D[58]*M[19] + D[61]*M[20] + D[62]*M[21] + D[65]*M[22] + D[66]*M[23] + D[67]*M[24] + D[70]*M[25] + D[71]*M[26] + D[72]*M[27] + D[73]*M[28] + D[76]*M[29] + D[77]*M[30] + D[78]*M[31] + D[79]*M[32] + D[80]*M[33];
#pragma omp atomic
L[8] += D[13]*M[0] + D[16]*M[1] + D[17]*M[2] + D[23]*M[3] + D[26]*M[4] + D[27]*M[5] + D[30]*M[6] + D[31]*M[7] + D[32]*M[8] + D[38]*M[9] + D[41]*M[10] + D[42]*M[11] + D[45]*M[12] + D[46]*M[13] + D[47]*M[14] + D[50]*M[15] + D[51]*M[16] + D[52]*M[17] + D[53]*M[18] + D[59]*M[19] + D[62]*M[20] + D[63]*M[21] + D[66]*M[22] + D[67]*M[23] + D[68]*M[24] + D[71]*M[25] + D[72]*M[26] + D[73]*M[27] + D[74]*M[28] + D[77]*M[29] + D[78]*M[30] + D[79]*M[31] + D[80]*M[32] + D[81]*M[33];
#pragma omp atomic
L[9] += D[14]*M[0] + D[17]*M[1] + D[18]*M[2] + D[24]*M[3] + D[27]*M[4] + D[28]*M[5] + D[31]*M[6] + D[32]*M[7] + D[33]*M[8] + D[39]*M[9] + D[42]*M[10] + D[43]*M[11] + D[46]*M[12] + D[47]*M[13] + D[48]*M[14] + D[51]*M[15] + D[52]*M[16] + D[53]*M[17] + D[54]*M[18] + D[60]*M[19] + D[63]*M[20] + D[64]*M[21] + D[67]*M[22] + D[68]*M[23] + D[69]*M[24] + D[72]*M[25] + D[73]*M[26] + D[74]*M[27] + D[75]*M[28] + D[78]*M[29] + D[79]*M[30] + D[80]*M[31] + D[81]*M[32] + D[82]*M[33];
#pragma omp atomic
L[10] += D[19]*M[0] + D[20]*M[1] + D[21]*M[2] + D[34]*M[3] + D[35]*M[4] + D[36]*M[5] + D[37]*M[6] + D[38]*M[7] + D[39]*M[8] + D[55]*M[9] + D[56]*M[10] + D[57]*M[11] + D[58]*M[12] + D[59]*M[13] + D[60]*M[14] + D[61]*M[15] + D[62]*M[16] + D[63]*M[17] + D[64]*M[18];
#pragma omp atomic
L[11] += D[20]*M[0] + D[22]*M[1] + D[23]*M[2] + D[35]*M[3] + D[37]*M[4] + D[38]*M[5] + D[40]*M[6] + D[41]*M[7] + D[42]*M[8] + D[56]*M[9] + D[58]*M[10] + D[59]*M[11] + D[61]*M[12] + D[62]*M[13] + D[63]*M[14] + D[65]*M[15] + D[66]*M[16] + D[67]*M[17] + D[68]*M[18];
#pragma omp atomic
L[12] += D[21]*M[0] + D[23]*M[1] + D[24]*M[2] + D[36]*M[3] + D[38]*M[4] + D[39]*M[5] + D[41]*M[6] + D[42]*M[7] + D[43]*M[8] + D[57]*M[9] + D[59]*M[10] + D[60]*M[11] + D[62]*M[12] + D[63]*M[13] + D[64]*M[14] + D[66]*M[15] + D[67]*M[16] + D[68]*M[17] + D[69]*M[18];
#pragma omp atomic
L[13] += D[22]*M[0] + D[25]*M[1] + D[26]*M[2] + D[37]*M[3] + D[40]*M[4] + D[41]*M[5] + D[44]*M[6] + D[45]*M[7] + D[46]*M[8] + D[58]*M[9] + D[61]*M[10] + D[62]*M[11] + D[65]*M[12] + D[66]*M[13] + D[67]*M[14] + D[70]*M[15] + D[71]*M[16] + D[72]*M[17] + D[73]*M[18];
#pragma omp atomic
L[14] += D[23]*M[0] + D[26]*M[1] + D[27]*M[2] + D[38]*M[3] + D[41]*M[4] + D[42]*M[5] + D[45]*M[6] + D[46]*M[7] + D[47]*M[8] + D[59]*M[9] + D[62]*M[10] + D[63]*M[11] + D[66]*M[12] + D[67]*M[13] + D[68]*M[14] + D[71]*M[15] + D[72]*M[16] + D[73]*M[17] + D[74]*M[18];
#pragma omp atomic
L[15] += D[24]*M[0] + D[27]*M[1] + D[28]*M[2] + D[39]*M[3] + D[42]*M[4] + D[43]*M[5] + D[46]*M[6] + D[47]*M[7] + D[48]*M[8] + D[60]*M[9] + D[63]*M[10] + D[64]*M[11] + D[67]*M[12] + D[68]*M[13] + D[69]*M[14] + D[72]*M[15] + D[73]*M[16] + D[74]*M[17] + D[75]*M[18];
#pragma omp atomic
L[16] += D[25]*M[0] + D[29]*M[1] + D[30]*M[2] + D[40]*M[3] + D[44]*M[4] + D[45]*M[5] + D[49]*M[6] + D[50]*M[7] + D[51]*M[8] + D[61]*M[9] + D[65]*M[10] + D[66]*M[11] + D[70]*M[12] + D[71]*M[13] + D[72]*M[14] + D[76]*M[15] + D[77]*M[16] + D[78]*M[17] + D[79]*M[18];
#pragma omp atomic
L[17] += D[26]*M[0] + D[30]*M[1] + D[31]*M[2] + D[41]*M[3] + D[45]*M[4] + D[46]*M[5] + D[50]*M[6] + D[51]*M[7] + D[52]*M[8] + D[62]*M[9] + D[66]*M[10] + D[67]*M[11] + D[71]*M[12] + D[72]*M[13] + D[73]*M[14] + D[77]*M[15] + D[78]*M[16] + D[79]*M[17] + D[80]*M[18];
#pragma omp atomic
L[18] += D[27]*M[0] + D[31]*M[1] + D[32]*M[2] + D[42]*M[3] + D[46]*M[4] + D[47]*M[5] + D[51]*M[6] + D[52]*M[7] + D[53]*M[8] + D[63]*M[9] + D[67]*M[10] + D[68]*M[11] + D[72]*M[12] + D[73]*M[13] + D[74]*M[14] + D[78]*M[15] + D[79]*M[16] + D[80]*M[17] + D[81]*M[18];
#pragma omp atomic
L[19] += D[28]*M[0] + D[32]*M[1] + D[33]*M[2] + D[43]*M[3] + D[47]*M[4] + D[48]*M[5] + D[52]*M[6] + D[53]*M[7] + D[54]*M[8] + D[64]*M[9] + D[68]*M[10] + D[69]*M[11] + D[73]*M[12] + D[74]*M[13] + D[75]*M[14] + D[79]*M[15] + D[80]*M[16] + D[81]*M[17] + D[82]*M[18];
#pragma omp atomic
L[20] += D[34]*M[0] + D[35]*M[1] + D[36]*M[2] + D[55]*M[3] + D[56]*M[4] + D[57]*M[5] + D[58]*M[6] + D[59]*M[7] + D[60]*M[8];
#pragma omp atomic
L[21] += D[35]*M[0] + D[37]*M[1] + D[38]*M[2] + D[56]*M[3] + D[58]*M[4] + D[59]*M[5] + D[61]*M[6] + D[62]*M[7] + D[63]*M[8];
#pragma omp atomic
L[22] += D[36]*M[0] + D[38]*M[1] + D[39]*M[2] + D[57]*M[3] + D[59]*M[4] + D[60]*M[5] + D[62]*M[6] + D[63]*M[7] + D[64]*M[8];
#pragma omp atomic
L[23] += D[37]*M[0] + D[40]*M[1] + D[41]*M[2] + D[58]*M[3] + D[61]*M[4] + D[62]*M[5] + D[65]*M[6] + D[66]*M[7] + D[67]*M[8];
#pragma omp atomic
L[24] += D[38]*M[0] + D[41]*M[1] + D[42]*M[2] + D[59]*M[3] + D[62]*M[4] + D[63]*M[5] + D[66]*M[6] + D[67]*M[7] + D[68]*M[8];
#pragma omp atomic
L[25] += D[39]*M[0] + D[42]*M[1] + D[43]*M[2] + D[60]*M[3] + D[63]*M[4] + D[64]*M[5] + D[67]*M[6] + D[68]*M[7] + D[69]*M[8];
#pragma omp atomic
L[26] += D[40]*M[0] + D[44]*M[1] + D[45]*M[2] + D[61]*M[3] + D[65]*M[4] + D[66]*M[5] + D[70]*M[6] + D[71]*M[7] + D[72]*M[8];
#pragma omp atomic
L[27] += D[41]*M[0] + D[45]*M[1] + D[46]*M[2] + D[62]*M[3] + D[66]*M[4] + D[67]*M[5] + D[71]*M[6] + D[72]*M[7] + D[73]*M[8];
#pragma omp atomic
L[28] += D[42]*M[0] + D[46]*M[1] + D[47]*M[2] + D[63]*M[3] + D[67]*M[4] + D[68]*M[5] + D[72]*M[6] + D[73]*M[7] + D[74]*M[8];
#pragma omp atomic
L[29] += D[43]*M[0] + D[47]*M[1] + D[48]*M[2] + D[64]*M[3] + D[68]*M[4] + D[69]*M[5] + D[73]*M[6] + D[74]*M[7] + D[75]*M[8];
#pragma omp atomic
L[30] += D[44]*M[0] + D[49]*M[1] + D[50]*M[2] + D[65]*M[3] + D[70]*M[4] + D[71]*M[5] + D[76]*M[6] + D[77]*M[7] + D[78]*M[8];
#pragma omp atomic
L[31] += D[45]*M[0] + D[50]*M[1] + D[51]*M[2] + D[66]*M[3] + D[71]*M[4] + D[72]*M[5] + D[77]*M[6] + D[78]*M[7] + D[79]*M[8];
#pragma omp atomic
L[32] += D[46]*M[0] + D[51]*M[1] + D[52]*M[2] + D[67]*M[3] + D[72]*M[4] + D[73]*M[5] + D[78]*M[6] + D[79]*M[7] + D[80]*M[8];
#pragma omp atomic
L[33] += D[47]*M[0] + D[52]*M[1] + D[53]*M[2] + D[68]*M[3] + D[73]*M[4] + D[74]*M[5] + D[79]*M[6] + D[80]*M[7] + D[81]*M[8];
#pragma omp atomic
L[34] += D[48]*M[0] + D[53]*M[1] + D[54]*M[2] + D[69]*M[3] + D[74]*M[4] + D[75]*M[5] + D[80]*M[6] + D[81]*M[7] + D[82]*M[8];
#pragma omp atomic
L[35] += D[55]*M[0] + D[56]*M[1] + D[57]*M[2];
#pragma omp atomic
L[36] += D[56]*M[0] + D[58]*M[1] + D[59]*M[2];
#pragma omp atomic
L[37] += D[57]*M[0] + D[59]*M[1] + D[60]*M[2];
#pragma omp atomic
L[38] += D[58]*M[0] + D[61]*M[1] + D[62]*M[2];
#pragma omp atomic
L[39] += D[59]*M[0] + D[62]*M[1] + D[63]*M[2];
#pragma omp atomic
L[40] += D[60]*M[0] + D[63]*M[1] + D[64]*M[2];
#pragma omp atomic
L[41] += D[61]*M[0] + D[65]*M[1] + D[66]*M[2];
#pragma omp atomic
L[42] += D[62]*M[0] + D[66]*M[1] + D[67]*M[2];
#pragma omp atomic
L[43] += D[63]*M[0] + D[67]*M[1] + D[68]*M[2];
#pragma omp atomic
L[44] += D[64]*M[0] + D[68]*M[1] + D[69]*M[2];
#pragma omp atomic
L[45] += D[65]*M[0] + D[70]*M[1] + D[71]*M[2];
#pragma omp atomic
L[46] += D[66]*M[0] + D[71]*M[1] + D[72]*M[2];
#pragma omp atomic
L[47] += D[67]*M[0] + D[72]*M[1] + D[73]*M[2];
#pragma omp atomic
L[48] += D[68]*M[0] + D[73]*M[1] + D[74]*M[2];
#pragma omp atomic
L[49] += D[69]*M[0] + D[74]*M[1] + D[75]*M[2];
#pragma omp atomic
L[50] += D[70]*M[0] + D[76]*M[1] + D[77]*M[2];
#pragma omp atomic
L[51] += D[71]*M[0] + D[77]*M[1] + D[78]*M[2];
#pragma omp atomic
L[52] += D[72]*M[0] + D[78]*M[1] + D[79]*M[2];
#pragma omp atomic
L[53] += D[73]*M[0] + D[79]*M[1] + D[80]*M[2];
#pragma omp atomic
L[54] += D[74]*M[0] + D[80]*M[1] + D[81]*M[2];
#pragma omp atomic
L[55] += D[75]*M[0] + D[81]*M[1] + D[82]*M[2];

}

void L2L_6(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[14];
double Lstmp4 = Lstmp3*y;
double Lstmp5 = pow(x, 2);
double Lstmp6 = (1.0/2.0)*Lstmp5;
double Lstmp7 = pow(x, 3);
double Lstmp8 = (1.0/6.0)*Lstmp7;
double Lstmp9 = (1.0/24.0)*pow(x, 4);
double Lstmp10 = pow(y, 2);
double Lstmp11 = (1.0/2.0)*Lstmp10;
double Lstmp12 = pow(y, 3);
double Lstmp13 = (1.0/6.0)*Lstmp12;
double Lstmp14 = (1.0/24.0)*pow(y, 4);
double Lstmp15 = pow(z, 2);
double Lstmp16 = (1.0/2.0)*Lstmp15;
double Lstmp17 = pow(z, 3);
double Lstmp18 = (1.0/6.0)*Lstmp17;
double Lstmp19 = (1.0/24.0)*pow(z, 4);
double Lstmp20 = x*L[13];
double Lstmp21 = x*L[26];
double Lstmp22 = x*L[45];
double Lstmp23 = x*L[15];
double Lstmp24 = x*L[29];
double Lstmp25 = x*L[49];
double Lstmp26 = y*L[11];
double Lstmp27 = z*L[12];
double Lstmp28 = y*L[21];
double Lstmp29 = z*L[22];
double Lstmp30 = y*L[36];
double Lstmp31 = z*L[37];
double Lstmp32 = y*L[18];
double Lstmp33 = y*L[33];
double Lstmp34 = y*L[54];
double Lstmp35 = z*L[17];
double Lstmp36 = z*L[31];
double Lstmp37 = z*L[51];
double Lstmp38 = y*L[28];
double Lstmp39 = Lstmp38*x;
double Lstmp40 = y*L[48];
double Lstmp41 = Lstmp40*x;
double Lstmp42 = z*L[27];
double Lstmp43 = Lstmp42*x;
double Lstmp44 = z*L[46];
double Lstmp45 = Lstmp44*x;
double Lstmp46 = z*L[24];
double Lstmp47 = Lstmp46*y;
double Lstmp48 = z*L[39];
double Lstmp49 = Lstmp48*y;
double Lstmp50 = (1.0/4.0)*Lstmp5;
double Lstmp51 = Lstmp10*Lstmp50;
double Lstmp52 = (1.0/12.0)*Lstmp5;
double Lstmp53 = Lstmp15*Lstmp50;
double Lstmp54 = (1.0/12.0)*Lstmp7;
double Lstmp55 = (1.0/4.0)*Lstmp10*Lstmp15;
double Lstmp56 = x*L[47];
double Lstmp57 = y*L[43];
double Lstmp58 = z*L[42];
double Lstmp59 = x*L[23];
double Lstmp60 = x*L[41];
double Lstmp61 = x*L[25];
double Lstmp62 = x*L[44];
double Lstmp63 = Lstmp57*x;
double Lstmp64 = Lstmp58*x;
double Lstmp65 = y*L[13];
double Lstmp66 = Lstmp42*y;
double Lstmp67 = x*L[28];
double Lstmp68 = x*L[48];
double Lstmp69 = y*L[23];
double Lstmp70 = y*L[38];
double Lstmp71 = y*L[32];
double Lstmp72 = y*L[53];
double Lstmp73 = y*L[47];
double Lstmp74 = Lstmp73*x;
double Lstmp75 = Lstmp58*y;
double Lstmp76 = y*L[14];
double Lstmp77 = z*L[15];
double Lstmp78 = z*L[18];
double Lstmp79 = z*L[28];
double Lstmp80 = Lstmp79*y;
double Lstmp81 = x*L[27];
double Lstmp82 = x*L[46];
double Lstmp83 = y*L[24];
double Lstmp84 = z*L[25];
double Lstmp85 = y*L[39];
double Lstmp86 = z*L[40];
double Lstmp87 = z*L[32];
double Lstmp88 = z*L[52];
double Lstmp89 = z*L[47];
double Lstmp90 = Lstmp89*x;
double Lstmp91 = z*L[43];
double Lstmp92 = Lstmp91*y;
double Lstmp93 = x*L[38];
double Lstmp94 = x*L[40];
double Lstmp95 = x*L[43];
double Lstmp96 = x*L[42];
double Lstmp97 = y*L[26];
double Lstmp98 = Lstmp44*y;
double Lstmp99 = y*L[41];
double Lstmp100 = y*L[52];
double Lstmp101 = y*L[27];
double Lstmp102 = Lstmp89*y;
double Lstmp103 = y*L[42];
double Lstmp104 = z*L[29];
double Lstmp105 = z*L[33];
double Lstmp106 = z*L[48];
double Lstmp107 = Lstmp106*y;
double Lstmp108 = z*L[44];
double Lstmp109 = z*L[53];
double Lstmp110 = y*L[45];
double Lstmp111 = y*L[46];
double Lstmp112 = z*L[49];
double Lstmp113 = z*L[54];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + (1.0/12.0)*Lstmp10*Lstmp17*L[53] + Lstmp10*Lstmp54*L[38] + Lstmp11*Lstmp20 + Lstmp11*Lstmp35 + Lstmp11*Lstmp43 + Lstmp11*L[7] + (1.0/12.0)*Lstmp12*Lstmp15*L[52] + Lstmp12*Lstmp52*L[41] + Lstmp13*Lstmp21 + Lstmp13*Lstmp36 + Lstmp13*Lstmp45 + Lstmp13*L[16] + Lstmp14*Lstmp22 + Lstmp14*Lstmp37 + Lstmp14*L[30] + Lstmp15*Lstmp54*L[40] + Lstmp16*Lstmp23 + Lstmp16*Lstmp32 + Lstmp16*Lstmp39 + Lstmp16*L[9] + Lstmp17*Lstmp52*L[44] + Lstmp18*Lstmp24 + Lstmp18*Lstmp33 + Lstmp18*Lstmp41 + Lstmp18*L[19] + Lstmp19*Lstmp25 + Lstmp19*Lstmp34 + Lstmp19*L[34] + Lstmp2*y + Lstmp26*Lstmp6 + Lstmp27*Lstmp6 + Lstmp28*Lstmp8 + Lstmp29*Lstmp8 + Lstmp30*Lstmp9 + Lstmp31*Lstmp9 + Lstmp4*x + Lstmp47*Lstmp6 + Lstmp49*Lstmp8 + Lstmp51*Lstmp58 + Lstmp51*L[23] + Lstmp53*Lstmp57 + Lstmp53*L[25] + Lstmp55*Lstmp56 + Lstmp55*L[32] + Lstmp6*L[4] + Lstmp8*L[10] + Lstmp9*L[20] + (1.0/120.0)*pow(x, 5)*L[35] + x*L[1] + (1.0/120.0)*pow(y, 5)*L[50] + y*L[2] + (1.0/120.0)*pow(z, 5)*L[55] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + Lstmp11*Lstmp42 + Lstmp11*Lstmp59 + Lstmp11*Lstmp64 + Lstmp11*L[13] + Lstmp13*Lstmp44 + Lstmp13*Lstmp60 + Lstmp13*L[26] + Lstmp14*L[45] + Lstmp16*Lstmp38 + Lstmp16*Lstmp61 + Lstmp16*Lstmp63 + Lstmp16*L[15] + Lstmp18*Lstmp40 + Lstmp18*Lstmp62 + Lstmp18*L[29] + Lstmp19*L[49] + Lstmp26*x + Lstmp27*x + Lstmp28*Lstmp6 + Lstmp29*Lstmp6 + Lstmp30*Lstmp8 + Lstmp31*Lstmp8 + Lstmp4 + Lstmp47*x + Lstmp49*Lstmp6 + Lstmp51*L[38] + Lstmp53*L[40] + Lstmp55*L[47] + Lstmp6*L[10] + Lstmp8*L[20] + Lstmp9*L[35] + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp11*Lstmp21 + Lstmp11*Lstmp36 + Lstmp11*Lstmp45 + Lstmp11*L[16] + Lstmp13*Lstmp22 + Lstmp13*Lstmp37 + Lstmp13*L[30] + Lstmp14*L[50] + Lstmp16*Lstmp67 + Lstmp16*Lstmp71 + Lstmp16*Lstmp74 + Lstmp16*L[18] + Lstmp18*Lstmp68 + Lstmp18*Lstmp72 + Lstmp18*L[33] + Lstmp19*L[54] + Lstmp2 + Lstmp3*x + Lstmp35*y + Lstmp46*Lstmp6 + Lstmp48*Lstmp8 + Lstmp51*L[41] + Lstmp53*L[43] + Lstmp55*L[52] + Lstmp6*Lstmp69 + Lstmp6*Lstmp75 + Lstmp6*L[11] + Lstmp65*x + Lstmp66*x + Lstmp70*Lstmp8 + Lstmp8*L[21] + Lstmp9*L[36] + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += Lstmp11*Lstmp81 + Lstmp11*Lstmp87 + Lstmp11*Lstmp90 + Lstmp11*L[17] + Lstmp13*Lstmp82 + Lstmp13*Lstmp88 + Lstmp13*L[31] + Lstmp14*L[51] + Lstmp16*Lstmp24 + Lstmp16*Lstmp33 + Lstmp16*Lstmp41 + Lstmp16*L[19] + Lstmp18*Lstmp25 + Lstmp18*Lstmp34 + Lstmp18*L[34] + Lstmp19*L[55] + Lstmp51*L[42] + Lstmp53*L[44] + Lstmp55*L[53] + Lstmp6*Lstmp83 + Lstmp6*Lstmp84 + Lstmp6*Lstmp92 + Lstmp6*L[12] + Lstmp76*x + Lstmp77*x + Lstmp78*y + Lstmp8*Lstmp85 + Lstmp8*Lstmp86 + Lstmp8*L[22] + Lstmp80*x + Lstmp9*L[37] + x*L[6] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += Lstmp11*Lstmp58 + Lstmp11*Lstmp93 + Lstmp11*L[23] + Lstmp13*L[41] + Lstmp16*Lstmp57 + Lstmp16*Lstmp94 + Lstmp16*L[25] + Lstmp18*L[44] + Lstmp26 + Lstmp27 + Lstmp28*x + Lstmp29*x + Lstmp30*Lstmp6 + Lstmp31*Lstmp6 + Lstmp47 + Lstmp49*x + Lstmp6*L[20] + Lstmp8*L[35] + x*L[10] + L[4];
#pragma omp atomic
Ls[5] += Lstmp11*Lstmp44 + Lstmp11*Lstmp60 + Lstmp11*L[26] + Lstmp13*L[45] + Lstmp16*Lstmp73 + Lstmp16*Lstmp95 + Lstmp16*L[28] + Lstmp18*L[48] + Lstmp3 + Lstmp46*x + Lstmp48*Lstmp6 + Lstmp6*Lstmp70 + Lstmp6*L[21] + Lstmp65 + Lstmp66 + Lstmp69*x + Lstmp75*x + Lstmp8*L[36] + x*L[11] + L[5];
#pragma omp atomic
Ls[6] += Lstmp11*Lstmp89 + Lstmp11*Lstmp96 + Lstmp11*L[27] + Lstmp13*L[46] + Lstmp16*Lstmp40 + Lstmp16*Lstmp62 + Lstmp16*L[29] + Lstmp18*L[49] + Lstmp6*Lstmp85 + Lstmp6*Lstmp86 + Lstmp6*L[22] + Lstmp76 + Lstmp77 + Lstmp8*L[37] + Lstmp80 + Lstmp83*x + Lstmp84*x + Lstmp92*x + x*L[12] + L[6];
#pragma omp atomic
Ls[7] += Lstmp100*Lstmp16 + Lstmp11*Lstmp22 + Lstmp11*Lstmp37 + Lstmp11*L[30] + Lstmp13*L[50] + Lstmp16*Lstmp56 + Lstmp16*L[32] + Lstmp18*L[53] + Lstmp20 + Lstmp35 + Lstmp36*y + Lstmp43 + Lstmp58*Lstmp6 + Lstmp6*Lstmp99 + Lstmp6*L[23] + Lstmp8*L[38] + Lstmp97*x + Lstmp98*x + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += Lstmp101*x + Lstmp102*x + Lstmp103*Lstmp6 + Lstmp11*Lstmp82 + Lstmp11*Lstmp88 + Lstmp11*L[31] + Lstmp13*L[51] + Lstmp16*Lstmp68 + Lstmp16*Lstmp72 + Lstmp16*L[33] + Lstmp18*L[54] + Lstmp6*Lstmp91 + Lstmp6*L[24] + Lstmp78 + Lstmp79*x + Lstmp8*L[39] + Lstmp87*y + x*L[14] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += Lstmp104*x + Lstmp105*y + Lstmp107*x + Lstmp108*Lstmp6 + Lstmp109*Lstmp11 + Lstmp11*Lstmp56 + Lstmp11*L[32] + Lstmp13*L[52] + Lstmp16*Lstmp25 + Lstmp16*Lstmp34 + Lstmp16*L[34] + Lstmp18*L[55] + Lstmp23 + Lstmp32 + Lstmp39 + Lstmp57*Lstmp6 + Lstmp6*L[25] + Lstmp8*L[40] + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += Lstmp11*L[38] + Lstmp16*L[40] + Lstmp28 + Lstmp29 + Lstmp30*x + Lstmp31*x + Lstmp49 + Lstmp6*L[35] + x*L[20] + L[10];
#pragma omp atomic
Ls[11] += Lstmp11*L[41] + Lstmp16*L[43] + Lstmp46 + Lstmp48*x + Lstmp6*L[36] + Lstmp69 + Lstmp70*x + Lstmp75 + x*L[21] + L[11];
#pragma omp atomic
Ls[12] += Lstmp11*L[42] + Lstmp16*L[44] + Lstmp6*L[37] + Lstmp83 + Lstmp84 + Lstmp85*x + Lstmp86*x + Lstmp92 + x*L[22] + L[12];
#pragma omp atomic
Ls[13] += Lstmp11*L[45] + Lstmp16*L[47] + Lstmp42 + Lstmp59 + Lstmp6*L[38] + Lstmp64 + Lstmp97 + Lstmp98 + Lstmp99*x + L[13];
#pragma omp atomic
Ls[14] += Lstmp101 + Lstmp102 + Lstmp103*x + Lstmp11*L[46] + Lstmp16*L[48] + Lstmp6*L[39] + Lstmp79 + Lstmp91*x + x*L[24] + L[14];
#pragma omp atomic
Ls[15] += Lstmp104 + Lstmp107 + Lstmp108*x + Lstmp11*L[47] + Lstmp16*L[49] + Lstmp38 + Lstmp6*L[40] + Lstmp61 + Lstmp63 + L[15];
#pragma omp atomic
Ls[16] += Lstmp11*L[50] + Lstmp110*x + Lstmp16*L[52] + Lstmp21 + Lstmp36 + Lstmp37*y + Lstmp45 + Lstmp6*L[41] + y*L[30] + L[16];
#pragma omp atomic
Ls[17] += Lstmp11*L[51] + Lstmp111*x + Lstmp16*L[53] + Lstmp6*L[42] + Lstmp81 + Lstmp87 + Lstmp88*y + Lstmp90 + y*L[31] + L[17];
#pragma omp atomic
Ls[18] += Lstmp105 + Lstmp106*x + Lstmp109*y + Lstmp11*L[52] + Lstmp16*L[54] + Lstmp6*L[43] + Lstmp67 + Lstmp71 + Lstmp74 + L[18];
#pragma omp atomic
Ls[19] += Lstmp11*L[53] + Lstmp112*x + Lstmp113*y + Lstmp16*L[55] + Lstmp24 + Lstmp33 + Lstmp41 + Lstmp6*L[44] + z*L[34] + L[19];
#pragma omp atomic
Ls[20] += Lstmp30 + Lstmp31 + x*L[35] + L[20];
#pragma omp atomic
Ls[21] += Lstmp48 + Lstmp70 + x*L[36] + L[21];
#pragma omp atomic
Ls[22] += Lstmp85 + Lstmp86 + x*L[37] + L[22];
#pragma omp atomic
Ls[23] += Lstmp58 + Lstmp93 + Lstmp99 + L[23];
#pragma omp atomic
Ls[24] += Lstmp103 + Lstmp91 + x*L[39] + L[24];
#pragma omp atomic
Ls[25] += Lstmp108 + Lstmp57 + Lstmp94 + L[25];
#pragma omp atomic
Ls[26] += Lstmp110 + Lstmp44 + Lstmp60 + L[26];
#pragma omp atomic
Ls[27] += Lstmp111 + Lstmp89 + Lstmp96 + L[27];
#pragma omp atomic
Ls[28] += Lstmp106 + Lstmp73 + Lstmp95 + L[28];
#pragma omp atomic
Ls[29] += Lstmp112 + Lstmp40 + Lstmp62 + L[29];
#pragma omp atomic
Ls[30] += Lstmp22 + Lstmp37 + y*L[50] + L[30];
#pragma omp atomic
Ls[31] += Lstmp82 + Lstmp88 + y*L[51] + L[31];
#pragma omp atomic
Ls[32] += Lstmp100 + Lstmp109 + Lstmp56 + L[32];
#pragma omp atomic
Ls[33] += Lstmp113 + Lstmp68 + Lstmp72 + L[33];
#pragma omp atomic
Ls[34] += Lstmp25 + Lstmp34 + z*L[55] + L[34];
#pragma omp atomic
Ls[35] += L[35];
#pragma omp atomic
Ls[36] += L[36];
#pragma omp atomic
Ls[37] += L[37];
#pragma omp atomic
Ls[38] += L[38];
#pragma omp atomic
Ls[39] += L[39];
#pragma omp atomic
Ls[40] += L[40];
#pragma omp atomic
Ls[41] += L[41];
#pragma omp atomic
Ls[42] += L[42];
#pragma omp atomic
Ls[43] += L[43];
#pragma omp atomic
Ls[44] += L[44];
#pragma omp atomic
Ls[45] += L[45];
#pragma omp atomic
Ls[46] += L[46];
#pragma omp atomic
Ls[47] += L[47];
#pragma omp atomic
Ls[48] += L[48];
#pragma omp atomic
Ls[49] += L[49];
#pragma omp atomic
Ls[50] += L[50];
#pragma omp atomic
Ls[51] += L[51];
#pragma omp atomic
Ls[52] += L[52];
#pragma omp atomic
Ls[53] += L[53];
#pragma omp atomic
Ls[54] += L[54];
#pragma omp atomic
Ls[55] += L[55];

}

void L2P_6(double x, double y, double z, double * L, double * F) {
double Ftmp0 = x*y;
double Ftmp1 = x*z;
double Ftmp2 = y*z;
double Ftmp3 = Ftmp0*z;
double Ftmp4 = pow(x, 2);
double Ftmp5 = (1.0/2.0)*Ftmp4;
double Ftmp6 = (1.0/6.0)*pow(x, 3);
double Ftmp7 = (1.0/24.0)*pow(x, 4);
double Ftmp8 = pow(y, 2);
double Ftmp9 = (1.0/2.0)*Ftmp8;
double Ftmp10 = (1.0/6.0)*pow(y, 3);
double Ftmp11 = (1.0/24.0)*pow(y, 4);
double Ftmp12 = pow(z, 2);
double Ftmp13 = (1.0/2.0)*Ftmp12;
double Ftmp14 = (1.0/6.0)*pow(z, 3);
double Ftmp15 = (1.0/24.0)*pow(z, 4);
double Ftmp16 = Ftmp9*x;
double Ftmp17 = Ftmp10*x;
double Ftmp18 = Ftmp13*x;
double Ftmp19 = Ftmp14*x;
double Ftmp20 = Ftmp5*y;
double Ftmp21 = Ftmp5*z;
double Ftmp22 = Ftmp6*y;
double Ftmp23 = Ftmp6*z;
double Ftmp24 = Ftmp13*y;
double Ftmp25 = Ftmp14*y;
double Ftmp26 = Ftmp9*z;
double Ftmp27 = Ftmp10*z;
double Ftmp28 = Ftmp0*Ftmp13;
double Ftmp29 = Ftmp1*Ftmp9;
double Ftmp30 = Ftmp2*Ftmp5;
double Ftmp31 = (1.0/4.0)*Ftmp4;
double Ftmp32 = Ftmp31*Ftmp8;
double Ftmp33 = Ftmp12*Ftmp31;
double Ftmp34 = (1.0/4.0)*Ftmp12*Ftmp8;
#pragma omp atomic
F[0] += -Ftmp0*L[11] - Ftmp1*L[12] - Ftmp10*L[26] - Ftmp11*L[45] - Ftmp13*L[15] - Ftmp14*L[29] - Ftmp15*L[49] - Ftmp16*L[23] - Ftmp17*L[41] - Ftmp18*L[25] - Ftmp19*L[44] - Ftmp2*L[14] - Ftmp20*L[21] - Ftmp21*L[22] - Ftmp22*L[36] - Ftmp23*L[37] - Ftmp24*L[28] - Ftmp25*L[48] - Ftmp26*L[27] - Ftmp27*L[46] - Ftmp28*L[43] - Ftmp29*L[42] - Ftmp3*L[24] - Ftmp30*L[39] - Ftmp32*L[38] - Ftmp33*L[40] - Ftmp34*L[47] - Ftmp5*L[10] - Ftmp6*L[20] - Ftmp7*L[35] - Ftmp9*L[13] - x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -Ftmp0*L[13] - Ftmp1*L[14] - Ftmp10*L[30] - Ftmp11*L[50] - Ftmp13*L[18] - Ftmp14*L[33] - Ftmp15*L[54] - Ftmp16*L[26] - Ftmp17*L[45] - Ftmp18*L[28] - Ftmp19*L[48] - Ftmp2*L[17] - Ftmp20*L[23] - Ftmp21*L[24] - Ftmp22*L[38] - Ftmp23*L[39] - Ftmp24*L[32] - Ftmp25*L[53] - Ftmp26*L[31] - Ftmp27*L[51] - Ftmp28*L[47] - Ftmp29*L[46] - Ftmp3*L[27] - Ftmp30*L[42] - Ftmp32*L[41] - Ftmp33*L[43] - Ftmp34*L[52] - Ftmp5*L[11] - Ftmp6*L[21] - Ftmp7*L[36] - Ftmp9*L[16] - x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -Ftmp0*L[14] - Ftmp1*L[15] - Ftmp10*L[31] - Ftmp11*L[51] - Ftmp13*L[19] - Ftmp14*L[34] - Ftmp15*L[55] - Ftmp16*L[27] - Ftmp17*L[46] - Ftmp18*L[29] - Ftmp19*L[49] - Ftmp2*L[18] - Ftmp20*L[24] - Ftmp21*L[25] - Ftmp22*L[39] - Ftmp23*L[40] - Ftmp24*L[33] - Ftmp25*L[54] - Ftmp26*L[32] - Ftmp27*L[52] - Ftmp28*L[48] - Ftmp29*L[47] - Ftmp3*L[28] - Ftmp30*L[43] - Ftmp32*L[42] - Ftmp33*L[44] - Ftmp34*L[53] - Ftmp5*L[12] - Ftmp6*L[22] - Ftmp7*L[37] - Ftmp9*L[17] - x*L[6] - y*L[8] - z*L[9] - L[3];

}

void M2P_6(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(R, -3);
double Ftmp1 = pow(R, -2);
double Ftmp2 = 3.0*Ftmp1;
double Ftmp3 = y*M[4];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = pow(R, -4);
double Ftmp6 = Ftmp5*z;
double Ftmp7 = 15.0*Ftmp6;
double Ftmp8 = y*M[13];
double Ftmp9 = Ftmp2*x;
double Ftmp10 = Ftmp9*y;
double Ftmp11 = Ftmp4*M[2];
double Ftmp12 = pow(x, 2);
double Ftmp13 = Ftmp1*Ftmp12;
double Ftmp14 = y*M[7];
double Ftmp15 = Ftmp7*x;
double Ftmp16 = Ftmp12*Ftmp5;
double Ftmp17 = 15.0*Ftmp16;
double Ftmp18 = pow(R, -6);
double Ftmp19 = Ftmp18*y;
double Ftmp20 = Ftmp19*z;
double Ftmp21 = 105.0*M[13];
double Ftmp22 = pow(y, 2);
double Ftmp23 = 15.0*Ftmp1;
double Ftmp24 = -Ftmp22*Ftmp23;
double Ftmp25 = Ftmp1*(Ftmp24 + 3.0);
double Ftmp26 = pow(z, 2);
double Ftmp27 = -Ftmp23*Ftmp26;
double Ftmp28 = Ftmp1*(Ftmp27 + 3.0);
double Ftmp29 = -15.0*Ftmp13;
double Ftmp30 = Ftmp1*(Ftmp29 + 9.0);
double Ftmp31 = -105.0*Ftmp13;
double Ftmp32 = Ftmp31 + 45.0;
double Ftmp33 = Ftmp32*Ftmp5;
double Ftmp34 = y*M[20];
double Ftmp35 = Ftmp33*M[21];
double Ftmp36 = Ftmp5*y;
double Ftmp37 = Ftmp1*Ftmp26;
double Ftmp38 = 3.0*M[27];
double Ftmp39 = Ftmp38*(5.0 - 35.0*Ftmp37);
double Ftmp40 = 105.0*Ftmp1;
double Ftmp41 = -Ftmp22*Ftmp40;
double Ftmp42 = Ftmp41 + 45.0;
double Ftmp43 = Ftmp36*Ftmp42;
double Ftmp44 = 1.0*M[25];
double Ftmp45 = 1.0*Ftmp6;
double Ftmp46 = Ftmp41 + 15.0;
double Ftmp47 = Ftmp46*M[26];
double Ftmp48 = -Ftmp26*Ftmp40;
double Ftmp49 = Ftmp48 + 45.0;
double Ftmp50 = Ftmp45*Ftmp49;
double Ftmp51 = Ftmp25*M[6];
double Ftmp52 = Ftmp28*M[8];
double Ftmp53 = Ftmp33*x;
double Ftmp54 = Ftmp53*y;
double Ftmp55 = Ftmp43*x;
double Ftmp56 = Ftmp46*M[16];
double Ftmp57 = Ftmp6*x;
double Ftmp58 = Ftmp53*z;
double Ftmp59 = Ftmp49*M[18];
double Ftmp60 = 315.0*Ftmp1;
double Ftmp61 = -Ftmp26*Ftmp60;
double Ftmp62 = Ftmp61 + 105.0;
double Ftmp63 = 3.0*M[47];
double Ftmp64 = Ftmp62*Ftmp63;
double Ftmp65 = -945.0*Ftmp13;
double Ftmp66 = Ftmp65 + 315.0;
double Ftmp67 = Ftmp66*M[38];
double Ftmp68 = 1.0*Ftmp19;
double Ftmp69 = Ftmp1*Ftmp22;
double Ftmp70 = -945.0*Ftmp69;
double Ftmp71 = Ftmp70 + 315.0;
double Ftmp72 = Ftmp71*M[45];
double Ftmp73 = Ftmp48 + 15.0;
double Ftmp74 = 1.0*Ftmp73*M[17];
double Ftmp75 = 1.0*Ftmp16;
double Ftmp76 = Ftmp46*M[12];
double Ftmp77 = Ftmp73*M[14];
double Ftmp78 = Ftmp18*x;
double Ftmp79 = Ftmp78*y;
double Ftmp80 = Ftmp79*z;
double Ftmp81 = Ftmp66*Ftmp80;
double Ftmp82 = Ftmp71*M[30];
double Ftmp83 = -945.0*Ftmp37;
double Ftmp84 = Ftmp83 + 315.0;
double Ftmp85 = Ftmp84*M[32];
double Ftmp86 = 1.0*Ftmp80;
double Ftmp87 = Ftmp12*Ftmp19;
double Ftmp88 = Ftmp38*(Ftmp61 + 35.0);
double Ftmp89 = Ftmp44*Ftmp71;
double Ftmp90 = Ftmp65 + 525.0;
double Ftmp91 = Ftmp12*Ftmp18;
double Ftmp92 = Ftmp18*z;
double Ftmp93 = 1.0*Ftmp92;
double Ftmp94 = Ftmp12*Ftmp93;
double Ftmp95 = Ftmp70 + 105.0;
double Ftmp96 = Ftmp95*M[26];
double Ftmp97 = Ftmp84*M[28];
double Ftmp98 = -10395.0*Ftmp13;
double Ftmp99 = Ftmp98 + 4725.0;
double Ftmp100 = pow(R, -8);
double Ftmp101 = Ftmp100*Ftmp12;
double Ftmp102 = Ftmp101*y;
double Ftmp103 = z*M[38];
double Ftmp104 = -3465.0*Ftmp37;
double Ftmp105 = Ftmp63*z*(Ftmp104 + 945.0);
double Ftmp106 = -10395.0*Ftmp69;
double Ftmp107 = Ftmp106 + 2835.0;
double Ftmp108 = 1.0*Ftmp102;
double Ftmp109 = z*M[45];
double Ftmp110 = pow(y, 4);
double Ftmp111 = 945.0*Ftmp5;
double Ftmp112 = Ftmp110*Ftmp111;
double Ftmp113 = 630.0*Ftmp1;
double Ftmp114 = Ftmp5*(Ftmp112 - Ftmp113*Ftmp22 + 45.0);
double Ftmp115 = pow(z, 4);
double Ftmp116 = Ftmp111*Ftmp115;
double Ftmp117 = Ftmp5*(-Ftmp113*Ftmp26 + Ftmp116 + 45.0);
double Ftmp118 = pow(x, 4);
double Ftmp119 = Ftmp111*Ftmp118;
double Ftmp120 = Ftmp5*(Ftmp119 - 1050.0*Ftmp13 + 225.0);
double Ftmp121 = Ftmp114*M[29];
double Ftmp122 = Ftmp117*M[33];
double Ftmp123 = Ftmp115*Ftmp5;
double Ftmp124 = 3.0*M[74];
double Ftmp125 = Ftmp124*(3465.0*Ftmp123 - 1890.0*Ftmp37 + 105.0);
double Ftmp126 = Ftmp118*Ftmp5;
double Ftmp127 = 10395.0*Ftmp126;
double Ftmp128 = Ftmp127 - 9450.0*Ftmp13 + 1575.0;
double Ftmp129 = Ftmp128*M[56];
double Ftmp130 = 10395.0*Ftmp5;
double Ftmp131 = Ftmp110*Ftmp130;
double Ftmp132 = Ftmp131 - 9450.0*Ftmp69 + 1575.0;
double Ftmp133 = Ftmp132*M[70];
double Ftmp134 = 5670.0*Ftmp1;
double Ftmp135 = Ftmp131 - Ftmp134*Ftmp22 + 315.0;
double Ftmp136 = Ftmp135*M[71];
double Ftmp137 = Ftmp128*M[57];
double Ftmp138 = Ftmp115*Ftmp130;
double Ftmp139 = Ftmp138 - 9450.0*Ftmp37 + 1575.0;
double Ftmp140 = Ftmp139*Ftmp93;
double Ftmp141 = -Ftmp134*Ftmp26 + Ftmp138 + 315.0;
double Ftmp142 = 1.0*Ftmp141*M[53];
double Ftmp143 = Ftmp128*Ftmp79;
double Ftmp144 = Ftmp132*M[49];
double Ftmp145 = Ftmp135*M[50];
double Ftmp146 = Ftmp78*z;
double Ftmp147 = Ftmp128*Ftmp146;
double Ftmp148 = Ftmp139*M[54];
double Ftmp149 = 135135.0*Ftmp126;
double Ftmp150 = -103950.0*Ftmp13 + Ftmp149 + 14175.0;
double Ftmp151 = Ftmp100*x*y*z;
double Ftmp152 = Ftmp110*Ftmp5;
double Ftmp153 = 135135.0*Ftmp152;
double Ftmp154 = Ftmp153 - 103950.0*Ftmp69 + 14175.0;
double Ftmp155 = Ftmp154*M[77];
double Ftmp156 = 1.0*Ftmp91;
double Ftmp157 = Ftmp135*M[44];
double Ftmp158 = Ftmp141*M[48];
double Ftmp159 = -145530.0*Ftmp13 + Ftmp149 + 33075.0;
double Ftmp160 = Ftmp101*z;
double Ftmp161 = 135135.0*Ftmp123;
double Ftmp162 = Ftmp161 - 103950.0*Ftmp37 + 14175.0;
double Ftmp163 = 1.0*M[81];
double Ftmp164 = Ftmp162*Ftmp163;
double Ftmp165 = 45045.0*Ftmp123;
double Ftmp166 = Ftmp124*(Ftmp165 - 20790.0*Ftmp37 + 945.0);
double Ftmp167 = Ftmp154*M[70];
double Ftmp168 = 1.0*Ftmp160;
double Ftmp169 = (Ftmp153 - 62370.0*Ftmp69 + 2835.0)*M[71];
double Ftmp170 = Ftmp162*M[75];
double Ftmp171 = 135135.0*Ftmp18;
double Ftmp172 = -Ftmp171*pow(y, 6);
double Ftmp173 = 155925.0*Ftmp5;
double Ftmp174 = 42525.0*Ftmp1;
double Ftmp175 = (Ftmp110*Ftmp173 + Ftmp172 - Ftmp174*Ftmp22 + 1575.0)*M[76];
double Ftmp176 = -Ftmp171*pow(z, 6);
double Ftmp177 = (Ftmp115*Ftmp173 - Ftmp174*Ftmp26 + Ftmp176 + 1575.0)*M[82];
double Ftmp178 = -Ftmp171*pow(x, 6);
double Ftmp179 = Ftmp111*Ftmp22;
double Ftmp180 = Ftmp179*Ftmp26;
double Ftmp181 = Ftmp5*(Ftmp180 + Ftmp46 + Ftmp48);
double Ftmp182 = -Ftmp22*Ftmp60;
double Ftmp183 = Ftmp12*Ftmp179;
double Ftmp184 = Ftmp5*(Ftmp182 + Ftmp183 + Ftmp32);
double Ftmp185 = Ftmp12*Ftmp26;
double Ftmp186 = Ftmp111*Ftmp185;
double Ftmp187 = Ftmp5*(Ftmp186 + Ftmp32 + Ftmp61);
double Ftmp188 = 2835.0*Ftmp1;
double Ftmp189 = -Ftmp188*Ftmp26;
double Ftmp190 = Ftmp130*Ftmp185;
double Ftmp191 = Ftmp189 + Ftmp190;
double Ftmp192 = Ftmp191 + Ftmp66;
double Ftmp193 = Ftmp192*M[63];
double Ftmp194 = Ftmp130*Ftmp22;
double Ftmp195 = Ftmp194*Ftmp26;
double Ftmp196 = Ftmp189 + Ftmp195;
double Ftmp197 = Ftmp196 + Ftmp71;
double Ftmp198 = Ftmp197*M[72];
double Ftmp199 = -Ftmp188*Ftmp22;
double Ftmp200 = Ftmp12*Ftmp194;
double Ftmp201 = Ftmp199 + Ftmp200;
double Ftmp202 = -2835.0*Ftmp13;
double Ftmp203 = Ftmp202 + 945.0;
double Ftmp204 = Ftmp201 + Ftmp203;
double Ftmp205 = Ftmp204*M[61];
double Ftmp206 = Ftmp201 + Ftmp66;
double Ftmp207 = Ftmp206*M[62];
double Ftmp208 = Ftmp195 + Ftmp199 + Ftmp84;
double Ftmp209 = Ftmp208*M[73];
double Ftmp210 = Ftmp191 + Ftmp203;
double Ftmp211 = Ftmp210*M[64];
double Ftmp212 = Ftmp192*Ftmp79;
double Ftmp213 = Ftmp197*Ftmp79;
double Ftmp214 = Ftmp204*Ftmp79;
double Ftmp215 = Ftmp146*Ftmp206;
double Ftmp216 = Ftmp146*Ftmp208;
double Ftmp217 = Ftmp146*Ftmp210;
double Ftmp218 = -31185.0*Ftmp13;
double Ftmp219 = Ftmp218 + 8505.0;
double Ftmp220 = Ftmp16*Ftmp22;
double Ftmp221 = 135135.0*Ftmp220;
double Ftmp222 = -31185.0*Ftmp69;
double Ftmp223 = Ftmp221 + Ftmp222;
double Ftmp224 = Ftmp151*(Ftmp219 + Ftmp223);
double Ftmp225 = Ftmp185*Ftmp5;
double Ftmp226 = 135135.0*Ftmp225;
double Ftmp227 = -31185.0*Ftmp37;
double Ftmp228 = Ftmp226 + Ftmp227;
double Ftmp229 = Ftmp151*(Ftmp219 + Ftmp228);
double Ftmp230 = Ftmp26*Ftmp5;
double Ftmp231 = Ftmp22*Ftmp230;
double Ftmp232 = 135135.0*Ftmp231;
double Ftmp233 = Ftmp227 + Ftmp232;
double Ftmp234 = Ftmp151*(Ftmp222 + Ftmp233 + 8505.0);
double Ftmp235 = 4725.0*Ftmp1;
double Ftmp236 = -Ftmp22*Ftmp235;
double Ftmp237 = -Ftmp235*Ftmp26;
double Ftmp238 = -51975.0*Ftmp37;
double Ftmp239 = Ftmp226 + Ftmp238;
double Ftmp240 = -51975.0*Ftmp69;
double Ftmp241 = Ftmp221 + Ftmp240;
double Ftmp242 = Ftmp218 + 14175.0;
double Ftmp243 = -10395.0*Ftmp37;
double Ftmp244 = Ftmp243 + 2835.0;
double Ftmp245 = Ftmp222 + Ftmp232;
double Ftmp246 = 62370.0*Ftmp22;
double Ftmp247 = Ftmp230*Ftmp246;
double Ftmp248 = Ftmp110*Ftmp171;
double Ftmp249 = -Ftmp248*Ftmp26;
double Ftmp250 = Ftmp189 + Ftmp247 + Ftmp249;
double Ftmp251 = Ftmp115*Ftmp171;
double Ftmp252 = -Ftmp22*Ftmp251;
double Ftmp253 = Ftmp247 + Ftmp252;
double Ftmp254 = Ftmp16*Ftmp246;
double Ftmp255 = -Ftmp12*Ftmp248;
double Ftmp256 = Ftmp254 + Ftmp255;
double Ftmp257 = 31185.0*Ftmp5;
double Ftmp258 = 17010.0*Ftmp1;
double Ftmp259 = Ftmp110*Ftmp257 - Ftmp22*Ftmp258;
double Ftmp260 = 62370.0*Ftmp225;
double Ftmp261 = -Ftmp12*Ftmp251;
double Ftmp262 = Ftmp260 + Ftmp261;
double Ftmp263 = Ftmp115*Ftmp257 - Ftmp258*Ftmp26;
double Ftmp264 = 14175.0*Ftmp1;
double Ftmp265 = -Ftmp22*Ftmp264;
double Ftmp266 = 103950.0*Ftmp220;
double Ftmp267 = Ftmp118*Ftmp171;
double Ftmp268 = -Ftmp22*Ftmp267;
double Ftmp269 = -Ftmp26*Ftmp264;
double Ftmp270 = 103950.0*Ftmp225;
double Ftmp271 = -Ftmp26*Ftmp267;
double Ftmp272 = Ftmp22*Ftmp257;
double Ftmp273 = -Ftmp171*Ftmp185*Ftmp22;
double Ftmp274 = x*M[13];
double Ftmp275 = Ftmp22*Ftmp5;
double Ftmp276 = 15.0*x;
double Ftmp277 = Ftmp1*(Ftmp29 + 3.0);
double Ftmp278 = Ftmp1*(Ftmp24 + 9.0);
double Ftmp279 = Ftmp31 + 15.0;
double Ftmp280 = Ftmp279*M[23];
double Ftmp281 = Ftmp42*M[30];
double Ftmp282 = Ftmp5*x;
double Ftmp283 = Ftmp277*M[3];
double Ftmp284 = Ftmp279*M[11];
double Ftmp285 = Ftmp6*y;
double Ftmp286 = Ftmp285*Ftmp42;
double Ftmp287 = Ftmp279*M[10];
double Ftmp288 = 1.0*Ftmp146;
double Ftmp289 = 1.0*Ftmp36;
double Ftmp290 = Ftmp22*Ftmp78;
double Ftmp291 = Ftmp66*Ftmp78;
double Ftmp292 = Ftmp70 + 525.0;
double Ftmp293 = Ftmp22*Ftmp92;
double Ftmp294 = Ftmp65 + 105.0;
double Ftmp295 = Ftmp294*M[23];
double Ftmp296 = Ftmp98 + 2835.0;
double Ftmp297 = Ftmp100*Ftmp22;
double Ftmp298 = Ftmp297*x;
double Ftmp299 = Ftmp106 + 4725.0;
double Ftmp300 = 1.0*Ftmp298;
double Ftmp301 = Ftmp5*(-Ftmp113*Ftmp12 + Ftmp119 + 45.0);
double Ftmp302 = Ftmp5*(Ftmp112 - 1050.0*Ftmp69 + 225.0);
double Ftmp303 = Ftmp301*M[19];
double Ftmp304 = 1.0*Ftmp78;
double Ftmp305 = Ftmp127 - 5670.0*Ftmp13 + 315.0;
double Ftmp306 = Ftmp305*M[59];
double Ftmp307 = Ftmp132*M[77];
double Ftmp308 = 1.0*Ftmp79;
double Ftmp309 = Ftmp305*M[36];
double Ftmp310 = Ftmp150*M[57];
double Ftmp311 = Ftmp18*Ftmp22;
double Ftmp312 = Ftmp305*M[35];
double Ftmp313 = Ftmp150*M[56];
double Ftmp314 = Ftmp297*z;
double Ftmp315 = (-62370.0*Ftmp13 + Ftmp149 + 2835.0)*M[59];
double Ftmp316 = Ftmp153 - 145530.0*Ftmp69 + 33075.0;
double Ftmp317 = 1.0*Ftmp151;
double Ftmp318 = (Ftmp118*Ftmp173 - 42525.0*Ftmp13 + Ftmp178 + 1575.0)*M[55];
double Ftmp319 = Ftmp5*(Ftmp186 + Ftmp31 + Ftmp73);
double Ftmp320 = -315.0*Ftmp13;
double Ftmp321 = Ftmp5*(Ftmp183 + Ftmp320 + Ftmp42);
double Ftmp322 = Ftmp5*(Ftmp180 + Ftmp42 + Ftmp61);
double Ftmp323 = Ftmp200 + Ftmp202;
double Ftmp324 = Ftmp323 + Ftmp71;
double Ftmp325 = Ftmp324*M[66];
double Ftmp326 = Ftmp190 + Ftmp202;
double Ftmp327 = Ftmp326 + Ftmp84;
double Ftmp328 = Ftmp327*M[68];
double Ftmp329 = Ftmp199 + 945.0;
double Ftmp330 = Ftmp196 + Ftmp329;
double Ftmp331 = Ftmp330*M[79];
double Ftmp332 = 1.0*M[46];
double Ftmp333 = Ftmp20*Ftmp324;
double Ftmp334 = Ftmp20*Ftmp327;
double Ftmp335 = Ftmp20*Ftmp330;
double Ftmp336 = -4725.0*Ftmp13;
double Ftmp337 = -51975.0*Ftmp13;
double Ftmp338 = Ftmp337 + 14175.0;
double Ftmp339 = 1.0*Ftmp234;
double Ftmp340 = Ftmp189 + Ftmp260 + Ftmp271;
double Ftmp341 = Ftmp254 + Ftmp268;
double Ftmp342 = 31185.0*Ftmp126 - 17010.0*Ftmp13;
double Ftmp343 = -14175.0*Ftmp13;
double Ftmp344 = 103950.0*Ftmp231;
double Ftmp345 = Ftmp1*(Ftmp27 + 9.0);
double Ftmp346 = 1.0*Ftmp282;
double Ftmp347 = Ftmp26*Ftmp304;
double Ftmp348 = Ftmp83 + 525.0;
double Ftmp349 = Ftmp19*Ftmp26;
double Ftmp350 = Ftmp100*Ftmp26;
double Ftmp351 = Ftmp350*x;
double Ftmp352 = Ftmp351*y;
double Ftmp353 = 1.0*Ftmp351;
double Ftmp354 = Ftmp5*(Ftmp116 - 1050.0*Ftmp37 + 225.0);
double Ftmp355 = Ftmp139*Ftmp68;
double Ftmp356 = Ftmp18*Ftmp26;
double Ftmp357 = Ftmp350*y;
double Ftmp358 = Ftmp161 - 145530.0*Ftmp37 + 33075.0;
double Ftmp359 = Ftmp5*(Ftmp183 + Ftmp31 + Ftmp46);
double Ftmp360 = Ftmp5*(Ftmp186 + Ftmp320 + Ftmp49);
double Ftmp361 = Ftmp5*(Ftmp180 + Ftmp182 + Ftmp49);
double Ftmp362 = Ftmp243 + 4725.0;
#pragma omp atomic
F[0] += Ftmp0*(-Ftmp10*M[1] + Ftmp102*Ftmp103*Ftmp99 + Ftmp102*Ftmp105 + Ftmp102*Ftmp159*M[56] + Ftmp102*Ftmp166 + Ftmp102*(Ftmp239 + Ftmp99)*M[63] + Ftmp102*(Ftmp241 + Ftmp242)*M[61] + Ftmp107*Ftmp108*Ftmp109 + Ftmp108*Ftmp167 + Ftmp108*(Ftmp107 + Ftmp233)*M[72] - Ftmp11*x + Ftmp114*M[44] + Ftmp117*M[48] - Ftmp12*Ftmp20*Ftmp21 - Ftmp12*Ftmp90*Ftmp92*M[21] + Ftmp120*x*M[19] + Ftmp120*M[34] + Ftmp121*x + Ftmp122*x - Ftmp125*Ftmp19 - Ftmp129*Ftmp19 - 3.0*Ftmp13*M[0] - Ftmp133*Ftmp68 - Ftmp136*Ftmp93 - Ftmp137*Ftmp92 + Ftmp14*Ftmp15 - Ftmp140*M[75] - Ftmp142*Ftmp79 - Ftmp143*M[35] - Ftmp144*Ftmp79 - Ftmp145*Ftmp146 - Ftmp146*Ftmp148 - Ftmp147*M[36] + Ftmp150*Ftmp151*M[59] + Ftmp151*Ftmp155 + Ftmp151*Ftmp164 - Ftmp156*Ftmp157 - Ftmp156*Ftmp158 - Ftmp156*(Ftmp195 + Ftmp83 + Ftmp95)*M[46] + Ftmp159*Ftmp160*M[57] + Ftmp16*(Ftmp31 + 75.0)*M[9] + Ftmp160*(Ftmp239 + Ftmp242)*M[64] + Ftmp160*(Ftmp241 + Ftmp99)*M[62] + Ftmp168*Ftmp169 + Ftmp168*Ftmp170 + Ftmp168*(Ftmp244 + Ftmp245)*M[73] + Ftmp17*Ftmp3 + Ftmp17*z*M[5] - Ftmp175*Ftmp78 - Ftmp177*Ftmp78 + Ftmp181*x*M[31] + Ftmp181*M[46] + Ftmp184*x*M[22] + Ftmp184*M[37] + Ftmp187*x*M[24] + Ftmp187*M[39] - Ftmp19*Ftmp193 - Ftmp19*Ftmp205 - Ftmp198*Ftmp68 - Ftmp2*Ftmp3 - Ftmp20*Ftmp64 - Ftmp20*Ftmp67 - Ftmp207*Ftmp92 - Ftmp209*Ftmp93 - Ftmp211*Ftmp92 - Ftmp212*M[42] - Ftmp213*M[51] - Ftmp214*M[40] - Ftmp215*M[41] - Ftmp216*M[52] - Ftmp217*M[43] + Ftmp224*M[66] + Ftmp229*M[68] + Ftmp234*M[79] - Ftmp25*M[12] - Ftmp28*M[14] - Ftmp30*x*M[3] - Ftmp30*M[9] + Ftmp33*Ftmp34 - Ftmp34*Ftmp90*Ftmp91 + Ftmp35*z + Ftmp36*Ftmp39 + Ftmp36*Ftmp74*x - Ftmp4*M[5] + Ftmp43*Ftmp44 + Ftmp45*Ftmp47 + Ftmp50*M[28] - Ftmp51*x - Ftmp52*x + Ftmp54*M[10] + Ftmp55*M[15] + Ftmp56*Ftmp57 + Ftmp57*Ftmp59 + Ftmp58*M[11] - Ftmp68*Ftmp72*z + Ftmp7*Ftmp8 + Ftmp75*Ftmp76 + Ftmp75*Ftmp77 - Ftmp78*(Ftmp135 + Ftmp250)*M[78] - Ftmp78*(Ftmp141 + Ftmp199 + Ftmp253)*M[80] - Ftmp78*(Ftmp203 + Ftmp256 + Ftmp259)*M[65] - Ftmp78*(Ftmp203 + Ftmp262 + Ftmp263)*M[69] - Ftmp78*(218295.0*Ftmp126 - 99225.0*Ftmp13 + Ftmp178 + 11025.0)*M[55] - Ftmp78*(Ftmp128 + Ftmp265 + Ftmp266 + Ftmp268)*M[58] - Ftmp78*(Ftmp128 + Ftmp269 + Ftmp270 + Ftmp271)*M[60] - Ftmp78*(Ftmp192 + Ftmp201 + Ftmp26*Ftmp272 + Ftmp273)*M[67] - Ftmp80*Ftmp82 - Ftmp81*M[23] - Ftmp85*Ftmp86 - Ftmp87*Ftmp88 - Ftmp87*Ftmp89 - Ftmp91*(Ftmp127 - 13230.0*Ftmp13 + 3675.0)*M[34] - Ftmp91*(Ftmp190 + Ftmp237 + Ftmp90)*M[39] - Ftmp91*(Ftmp200 + Ftmp236 + Ftmp90)*M[37] - Ftmp94*Ftmp96 - Ftmp94*Ftmp97 + 1.0*M[0]);
#pragma omp atomic
F[1] += Ftmp0*(-Ftmp10*M[0] + Ftmp103*Ftmp296*Ftmp298 + Ftmp105*Ftmp298 + Ftmp109*Ftmp299*Ftmp300 - Ftmp11*y + Ftmp117*M[53] + Ftmp122*y - Ftmp125*Ftmp78 - Ftmp129*Ftmp78 - Ftmp132*Ftmp20*M[50] - Ftmp132*Ftmp308*M[44] - Ftmp133*Ftmp304 - Ftmp140*M[81] - Ftmp142*Ftmp311 - Ftmp143*M[34] - Ftmp146*Ftmp21*Ftmp22 - Ftmp146*Ftmp64 - Ftmp146*Ftmp67 - Ftmp148*Ftmp20 + Ftmp15*y*M[5] + Ftmp151*Ftmp310 + Ftmp154*Ftmp317*M[71] - Ftmp158*Ftmp308 + Ftmp164*Ftmp314 + Ftmp166*Ftmp298 + Ftmp170*Ftmp317 - Ftmp177*Ftmp19 - Ftmp19*Ftmp318 - Ftmp19*(Ftmp305 + Ftmp340)*M[60] - Ftmp19*(Ftmp141 + Ftmp202 + Ftmp262)*M[69] - Ftmp19*(Ftmp253 + Ftmp263 + Ftmp329)*M[80] - Ftmp19*(Ftmp329 + Ftmp341 + Ftmp342)*M[58] - Ftmp19*(Ftmp132 + Ftmp249 + Ftmp269 + Ftmp344)*M[78] - Ftmp19*(Ftmp132 + Ftmp255 + Ftmp266 + Ftmp343)*M[65] - Ftmp19*(218295.0*Ftmp152 + Ftmp172 - 99225.0*Ftmp69 + 11025.0)*M[76] - Ftmp19*(Ftmp185*Ftmp257 + Ftmp197 + Ftmp273 + Ftmp323)*M[67] - Ftmp193*Ftmp78 - Ftmp198*Ftmp304 - Ftmp20*Ftmp309 - Ftmp205*Ftmp78 - Ftmp212*M[39] - Ftmp213*Ftmp332 - Ftmp214*M[37] - Ftmp22*Ftmp291*M[20] + Ftmp22*Ftmp7*M[7] - Ftmp22*Ftmp85*Ftmp93 + Ftmp224*M[62] + Ftmp229*M[64] + Ftmp274*Ftmp7 + Ftmp275*Ftmp276*M[4] + Ftmp275*Ftmp287 + Ftmp275*Ftmp74 + Ftmp275*(Ftmp41 + 75.0)*M[15] - Ftmp277*M[10] - Ftmp278*y*M[6] - Ftmp278*M[15] - Ftmp28*M[17] + Ftmp280*Ftmp6 + Ftmp281*Ftmp6 + Ftmp282*Ftmp39 + Ftmp282*Ftmp42*Ftmp44 - Ftmp283*y + Ftmp284*Ftmp285 + Ftmp285*Ftmp59 + Ftmp286*M[16] - Ftmp288*Ftmp72 + Ftmp289*Ftmp77*x - Ftmp290*Ftmp292*Ftmp44 - Ftmp290*Ftmp88 - Ftmp292*Ftmp293*M[30] - Ftmp293*Ftmp295 + Ftmp298*Ftmp313 + Ftmp298*(Ftmp223 + Ftmp338)*M[61] + Ftmp298*(Ftmp228 + Ftmp296)*M[63] + Ftmp300*Ftmp316*M[70] + Ftmp300*(Ftmp232 + Ftmp238 + Ftmp299)*M[72] + Ftmp301*M[35] + Ftmp302*y*M[29] + Ftmp302*M[49] + Ftmp303*y - Ftmp306*Ftmp92 - Ftmp307*Ftmp92 - Ftmp311*Ftmp312 - Ftmp311*(Ftmp131 - 13230.0*Ftmp69 + 3675.0)*M[49] - Ftmp311*(Ftmp190 + Ftmp294 + Ftmp83)*M[42] - Ftmp311*(Ftmp195 + Ftmp237 + Ftmp292)*M[51] - Ftmp311*(Ftmp200 + Ftmp292 + Ftmp336)*M[40] + Ftmp314*Ftmp315 + Ftmp314*Ftmp316*M[77] + Ftmp314*(Ftmp218 + Ftmp226 + Ftmp244)*M[68] + Ftmp314*(Ftmp221 + Ftmp299 + Ftmp337)*M[66] + Ftmp314*(Ftmp238 + Ftmp245 + 14175.0)*M[79] + Ftmp319*y*M[24] + Ftmp319*M[42] + Ftmp321*y*M[22] + Ftmp321*M[40] + Ftmp322*y*M[31] + Ftmp322*M[51] - Ftmp325*Ftmp92 - Ftmp328*Ftmp92 - Ftmp331*Ftmp92 - Ftmp333*M[41] - Ftmp334*M[43] - Ftmp335*M[52] + Ftmp339*M[73] - Ftmp4*M[7] + Ftmp50*M[32] - Ftmp52*y + Ftmp53*M[20] + Ftmp54*M[9] + 1.0*Ftmp55*M[12] - 3.0*Ftmp69*M[1] - Ftmp71*Ftmp86*M[26] - Ftmp81*M[21] - Ftmp86*Ftmp97 - Ftmp9*M[4] + 1.0*M[1]);
#pragma omp atomic
F[2] += Ftmp0*(Ftmp107*Ftmp353*y*M[45] + Ftmp114*M[50] + Ftmp121*z + Ftmp124*Ftmp151*(Ftmp165 - 34650.0*Ftmp37 + 4725.0) - Ftmp136*Ftmp304 - Ftmp137*Ftmp78 - Ftmp139*Ftmp288*M[48] - Ftmp139*Ftmp304*M[75] - Ftmp14*Ftmp2 + 15.0*Ftmp14*Ftmp230 - Ftmp144*Ftmp20 - Ftmp145*Ftmp356 - Ftmp146*Ftmp34*Ftmp66 - Ftmp147*M[34] + Ftmp15*Ftmp3 + Ftmp151*Ftmp313 + Ftmp155*Ftmp357 - Ftmp157*Ftmp288 + Ftmp163*Ftmp357*Ftmp358 + Ftmp167*Ftmp317 + Ftmp169*Ftmp353 - Ftmp175*Ftmp92 - Ftmp19*Ftmp306 - Ftmp19*Ftmp307 - Ftmp19*Ftmp325 - Ftmp19*Ftmp328 - Ftmp19*Ftmp331 - Ftmp20*Ftmp312 - Ftmp207*Ftmp78 - Ftmp209*Ftmp304 - Ftmp211*Ftmp78 - Ftmp215*M[37] - Ftmp216*Ftmp332 - Ftmp217*M[39] + Ftmp224*M[61] + Ftmp229*M[63] + Ftmp230*Ftmp276*M[5] + Ftmp230*Ftmp284 + Ftmp230*Ftmp56 + Ftmp230*(Ftmp48 + 75.0)*M[18] - Ftmp25*M[16] - Ftmp26*Ftmp291*M[21] - Ftmp26*Ftmp348*Ftmp68*M[32] - 105.0*Ftmp26*Ftmp78*Ftmp8 + 15.0*Ftmp274*Ftmp36 - Ftmp277*M[11] + Ftmp280*Ftmp36 + Ftmp281*Ftmp36 - Ftmp283*z + Ftmp285*Ftmp287 + Ftmp286*M[15] + Ftmp289*Ftmp49*M[32] - Ftmp295*Ftmp349 + Ftmp296*Ftmp352*M[38] + Ftmp301*M[36] + Ftmp303*z - Ftmp308*Ftmp72 - Ftmp309*Ftmp356 + Ftmp310*Ftmp351 + Ftmp315*Ftmp357 - Ftmp318*Ftmp92 - Ftmp333*M[40] - Ftmp334*M[42] - Ftmp335*M[51] + Ftmp339*M[72] - Ftmp345*z*M[8] - Ftmp345*M[18] + Ftmp346*Ftmp47 + Ftmp346*Ftmp49*M[28] - Ftmp347*Ftmp348*M[28] - Ftmp347*Ftmp96 - Ftmp349*Ftmp82 + Ftmp35*x + Ftmp351*(Ftmp223 + Ftmp296)*M[62] + Ftmp351*(Ftmp228 + Ftmp338)*M[64] + Ftmp352*Ftmp63*(Ftmp104 + 1575.0) + Ftmp353*Ftmp358*M[75] + Ftmp353*(Ftmp232 + Ftmp240 + Ftmp362)*M[73] + Ftmp354*z*M[33] + Ftmp354*M[54] - Ftmp355*z*M[53] - Ftmp355*M[81] - Ftmp356*(Ftmp138 - 13230.0*Ftmp37 + 3675.0)*M[54] - Ftmp356*(Ftmp190 + Ftmp336 + Ftmp348)*M[43] - Ftmp356*(Ftmp195 + Ftmp236 + Ftmp348)*M[52] - Ftmp356*(Ftmp200 + Ftmp65 + Ftmp95)*M[41] + Ftmp357*(Ftmp107 + Ftmp218 + Ftmp221)*M[66] + Ftmp357*(Ftmp226 + Ftmp337 + Ftmp362)*M[68] + Ftmp357*(Ftmp233 + Ftmp240 + 14175.0)*M[79] + Ftmp359*z*M[22] + Ftmp359*M[41] + Ftmp360*z*M[24] + Ftmp360*M[43] + Ftmp361*z*M[31] + Ftmp361*M[52] - 3.0*Ftmp37*M[2] - Ftmp38*Ftmp62*Ftmp80 - Ftmp4*x*M[0] - Ftmp4*y*M[1] + Ftmp45*Ftmp76*x + Ftmp50*x*M[14] + Ftmp50*y*M[17] - Ftmp51*z + Ftmp58*M[9] - Ftmp64*Ftmp79 - Ftmp67*Ftmp79 - Ftmp80*Ftmp89 - Ftmp9*M[5] - Ftmp92*(Ftmp135 + Ftmp202 + Ftmp256)*M[65] - Ftmp92*(Ftmp199 + Ftmp305 + Ftmp341)*M[58] - Ftmp92*(Ftmp250 + Ftmp259 + 945.0)*M[78] - Ftmp92*(Ftmp340 + Ftmp342 + 945.0)*M[60] - Ftmp92*(218295.0*Ftmp123 + Ftmp176 - 99225.0*Ftmp37 + 11025.0)*M[82] - Ftmp92*(Ftmp139 + Ftmp252 + Ftmp265 + Ftmp344)*M[80] - Ftmp92*(Ftmp139 + Ftmp261 + Ftmp270 + Ftmp343)*M[69] - Ftmp92*(Ftmp12*Ftmp272 + Ftmp208 + Ftmp273 + Ftmp326)*M[67] + 1.0*M[2]);

}

void P2M_7(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = q*z;
double Mtmp3 = pow(x, 2);
double Mtmp4 = (1.0/2.0)*q;
double Mtmp5 = Mtmp0*y;
double Mtmp6 = Mtmp0*z;
double Mtmp7 = pow(y, 2);
double Mtmp8 = Mtmp1*z;
double Mtmp9 = pow(z, 2);
double Mtmp10 = pow(x, 3);
double Mtmp11 = (1.0/6.0)*q;
double Mtmp12 = (1.0/2.0)*Mtmp3;
double Mtmp13 = (1.0/2.0)*Mtmp0;
double Mtmp14 = pow(y, 3);
double Mtmp15 = (1.0/2.0)*Mtmp7;
double Mtmp16 = (1.0/2.0)*Mtmp9;
double Mtmp17 = pow(z, 3);
double Mtmp18 = pow(x, 4);
double Mtmp19 = (1.0/24.0)*q;
double Mtmp20 = (1.0/6.0)*Mtmp10;
double Mtmp21 = Mtmp7*q;
double Mtmp22 = (1.0/4.0)*Mtmp3;
double Mtmp23 = Mtmp9*q;
double Mtmp24 = (1.0/6.0)*Mtmp0;
double Mtmp25 = pow(y, 4);
double Mtmp26 = (1.0/6.0)*Mtmp14;
double Mtmp27 = (1.0/4.0)*Mtmp9;
double Mtmp28 = (1.0/6.0)*Mtmp17;
double Mtmp29 = pow(z, 4);
double Mtmp30 = pow(x, 5);
double Mtmp31 = (1.0/120.0)*q;
double Mtmp32 = (1.0/24.0)*Mtmp18;
double Mtmp33 = (1.0/12.0)*Mtmp10;
double Mtmp34 = (1.0/12.0)*Mtmp14;
double Mtmp35 = Mtmp3*q;
double Mtmp36 = Mtmp2*Mtmp7;
double Mtmp37 = Mtmp1*Mtmp9;
double Mtmp38 = (1.0/12.0)*Mtmp17;
double Mtmp39 = (1.0/24.0)*Mtmp0;
double Mtmp40 = Mtmp0*Mtmp7;
double Mtmp41 = pow(y, 5);
double Mtmp42 = (1.0/24.0)*Mtmp25;
double Mtmp43 = (1.0/24.0)*Mtmp29;
double Mtmp44 = pow(z, 5);
double Mtmp45 = pow(x, 6);
double Mtmp46 = (1.0/720.0)*q;
double Mtmp47 = (1.0/120.0)*Mtmp30;
double Mtmp48 = (1.0/48.0)*Mtmp18;
double Mtmp49 = Mtmp14*q;
double Mtmp50 = (1.0/36.0)*Mtmp10;
double Mtmp51 = Mtmp17*q;
double Mtmp52 = (1.0/48.0)*Mtmp35;
double Mtmp53 = Mtmp2*Mtmp3;
double Mtmp54 = Mtmp3*Mtmp9;
double Mtmp55 = Mtmp1*Mtmp3;
double Mtmp56 = (1.0/120.0)*Mtmp0;
double Mtmp57 = Mtmp0*Mtmp9;
double Mtmp58 = pow(y, 6);
double Mtmp59 = (1.0/120.0)*Mtmp41;
double Mtmp60 = (1.0/48.0)*Mtmp25;
double Mtmp61 = (1.0/36.0)*Mtmp17;
double Mtmp62 = (1.0/48.0)*Mtmp29;
double Mtmp63 = (1.0/120.0)*Mtmp44;
double Mtmp64 = pow(z, 6);
double Mtmp65 = (1.0/5040.0)*q;
double Mtmp66 = (1.0/720.0)*Mtmp45;
double Mtmp67 = (1.0/240.0)*Mtmp30;
double Mtmp68 = (1.0/144.0)*Mtmp18;
double Mtmp69 = (1.0/144.0)*Mtmp25;
double Mtmp70 = Mtmp10*q;
double Mtmp71 = Mtmp19*Mtmp7;
double Mtmp72 = (1.0/144.0)*Mtmp29;
double Mtmp73 = (1.0/240.0)*Mtmp35;
double Mtmp74 = (1.0/720.0)*Mtmp0;
M[0] += -Mtmp0;
M[1] += -Mtmp1;
M[2] += -Mtmp2;
M[3] += Mtmp3*Mtmp4;
M[4] += Mtmp5;
M[5] += Mtmp6;
M[6] += Mtmp4*Mtmp7;
M[7] += Mtmp8;
M[8] += Mtmp4*Mtmp9;
M[9] += -Mtmp10*Mtmp11;
M[10] += -Mtmp1*Mtmp12;
M[11] += -Mtmp12*Mtmp2;
M[12] += -Mtmp13*Mtmp7;
M[13] += -Mtmp5*z;
M[14] += -Mtmp13*Mtmp9;
M[15] += -Mtmp11*Mtmp14;
M[16] += -Mtmp15*Mtmp2;
M[17] += -Mtmp1*Mtmp16;
M[18] += -Mtmp11*Mtmp17;
M[19] += Mtmp18*Mtmp19;
M[20] += Mtmp1*Mtmp20;
M[21] += Mtmp2*Mtmp20;
M[22] += Mtmp21*Mtmp22;
M[23] += Mtmp12*Mtmp8;
M[24] += Mtmp22*Mtmp23;
M[25] += Mtmp14*Mtmp24;
M[26] += Mtmp15*Mtmp6;
M[27] += Mtmp16*Mtmp5;
M[28] += Mtmp17*Mtmp24;
M[29] += Mtmp19*Mtmp25;
M[30] += Mtmp2*Mtmp26;
M[31] += Mtmp21*Mtmp27;
M[32] += Mtmp1*Mtmp28;
M[33] += Mtmp19*Mtmp29;
M[34] += -Mtmp30*Mtmp31;
M[35] += -Mtmp1*Mtmp32;
M[36] += -Mtmp2*Mtmp32;
M[37] += -Mtmp21*Mtmp33;
M[38] += -Mtmp20*Mtmp8;
M[39] += -Mtmp23*Mtmp33;
M[40] += -Mtmp34*Mtmp35;
M[41] += -Mtmp22*Mtmp36;
M[42] += -Mtmp22*Mtmp37;
M[43] += -Mtmp35*Mtmp38;
M[44] += -Mtmp25*Mtmp39;
M[45] += -Mtmp26*Mtmp6;
M[46] += -Mtmp27*Mtmp40;
M[47] += -Mtmp28*Mtmp5;
M[48] += -Mtmp29*Mtmp39;
M[49] += -Mtmp31*Mtmp41;
M[50] += -Mtmp2*Mtmp42;
M[51] += -Mtmp23*Mtmp34;
M[52] += -Mtmp21*Mtmp38;
M[53] += -Mtmp1*Mtmp43;
M[54] += -Mtmp31*Mtmp44;
M[55] += Mtmp45*Mtmp46;
M[56] += Mtmp1*Mtmp47;
M[57] += Mtmp2*Mtmp47;
M[58] += Mtmp21*Mtmp48;
M[59] += Mtmp32*Mtmp8;
M[60] += Mtmp23*Mtmp48;
M[61] += Mtmp49*Mtmp50;
M[62] += Mtmp33*Mtmp36;
M[63] += Mtmp33*Mtmp37;
M[64] += Mtmp50*Mtmp51;
M[65] += Mtmp25*Mtmp52;
M[66] += Mtmp34*Mtmp53;
M[67] += (1.0/8.0)*Mtmp21*Mtmp54;
M[68] += Mtmp38*Mtmp55;
M[69] += Mtmp29*Mtmp52;
M[70] += Mtmp41*Mtmp56;
M[71] += Mtmp42*Mtmp6;
M[72] += Mtmp34*Mtmp57;
M[73] += Mtmp38*Mtmp40;
M[74] += Mtmp43*Mtmp5;
M[75] += Mtmp44*Mtmp56;
M[76] += Mtmp46*Mtmp58;
M[77] += Mtmp2*Mtmp59;
M[78] += Mtmp23*Mtmp60;
M[79] += Mtmp49*Mtmp61;
M[80] += Mtmp21*Mtmp62;
M[81] += Mtmp1*Mtmp63;
M[82] += Mtmp46*Mtmp64;
M[83] += -Mtmp65*pow(x, 7);
M[84] += -Mtmp1*Mtmp66;
M[85] += -Mtmp2*Mtmp66;
M[86] += -Mtmp21*Mtmp67;
M[87] += -Mtmp47*Mtmp8;
M[88] += -Mtmp23*Mtmp67;
M[89] += -Mtmp49*Mtmp68;
M[90] += -Mtmp36*Mtmp48;
M[91] += -Mtmp37*Mtmp48;
M[92] += -Mtmp51*Mtmp68;
M[93] += -Mtmp69*Mtmp70;
M[94] += -Mtmp14*Mtmp2*Mtmp50;
M[95] += -Mtmp10*Mtmp71*Mtmp9;
M[96] += -Mtmp1*Mtmp17*Mtmp50;
M[97] += -Mtmp70*Mtmp72;
M[98] += -Mtmp41*Mtmp73;
M[99] += -Mtmp53*Mtmp60;
M[100] += -Mtmp14*Mtmp19*Mtmp54;
M[101] += -Mtmp17*Mtmp3*Mtmp71;
M[102] += -Mtmp55*Mtmp62;
M[103] += -Mtmp44*Mtmp73;
M[104] += -Mtmp58*Mtmp74;
M[105] += -Mtmp59*Mtmp6;
M[106] += -Mtmp57*Mtmp60;
M[107] += -Mtmp0*Mtmp14*Mtmp61;
M[108] += -Mtmp40*Mtmp62;
M[109] += -Mtmp5*Mtmp63;
M[110] += -Mtmp64*Mtmp74;
M[111] += -Mtmp65*pow(y, 7);
M[112] += -1.0/720.0*Mtmp2*Mtmp58;
M[113] += -1.0/240.0*Mtmp23*Mtmp41;
M[114] += -Mtmp51*Mtmp69;
M[115] += -Mtmp49*Mtmp72;
M[116] += -1.0/240.0*Mtmp21*Mtmp44;
M[117] += -1.0/720.0*Mtmp1*Mtmp64;
M[118] += -Mtmp65*pow(z, 7);

}
void M2M_7(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = x*M[1];
double Mstmp2 = y*M[0];
double Mstmp3 = x*M[2];
double Mstmp4 = z*M[0];
double Mstmp5 = y*M[1];
double Mstmp6 = y*M[2];
double Mstmp7 = z*M[1];
double Mstmp8 = z*M[2];
double Mstmp9 = x*M[3];
double Mstmp10 = pow(x, 2);
double Mstmp11 = (1.0/2.0)*Mstmp10;
double Mstmp12 = x*M[4];
double Mstmp13 = y*M[3];
double Mstmp14 = Mstmp0*y;
double Mstmp15 = x*M[5];
double Mstmp16 = z*M[3];
double Mstmp17 = Mstmp0*z;
double Mstmp18 = x*M[6];
double Mstmp19 = y*M[4];
double Mstmp20 = Mstmp1*y;
double Mstmp21 = pow(y, 2);
double Mstmp22 = (1.0/2.0)*M[0];
double Mstmp23 = x*M[7];
double Mstmp24 = y*M[5];
double Mstmp25 = z*M[4];
double Mstmp26 = Mstmp3*y;
double Mstmp27 = Mstmp1*z;
double Mstmp28 = Mstmp2*z;
double Mstmp29 = x*M[8];
double Mstmp30 = z*M[5];
double Mstmp31 = Mstmp3*z;
double Mstmp32 = pow(z, 2);
double Mstmp33 = y*M[6];
double Mstmp34 = (1.0/2.0)*Mstmp21;
double Mstmp35 = y*M[7];
double Mstmp36 = z*M[6];
double Mstmp37 = Mstmp5*z;
double Mstmp38 = y*M[8];
double Mstmp39 = z*M[7];
double Mstmp40 = Mstmp6*z;
double Mstmp41 = (1.0/2.0)*Mstmp32;
double Mstmp42 = z*M[8];
double Mstmp43 = x*M[9];
double Mstmp44 = pow(x, 3);
double Mstmp45 = (1.0/6.0)*Mstmp44;
double Mstmp46 = x*M[10];
double Mstmp47 = y*M[9];
double Mstmp48 = Mstmp9*y;
double Mstmp49 = x*M[11];
double Mstmp50 = z*M[9];
double Mstmp51 = Mstmp9*z;
double Mstmp52 = x*M[12];
double Mstmp53 = y*M[10];
double Mstmp54 = Mstmp12*y;
double Mstmp55 = x*M[13];
double Mstmp56 = y*M[11];
double Mstmp57 = z*M[10];
double Mstmp58 = Mstmp15*y;
double Mstmp59 = Mstmp12*z;
double Mstmp60 = Mstmp13*z;
double Mstmp61 = x*M[14];
double Mstmp62 = z*M[11];
double Mstmp63 = Mstmp15*z;
double Mstmp64 = x*M[15];
double Mstmp65 = y*M[12];
double Mstmp66 = Mstmp18*y;
double Mstmp67 = pow(y, 3);
double Mstmp68 = (1.0/6.0)*M[0];
double Mstmp69 = x*M[16];
double Mstmp70 = y*M[13];
double Mstmp71 = z*M[12];
double Mstmp72 = Mstmp23*y;
double Mstmp73 = Mstmp18*z;
double Mstmp74 = Mstmp19*z;
double Mstmp75 = x*M[17];
double Mstmp76 = y*M[14];
double Mstmp77 = z*M[13];
double Mstmp78 = Mstmp29*y;
double Mstmp79 = Mstmp23*z;
double Mstmp80 = Mstmp24*z;
double Mstmp81 = x*M[18];
double Mstmp82 = z*M[14];
double Mstmp83 = Mstmp29*z;
double Mstmp84 = pow(z, 3);
double Mstmp85 = y*M[15];
double Mstmp86 = (1.0/6.0)*Mstmp67;
double Mstmp87 = y*M[16];
double Mstmp88 = z*M[15];
double Mstmp89 = Mstmp33*z;
double Mstmp90 = y*M[17];
double Mstmp91 = z*M[16];
double Mstmp92 = Mstmp35*z;
double Mstmp93 = y*M[18];
double Mstmp94 = z*M[17];
double Mstmp95 = Mstmp38*z;
double Mstmp96 = (1.0/6.0)*Mstmp84;
double Mstmp97 = z*M[18];
double Mstmp98 = x*M[19];
double Mstmp99 = pow(x, 4);
double Mstmp100 = (1.0/24.0)*Mstmp99;
double Mstmp101 = x*M[20];
double Mstmp102 = y*M[19];
double Mstmp103 = Mstmp43*y;
double Mstmp104 = x*M[21];
double Mstmp105 = z*M[19];
double Mstmp106 = Mstmp43*z;
double Mstmp107 = x*M[22];
double Mstmp108 = y*M[20];
double Mstmp109 = Mstmp46*y;
double Mstmp110 = (1.0/4.0)*Mstmp10;
double Mstmp111 = Mstmp21*M[0];
double Mstmp112 = x*M[23];
double Mstmp113 = y*M[21];
double Mstmp114 = z*M[20];
double Mstmp115 = Mstmp49*y;
double Mstmp116 = Mstmp46*z;
double Mstmp117 = Mstmp47*z;
double Mstmp118 = x*M[24];
double Mstmp119 = z*M[21];
double Mstmp120 = Mstmp49*z;
double Mstmp121 = Mstmp110*Mstmp32;
double Mstmp122 = x*M[25];
double Mstmp123 = y*M[22];
double Mstmp124 = Mstmp52*y;
double Mstmp125 = Mstmp110*Mstmp21;
double Mstmp126 = x*M[26];
double Mstmp127 = y*M[23];
double Mstmp128 = z*M[22];
double Mstmp129 = Mstmp55*y;
double Mstmp130 = Mstmp52*z;
double Mstmp131 = Mstmp53*z;
double Mstmp132 = x*M[27];
double Mstmp133 = y*M[24];
double Mstmp134 = z*M[23];
double Mstmp135 = Mstmp61*y;
double Mstmp136 = Mstmp55*z;
double Mstmp137 = Mstmp56*z;
double Mstmp138 = x*M[28];
double Mstmp139 = z*M[24];
double Mstmp140 = Mstmp61*z;
double Mstmp141 = x*M[29];
double Mstmp142 = y*M[25];
double Mstmp143 = Mstmp64*y;
double Mstmp144 = pow(y, 4);
double Mstmp145 = (1.0/24.0)*M[0];
double Mstmp146 = x*M[30];
double Mstmp147 = y*M[26];
double Mstmp148 = z*M[25];
double Mstmp149 = Mstmp69*y;
double Mstmp150 = Mstmp64*z;
double Mstmp151 = Mstmp65*z;
double Mstmp152 = x*M[31];
double Mstmp153 = y*M[27];
double Mstmp154 = z*M[26];
double Mstmp155 = Mstmp75*y;
double Mstmp156 = Mstmp69*z;
double Mstmp157 = Mstmp70*z;
double Mstmp158 = (1.0/4.0)*Mstmp32;
double Mstmp159 = x*M[32];
double Mstmp160 = y*M[28];
double Mstmp161 = z*M[27];
double Mstmp162 = Mstmp81*y;
double Mstmp163 = Mstmp75*z;
double Mstmp164 = Mstmp76*z;
double Mstmp165 = x*M[33];
double Mstmp166 = z*M[28];
double Mstmp167 = Mstmp81*z;
double Mstmp168 = pow(z, 4);
double Mstmp169 = y*M[29];
double Mstmp170 = (1.0/24.0)*Mstmp144;
double Mstmp171 = y*M[30];
double Mstmp172 = z*M[29];
double Mstmp173 = Mstmp85*z;
double Mstmp174 = y*M[31];
double Mstmp175 = z*M[30];
double Mstmp176 = Mstmp87*z;
double Mstmp177 = Mstmp158*Mstmp21;
double Mstmp178 = y*M[32];
double Mstmp179 = z*M[31];
double Mstmp180 = Mstmp90*z;
double Mstmp181 = y*M[33];
double Mstmp182 = z*M[32];
double Mstmp183 = Mstmp93*z;
double Mstmp184 = (1.0/24.0)*Mstmp168;
double Mstmp185 = z*M[33];
double Mstmp186 = x*M[34];
double Mstmp187 = (1.0/120.0)*pow(x, 5);
double Mstmp188 = x*M[35];
double Mstmp189 = y*M[34];
double Mstmp190 = Mstmp98*y;
double Mstmp191 = x*M[36];
double Mstmp192 = x*M[37];
double Mstmp193 = y*M[35];
double Mstmp194 = Mstmp101*y;
double Mstmp195 = (1.0/12.0)*Mstmp44;
double Mstmp196 = x*M[38];
double Mstmp197 = y*M[36];
double Mstmp198 = Mstmp104*y;
double Mstmp199 = x*M[39];
double Mstmp200 = Mstmp195*Mstmp32;
double Mstmp201 = x*M[40];
double Mstmp202 = y*M[37];
double Mstmp203 = Mstmp107*y;
double Mstmp204 = (1.0/12.0)*Mstmp10;
double Mstmp205 = Mstmp67*M[0];
double Mstmp206 = Mstmp195*Mstmp21;
double Mstmp207 = x*M[41];
double Mstmp208 = y*M[38];
double Mstmp209 = Mstmp112*y;
double Mstmp210 = x*M[42];
double Mstmp211 = y*M[39];
double Mstmp212 = Mstmp118*y;
double Mstmp213 = x*M[43];
double Mstmp214 = Mstmp204*Mstmp84;
double Mstmp215 = x*M[44];
double Mstmp216 = y*M[40];
double Mstmp217 = Mstmp122*y;
double Mstmp218 = Mstmp204*Mstmp67;
double Mstmp219 = x*M[45];
double Mstmp220 = y*M[41];
double Mstmp221 = Mstmp126*y;
double Mstmp222 = x*M[46];
double Mstmp223 = y*M[42];
double Mstmp224 = Mstmp132*y;
double Mstmp225 = x*M[47];
double Mstmp226 = y*M[43];
double Mstmp227 = Mstmp138*y;
double Mstmp228 = x*M[48];
double Mstmp229 = x*M[49];
double Mstmp230 = y*M[44];
double Mstmp231 = Mstmp141*y;
double Mstmp232 = pow(y, 5);
double Mstmp233 = (1.0/120.0)*M[0];
double Mstmp234 = x*M[50];
double Mstmp235 = y*M[45];
double Mstmp236 = Mstmp146*y;
double Mstmp237 = x*M[51];
double Mstmp238 = y*M[46];
double Mstmp239 = Mstmp152*y;
double Mstmp240 = (1.0/12.0)*Mstmp32;
double Mstmp241 = x*M[52];
double Mstmp242 = y*M[47];
double Mstmp243 = Mstmp159*y;
double Mstmp244 = (1.0/12.0)*Mstmp84;
double Mstmp245 = x*M[53];
double Mstmp246 = y*M[48];
double Mstmp247 = Mstmp165*y;
double Mstmp248 = x*M[54];
double Mstmp249 = pow(z, 5);
double Mstmp250 = y*M[49];
double Mstmp251 = (1.0/120.0)*Mstmp232;
double Mstmp252 = y*M[50];
double Mstmp253 = y*M[51];
double Mstmp254 = Mstmp240*Mstmp67;
double Mstmp255 = y*M[52];
double Mstmp256 = Mstmp21*Mstmp244;
double Mstmp257 = y*M[53];
double Mstmp258 = y*M[54];
double Mstmp259 = (1.0/120.0)*Mstmp249;
double Mstmp260 = (1.0/720.0)*pow(x, 6);
double Mstmp261 = (1.0/48.0)*Mstmp99;
double Mstmp262 = Mstmp261*Mstmp32;
double Mstmp263 = (1.0/36.0)*Mstmp44;
double Mstmp264 = Mstmp21*Mstmp261;
double Mstmp265 = Mstmp263*Mstmp84;
double Mstmp266 = (1.0/48.0)*Mstmp10;
double Mstmp267 = Mstmp266*M[0];
double Mstmp268 = Mstmp263*Mstmp67;
double Mstmp269 = (1.0/8.0)*Mstmp10*Mstmp32;
double Mstmp270 = Mstmp144*Mstmp266;
double Mstmp271 = Mstmp21*Mstmp269;
double Mstmp272 = Mstmp168*Mstmp266;
double Mstmp273 = pow(y, 6);
double Mstmp274 = (1.0/720.0)*M[0];
double Mstmp275 = (1.0/48.0)*Mstmp144*Mstmp32;
double Mstmp276 = (1.0/36.0)*Mstmp84;
double Mstmp277 = (1.0/48.0)*Mstmp168;
double Mstmp278 = pow(z, 6);
double Mstmp279 = (1.0/720.0)*Mstmp273;
double Mstmp280 = Mstmp276*Mstmp67;
double Mstmp281 = Mstmp21*Mstmp277;
double Mstmp282 = (1.0/720.0)*Mstmp278;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += Mstmp0 + M[3];
#pragma omp atomic
Ms[4] += Mstmp1 + Mstmp2 + M[4];
#pragma omp atomic
Ms[5] += Mstmp3 + Mstmp4 + M[5];
#pragma omp atomic
Ms[6] += Mstmp5 + M[6];
#pragma omp atomic
Ms[7] += Mstmp6 + Mstmp7 + M[7];
#pragma omp atomic
Ms[8] += Mstmp8 + M[8];
#pragma omp atomic
Ms[9] += Mstmp11*M[0] + Mstmp9 + M[9];
#pragma omp atomic
Ms[10] += Mstmp11*M[1] + Mstmp12 + Mstmp13 + Mstmp14 + M[10];
#pragma omp atomic
Ms[11] += Mstmp11*M[2] + Mstmp15 + Mstmp16 + Mstmp17 + M[11];
#pragma omp atomic
Ms[12] += Mstmp18 + Mstmp19 + Mstmp20 + Mstmp21*Mstmp22 + M[12];
#pragma omp atomic
Ms[13] += Mstmp23 + Mstmp24 + Mstmp25 + Mstmp26 + Mstmp27 + Mstmp28 + M[13];
#pragma omp atomic
Ms[14] += Mstmp22*Mstmp32 + Mstmp29 + Mstmp30 + Mstmp31 + M[14];
#pragma omp atomic
Ms[15] += Mstmp33 + Mstmp34*M[1] + M[15];
#pragma omp atomic
Ms[16] += Mstmp34*M[2] + Mstmp35 + Mstmp36 + Mstmp37 + M[16];
#pragma omp atomic
Ms[17] += Mstmp38 + Mstmp39 + Mstmp40 + Mstmp41*M[1] + M[17];
#pragma omp atomic
Ms[18] += Mstmp41*M[2] + Mstmp42 + M[18];
#pragma omp atomic
Ms[19] += Mstmp11*M[3] + Mstmp43 + Mstmp45*M[0] + M[19];
#pragma omp atomic
Ms[20] += Mstmp11*Mstmp2 + Mstmp11*M[4] + Mstmp45*M[1] + Mstmp46 + Mstmp47 + Mstmp48 + M[20];
#pragma omp atomic
Ms[21] += Mstmp11*Mstmp4 + Mstmp11*M[5] + Mstmp45*M[2] + Mstmp49 + Mstmp50 + Mstmp51 + M[21];
#pragma omp atomic
Ms[22] += Mstmp0*Mstmp34 + Mstmp11*Mstmp5 + Mstmp11*M[6] + Mstmp34*M[3] + Mstmp52 + Mstmp53 + Mstmp54 + M[22];
#pragma omp atomic
Ms[23] += Mstmp11*Mstmp6 + Mstmp11*Mstmp7 + Mstmp11*M[7] + Mstmp14*z + Mstmp55 + Mstmp56 + Mstmp57 + Mstmp58 + Mstmp59 + Mstmp60 + M[23];
#pragma omp atomic
Ms[24] += Mstmp0*Mstmp41 + Mstmp11*Mstmp8 + Mstmp11*M[8] + Mstmp41*M[3] + Mstmp61 + Mstmp62 + Mstmp63 + M[24];
#pragma omp atomic
Ms[25] += Mstmp1*Mstmp34 + Mstmp34*M[4] + Mstmp64 + Mstmp65 + Mstmp66 + Mstmp67*Mstmp68 + M[25];
#pragma omp atomic
Ms[26] += Mstmp20*z + Mstmp3*Mstmp34 + Mstmp34*Mstmp4 + Mstmp34*M[5] + Mstmp69 + Mstmp70 + Mstmp71 + Mstmp72 + Mstmp73 + Mstmp74 + M[26];
#pragma omp atomic
Ms[27] += Mstmp1*Mstmp41 + Mstmp2*Mstmp41 + Mstmp26*z + Mstmp41*M[4] + Mstmp75 + Mstmp76 + Mstmp77 + Mstmp78 + Mstmp79 + Mstmp80 + M[27];
#pragma omp atomic
Ms[28] += Mstmp3*Mstmp41 + Mstmp41*M[5] + Mstmp68*Mstmp84 + Mstmp81 + Mstmp82 + Mstmp83 + M[28];
#pragma omp atomic
Ms[29] += Mstmp34*M[6] + Mstmp85 + Mstmp86*M[1] + M[29];
#pragma omp atomic
Ms[30] += Mstmp34*Mstmp7 + Mstmp34*M[7] + Mstmp86*M[2] + Mstmp87 + Mstmp88 + Mstmp89 + M[30];
#pragma omp atomic
Ms[31] += Mstmp34*Mstmp8 + Mstmp34*M[8] + Mstmp41*Mstmp5 + Mstmp41*M[6] + Mstmp90 + Mstmp91 + Mstmp92 + M[31];
#pragma omp atomic
Ms[32] += Mstmp41*Mstmp6 + Mstmp41*M[7] + Mstmp93 + Mstmp94 + Mstmp95 + Mstmp96*M[1] + M[32];
#pragma omp atomic
Ms[33] += Mstmp41*M[8] + Mstmp96*M[2] + Mstmp97 + M[33];
#pragma omp atomic
Ms[34] += Mstmp100*M[0] + Mstmp11*M[9] + Mstmp45*M[3] + Mstmp98 + M[34];
#pragma omp atomic
Ms[35] += Mstmp100*M[1] + Mstmp101 + Mstmp102 + Mstmp103 + Mstmp11*Mstmp13 + Mstmp11*M[10] + Mstmp2*Mstmp45 + Mstmp45*M[4] + M[35];
#pragma omp atomic
Ms[36] += Mstmp100*M[2] + Mstmp104 + Mstmp105 + Mstmp106 + Mstmp11*Mstmp16 + Mstmp11*M[11] + Mstmp4*Mstmp45 + Mstmp45*M[5] + M[36];
#pragma omp atomic
Ms[37] += Mstmp107 + Mstmp108 + Mstmp109 + Mstmp11*Mstmp19 + Mstmp11*M[12] + Mstmp110*Mstmp111 + Mstmp34*Mstmp9 + Mstmp34*M[9] + Mstmp45*Mstmp5 + Mstmp45*M[6] + M[37];
#pragma omp atomic
Ms[38] += Mstmp11*Mstmp24 + Mstmp11*Mstmp25 + Mstmp11*Mstmp28 + Mstmp11*M[13] + Mstmp112 + Mstmp113 + Mstmp114 + Mstmp115 + Mstmp116 + Mstmp117 + Mstmp45*Mstmp6 + Mstmp45*Mstmp7 + Mstmp45*M[7] + Mstmp48*z + M[38];
#pragma omp atomic
Ms[39] += Mstmp11*Mstmp30 + Mstmp11*M[14] + Mstmp118 + Mstmp119 + Mstmp120 + Mstmp121*M[0] + Mstmp41*Mstmp9 + Mstmp41*M[9] + Mstmp45*Mstmp8 + Mstmp45*M[8] + M[39];
#pragma omp atomic
Ms[40] += Mstmp0*Mstmp86 + Mstmp11*Mstmp33 + Mstmp11*M[15] + Mstmp12*Mstmp34 + Mstmp122 + Mstmp123 + Mstmp124 + Mstmp125*M[1] + Mstmp34*M[10] + Mstmp86*M[3] + M[40];
#pragma omp atomic
Ms[41] += Mstmp11*Mstmp35 + Mstmp11*Mstmp36 + Mstmp11*Mstmp37 + Mstmp11*M[16] + Mstmp125*M[2] + Mstmp126 + Mstmp127 + Mstmp128 + Mstmp129 + Mstmp130 + Mstmp131 + Mstmp15*Mstmp34 + Mstmp16*Mstmp34 + Mstmp17*Mstmp34 + Mstmp34*M[11] + Mstmp54*z + M[41];
#pragma omp atomic
Ms[42] += Mstmp11*Mstmp38 + Mstmp11*Mstmp39 + Mstmp11*Mstmp40 + Mstmp11*M[17] + Mstmp12*Mstmp41 + Mstmp121*M[1] + Mstmp13*Mstmp41 + Mstmp132 + Mstmp133 + Mstmp134 + Mstmp135 + Mstmp136 + Mstmp137 + Mstmp14*Mstmp41 + Mstmp41*M[10] + Mstmp58*z + M[42];
#pragma omp atomic
Ms[43] += Mstmp0*Mstmp96 + Mstmp11*Mstmp42 + Mstmp11*M[18] + Mstmp121*M[2] + Mstmp138 + Mstmp139 + Mstmp140 + Mstmp15*Mstmp41 + Mstmp41*M[11] + Mstmp96*M[3] + M[43];
#pragma omp atomic
Ms[44] += Mstmp1*Mstmp86 + Mstmp141 + Mstmp142 + Mstmp143 + Mstmp144*Mstmp145 + Mstmp18*Mstmp34 + Mstmp34*M[12] + Mstmp86*M[4] + M[44];
#pragma omp atomic
Ms[45] += Mstmp146 + Mstmp147 + Mstmp148 + Mstmp149 + Mstmp150 + Mstmp151 + Mstmp23*Mstmp34 + Mstmp25*Mstmp34 + Mstmp27*Mstmp34 + Mstmp3*Mstmp86 + Mstmp34*M[13] + Mstmp4*Mstmp86 + Mstmp66*z + Mstmp86*M[5] + M[45];
#pragma omp atomic
Ms[46] += Mstmp111*Mstmp158 + Mstmp152 + Mstmp153 + Mstmp154 + Mstmp155 + Mstmp156 + Mstmp157 + Mstmp18*Mstmp41 + Mstmp19*Mstmp41 + Mstmp20*Mstmp41 + Mstmp29*Mstmp34 + Mstmp30*Mstmp34 + Mstmp31*Mstmp34 + Mstmp34*M[14] + Mstmp41*M[12] + Mstmp72*z + M[46];
#pragma omp atomic
Ms[47] += Mstmp1*Mstmp96 + Mstmp159 + Mstmp160 + Mstmp161 + Mstmp162 + Mstmp163 + Mstmp164 + Mstmp2*Mstmp96 + Mstmp23*Mstmp41 + Mstmp24*Mstmp41 + Mstmp26*Mstmp41 + Mstmp41*M[13] + Mstmp78*z + Mstmp96*M[4] + M[47];
#pragma omp atomic
Ms[48] += Mstmp145*Mstmp168 + Mstmp165 + Mstmp166 + Mstmp167 + Mstmp29*Mstmp41 + Mstmp3*Mstmp96 + Mstmp41*M[14] + Mstmp96*M[5] + M[48];
#pragma omp atomic
Ms[49] += Mstmp169 + Mstmp170*M[1] + Mstmp34*M[15] + Mstmp86*M[6] + M[49];
#pragma omp atomic
Ms[50] += Mstmp170*M[2] + Mstmp171 + Mstmp172 + Mstmp173 + Mstmp34*Mstmp36 + Mstmp34*M[16] + Mstmp7*Mstmp86 + Mstmp86*M[7] + M[50];
#pragma omp atomic
Ms[51] += Mstmp174 + Mstmp175 + Mstmp176 + Mstmp177*M[1] + Mstmp33*Mstmp41 + Mstmp34*Mstmp39 + Mstmp34*M[17] + Mstmp41*M[15] + Mstmp8*Mstmp86 + Mstmp86*M[8] + M[51];
#pragma omp atomic
Ms[52] += Mstmp177*M[2] + Mstmp178 + Mstmp179 + Mstmp180 + Mstmp34*Mstmp42 + Mstmp34*M[18] + Mstmp35*Mstmp41 + Mstmp41*M[16] + Mstmp5*Mstmp96 + Mstmp96*M[6] + M[52];
#pragma omp atomic
Ms[53] += Mstmp181 + Mstmp182 + Mstmp183 + Mstmp184*M[1] + Mstmp38*Mstmp41 + Mstmp41*M[17] + Mstmp6*Mstmp96 + Mstmp96*M[7] + M[53];
#pragma omp atomic
Ms[54] += Mstmp184*M[2] + Mstmp185 + Mstmp41*M[18] + Mstmp96*M[8] + M[54];
#pragma omp atomic
Ms[55] += Mstmp100*M[3] + Mstmp11*M[19] + Mstmp186 + Mstmp187*M[0] + Mstmp45*M[9] + M[55];
#pragma omp atomic
Ms[56] += Mstmp100*Mstmp2 + Mstmp100*M[4] + Mstmp11*Mstmp47 + Mstmp11*M[20] + Mstmp13*Mstmp45 + Mstmp187*M[1] + Mstmp188 + Mstmp189 + Mstmp190 + Mstmp45*M[10] + M[56];
#pragma omp atomic
Ms[57] += Mstmp100*Mstmp4 + Mstmp100*M[5] + Mstmp11*Mstmp50 + Mstmp11*M[21] + Mstmp16*Mstmp45 + Mstmp187*M[2] + Mstmp191 + Mstmp45*M[11] + Mstmp98*z + z*M[34] + M[57];
#pragma omp atomic
Ms[58] += Mstmp100*Mstmp5 + Mstmp100*M[6] + Mstmp11*Mstmp53 + Mstmp11*M[22] + Mstmp111*Mstmp195 + Mstmp125*M[3] + Mstmp19*Mstmp45 + Mstmp192 + Mstmp193 + Mstmp194 + Mstmp34*Mstmp43 + Mstmp34*M[19] + Mstmp45*M[12] + M[58];
#pragma omp atomic
Ms[59] += Mstmp100*Mstmp6 + Mstmp100*Mstmp7 + Mstmp100*M[7] + Mstmp101*z + Mstmp102*z + Mstmp103*z + Mstmp11*Mstmp56 + Mstmp11*Mstmp57 + Mstmp11*Mstmp60 + Mstmp11*M[23] + Mstmp196 + Mstmp197 + Mstmp198 + Mstmp24*Mstmp45 + Mstmp25*Mstmp45 + Mstmp28*Mstmp45 + Mstmp45*M[13] + z*M[35] + M[59];
#pragma omp atomic
Ms[60] += Mstmp100*Mstmp8 + Mstmp100*M[8] + Mstmp104*z + Mstmp11*Mstmp62 + Mstmp11*M[24] + Mstmp121*M[3] + Mstmp199 + Mstmp200*M[0] + Mstmp30*Mstmp45 + Mstmp41*Mstmp43 + Mstmp41*M[19] + Mstmp45*M[14] + z*M[36] + M[60];
#pragma omp atomic
Ms[61] += Mstmp11*Mstmp65 + Mstmp11*M[25] + Mstmp125*M[4] + Mstmp201 + Mstmp202 + Mstmp203 + Mstmp204*Mstmp205 + Mstmp206*M[1] + Mstmp33*Mstmp45 + Mstmp34*Mstmp46 + Mstmp34*M[20] + Mstmp45*M[15] + Mstmp86*Mstmp9 + Mstmp86*M[9] + M[61];
#pragma omp atomic
Ms[62] += Mstmp107*z + Mstmp108*z + Mstmp109*z + Mstmp11*Mstmp70 + Mstmp11*Mstmp71 + Mstmp11*Mstmp74 + Mstmp11*M[26] + Mstmp125*Mstmp4 + Mstmp125*M[5] + Mstmp206*M[2] + Mstmp207 + Mstmp208 + Mstmp209 + Mstmp34*Mstmp49 + Mstmp34*Mstmp50 + Mstmp34*Mstmp51 + Mstmp34*M[21] + Mstmp35*Mstmp45 + Mstmp36*Mstmp45 + Mstmp37*Mstmp45 + Mstmp45*M[16] + z*M[37] + M[62];
#pragma omp atomic
Ms[63] += Mstmp11*Mstmp76 + Mstmp11*Mstmp77 + Mstmp11*Mstmp80 + Mstmp11*M[27] + Mstmp112*z + Mstmp113*z + Mstmp115*z + Mstmp121*Mstmp2 + Mstmp121*M[4] + Mstmp200*M[1] + Mstmp210 + Mstmp211 + Mstmp212 + Mstmp38*Mstmp45 + Mstmp39*Mstmp45 + Mstmp40*Mstmp45 + Mstmp41*Mstmp46 + Mstmp41*Mstmp47 + Mstmp41*Mstmp48 + Mstmp41*M[20] + Mstmp45*M[17] + z*M[38] + M[63];
#pragma omp atomic
Ms[64] += Mstmp11*Mstmp82 + Mstmp11*M[28] + Mstmp118*z + Mstmp121*M[5] + Mstmp200*M[2] + Mstmp213 + Mstmp214*M[0] + Mstmp41*Mstmp49 + Mstmp41*M[21] + Mstmp42*Mstmp45 + Mstmp45*M[18] + Mstmp9*Mstmp96 + Mstmp96*M[9] + z*M[39] + M[64];
#pragma omp atomic
Ms[65] += Mstmp0*Mstmp170 + Mstmp11*Mstmp85 + Mstmp11*M[29] + Mstmp12*Mstmp86 + Mstmp125*M[6] + Mstmp170*M[3] + Mstmp215 + Mstmp216 + Mstmp217 + Mstmp218*M[1] + Mstmp34*Mstmp52 + Mstmp34*M[22] + Mstmp86*M[10] + M[65];
#pragma omp atomic
Ms[66] += Mstmp11*Mstmp87 + Mstmp11*Mstmp88 + Mstmp11*Mstmp89 + Mstmp11*M[30] + Mstmp122*z + Mstmp123*z + Mstmp124*z + Mstmp125*Mstmp7 + Mstmp125*M[7] + Mstmp15*Mstmp86 + Mstmp16*Mstmp86 + Mstmp17*Mstmp86 + Mstmp218*M[2] + Mstmp219 + Mstmp220 + Mstmp221 + Mstmp34*Mstmp55 + Mstmp34*Mstmp57 + Mstmp34*Mstmp59 + Mstmp34*M[23] + Mstmp86*M[11] + z*M[40] + M[66];
#pragma omp atomic
Ms[67] += Mstmp0*Mstmp177 + Mstmp11*Mstmp90 + Mstmp11*Mstmp91 + Mstmp11*Mstmp92 + Mstmp11*M[31] + Mstmp121*Mstmp5 + Mstmp121*M[6] + Mstmp125*Mstmp8 + Mstmp125*M[8] + Mstmp126*z + Mstmp127*z + Mstmp129*z + Mstmp177*M[3] + Mstmp222 + Mstmp223 + Mstmp224 + Mstmp34*Mstmp61 + Mstmp34*Mstmp62 + Mstmp34*Mstmp63 + Mstmp34*M[24] + Mstmp41*Mstmp52 + Mstmp41*Mstmp53 + Mstmp41*Mstmp54 + Mstmp41*M[22] + z*M[41] + M[67];
#pragma omp atomic
Ms[68] += Mstmp11*Mstmp93 + Mstmp11*Mstmp94 + Mstmp11*Mstmp95 + Mstmp11*M[32] + Mstmp12*Mstmp96 + Mstmp121*Mstmp6 + Mstmp121*M[7] + Mstmp13*Mstmp96 + Mstmp132*z + Mstmp133*z + Mstmp135*z + Mstmp14*Mstmp96 + Mstmp214*M[1] + Mstmp225 + Mstmp226 + Mstmp227 + Mstmp41*Mstmp55 + Mstmp41*Mstmp56 + Mstmp41*Mstmp58 + Mstmp41*M[23] + Mstmp96*M[10] + z*M[42] + M[68];
#pragma omp atomic
Ms[69] += Mstmp0*Mstmp184 + Mstmp11*Mstmp97 + Mstmp11*M[33] + Mstmp121*M[8] + Mstmp138*z + Mstmp15*Mstmp96 + Mstmp184*M[3] + Mstmp214*M[2] + Mstmp228 + Mstmp41*Mstmp61 + Mstmp41*M[24] + Mstmp96*M[11] + z*M[43] + M[69];
#pragma omp atomic
Ms[70] += Mstmp1*Mstmp170 + Mstmp170*M[4] + Mstmp18*Mstmp86 + Mstmp229 + Mstmp230 + Mstmp231 + Mstmp232*Mstmp233 + Mstmp34*Mstmp64 + Mstmp34*M[25] + Mstmp86*M[12] + M[70];
#pragma omp atomic
Ms[71] += Mstmp141*z + Mstmp142*z + Mstmp143*z + Mstmp170*Mstmp3 + Mstmp170*Mstmp4 + Mstmp170*M[5] + Mstmp23*Mstmp86 + Mstmp234 + Mstmp235 + Mstmp236 + Mstmp25*Mstmp86 + Mstmp27*Mstmp86 + Mstmp34*Mstmp69 + Mstmp34*Mstmp71 + Mstmp34*Mstmp73 + Mstmp34*M[26] + Mstmp86*M[13] + z*M[44] + M[71];
#pragma omp atomic
Ms[72] += Mstmp1*Mstmp177 + Mstmp146*z + Mstmp147*z + Mstmp149*z + Mstmp177*M[4] + Mstmp205*Mstmp240 + Mstmp237 + Mstmp238 + Mstmp239 + Mstmp29*Mstmp86 + Mstmp30*Mstmp86 + Mstmp31*Mstmp86 + Mstmp34*Mstmp75 + Mstmp34*Mstmp77 + Mstmp34*Mstmp79 + Mstmp34*M[27] + Mstmp41*Mstmp64 + Mstmp41*Mstmp65 + Mstmp41*Mstmp66 + Mstmp41*M[25] + Mstmp86*M[14] + z*M[45] + M[72];
#pragma omp atomic
Ms[73] += Mstmp111*Mstmp244 + Mstmp152*z + Mstmp153*z + Mstmp155*z + Mstmp177*Mstmp3 + Mstmp177*M[5] + Mstmp18*Mstmp96 + Mstmp19*Mstmp96 + Mstmp20*Mstmp96 + Mstmp241 + Mstmp242 + Mstmp243 + Mstmp34*Mstmp81 + Mstmp34*Mstmp82 + Mstmp34*Mstmp83 + Mstmp34*M[28] + Mstmp41*Mstmp69 + Mstmp41*Mstmp70 + Mstmp41*Mstmp72 + Mstmp41*M[26] + Mstmp96*M[12] + z*M[46] + M[73];
#pragma omp atomic
Ms[74] += Mstmp1*Mstmp184 + Mstmp159*z + Mstmp160*z + Mstmp162*z + Mstmp184*Mstmp2 + Mstmp184*M[4] + Mstmp23*Mstmp96 + Mstmp24*Mstmp96 + Mstmp245 + Mstmp246 + Mstmp247 + Mstmp26*Mstmp96 + Mstmp41*Mstmp75 + Mstmp41*Mstmp76 + Mstmp41*Mstmp78 + Mstmp41*M[27] + Mstmp96*M[13] + z*M[47] + M[74];
#pragma omp atomic
Ms[75] += Mstmp165*z + Mstmp184*Mstmp3 + Mstmp184*M[5] + Mstmp233*Mstmp249 + Mstmp248 + Mstmp29*Mstmp96 + Mstmp41*Mstmp81 + Mstmp41*M[28] + Mstmp96*M[14] + z*M[48] + M[75];
#pragma omp atomic
Ms[76] += Mstmp170*M[6] + Mstmp250 + Mstmp251*M[1] + Mstmp34*M[29] + Mstmp86*M[15] + M[76];
#pragma omp atomic
Ms[77] += Mstmp169*z + Mstmp170*Mstmp7 + Mstmp170*M[7] + Mstmp251*M[2] + Mstmp252 + Mstmp34*Mstmp88 + Mstmp34*M[30] + Mstmp36*Mstmp86 + Mstmp86*M[16] + z*M[49] + M[77];
#pragma omp atomic
Ms[78] += Mstmp170*Mstmp8 + Mstmp170*M[8] + Mstmp171*z + Mstmp177*M[6] + Mstmp253 + Mstmp254*M[1] + Mstmp34*Mstmp91 + Mstmp34*M[31] + Mstmp39*Mstmp86 + Mstmp41*Mstmp85 + Mstmp41*M[29] + Mstmp86*M[17] + z*M[50] + M[78];
#pragma omp atomic
Ms[79] += Mstmp174*z + Mstmp177*M[7] + Mstmp254*M[2] + Mstmp255 + Mstmp256*M[1] + Mstmp33*Mstmp96 + Mstmp34*Mstmp94 + Mstmp34*M[32] + Mstmp41*Mstmp87 + Mstmp41*M[30] + Mstmp42*Mstmp86 + Mstmp86*M[18] + Mstmp96*M[15] + z*M[51] + M[79];
#pragma omp atomic
Ms[80] += Mstmp177*M[8] + Mstmp178*z + Mstmp184*Mstmp5 + Mstmp184*M[6] + Mstmp256*M[2] + Mstmp257 + Mstmp34*Mstmp97 + Mstmp34*M[33] + Mstmp35*Mstmp96 + Mstmp41*Mstmp90 + Mstmp41*M[31] + Mstmp96*M[16] + z*M[52] + M[80];
#pragma omp atomic
Ms[81] += Mstmp181*z + Mstmp184*Mstmp6 + Mstmp184*M[7] + Mstmp258 + Mstmp259*M[1] + Mstmp38*Mstmp96 + Mstmp41*Mstmp93 + Mstmp41*M[32] + Mstmp96*M[17] + z*M[53] + M[81];
#pragma omp atomic
Ms[82] += Mstmp184*M[8] + Mstmp259*M[2] + Mstmp41*M[33] + Mstmp96*M[18] + z*M[54] + M[82];
#pragma omp atomic
Ms[83] += Mstmp100*M[9] + Mstmp11*M[34] + Mstmp187*M[3] + Mstmp260*M[0] + Mstmp45*M[19] + x*M[55] + M[83];
#pragma omp atomic
Ms[84] += Mstmp100*Mstmp13 + Mstmp100*M[10] + Mstmp102*Mstmp11 + Mstmp11*M[35] + Mstmp186*y + Mstmp187*Mstmp2 + Mstmp187*M[4] + Mstmp260*M[1] + Mstmp45*Mstmp47 + Mstmp45*M[20] + x*M[56] + y*M[55] + M[84];
#pragma omp atomic
Ms[85] += Mstmp100*Mstmp16 + Mstmp100*M[11] + Mstmp105*Mstmp11 + Mstmp11*M[36] + Mstmp186*z + Mstmp187*Mstmp4 + Mstmp187*M[5] + Mstmp260*M[2] + Mstmp45*Mstmp50 + Mstmp45*M[21] + x*M[57] + z*M[55] + M[85];
#pragma omp atomic
Ms[86] += Mstmp100*Mstmp19 + Mstmp100*M[12] + Mstmp108*Mstmp11 + Mstmp11*M[37] + Mstmp111*Mstmp261 + Mstmp125*M[9] + Mstmp187*Mstmp5 + Mstmp187*M[6] + Mstmp188*y + Mstmp206*M[3] + Mstmp34*Mstmp98 + Mstmp34*M[34] + Mstmp45*Mstmp53 + Mstmp45*M[22] + x*M[58] + y*M[56] + M[86];
#pragma omp atomic
Ms[87] += Mstmp100*Mstmp24 + Mstmp100*Mstmp25 + Mstmp100*Mstmp28 + Mstmp100*M[13] + Mstmp11*Mstmp113 + Mstmp11*Mstmp114 + Mstmp11*Mstmp117 + Mstmp11*M[38] + Mstmp187*Mstmp6 + Mstmp187*Mstmp7 + Mstmp187*M[7] + Mstmp188*z + Mstmp189*z + Mstmp190*z + Mstmp191*y + Mstmp45*Mstmp56 + Mstmp45*Mstmp57 + Mstmp45*Mstmp60 + Mstmp45*M[23] + x*M[59] + y*M[57] + z*M[56] + M[87];
#pragma omp atomic
Ms[88] += Mstmp100*Mstmp30 + Mstmp100*M[14] + Mstmp11*Mstmp119 + Mstmp11*M[39] + Mstmp121*M[9] + Mstmp187*Mstmp8 + Mstmp187*M[8] + Mstmp191*z + Mstmp200*M[3] + Mstmp262*M[0] + Mstmp41*Mstmp98 + Mstmp41*M[34] + Mstmp45*Mstmp62 + Mstmp45*M[24] + x*M[60] + z*M[57] + M[88];
#pragma omp atomic
Ms[89] += Mstmp100*Mstmp33 + Mstmp100*M[15] + Mstmp101*Mstmp34 + Mstmp11*Mstmp123 + Mstmp11*M[40] + Mstmp125*M[10] + Mstmp192*y + Mstmp205*Mstmp263 + Mstmp206*M[4] + Mstmp218*M[3] + Mstmp264*M[1] + Mstmp34*M[35] + Mstmp43*Mstmp86 + Mstmp45*Mstmp65 + Mstmp45*M[25] + Mstmp86*M[19] + x*M[61] + y*M[58] + M[89];
#pragma omp atomic
Ms[90] += Mstmp100*Mstmp35 + Mstmp100*Mstmp36 + Mstmp100*Mstmp37 + Mstmp100*M[16] + Mstmp104*Mstmp34 + Mstmp105*Mstmp34 + Mstmp106*Mstmp34 + Mstmp11*Mstmp127 + Mstmp11*Mstmp128 + Mstmp11*Mstmp131 + Mstmp11*M[41] + Mstmp125*Mstmp16 + Mstmp125*M[11] + Mstmp192*z + Mstmp193*z + Mstmp194*z + Mstmp196*y + Mstmp206*Mstmp4 + Mstmp206*M[5] + Mstmp264*M[2] + Mstmp34*M[36] + Mstmp45*Mstmp70 + Mstmp45*Mstmp71 + Mstmp45*Mstmp74 + Mstmp45*M[26] + x*M[62] + y*M[59] + z*M[58] + M[90];
#pragma omp atomic
Ms[91] += Mstmp100*Mstmp38 + Mstmp100*Mstmp39 + Mstmp100*Mstmp40 + Mstmp100*M[17] + Mstmp101*Mstmp41 + Mstmp102*Mstmp41 + Mstmp103*Mstmp41 + Mstmp11*Mstmp133 + Mstmp11*Mstmp134 + Mstmp11*Mstmp137 + Mstmp11*M[42] + Mstmp121*Mstmp13 + Mstmp121*M[10] + Mstmp196*z + Mstmp197*z + Mstmp198*z + Mstmp199*y + Mstmp2*Mstmp200 + Mstmp200*M[4] + Mstmp262*M[1] + Mstmp41*M[35] + Mstmp45*Mstmp76 + Mstmp45*Mstmp77 + Mstmp45*Mstmp80 + Mstmp45*M[27] + x*M[63] + y*M[60] + z*M[59] + M[91];
#pragma omp atomic
Ms[92] += Mstmp100*Mstmp42 + Mstmp100*M[18] + Mstmp104*Mstmp41 + Mstmp11*Mstmp139 + Mstmp11*M[43] + Mstmp121*M[11] + Mstmp199*z + Mstmp200*M[5] + Mstmp214*M[3] + Mstmp262*M[2] + Mstmp265*M[0] + Mstmp41*M[36] + Mstmp43*Mstmp96 + Mstmp45*Mstmp82 + Mstmp45*M[28] + Mstmp96*M[19] + x*M[64] + z*M[60] + M[92];
#pragma omp atomic
Ms[93] += Mstmp107*Mstmp34 + Mstmp11*Mstmp142 + Mstmp11*M[44] + Mstmp125*M[12] + Mstmp144*Mstmp267 + Mstmp170*Mstmp9 + Mstmp170*M[9] + Mstmp201*y + Mstmp206*M[6] + Mstmp218*M[4] + Mstmp268*M[1] + Mstmp34*M[37] + Mstmp45*Mstmp85 + Mstmp45*M[29] + Mstmp46*Mstmp86 + Mstmp86*M[20] + x*M[65] + y*M[61] + M[93];
#pragma omp atomic
Ms[94] += Mstmp11*Mstmp147 + Mstmp11*Mstmp148 + Mstmp11*Mstmp151 + Mstmp11*M[45] + Mstmp112*Mstmp34 + Mstmp114*Mstmp34 + Mstmp116*Mstmp34 + Mstmp125*Mstmp25 + Mstmp125*M[13] + Mstmp201*z + Mstmp202*z + Mstmp203*z + Mstmp206*Mstmp7 + Mstmp206*M[7] + Mstmp207*y + Mstmp218*Mstmp4 + Mstmp218*M[5] + Mstmp268*M[2] + Mstmp34*M[38] + Mstmp45*Mstmp87 + Mstmp45*Mstmp88 + Mstmp45*Mstmp89 + Mstmp45*M[30] + Mstmp49*Mstmp86 + Mstmp50*Mstmp86 + Mstmp51*Mstmp86 + Mstmp86*M[21] + x*M[66] + y*M[62] + z*M[61] + M[94];
#pragma omp atomic
Ms[95] += Mstmp107*Mstmp41 + Mstmp108*Mstmp41 + Mstmp109*Mstmp41 + Mstmp11*Mstmp153 + Mstmp11*Mstmp154 + Mstmp11*Mstmp157 + Mstmp11*M[46] + Mstmp111*Mstmp269 + Mstmp118*Mstmp34 + Mstmp119*Mstmp34 + Mstmp120*Mstmp34 + Mstmp121*Mstmp19 + Mstmp121*M[12] + Mstmp125*Mstmp30 + Mstmp125*M[14] + Mstmp177*Mstmp9 + Mstmp177*M[9] + Mstmp200*Mstmp5 + Mstmp200*M[6] + Mstmp206*Mstmp8 + Mstmp206*M[8] + Mstmp207*z + Mstmp208*z + Mstmp209*z + Mstmp210*y + Mstmp34*M[39] + Mstmp41*M[37] + Mstmp45*Mstmp90 + Mstmp45*Mstmp91 + Mstmp45*Mstmp92 + Mstmp45*M[31] + x*M[67] + y*M[63] + z*M[62] + M[95];
#pragma omp atomic
Ms[96] += Mstmp11*Mstmp160 + Mstmp11*Mstmp161 + Mstmp11*Mstmp164 + Mstmp11*M[47] + Mstmp112*Mstmp41 + Mstmp113*Mstmp41 + Mstmp115*Mstmp41 + Mstmp121*Mstmp24 + Mstmp121*M[13] + Mstmp2*Mstmp214 + Mstmp200*Mstmp6 + Mstmp200*M[7] + Mstmp210*z + Mstmp211*z + Mstmp212*z + Mstmp213*y + Mstmp214*M[4] + Mstmp265*M[1] + Mstmp41*M[38] + Mstmp45*Mstmp93 + Mstmp45*Mstmp94 + Mstmp45*Mstmp95 + Mstmp45*M[32] + Mstmp46*Mstmp96 + Mstmp47*Mstmp96 + Mstmp48*Mstmp96 + Mstmp96*M[20] + x*M[68] + y*M[64] + z*M[63] + M[96];
#pragma omp atomic
Ms[97] += Mstmp11*Mstmp166 + Mstmp11*M[48] + Mstmp118*Mstmp41 + Mstmp121*M[14] + Mstmp168*Mstmp267 + Mstmp184*Mstmp9 + Mstmp184*M[9] + Mstmp200*M[8] + Mstmp213*z + Mstmp214*M[5] + Mstmp265*M[2] + Mstmp41*M[39] + Mstmp45*Mstmp97 + Mstmp45*M[33] + Mstmp49*Mstmp96 + Mstmp96*M[21] + x*M[69] + z*M[64] + M[97];
#pragma omp atomic
Ms[98] += Mstmp0*Mstmp251 + Mstmp11*Mstmp169 + Mstmp11*M[49] + Mstmp12*Mstmp170 + Mstmp122*Mstmp34 + Mstmp125*M[15] + Mstmp170*M[10] + Mstmp215*y + Mstmp218*M[6] + Mstmp251*M[3] + Mstmp270*M[1] + Mstmp34*M[40] + Mstmp52*Mstmp86 + Mstmp86*M[22] + x*M[70] + y*M[65] + M[98];
#pragma omp atomic
Ms[99] += Mstmp11*Mstmp171 + Mstmp11*Mstmp172 + Mstmp11*Mstmp173 + Mstmp11*M[50] + Mstmp125*Mstmp36 + Mstmp125*M[16] + Mstmp126*Mstmp34 + Mstmp128*Mstmp34 + Mstmp130*Mstmp34 + Mstmp15*Mstmp170 + Mstmp16*Mstmp170 + Mstmp17*Mstmp170 + Mstmp170*M[11] + Mstmp215*z + Mstmp216*z + Mstmp217*z + Mstmp218*Mstmp7 + Mstmp218*M[7] + Mstmp219*y + Mstmp270*M[2] + Mstmp34*M[41] + Mstmp55*Mstmp86 + Mstmp57*Mstmp86 + Mstmp59*Mstmp86 + Mstmp86*M[23] + x*M[71] + y*M[66] + z*M[65] + M[99];
#pragma omp atomic
Ms[100] += Mstmp0*Mstmp254 + Mstmp11*Mstmp174 + Mstmp11*Mstmp175 + Mstmp11*Mstmp176 + Mstmp11*M[51] + Mstmp12*Mstmp177 + Mstmp121*Mstmp33 + Mstmp121*M[15] + Mstmp122*Mstmp41 + Mstmp123*Mstmp41 + Mstmp124*Mstmp41 + Mstmp125*Mstmp39 + Mstmp125*M[17] + Mstmp132*Mstmp34 + Mstmp134*Mstmp34 + Mstmp136*Mstmp34 + Mstmp177*M[10] + Mstmp218*Mstmp8 + Mstmp218*M[8] + Mstmp219*z + Mstmp220*z + Mstmp221*z + Mstmp222*y + Mstmp254*M[3] + Mstmp271*M[1] + Mstmp34*M[42] + Mstmp41*M[40] + Mstmp61*Mstmp86 + Mstmp62*Mstmp86 + Mstmp63*Mstmp86 + Mstmp86*M[24] + x*M[72] + y*M[67] + z*M[66] + M[100];
#pragma omp atomic
Ms[101] += Mstmp0*Mstmp256 + Mstmp11*Mstmp178 + Mstmp11*Mstmp179 + Mstmp11*Mstmp180 + Mstmp11*M[52] + Mstmp121*Mstmp35 + Mstmp121*M[16] + Mstmp125*Mstmp42 + Mstmp125*M[18] + Mstmp126*Mstmp41 + Mstmp127*Mstmp41 + Mstmp129*Mstmp41 + Mstmp138*Mstmp34 + Mstmp139*Mstmp34 + Mstmp140*Mstmp34 + Mstmp15*Mstmp177 + Mstmp177*M[11] + Mstmp214*Mstmp5 + Mstmp214*M[6] + Mstmp222*z + Mstmp223*z + Mstmp224*z + Mstmp225*y + Mstmp256*M[3] + Mstmp271*M[2] + Mstmp34*M[43] + Mstmp41*M[41] + Mstmp52*Mstmp96 + Mstmp53*Mstmp96 + Mstmp54*Mstmp96 + Mstmp96*M[22] + x*M[73] + y*M[68] + z*M[67] + M[101];
#pragma omp atomic
Ms[102] += Mstmp11*Mstmp181 + Mstmp11*Mstmp182 + Mstmp11*Mstmp183 + Mstmp11*M[53] + Mstmp12*Mstmp184 + Mstmp121*Mstmp38 + Mstmp121*M[17] + Mstmp13*Mstmp184 + Mstmp132*Mstmp41 + Mstmp133*Mstmp41 + Mstmp135*Mstmp41 + Mstmp14*Mstmp184 + Mstmp184*M[10] + Mstmp214*Mstmp6 + Mstmp214*M[7] + Mstmp225*z + Mstmp226*z + Mstmp227*z + Mstmp228*y + Mstmp272*M[1] + Mstmp41*M[42] + Mstmp55*Mstmp96 + Mstmp56*Mstmp96 + Mstmp58*Mstmp96 + Mstmp96*M[23] + x*M[74] + y*M[69] + z*M[68] + M[102];
#pragma omp atomic
Ms[103] += Mstmp0*Mstmp259 + Mstmp11*Mstmp185 + Mstmp11*M[54] + Mstmp121*M[18] + Mstmp138*Mstmp41 + Mstmp15*Mstmp184 + Mstmp184*M[11] + Mstmp214*M[8] + Mstmp228*z + Mstmp259*M[3] + Mstmp272*M[2] + Mstmp41*M[43] + Mstmp61*Mstmp96 + Mstmp96*M[24] + x*M[75] + z*M[69] + M[103];
#pragma omp atomic
Ms[104] += Mstmp1*Mstmp251 + Mstmp141*Mstmp34 + Mstmp170*Mstmp18 + Mstmp170*M[12] + Mstmp229*y + Mstmp251*M[4] + Mstmp273*Mstmp274 + Mstmp34*M[44] + Mstmp64*Mstmp86 + Mstmp86*M[25] + x*M[76] + y*M[70] + M[104];
#pragma omp atomic
Ms[105] += Mstmp146*Mstmp34 + Mstmp148*Mstmp34 + Mstmp150*Mstmp34 + Mstmp170*Mstmp23 + Mstmp170*Mstmp25 + Mstmp170*Mstmp27 + Mstmp170*M[13] + Mstmp229*z + Mstmp230*z + Mstmp231*z + Mstmp234*y + Mstmp251*Mstmp3 + Mstmp251*Mstmp4 + Mstmp251*M[5] + Mstmp34*M[45] + Mstmp69*Mstmp86 + Mstmp71*Mstmp86 + Mstmp73*Mstmp86 + Mstmp86*M[26] + x*M[77] + y*M[71] + z*M[70] + M[105];
#pragma omp atomic
Ms[106] += Mstmp1*Mstmp254 + Mstmp141*Mstmp41 + Mstmp142*Mstmp41 + Mstmp143*Mstmp41 + Mstmp152*Mstmp34 + Mstmp154*Mstmp34 + Mstmp156*Mstmp34 + Mstmp170*Mstmp29 + Mstmp170*Mstmp30 + Mstmp170*Mstmp31 + Mstmp170*M[14] + Mstmp177*Mstmp18 + Mstmp177*M[12] + Mstmp234*z + Mstmp235*z + Mstmp236*z + Mstmp237*y + Mstmp254*M[4] + Mstmp275*M[0] + Mstmp34*M[46] + Mstmp41*M[44] + Mstmp75*Mstmp86 + Mstmp77*Mstmp86 + Mstmp79*Mstmp86 + Mstmp86*M[27] + x*M[78] + y*M[72] + z*M[71] + M[106];
#pragma omp atomic
Ms[107] += Mstmp1*Mstmp256 + Mstmp146*Mstmp41 + Mstmp147*Mstmp41 + Mstmp149*Mstmp41 + Mstmp159*Mstmp34 + Mstmp161*Mstmp34 + Mstmp163*Mstmp34 + Mstmp177*Mstmp23 + Mstmp177*M[13] + Mstmp205*Mstmp276 + Mstmp237*z + Mstmp238*z + Mstmp239*z + Mstmp241*y + Mstmp254*Mstmp3 + Mstmp254*M[5] + Mstmp256*M[4] + Mstmp34*M[47] + Mstmp41*M[45] + Mstmp64*Mstmp96 + Mstmp65*Mstmp96 + Mstmp66*Mstmp96 + Mstmp81*Mstmp86 + Mstmp82*Mstmp86 + Mstmp83*Mstmp86 + Mstmp86*M[28] + Mstmp96*M[25] + x*M[79] + y*M[73] + z*M[72] + M[107];
#pragma omp atomic
Ms[108] += Mstmp111*Mstmp277 + Mstmp152*Mstmp41 + Mstmp153*Mstmp41 + Mstmp155*Mstmp41 + Mstmp165*Mstmp34 + Mstmp166*Mstmp34 + Mstmp167*Mstmp34 + Mstmp177*Mstmp29 + Mstmp177*M[14] + Mstmp18*Mstmp184 + Mstmp184*Mstmp19 + Mstmp184*Mstmp20 + Mstmp184*M[12] + Mstmp241*z + Mstmp242*z + Mstmp243*z + Mstmp245*y + Mstmp256*Mstmp3 + Mstmp256*M[5] + Mstmp34*M[48] + Mstmp41*M[46] + Mstmp69*Mstmp96 + Mstmp70*Mstmp96 + Mstmp72*Mstmp96 + Mstmp96*M[26] + x*M[80] + y*M[74] + z*M[73] + M[108];
#pragma omp atomic
Ms[109] += Mstmp1*Mstmp259 + Mstmp159*Mstmp41 + Mstmp160*Mstmp41 + Mstmp162*Mstmp41 + Mstmp184*Mstmp23 + Mstmp184*Mstmp24 + Mstmp184*Mstmp26 + Mstmp184*M[13] + Mstmp2*Mstmp259 + Mstmp245*z + Mstmp246*z + Mstmp247*z + Mstmp248*y + Mstmp259*M[4] + Mstmp41*M[47] + Mstmp75*Mstmp96 + Mstmp76*Mstmp96 + Mstmp78*Mstmp96 + Mstmp96*M[27] + x*M[81] + y*M[75] + z*M[74] + M[109];
#pragma omp atomic
Ms[110] += Mstmp165*Mstmp41 + Mstmp184*Mstmp29 + Mstmp184*M[14] + Mstmp248*z + Mstmp259*Mstmp3 + Mstmp259*M[5] + Mstmp274*Mstmp278 + Mstmp41*M[48] + Mstmp81*Mstmp96 + Mstmp96*M[28] + x*M[82] + z*M[75] + M[110];
#pragma omp atomic
Ms[111] += Mstmp170*M[15] + Mstmp251*M[6] + Mstmp279*M[1] + Mstmp34*M[49] + Mstmp86*M[29] + y*M[76] + M[111];
#pragma omp atomic
Ms[112] += Mstmp170*Mstmp36 + Mstmp170*M[16] + Mstmp172*Mstmp34 + Mstmp250*z + Mstmp251*Mstmp7 + Mstmp251*M[7] + Mstmp279*M[2] + Mstmp34*M[50] + Mstmp86*Mstmp88 + Mstmp86*M[30] + y*M[77] + z*M[76] + M[112];
#pragma omp atomic
Ms[113] += Mstmp169*Mstmp41 + Mstmp170*Mstmp39 + Mstmp170*M[17] + Mstmp175*Mstmp34 + Mstmp177*M[15] + Mstmp251*Mstmp8 + Mstmp251*M[8] + Mstmp252*z + Mstmp254*M[6] + Mstmp275*M[1] + Mstmp34*M[51] + Mstmp41*M[49] + Mstmp86*Mstmp91 + Mstmp86*M[31] + y*M[78] + z*M[77] + M[113];
#pragma omp atomic
Ms[114] += Mstmp170*Mstmp42 + Mstmp170*M[18] + Mstmp171*Mstmp41 + Mstmp177*M[16] + Mstmp179*Mstmp34 + Mstmp253*z + Mstmp254*M[7] + Mstmp256*M[6] + Mstmp275*M[2] + Mstmp280*M[1] + Mstmp34*M[52] + Mstmp41*M[50] + Mstmp85*Mstmp96 + Mstmp86*Mstmp94 + Mstmp86*M[32] + Mstmp96*M[29] + y*M[79] + z*M[78] + M[114];
#pragma omp atomic
Ms[115] += Mstmp174*Mstmp41 + Mstmp177*M[17] + Mstmp182*Mstmp34 + Mstmp184*Mstmp33 + Mstmp184*M[15] + Mstmp254*M[8] + Mstmp255*z + Mstmp256*M[7] + Mstmp280*M[2] + Mstmp281*M[1] + Mstmp34*M[53] + Mstmp41*M[51] + Mstmp86*Mstmp97 + Mstmp86*M[33] + Mstmp87*Mstmp96 + Mstmp96*M[30] + y*M[80] + z*M[79] + M[115];
#pragma omp atomic
Ms[116] += Mstmp177*M[18] + Mstmp178*Mstmp41 + Mstmp184*Mstmp35 + Mstmp184*M[16] + Mstmp185*Mstmp34 + Mstmp256*M[8] + Mstmp257*z + Mstmp259*Mstmp5 + Mstmp259*M[6] + Mstmp281*M[2] + Mstmp34*M[54] + Mstmp41*M[52] + Mstmp90*Mstmp96 + Mstmp96*M[31] + y*M[81] + z*M[80] + M[116];
#pragma omp atomic
Ms[117] += Mstmp181*Mstmp41 + Mstmp184*Mstmp38 + Mstmp184*M[17] + Mstmp258*z + Mstmp259*Mstmp6 + Mstmp259*M[7] + Mstmp282*M[1] + Mstmp41*M[53] + Mstmp93*Mstmp96 + Mstmp96*M[32] + y*M[82] + z*M[81] + M[117];
#pragma omp atomic
Ms[118] += Mstmp184*M[18] + Mstmp259*M[8] + Mstmp282*M[2] + Mstmp41*M[54] + Mstmp96*M[33] + z*M[82] + M[118];

}

void M2L_7(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[119];
double Dtmp0 = pow(R, -3);
double Dtmp1 = 1.0*Dtmp0;
double Dtmp2 = pow(x, 2);
double Dtmp3 = pow(R, -2);
double Dtmp4 = 3.0*Dtmp3;
double Dtmp5 = pow(R, -5);
double Dtmp6 = Dtmp5*x;
double Dtmp7 = 3.0*Dtmp6;
double Dtmp8 = pow(y, 2);
double Dtmp9 = Dtmp5*y;
double Dtmp10 = 15.0*Dtmp3;
double Dtmp11 = -Dtmp10*Dtmp2;
double Dtmp12 = Dtmp11 + 3.0;
double Dtmp13 = Dtmp12*Dtmp5;
double Dtmp14 = -Dtmp10*Dtmp8;
double Dtmp15 = Dtmp14 + 3.0;
double Dtmp16 = pow(R, -7);
double Dtmp17 = Dtmp16*x;
double Dtmp18 = y*z;
double Dtmp19 = pow(x, 4);
double Dtmp20 = pow(R, -4);
double Dtmp21 = 105.0*Dtmp20;
double Dtmp22 = Dtmp2*Dtmp3;
double Dtmp23 = -105.0*Dtmp22;
double Dtmp24 = Dtmp23 + 45.0;
double Dtmp25 = Dtmp17*Dtmp24;
double Dtmp26 = Dtmp2*Dtmp8;
double Dtmp27 = Dtmp23 + 15.0;
double Dtmp28 = Dtmp16*y;
double Dtmp29 = Dtmp28*z;
double Dtmp30 = Dtmp3*Dtmp8;
double Dtmp31 = -105.0*Dtmp30;
double Dtmp32 = Dtmp31 + 45.0;
double Dtmp33 = 1.0*Dtmp17;
double Dtmp34 = pow(y, 4);
double Dtmp35 = 945.0*Dtmp20;
double Dtmp36 = Dtmp19*Dtmp35;
double Dtmp37 = 630.0*Dtmp22;
double Dtmp38 = Dtmp16*(Dtmp36 - Dtmp37 + 45.0);
double Dtmp39 = 315.0*Dtmp30;
double Dtmp40 = Dtmp26*Dtmp35;
double Dtmp41 = 315.0 - 945.0*Dtmp22;
double Dtmp42 = pow(R, -9);
double Dtmp43 = Dtmp42*x;
double Dtmp44 = Dtmp43*y;
double Dtmp45 = Dtmp44*z;
double Dtmp46 = 315.0*Dtmp22;
double Dtmp47 = Dtmp16*z;
double Dtmp48 = Dtmp34*Dtmp35;
double Dtmp49 = 630.0*Dtmp30;
double Dtmp50 = Dtmp48 - Dtmp49 + 45.0;
double Dtmp51 = 315.0 - 945.0*Dtmp30;
double Dtmp52 = pow(x, 6);
double Dtmp53 = pow(R, -6);
double Dtmp54 = 10395.0*Dtmp53;
double Dtmp55 = Dtmp19*Dtmp20;
double Dtmp56 = 10395.0*Dtmp55;
double Dtmp57 = -9450.0*Dtmp22 + Dtmp56 + 1575.0;
double Dtmp58 = Dtmp43*Dtmp57;
double Dtmp59 = Dtmp19*Dtmp8;
double Dtmp60 = Dtmp20*Dtmp26;
double Dtmp61 = -5670.0*Dtmp60 - 45.0;
double Dtmp62 = -5670.0*Dtmp22 + Dtmp56 + 315.0;
double Dtmp63 = Dtmp42*y;
double Dtmp64 = Dtmp63*z;
double Dtmp65 = -2835.0*Dtmp30;
double Dtmp66 = 10395.0*Dtmp60;
double Dtmp67 = Dtmp65 + Dtmp66;
double Dtmp68 = -2835.0*Dtmp22;
double Dtmp69 = Dtmp68 + 945.0;
double Dtmp70 = Dtmp43*z;
double Dtmp71 = Dtmp2*Dtmp34;
double Dtmp72 = Dtmp20*Dtmp34;
double Dtmp73 = 10395.0*Dtmp72;
double Dtmp74 = -9450.0*Dtmp30 + Dtmp73 + 1575.0;
double Dtmp75 = -5670.0*Dtmp30 + Dtmp73 + 315.0;
double Dtmp76 = pow(y, 6);
double Dtmp77 = 135135.0*Dtmp53;
double Dtmp78 = -Dtmp52*Dtmp77;
double Dtmp79 = Dtmp42*(-42525.0*Dtmp22 + 155925.0*Dtmp55 + Dtmp78 + 1575.0);
double Dtmp80 = 103950.0*Dtmp60;
double Dtmp81 = -Dtmp59*Dtmp77;
double Dtmp82 = Dtmp18*x/pow(R, 11);
double Dtmp83 = 62370.0*Dtmp60;
double Dtmp84 = Dtmp65 + Dtmp81 + Dtmp83;
double Dtmp85 = Dtmp42*z;
double Dtmp86 = -Dtmp71*Dtmp77;
double Dtmp87 = Dtmp83 + Dtmp86;
double Dtmp88 = -Dtmp76*Dtmp77;
double Dtmp89 = -42525.0*Dtmp30 + 155925.0*Dtmp72 + Dtmp88 + 1575.0;
D[0] = -Dtmp1*x;
D[1] = -Dtmp1*y;
D[2] = -Dtmp1*z;
D[3] = Dtmp0*(Dtmp2*Dtmp4 - 1.0);
D[4] = Dtmp7*y;
D[5] = Dtmp7*z;
D[6] = Dtmp0*(Dtmp4*Dtmp8 - 1.0);
D[7] = 3.0*Dtmp9*z;
D[8] = -D[3] - D[6];
D[9] = Dtmp6*(Dtmp11 + 9.0);
D[10] = Dtmp13*y;
D[11] = Dtmp13*z;
D[12] = 1.0*Dtmp15*Dtmp6;
D[13] = -15.0*Dtmp17*Dtmp18;
D[14] = -D[9] - D[12];
D[15] = Dtmp9*(Dtmp14 + 9.0);
D[16] = Dtmp15*Dtmp5*z;
D[17] = -D[10] - D[15];
D[18] = -D[11] - D[16];
D[19] = Dtmp5*(Dtmp19*Dtmp21 - 90.0*Dtmp22 + 9.0);
D[20] = -Dtmp25*y;
D[21] = -Dtmp25*z;
D[22] = Dtmp5*(Dtmp12 + Dtmp14 + Dtmp21*Dtmp26);
D[23] = -Dtmp27*Dtmp29;
D[24] = -D[19] - D[22];
D[25] = -Dtmp32*Dtmp33*y;
D[26] = -Dtmp33*z*(Dtmp31 + 15.0);
D[27] = -D[20] - D[25];
D[28] = -D[21] - D[26];
D[29] = Dtmp5*(Dtmp21*Dtmp34 - 90.0*Dtmp30 + 9.0);
D[30] = -Dtmp29*Dtmp32;
D[31] = -D[22] - D[29];
D[32] = -D[23] - D[30];
D[33] = -D[24] - D[31];
D[34] = -Dtmp17*(-1050.0*Dtmp22 + Dtmp36 + 225.0);
D[35] = -Dtmp38*y;
D[36] = -Dtmp38*z;
D[37] = -Dtmp17*(Dtmp24 - Dtmp39 + Dtmp40);
D[38] = Dtmp41*Dtmp45;
D[39] = -D[34] - D[37];
D[40] = -Dtmp28*(Dtmp32 + Dtmp40 - Dtmp46);
D[41] = -Dtmp47*(Dtmp27 + Dtmp31 + Dtmp40);
D[42] = -D[35] - D[40];
D[43] = -D[36] - D[41];
D[44] = -Dtmp33*Dtmp50;
D[45] = 1.0*Dtmp45*Dtmp51;
D[46] = -D[37] - D[44];
D[47] = -D[38] - D[45];
D[48] = -D[39] - D[46];
D[49] = -Dtmp28*(-1050.0*Dtmp30 + Dtmp48 + 225.0);
D[50] = -Dtmp47*Dtmp50;
D[51] = -D[40] - D[49];
D[52] = -D[41] - D[50];
D[53] = -D[42] - D[51];
D[54] = -D[43] - D[52];
D[55] = Dtmp16*(4725.0*Dtmp22 + Dtmp52*Dtmp54 - 14175.0*Dtmp55 - 225.0);
D[56] = Dtmp58*y;
D[57] = Dtmp58*z;
D[58] = Dtmp16*(-Dtmp36 + Dtmp37 + Dtmp39 + Dtmp54*Dtmp59 + Dtmp61);
D[59] = Dtmp62*Dtmp64;
D[60] = -D[55] - D[58];
D[61] = Dtmp44*(Dtmp67 + Dtmp69);
D[62] = Dtmp70*(Dtmp41 + Dtmp67);
D[63] = -D[56] - D[61];
D[64] = -D[57] - D[62];
D[65] = Dtmp16*(Dtmp46 - Dtmp48 + Dtmp49 + Dtmp54*Dtmp71 + Dtmp61);
D[66] = Dtmp64*(Dtmp51 + Dtmp66 + Dtmp68);
D[67] = -D[58] - D[65];
D[68] = -D[59] - D[66];
D[69] = -D[60] - D[67];
D[70] = 1.0*Dtmp44*Dtmp74;
D[71] = 1.0*Dtmp70*Dtmp75;
D[72] = -D[61] - D[70];
D[73] = -D[62] - D[71];
D[74] = -D[63] - D[72];
D[75] = -D[64] - D[73];
D[76] = Dtmp16*(4725.0*Dtmp30 + Dtmp54*Dtmp76 - 14175.0*Dtmp72 - 225.0);
D[77] = Dtmp64*Dtmp74;
D[78] = -D[65] - D[76];
D[79] = -D[66] - D[77];
D[80] = -D[67] - D[78];
D[81] = -D[68] - D[79];
D[82] = -D[69] - D[80];
D[83] = Dtmp43*(-99225.0*Dtmp22 + 218295.0*Dtmp55 + Dtmp78 + 11025.0);
D[84] = Dtmp79*y;
D[85] = Dtmp79*z;
D[86] = Dtmp43*(-14175.0*Dtmp30 + Dtmp57 + Dtmp80 + Dtmp81);
D[87] = -Dtmp82*(-103950.0*Dtmp22 + 135135.0*Dtmp55 + 14175.0);
D[88] = -D[83] - D[86];
D[89] = Dtmp63*(-17010.0*Dtmp22 + 31185.0*Dtmp55 + Dtmp84 + 945.0);
D[90] = Dtmp85*(Dtmp62 + Dtmp84);
D[91] = -D[84] - D[89];
D[92] = -D[85] - D[90];
D[93] = Dtmp43*(-17010.0*Dtmp30 + Dtmp69 + 31185.0*Dtmp72 + Dtmp87);
D[94] = -Dtmp82*(-31185.0*Dtmp22 - 31185.0*Dtmp30 + 135135.0*Dtmp60 + 8505.0);
D[95] = -D[86] - D[93];
D[96] = -D[87] - D[94];
D[97] = -D[88] - D[95];
D[98] = Dtmp63*(-14175.0*Dtmp22 + Dtmp74 + Dtmp80 + Dtmp86);
D[99] = Dtmp85*(Dtmp68 + Dtmp75 + Dtmp87);
D[100] = -D[89] - D[98];
D[101] = -D[90] - D[99];
D[102] = -D[91] - D[100];
D[103] = -D[92] - D[101];
D[104] = 1.0*Dtmp43*Dtmp89;
D[105] = -1.0*Dtmp82*(-103950.0*Dtmp30 + 135135.0*Dtmp72 + 14175.0);
D[106] = -D[93] - D[104];
D[107] = -D[94] - D[105];
D[108] = -D[95] - D[106];
D[109] = -D[96] - D[107];
D[110] = -D[97] - D[108];
D[111] = Dtmp63*(-99225.0*Dtmp30 + 218295.0*Dtmp72 + Dtmp88 + 11025.0);
D[112] = Dtmp85*Dtmp89;
D[113] = -D[98] - D[111];
D[114] = -D[99] - D[112];
D[115] = -D[100] - D[113];
D[116] = -D[101] - D[114];
D[117] = -D[102] - D[115];
D[118] = -D[103] - D[116];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8] + D[9]*M[9] + D[10]*M[10] + D[11]*M[11] + D[12]*M[12] + D[13]*M[13] + D[14]*M[14] + D[15]*M[15] + D[16]*M[16] + D[17]*M[17] + D[18]*M[18] + D[19]*M[19] + D[20]*M[20] + D[21]*M[21] + D[22]*M[22] + D[23]*M[23] + D[24]*M[24] + D[25]*M[25] + D[26]*M[26] + D[27]*M[27] + D[28]*M[28] + D[29]*M[29] + D[30]*M[30] + D[31]*M[31] + D[32]*M[32] + D[33]*M[33] + D[34]*M[34] + D[35]*M[35] + D[36]*M[36] + D[37]*M[37] + D[38]*M[38] + D[39]*M[39] + D[40]*M[40] + D[41]*M[41] + D[42]*M[42] + D[43]*M[43] + D[44]*M[44] + D[45]*M[45] + D[46]*M[46] + D[47]*M[47] + D[48]*M[48] + D[49]*M[49] + D[50]*M[50] + D[51]*M[51] + D[52]*M[52] + D[53]*M[53] + D[54]*M[54] + D[55]*M[55] + D[56]*M[56] + D[57]*M[57] + D[58]*M[58] + D[59]*M[59] + D[60]*M[60] + D[61]*M[61] + D[62]*M[62] + D[63]*M[63] + D[64]*M[64] + D[65]*M[65] + D[66]*M[66] + D[67]*M[67] + D[68]*M[68] + D[69]*M[69] + D[70]*M[70] + D[71]*M[71] + D[72]*M[72] + D[73]*M[73] + D[74]*M[74] + D[75]*M[75] + D[76]*M[76] + D[77]*M[77] + D[78]*M[78] + D[79]*M[79] + D[80]*M[80] + D[81]*M[81] + D[82]*M[82] + D[83]*M[83] + D[84]*M[84] + D[85]*M[85] + D[86]*M[86] + D[87]*M[87] + D[88]*M[88] + D[89]*M[89] + D[90]*M[90] + D[91]*M[91] + D[92]*M[92] + D[93]*M[93] + D[94]*M[94] + D[95]*M[95] + D[96]*M[96] + D[97]*M[97] + D[98]*M[98] + D[99]*M[99] + D[100]*M[100] + D[101]*M[101] + D[102]*M[102] + D[103]*M[103] + D[104]*M[104] + D[105]*M[105] + D[106]*M[106] + D[107]*M[107] + D[108]*M[108] + D[109]*M[109] + D[110]*M[110] + D[111]*M[111] + D[112]*M[112] + D[113]*M[113] + D[114]*M[114] + D[115]*M[115] + D[116]*M[116] + D[117]*M[117] + D[118]*M[118];
#pragma omp atomic
L[1] += D[3]*M[0] + D[4]*M[1] + D[5]*M[2] + D[9]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[14]*M[8] + D[19]*M[9] + D[20]*M[10] + D[21]*M[11] + D[22]*M[12] + D[23]*M[13] + D[24]*M[14] + D[25]*M[15] + D[26]*M[16] + D[27]*M[17] + D[28]*M[18] + D[34]*M[19] + D[35]*M[20] + D[36]*M[21] + D[37]*M[22] + D[38]*M[23] + D[39]*M[24] + D[40]*M[25] + D[41]*M[26] + D[42]*M[27] + D[43]*M[28] + D[44]*M[29] + D[45]*M[30] + D[46]*M[31] + D[47]*M[32] + D[48]*M[33] + D[55]*M[34] + D[56]*M[35] + D[57]*M[36] + D[58]*M[37] + D[59]*M[38] + D[60]*M[39] + D[61]*M[40] + D[62]*M[41] + D[63]*M[42] + D[64]*M[43] + D[65]*M[44] + D[66]*M[45] + D[67]*M[46] + D[68]*M[47] + D[69]*M[48] + D[70]*M[49] + D[71]*M[50] + D[72]*M[51] + D[73]*M[52] + D[74]*M[53] + D[75]*M[54] + D[83]*M[55] + D[84]*M[56] + D[85]*M[57] + D[86]*M[58] + D[87]*M[59] + D[88]*M[60] + D[89]*M[61] + D[90]*M[62] + D[91]*M[63] + D[92]*M[64] + D[93]*M[65] + D[94]*M[66] + D[95]*M[67] + D[96]*M[68] + D[97]*M[69] + D[98]*M[70] + D[99]*M[71] + D[100]*M[72] + D[101]*M[73] + D[102]*M[74] + D[103]*M[75] + D[104]*M[76] + D[105]*M[77] + D[106]*M[78] + D[107]*M[79] + D[108]*M[80] + D[109]*M[81] + D[110]*M[82];
#pragma omp atomic
L[2] += D[4]*M[0] + D[6]*M[1] + D[7]*M[2] + D[10]*M[3] + D[12]*M[4] + D[13]*M[5] + D[15]*M[6] + D[16]*M[7] + D[17]*M[8] + D[20]*M[9] + D[22]*M[10] + D[23]*M[11] + D[25]*M[12] + D[26]*M[13] + D[27]*M[14] + D[29]*M[15] + D[30]*M[16] + D[31]*M[17] + D[32]*M[18] + D[35]*M[19] + D[37]*M[20] + D[38]*M[21] + D[40]*M[22] + D[41]*M[23] + D[42]*M[24] + D[44]*M[25] + D[45]*M[26] + D[46]*M[27] + D[47]*M[28] + D[49]*M[29] + D[50]*M[30] + D[51]*M[31] + D[52]*M[32] + D[53]*M[33] + D[56]*M[34] + D[58]*M[35] + D[59]*M[36] + D[61]*M[37] + D[62]*M[38] + D[63]*M[39] + D[65]*M[40] + D[66]*M[41] + D[67]*M[42] + D[68]*M[43] + D[70]*M[44] + D[71]*M[45] + D[72]*M[46] + D[73]*M[47] + D[74]*M[48] + D[76]*M[49] + D[77]*M[50] + D[78]*M[51] + D[79]*M[52] + D[80]*M[53] + D[81]*M[54] + D[84]*M[55] + D[86]*M[56] + D[87]*M[57] + D[89]*M[58] + D[90]*M[59] + D[91]*M[60] + D[93]*M[61] + D[94]*M[62] + D[95]*M[63] + D[96]*M[64] + D[98]*M[65] + D[99]*M[66] + D[100]*M[67] + D[101]*M[68] + D[102]*M[69] + D[104]*M[70] + D[105]*M[71] + D[106]*M[72] + D[107]*M[73] + D[108]*M[74] + D[109]*M[75] + D[111]*M[76] + D[112]*M[77] + D[113]*M[78] + D[114]*M[79] + D[115]*M[80] + D[116]*M[81] + D[117]*M[82];
#pragma omp atomic
L[3] += D[5]*M[0] + D[7]*M[1] + D[8]*M[2] + D[11]*M[3] + D[13]*M[4] + D[14]*M[5] + D[16]*M[6] + D[17]*M[7] + D[18]*M[8] + D[21]*M[9] + D[23]*M[10] + D[24]*M[11] + D[26]*M[12] + D[27]*M[13] + D[28]*M[14] + D[30]*M[15] + D[31]*M[16] + D[32]*M[17] + D[33]*M[18] + D[36]*M[19] + D[38]*M[20] + D[39]*M[21] + D[41]*M[22] + D[42]*M[23] + D[43]*M[24] + D[45]*M[25] + D[46]*M[26] + D[47]*M[27] + D[48]*M[28] + D[50]*M[29] + D[51]*M[30] + D[52]*M[31] + D[53]*M[32] + D[54]*M[33] + D[57]*M[34] + D[59]*M[35] + D[60]*M[36] + D[62]*M[37] + D[63]*M[38] + D[64]*M[39] + D[66]*M[40] + D[67]*M[41] + D[68]*M[42] + D[69]*M[43] + D[71]*M[44] + D[72]*M[45] + D[73]*M[46] + D[74]*M[47] + D[75]*M[48] + D[77]*M[49] + D[78]*M[50] + D[79]*M[51] + D[80]*M[52] + D[81]*M[53] + D[82]*M[54] + D[85]*M[55] + D[87]*M[56] + D[88]*M[57] + D[90]*M[58] + D[91]*M[59] + D[92]*M[60] + D[94]*M[61] + D[95]*M[62] + D[96]*M[63] + D[97]*M[64] + D[99]*M[65] + D[100]*M[66] + D[101]*M[67] + D[102]*M[68] + D[103]*M[69] + D[105]*M[70] + D[106]*M[71] + D[107]*M[72] + D[108]*M[73] + D[109]*M[74] + D[110]*M[75] + D[112]*M[76] + D[113]*M[77] + D[114]*M[78] + D[115]*M[79] + D[116]*M[80] + D[117]*M[81] + D[118]*M[82];
#pragma omp atomic
L[4] += D[9]*M[0] + D[10]*M[1] + D[11]*M[2] + D[19]*M[3] + D[20]*M[4] + D[21]*M[5] + D[22]*M[6] + D[23]*M[7] + D[24]*M[8] + D[34]*M[9] + D[35]*M[10] + D[36]*M[11] + D[37]*M[12] + D[38]*M[13] + D[39]*M[14] + D[40]*M[15] + D[41]*M[16] + D[42]*M[17] + D[43]*M[18] + D[55]*M[19] + D[56]*M[20] + D[57]*M[21] + D[58]*M[22] + D[59]*M[23] + D[60]*M[24] + D[61]*M[25] + D[62]*M[26] + D[63]*M[27] + D[64]*M[28] + D[65]*M[29] + D[66]*M[30] + D[67]*M[31] + D[68]*M[32] + D[69]*M[33] + D[83]*M[34] + D[84]*M[35] + D[85]*M[36] + D[86]*M[37] + D[87]*M[38] + D[88]*M[39] + D[89]*M[40] + D[90]*M[41] + D[91]*M[42] + D[92]*M[43] + D[93]*M[44] + D[94]*M[45] + D[95]*M[46] + D[96]*M[47] + D[97]*M[48] + D[98]*M[49] + D[99]*M[50] + D[100]*M[51] + D[101]*M[52] + D[102]*M[53] + D[103]*M[54];
#pragma omp atomic
L[5] += D[10]*M[0] + D[12]*M[1] + D[13]*M[2] + D[20]*M[3] + D[22]*M[4] + D[23]*M[5] + D[25]*M[6] + D[26]*M[7] + D[27]*M[8] + D[35]*M[9] + D[37]*M[10] + D[38]*M[11] + D[40]*M[12] + D[41]*M[13] + D[42]*M[14] + D[44]*M[15] + D[45]*M[16] + D[46]*M[17] + D[47]*M[18] + D[56]*M[19] + D[58]*M[20] + D[59]*M[21] + D[61]*M[22] + D[62]*M[23] + D[63]*M[24] + D[65]*M[25] + D[66]*M[26] + D[67]*M[27] + D[68]*M[28] + D[70]*M[29] + D[71]*M[30] + D[72]*M[31] + D[73]*M[32] + D[74]*M[33] + D[84]*M[34] + D[86]*M[35] + D[87]*M[36] + D[89]*M[37] + D[90]*M[38] + D[91]*M[39] + D[93]*M[40] + D[94]*M[41] + D[95]*M[42] + D[96]*M[43] + D[98]*M[44] + D[99]*M[45] + D[100]*M[46] + D[101]*M[47] + D[102]*M[48] + D[104]*M[49] + D[105]*M[50] + D[106]*M[51] + D[107]*M[52] + D[108]*M[53] + D[109]*M[54];
#pragma omp atomic
L[6] += D[11]*M[0] + D[13]*M[1] + D[14]*M[2] + D[21]*M[3] + D[23]*M[4] + D[24]*M[5] + D[26]*M[6] + D[27]*M[7] + D[28]*M[8] + D[36]*M[9] + D[38]*M[10] + D[39]*M[11] + D[41]*M[12] + D[42]*M[13] + D[43]*M[14] + D[45]*M[15] + D[46]*M[16] + D[47]*M[17] + D[48]*M[18] + D[57]*M[19] + D[59]*M[20] + D[60]*M[21] + D[62]*M[22] + D[63]*M[23] + D[64]*M[24] + D[66]*M[25] + D[67]*M[26] + D[68]*M[27] + D[69]*M[28] + D[71]*M[29] + D[72]*M[30] + D[73]*M[31] + D[74]*M[32] + D[75]*M[33] + D[85]*M[34] + D[87]*M[35] + D[88]*M[36] + D[90]*M[37] + D[91]*M[38] + D[92]*M[39] + D[94]*M[40] + D[95]*M[41] + D[96]*M[42] + D[97]*M[43] + D[99]*M[44] + D[100]*M[45] + D[101]*M[46] + D[102]*M[47] + D[103]*M[48] + D[105]*M[49] + D[106]*M[50] + D[107]*M[51] + D[108]*M[52] + D[109]*M[53] + D[110]*M[54];
#pragma omp atomic
L[7] += D[12]*M[0] + D[15]*M[1] + D[16]*M[2] + D[22]*M[3] + D[25]*M[4] + D[26]*M[5] + D[29]*M[6] + D[30]*M[7] + D[31]*M[8] + D[37]*M[9] + D[40]*M[10] + D[41]*M[11] + D[44]*M[12] + D[45]*M[13] + D[46]*M[14] + D[49]*M[15] + D[50]*M[16] + D[51]*M[17] + D[52]*M[18] + D[58]*M[19] + D[61]*M[20] + D[62]*M[21] + D[65]*M[22] + D[66]*M[23] + D[67]*M[24] + D[70]*M[25] + D[71]*M[26] + D[72]*M[27] + D[73]*M[28] + D[76]*M[29] + D[77]*M[30] + D[78]*M[31] + D[79]*M[32] + D[80]*M[33] + D[86]*M[34] + D[89]*M[35] + D[90]*M[36] + D[93]*M[37] + D[94]*M[38] + D[95]*M[39] + D[98]*M[40] + D[99]*M[41] + D[100]*M[42] + D[101]*M[43] + D[104]*M[44] + D[105]*M[45] + D[106]*M[46] + D[107]*M[47] + D[108]*M[48] + D[111]*M[49] + D[112]*M[50] + D[113]*M[51] + D[114]*M[52] + D[115]*M[53] + D[116]*M[54];
#pragma omp atomic
L[8] += D[13]*M[0] + D[16]*M[1] + D[17]*M[2] + D[23]*M[3] + D[26]*M[4] + D[27]*M[5] + D[30]*M[6] + D[31]*M[7] + D[32]*M[8] + D[38]*M[9] + D[41]*M[10] + D[42]*M[11] + D[45]*M[12] + D[46]*M[13] + D[47]*M[14] + D[50]*M[15] + D[51]*M[16] + D[52]*M[17] + D[53]*M[18] + D[59]*M[19] + D[62]*M[20] + D[63]*M[21] + D[66]*M[22] + D[67]*M[23] + D[68]*M[24] + D[71]*M[25] + D[72]*M[26] + D[73]*M[27] + D[74]*M[28] + D[77]*M[29] + D[78]*M[30] + D[79]*M[31] + D[80]*M[32] + D[81]*M[33] + D[87]*M[34] + D[90]*M[35] + D[91]*M[36] + D[94]*M[37] + D[95]*M[38] + D[96]*M[39] + D[99]*M[40] + D[100]*M[41] + D[101]*M[42] + D[102]*M[43] + D[105]*M[44] + D[106]*M[45] + D[107]*M[46] + D[108]*M[47] + D[109]*M[48] + D[112]*M[49] + D[113]*M[50] + D[114]*M[51] + D[115]*M[52] + D[116]*M[53] + D[117]*M[54];
#pragma omp atomic
L[9] += D[14]*M[0] + D[17]*M[1] + D[18]*M[2] + D[24]*M[3] + D[27]*M[4] + D[28]*M[5] + D[31]*M[6] + D[32]*M[7] + D[33]*M[8] + D[39]*M[9] + D[42]*M[10] + D[43]*M[11] + D[46]*M[12] + D[47]*M[13] + D[48]*M[14] + D[51]*M[15] + D[52]*M[16] + D[53]*M[17] + D[54]*M[18] + D[60]*M[19] + D[63]*M[20] + D[64]*M[21] + D[67]*M[22] + D[68]*M[23] + D[69]*M[24] + D[72]*M[25] + D[73]*M[26] + D[74]*M[27] + D[75]*M[28] + D[78]*M[29] + D[79]*M[30] + D[80]*M[31] + D[81]*M[32] + D[82]*M[33] + D[88]*M[34] + D[91]*M[35] + D[92]*M[36] + D[95]*M[37] + D[96]*M[38] + D[97]*M[39] + D[100]*M[40] + D[101]*M[41] + D[102]*M[42] + D[103]*M[43] + D[106]*M[44] + D[107]*M[45] + D[108]*M[46] + D[109]*M[47] + D[110]*M[48] + D[113]*M[49] + D[114]*M[50] + D[115]*M[51] + D[116]*M[52] + D[117]*M[53] + D[118]*M[54];
#pragma omp atomic
L[10] += D[19]*M[0] + D[20]*M[1] + D[21]*M[2] + D[34]*M[3] + D[35]*M[4] + D[36]*M[5] + D[37]*M[6] + D[38]*M[7] + D[39]*M[8] + D[55]*M[9] + D[56]*M[10] + D[57]*M[11] + D[58]*M[12] + D[59]*M[13] + D[60]*M[14] + D[61]*M[15] + D[62]*M[16] + D[63]*M[17] + D[64]*M[18] + D[83]*M[19] + D[84]*M[20] + D[85]*M[21] + D[86]*M[22] + D[87]*M[23] + D[88]*M[24] + D[89]*M[25] + D[90]*M[26] + D[91]*M[27] + D[92]*M[28] + D[93]*M[29] + D[94]*M[30] + D[95]*M[31] + D[96]*M[32] + D[97]*M[33];
#pragma omp atomic
L[11] += D[20]*M[0] + D[22]*M[1] + D[23]*M[2] + D[35]*M[3] + D[37]*M[4] + D[38]*M[5] + D[40]*M[6] + D[41]*M[7] + D[42]*M[8] + D[56]*M[9] + D[58]*M[10] + D[59]*M[11] + D[61]*M[12] + D[62]*M[13] + D[63]*M[14] + D[65]*M[15] + D[66]*M[16] + D[67]*M[17] + D[68]*M[18] + D[84]*M[19] + D[86]*M[20] + D[87]*M[21] + D[89]*M[22] + D[90]*M[23] + D[91]*M[24] + D[93]*M[25] + D[94]*M[26] + D[95]*M[27] + D[96]*M[28] + D[98]*M[29] + D[99]*M[30] + D[100]*M[31] + D[101]*M[32] + D[102]*M[33];
#pragma omp atomic
L[12] += D[21]*M[0] + D[23]*M[1] + D[24]*M[2] + D[36]*M[3] + D[38]*M[4] + D[39]*M[5] + D[41]*M[6] + D[42]*M[7] + D[43]*M[8] + D[57]*M[9] + D[59]*M[10] + D[60]*M[11] + D[62]*M[12] + D[63]*M[13] + D[64]*M[14] + D[66]*M[15] + D[67]*M[16] + D[68]*M[17] + D[69]*M[18] + D[85]*M[19] + D[87]*M[20] + D[88]*M[21] + D[90]*M[22] + D[91]*M[23] + D[92]*M[24] + D[94]*M[25] + D[95]*M[26] + D[96]*M[27] + D[97]*M[28] + D[99]*M[29] + D[100]*M[30] + D[101]*M[31] + D[102]*M[32] + D[103]*M[33];
#pragma omp atomic
L[13] += D[22]*M[0] + D[25]*M[1] + D[26]*M[2] + D[37]*M[3] + D[40]*M[4] + D[41]*M[5] + D[44]*M[6] + D[45]*M[7] + D[46]*M[8] + D[58]*M[9] + D[61]*M[10] + D[62]*M[11] + D[65]*M[12] + D[66]*M[13] + D[67]*M[14] + D[70]*M[15] + D[71]*M[16] + D[72]*M[17] + D[73]*M[18] + D[86]*M[19] + D[89]*M[20] + D[90]*M[21] + D[93]*M[22] + D[94]*M[23] + D[95]*M[24] + D[98]*M[25] + D[99]*M[26] + D[100]*M[27] + D[101]*M[28] + D[104]*M[29] + D[105]*M[30] + D[106]*M[31] + D[107]*M[32] + D[108]*M[33];
#pragma omp atomic
L[14] += D[23]*M[0] + D[26]*M[1] + D[27]*M[2] + D[38]*M[3] + D[41]*M[4] + D[42]*M[5] + D[45]*M[6] + D[46]*M[7] + D[47]*M[8] + D[59]*M[9] + D[62]*M[10] + D[63]*M[11] + D[66]*M[12] + D[67]*M[13] + D[68]*M[14] + D[71]*M[15] + D[72]*M[16] + D[73]*M[17] + D[74]*M[18] + D[87]*M[19] + D[90]*M[20] + D[91]*M[21] + D[94]*M[22] + D[95]*M[23] + D[96]*M[24] + D[99]*M[25] + D[100]*M[26] + D[101]*M[27] + D[102]*M[28] + D[105]*M[29] + D[106]*M[30] + D[107]*M[31] + D[108]*M[32] + D[109]*M[33];
#pragma omp atomic
L[15] += D[24]*M[0] + D[27]*M[1] + D[28]*M[2] + D[39]*M[3] + D[42]*M[4] + D[43]*M[5] + D[46]*M[6] + D[47]*M[7] + D[48]*M[8] + D[60]*M[9] + D[63]*M[10] + D[64]*M[11] + D[67]*M[12] + D[68]*M[13] + D[69]*M[14] + D[72]*M[15] + D[73]*M[16] + D[74]*M[17] + D[75]*M[18] + D[88]*M[19] + D[91]*M[20] + D[92]*M[21] + D[95]*M[22] + D[96]*M[23] + D[97]*M[24] + D[100]*M[25] + D[101]*M[26] + D[102]*M[27] + D[103]*M[28] + D[106]*M[29] + D[107]*M[30] + D[108]*M[31] + D[109]*M[32] + D[110]*M[33];
#pragma omp atomic
L[16] += D[25]*M[0] + D[29]*M[1] + D[30]*M[2] + D[40]*M[3] + D[44]*M[4] + D[45]*M[5] + D[49]*M[6] + D[50]*M[7] + D[51]*M[8] + D[61]*M[9] + D[65]*M[10] + D[66]*M[11] + D[70]*M[12] + D[71]*M[13] + D[72]*M[14] + D[76]*M[15] + D[77]*M[16] + D[78]*M[17] + D[79]*M[18] + D[89]*M[19] + D[93]*M[20] + D[94]*M[21] + D[98]*M[22] + D[99]*M[23] + D[100]*M[24] + D[104]*M[25] + D[105]*M[26] + D[106]*M[27] + D[107]*M[28] + D[111]*M[29] + D[112]*M[30] + D[113]*M[31] + D[114]*M[32] + D[115]*M[33];
#pragma omp atomic
L[17] += D[26]*M[0] + D[30]*M[1] + D[31]*M[2] + D[41]*M[3] + D[45]*M[4] + D[46]*M[5] + D[50]*M[6] + D[51]*M[7] + D[52]*M[8] + D[62]*M[9] + D[66]*M[10] + D[67]*M[11] + D[71]*M[12] + D[72]*M[13] + D[73]*M[14] + D[77]*M[15] + D[78]*M[16] + D[79]*M[17] + D[80]*M[18] + D[90]*M[19] + D[94]*M[20] + D[95]*M[21] + D[99]*M[22] + D[100]*M[23] + D[101]*M[24] + D[105]*M[25] + D[106]*M[26] + D[107]*M[27] + D[108]*M[28] + D[112]*M[29] + D[113]*M[30] + D[114]*M[31] + D[115]*M[32] + D[116]*M[33];
#pragma omp atomic
L[18] += D[27]*M[0] + D[31]*M[1] + D[32]*M[2] + D[42]*M[3] + D[46]*M[4] + D[47]*M[5] + D[51]*M[6] + D[52]*M[7] + D[53]*M[8] + D[63]*M[9] + D[67]*M[10] + D[68]*M[11] + D[72]*M[12] + D[73]*M[13] + D[74]*M[14] + D[78]*M[15] + D[79]*M[16] + D[80]*M[17] + D[81]*M[18] + D[91]*M[19] + D[95]*M[20] + D[96]*M[21] + D[100]*M[22] + D[101]*M[23] + D[102]*M[24] + D[106]*M[25] + D[107]*M[26] + D[108]*M[27] + D[109]*M[28] + D[113]*M[29] + D[114]*M[30] + D[115]*M[31] + D[116]*M[32] + D[117]*M[33];
#pragma omp atomic
L[19] += D[28]*M[0] + D[32]*M[1] + D[33]*M[2] + D[43]*M[3] + D[47]*M[4] + D[48]*M[5] + D[52]*M[6] + D[53]*M[7] + D[54]*M[8] + D[64]*M[9] + D[68]*M[10] + D[69]*M[11] + D[73]*M[12] + D[74]*M[13] + D[75]*M[14] + D[79]*M[15] + D[80]*M[16] + D[81]*M[17] + D[82]*M[18] + D[92]*M[19] + D[96]*M[20] + D[97]*M[21] + D[101]*M[22] + D[102]*M[23] + D[103]*M[24] + D[107]*M[25] + D[108]*M[26] + D[109]*M[27] + D[110]*M[28] + D[114]*M[29] + D[115]*M[30] + D[116]*M[31] + D[117]*M[32] + D[118]*M[33];
#pragma omp atomic
L[20] += D[34]*M[0] + D[35]*M[1] + D[36]*M[2] + D[55]*M[3] + D[56]*M[4] + D[57]*M[5] + D[58]*M[6] + D[59]*M[7] + D[60]*M[8] + D[83]*M[9] + D[84]*M[10] + D[85]*M[11] + D[86]*M[12] + D[87]*M[13] + D[88]*M[14] + D[89]*M[15] + D[90]*M[16] + D[91]*M[17] + D[92]*M[18];
#pragma omp atomic
L[21] += D[35]*M[0] + D[37]*M[1] + D[38]*M[2] + D[56]*M[3] + D[58]*M[4] + D[59]*M[5] + D[61]*M[6] + D[62]*M[7] + D[63]*M[8] + D[84]*M[9] + D[86]*M[10] + D[87]*M[11] + D[89]*M[12] + D[90]*M[13] + D[91]*M[14] + D[93]*M[15] + D[94]*M[16] + D[95]*M[17] + D[96]*M[18];
#pragma omp atomic
L[22] += D[36]*M[0] + D[38]*M[1] + D[39]*M[2] + D[57]*M[3] + D[59]*M[4] + D[60]*M[5] + D[62]*M[6] + D[63]*M[7] + D[64]*M[8] + D[85]*M[9] + D[87]*M[10] + D[88]*M[11] + D[90]*M[12] + D[91]*M[13] + D[92]*M[14] + D[94]*M[15] + D[95]*M[16] + D[96]*M[17] + D[97]*M[18];
#pragma omp atomic
L[23] += D[37]*M[0] + D[40]*M[1] + D[41]*M[2] + D[58]*M[3] + D[61]*M[4] + D[62]*M[5] + D[65]*M[6] + D[66]*M[7] + D[67]*M[8] + D[86]*M[9] + D[89]*M[10] + D[90]*M[11] + D[93]*M[12] + D[94]*M[13] + D[95]*M[14] + D[98]*M[15] + D[99]*M[16] + D[100]*M[17] + D[101]*M[18];
#pragma omp atomic
L[24] += D[38]*M[0] + D[41]*M[1] + D[42]*M[2] + D[59]*M[3] + D[62]*M[4] + D[63]*M[5] + D[66]*M[6] + D[67]*M[7] + D[68]*M[8] + D[87]*M[9] + D[90]*M[10] + D[91]*M[11] + D[94]*M[12] + D[95]*M[13] + D[96]*M[14] + D[99]*M[15] + D[100]*M[16] + D[101]*M[17] + D[102]*M[18];
#pragma omp atomic
L[25] += D[39]*M[0] + D[42]*M[1] + D[43]*M[2] + D[60]*M[3] + D[63]*M[4] + D[64]*M[5] + D[67]*M[6] + D[68]*M[7] + D[69]*M[8] + D[88]*M[9] + D[91]*M[10] + D[92]*M[11] + D[95]*M[12] + D[96]*M[13] + D[97]*M[14] + D[100]*M[15] + D[101]*M[16] + D[102]*M[17] + D[103]*M[18];
#pragma omp atomic
L[26] += D[40]*M[0] + D[44]*M[1] + D[45]*M[2] + D[61]*M[3] + D[65]*M[4] + D[66]*M[5] + D[70]*M[6] + D[71]*M[7] + D[72]*M[8] + D[89]*M[9] + D[93]*M[10] + D[94]*M[11] + D[98]*M[12] + D[99]*M[13] + D[100]*M[14] + D[104]*M[15] + D[105]*M[16] + D[106]*M[17] + D[107]*M[18];
#pragma omp atomic
L[27] += D[41]*M[0] + D[45]*M[1] + D[46]*M[2] + D[62]*M[3] + D[66]*M[4] + D[67]*M[5] + D[71]*M[6] + D[72]*M[7] + D[73]*M[8] + D[90]*M[9] + D[94]*M[10] + D[95]*M[11] + D[99]*M[12] + D[100]*M[13] + D[101]*M[14] + D[105]*M[15] + D[106]*M[16] + D[107]*M[17] + D[108]*M[18];
#pragma omp atomic
L[28] += D[42]*M[0] + D[46]*M[1] + D[47]*M[2] + D[63]*M[3] + D[67]*M[4] + D[68]*M[5] + D[72]*M[6] + D[73]*M[7] + D[74]*M[8] + D[91]*M[9] + D[95]*M[10] + D[96]*M[11] + D[100]*M[12] + D[101]*M[13] + D[102]*M[14] + D[106]*M[15] + D[107]*M[16] + D[108]*M[17] + D[109]*M[18];
#pragma omp atomic
L[29] += D[43]*M[0] + D[47]*M[1] + D[48]*M[2] + D[64]*M[3] + D[68]*M[4] + D[69]*M[5] + D[73]*M[6] + D[74]*M[7] + D[75]*M[8] + D[92]*M[9] + D[96]*M[10] + D[97]*M[11] + D[101]*M[12] + D[102]*M[13] + D[103]*M[14] + D[107]*M[15] + D[108]*M[16] + D[109]*M[17] + D[110]*M[18];
#pragma omp atomic
L[30] += D[44]*M[0] + D[49]*M[1] + D[50]*M[2] + D[65]*M[3] + D[70]*M[4] + D[71]*M[5] + D[76]*M[6] + D[77]*M[7] + D[78]*M[8] + D[93]*M[9] + D[98]*M[10] + D[99]*M[11] + D[104]*M[12] + D[105]*M[13] + D[106]*M[14] + D[111]*M[15] + D[112]*M[16] + D[113]*M[17] + D[114]*M[18];
#pragma omp atomic
L[31] += D[45]*M[0] + D[50]*M[1] + D[51]*M[2] + D[66]*M[3] + D[71]*M[4] + D[72]*M[5] + D[77]*M[6] + D[78]*M[7] + D[79]*M[8] + D[94]*M[9] + D[99]*M[10] + D[100]*M[11] + D[105]*M[12] + D[106]*M[13] + D[107]*M[14] + D[112]*M[15] + D[113]*M[16] + D[114]*M[17] + D[115]*M[18];
#pragma omp atomic
L[32] += D[46]*M[0] + D[51]*M[1] + D[52]*M[2] + D[67]*M[3] + D[72]*M[4] + D[73]*M[5] + D[78]*M[6] + D[79]*M[7] + D[80]*M[8] + D[95]*M[9] + D[100]*M[10] + D[101]*M[11] + D[106]*M[12] + D[107]*M[13] + D[108]*M[14] + D[113]*M[15] + D[114]*M[16] + D[115]*M[17] + D[116]*M[18];
#pragma omp atomic
L[33] += D[47]*M[0] + D[52]*M[1] + D[53]*M[2] + D[68]*M[3] + D[73]*M[4] + D[74]*M[5] + D[79]*M[6] + D[80]*M[7] + D[81]*M[8] + D[96]*M[9] + D[101]*M[10] + D[102]*M[11] + D[107]*M[12] + D[108]*M[13] + D[109]*M[14] + D[114]*M[15] + D[115]*M[16] + D[116]*M[17] + D[117]*M[18];
#pragma omp atomic
L[34] += D[48]*M[0] + D[53]*M[1] + D[54]*M[2] + D[69]*M[3] + D[74]*M[4] + D[75]*M[5] + D[80]*M[6] + D[81]*M[7] + D[82]*M[8] + D[97]*M[9] + D[102]*M[10] + D[103]*M[11] + D[108]*M[12] + D[109]*M[13] + D[110]*M[14] + D[115]*M[15] + D[116]*M[16] + D[117]*M[17] + D[118]*M[18];
#pragma omp atomic
L[35] += D[55]*M[0] + D[56]*M[1] + D[57]*M[2] + D[83]*M[3] + D[84]*M[4] + D[85]*M[5] + D[86]*M[6] + D[87]*M[7] + D[88]*M[8];
#pragma omp atomic
L[36] += D[56]*M[0] + D[58]*M[1] + D[59]*M[2] + D[84]*M[3] + D[86]*M[4] + D[87]*M[5] + D[89]*M[6] + D[90]*M[7] + D[91]*M[8];
#pragma omp atomic
L[37] += D[57]*M[0] + D[59]*M[1] + D[60]*M[2] + D[85]*M[3] + D[87]*M[4] + D[88]*M[5] + D[90]*M[6] + D[91]*M[7] + D[92]*M[8];
#pragma omp atomic
L[38] += D[58]*M[0] + D[61]*M[1] + D[62]*M[2] + D[86]*M[3] + D[89]*M[4] + D[90]*M[5] + D[93]*M[6] + D[94]*M[7] + D[95]*M[8];
#pragma omp atomic
L[39] += D[59]*M[0] + D[62]*M[1] + D[63]*M[2] + D[87]*M[3] + D[90]*M[4] + D[91]*M[5] + D[94]*M[6] + D[95]*M[7] + D[96]*M[8];
#pragma omp atomic
L[40] += D[60]*M[0] + D[63]*M[1] + D[64]*M[2] + D[88]*M[3] + D[91]*M[4] + D[92]*M[5] + D[95]*M[6] + D[96]*M[7] + D[97]*M[8];
#pragma omp atomic
L[41] += D[61]*M[0] + D[65]*M[1] + D[66]*M[2] + D[89]*M[3] + D[93]*M[4] + D[94]*M[5] + D[98]*M[6] + D[99]*M[7] + D[100]*M[8];
#pragma omp atomic
L[42] += D[62]*M[0] + D[66]*M[1] + D[67]*M[2] + D[90]*M[3] + D[94]*M[4] + D[95]*M[5] + D[99]*M[6] + D[100]*M[7] + D[101]*M[8];
#pragma omp atomic
L[43] += D[63]*M[0] + D[67]*M[1] + D[68]*M[2] + D[91]*M[3] + D[95]*M[4] + D[96]*M[5] + D[100]*M[6] + D[101]*M[7] + D[102]*M[8];
#pragma omp atomic
L[44] += D[64]*M[0] + D[68]*M[1] + D[69]*M[2] + D[92]*M[3] + D[96]*M[4] + D[97]*M[5] + D[101]*M[6] + D[102]*M[7] + D[103]*M[8];
#pragma omp atomic
L[45] += D[65]*M[0] + D[70]*M[1] + D[71]*M[2] + D[93]*M[3] + D[98]*M[4] + D[99]*M[5] + D[104]*M[6] + D[105]*M[7] + D[106]*M[8];
#pragma omp atomic
L[46] += D[66]*M[0] + D[71]*M[1] + D[72]*M[2] + D[94]*M[3] + D[99]*M[4] + D[100]*M[5] + D[105]*M[6] + D[106]*M[7] + D[107]*M[8];
#pragma omp atomic
L[47] += D[67]*M[0] + D[72]*M[1] + D[73]*M[2] + D[95]*M[3] + D[100]*M[4] + D[101]*M[5] + D[106]*M[6] + D[107]*M[7] + D[108]*M[8];
#pragma omp atomic
L[48] += D[68]*M[0] + D[73]*M[1] + D[74]*M[2] + D[96]*M[3] + D[101]*M[4] + D[102]*M[5] + D[107]*M[6] + D[108]*M[7] + D[109]*M[8];
#pragma omp atomic
L[49] += D[69]*M[0] + D[74]*M[1] + D[75]*M[2] + D[97]*M[3] + D[102]*M[4] + D[103]*M[5] + D[108]*M[6] + D[109]*M[7] + D[110]*M[8];
#pragma omp atomic
L[50] += D[70]*M[0] + D[76]*M[1] + D[77]*M[2] + D[98]*M[3] + D[104]*M[4] + D[105]*M[5] + D[111]*M[6] + D[112]*M[7] + D[113]*M[8];
#pragma omp atomic
L[51] += D[71]*M[0] + D[77]*M[1] + D[78]*M[2] + D[99]*M[3] + D[105]*M[4] + D[106]*M[5] + D[112]*M[6] + D[113]*M[7] + D[114]*M[8];
#pragma omp atomic
L[52] += D[72]*M[0] + D[78]*M[1] + D[79]*M[2] + D[100]*M[3] + D[106]*M[4] + D[107]*M[5] + D[113]*M[6] + D[114]*M[7] + D[115]*M[8];
#pragma omp atomic
L[53] += D[73]*M[0] + D[79]*M[1] + D[80]*M[2] + D[101]*M[3] + D[107]*M[4] + D[108]*M[5] + D[114]*M[6] + D[115]*M[7] + D[116]*M[8];
#pragma omp atomic
L[54] += D[74]*M[0] + D[80]*M[1] + D[81]*M[2] + D[102]*M[3] + D[108]*M[4] + D[109]*M[5] + D[115]*M[6] + D[116]*M[7] + D[117]*M[8];
#pragma omp atomic
L[55] += D[75]*M[0] + D[81]*M[1] + D[82]*M[2] + D[103]*M[3] + D[109]*M[4] + D[110]*M[5] + D[116]*M[6] + D[117]*M[7] + D[118]*M[8];
#pragma omp atomic
L[56] += D[83]*M[0] + D[84]*M[1] + D[85]*M[2];
#pragma omp atomic
L[57] += D[84]*M[0] + D[86]*M[1] + D[87]*M[2];
#pragma omp atomic
L[58] += D[85]*M[0] + D[87]*M[1] + D[88]*M[2];
#pragma omp atomic
L[59] += D[86]*M[0] + D[89]*M[1] + D[90]*M[2];
#pragma omp atomic
L[60] += D[87]*M[0] + D[90]*M[1] + D[91]*M[2];
#pragma omp atomic
L[61] += D[88]*M[0] + D[91]*M[1] + D[92]*M[2];
#pragma omp atomic
L[62] += D[89]*M[0] + D[93]*M[1] + D[94]*M[2];
#pragma omp atomic
L[63] += D[90]*M[0] + D[94]*M[1] + D[95]*M[2];
#pragma omp atomic
L[64] += D[91]*M[0] + D[95]*M[1] + D[96]*M[2];
#pragma omp atomic
L[65] += D[92]*M[0] + D[96]*M[1] + D[97]*M[2];
#pragma omp atomic
L[66] += D[93]*M[0] + D[98]*M[1] + D[99]*M[2];
#pragma omp atomic
L[67] += D[94]*M[0] + D[99]*M[1] + D[100]*M[2];
#pragma omp atomic
L[68] += D[95]*M[0] + D[100]*M[1] + D[101]*M[2];
#pragma omp atomic
L[69] += D[96]*M[0] + D[101]*M[1] + D[102]*M[2];
#pragma omp atomic
L[70] += D[97]*M[0] + D[102]*M[1] + D[103]*M[2];
#pragma omp atomic
L[71] += D[98]*M[0] + D[104]*M[1] + D[105]*M[2];
#pragma omp atomic
L[72] += D[99]*M[0] + D[105]*M[1] + D[106]*M[2];
#pragma omp atomic
L[73] += D[100]*M[0] + D[106]*M[1] + D[107]*M[2];
#pragma omp atomic
L[74] += D[101]*M[0] + D[107]*M[1] + D[108]*M[2];
#pragma omp atomic
L[75] += D[102]*M[0] + D[108]*M[1] + D[109]*M[2];
#pragma omp atomic
L[76] += D[103]*M[0] + D[109]*M[1] + D[110]*M[2];
#pragma omp atomic
L[77] += D[104]*M[0] + D[111]*M[1] + D[112]*M[2];
#pragma omp atomic
L[78] += D[105]*M[0] + D[112]*M[1] + D[113]*M[2];
#pragma omp atomic
L[79] += D[106]*M[0] + D[113]*M[1] + D[114]*M[2];
#pragma omp atomic
L[80] += D[107]*M[0] + D[114]*M[1] + D[115]*M[2];
#pragma omp atomic
L[81] += D[108]*M[0] + D[115]*M[1] + D[116]*M[2];
#pragma omp atomic
L[82] += D[109]*M[0] + D[116]*M[1] + D[117]*M[2];
#pragma omp atomic
L[83] += D[110]*M[0] + D[117]*M[1] + D[118]*M[2];

}

void L2L_7(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[14];
double Lstmp4 = Lstmp3*y;
double Lstmp5 = pow(x, 2);
double Lstmp6 = (1.0/2.0)*Lstmp5;
double Lstmp7 = pow(x, 3);
double Lstmp8 = (1.0/6.0)*Lstmp7;
double Lstmp9 = pow(x, 4);
double Lstmp10 = (1.0/24.0)*Lstmp9;
double Lstmp11 = (1.0/120.0)*pow(x, 5);
double Lstmp12 = pow(y, 2);
double Lstmp13 = (1.0/2.0)*Lstmp12;
double Lstmp14 = pow(y, 3);
double Lstmp15 = (1.0/6.0)*Lstmp14;
double Lstmp16 = pow(y, 4);
double Lstmp17 = (1.0/24.0)*Lstmp16;
double Lstmp18 = (1.0/120.0)*pow(y, 5);
double Lstmp19 = pow(z, 2);
double Lstmp20 = (1.0/2.0)*Lstmp19;
double Lstmp21 = pow(z, 3);
double Lstmp22 = (1.0/6.0)*Lstmp21;
double Lstmp23 = pow(z, 4);
double Lstmp24 = (1.0/24.0)*Lstmp23;
double Lstmp25 = (1.0/120.0)*pow(z, 5);
double Lstmp26 = x*L[13];
double Lstmp27 = x*L[26];
double Lstmp28 = x*L[45];
double Lstmp29 = x*L[71];
double Lstmp30 = x*L[15];
double Lstmp31 = x*L[29];
double Lstmp32 = x*L[49];
double Lstmp33 = x*L[76];
double Lstmp34 = y*L[11];
double Lstmp35 = z*L[12];
double Lstmp36 = y*L[21];
double Lstmp37 = z*L[22];
double Lstmp38 = y*L[36];
double Lstmp39 = z*L[37];
double Lstmp40 = y*L[57];
double Lstmp41 = z*L[58];
double Lstmp42 = y*L[18];
double Lstmp43 = y*L[33];
double Lstmp44 = y*L[54];
double Lstmp45 = y*L[82];
double Lstmp46 = z*L[17];
double Lstmp47 = z*L[31];
double Lstmp48 = z*L[51];
double Lstmp49 = z*L[78];
double Lstmp50 = y*L[28];
double Lstmp51 = Lstmp50*x;
double Lstmp52 = y*L[48];
double Lstmp53 = Lstmp52*x;
double Lstmp54 = y*L[75];
double Lstmp55 = Lstmp54*x;
double Lstmp56 = z*L[27];
double Lstmp57 = Lstmp56*x;
double Lstmp58 = z*L[46];
double Lstmp59 = Lstmp58*x;
double Lstmp60 = z*L[72];
double Lstmp61 = Lstmp60*x;
double Lstmp62 = z*L[24];
double Lstmp63 = Lstmp62*y;
double Lstmp64 = z*L[39];
double Lstmp65 = Lstmp64*y;
double Lstmp66 = z*L[60];
double Lstmp67 = Lstmp66*y;
double Lstmp68 = (1.0/4.0)*Lstmp5;
double Lstmp69 = Lstmp12*Lstmp68;
double Lstmp70 = (1.0/12.0)*Lstmp5;
double Lstmp71 = Lstmp14*Lstmp70;
double Lstmp72 = (1.0/48.0)*Lstmp5;
double Lstmp73 = Lstmp19*Lstmp68;
double Lstmp74 = Lstmp21*Lstmp70;
double Lstmp75 = (1.0/12.0)*Lstmp7;
double Lstmp76 = Lstmp12*Lstmp75;
double Lstmp77 = (1.0/36.0)*Lstmp7;
double Lstmp78 = Lstmp19*Lstmp75;
double Lstmp79 = (1.0/48.0)*Lstmp9;
double Lstmp80 = Lstmp12*Lstmp19;
double Lstmp81 = (1.0/4.0)*Lstmp80;
double Lstmp82 = (1.0/12.0)*Lstmp12*Lstmp21;
double Lstmp83 = (1.0/12.0)*Lstmp14*Lstmp19;
double Lstmp84 = x*L[47];
double Lstmp85 = x*L[74];
double Lstmp86 = x*L[73];
double Lstmp87 = y*L[43];
double Lstmp88 = y*L[69];
double Lstmp89 = z*L[42];
double Lstmp90 = z*L[67];
double Lstmp91 = y*L[64];
double Lstmp92 = z*L[63];
double Lstmp93 = x*L[23];
double Lstmp94 = x*L[41];
double Lstmp95 = x*L[66];
double Lstmp96 = x*L[25];
double Lstmp97 = x*L[44];
double Lstmp98 = x*L[70];
double Lstmp99 = Lstmp87*x;
double Lstmp100 = Lstmp88*x;
double Lstmp101 = Lstmp89*x;
double Lstmp102 = Lstmp90*x;
double Lstmp103 = x*L[68];
double Lstmp104 = y*L[13];
double Lstmp105 = Lstmp56*y;
double Lstmp106 = x*L[28];
double Lstmp107 = x*L[48];
double Lstmp108 = x*L[75];
double Lstmp109 = y*L[23];
double Lstmp110 = y*L[38];
double Lstmp111 = y*L[59];
double Lstmp112 = y*L[32];
double Lstmp113 = y*L[53];
double Lstmp114 = y*L[81];
double Lstmp115 = y*L[47];
double Lstmp116 = Lstmp115*x;
double Lstmp117 = y*L[74];
double Lstmp118 = Lstmp117*x;
double Lstmp119 = Lstmp89*y;
double Lstmp120 = Lstmp92*y;
double Lstmp121 = y*L[68];
double Lstmp122 = y*L[14];
double Lstmp123 = z*L[15];
double Lstmp124 = z*L[18];
double Lstmp125 = z*L[28];
double Lstmp126 = Lstmp125*y;
double Lstmp127 = x*L[27];
double Lstmp128 = x*L[46];
double Lstmp129 = x*L[72];
double Lstmp130 = y*L[24];
double Lstmp131 = z*L[25];
double Lstmp132 = y*L[39];
double Lstmp133 = z*L[40];
double Lstmp134 = y*L[60];
double Lstmp135 = z*L[61];
double Lstmp136 = z*L[32];
double Lstmp137 = z*L[52];
double Lstmp138 = z*L[79];
double Lstmp139 = z*L[47];
double Lstmp140 = Lstmp139*x;
double Lstmp141 = z*L[73];
double Lstmp142 = Lstmp141*x;
double Lstmp143 = z*L[43];
double Lstmp144 = Lstmp143*y;
double Lstmp145 = z*L[64];
double Lstmp146 = Lstmp145*y;
double Lstmp147 = z*L[68];
double Lstmp148 = x*L[38];
double Lstmp149 = x*L[62];
double Lstmp150 = x*L[40];
double Lstmp151 = x*L[65];
double Lstmp152 = Lstmp91*x;
double Lstmp153 = Lstmp92*x;
double Lstmp154 = x*L[43];
double Lstmp155 = x*L[69];
double Lstmp156 = Lstmp121*x;
double Lstmp157 = x*L[42];
double Lstmp158 = x*L[67];
double Lstmp159 = Lstmp147*x;
double Lstmp160 = y*L[26];
double Lstmp161 = Lstmp58*y;
double Lstmp162 = y*L[41];
double Lstmp163 = y*L[62];
double Lstmp164 = y*L[52];
double Lstmp165 = y*L[80];
double Lstmp166 = y*L[73];
double Lstmp167 = Lstmp166*x;
double Lstmp168 = Lstmp90*y;
double Lstmp169 = y*L[27];
double Lstmp170 = Lstmp139*y;
double Lstmp171 = y*L[42];
double Lstmp172 = y*L[63];
double Lstmp173 = Lstmp147*y;
double Lstmp174 = z*L[29];
double Lstmp175 = z*L[33];
double Lstmp176 = z*L[48];
double Lstmp177 = Lstmp176*y;
double Lstmp178 = z*L[44];
double Lstmp179 = z*L[65];
double Lstmp180 = z*L[53];
double Lstmp181 = z*L[80];
double Lstmp182 = z*L[74];
double Lstmp183 = Lstmp182*x;
double Lstmp184 = z*L[69];
double Lstmp185 = Lstmp184*y;
double Lstmp186 = x*L[59];
double Lstmp187 = x*L[61];
double Lstmp188 = x*L[64];
double Lstmp189 = x*L[63];
double Lstmp190 = y*L[45];
double Lstmp191 = Lstmp60*y;
double Lstmp192 = y*L[66];
double Lstmp193 = y*L[79];
double Lstmp194 = y*L[46];
double Lstmp195 = Lstmp141*y;
double Lstmp196 = y*L[67];
double Lstmp197 = Lstmp182*y;
double Lstmp198 = z*L[49];
double Lstmp199 = z*L[54];
double Lstmp200 = z*L[75];
double Lstmp201 = Lstmp200*y;
double Lstmp202 = z*L[70];
double Lstmp203 = z*L[81];
double Lstmp204 = y*L[71];
double Lstmp205 = y*L[72];
double Lstmp206 = z*L[76];
double Lstmp207 = z*L[82];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp10*Lstmp38 + Lstmp10*Lstmp39 + Lstmp10*Lstmp67 + Lstmp10*L[20] + Lstmp11*Lstmp40 + Lstmp11*Lstmp41 + Lstmp11*L[35] + (1.0/48.0)*Lstmp12*Lstmp23*L[81] + Lstmp12*Lstmp79*L[59] + Lstmp13*Lstmp26 + Lstmp13*Lstmp46 + Lstmp13*Lstmp57 + Lstmp13*L[7] + (1.0/36.0)*Lstmp14*Lstmp21*L[80] + Lstmp14*Lstmp77*L[62] + Lstmp15*Lstmp27 + Lstmp15*Lstmp47 + Lstmp15*Lstmp59 + Lstmp15*L[16] + (1.0/48.0)*Lstmp16*Lstmp19*L[79] + Lstmp16*Lstmp72*L[66] + Lstmp17*Lstmp28 + Lstmp17*Lstmp48 + Lstmp17*Lstmp61 + Lstmp17*L[30] + Lstmp18*Lstmp29 + Lstmp18*Lstmp49 + Lstmp18*L[50] + Lstmp19*Lstmp79*L[61] + Lstmp2*y + Lstmp20*Lstmp30 + Lstmp20*Lstmp42 + Lstmp20*Lstmp51 + Lstmp20*L[9] + Lstmp21*Lstmp77*L[65] + Lstmp22*Lstmp31 + Lstmp22*Lstmp43 + Lstmp22*Lstmp53 + Lstmp22*L[19] + Lstmp23*Lstmp72*L[70] + Lstmp24*Lstmp32 + Lstmp24*Lstmp44 + Lstmp24*Lstmp55 + Lstmp24*L[34] + Lstmp25*Lstmp33 + Lstmp25*Lstmp45 + Lstmp25*L[55] + Lstmp34*Lstmp6 + Lstmp35*Lstmp6 + Lstmp36*Lstmp8 + Lstmp37*Lstmp8 + Lstmp4*x + (1.0/8.0)*Lstmp5*Lstmp80*L[68] + Lstmp6*Lstmp63 + Lstmp6*L[4] + Lstmp65*Lstmp8 + Lstmp69*Lstmp89 + Lstmp69*L[23] + Lstmp71*Lstmp90 + Lstmp71*L[41] + Lstmp73*Lstmp87 + Lstmp73*L[25] + Lstmp74*Lstmp88 + Lstmp74*L[44] + Lstmp76*Lstmp92 + Lstmp76*L[38] + Lstmp78*Lstmp91 + Lstmp78*L[40] + Lstmp8*L[10] + Lstmp81*Lstmp84 + Lstmp81*L[32] + Lstmp82*Lstmp85 + Lstmp82*L[53] + Lstmp83*Lstmp86 + Lstmp83*L[52] + (1.0/720.0)*pow(x, 6)*L[56] + x*L[1] + (1.0/720.0)*pow(y, 6)*L[77] + y*L[2] + (1.0/720.0)*pow(z, 6)*L[83] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + Lstmp10*Lstmp40 + Lstmp10*Lstmp41 + Lstmp10*L[35] + Lstmp100*Lstmp22 + Lstmp101*Lstmp13 + Lstmp102*Lstmp15 + Lstmp103*Lstmp81 + Lstmp11*L[56] + Lstmp13*Lstmp56 + Lstmp13*Lstmp93 + Lstmp13*L[13] + Lstmp15*Lstmp58 + Lstmp15*Lstmp94 + Lstmp15*L[26] + Lstmp17*Lstmp60 + Lstmp17*Lstmp95 + Lstmp17*L[45] + Lstmp18*L[71] + Lstmp20*Lstmp50 + Lstmp20*Lstmp96 + Lstmp20*Lstmp99 + Lstmp20*L[15] + Lstmp22*Lstmp52 + Lstmp22*Lstmp97 + Lstmp22*L[29] + Lstmp24*Lstmp54 + Lstmp24*Lstmp98 + Lstmp24*L[49] + Lstmp25*L[76] + Lstmp34*x + Lstmp35*x + Lstmp36*Lstmp6 + Lstmp37*Lstmp6 + Lstmp38*Lstmp8 + Lstmp39*Lstmp8 + Lstmp4 + Lstmp6*Lstmp65 + Lstmp6*L[10] + Lstmp63*x + Lstmp67*Lstmp8 + Lstmp69*Lstmp92 + Lstmp69*L[38] + Lstmp71*L[62] + Lstmp73*Lstmp91 + Lstmp73*L[40] + Lstmp74*L[65] + Lstmp76*L[59] + Lstmp78*L[61] + Lstmp8*L[20] + Lstmp81*L[47] + Lstmp82*L[74] + Lstmp83*L[73] + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp10*Lstmp111 + Lstmp10*Lstmp66 + Lstmp10*L[36] + Lstmp104*x + Lstmp105*x + Lstmp106*Lstmp20 + Lstmp107*Lstmp22 + Lstmp108*Lstmp24 + Lstmp109*Lstmp6 + Lstmp11*L[57] + Lstmp110*Lstmp8 + Lstmp112*Lstmp20 + Lstmp113*Lstmp22 + Lstmp114*Lstmp24 + Lstmp116*Lstmp20 + Lstmp118*Lstmp22 + Lstmp119*Lstmp6 + Lstmp120*Lstmp8 + Lstmp121*Lstmp73 + Lstmp13*Lstmp27 + Lstmp13*Lstmp47 + Lstmp13*Lstmp59 + Lstmp13*L[16] + Lstmp15*Lstmp28 + Lstmp15*Lstmp48 + Lstmp15*Lstmp61 + Lstmp15*L[30] + Lstmp17*Lstmp29 + Lstmp17*Lstmp49 + Lstmp17*L[50] + Lstmp18*L[77] + Lstmp2 + Lstmp20*L[18] + Lstmp22*L[33] + Lstmp24*L[54] + Lstmp25*L[82] + Lstmp3*x + Lstmp46*y + Lstmp6*Lstmp62 + Lstmp6*L[11] + Lstmp64*Lstmp8 + Lstmp69*Lstmp90 + Lstmp69*L[41] + Lstmp71*L[66] + Lstmp73*L[43] + Lstmp74*L[69] + Lstmp76*L[62] + Lstmp78*L[64] + Lstmp8*L[21] + Lstmp81*Lstmp86 + Lstmp81*L[52] + Lstmp82*L[80] + Lstmp83*L[79] + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += Lstmp10*Lstmp134 + Lstmp10*Lstmp135 + Lstmp10*L[37] + Lstmp11*L[58] + Lstmp122*x + Lstmp123*x + Lstmp124*y + Lstmp126*x + Lstmp127*Lstmp13 + Lstmp128*Lstmp15 + Lstmp129*Lstmp17 + Lstmp13*Lstmp136 + Lstmp13*Lstmp140 + Lstmp13*L[17] + Lstmp130*Lstmp6 + Lstmp131*Lstmp6 + Lstmp132*Lstmp8 + Lstmp133*Lstmp8 + Lstmp137*Lstmp15 + Lstmp138*Lstmp17 + Lstmp142*Lstmp15 + Lstmp144*Lstmp6 + Lstmp146*Lstmp8 + Lstmp147*Lstmp69 + Lstmp15*L[31] + Lstmp17*L[51] + Lstmp18*L[78] + Lstmp20*Lstmp31 + Lstmp20*Lstmp43 + Lstmp20*Lstmp53 + Lstmp20*L[19] + Lstmp22*Lstmp32 + Lstmp22*Lstmp44 + Lstmp22*Lstmp55 + Lstmp22*L[34] + Lstmp24*Lstmp33 + Lstmp24*Lstmp45 + Lstmp24*L[55] + Lstmp25*L[83] + Lstmp6*L[12] + Lstmp69*L[42] + Lstmp71*L[67] + Lstmp73*Lstmp88 + Lstmp73*L[44] + Lstmp74*L[70] + Lstmp76*L[63] + Lstmp78*L[65] + Lstmp8*L[22] + Lstmp81*Lstmp85 + Lstmp81*L[53] + Lstmp82*L[81] + Lstmp83*L[80] + x*L[6] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += Lstmp10*L[56] + Lstmp13*Lstmp148 + Lstmp13*Lstmp153 + Lstmp13*Lstmp89 + Lstmp13*L[23] + Lstmp149*Lstmp15 + Lstmp15*Lstmp90 + Lstmp15*L[41] + Lstmp150*Lstmp20 + Lstmp151*Lstmp22 + Lstmp152*Lstmp20 + Lstmp17*L[66] + Lstmp20*Lstmp87 + Lstmp20*L[25] + Lstmp22*Lstmp88 + Lstmp22*L[44] + Lstmp24*L[70] + Lstmp34 + Lstmp35 + Lstmp36*x + Lstmp37*x + Lstmp38*Lstmp6 + Lstmp39*Lstmp6 + Lstmp40*Lstmp8 + Lstmp41*Lstmp8 + Lstmp6*Lstmp67 + Lstmp6*L[20] + Lstmp63 + Lstmp65*x + Lstmp69*L[59] + Lstmp73*L[61] + Lstmp8*L[35] + Lstmp81*L[68] + x*L[10] + L[4];
#pragma omp atomic
Ls[5] += Lstmp10*L[57] + Lstmp102*Lstmp13 + Lstmp104 + Lstmp105 + Lstmp109*x + Lstmp110*Lstmp6 + Lstmp111*Lstmp8 + Lstmp115*Lstmp20 + Lstmp117*Lstmp22 + Lstmp119*x + Lstmp120*Lstmp6 + Lstmp13*Lstmp58 + Lstmp13*Lstmp94 + Lstmp13*L[26] + Lstmp15*Lstmp60 + Lstmp15*Lstmp95 + Lstmp15*L[45] + Lstmp154*Lstmp20 + Lstmp155*Lstmp22 + Lstmp156*Lstmp20 + Lstmp17*L[71] + Lstmp20*L[28] + Lstmp22*L[48] + Lstmp24*L[75] + Lstmp3 + Lstmp6*Lstmp64 + Lstmp6*L[21] + Lstmp62*x + Lstmp66*Lstmp8 + Lstmp69*L[62] + Lstmp73*L[64] + Lstmp8*L[36] + Lstmp81*L[73] + x*L[11] + L[5];
#pragma omp atomic
Ls[6] += Lstmp10*L[58] + Lstmp100*Lstmp20 + Lstmp122 + Lstmp123 + Lstmp126 + Lstmp13*Lstmp139 + Lstmp13*Lstmp157 + Lstmp13*Lstmp159 + Lstmp13*L[27] + Lstmp130*x + Lstmp131*x + Lstmp132*Lstmp6 + Lstmp133*Lstmp6 + Lstmp134*Lstmp8 + Lstmp135*Lstmp8 + Lstmp141*Lstmp15 + Lstmp144*x + Lstmp146*Lstmp6 + Lstmp15*Lstmp158 + Lstmp15*L[46] + Lstmp17*L[72] + Lstmp20*Lstmp52 + Lstmp20*Lstmp97 + Lstmp20*L[29] + Lstmp22*Lstmp54 + Lstmp22*Lstmp98 + Lstmp22*L[49] + Lstmp24*L[76] + Lstmp6*L[22] + Lstmp69*L[63] + Lstmp73*L[65] + Lstmp8*L[37] + Lstmp81*L[74] + x*L[12] + L[6];
#pragma omp atomic
Ls[7] += Lstmp10*L[59] + Lstmp13*Lstmp28 + Lstmp13*Lstmp48 + Lstmp13*Lstmp61 + Lstmp13*L[30] + Lstmp15*Lstmp29 + Lstmp15*Lstmp49 + Lstmp15*L[50] + Lstmp160*x + Lstmp161*x + Lstmp162*Lstmp6 + Lstmp163*Lstmp8 + Lstmp164*Lstmp20 + Lstmp165*Lstmp22 + Lstmp167*Lstmp20 + Lstmp168*Lstmp6 + Lstmp17*L[77] + Lstmp20*Lstmp84 + Lstmp20*L[32] + Lstmp22*Lstmp85 + Lstmp22*L[53] + Lstmp24*L[81] + Lstmp26 + Lstmp46 + Lstmp47*y + Lstmp57 + Lstmp6*Lstmp89 + Lstmp6*L[23] + Lstmp69*L[66] + Lstmp73*L[68] + Lstmp8*Lstmp92 + Lstmp8*L[38] + Lstmp81*L[79] + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += Lstmp10*L[60] + Lstmp107*Lstmp20 + Lstmp108*Lstmp22 + Lstmp113*Lstmp20 + Lstmp114*Lstmp22 + Lstmp118*Lstmp20 + Lstmp124 + Lstmp125*x + Lstmp128*Lstmp13 + Lstmp129*Lstmp15 + Lstmp13*Lstmp137 + Lstmp13*Lstmp142 + Lstmp13*L[31] + Lstmp136*y + Lstmp138*Lstmp15 + Lstmp143*Lstmp6 + Lstmp145*Lstmp8 + Lstmp15*L[51] + Lstmp169*x + Lstmp17*L[78] + Lstmp170*x + Lstmp171*Lstmp6 + Lstmp172*Lstmp8 + Lstmp173*Lstmp6 + Lstmp20*L[33] + Lstmp22*L[54] + Lstmp24*L[82] + Lstmp6*L[24] + Lstmp69*L[67] + Lstmp73*L[69] + Lstmp8*L[39] + Lstmp81*L[80] + x*L[14] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += Lstmp10*L[61] + Lstmp13*Lstmp180 + Lstmp13*Lstmp183 + Lstmp13*Lstmp84 + Lstmp13*L[32] + Lstmp15*Lstmp181 + Lstmp15*Lstmp86 + Lstmp15*L[52] + Lstmp17*L[79] + Lstmp174*x + Lstmp175*y + Lstmp177*x + Lstmp178*Lstmp6 + Lstmp179*Lstmp8 + Lstmp185*Lstmp6 + Lstmp20*Lstmp32 + Lstmp20*Lstmp44 + Lstmp20*Lstmp55 + Lstmp20*L[34] + Lstmp22*Lstmp33 + Lstmp22*Lstmp45 + Lstmp22*L[55] + Lstmp24*L[83] + Lstmp30 + Lstmp42 + Lstmp51 + Lstmp6*Lstmp87 + Lstmp6*L[25] + Lstmp69*L[68] + Lstmp73*L[70] + Lstmp8*Lstmp91 + Lstmp8*L[40] + Lstmp81*L[81] + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += Lstmp13*Lstmp186 + Lstmp13*Lstmp92 + Lstmp13*L[38] + Lstmp15*L[62] + Lstmp187*Lstmp20 + Lstmp20*Lstmp91 + Lstmp20*L[40] + Lstmp22*L[65] + Lstmp36 + Lstmp37 + Lstmp38*x + Lstmp39*x + Lstmp40*Lstmp6 + Lstmp41*Lstmp6 + Lstmp6*L[35] + Lstmp65 + Lstmp67*x + Lstmp8*L[56] + x*L[20] + L[10];
#pragma omp atomic
Ls[11] += Lstmp109 + Lstmp110*x + Lstmp111*Lstmp6 + Lstmp119 + Lstmp120*x + Lstmp121*Lstmp20 + Lstmp13*Lstmp149 + Lstmp13*Lstmp90 + Lstmp13*L[41] + Lstmp15*L[66] + Lstmp188*Lstmp20 + Lstmp20*L[43] + Lstmp22*L[69] + Lstmp6*Lstmp66 + Lstmp6*L[36] + Lstmp62 + Lstmp64*x + Lstmp8*L[57] + x*L[21] + L[11];
#pragma omp atomic
Ls[12] += Lstmp13*Lstmp147 + Lstmp13*Lstmp189 + Lstmp13*L[42] + Lstmp130 + Lstmp131 + Lstmp132*x + Lstmp133*x + Lstmp134*Lstmp6 + Lstmp135*Lstmp6 + Lstmp144 + Lstmp146*x + Lstmp15*L[67] + Lstmp151*Lstmp20 + Lstmp20*Lstmp88 + Lstmp20*L[44] + Lstmp22*L[70] + Lstmp6*L[37] + Lstmp8*L[58] + x*L[22] + L[12];
#pragma omp atomic
Ls[13] += Lstmp101 + Lstmp103*Lstmp20 + Lstmp13*Lstmp60 + Lstmp13*Lstmp95 + Lstmp13*L[45] + Lstmp15*L[71] + Lstmp160 + Lstmp161 + Lstmp162*x + Lstmp163*Lstmp6 + Lstmp166*Lstmp20 + Lstmp168*x + Lstmp20*L[47] + Lstmp22*L[74] + Lstmp56 + Lstmp6*Lstmp92 + Lstmp6*L[38] + Lstmp8*L[59] + Lstmp93 + L[13];
#pragma omp atomic
Ls[14] += Lstmp117*Lstmp20 + Lstmp125 + Lstmp13*Lstmp141 + Lstmp13*Lstmp158 + Lstmp13*L[46] + Lstmp143*x + Lstmp145*Lstmp6 + Lstmp15*L[72] + Lstmp155*Lstmp20 + Lstmp169 + Lstmp170 + Lstmp171*x + Lstmp172*Lstmp6 + Lstmp173*x + Lstmp20*L[48] + Lstmp22*L[75] + Lstmp6*L[39] + Lstmp8*L[60] + x*L[24] + L[14];
#pragma omp atomic
Ls[15] += Lstmp103*Lstmp13 + Lstmp13*Lstmp182 + Lstmp13*L[47] + Lstmp15*L[73] + Lstmp174 + Lstmp177 + Lstmp178*x + Lstmp179*Lstmp6 + Lstmp185*x + Lstmp20*Lstmp54 + Lstmp20*Lstmp98 + Lstmp20*L[49] + Lstmp22*L[76] + Lstmp50 + Lstmp6*Lstmp91 + Lstmp6*L[40] + Lstmp8*L[61] + Lstmp96 + Lstmp99 + L[15];
#pragma omp atomic
Ls[16] += Lstmp13*Lstmp29 + Lstmp13*Lstmp49 + Lstmp13*L[50] + Lstmp15*L[77] + Lstmp190*x + Lstmp191*x + Lstmp192*Lstmp6 + Lstmp193*Lstmp20 + Lstmp20*Lstmp86 + Lstmp20*L[52] + Lstmp22*L[80] + Lstmp27 + Lstmp47 + Lstmp48*y + Lstmp59 + Lstmp6*Lstmp90 + Lstmp6*L[41] + Lstmp8*L[62] + y*L[30] + L[16];
#pragma omp atomic
Ls[17] += Lstmp127 + Lstmp129*Lstmp13 + Lstmp13*Lstmp138 + Lstmp13*L[51] + Lstmp136 + Lstmp137*y + Lstmp140 + Lstmp147*Lstmp6 + Lstmp15*L[78] + Lstmp165*Lstmp20 + Lstmp194*x + Lstmp195*x + Lstmp196*Lstmp6 + Lstmp20*Lstmp85 + Lstmp20*L[53] + Lstmp22*L[81] + Lstmp6*L[42] + Lstmp8*L[63] + y*L[31] + L[17];
#pragma omp atomic
Ls[18] += Lstmp106 + Lstmp108*Lstmp20 + Lstmp112 + Lstmp114*Lstmp20 + Lstmp116 + Lstmp121*Lstmp6 + Lstmp13*Lstmp181 + Lstmp13*Lstmp86 + Lstmp13*L[52] + Lstmp15*L[79] + Lstmp175 + Lstmp176*x + Lstmp180*y + Lstmp184*Lstmp6 + Lstmp197*x + Lstmp20*L[54] + Lstmp22*L[82] + Lstmp6*L[43] + Lstmp8*L[64] + L[18];
#pragma omp atomic
Ls[19] += Lstmp13*Lstmp203 + Lstmp13*Lstmp85 + Lstmp13*L[53] + Lstmp15*L[80] + Lstmp198*x + Lstmp199*y + Lstmp20*Lstmp33 + Lstmp20*Lstmp45 + Lstmp20*L[55] + Lstmp201*x + Lstmp202*Lstmp6 + Lstmp22*L[83] + Lstmp31 + Lstmp43 + Lstmp53 + Lstmp6*Lstmp88 + Lstmp6*L[44] + Lstmp8*L[65] + z*L[34] + L[19];
#pragma omp atomic
Ls[20] += Lstmp13*L[59] + Lstmp20*L[61] + Lstmp38 + Lstmp39 + Lstmp40*x + Lstmp41*x + Lstmp6*L[56] + Lstmp67 + x*L[35] + L[20];
#pragma omp atomic
Ls[21] += Lstmp110 + Lstmp111*x + Lstmp120 + Lstmp13*L[62] + Lstmp20*L[64] + Lstmp6*L[57] + Lstmp64 + Lstmp66*x + x*L[36] + L[21];
#pragma omp atomic
Ls[22] += Lstmp13*L[63] + Lstmp132 + Lstmp133 + Lstmp134*x + Lstmp135*x + Lstmp146 + Lstmp20*L[65] + Lstmp6*L[58] + x*L[37] + L[22];
#pragma omp atomic
Ls[23] += Lstmp13*L[66] + Lstmp148 + Lstmp153 + Lstmp162 + Lstmp163*x + Lstmp168 + Lstmp20*L[68] + Lstmp6*L[59] + Lstmp89 + L[23];
#pragma omp atomic
Ls[24] += Lstmp13*L[67] + Lstmp143 + Lstmp145*x + Lstmp171 + Lstmp172*x + Lstmp173 + Lstmp20*L[69] + Lstmp6*L[60] + x*L[39] + L[24];
#pragma omp atomic
Ls[25] += Lstmp13*L[68] + Lstmp150 + Lstmp152 + Lstmp178 + Lstmp179*x + Lstmp185 + Lstmp20*L[70] + Lstmp6*L[61] + Lstmp87 + L[25];
#pragma omp atomic
Ls[26] += Lstmp102 + Lstmp13*L[71] + Lstmp190 + Lstmp191 + Lstmp192*x + Lstmp20*L[73] + Lstmp58 + Lstmp6*L[62] + Lstmp94 + L[26];
#pragma omp atomic
Ls[27] += Lstmp13*L[72] + Lstmp139 + Lstmp157 + Lstmp159 + Lstmp194 + Lstmp195 + Lstmp196*x + Lstmp20*L[74] + Lstmp6*L[63] + L[27];
#pragma omp atomic
Ls[28] += Lstmp115 + Lstmp13*L[73] + Lstmp154 + Lstmp156 + Lstmp176 + Lstmp184*x + Lstmp197 + Lstmp20*L[75] + Lstmp6*L[64] + L[28];
#pragma omp atomic
Ls[29] += Lstmp100 + Lstmp13*L[74] + Lstmp198 + Lstmp20*L[76] + Lstmp201 + Lstmp202*x + Lstmp52 + Lstmp6*L[65] + Lstmp97 + L[29];
#pragma omp atomic
Ls[30] += Lstmp13*L[77] + Lstmp20*L[79] + Lstmp204*x + Lstmp28 + Lstmp48 + Lstmp49*y + Lstmp6*L[66] + Lstmp61 + y*L[50] + L[30];
#pragma omp atomic
Ls[31] += Lstmp128 + Lstmp13*L[78] + Lstmp137 + Lstmp138*y + Lstmp142 + Lstmp20*L[80] + Lstmp205*x + Lstmp6*L[67] + y*L[51] + L[31];
#pragma omp atomic
Ls[32] += Lstmp13*L[79] + Lstmp164 + Lstmp167 + Lstmp180 + Lstmp181*y + Lstmp183 + Lstmp20*L[81] + Lstmp6*L[68] + Lstmp84 + L[32];
#pragma omp atomic
Ls[33] += Lstmp107 + Lstmp113 + Lstmp118 + Lstmp13*L[80] + Lstmp199 + Lstmp20*L[82] + Lstmp200*x + Lstmp203*y + Lstmp6*L[69] + L[33];
#pragma omp atomic
Ls[34] += Lstmp13*L[81] + Lstmp20*L[83] + Lstmp206*x + Lstmp207*y + Lstmp32 + Lstmp44 + Lstmp55 + Lstmp6*L[70] + z*L[55] + L[34];
#pragma omp atomic
Ls[35] += Lstmp40 + Lstmp41 + x*L[56] + L[35];
#pragma omp atomic
Ls[36] += Lstmp111 + Lstmp66 + x*L[57] + L[36];
#pragma omp atomic
Ls[37] += Lstmp134 + Lstmp135 + x*L[58] + L[37];
#pragma omp atomic
Ls[38] += Lstmp163 + Lstmp186 + Lstmp92 + L[38];
#pragma omp atomic
Ls[39] += Lstmp145 + Lstmp172 + x*L[60] + L[39];
#pragma omp atomic
Ls[40] += Lstmp179 + Lstmp187 + Lstmp91 + L[40];
#pragma omp atomic
Ls[41] += Lstmp149 + Lstmp192 + Lstmp90 + L[41];
#pragma omp atomic
Ls[42] += Lstmp147 + Lstmp189 + Lstmp196 + L[42];
#pragma omp atomic
Ls[43] += Lstmp121 + Lstmp184 + Lstmp188 + L[43];
#pragma omp atomic
Ls[44] += Lstmp151 + Lstmp202 + Lstmp88 + L[44];
#pragma omp atomic
Ls[45] += Lstmp204 + Lstmp60 + Lstmp95 + L[45];
#pragma omp atomic
Ls[46] += Lstmp141 + Lstmp158 + Lstmp205 + L[46];
#pragma omp atomic
Ls[47] += Lstmp103 + Lstmp166 + Lstmp182 + L[47];
#pragma omp atomic
Ls[48] += Lstmp117 + Lstmp155 + Lstmp200 + L[48];
#pragma omp atomic
Ls[49] += Lstmp206 + Lstmp54 + Lstmp98 + L[49];
#pragma omp atomic
Ls[50] += Lstmp29 + Lstmp49 + y*L[77] + L[50];
#pragma omp atomic
Ls[51] += Lstmp129 + Lstmp138 + y*L[78] + L[51];
#pragma omp atomic
Ls[52] += Lstmp181 + Lstmp193 + Lstmp86 + L[52];
#pragma omp atomic
Ls[53] += Lstmp165 + Lstmp203 + Lstmp85 + L[53];
#pragma omp atomic
Ls[54] += Lstmp108 + Lstmp114 + Lstmp207 + L[54];
#pragma omp atomic
Ls[55] += Lstmp33 + Lstmp45 + z*L[83] + L[55];
#pragma omp atomic
Ls[56] += L[56];
#pragma omp atomic
Ls[57] += L[57];
#pragma omp atomic
Ls[58] += L[58];
#pragma omp atomic
Ls[59] += L[59];
#pragma omp atomic
Ls[60] += L[60];
#pragma omp atomic
Ls[61] += L[61];
#pragma omp atomic
Ls[62] += L[62];
#pragma omp atomic
Ls[63] += L[63];
#pragma omp atomic
Ls[64] += L[64];
#pragma omp atomic
Ls[65] += L[65];
#pragma omp atomic
Ls[66] += L[66];
#pragma omp atomic
Ls[67] += L[67];
#pragma omp atomic
Ls[68] += L[68];
#pragma omp atomic
Ls[69] += L[69];
#pragma omp atomic
Ls[70] += L[70];
#pragma omp atomic
Ls[71] += L[71];
#pragma omp atomic
Ls[72] += L[72];
#pragma omp atomic
Ls[73] += L[73];
#pragma omp atomic
Ls[74] += L[74];
#pragma omp atomic
Ls[75] += L[75];
#pragma omp atomic
Ls[76] += L[76];
#pragma omp atomic
Ls[77] += L[77];
#pragma omp atomic
Ls[78] += L[78];
#pragma omp atomic
Ls[79] += L[79];
#pragma omp atomic
Ls[80] += L[80];
#pragma omp atomic
Ls[81] += L[81];
#pragma omp atomic
Ls[82] += L[82];
#pragma omp atomic
Ls[83] += L[83];

}

void L2P_7(double x, double y, double z, double * L, double * F) {
double Ftmp0 = x*y;
double Ftmp1 = x*z;
double Ftmp2 = y*z;
double Ftmp3 = Ftmp0*z;
double Ftmp4 = pow(x, 2);
double Ftmp5 = (1.0/2.0)*Ftmp4;
double Ftmp6 = pow(x, 3);
double Ftmp7 = (1.0/6.0)*Ftmp6;
double Ftmp8 = (1.0/24.0)*pow(x, 4);
double Ftmp9 = (1.0/120.0)*pow(x, 5);
double Ftmp10 = pow(y, 2);
double Ftmp11 = (1.0/2.0)*Ftmp10;
double Ftmp12 = pow(y, 3);
double Ftmp13 = (1.0/6.0)*Ftmp12;
double Ftmp14 = (1.0/24.0)*pow(y, 4);
double Ftmp15 = (1.0/120.0)*pow(y, 5);
double Ftmp16 = pow(z, 2);
double Ftmp17 = (1.0/2.0)*Ftmp16;
double Ftmp18 = pow(z, 3);
double Ftmp19 = (1.0/6.0)*Ftmp18;
double Ftmp20 = (1.0/24.0)*pow(z, 4);
double Ftmp21 = (1.0/120.0)*pow(z, 5);
double Ftmp22 = Ftmp11*x;
double Ftmp23 = Ftmp13*x;
double Ftmp24 = Ftmp14*x;
double Ftmp25 = Ftmp17*x;
double Ftmp26 = Ftmp19*x;
double Ftmp27 = Ftmp20*x;
double Ftmp28 = Ftmp5*y;
double Ftmp29 = Ftmp5*z;
double Ftmp30 = Ftmp7*y;
double Ftmp31 = Ftmp7*z;
double Ftmp32 = Ftmp8*y;
double Ftmp33 = Ftmp8*z;
double Ftmp34 = Ftmp17*y;
double Ftmp35 = Ftmp19*y;
double Ftmp36 = Ftmp20*y;
double Ftmp37 = Ftmp11*z;
double Ftmp38 = Ftmp13*z;
double Ftmp39 = Ftmp14*z;
double Ftmp40 = Ftmp0*Ftmp17;
double Ftmp41 = Ftmp0*Ftmp19;
double Ftmp42 = Ftmp1*Ftmp11;
double Ftmp43 = Ftmp1*Ftmp13;
double Ftmp44 = Ftmp2*Ftmp5;
double Ftmp45 = Ftmp2*Ftmp7;
double Ftmp46 = (1.0/4.0)*Ftmp4;
double Ftmp47 = Ftmp10*Ftmp46;
double Ftmp48 = (1.0/12.0)*Ftmp4;
double Ftmp49 = Ftmp12*Ftmp48;
double Ftmp50 = Ftmp16*Ftmp46;
double Ftmp51 = Ftmp18*Ftmp48;
double Ftmp52 = (1.0/12.0)*Ftmp6;
double Ftmp53 = Ftmp10*Ftmp52;
double Ftmp54 = Ftmp16*Ftmp52;
double Ftmp55 = (1.0/4.0)*Ftmp10*Ftmp16;
double Ftmp56 = (1.0/12.0)*Ftmp10*Ftmp18;
double Ftmp57 = (1.0/12.0)*Ftmp12*Ftmp16;
double Ftmp58 = Ftmp55*x;
double Ftmp59 = Ftmp50*y;
double Ftmp60 = Ftmp47*z;
#pragma omp atomic
F[0] += -Ftmp0*L[11] - Ftmp1*L[12] - Ftmp11*L[13] - Ftmp13*L[26] - Ftmp14*L[45] - Ftmp15*L[71] - Ftmp17*L[15] - Ftmp19*L[29] - Ftmp2*L[14] - Ftmp20*L[49] - Ftmp21*L[76] - Ftmp22*L[23] - Ftmp23*L[41] - Ftmp24*L[66] - Ftmp25*L[25] - Ftmp26*L[44] - Ftmp27*L[70] - Ftmp28*L[21] - Ftmp29*L[22] - Ftmp3*L[24] - Ftmp30*L[36] - Ftmp31*L[37] - Ftmp32*L[57] - Ftmp33*L[58] - Ftmp34*L[28] - Ftmp35*L[48] - Ftmp36*L[75] - Ftmp37*L[27] - Ftmp38*L[46] - Ftmp39*L[72] - Ftmp40*L[43] - Ftmp41*L[69] - Ftmp42*L[42] - Ftmp43*L[67] - Ftmp44*L[39] - Ftmp45*L[60] - Ftmp47*L[38] - Ftmp49*L[62] - Ftmp5*L[10] - Ftmp50*L[40] - Ftmp51*L[65] - Ftmp53*L[59] - Ftmp54*L[61] - Ftmp55*L[47] - Ftmp56*L[74] - Ftmp57*L[73] - Ftmp58*L[68] - Ftmp59*L[64] - Ftmp60*L[63] - Ftmp7*L[20] - Ftmp8*L[35] - Ftmp9*L[56] - x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -Ftmp0*L[13] - Ftmp1*L[14] - Ftmp11*L[16] - Ftmp13*L[30] - Ftmp14*L[50] - Ftmp15*L[77] - Ftmp17*L[18] - Ftmp19*L[33] - Ftmp2*L[17] - Ftmp20*L[54] - Ftmp21*L[82] - Ftmp22*L[26] - Ftmp23*L[45] - Ftmp24*L[71] - Ftmp25*L[28] - Ftmp26*L[48] - Ftmp27*L[75] - Ftmp28*L[23] - Ftmp29*L[24] - Ftmp3*L[27] - Ftmp30*L[38] - Ftmp31*L[39] - Ftmp32*L[59] - Ftmp33*L[60] - Ftmp34*L[32] - Ftmp35*L[53] - Ftmp36*L[81] - Ftmp37*L[31] - Ftmp38*L[51] - Ftmp39*L[78] - Ftmp40*L[47] - Ftmp41*L[74] - Ftmp42*L[46] - Ftmp43*L[72] - Ftmp44*L[42] - Ftmp45*L[63] - Ftmp47*L[41] - Ftmp49*L[66] - Ftmp5*L[11] - Ftmp50*L[43] - Ftmp51*L[69] - Ftmp53*L[62] - Ftmp54*L[64] - Ftmp55*L[52] - Ftmp56*L[80] - Ftmp57*L[79] - Ftmp58*L[73] - Ftmp59*L[68] - Ftmp60*L[67] - Ftmp7*L[21] - Ftmp8*L[36] - Ftmp9*L[57] - x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -Ftmp0*L[14] - Ftmp1*L[15] - Ftmp11*L[17] - Ftmp13*L[31] - Ftmp14*L[51] - Ftmp15*L[78] - Ftmp17*L[19] - Ftmp19*L[34] - Ftmp2*L[18] - Ftmp20*L[55] - Ftmp21*L[83] - Ftmp22*L[27] - Ftmp23*L[46] - Ftmp24*L[72] - Ftmp25*L[29] - Ftmp26*L[49] - Ftmp27*L[76] - Ftmp28*L[24] - Ftmp29*L[25] - Ftmp3*L[28] - Ftmp30*L[39] - Ftmp31*L[40] - Ftmp32*L[60] - Ftmp33*L[61] - Ftmp34*L[33] - Ftmp35*L[54] - Ftmp36*L[82] - Ftmp37*L[32] - Ftmp38*L[52] - Ftmp39*L[79] - Ftmp40*L[48] - Ftmp41*L[75] - Ftmp42*L[47] - Ftmp43*L[73] - Ftmp44*L[43] - Ftmp45*L[64] - Ftmp47*L[42] - Ftmp49*L[67] - Ftmp5*L[12] - Ftmp50*L[44] - Ftmp51*L[70] - Ftmp53*L[63] - Ftmp54*L[65] - Ftmp55*L[53] - Ftmp56*L[81] - Ftmp57*L[80] - Ftmp58*L[74] - Ftmp59*L[69] - Ftmp60*L[68] - Ftmp7*L[22] - Ftmp8*L[37] - Ftmp9*L[58] - x*L[6] - y*L[8] - z*L[9] - L[3];

}

void M2P_7(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(R, -3);
double Ftmp1 = pow(R, -2);
double Ftmp2 = 3.0*Ftmp1;
double Ftmp3 = y*M[4];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = pow(R, -4);
double Ftmp6 = Ftmp5*z;
double Ftmp7 = 15.0*Ftmp6;
double Ftmp8 = Ftmp7*M[13];
double Ftmp9 = Ftmp2*x;
double Ftmp10 = Ftmp9*y;
double Ftmp11 = Ftmp4*M[2];
double Ftmp12 = pow(x, 2);
double Ftmp13 = Ftmp1*Ftmp12;
double Ftmp14 = y*M[7];
double Ftmp15 = Ftmp7*x;
double Ftmp16 = Ftmp12*Ftmp5;
double Ftmp17 = 15.0*Ftmp16;
double Ftmp18 = pow(R, -6);
double Ftmp19 = Ftmp18*y;
double Ftmp20 = Ftmp19*z;
double Ftmp21 = 105.0*M[13];
double Ftmp22 = pow(y, 2);
double Ftmp23 = 15.0*Ftmp1;
double Ftmp24 = -Ftmp22*Ftmp23;
double Ftmp25 = Ftmp1*(Ftmp24 + 3.0);
double Ftmp26 = pow(z, 2);
double Ftmp27 = -Ftmp23*Ftmp26;
double Ftmp28 = Ftmp1*(Ftmp27 + 3.0);
double Ftmp29 = -15.0*Ftmp13;
double Ftmp30 = Ftmp1*(Ftmp29 + 9.0);
double Ftmp31 = -105.0*Ftmp13;
double Ftmp32 = Ftmp31 + 45.0;
double Ftmp33 = Ftmp32*Ftmp5;
double Ftmp34 = y*M[20];
double Ftmp35 = Ftmp33*M[21];
double Ftmp36 = Ftmp5*y;
double Ftmp37 = Ftmp1*Ftmp26;
double Ftmp38 = 3.0*M[27];
double Ftmp39 = Ftmp38*(5.0 - 35.0*Ftmp37);
double Ftmp40 = 105.0*Ftmp1;
double Ftmp41 = -Ftmp22*Ftmp40;
double Ftmp42 = Ftmp41 + 45.0;
double Ftmp43 = Ftmp36*Ftmp42;
double Ftmp44 = 1.0*M[25];
double Ftmp45 = 1.0*Ftmp6;
double Ftmp46 = Ftmp41 + 15.0;
double Ftmp47 = Ftmp46*M[26];
double Ftmp48 = -Ftmp26*Ftmp40;
double Ftmp49 = Ftmp48 + 45.0;
double Ftmp50 = Ftmp45*Ftmp49;
double Ftmp51 = Ftmp25*M[6];
double Ftmp52 = Ftmp28*M[8];
double Ftmp53 = Ftmp33*x;
double Ftmp54 = Ftmp53*y;
double Ftmp55 = Ftmp43*x;
double Ftmp56 = Ftmp46*M[16];
double Ftmp57 = Ftmp6*x;
double Ftmp58 = Ftmp53*z;
double Ftmp59 = Ftmp49*M[18];
double Ftmp60 = 315.0*Ftmp1;
double Ftmp61 = -Ftmp26*Ftmp60;
double Ftmp62 = Ftmp61 + 105.0;
double Ftmp63 = 3.0*M[47];
double Ftmp64 = Ftmp62*Ftmp63;
double Ftmp65 = -945.0*Ftmp13;
double Ftmp66 = Ftmp65 + 315.0;
double Ftmp67 = Ftmp66*M[38];
double Ftmp68 = 1.0*Ftmp19;
double Ftmp69 = Ftmp1*Ftmp22;
double Ftmp70 = -945.0*Ftmp69;
double Ftmp71 = Ftmp70 + 315.0;
double Ftmp72 = Ftmp71*M[45];
double Ftmp73 = Ftmp48 + 15.0;
double Ftmp74 = 1.0*Ftmp73*M[17];
double Ftmp75 = 1.0*Ftmp16;
double Ftmp76 = Ftmp46*M[12];
double Ftmp77 = Ftmp73*M[14];
double Ftmp78 = Ftmp19*x;
double Ftmp79 = Ftmp78*z;
double Ftmp80 = Ftmp66*Ftmp79;
double Ftmp81 = Ftmp71*M[30];
double Ftmp82 = -945.0*Ftmp37;
double Ftmp83 = Ftmp82 + 315.0;
double Ftmp84 = Ftmp83*M[32];
double Ftmp85 = Ftmp68*x;
double Ftmp86 = Ftmp85*z;
double Ftmp87 = Ftmp12*Ftmp19;
double Ftmp88 = Ftmp38*(Ftmp61 + 35.0);
double Ftmp89 = Ftmp44*Ftmp71;
double Ftmp90 = Ftmp65 + 525.0;
double Ftmp91 = Ftmp12*Ftmp18;
double Ftmp92 = Ftmp18*z;
double Ftmp93 = 1.0*Ftmp92;
double Ftmp94 = Ftmp12*Ftmp93;
double Ftmp95 = Ftmp70 + 105.0;
double Ftmp96 = Ftmp95*M[26];
double Ftmp97 = Ftmp83*M[28];
double Ftmp98 = -10395.0*Ftmp13;
double Ftmp99 = Ftmp98 + 4725.0;
double Ftmp100 = pow(R, -8);
double Ftmp101 = Ftmp100*Ftmp12;
double Ftmp102 = y*z;
double Ftmp103 = Ftmp101*Ftmp102;
double Ftmp104 = -3465.0*Ftmp37;
double Ftmp105 = Ftmp63*(Ftmp104 + 945.0);
double Ftmp106 = 1.0*Ftmp101;
double Ftmp107 = -10395.0*Ftmp69;
double Ftmp108 = Ftmp107 + 2835.0;
double Ftmp109 = Ftmp108*M[45];
double Ftmp110 = pow(y, 4);
double Ftmp111 = 945.0*Ftmp5;
double Ftmp112 = Ftmp110*Ftmp111;
double Ftmp113 = 630.0*Ftmp1;
double Ftmp114 = Ftmp5*(Ftmp112 - Ftmp113*Ftmp22 + 45.0);
double Ftmp115 = pow(z, 4);
double Ftmp116 = Ftmp111*Ftmp115;
double Ftmp117 = Ftmp5*(-Ftmp113*Ftmp26 + Ftmp116 + 45.0);
double Ftmp118 = pow(x, 4);
double Ftmp119 = Ftmp111*Ftmp118;
double Ftmp120 = Ftmp5*(Ftmp119 - 1050.0*Ftmp13 + 225.0);
double Ftmp121 = Ftmp114*M[29];
double Ftmp122 = Ftmp117*M[33];
double Ftmp123 = Ftmp115*Ftmp5;
double Ftmp124 = 3.0*M[74];
double Ftmp125 = Ftmp124*(3465.0*Ftmp123 - 1890.0*Ftmp37 + 105.0);
double Ftmp126 = Ftmp118*Ftmp5;
double Ftmp127 = 10395.0*Ftmp126;
double Ftmp128 = Ftmp127 - 9450.0*Ftmp13 + 1575.0;
double Ftmp129 = Ftmp128*M[56];
double Ftmp130 = 10395.0*Ftmp5;
double Ftmp131 = Ftmp110*Ftmp130;
double Ftmp132 = Ftmp131 - 9450.0*Ftmp69 + 1575.0;
double Ftmp133 = Ftmp132*M[70];
double Ftmp134 = 5670.0*Ftmp1;
double Ftmp135 = Ftmp131 - Ftmp134*Ftmp22 + 315.0;
double Ftmp136 = Ftmp135*M[71];
double Ftmp137 = Ftmp128*M[57];
double Ftmp138 = Ftmp115*Ftmp130;
double Ftmp139 = Ftmp138 - 9450.0*Ftmp37 + 1575.0;
double Ftmp140 = Ftmp139*Ftmp93;
double Ftmp141 = 135135.0*Ftmp126;
double Ftmp142 = -103950.0*Ftmp13 + Ftmp141 + 14175.0;
double Ftmp143 = Ftmp142*M[87];
double Ftmp144 = Ftmp100*Ftmp102;
double Ftmp145 = 45045.0*Ftmp123;
double Ftmp146 = Ftmp145 - 34650.0*Ftmp37 + 4725.0;
double Ftmp147 = 3.0*M[109];
double Ftmp148 = Ftmp146*Ftmp147;
double Ftmp149 = 1.0*Ftmp144;
double Ftmp150 = 135135.0*Ftmp5;
double Ftmp151 = Ftmp110*Ftmp150;
double Ftmp152 = Ftmp151 - 103950.0*Ftmp69 + 14175.0;
double Ftmp153 = Ftmp152*M[105];
double Ftmp154 = -Ftmp134*Ftmp26 + Ftmp138 + 315.0;
double Ftmp155 = Ftmp154*M[53];
double Ftmp156 = Ftmp128*Ftmp78;
double Ftmp157 = Ftmp132*M[49];
double Ftmp158 = Ftmp135*M[50];
double Ftmp159 = Ftmp92*x;
double Ftmp160 = Ftmp128*Ftmp159;
double Ftmp161 = Ftmp139*M[54];
double Ftmp162 = Ftmp100*x;
double Ftmp163 = Ftmp162*y;
double Ftmp164 = Ftmp163*z;
double Ftmp165 = Ftmp152*M[77];
double Ftmp166 = 1.0*Ftmp91;
double Ftmp167 = Ftmp135*M[44];
double Ftmp168 = Ftmp154*M[48];
double Ftmp169 = -145530.0*Ftmp13 + Ftmp141 + 33075.0;
double Ftmp170 = Ftmp101*y;
double Ftmp171 = Ftmp101*z;
double Ftmp172 = 1.0*Ftmp163;
double Ftmp173 = Ftmp115*Ftmp150;
double Ftmp174 = Ftmp173 - 103950.0*Ftmp37 + 14175.0;
double Ftmp175 = Ftmp174*z*M[81];
double Ftmp176 = Ftmp124*(Ftmp145 - 20790.0*Ftmp37 + 945.0);
double Ftmp177 = Ftmp106*y;
double Ftmp178 = Ftmp152*M[70];
double Ftmp179 = Ftmp106*z;
double Ftmp180 = 62370.0*Ftmp1;
double Ftmp181 = Ftmp151 - Ftmp180*Ftmp22 + 2835.0;
double Ftmp182 = Ftmp181*M[71];
double Ftmp183 = Ftmp174*M[75];
double Ftmp184 = pow(R, -10);
double Ftmp185 = Ftmp102*Ftmp12;
double Ftmp186 = Ftmp184*Ftmp185;
double Ftmp187 = 675675.0*Ftmp5;
double Ftmp188 = Ftmp115*Ftmp187;
double Ftmp189 = Ftmp147*(Ftmp188 - 450450.0*Ftmp37 + 51975.0);
double Ftmp190 = 1.0*Ftmp186;
double Ftmp191 = 2027025.0*Ftmp5;
double Ftmp192 = Ftmp110*Ftmp191;
double Ftmp193 = (Ftmp192 - 1351350.0*Ftmp69 + 155925.0)*M[105];
double Ftmp194 = 2027025.0*Ftmp126;
double Ftmp195 = Ftmp184*M[87];
double Ftmp196 = pow(y, 6);
double Ftmp197 = 135135.0*Ftmp18;
double Ftmp198 = -Ftmp196*Ftmp197;
double Ftmp199 = 155925.0*Ftmp5;
double Ftmp200 = 42525.0*Ftmp1;
double Ftmp201 = Ftmp18*(Ftmp110*Ftmp199 + Ftmp198 - Ftmp200*Ftmp22 + 1575.0);
double Ftmp202 = pow(z, 6);
double Ftmp203 = -Ftmp197*Ftmp202;
double Ftmp204 = Ftmp18*(Ftmp115*Ftmp199 - Ftmp200*Ftmp26 + Ftmp203 + 1575.0);
double Ftmp205 = pow(x, 6);
double Ftmp206 = -Ftmp197*Ftmp205;
double Ftmp207 = Ftmp18*(218295.0*Ftmp126 - 99225.0*Ftmp13 + Ftmp206 + 11025.0);
double Ftmp208 = Ftmp201*M[76];
double Ftmp209 = Ftmp204*M[82];
double Ftmp210 = 2027025.0*Ftmp18;
double Ftmp211 = -Ftmp205*Ftmp210;
double Ftmp212 = 2837835.0*Ftmp126 - 1091475.0*Ftmp13 + Ftmp211 + 99225.0;
double Ftmp213 = Ftmp163*Ftmp212;
double Ftmp214 = -Ftmp196*Ftmp210;
double Ftmp215 = Ftmp110*Ftmp5;
double Ftmp216 = Ftmp214 + 2837835.0*Ftmp215 - 1091475.0*Ftmp69 + 99225.0;
double Ftmp217 = Ftmp216*M[111];
double Ftmp218 = 467775.0*Ftmp1;
double Ftmp219 = Ftmp192 + Ftmp214 - Ftmp218*Ftmp22 + 14175.0;
double Ftmp220 = Ftmp219*M[112];
double Ftmp221 = Ftmp162*z;
double Ftmp222 = Ftmp212*Ftmp221;
double Ftmp223 = -Ftmp202*Ftmp210;
double Ftmp224 = 2837835.0*Ftmp123 + Ftmp223 - 1091475.0*Ftmp37 + 99225.0;
double Ftmp225 = Ftmp224*M[118];
double Ftmp226 = Ftmp111*Ftmp22;
double Ftmp227 = Ftmp226*Ftmp26;
double Ftmp228 = Ftmp5*(Ftmp227 + Ftmp46 + Ftmp48);
double Ftmp229 = -Ftmp22*Ftmp60;
double Ftmp230 = Ftmp12*Ftmp226;
double Ftmp231 = Ftmp5*(Ftmp229 + Ftmp230 + Ftmp32);
double Ftmp232 = Ftmp12*Ftmp26;
double Ftmp233 = Ftmp111*Ftmp232;
double Ftmp234 = Ftmp5*(Ftmp233 + Ftmp32 + Ftmp61);
double Ftmp235 = Ftmp115*Ftmp191 - Ftmp218*Ftmp26 + Ftmp223 + 14175.0;
double Ftmp236 = Ftmp235*M[117];
double Ftmp237 = Ftmp219*M[104];
double Ftmp238 = Ftmp235*M[110];
double Ftmp239 = 2835.0*Ftmp1;
double Ftmp240 = -Ftmp239*Ftmp26;
double Ftmp241 = Ftmp130*Ftmp232;
double Ftmp242 = Ftmp240 + Ftmp241;
double Ftmp243 = Ftmp242 + Ftmp66;
double Ftmp244 = Ftmp243*M[63];
double Ftmp245 = Ftmp22*Ftmp26;
double Ftmp246 = Ftmp130*Ftmp245;
double Ftmp247 = Ftmp240 + Ftmp246;
double Ftmp248 = Ftmp247 + Ftmp71;
double Ftmp249 = Ftmp248*M[72];
double Ftmp250 = -Ftmp22*Ftmp239;
double Ftmp251 = Ftmp12*Ftmp22;
double Ftmp252 = Ftmp130*Ftmp251;
double Ftmp253 = Ftmp250 + Ftmp252;
double Ftmp254 = -2835.0*Ftmp13;
double Ftmp255 = Ftmp254 + 945.0;
double Ftmp256 = Ftmp253 + Ftmp255;
double Ftmp257 = Ftmp256*M[61];
double Ftmp258 = Ftmp253 + Ftmp66;
double Ftmp259 = Ftmp258*M[62];
double Ftmp260 = Ftmp246 + Ftmp250 + Ftmp83;
double Ftmp261 = Ftmp260*M[73];
double Ftmp262 = Ftmp242 + Ftmp255;
double Ftmp263 = Ftmp262*M[64];
double Ftmp264 = -31185.0*Ftmp13;
double Ftmp265 = Ftmp264 + 8505.0;
double Ftmp266 = 31185.0*Ftmp1;
double Ftmp267 = -Ftmp22*Ftmp266;
double Ftmp268 = Ftmp150*Ftmp251;
double Ftmp269 = Ftmp267 + Ftmp268;
double Ftmp270 = Ftmp265 + Ftmp269;
double Ftmp271 = Ftmp270*M[94];
double Ftmp272 = -Ftmp26*Ftmp266;
double Ftmp273 = Ftmp150*Ftmp232;
double Ftmp274 = Ftmp272 + Ftmp273;
double Ftmp275 = Ftmp265 + Ftmp274;
double Ftmp276 = Ftmp275*M[96];
double Ftmp277 = Ftmp150*Ftmp245;
double Ftmp278 = Ftmp272 + Ftmp277;
double Ftmp279 = Ftmp267 + 8505.0;
double Ftmp280 = Ftmp278 + Ftmp279;
double Ftmp281 = Ftmp280*M[107];
double Ftmp282 = Ftmp243*Ftmp78;
double Ftmp283 = Ftmp256*Ftmp78;
double Ftmp284 = Ftmp159*Ftmp258;
double Ftmp285 = Ftmp159*Ftmp262;
double Ftmp286 = Ftmp164*Ftmp270;
double Ftmp287 = Ftmp164*Ftmp275;
double Ftmp288 = 4725.0*Ftmp1;
double Ftmp289 = -Ftmp22*Ftmp288;
double Ftmp290 = -Ftmp26*Ftmp288;
double Ftmp291 = 51975.0*Ftmp1;
double Ftmp292 = -Ftmp26*Ftmp291;
double Ftmp293 = Ftmp273 + Ftmp292;
double Ftmp294 = Ftmp293 + Ftmp99;
double Ftmp295 = -Ftmp22*Ftmp291;
double Ftmp296 = Ftmp268 + Ftmp295;
double Ftmp297 = Ftmp264 + 14175.0;
double Ftmp298 = -10395.0*Ftmp37;
double Ftmp299 = Ftmp298 + 2835.0;
double Ftmp300 = Ftmp267 + Ftmp277;
double Ftmp301 = -405405.0*Ftmp37;
double Ftmp302 = Ftmp301 + 93555.0;
double Ftmp303 = -405405.0*Ftmp69;
double Ftmp304 = Ftmp191*Ftmp245;
double Ftmp305 = Ftmp303 + Ftmp304;
double Ftmp306 = -675675.0*Ftmp69;
double Ftmp307 = Ftmp191*Ftmp251;
double Ftmp308 = -405405.0*Ftmp13;
double Ftmp309 = Ftmp308 + 155925.0;
double Ftmp310 = -675675.0*Ftmp37;
double Ftmp311 = Ftmp191*Ftmp232;
double Ftmp312 = 62370.0*Ftmp5;
double Ftmp313 = Ftmp245*Ftmp312;
double Ftmp314 = Ftmp110*Ftmp197;
double Ftmp315 = -Ftmp26*Ftmp314;
double Ftmp316 = Ftmp240 + Ftmp313 + Ftmp315;
double Ftmp317 = Ftmp18*(Ftmp135 + Ftmp316);
double Ftmp318 = Ftmp115*Ftmp197;
double Ftmp319 = -Ftmp22*Ftmp318;
double Ftmp320 = Ftmp313 + Ftmp319;
double Ftmp321 = Ftmp18*(Ftmp154 + Ftmp250 + Ftmp320);
double Ftmp322 = Ftmp251*Ftmp312;
double Ftmp323 = -Ftmp12*Ftmp314;
double Ftmp324 = Ftmp322 + Ftmp323;
double Ftmp325 = 31185.0*Ftmp5;
double Ftmp326 = 17010.0*Ftmp1;
double Ftmp327 = Ftmp110*Ftmp325 - Ftmp22*Ftmp326;
double Ftmp328 = Ftmp18*(Ftmp255 + Ftmp324 + Ftmp327);
double Ftmp329 = Ftmp232*Ftmp312;
double Ftmp330 = -Ftmp12*Ftmp318;
double Ftmp331 = Ftmp329 + Ftmp330;
double Ftmp332 = Ftmp115*Ftmp325 - Ftmp26*Ftmp326;
double Ftmp333 = Ftmp18*(Ftmp255 + Ftmp331 + Ftmp332);
double Ftmp334 = 14175.0*Ftmp1;
double Ftmp335 = -Ftmp22*Ftmp334;
double Ftmp336 = 103950.0*Ftmp5;
double Ftmp337 = Ftmp251*Ftmp336;
double Ftmp338 = Ftmp118*Ftmp197;
double Ftmp339 = -Ftmp22*Ftmp338;
double Ftmp340 = Ftmp18*(Ftmp128 + Ftmp335 + Ftmp337 + Ftmp339);
double Ftmp341 = -Ftmp26*Ftmp334;
double Ftmp342 = Ftmp232*Ftmp336;
double Ftmp343 = -Ftmp26*Ftmp338;
double Ftmp344 = Ftmp18*(Ftmp128 + Ftmp341 + Ftmp342 + Ftmp343);
double Ftmp345 = 810810.0*Ftmp5;
double Ftmp346 = Ftmp232*Ftmp345;
double Ftmp347 = Ftmp12*Ftmp210;
double Ftmp348 = -Ftmp115*Ftmp347;
double Ftmp349 = Ftmp346 + Ftmp348;
double Ftmp350 = 405405.0*Ftmp123;
double Ftmp351 = Ftmp350 - 187110.0*Ftmp37;
double Ftmp352 = Ftmp163*(Ftmp265 + Ftmp349 + Ftmp351);
double Ftmp353 = Ftmp245*Ftmp345;
double Ftmp354 = Ftmp210*Ftmp22;
double Ftmp355 = -Ftmp115*Ftmp354;
double Ftmp356 = Ftmp353 + Ftmp355;
double Ftmp357 = Ftmp163*(Ftmp279 + Ftmp351 + Ftmp356);
double Ftmp358 = Ftmp118*Ftmp210;
double Ftmp359 = -Ftmp26*Ftmp358;
double Ftmp360 = -155925.0*Ftmp37;
double Ftmp361 = 1351350.0*Ftmp5;
double Ftmp362 = Ftmp232*Ftmp361;
double Ftmp363 = Ftmp359 + Ftmp360 + Ftmp362;
double Ftmp364 = Ftmp163*(Ftmp142 + Ftmp363);
double Ftmp365 = -Ftmp110*Ftmp210*Ftmp26;
double Ftmp366 = Ftmp245*Ftmp361;
double Ftmp367 = Ftmp360 + Ftmp365 + Ftmp366;
double Ftmp368 = Ftmp163*(Ftmp152 + Ftmp367);
double Ftmp369 = 405405.0*Ftmp126;
double Ftmp370 = -155925.0*Ftmp69;
double Ftmp371 = Ftmp251*Ftmp361;
double Ftmp372 = Ftmp371 + 42525.0;
double Ftmp373 = -Ftmp22*Ftmp358;
double Ftmp374 = -311850.0*Ftmp13;
double Ftmp375 = Ftmp373 + Ftmp374;
double Ftmp376 = Ftmp163*(Ftmp369 + Ftmp370 + Ftmp372 + Ftmp375);
double Ftmp377 = 405405.0*Ftmp215;
double Ftmp378 = -155925.0*Ftmp13;
double Ftmp379 = 311850.0*Ftmp1;
double Ftmp380 = -Ftmp22*Ftmp379;
double Ftmp381 = -Ftmp110*Ftmp347;
double Ftmp382 = Ftmp380 + Ftmp381;
double Ftmp383 = Ftmp163*(Ftmp372 + Ftmp377 + Ftmp378 + Ftmp382);
double Ftmp384 = Ftmp377 - 187110.0*Ftmp69;
double Ftmp385 = Ftmp251*Ftmp345;
double Ftmp386 = Ftmp381 + Ftmp385;
double Ftmp387 = Ftmp221*(Ftmp265 + Ftmp384 + Ftmp386);
double Ftmp388 = Ftmp272 + Ftmp353 + Ftmp365;
double Ftmp389 = Ftmp221*(Ftmp384 + Ftmp388 + 8505.0);
double Ftmp390 = Ftmp221*(Ftmp142 + Ftmp370 + Ftmp371 + Ftmp373);
double Ftmp391 = Ftmp355 + Ftmp366 + Ftmp370;
double Ftmp392 = Ftmp221*(Ftmp174 + Ftmp391);
double Ftmp393 = Ftmp221*(Ftmp363 + Ftmp369 + Ftmp374 + 42525.0);
double Ftmp394 = Ftmp348 + Ftmp362 + Ftmp378;
double Ftmp395 = -Ftmp26*Ftmp379;
double Ftmp396 = Ftmp350 + Ftmp395 + 42525.0;
double Ftmp397 = Ftmp221*(Ftmp394 + Ftmp396);
double Ftmp398 = Ftmp110*Ftmp187;
double Ftmp399 = Ftmp188 + Ftmp395;
double Ftmp400 = 363825.0*Ftmp1;
double Ftmp401 = -Ftmp22*Ftmp400;
double Ftmp402 = 1891890.0*Ftmp5;
double Ftmp403 = Ftmp251*Ftmp402;
double Ftmp404 = -Ftmp26*Ftmp400;
double Ftmp405 = Ftmp232*Ftmp402;
double Ftmp406 = Ftmp173 - Ftmp180*Ftmp26 + 2835.0;
double Ftmp407 = Ftmp267 + Ftmp356;
double Ftmp408 = -Ftmp197*Ftmp22*Ftmp232;
double Ftmp409 = Ftmp18*(Ftmp243 + Ftmp245*Ftmp325 + Ftmp253 + Ftmp408);
double Ftmp410 = 405405.0*Ftmp5;
double Ftmp411 = Ftmp232*Ftmp410;
double Ftmp412 = -Ftmp232*Ftmp354;
double Ftmp413 = Ftmp245*Ftmp410 + Ftmp412;
double Ftmp414 = Ftmp163*(Ftmp270 - 93555.0*Ftmp37 + Ftmp411 + Ftmp413);
double Ftmp415 = Ftmp251*Ftmp410;
double Ftmp416 = Ftmp221*(Ftmp275 + Ftmp413 + Ftmp415 - 93555.0*Ftmp69);
double Ftmp417 = Ftmp22*Ftmp5;
double Ftmp418 = 15.0*x;
double Ftmp419 = Ftmp1*(Ftmp29 + 3.0);
double Ftmp420 = Ftmp1*(Ftmp24 + 9.0);
double Ftmp421 = Ftmp31 + 15.0;
double Ftmp422 = Ftmp421*M[23];
double Ftmp423 = Ftmp42*M[30];
double Ftmp424 = Ftmp5*x;
double Ftmp425 = Ftmp419*M[3];
double Ftmp426 = Ftmp421*M[11];
double Ftmp427 = Ftmp6*y;
double Ftmp428 = Ftmp42*Ftmp427;
double Ftmp429 = Ftmp421*M[10];
double Ftmp430 = Ftmp93*x;
double Ftmp431 = 1.0*Ftmp36;
double Ftmp432 = Ftmp18*x;
double Ftmp433 = Ftmp22*Ftmp432;
double Ftmp434 = Ftmp432*Ftmp66;
double Ftmp435 = Ftmp70 + 525.0;
double Ftmp436 = Ftmp22*Ftmp92;
double Ftmp437 = Ftmp65 + 105.0;
double Ftmp438 = Ftmp437*M[23];
double Ftmp439 = Ftmp98 + 2835.0;
double Ftmp440 = Ftmp439*M[38];
double Ftmp441 = Ftmp100*Ftmp22;
double Ftmp442 = Ftmp441*x;
double Ftmp443 = Ftmp442*z;
double Ftmp444 = Ftmp107 + 4725.0;
double Ftmp445 = 1.0*Ftmp441;
double Ftmp446 = Ftmp445*x;
double Ftmp447 = Ftmp5*(-Ftmp113*Ftmp12 + Ftmp119 + 45.0);
double Ftmp448 = Ftmp5*(Ftmp112 - 1050.0*Ftmp69 + 225.0);
double Ftmp449 = Ftmp447*M[19];
double Ftmp450 = 1.0*Ftmp432;
double Ftmp451 = Ftmp127 - 5670.0*Ftmp13 + 315.0;
double Ftmp452 = Ftmp451*M[59];
double Ftmp453 = Ftmp132*M[77];
double Ftmp454 = 1.0*Ftmp221;
double Ftmp455 = Ftmp451*M[36];
double Ftmp456 = Ftmp142*M[57];
double Ftmp457 = Ftmp18*Ftmp22;
double Ftmp458 = Ftmp451*M[35];
double Ftmp459 = Ftmp142*M[56];
double Ftmp460 = Ftmp441*z;
double Ftmp461 = -62370.0*Ftmp13 + Ftmp141 + 2835.0;
double Ftmp462 = Ftmp461*M[59];
double Ftmp463 = Ftmp151 - 145530.0*Ftmp69 + 33075.0;
double Ftmp464 = Ftmp172*z;
double Ftmp465 = Ftmp22*x*z;
double Ftmp466 = Ftmp184*Ftmp465;
double Ftmp467 = Ftmp195*(-1351350.0*Ftmp13 + Ftmp194 + 155925.0);
double Ftmp468 = 1.0*Ftmp466;
double Ftmp469 = Ftmp18*(Ftmp118*Ftmp199 - 42525.0*Ftmp13 + Ftmp206 + 1575.0);
double Ftmp470 = Ftmp18*(Ftmp198 + 218295.0*Ftmp215 - 99225.0*Ftmp69 + 11025.0);
double Ftmp471 = Ftmp469*M[55];
double Ftmp472 = -467775.0*Ftmp13 + Ftmp194 + Ftmp211 + 14175.0;
double Ftmp473 = Ftmp472*M[85];
double Ftmp474 = Ftmp5*(Ftmp233 + Ftmp31 + Ftmp73);
double Ftmp475 = -315.0*Ftmp13;
double Ftmp476 = Ftmp5*(Ftmp230 + Ftmp42 + Ftmp475);
double Ftmp477 = Ftmp5*(Ftmp227 + Ftmp42 + Ftmp61);
double Ftmp478 = Ftmp472*M[84];
double Ftmp479 = Ftmp252 + Ftmp254;
double Ftmp480 = Ftmp479 + Ftmp71;
double Ftmp481 = Ftmp480*M[66];
double Ftmp482 = Ftmp241 + Ftmp254;
double Ftmp483 = Ftmp482 + Ftmp83;
double Ftmp484 = Ftmp483*M[68];
double Ftmp485 = Ftmp250 + 945.0;
double Ftmp486 = Ftmp247 + Ftmp485;
double Ftmp487 = Ftmp486*M[79];
double Ftmp488 = Ftmp20*Ftmp480;
double Ftmp489 = Ftmp20*Ftmp483;
double Ftmp490 = Ftmp20*Ftmp486;
double Ftmp491 = -4725.0*Ftmp13;
double Ftmp492 = -51975.0*Ftmp13;
double Ftmp493 = Ftmp492 + 14175.0;
double Ftmp494 = Ftmp268 + Ftmp444 + Ftmp492;
double Ftmp495 = Ftmp280*Ftmp464;
double Ftmp496 = Ftmp277 + Ftmp292;
double Ftmp497 = Ftmp303 + Ftmp307;
double Ftmp498 = 155925.0 - 675675.0*Ftmp13;
double Ftmp499 = Ftmp240 + Ftmp329 + Ftmp343;
double Ftmp500 = Ftmp18*(Ftmp451 + Ftmp499);
double Ftmp501 = Ftmp18*(Ftmp154 + Ftmp254 + Ftmp331);
double Ftmp502 = Ftmp322 + Ftmp339;
double Ftmp503 = 31185.0*Ftmp126 - 17010.0*Ftmp13;
double Ftmp504 = Ftmp18*(Ftmp485 + Ftmp502 + Ftmp503);
double Ftmp505 = Ftmp18*(Ftmp320 + Ftmp332 + Ftmp485);
double Ftmp506 = -14175.0*Ftmp13;
double Ftmp507 = Ftmp18*(Ftmp132 + Ftmp323 + Ftmp337 + Ftmp506);
double Ftmp508 = Ftmp245*Ftmp336;
double Ftmp509 = Ftmp18*(Ftmp132 + Ftmp315 + Ftmp341 + Ftmp508);
double Ftmp510 = -187110.0*Ftmp13 + Ftmp369;
double Ftmp511 = Ftmp144*(Ftmp279 + Ftmp373 + Ftmp385 + Ftmp510);
double Ftmp512 = Ftmp272 + Ftmp346 + Ftmp359;
double Ftmp513 = Ftmp144*(Ftmp510 + Ftmp512 + 8505.0);
double Ftmp514 = Ftmp144*(Ftmp152 + Ftmp371 + Ftmp378 + Ftmp381);
double Ftmp515 = Ftmp144*(Ftmp174 + Ftmp394);
double Ftmp516 = Ftmp144*(Ftmp367 + Ftmp377 + Ftmp380 + 42525.0);
double Ftmp517 = Ftmp144*(Ftmp391 + Ftmp396);
double Ftmp518 = Ftmp267 + Ftmp385;
double Ftmp519 = 675675.0*Ftmp126 + 14175.0;
double Ftmp520 = -363825.0*Ftmp13;
double Ftmp521 = Ftmp245*Ftmp402;
double Ftmp522 = 1.0*M[108];
double Ftmp523 = 1.0*M[106];
double Ftmp524 = Ftmp18*(Ftmp232*Ftmp325 + Ftmp248 + Ftmp408 + Ftmp479);
double Ftmp525 = Ftmp144*(-93555.0*Ftmp13 + Ftmp280 + Ftmp411 + Ftmp412 + Ftmp415);
double Ftmp526 = Ftmp26*Ftmp5;
double Ftmp527 = Ftmp1*(Ftmp27 + 9.0);
double Ftmp528 = 1.0*Ftmp424;
double Ftmp529 = Ftmp26*Ftmp450;
double Ftmp530 = Ftmp82 + 525.0;
double Ftmp531 = Ftmp19*Ftmp26;
double Ftmp532 = Ftmp100*Ftmp26;
double Ftmp533 = Ftmp532*x;
double Ftmp534 = Ftmp533*y;
double Ftmp535 = 1.0*Ftmp533;
double Ftmp536 = Ftmp5*(Ftmp116 - 1050.0*Ftmp37 + 225.0);
double Ftmp537 = Ftmp139*Ftmp68;
double Ftmp538 = Ftmp18*Ftmp26;
double Ftmp539 = Ftmp532*y;
double Ftmp540 = Ftmp173 - 145530.0*Ftmp37 + 33075.0;
double Ftmp541 = Ftmp26*x*y;
double Ftmp542 = Ftmp184*Ftmp541;
double Ftmp543 = 1.0*Ftmp542;
double Ftmp544 = Ftmp18*(218295.0*Ftmp123 + Ftmp203 - 99225.0*Ftmp37 + 11025.0);
double Ftmp545 = Ftmp5*(Ftmp230 + Ftmp31 + Ftmp46);
double Ftmp546 = Ftmp5*(Ftmp233 + Ftmp475 + Ftmp49);
double Ftmp547 = Ftmp5*(Ftmp227 + Ftmp229 + Ftmp49);
double Ftmp548 = Ftmp298 + 4725.0;
double Ftmp549 = Ftmp273 + Ftmp492 + Ftmp548;
double Ftmp550 = Ftmp277 + Ftmp295;
double Ftmp551 = Ftmp18*(Ftmp250 + Ftmp451 + Ftmp502);
double Ftmp552 = Ftmp18*(Ftmp135 + Ftmp254 + Ftmp324);
double Ftmp553 = Ftmp18*(Ftmp499 + Ftmp503 + 945.0);
double Ftmp554 = Ftmp18*(Ftmp316 + Ftmp327 + 945.0);
double Ftmp555 = Ftmp18*(Ftmp139 + Ftmp330 + Ftmp342 + Ftmp506);
double Ftmp556 = Ftmp18*(Ftmp139 + Ftmp319 + Ftmp335 + Ftmp508);
double Ftmp557 = Ftmp18*(Ftmp251*Ftmp325 + Ftmp260 + Ftmp408 + Ftmp482);
#pragma omp atomic
F[0] += Ftmp0*(-Ftmp10*M[1] + Ftmp101*(Ftmp297 + Ftmp349 + Ftmp399)*M[97] + Ftmp101*(3648645.0*Ftmp126 - 1964655.0*Ftmp13 + Ftmp211 + 297675.0)*M[83] + Ftmp101*(Ftmp169 + Ftmp359 + Ftmp404 + Ftmp405)*M[88] + Ftmp101*(Ftmp169 + Ftmp373 + Ftmp401 + Ftmp403)*M[86] + Ftmp101*(Ftmp297 + Ftmp382 + Ftmp385 + Ftmp398)*M[93] + Ftmp101*(Ftmp187*Ftmp245 + Ftmp294 + Ftmp296 + Ftmp412)*M[95] + Ftmp102*Ftmp106*Ftmp109 + Ftmp103*Ftmp105 + Ftmp103*Ftmp99*M[38] + Ftmp106*Ftmp237 + Ftmp106*Ftmp238 + Ftmp106*(Ftmp181 + Ftmp388)*M[106] + Ftmp106*(Ftmp406 + Ftmp407)*M[108] - Ftmp11*x + Ftmp114*M[44] + Ftmp117*M[48] - Ftmp12*Ftmp20*Ftmp21 - Ftmp12*Ftmp90*Ftmp92*M[21] + Ftmp120*x*M[19] + Ftmp120*M[34] + Ftmp121*x + Ftmp122*x - Ftmp125*Ftmp19 - Ftmp129*Ftmp19 - 3.0*Ftmp13*M[0] - Ftmp133*Ftmp68 - Ftmp136*Ftmp93 - Ftmp137*Ftmp92 + Ftmp14*Ftmp15 - Ftmp140*M[75] + Ftmp142*Ftmp164*M[59] + Ftmp143*Ftmp144 + Ftmp144*Ftmp148 + Ftmp144*Ftmp271 + Ftmp144*Ftmp276 + Ftmp149*Ftmp153 + Ftmp149*Ftmp281 - Ftmp155*Ftmp85 - Ftmp156*M[35] - Ftmp157*Ftmp78 - Ftmp158*Ftmp159 - Ftmp159*Ftmp161 - Ftmp159*Ftmp260*M[52] + Ftmp16*(Ftmp31 + 75.0)*M[9] - Ftmp160*M[36] + Ftmp163*Ftmp217 + Ftmp164*Ftmp165 + Ftmp164*Ftmp280*M[79] - Ftmp166*Ftmp167 - Ftmp166*Ftmp168 - Ftmp166*(Ftmp246 + Ftmp82 + Ftmp95)*M[46] + Ftmp169*Ftmp170*M[56] + Ftmp169*Ftmp171*M[57] + Ftmp17*Ftmp3 + Ftmp17*z*M[5] + Ftmp170*Ftmp176 + Ftmp170*Ftmp294*M[63] + Ftmp170*(Ftmp296 + Ftmp297)*M[61] + Ftmp171*(Ftmp293 + Ftmp297)*M[64] + Ftmp171*(Ftmp296 + Ftmp99)*M[62] + Ftmp172*Ftmp175 + Ftmp172*Ftmp236 + Ftmp177*Ftmp178 + Ftmp177*(Ftmp108 + Ftmp278)*M[72] + Ftmp179*Ftmp182 + Ftmp179*Ftmp183 + Ftmp179*(Ftmp299 + Ftmp300)*M[73] - Ftmp185*Ftmp195*(-1891890.0*Ftmp13 + Ftmp194 + 363825.0) - Ftmp186*Ftmp189 - Ftmp186*(Ftmp306 + Ftmp307 + Ftmp309)*M[94] - Ftmp186*(Ftmp309 + Ftmp310 + Ftmp311)*M[96] - Ftmp19*Ftmp244 - Ftmp19*Ftmp257 - Ftmp190*Ftmp193 - Ftmp190*(Ftmp302 + Ftmp305)*M[107] - Ftmp2*Ftmp3 - Ftmp20*Ftmp64 - Ftmp20*Ftmp67 - Ftmp201*M[104] - Ftmp204*M[110] - Ftmp207*x*M[55] - Ftmp207*M[83] - Ftmp208*x - Ftmp209*x + Ftmp213*M[84] + Ftmp220*Ftmp221 + Ftmp221*Ftmp225 + Ftmp222*M[85] + Ftmp228*x*M[31] + Ftmp228*M[46] + Ftmp231*x*M[22] + Ftmp231*M[37] + Ftmp234*x*M[24] + Ftmp234*M[39] - Ftmp248*Ftmp78*M[51] - Ftmp249*Ftmp68 - Ftmp25*M[12] - Ftmp259*Ftmp92 - Ftmp261*Ftmp93 - Ftmp263*Ftmp92 - Ftmp28*M[14] - Ftmp282*M[42] - Ftmp283*M[40] - Ftmp284*M[41] - Ftmp285*M[43] + Ftmp286*M[66] + Ftmp287*M[68] - Ftmp30*x*M[3] - Ftmp30*M[9] - Ftmp317*x*M[78] - Ftmp317*M[106] - Ftmp321*x*M[80] - Ftmp321*M[108] - Ftmp328*x*M[65] - Ftmp328*M[93] + Ftmp33*Ftmp34 - Ftmp333*x*M[69] - Ftmp333*M[97] - Ftmp34*Ftmp90*Ftmp91 - Ftmp340*x*M[58] - Ftmp340*M[86] - Ftmp344*x*M[60] - Ftmp344*M[88] + Ftmp35*z + Ftmp352*M[102] + Ftmp357*M[115] + Ftmp36*Ftmp39 + Ftmp36*Ftmp74*x + Ftmp364*M[91] + Ftmp368*M[113] + Ftmp376*M[89] + Ftmp383*M[98] + Ftmp387*M[99] + Ftmp389*M[114] + Ftmp390*M[90] + Ftmp392*M[116] + Ftmp393*M[92] + Ftmp397*M[103] - Ftmp4*M[5] - Ftmp409*x*M[67] - Ftmp409*M[95] + Ftmp414*M[100] + Ftmp416*M[101] + Ftmp43*Ftmp44 + Ftmp45*Ftmp47 + Ftmp50*M[28] - Ftmp51*x - Ftmp52*x + Ftmp54*M[10] + Ftmp55*M[15] + Ftmp56*Ftmp57 + Ftmp57*Ftmp59 + Ftmp58*M[11] - Ftmp68*Ftmp72*z + Ftmp75*Ftmp76 + Ftmp75*Ftmp77 - Ftmp79*Ftmp81 + Ftmp8*y - Ftmp80*M[23] - Ftmp84*Ftmp86 - Ftmp87*Ftmp88 - Ftmp87*Ftmp89 - Ftmp91*(Ftmp127 - 13230.0*Ftmp13 + 3675.0)*M[34] - Ftmp91*(Ftmp241 + Ftmp290 + Ftmp90)*M[39] - Ftmp91*(Ftmp252 + Ftmp289 + Ftmp90)*M[37] - Ftmp94*Ftmp96 - Ftmp94*Ftmp97 + 1.0*M[0]);
#pragma omp atomic
F[1] += Ftmp0*(-Ftmp10*M[0] + Ftmp105*Ftmp443 - Ftmp11*y + Ftmp117*M[53] + Ftmp122*y - Ftmp125*Ftmp432 - Ftmp129*Ftmp432 - Ftmp132*Ftmp20*M[50] - Ftmp132*Ftmp85*M[44] - Ftmp133*Ftmp450 - Ftmp140*M[81] + Ftmp143*Ftmp221 + Ftmp144*Ftmp216*M[112] + Ftmp144*Ftmp225 + Ftmp144*Ftmp473 + Ftmp148*Ftmp221 + Ftmp15*y*M[5] + Ftmp152*Ftmp464*M[71] + Ftmp153*Ftmp454 - 1.0*Ftmp155*Ftmp457 - Ftmp156*M[34] - Ftmp159*Ftmp21*Ftmp22 - Ftmp159*Ftmp64 - Ftmp159*Ftmp67 - Ftmp161*Ftmp20 + Ftmp164*Ftmp456 - Ftmp168*Ftmp85 + Ftmp172*Ftmp216*M[104] + Ftmp172*Ftmp238 + Ftmp175*Ftmp445 + Ftmp176*Ftmp442 + Ftmp183*Ftmp464 - Ftmp189*Ftmp466 - Ftmp20*Ftmp455 - Ftmp204*M[117] - Ftmp209*y + Ftmp213*M[83] - Ftmp22*Ftmp434*M[20] + Ftmp22*Ftmp7*M[7] - Ftmp22*Ftmp84*Ftmp93 + Ftmp221*Ftmp271 + Ftmp221*Ftmp276 + Ftmp236*Ftmp445 - Ftmp244*Ftmp432 - Ftmp248*Ftmp85*M[46] - Ftmp249*Ftmp450 - Ftmp257*Ftmp432 - Ftmp28*M[17] + Ftmp281*Ftmp454 - Ftmp282*M[39] - Ftmp283*M[37] + Ftmp286*M[62] + Ftmp287*M[64] + Ftmp352*M[97] + Ftmp357*Ftmp522 + Ftmp364*M[88] + Ftmp368*Ftmp523 + Ftmp376*M[86] + Ftmp383*M[93] + Ftmp39*Ftmp424 - Ftmp4*M[7] + Ftmp414*M[95] + Ftmp417*Ftmp418*M[4] + Ftmp417*Ftmp429 + Ftmp417*Ftmp74 + Ftmp417*(Ftmp41 + 75.0)*M[15] - Ftmp419*M[10] + Ftmp42*Ftmp424*Ftmp44 - Ftmp420*y*M[6] - Ftmp420*M[15] + Ftmp422*Ftmp6 + Ftmp423*Ftmp6 - Ftmp425*y + Ftmp426*Ftmp427 + Ftmp427*Ftmp59 + Ftmp428*M[16] - Ftmp430*Ftmp72 + Ftmp431*Ftmp77*x - Ftmp433*Ftmp435*Ftmp44 - Ftmp433*Ftmp88 - Ftmp435*Ftmp436*M[30] - Ftmp436*Ftmp438 + Ftmp440*Ftmp443 + Ftmp441*Ftmp478 + Ftmp441*(Ftmp461 + Ftmp512)*M[91] + Ftmp441*(Ftmp264 + Ftmp349 + Ftmp406)*M[102] + Ftmp441*(Ftmp375 + Ftmp518 + Ftmp519)*M[89] + Ftmp441*(Ftmp399 + Ftmp407 + 14175.0)*M[115] + Ftmp441*(Ftmp214 + 3648645.0*Ftmp215 - 1964655.0*Ftmp69 + 297675.0)*M[111] + Ftmp441*(Ftmp365 + Ftmp404 + Ftmp463 + Ftmp521)*M[113] + Ftmp441*(Ftmp381 + Ftmp403 + Ftmp463 + Ftmp520)*M[98] + Ftmp441*(Ftmp187*Ftmp232 + Ftmp412 + Ftmp494 + Ftmp496)*M[100] + Ftmp442*Ftmp459 + Ftmp442*(Ftmp269 + Ftmp493)*M[61] + Ftmp442*(Ftmp274 + Ftmp439)*M[63] + Ftmp444*Ftmp446*z*M[45] + Ftmp446*Ftmp463*M[70] + Ftmp446*(Ftmp444 + Ftmp496)*M[72] + Ftmp447*M[35] + Ftmp448*y*M[29] + Ftmp448*M[49] + Ftmp449*y - Ftmp452*Ftmp92 - Ftmp453*Ftmp92 - Ftmp457*Ftmp458 - Ftmp457*(Ftmp131 - 13230.0*Ftmp69 + 3675.0)*M[49] - Ftmp457*(Ftmp241 + Ftmp437 + Ftmp82)*M[42] - Ftmp457*(Ftmp246 + Ftmp290 + Ftmp435)*M[51] - Ftmp457*(Ftmp252 + Ftmp435 + Ftmp491)*M[40] + Ftmp460*Ftmp462 + Ftmp460*Ftmp463*M[77] + Ftmp460*Ftmp494*M[66] + Ftmp460*(Ftmp264 + Ftmp273 + Ftmp299)*M[68] + Ftmp460*(Ftmp292 + Ftmp300 + 14175.0)*M[79] - Ftmp465*Ftmp467 - Ftmp466*(Ftmp497 + Ftmp498)*M[94] - Ftmp466*(Ftmp302 + Ftmp308 + Ftmp311)*M[96] - Ftmp468*(Ftmp192 - 1891890.0*Ftmp69 + 363825.0)*M[105] - Ftmp468*(Ftmp305 + Ftmp310 + 155925.0)*M[107] - Ftmp469*M[84] - Ftmp470*y*M[76] - Ftmp470*M[111] - Ftmp471*y + Ftmp474*y*M[24] + Ftmp474*M[42] + Ftmp476*y*M[22] + Ftmp476*M[40] + Ftmp477*y*M[31] + Ftmp477*M[51] - Ftmp481*Ftmp92 - Ftmp484*Ftmp92 - Ftmp487*Ftmp92 - Ftmp488*M[41] - Ftmp489*M[43] - Ftmp490*M[52] + Ftmp495*M[73] + Ftmp50*M[32] - Ftmp500*y*M[60] - Ftmp500*M[91] - Ftmp501*y*M[69] - Ftmp501*M[102] - Ftmp504*y*M[58] - Ftmp504*M[89] - Ftmp505*y*M[80] - Ftmp505*M[115] - Ftmp507*y*M[65] - Ftmp507*M[98] - Ftmp509*y*M[78] - Ftmp509*M[113] + Ftmp511*M[90] + Ftmp513*M[92] + Ftmp514*M[99] + Ftmp515*M[103] + Ftmp516*M[114] + Ftmp517*M[116] - Ftmp52*y - Ftmp524*y*M[67] - Ftmp524*M[100] + Ftmp525*M[101] + Ftmp53*M[20] + Ftmp54*M[9] + 1.0*Ftmp55*M[12] - 3.0*Ftmp69*M[1] - Ftmp71*Ftmp86*M[26] + Ftmp8*x - Ftmp80*M[21] - Ftmp86*Ftmp97 - Ftmp9*M[4] + 1.0*M[1]);
#pragma omp atomic
F[2] += Ftmp0*(Ftmp109*Ftmp535*y + Ftmp114*M[50] + Ftmp121*z + Ftmp124*Ftmp146*Ftmp164 - Ftmp136*Ftmp450 - Ftmp137*Ftmp432 - Ftmp139*Ftmp450*M[75] - Ftmp14*Ftmp2 + 15.0*Ftmp14*Ftmp526 - Ftmp140*x*M[48] + Ftmp143*Ftmp163 + Ftmp144*Ftmp217 + Ftmp144*Ftmp478 - Ftmp147*Ftmp542*(Ftmp188 - 630630.0*Ftmp37 + 121275.0) + Ftmp148*Ftmp163 + Ftmp149*Ftmp224*M[117] + Ftmp15*Ftmp3 + Ftmp153*Ftmp172 - Ftmp157*Ftmp20 - Ftmp158*Ftmp538 - Ftmp159*Ftmp34*Ftmp66 - Ftmp160*M[34] + Ftmp163*Ftmp271 + Ftmp163*Ftmp276 + Ftmp164*Ftmp459 + Ftmp165*Ftmp539 - Ftmp167*Ftmp430 + Ftmp172*Ftmp281 + Ftmp178*Ftmp464 + Ftmp182*Ftmp535 - Ftmp19*Ftmp452 - Ftmp19*Ftmp453 - Ftmp19*Ftmp481 - Ftmp19*Ftmp484 - Ftmp19*Ftmp487 - Ftmp193*Ftmp543 - Ftmp20*Ftmp458 - Ftmp201*M[112] - Ftmp208*z - Ftmp21*Ftmp26*Ftmp78 + Ftmp220*Ftmp532 + Ftmp222*M[83] + Ftmp224*Ftmp454*M[110] + Ftmp237*Ftmp454 - Ftmp25*M[16] - Ftmp259*Ftmp432 - Ftmp26*Ftmp434*M[21] - Ftmp26*Ftmp530*Ftmp68*M[32] - Ftmp260*Ftmp430*M[46] - Ftmp261*Ftmp450 - Ftmp263*Ftmp432 - Ftmp284*M[37] - Ftmp285*M[39] + Ftmp286*M[61] + Ftmp287*M[63] + Ftmp35*x + Ftmp36*Ftmp418*M[13] + Ftmp36*Ftmp422 + Ftmp36*Ftmp423 - 3.0*Ftmp37*M[2] - Ftmp38*Ftmp62*Ftmp79 + Ftmp387*M[93] + Ftmp389*Ftmp523 + Ftmp390*M[86] + Ftmp392*Ftmp522 + Ftmp393*M[88] + Ftmp397*M[97] - Ftmp4*x*M[0] - Ftmp4*y*M[1] + Ftmp416*M[95] + Ftmp418*Ftmp526*M[5] - Ftmp419*M[11] - Ftmp425*z + Ftmp426*Ftmp526 + Ftmp427*Ftmp429 + Ftmp428*M[15] + Ftmp431*Ftmp49*M[32] - Ftmp438*Ftmp531 + Ftmp440*Ftmp534 + Ftmp447*M[36] + Ftmp449*z + Ftmp45*Ftmp76*x - Ftmp455*Ftmp538 + Ftmp456*Ftmp533 + Ftmp462*Ftmp539 - Ftmp467*Ftmp541 - Ftmp469*M[85] + Ftmp47*Ftmp528 - Ftmp471*z + Ftmp473*Ftmp532 - Ftmp488*M[40] - Ftmp489*M[42] + Ftmp49*Ftmp528*M[28] - Ftmp490*M[51] + Ftmp495*M[72] + Ftmp50*x*M[14] + Ftmp50*y*M[17] - Ftmp51*z + Ftmp511*M[89] + Ftmp513*M[91] + Ftmp514*M[98] + Ftmp515*M[102] + Ftmp516*M[113] + Ftmp517*M[115] + Ftmp525*M[100] + Ftmp526*Ftmp56 + Ftmp526*(Ftmp48 + 75.0)*M[18] - Ftmp527*z*M[8] - Ftmp527*M[18] - Ftmp529*Ftmp530*M[28] - Ftmp529*Ftmp96 - Ftmp531*Ftmp81 + Ftmp532*(Ftmp181 + Ftmp264 + Ftmp386)*M[99] + Ftmp532*(Ftmp373 + Ftmp461 + Ftmp518)*M[90] + Ftmp532*(Ftmp374 + Ftmp512 + Ftmp519)*M[92] + Ftmp532*(3648645.0*Ftmp123 + Ftmp223 - 1964655.0*Ftmp37 + 297675.0)*M[118] + Ftmp532*(Ftmp348 + Ftmp405 + Ftmp520 + Ftmp540)*M[103] + Ftmp532*(Ftmp355 + Ftmp401 + Ftmp521 + Ftmp540)*M[116] + Ftmp532*(Ftmp380 + Ftmp388 + Ftmp398 + 14175.0)*M[114] + Ftmp532*(Ftmp187*Ftmp251 + Ftmp412 + Ftmp549 + Ftmp550)*M[101] + Ftmp533*(Ftmp269 + Ftmp439)*M[62] + Ftmp533*(Ftmp274 + Ftmp493)*M[64] + Ftmp534*Ftmp63*(Ftmp104 + 1575.0) + Ftmp535*Ftmp540*M[75] + Ftmp535*(Ftmp548 + Ftmp550)*M[73] + Ftmp536*z*M[33] + Ftmp536*M[54] - Ftmp537*z*M[53] - Ftmp537*M[81] - Ftmp538*(Ftmp138 - 13230.0*Ftmp37 + 3675.0)*M[54] - Ftmp538*(Ftmp241 + Ftmp491 + Ftmp530)*M[43] - Ftmp538*(Ftmp246 + Ftmp289 + Ftmp530)*M[52] - Ftmp538*(Ftmp252 + Ftmp65 + Ftmp95)*M[41] + 1.0*Ftmp539*Ftmp540*M[81] + Ftmp539*Ftmp549*M[68] + Ftmp539*(Ftmp108 + Ftmp264 + Ftmp268)*M[66] + Ftmp539*(Ftmp278 + Ftmp295 + 14175.0)*M[79] - Ftmp542*(Ftmp301 + Ftmp311 + Ftmp498)*M[96] - Ftmp542*(Ftmp308 + Ftmp497 + 93555.0)*M[94] - Ftmp543*(Ftmp301 + Ftmp304 + Ftmp306 + 155925.0)*M[107] - Ftmp544*z*M[82] - Ftmp544*M[118] + Ftmp545*z*M[22] + Ftmp545*M[41] + Ftmp546*z*M[24] + Ftmp546*M[43] + Ftmp547*z*M[31] + Ftmp547*M[52] - Ftmp551*z*M[58] - Ftmp551*M[90] - Ftmp552*z*M[65] - Ftmp552*M[99] - Ftmp553*z*M[60] - Ftmp553*M[92] - Ftmp554*z*M[78] - Ftmp554*M[114] - Ftmp555*z*M[69] - Ftmp555*M[103] - Ftmp556*z*M[80] - Ftmp556*M[116] - Ftmp557*z*M[67] - Ftmp557*M[101] + Ftmp58*M[9] - Ftmp64*Ftmp78 - Ftmp67*Ftmp78 - Ftmp72*Ftmp85 - Ftmp79*Ftmp89 - Ftmp9*M[5] + 1.0*M[2]);

}

void P2M(double x, double y, double z, double q, double * M, int order) {
switch (order) {
  case 2:
    P2M_2(x, y, z, q, M);
    break;
  case 3:
    P2M_3(x, y, z, q, M);
    break;
  case 4:
    P2M_4(x, y, z, q, M);
    break;
  case 5:
    P2M_5(x, y, z, q, M);
    break;
  case 6:
    P2M_6(x, y, z, q, M);
    break;
  case 7:
    P2M_7(x, y, z, q, M);
    break;
  }
}
void M2M(double x, double y, double z, double * M, double * Ms, int order) {
switch (order) {
  case 2:
    M2M_2(x, y, z, M, Ms);
    break;
  case 3:
    M2M_3(x, y, z, M, Ms);
    break;
  case 4:
    M2M_4(x, y, z, M, Ms);
    break;
  case 5:
    M2M_5(x, y, z, M, Ms);
    break;
  case 6:
    M2M_6(x, y, z, M, Ms);
    break;
  case 7:
    M2M_7(x, y, z, M, Ms);
    break;
  }
}
void M2L(double x, double y, double z, double * M, double * L, int order) {
switch (order) {
  case 2:
    M2L_2(x, y, z, M, L);
    break;
  case 3:
    M2L_3(x, y, z, M, L);
    break;
  case 4:
    M2L_4(x, y, z, M, L);
    break;
  case 5:
    M2L_5(x, y, z, M, L);
    break;
  case 6:
    M2L_6(x, y, z, M, L);
    break;
  case 7:
    M2L_7(x, y, z, M, L);
    break;
  }
}
void L2L(double x, double y, double z, double * L, double * Ls, int order) {
switch (order) {
  case 2:
    L2L_2(x, y, z, L, Ls);
    break;
  case 3:
    L2L_3(x, y, z, L, Ls);
    break;
  case 4:
    L2L_4(x, y, z, L, Ls);
    break;
  case 5:
    L2L_5(x, y, z, L, Ls);
    break;
  case 6:
    L2L_6(x, y, z, L, Ls);
    break;
  case 7:
    L2L_7(x, y, z, L, Ls);
    break;
  }
}
void L2P(double x, double y, double z, double * L, double * F, int order) {
switch (order) {
  case 2:
    L2P_2(x, y, z, L, F);
    break;
  case 3:
    L2P_3(x, y, z, L, F);
    break;
  case 4:
    L2P_4(x, y, z, L, F);
    break;
  case 5:
    L2P_5(x, y, z, L, F);
    break;
  case 6:
    L2P_6(x, y, z, L, F);
    break;
  case 7:
    L2P_7(x, y, z, L, F);
    break;
  }
}
void M2P(double x, double y, double z, double * M, double * F, int order) {
switch (order) {
  case 2:
    M2P_2(x, y, z, M, F);
    break;
  case 3:
    M2P_3(x, y, z, M, F);
    break;
  case 4:
    M2P_4(x, y, z, M, F);
    break;
  case 5:
    M2P_5(x, y, z, M, F);
    break;
  case 6:
    M2P_6(x, y, z, M, F);
    break;
  case 7:
    M2P_7(x, y, z, M, F);
    break;
  }
}
