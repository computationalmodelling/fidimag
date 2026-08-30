#include "operators.h"
#include "math.h"
void S2M_2(double x, double y, double z, double * S, double * M) {
#pragma omp atomic
M[0] += S[0];
#pragma omp atomic
M[1] += S[1];
#pragma omp atomic
M[2] += S[2];
#pragma omp atomic
M[3] += x*S[0];
#pragma omp atomic
M[4] += x*S[1] + y*S[0];
#pragma omp atomic
M[5] += x*S[2] + z*S[0];
#pragma omp atomic
M[6] += y*S[1];
#pragma omp atomic
M[7] += y*S[2] + z*S[1];
#pragma omp atomic
M[8] += z*S[2];

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
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double D[9];
double Dtmp0 = pow(Rinv, 3);
double Dtmp1 = 3*pow(Rinv, 2);
double Dtmp2 = pow(Rinv, 5);
double Dtmp3 = 3*Dtmp2*x;
D[0] = -Dtmp0*x;
D[1] = -Dtmp0*y;
D[2] = -Dtmp0*z;
D[3] = Dtmp0*(Dtmp1*pow(x, 2) - 1);
D[4] = Dtmp3*y;
D[5] = Dtmp3*z;
D[6] = Dtmp0*(Dtmp1*pow(y, 2) - 1);
D[7] = 3*Dtmp2*y*z;
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
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(Rinv, 3);
double Ftmp1 = pow(Rinv, 2);
double Ftmp2 = 3*Ftmp1;
double Ftmp3 = y*M[4];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = Ftmp2*x;
double Ftmp6 = Ftmp5*y;
double Ftmp7 = Ftmp4*M[2];
double Ftmp8 = pow(x, 2);
double Ftmp9 = Ftmp1*Ftmp8;
double Ftmp10 = 3*Ftmp9;
double Ftmp11 = z*M[7];
double Ftmp12 = pow(Rinv, 4);
double Ftmp13 = 15*Ftmp12;
double Ftmp14 = Ftmp13*x;
double Ftmp15 = Ftmp14*y;
double Ftmp16 = 6*x;
double Ftmp17 = Ftmp12*Ftmp16;
double Ftmp18 = pow(y, 2);
double Ftmp19 = Ftmp18*M[6];
double Ftmp20 = pow(z, 2);
double Ftmp21 = Ftmp20*M[8];
double Ftmp22 = Ftmp12*Ftmp8;
double Ftmp23 = 15*Ftmp22;
double Ftmp24 = z*M[5];
double Ftmp25 = (Ftmp10 - 1)*M[3];
double Ftmp26 = Ftmp1*Ftmp18;
double Ftmp27 = 3*Ftmp26;
double Ftmp28 = (Ftmp27 - 1)*M[6];
double Ftmp29 = Ftmp1*Ftmp20;
double Ftmp30 = 3*Ftmp29;
double Ftmp31 = (Ftmp30 - 1)*M[8];
double Ftmp32 = Ftmp13*Ftmp18;
double Ftmp33 = 6*y;
double Ftmp34 = Ftmp22*M[3];
double Ftmp35 = Ftmp2*y;
double Ftmp36 = Ftmp13*Ftmp20;
double Ftmp37 = 6*z;
#pragma omp atomic
F[0] += Ftmp0*(Ftmp1*Ftmp16*(Ftmp9 - 1)*M[3] - Ftmp10*M[0] + Ftmp11*Ftmp15 + Ftmp17*Ftmp19 + Ftmp17*Ftmp21 - Ftmp2*Ftmp3 + Ftmp23*Ftmp24 + Ftmp23*Ftmp3 + Ftmp25*Ftmp5 + Ftmp28*Ftmp5 + Ftmp31*Ftmp5 - Ftmp4*M[5] - Ftmp6*M[1] - Ftmp7*x + M[0]);
#pragma omp atomic
F[1] += Ftmp0*(Ftmp1*Ftmp33*(Ftmp26 - 1)*M[6] + Ftmp11*Ftmp32 + Ftmp12*Ftmp21*Ftmp33 + Ftmp15*Ftmp24 + Ftmp25*Ftmp35 - Ftmp27*M[1] + Ftmp28*Ftmp35 + Ftmp31*Ftmp35 + Ftmp32*x*M[4] + Ftmp33*Ftmp34 - Ftmp4*M[7] - Ftmp5*M[4] - Ftmp6*M[0] - Ftmp7*y + M[1]);
#pragma omp atomic
F[2] += Ftmp0*(Ftmp1*Ftmp37*(Ftmp29 - 1)*M[8] + Ftmp12*Ftmp19*Ftmp37 + Ftmp14*Ftmp3*z + Ftmp25*Ftmp4 + Ftmp28*Ftmp4 - Ftmp30*M[2] + Ftmp31*Ftmp4 + Ftmp34*Ftmp37 - Ftmp35*M[7] + Ftmp36*x*M[5] + Ftmp36*y*M[7] - Ftmp4*x*M[0] - Ftmp4*y*M[1] - Ftmp5*M[5] + M[2]);

}

void S2Mc_2(double x, double y, double z, double * S, double * M) {
double Mtmp0 = -z*S[2];
#pragma omp atomic
M[0] += S[0];
#pragma omp atomic
M[1] += S[1];
#pragma omp atomic
M[2] += S[2];
#pragma omp atomic
M[3] += Mtmp0 + x*S[0];
#pragma omp atomic
M[4] += x*S[1] + y*S[0];
#pragma omp atomic
M[5] += x*S[2] + z*S[0];
#pragma omp atomic
M[6] += Mtmp0 + y*S[1];
#pragma omp atomic
M[7] += y*S[2] + z*S[1];

}

void M2Mc_2(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = -z*M[2];
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += Mstmp0 + x*M[0] + M[3];
#pragma omp atomic
Ms[4] += x*M[1] + y*M[0] + M[4];
#pragma omp atomic
Ms[5] += x*M[2] + z*M[0] + M[5];
#pragma omp atomic
Ms[6] += Mstmp0 + y*M[1] + M[6];
#pragma omp atomic
Ms[7] += y*M[2] + z*M[1] + M[7];

}

void L2Lc_2(double x, double y, double z, double * L, double * Ls) {
#pragma omp atomic
Ls[0] += x*L[1] + y*L[2] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += L[1];
#pragma omp atomic
Ls[2] += L[2];
#pragma omp atomic
Ls[3] += L[3];

}

void L2Pc_2(double x, double y, double z, double * L, double * F) {
#pragma omp atomic
F[0] += -L[1];
#pragma omp atomic
F[1] += -L[2];
#pragma omp atomic
F[2] += -L[3];

}

void M2Pc_2(double x, double y, double z, double * M, double * F) {
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(Rinv, 3);
double Ftmp1 = pow(Rinv, 2);
double Ftmp2 = 3*Ftmp1;
double Ftmp3 = y*M[4];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = Ftmp2*x;
double Ftmp6 = Ftmp5*y;
double Ftmp7 = Ftmp4*M[2];
double Ftmp8 = pow(x, 2);
double Ftmp9 = Ftmp1*Ftmp8;
double Ftmp10 = 3*Ftmp9;
double Ftmp11 = z*M[7];
double Ftmp12 = pow(Rinv, 4);
double Ftmp13 = 15*Ftmp12*x;
double Ftmp14 = Ftmp13*y;
double Ftmp15 = 6*x;
double Ftmp16 = pow(y, 2);
double Ftmp17 = Ftmp12*Ftmp16;
double Ftmp18 = Ftmp17*M[6];
double Ftmp19 = Ftmp12*Ftmp8;
double Ftmp20 = 15*Ftmp19;
double Ftmp21 = z*M[5];
double Ftmp22 = (Ftmp10 - 1)*M[3];
double Ftmp23 = Ftmp1*Ftmp16;
double Ftmp24 = 3*Ftmp23;
double Ftmp25 = (Ftmp24 - 1)*M[6];
double Ftmp26 = 15*Ftmp17;
double Ftmp27 = 6*y;
double Ftmp28 = Ftmp19*M[3];
double Ftmp29 = Ftmp2*y;
double Ftmp30 = pow(z, 2);
double Ftmp31 = 15*Ftmp12*Ftmp30;
double Ftmp32 = 6*z;
#pragma omp atomic
F[0] += Ftmp0*(Ftmp1*Ftmp15*(Ftmp9 - 1)*M[3] - Ftmp10*M[0] + Ftmp11*Ftmp14 + Ftmp15*Ftmp18 - Ftmp2*Ftmp3 + Ftmp20*Ftmp21 + Ftmp20*Ftmp3 + Ftmp22*Ftmp5 + Ftmp25*Ftmp5 - Ftmp4*M[5] - Ftmp6*M[1] - Ftmp7*x + M[0]);
#pragma omp atomic
F[1] += Ftmp0*(Ftmp1*Ftmp27*(Ftmp23 - 1)*M[6] + Ftmp11*Ftmp26 + Ftmp14*Ftmp21 + Ftmp22*Ftmp29 - Ftmp24*M[1] + Ftmp25*Ftmp29 + Ftmp26*x*M[4] + Ftmp27*Ftmp28 - Ftmp4*M[7] - Ftmp5*M[4] - Ftmp6*M[0] - Ftmp7*y + M[1]);
#pragma omp atomic
F[2] += Ftmp0*(Ftmp13*Ftmp3*z + Ftmp18*Ftmp32 - Ftmp2*Ftmp30*M[2] + Ftmp22*Ftmp4 + Ftmp25*Ftmp4 + Ftmp28*Ftmp32 - Ftmp29*M[7] + Ftmp31*x*M[5] + Ftmp31*y*M[7] - Ftmp4*x*M[0] - Ftmp4*y*M[1] - Ftmp5*M[5] + M[2]);

}

void M2Lc_2(double x, double y, double z, double * M, double * L) {
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double D[9];
double Dtmp0 = pow(Rinv, 3);
double Dtmp1 = 3*pow(Rinv, 2);
double Dtmp2 = pow(Rinv, 5);
double Dtmp3 = 3*Dtmp2*x;
D[0] = -Dtmp0*x;
D[1] = -Dtmp0*y;
D[2] = -Dtmp0*z;
D[3] = Dtmp0*(Dtmp1*pow(x, 2) - 1);
D[4] = Dtmp3*y;
D[5] = Dtmp3*z;
D[6] = Dtmp0*(Dtmp1*pow(y, 2) - 1);
D[7] = 3*Dtmp2*y*z;
D[8] = -D[3] - D[6];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7];
#pragma omp atomic
L[1] += D[3]*M[0] + D[4]*M[1] + D[5]*M[2];
#pragma omp atomic
L[2] += D[4]*M[0] + D[6]*M[1] + D[7]*M[2];
#pragma omp atomic
L[3] += D[5]*M[0] + D[7]*M[1] + D[8]*M[2];

}

void P2P(double x, double y, double z, double * S, double * F) {
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(Rinv, 3);
double Ftmp1 = 3*pow(Rinv, 2);
double Ftmp2 = Ftmp1*S[1];
double Ftmp3 = x*y;
double Ftmp4 = Ftmp1*S[2];
double Ftmp5 = Ftmp4*z;
double Ftmp6 = Ftmp1*S[0];
#pragma omp atomic
F[0] += Ftmp0*(-Ftmp2*Ftmp3 - Ftmp5*x - Ftmp6*pow(x, 2) + S[0]);
#pragma omp atomic
F[1] += Ftmp0*(-Ftmp2*pow(y, 2) - Ftmp3*Ftmp6 - Ftmp5*y + S[1]);
#pragma omp atomic
F[2] += Ftmp0*(-Ftmp2*y*z - Ftmp4*pow(z, 2) - Ftmp6*x*z + S[2]);

}

void P2P_batch(double tx, double ty, double tz, const double * sx, const double * sy, const double * sz, const double * S, size_t begin, size_t end, double * F) {
double facc0 = 0.0;
double facc1 = 0.0;
double facc2 = 0.0;
#pragma omp simd reduction(+:facc0,facc1,facc2)
for (size_t u = begin; u < end; u++) {
double x = tx - sx[u];
double y = ty - sy[u];
double z = tz - sz[u];
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(Rinv, 3);
double Ftmp1 = 3*pow(Rinv, 2);
double Ftmp2 = Ftmp1*S[3*u + 1];
double Ftmp3 = x*y;
double Ftmp4 = Ftmp1*S[3*u + 2];
double Ftmp5 = Ftmp4*z;
double Ftmp6 = Ftmp1*S[3*u + 0];
facc0 += Ftmp0*(-Ftmp2*Ftmp3 - Ftmp5*x - Ftmp6*pow(x, 2) + S[3*u + 0]);
facc1 += Ftmp0*(-Ftmp2*pow(y, 2) - Ftmp3*Ftmp6 - Ftmp5*y + S[3*u + 1]);
facc2 += Ftmp0*(-Ftmp2*y*z - Ftmp4*pow(z, 2) - Ftmp6*x*z + S[3*u + 2]);
}
F[0] += facc0;
F[1] += facc1;
F[2] += facc2;
}

void S2M_3(double x, double y, double z, double * S, double * M) {
double Mtmp0 = x*S[1];
double Mtmp1 = y*S[0];
double Mtmp2 = x*S[2];
double Mtmp3 = z*S[0];
double Mtmp4 = y*S[2];
double Mtmp5 = z*S[1];
#pragma omp atomic
M[0] += S[0];
#pragma omp atomic
M[1] += S[1];
#pragma omp atomic
M[2] += S[2];
#pragma omp atomic
M[3] += x*S[0];
#pragma omp atomic
M[4] += Mtmp0 + Mtmp1;
#pragma omp atomic
M[5] += Mtmp2 + Mtmp3;
#pragma omp atomic
M[6] += y*S[1];
#pragma omp atomic
M[7] += Mtmp4 + Mtmp5;
#pragma omp atomic
M[8] += z*S[2];
#pragma omp atomic
M[9] += (1.0/2.0)*pow(x, 2)*S[0];
#pragma omp atomic
M[10] += x*((1.0/2.0)*Mtmp0 + Mtmp1);
#pragma omp atomic
M[11] += x*((1.0/2.0)*Mtmp2 + Mtmp3);
#pragma omp atomic
M[12] += y*(Mtmp0 + (1.0/2.0)*Mtmp1);
#pragma omp atomic
M[13] += Mtmp0*z + Mtmp1*z + Mtmp2*y;
#pragma omp atomic
M[14] += z*(Mtmp2 + (1.0/2.0)*Mtmp3);
#pragma omp atomic
M[15] += (1.0/2.0)*pow(y, 2)*S[1];
#pragma omp atomic
M[16] += y*((1.0/2.0)*Mtmp4 + Mtmp5);
#pragma omp atomic
M[17] += z*(Mtmp4 + (1.0/2.0)*Mtmp5);
#pragma omp atomic
M[18] += (1.0/2.0)*pow(z, 2)*S[2];

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
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double D[19];
double Dtmp0 = pow(Rinv, 3);
double Dtmp1 = pow(x, 2);
double Dtmp2 = pow(Rinv, 2);
double Dtmp3 = 3*Dtmp2;
double Dtmp4 = 3*pow(Rinv, 5);
double Dtmp5 = x*y;
double Dtmp6 = Dtmp4*z;
double Dtmp7 = pow(y, 2);
double Dtmp8 = 5*Dtmp2;
double Dtmp9 = Dtmp1*Dtmp8;
double Dtmp10 = Dtmp4*x;
double Dtmp11 = Dtmp9 - 1;
double Dtmp12 = Dtmp4*y;
double Dtmp13 = Dtmp7*Dtmp8;
double Dtmp14 = Dtmp13 - 1;
D[0] = -Dtmp0*x;
D[1] = -Dtmp0*y;
D[2] = -Dtmp0*z;
D[3] = Dtmp0*(Dtmp1*Dtmp3 - 1);
D[4] = Dtmp4*Dtmp5;
D[5] = Dtmp6*x;
D[6] = Dtmp0*(Dtmp3*Dtmp7 - 1);
D[7] = Dtmp6*y;
D[8] = -D[3] - D[6];
D[9] = -Dtmp10*(Dtmp9 - 3);
D[10] = -Dtmp11*Dtmp12;
D[11] = -Dtmp11*Dtmp6;
D[12] = -Dtmp10*Dtmp14;
D[13] = -15*Dtmp5*pow(Rinv, 7)*z;
D[14] = -D[9] - D[12];
D[15] = -Dtmp12*(Dtmp13 - 3);
D[16] = -Dtmp14*Dtmp6;
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
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(Rinv, 3);
double Ftmp1 = pow(Rinv, 2);
double Ftmp2 = 3*Ftmp1;
double Ftmp3 = Ftmp2*z;
double Ftmp4 = Ftmp2*x;
double Ftmp5 = Ftmp4*y;
double Ftmp6 = Ftmp3*M[2];
double Ftmp7 = pow(Rinv, 4);
double Ftmp8 = pow(x, 2);
double Ftmp9 = Ftmp1*Ftmp8;
double Ftmp10 = 3*Ftmp9;
double Ftmp11 = pow(y, 2);
double Ftmp12 = pow(z, 2);
double Ftmp13 = pow(Rinv, 6);
double Ftmp14 = 30*Ftmp13;
double Ftmp15 = Ftmp14*x;
double Ftmp16 = pow(y, 3)*M[15];
double Ftmp17 = pow(z, 3)*M[18];
double Ftmp18 = y*M[17];
double Ftmp19 = Ftmp12*Ftmp15;
double Ftmp20 = Ftmp11*M[16];
double Ftmp21 = Ftmp15*z;
double Ftmp22 = 105*Ftmp13*z*M[13];
double Ftmp23 = Ftmp14*Ftmp8;
double Ftmp24 = Ftmp11*M[12];
double Ftmp25 = Ftmp12*Ftmp23;
double Ftmp26 = 5*Ftmp9;
double Ftmp27 = Ftmp26 - 3;
double Ftmp28 = Ftmp1*Ftmp11;
double Ftmp29 = 5*Ftmp28;
double Ftmp30 = Ftmp29 - 1;
double Ftmp31 = Ftmp1*Ftmp12;
double Ftmp32 = 5*Ftmp31;
double Ftmp33 = Ftmp32 - 1;
double Ftmp34 = Ftmp9 - 1;
double Ftmp35 = Ftmp10 - 1;
double Ftmp36 = 3*Ftmp28;
double Ftmp37 = Ftmp36 - 1;
double Ftmp38 = 3*Ftmp31;
double Ftmp39 = Ftmp38 - 1;
double Ftmp40 = x*y;
double Ftmp41 = 30*Ftmp34;
double Ftmp42 = Ftmp41*Ftmp7;
double Ftmp43 = x*z;
double Ftmp44 = Ftmp7*Ftmp8;
double Ftmp45 = 15*Ftmp7;
double Ftmp46 = Ftmp45*y;
double Ftmp47 = Ftmp46*x;
double Ftmp48 = Ftmp26 - 1;
double Ftmp49 = Ftmp48*M[10];
double Ftmp50 = Ftmp29 - 3;
double Ftmp51 = Ftmp50*M[15];
double Ftmp52 = Ftmp33*M[17];
double Ftmp53 = Ftmp45*Ftmp48;
double Ftmp54 = Ftmp53*M[11];
double Ftmp55 = Ftmp30*M[16];
double Ftmp56 = Ftmp45*x*z;
double Ftmp57 = Ftmp32 - 3;
double Ftmp58 = Ftmp57*M[18];
double Ftmp59 = Ftmp27*M[9];
double Ftmp60 = 15*Ftmp44;
double Ftmp61 = Ftmp30*M[12];
double Ftmp62 = Ftmp33*M[14];
double Ftmp63 = Ftmp14*y;
double Ftmp64 = pow(x, 3)*M[9];
double Ftmp65 = y*z;
double Ftmp66 = Ftmp23*Ftmp65;
double Ftmp67 = Ftmp11*M[10];
double Ftmp68 = Ftmp11*M[17];
double Ftmp69 = Ftmp12*Ftmp14;
double Ftmp70 = Ftmp28 - 1;
double Ftmp71 = 30*Ftmp7;
double Ftmp72 = Ftmp70*Ftmp71;
double Ftmp73 = Ftmp46*z;
double Ftmp74 = Ftmp14*z;
double Ftmp75 = Ftmp31 - 1;
double Ftmp76 = Ftmp71*Ftmp75;
double Ftmp77 = Ftmp12*Ftmp45;
#pragma omp atomic
F[0] += Ftmp0*(3*Ftmp1*Ftmp27*M[9] + 3*Ftmp1*Ftmp30*M[12] + 3*Ftmp1*Ftmp33*M[14] + 6*Ftmp1*Ftmp34*x*M[3] + 3*Ftmp1*Ftmp35*x*M[3] + 3*Ftmp1*Ftmp37*x*M[6] + 3*Ftmp1*Ftmp39*x*M[8] - Ftmp10*M[0] + 6*Ftmp11*Ftmp7*x*M[6] + 6*Ftmp12*Ftmp7*x*M[8] - Ftmp15*Ftmp16 - Ftmp15*Ftmp17 - Ftmp18*Ftmp19 - Ftmp2*y*M[4] - Ftmp20*Ftmp21 - Ftmp22*Ftmp8*y - Ftmp23*Ftmp24 - Ftmp25*M[14] - Ftmp3*M[5] - Ftmp40*Ftmp42*M[10] - Ftmp41*Ftmp44*M[9] - Ftmp42*Ftmp43*M[11] - Ftmp43*Ftmp54 - Ftmp47*Ftmp49 - Ftmp47*Ftmp51 - Ftmp47*Ftmp52 - Ftmp5*M[1] - Ftmp55*Ftmp56 - Ftmp56*Ftmp58 - Ftmp59*Ftmp60 - Ftmp6*x - Ftmp60*Ftmp61 - Ftmp60*Ftmp62 + 15*Ftmp7*Ftmp8*y*M[4] + 15*Ftmp7*Ftmp8*z*M[5] + 15*Ftmp7*x*y*z*M[7] + 15*Ftmp7*y*z*M[13] + M[0]);
#pragma omp atomic
F[1] += Ftmp0*(3*Ftmp1*Ftmp33*M[17] + 3*Ftmp1*Ftmp35*y*M[3] + 3*Ftmp1*Ftmp37*y*M[6] + 3*Ftmp1*Ftmp39*y*M[8] + 3*Ftmp1*Ftmp48*M[10] + 3*Ftmp1*Ftmp50*M[15] + 6*Ftmp1*Ftmp70*y*M[6] - Ftmp11*Ftmp22*x - Ftmp11*Ftmp45*Ftmp51 + 15*Ftmp11*Ftmp7*x*M[4] + 15*Ftmp11*Ftmp7*z*M[7] - Ftmp11*Ftmp72*M[15] + 6*Ftmp12*Ftmp7*y*M[8] - Ftmp17*Ftmp63 - Ftmp19*y*M[14] - Ftmp23*Ftmp67 - Ftmp3*M[7] - Ftmp33*Ftmp45*Ftmp68 - Ftmp36*M[1] - Ftmp4*M[4] - Ftmp40*Ftmp72*M[12] - Ftmp47*Ftmp59 - Ftmp47*Ftmp61 - Ftmp47*Ftmp62 - Ftmp48*Ftmp73*M[11] - Ftmp5*M[0] - Ftmp53*Ftmp67 - Ftmp55*Ftmp73 - Ftmp58*Ftmp73 - Ftmp6*y - Ftmp63*Ftmp64 - Ftmp65*Ftmp72*M[16] - Ftmp66*M[11] - Ftmp68*Ftmp69 + 6*Ftmp7*Ftmp8*y*M[3] + 15*Ftmp7*x*y*z*M[5] + 15*Ftmp7*x*z*M[13] + M[1]);
#pragma omp atomic
F[2] += Ftmp0*(3*Ftmp1*Ftmp30*M[16] + 3*Ftmp1*Ftmp35*z*M[3] + 3*Ftmp1*Ftmp37*z*M[6] + 3*Ftmp1*Ftmp39*z*M[8] + 3*Ftmp1*Ftmp48*M[11] + 3*Ftmp1*Ftmp57*M[18] + 6*Ftmp1*Ftmp75*z*M[8] + 6*Ftmp11*Ftmp7*z*M[6] - 105*Ftmp12*Ftmp13*Ftmp40*M[13] - Ftmp12*Ftmp54 + 15*Ftmp12*Ftmp7*x*M[5] + 15*Ftmp12*Ftmp7*y*M[7] - Ftmp12*Ftmp76*M[18] - Ftmp16*Ftmp74 - Ftmp18*Ftmp76*z - Ftmp2*y*M[7] - Ftmp20*Ftmp69 - Ftmp21*Ftmp24 - Ftmp25*M[11] - Ftmp3*x*M[0] - Ftmp3*y*M[1] - Ftmp38*M[2] - Ftmp4*M[5] - Ftmp43*Ftmp76*M[14] - Ftmp49*Ftmp73 - Ftmp51*Ftmp73 - Ftmp52*Ftmp73 - Ftmp55*Ftmp77 - Ftmp56*Ftmp59 - Ftmp56*Ftmp61 - Ftmp56*Ftmp62 - Ftmp58*Ftmp77 - Ftmp64*Ftmp74 - Ftmp66*M[10] + 6*Ftmp7*Ftmp8*z*M[3] + 15*Ftmp7*x*y*z*M[4] + 15*Ftmp7*x*y*M[13] + M[2]);

}

void S2Mc_3(double x, double y, double z, double * S, double * M) {
double Mtmp0 = x*S[0];
double Mtmp1 = z*S[2];
double Mtmp2 = -Mtmp1;
double Mtmp3 = x*S[1];
double Mtmp4 = y*S[0];
double Mtmp5 = x*S[2];
double Mtmp6 = y*S[1];
double Mtmp7 = Mtmp1*x;
double Mtmp8 = pow(x, 2);
double Mtmp9 = pow(z, 2);
double Mtmp10 = (1.0/2.0)*S[0];
double Mtmp11 = Mtmp10*Mtmp9;
double Mtmp12 = Mtmp1*y;
double Mtmp13 = (1.0/2.0)*S[1];
double Mtmp14 = Mtmp13*Mtmp9;
double Mtmp15 = (1.0/2.0)*S[2];
double Mtmp16 = -Mtmp15*Mtmp9;
double Mtmp17 = pow(y, 2);
#pragma omp atomic
M[0] += S[0];
#pragma omp atomic
M[1] += S[1];
#pragma omp atomic
M[2] += S[2];
#pragma omp atomic
M[3] += Mtmp0 + Mtmp2;
#pragma omp atomic
M[4] += Mtmp3 + Mtmp4;
#pragma omp atomic
M[5] += Mtmp5 + z*S[0];
#pragma omp atomic
M[6] += Mtmp2 + Mtmp6;
#pragma omp atomic
M[7] += y*S[2] + z*S[1];
#pragma omp atomic
M[8] += -Mtmp11 - Mtmp7 + (1.0/2.0)*Mtmp8*S[0];
#pragma omp atomic
M[9] += Mtmp0*y - Mtmp12 + Mtmp13*Mtmp8 - Mtmp14;
#pragma omp atomic
M[10] += Mtmp0*z + Mtmp15*Mtmp8 + Mtmp16;
#pragma omp atomic
M[11] += Mtmp10*Mtmp17 - Mtmp11 + Mtmp3*y - Mtmp7;
#pragma omp atomic
M[12] += Mtmp3*z + Mtmp4*z + Mtmp5*y;
#pragma omp atomic
M[13] += -Mtmp12 - Mtmp14 + (1.0/2.0)*Mtmp17*S[1];
#pragma omp atomic
M[14] += Mtmp15*Mtmp17 + Mtmp16 + Mtmp6*z;

}

void M2Mc_3(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = z*M[2];
double Mstmp2 = -Mstmp1;
double Mstmp3 = x*M[1];
double Mstmp4 = y*M[0];
double Mstmp5 = x*M[2];
double Mstmp6 = y*M[1];
double Mstmp7 = pow(x, 2);
double Mstmp8 = (1.0/2.0)*M[0];
double Mstmp9 = pow(z, 2);
double Mstmp10 = -Mstmp1*x - Mstmp8*Mstmp9 - z*M[5];
double Mstmp11 = (1.0/2.0)*M[1];
double Mstmp12 = -Mstmp1*y - Mstmp11*Mstmp9 - z*M[7];
double Mstmp13 = (1.0/2.0)*M[2];
double Mstmp14 = -Mstmp13*Mstmp9;
double Mstmp15 = pow(y, 2);
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += Mstmp0 + Mstmp2 + M[3];
#pragma omp atomic
Ms[4] += Mstmp3 + Mstmp4 + M[4];
#pragma omp atomic
Ms[5] += Mstmp5 + z*M[0] + M[5];
#pragma omp atomic
Ms[6] += Mstmp2 + Mstmp6 + M[6];
#pragma omp atomic
Ms[7] += y*M[2] + z*M[1] + M[7];
#pragma omp atomic
Ms[8] += Mstmp10 + Mstmp7*Mstmp8 + x*M[3] + M[8];
#pragma omp atomic
Ms[9] += Mstmp0*y + Mstmp11*Mstmp7 + Mstmp12 + x*M[4] + y*M[3] + M[9];
#pragma omp atomic
Ms[10] += Mstmp0*z + Mstmp13*Mstmp7 + Mstmp14 + x*M[5] + z*M[3] + M[10];
#pragma omp atomic
Ms[11] += Mstmp10 + Mstmp15*Mstmp8 + Mstmp3*y + x*M[6] + y*M[4] + M[11];
#pragma omp atomic
Ms[12] += Mstmp3*z + Mstmp4*z + Mstmp5*y + x*M[7] + y*M[5] + z*M[4] + M[12];
#pragma omp atomic
Ms[13] += Mstmp11*Mstmp15 + Mstmp12 + y*M[6] + M[13];
#pragma omp atomic
Ms[14] += Mstmp13*Mstmp15 + Mstmp14 + Mstmp6*z + y*M[7] + z*M[6] + M[14];

}

void L2Lc_3(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = L[4] + L[7];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp2*y - 1.0/2.0*Lstmp3*pow(z, 2) + (1.0/2.0)*pow(x, 2)*L[4] + x*L[1] + (1.0/2.0)*pow(y, 2)*L[7] + y*L[2] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp2 + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += -Lstmp3*z + x*L[6] + y*L[8] + L[3];
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

}

void L2Pc_3(double x, double y, double z, double * L, double * F) {
#pragma omp atomic
F[0] += -x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -x*L[6] - y*L[8] + z*(L[4] + L[7]) - L[3];

}

void M2Pc_3(double x, double y, double z, double * M, double * F) {
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(Rinv, 3);
double Ftmp1 = pow(Rinv, 2);
double Ftmp2 = 3*Ftmp1;
double Ftmp3 = Ftmp2*z;
double Ftmp4 = Ftmp2*x;
double Ftmp5 = Ftmp4*y;
double Ftmp6 = Ftmp3*M[2];
double Ftmp7 = pow(Rinv, 4);
double Ftmp8 = pow(x, 2);
double Ftmp9 = Ftmp1*Ftmp8;
double Ftmp10 = 3*Ftmp9;
double Ftmp11 = pow(y, 2);
double Ftmp12 = pow(Rinv, 6);
double Ftmp13 = 30*Ftmp12;
double Ftmp14 = Ftmp13*x;
double Ftmp15 = pow(y, 3)*M[13];
double Ftmp16 = Ftmp11*z;
double Ftmp17 = y*z;
double Ftmp18 = 105*Ftmp12;
double Ftmp19 = Ftmp13*Ftmp8;
double Ftmp20 = Ftmp11*M[11];
double Ftmp21 = 5*Ftmp9;
double Ftmp22 = Ftmp21 - 3;
double Ftmp23 = Ftmp1*Ftmp11;
double Ftmp24 = 5*Ftmp23;
double Ftmp25 = Ftmp24 - 1;
double Ftmp26 = Ftmp9 - 1;
double Ftmp27 = Ftmp10 - 1;
double Ftmp28 = 3*Ftmp23;
double Ftmp29 = Ftmp28 - 1;
double Ftmp30 = 30*Ftmp26;
double Ftmp31 = Ftmp7*y;
double Ftmp32 = Ftmp31*x;
double Ftmp33 = Ftmp7*x;
double Ftmp34 = Ftmp33*M[10];
double Ftmp35 = Ftmp7*Ftmp8;
double Ftmp36 = Ftmp21 - 1;
double Ftmp37 = Ftmp36*M[9];
double Ftmp38 = 15*Ftmp31*x;
double Ftmp39 = Ftmp24 - 3;
double Ftmp40 = Ftmp39*M[13];
double Ftmp41 = 15*z;
double Ftmp42 = Ftmp36*Ftmp41;
double Ftmp43 = Ftmp33*Ftmp41;
double Ftmp44 = Ftmp25*M[14];
double Ftmp45 = Ftmp22*M[8];
double Ftmp46 = 15*Ftmp35;
double Ftmp47 = Ftmp25*M[11];
double Ftmp48 = Ftmp13*pow(x, 3)*M[8];
double Ftmp49 = Ftmp18*x*M[12];
double Ftmp50 = Ftmp19*M[9];
double Ftmp51 = Ftmp23 - 1;
double Ftmp52 = 30*Ftmp51;
double Ftmp53 = Ftmp11*Ftmp7;
double Ftmp54 = Ftmp31*Ftmp41;
double Ftmp55 = 15*Ftmp53;
double Ftmp56 = pow(z, 2);
double Ftmp57 = Ftmp56*M[10];
double Ftmp58 = Ftmp56*M[14];
double Ftmp59 = 15*Ftmp7;
#pragma omp atomic
F[0] += Ftmp0*(3*Ftmp1*Ftmp22*M[8] + 3*Ftmp1*Ftmp25*M[11] + 6*Ftmp1*Ftmp26*x*M[3] + 3*Ftmp1*Ftmp27*x*M[3] + 3*Ftmp1*Ftmp29*x*M[6] - Ftmp10*M[0] + 6*Ftmp11*Ftmp7*x*M[6] - Ftmp14*Ftmp15 - Ftmp14*Ftmp16*M[14] - Ftmp17*Ftmp18*Ftmp8*M[12] - Ftmp19*Ftmp20 - Ftmp2*y*M[4] - Ftmp3*M[5] - Ftmp30*Ftmp32*M[9] - Ftmp30*Ftmp34*z - Ftmp30*Ftmp35*M[8] - Ftmp34*Ftmp42 - Ftmp37*Ftmp38 - Ftmp38*Ftmp40 - Ftmp43*Ftmp44 - Ftmp45*Ftmp46 - Ftmp46*Ftmp47 - Ftmp5*M[1] - Ftmp6*x + 15*Ftmp7*Ftmp8*y*M[4] + 15*Ftmp7*Ftmp8*z*M[5] + 15*Ftmp7*x*y*z*M[7] + 15*Ftmp7*y*z*M[12] + M[0]);
#pragma omp atomic
F[1] += Ftmp0*(3*Ftmp1*Ftmp27*y*M[3] + 3*Ftmp1*Ftmp29*y*M[6] + 3*Ftmp1*Ftmp36*M[9] + 3*Ftmp1*Ftmp39*M[13] + 6*Ftmp1*Ftmp51*y*M[6] - Ftmp11*Ftmp50 + 15*Ftmp11*Ftmp7*x*M[4] + 15*Ftmp11*Ftmp7*z*M[7] - Ftmp16*Ftmp49 - Ftmp17*Ftmp19*M[10] - Ftmp28*M[1] - Ftmp3*M[7] - Ftmp31*Ftmp42*M[10] - Ftmp31*Ftmp52*z*M[14] - Ftmp32*Ftmp52*M[11] - Ftmp37*Ftmp55 - Ftmp38*Ftmp45 - Ftmp38*Ftmp47 - Ftmp4*M[4] - Ftmp40*Ftmp55 - Ftmp44*Ftmp54 - Ftmp48*y - Ftmp5*M[0] - Ftmp52*Ftmp53*M[13] - Ftmp6*y + 6*Ftmp7*Ftmp8*y*M[3] + 15*Ftmp7*x*y*z*M[5] + 15*Ftmp7*x*z*M[12] + M[1]);
#pragma omp atomic
F[2] += Ftmp0*(3*Ftmp1*Ftmp25*M[14] + 3*Ftmp1*Ftmp27*z*M[3] + 3*Ftmp1*Ftmp29*z*M[6] + 3*Ftmp1*Ftmp36*M[10] - Ftmp11*Ftmp13*Ftmp58 + 6*Ftmp11*Ftmp7*z*M[6] - Ftmp13*Ftmp15*z - Ftmp14*Ftmp20*z - Ftmp17*Ftmp50 - Ftmp19*Ftmp57 - Ftmp2*Ftmp56*M[2] - Ftmp2*y*M[7] - Ftmp25*Ftmp58*Ftmp59 - Ftmp3*x*M[0] - Ftmp3*y*M[1] - Ftmp36*Ftmp57*Ftmp59 - Ftmp37*Ftmp54 - Ftmp4*M[5] - Ftmp40*Ftmp54 - Ftmp43*Ftmp45 - Ftmp43*Ftmp47 - Ftmp48*z - Ftmp49*Ftmp56*y + 15*Ftmp56*Ftmp7*x*M[5] + 15*Ftmp56*Ftmp7*y*M[7] + 6*Ftmp7*Ftmp8*z*M[3] + 15*Ftmp7*x*y*z*M[4] + 15*Ftmp7*x*y*M[12] + M[2]);

}

void M2Lc_3(double x, double y, double z, double * M, double * L) {
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double D[18];
double Dtmp0 = pow(Rinv, 3);
double Dtmp1 = pow(x, 2);
double Dtmp2 = pow(Rinv, 2);
double Dtmp3 = 3*Dtmp2;
double Dtmp4 = 3*pow(Rinv, 5);
double Dtmp5 = x*y;
double Dtmp6 = Dtmp4*z;
double Dtmp7 = pow(y, 2);
double Dtmp8 = 5*Dtmp2;
double Dtmp9 = Dtmp1*Dtmp8;
double Dtmp10 = Dtmp4*x;
double Dtmp11 = Dtmp9 - 1;
double Dtmp12 = Dtmp4*y;
double Dtmp13 = Dtmp7*Dtmp8;
double Dtmp14 = Dtmp13 - 1;
D[0] = -Dtmp0*x;
D[1] = -Dtmp0*y;
D[2] = -Dtmp0*z;
D[3] = Dtmp0*(Dtmp1*Dtmp3 - 1);
D[4] = Dtmp4*Dtmp5;
D[5] = Dtmp6*x;
D[6] = Dtmp0*(Dtmp3*Dtmp7 - 1);
D[7] = Dtmp6*y;
D[8] = -D[3] - D[6];
D[9] = -Dtmp10*(Dtmp9 - 3);
D[10] = -Dtmp11*Dtmp12;
D[11] = -Dtmp11*Dtmp6;
D[12] = -Dtmp10*Dtmp14;
D[13] = -15*Dtmp5*pow(Rinv, 7)*z;
D[14] = -D[9] - D[12];
D[15] = -Dtmp12*(Dtmp13 - 3);
D[16] = -Dtmp14*Dtmp6;
D[17] = -D[10] - D[15];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[9]*M[8] + D[10]*M[9] + D[11]*M[10] + D[12]*M[11] + D[13]*M[12] + D[15]*M[13] + D[16]*M[14];
#pragma omp atomic
L[1] += D[3]*M[0] + D[4]*M[1] + D[5]*M[2] + D[9]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7];
#pragma omp atomic
L[2] += D[4]*M[0] + D[6]*M[1] + D[7]*M[2] + D[10]*M[3] + D[12]*M[4] + D[13]*M[5] + D[15]*M[6] + D[16]*M[7];
#pragma omp atomic
L[3] += D[5]*M[0] + D[7]*M[1] + D[8]*M[2] + D[11]*M[3] + D[13]*M[4] + D[14]*M[5] + D[16]*M[6] + D[17]*M[7];
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

}

void S2M_4(double x, double y, double z, double * S, double * M) {
double Mtmp0 = x*S[1];
double Mtmp1 = y*S[0];
double Mtmp2 = Mtmp0 + Mtmp1;
double Mtmp3 = x*S[2];
double Mtmp4 = z*S[0];
double Mtmp5 = Mtmp3 + Mtmp4;
double Mtmp6 = y*S[2];
double Mtmp7 = z*S[1];
double Mtmp8 = Mtmp6 + Mtmp7;
double Mtmp9 = pow(x, 2);
double Mtmp10 = (1.0/2.0)*Mtmp0;
double Mtmp11 = (1.0/2.0)*Mtmp3;
double Mtmp12 = (1.0/2.0)*Mtmp1;
double Mtmp13 = Mtmp3*y;
double Mtmp14 = Mtmp0*z;
double Mtmp15 = Mtmp1*z;
double Mtmp16 = pow(y, 2);
double Mtmp17 = pow(z, 2);
double Mtmp18 = (1.0/6.0)*Mtmp9;
double Mtmp19 = (1.0/2.0)*x;
double Mtmp20 = Mtmp11*y;
double Mtmp21 = Mtmp10*z;
double Mtmp22 = (1.0/6.0)*Mtmp16;
double Mtmp23 = Mtmp12*z;
double Mtmp24 = (1.0/6.0)*Mtmp17;
#pragma omp atomic
M[0] += S[0];
#pragma omp atomic
M[1] += S[1];
#pragma omp atomic
M[2] += S[2];
#pragma omp atomic
M[3] += x*S[0];
#pragma omp atomic
M[4] += Mtmp2;
#pragma omp atomic
M[5] += Mtmp5;
#pragma omp atomic
M[6] += y*S[1];
#pragma omp atomic
M[7] += Mtmp8;
#pragma omp atomic
M[8] += z*S[2];
#pragma omp atomic
M[9] += (1.0/2.0)*Mtmp9*S[0];
#pragma omp atomic
M[10] += x*(Mtmp1 + Mtmp10);
#pragma omp atomic
M[11] += x*(Mtmp11 + Mtmp4);
#pragma omp atomic
M[12] += y*(Mtmp0 + Mtmp12);
#pragma omp atomic
M[13] += Mtmp13 + Mtmp14 + Mtmp15;
#pragma omp atomic
M[14] += z*(Mtmp3 + (1.0/2.0)*Mtmp4);
#pragma omp atomic
M[15] += (1.0/2.0)*Mtmp16*S[1];
#pragma omp atomic
M[16] += y*((1.0/2.0)*Mtmp6 + Mtmp7);
#pragma omp atomic
M[17] += z*(Mtmp6 + (1.0/2.0)*Mtmp7);
#pragma omp atomic
M[18] += (1.0/2.0)*Mtmp17*S[2];
#pragma omp atomic
M[19] += (1.0/6.0)*pow(x, 3)*S[0];
#pragma omp atomic
M[20] += Mtmp18*(Mtmp0 + 3*Mtmp1);
#pragma omp atomic
M[21] += Mtmp18*(Mtmp3 + 3*Mtmp4);
#pragma omp atomic
M[22] += Mtmp19*Mtmp2*y;
#pragma omp atomic
M[23] += x*(Mtmp15 + Mtmp20 + Mtmp21);
#pragma omp atomic
M[24] += Mtmp19*Mtmp5*z;
#pragma omp atomic
M[25] += Mtmp22*(3*Mtmp0 + Mtmp1);
#pragma omp atomic
M[26] += y*(Mtmp14 + Mtmp20 + Mtmp23);
#pragma omp atomic
M[27] += z*(Mtmp13 + Mtmp21 + Mtmp23);
#pragma omp atomic
M[28] += Mtmp24*(3*Mtmp3 + Mtmp4);
#pragma omp atomic
M[29] += (1.0/6.0)*pow(y, 3)*S[1];
#pragma omp atomic
M[30] += Mtmp22*(Mtmp6 + 3*Mtmp7);
#pragma omp atomic
M[31] += (1.0/2.0)*Mtmp8*y*z;
#pragma omp atomic
M[32] += Mtmp24*(3*Mtmp6 + Mtmp7);
#pragma omp atomic
M[33] += (1.0/6.0)*pow(z, 3)*S[2];

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
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double D[34];
double Dtmp0 = pow(Rinv, 3);
double Dtmp1 = pow(x, 2);
double Dtmp2 = pow(Rinv, 2);
double Dtmp3 = 3*Dtmp2;
double Dtmp4 = 3*pow(Rinv, 5);
double Dtmp5 = Dtmp4*x;
double Dtmp6 = pow(y, 2);
double Dtmp7 = Dtmp4*y;
double Dtmp8 = 5*Dtmp2;
double Dtmp9 = Dtmp1*Dtmp8;
double Dtmp10 = Dtmp9 - 1;
double Dtmp11 = Dtmp4*z;
double Dtmp12 = Dtmp6*Dtmp8;
double Dtmp13 = Dtmp12 - 1;
double Dtmp14 = pow(Rinv, 7);
double Dtmp15 = 15*Dtmp14*x*y;
double Dtmp16 = Dtmp1*Dtmp2;
double Dtmp17 = 35*pow(Rinv, 4);
double Dtmp18 = 7*Dtmp16;
double Dtmp19 = Dtmp18 - 3;
double Dtmp20 = 15*Dtmp14*z;
double Dtmp21 = Dtmp20*x;
double Dtmp22 = Dtmp20*y;
double Dtmp23 = Dtmp2*Dtmp6;
double Dtmp24 = 7*Dtmp23;
double Dtmp25 = Dtmp24 - 3;
D[0] = -Dtmp0*x;
D[1] = -Dtmp0*y;
D[2] = -Dtmp0*z;
D[3] = Dtmp0*(Dtmp1*Dtmp3 - 1);
D[4] = Dtmp5*y;
D[5] = Dtmp5*z;
D[6] = Dtmp0*(Dtmp3*Dtmp6 - 1);
D[7] = Dtmp7*z;
D[8] = -D[3] - D[6];
D[9] = -Dtmp5*(Dtmp9 - 3);
D[10] = -Dtmp10*Dtmp7;
D[11] = -Dtmp10*Dtmp11;
D[12] = -Dtmp13*Dtmp5;
D[13] = -Dtmp15*z;
D[14] = -D[9] - D[12];
D[15] = -Dtmp7*(Dtmp12 - 3);
D[16] = -Dtmp11*Dtmp13;
D[17] = -D[10] - D[15];
D[18] = -D[11] - D[16];
D[19] = Dtmp4*(-30*Dtmp16 + Dtmp17*pow(x, 4) + 3);
D[20] = Dtmp15*Dtmp19;
D[21] = Dtmp19*Dtmp21;
D[22] = Dtmp4*(Dtmp1*Dtmp17*Dtmp6 - Dtmp12 - Dtmp9 + 1);
D[23] = Dtmp22*(Dtmp18 - 1);
D[24] = -D[19] - D[22];
D[25] = Dtmp15*Dtmp25;
D[26] = Dtmp21*(Dtmp24 - 1);
D[27] = -D[20] - D[25];
D[28] = -D[21] - D[26];
D[29] = Dtmp4*(Dtmp17*pow(y, 4) - 30*Dtmp23 + 3);
D[30] = Dtmp22*Dtmp25;
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
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(Rinv, 3);
double Ftmp1 = pow(Rinv, 2);
double Ftmp2 = 3*Ftmp1;
double Ftmp3 = y*M[4];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = Ftmp2*x;
double Ftmp6 = Ftmp5*y;
double Ftmp7 = Ftmp4*M[2];
double Ftmp8 = pow(Rinv, 4);
double Ftmp9 = 15*Ftmp8;
double Ftmp10 = Ftmp9*y;
double Ftmp11 = Ftmp10*M[13];
double Ftmp12 = pow(x, 2);
double Ftmp13 = Ftmp1*Ftmp12;
double Ftmp14 = 3*Ftmp13;
double Ftmp15 = z*M[7];
double Ftmp16 = Ftmp10*x;
double Ftmp17 = 6*x;
double Ftmp18 = pow(y, 2);
double Ftmp19 = Ftmp18*Ftmp8;
double Ftmp20 = Ftmp19*M[6];
double Ftmp21 = pow(z, 2);
double Ftmp22 = Ftmp21*Ftmp8;
double Ftmp23 = Ftmp22*M[8];
double Ftmp24 = Ftmp12*Ftmp8;
double Ftmp25 = 15*Ftmp24;
double Ftmp26 = z*M[5];
double Ftmp27 = pow(y, 3);
double Ftmp28 = Ftmp27*x;
double Ftmp29 = pow(Rinv, 6);
double Ftmp30 = 30*Ftmp29;
double Ftmp31 = Ftmp30*M[15];
double Ftmp32 = Ftmp30*x;
double Ftmp33 = pow(z, 3);
double Ftmp34 = Ftmp33*M[18];
double Ftmp35 = Ftmp21*y;
double Ftmp36 = Ftmp32*Ftmp35;
double Ftmp37 = Ftmp18*z;
double Ftmp38 = y*z;
double Ftmp39 = 105*Ftmp29;
double Ftmp40 = Ftmp12*Ftmp39;
double Ftmp41 = Ftmp33*M[32];
double Ftmp42 = 210*pow(Rinv, 8);
double Ftmp43 = Ftmp42*x;
double Ftmp44 = Ftmp43*y;
double Ftmp45 = Ftmp28*Ftmp42*z;
double Ftmp46 = Ftmp12*Ftmp30;
double Ftmp47 = Ftmp18*M[12];
double Ftmp48 = Ftmp21*Ftmp46;
double Ftmp49 = Ftmp12*Ftmp42;
double Ftmp50 = Ftmp33*M[28];
double Ftmp51 = y*M[27];
double Ftmp52 = Ftmp37*Ftmp49;
double Ftmp53 = 5*Ftmp13;
double Ftmp54 = (Ftmp53 - 3)*M[9];
double Ftmp55 = Ftmp1*Ftmp18;
double Ftmp56 = 5*Ftmp55;
double Ftmp57 = Ftmp56 - 1;
double Ftmp58 = Ftmp2*Ftmp57;
double Ftmp59 = Ftmp1*Ftmp21;
double Ftmp60 = 5*Ftmp59;
double Ftmp61 = Ftmp60 - 1;
double Ftmp62 = Ftmp2*Ftmp61;
double Ftmp63 = Ftmp13 - 1;
double Ftmp64 = (Ftmp14 - 1)*M[3];
double Ftmp65 = 3*Ftmp55;
double Ftmp66 = (Ftmp65 - 1)*M[6];
double Ftmp67 = 3*Ftmp59;
double Ftmp68 = (Ftmp67 - 1)*M[8];
double Ftmp69 = 7*Ftmp13;
double Ftmp70 = Ftmp69 - 3;
double Ftmp71 = Ftmp70*M[20];
double Ftmp72 = 7*Ftmp55;
double Ftmp73 = Ftmp72 - 3;
double Ftmp74 = Ftmp73*M[25];
double Ftmp75 = 7*Ftmp59;
double Ftmp76 = Ftmp75 - 1;
double Ftmp77 = Ftmp51*Ftmp76;
double Ftmp78 = Ftmp9*z;
double Ftmp79 = Ftmp70*M[21];
double Ftmp80 = (Ftmp72 - 1)*M[26];
double Ftmp81 = Ftmp75 - 3;
double Ftmp82 = Ftmp78*Ftmp81;
double Ftmp83 = 30*Ftmp63;
double Ftmp84 = Ftmp8*x;
double Ftmp85 = Ftmp83*Ftmp84;
double Ftmp86 = Ftmp53 - 1;
double Ftmp87 = Ftmp86*M[10];
double Ftmp88 = (Ftmp56 - 3)*M[15];
double Ftmp89 = Ftmp61*M[17];
double Ftmp90 = Ftmp78*x;
double Ftmp91 = Ftmp86*M[11];
double Ftmp92 = Ftmp57*M[16];
double Ftmp93 = (Ftmp60 - 3)*M[18];
double Ftmp94 = Ftmp29*x;
double Ftmp95 = Ftmp38*Ftmp94;
double Ftmp96 = Ftmp57*M[12];
double Ftmp97 = Ftmp61*M[14];
double Ftmp98 = Ftmp12*y;
double Ftmp99 = 210*Ftmp29;
double Ftmp100 = Ftmp63*Ftmp99;
double Ftmp101 = Ftmp12*z;
double Ftmp102 = (Ftmp69 - 1)*M[23];
double Ftmp103 = Ftmp39*x;
double Ftmp104 = Ftmp103*Ftmp38;
double Ftmp105 = Ftmp73*M[30];
double Ftmp106 = Ftmp81*M[32];
double Ftmp107 = Ftmp18*Ftmp94;
double Ftmp108 = 60*M[29];
double Ftmp109 = Ftmp108*Ftmp73;
double Ftmp110 = 60*M[33];
double Ftmp111 = Ftmp110*Ftmp81;
double Ftmp112 = Ftmp21*Ftmp94;
double Ftmp113 = Ftmp40*y;
double Ftmp114 = Ftmp40*z;
double Ftmp115 = Ftmp81*M[28];
double Ftmp116 = pow(x, 4);
double Ftmp117 = 7*Ftmp8;
double Ftmp118 = 60*M[19];
double Ftmp119 = Ftmp9*x;
double Ftmp120 = 35*Ftmp8;
double Ftmp121 = (Ftmp116*Ftmp120 - 30*Ftmp13 + 3)*M[19];
double Ftmp122 = pow(y, 4);
double Ftmp123 = (Ftmp120*Ftmp122 - 30*Ftmp55 + 3)*M[29];
double Ftmp124 = pow(z, 4);
double Ftmp125 = (Ftmp120*Ftmp124 - 30*Ftmp59 + 3)*M[33];
double Ftmp126 = 14*Ftmp21;
double Ftmp127 = -8*Ftmp55;
double Ftmp128 = 14*Ftmp18;
double Ftmp129 = Ftmp128*Ftmp24;
double Ftmp130 = 1 - Ftmp13;
double Ftmp131 = 30*Ftmp84;
double Ftmp132 = -Ftmp56;
double Ftmp133 = Ftmp12*Ftmp120;
double Ftmp134 = 1 - Ftmp53;
double Ftmp135 = (Ftmp132 + Ftmp133*Ftmp18 + Ftmp134)*M[22];
double Ftmp136 = -8*Ftmp59;
double Ftmp137 = Ftmp126*Ftmp24;
double Ftmp138 = -Ftmp60;
double Ftmp139 = (Ftmp133*Ftmp21 + Ftmp134 + Ftmp138)*M[24];
double Ftmp140 = Ftmp18*Ftmp21;
double Ftmp141 = (Ftmp120*Ftmp140 + Ftmp132 + Ftmp138 + 1)*M[31];
double Ftmp142 = 15*Ftmp19;
double Ftmp143 = 6*y;
double Ftmp144 = Ftmp24*M[3];
double Ftmp145 = Ftmp30*y;
double Ftmp146 = pow(x, 3);
double Ftmp147 = Ftmp146*M[9];
double Ftmp148 = Ftmp103*M[13];
double Ftmp149 = Ftmp146*M[21];
double Ftmp150 = Ftmp38*Ftmp42;
double Ftmp151 = Ftmp46*M[10];
double Ftmp152 = Ftmp140*Ftmp30;
double Ftmp153 = Ftmp18*Ftmp42;
double Ftmp154 = Ftmp146*M[20];
double Ftmp155 = Ftmp140*Ftmp43;
double Ftmp156 = Ftmp2*Ftmp86;
double Ftmp157 = Ftmp55 - 1;
double Ftmp158 = Ftmp2*y;
double Ftmp159 = Ftmp76*M[27];
double Ftmp160 = 30*Ftmp157;
double Ftmp161 = Ftmp8*y;
double Ftmp162 = Ftmp10*z;
double Ftmp163 = 210*Ftmp157;
double Ftmp164 = Ftmp103*Ftmp18;
double Ftmp165 = Ftmp118*Ftmp29*Ftmp70;
double Ftmp166 = Ftmp37*Ftmp39;
double Ftmp167 = -8*Ftmp13;
double Ftmp168 = 1 - Ftmp55;
double Ftmp169 = 30*Ftmp161;
double Ftmp170 = Ftmp126*Ftmp19;
double Ftmp171 = 6*z;
double Ftmp172 = Ftmp30*z;
double Ftmp173 = Ftmp21*Ftmp42;
double Ftmp174 = Ftmp59 - 1;
double Ftmp175 = Ftmp174*z;
double Ftmp176 = Ftmp21*Ftmp9;
double Ftmp177 = Ftmp103*Ftmp21;
double Ftmp178 = Ftmp35*Ftmp39;
double Ftmp179 = Ftmp8*z;
double Ftmp180 = 1 - Ftmp59;
double Ftmp181 = 30*Ftmp179;
#pragma omp atomic
F[0] += Ftmp0*(Ftmp1*Ftmp17*Ftmp63*M[3] - Ftmp10*Ftmp71 - Ftmp10*Ftmp74 + Ftmp100*Ftmp101*M[21] + Ftmp100*Ftmp98*M[20] + Ftmp102*Ftmp104 + Ftmp104*Ftmp105 + Ftmp104*Ftmp106 + Ftmp107*Ftmp109 + Ftmp11*z + Ftmp111*Ftmp112 + Ftmp113*Ftmp71 + Ftmp113*Ftmp74 + Ftmp114*Ftmp115 + Ftmp114*Ftmp79 + Ftmp114*Ftmp80 + Ftmp118*Ftmp84*(Ftmp116*Ftmp117 - 10*Ftmp13 + 3) + Ftmp119*Ftmp121 + Ftmp119*Ftmp123 + Ftmp119*Ftmp125 + Ftmp119*Ftmp135 + Ftmp119*Ftmp139 + Ftmp119*Ftmp141 + Ftmp131*(Ftmp127 + Ftmp129 + Ftmp130)*M[22] + Ftmp131*(Ftmp130 + Ftmp136 + Ftmp137)*M[24] - Ftmp14*M[0] + Ftmp15*Ftmp16 - Ftmp16*Ftmp87 - Ftmp16*Ftmp88 - Ftmp16*Ftmp89 + Ftmp17*Ftmp20 + Ftmp17*Ftmp23 - Ftmp2*Ftmp3 + Ftmp2*Ftmp54 + Ftmp21*Ftmp49*Ftmp51 - Ftmp24*Ftmp83*M[9] + Ftmp25*Ftmp26 + Ftmp25*Ftmp3 - Ftmp25*Ftmp54 - Ftmp25*Ftmp96 - Ftmp25*Ftmp97 + Ftmp27*Ftmp49*M[25] - Ftmp28*Ftmp31 - Ftmp32*Ftmp34 - Ftmp32*Ftmp37*M[16] - Ftmp32*(-Ftmp126*Ftmp55 + Ftmp18 + Ftmp21)*M[31] - Ftmp36*M[17] - Ftmp38*Ftmp40*M[13] - Ftmp4*M[5] + Ftmp40*Ftmp77 + Ftmp41*Ftmp44 + Ftmp45*M[30] - Ftmp46*Ftmp47 - Ftmp48*M[14] + Ftmp49*Ftmp50 + Ftmp5*Ftmp64 + Ftmp5*Ftmp66 + Ftmp5*Ftmp68 + Ftmp52*M[26] + Ftmp58*M[12] - Ftmp6*M[1] + Ftmp62*M[14] + 210*Ftmp63*Ftmp95*M[23] - Ftmp7*x - Ftmp77*Ftmp9 - Ftmp78*Ftmp79 - Ftmp78*Ftmp80 - Ftmp82*M[28] - Ftmp85*y*M[10] - Ftmp85*z*M[11] - Ftmp90*Ftmp91 - Ftmp90*Ftmp92 - Ftmp90*Ftmp93 + M[0]);
#pragma omp atomic
F[1] += Ftmp0*(Ftmp1*Ftmp143*Ftmp157*M[6] + Ftmp10*Ftmp121 + Ftmp10*Ftmp123 + Ftmp10*Ftmp125 + Ftmp10*Ftmp135 + Ftmp10*Ftmp139 + Ftmp10*Ftmp141 + Ftmp102*Ftmp166 - Ftmp102*Ftmp78 + Ftmp104*Ftmp115 + Ftmp104*Ftmp79 + Ftmp104*Ftmp80 + Ftmp105*Ftmp166 - Ftmp105*Ftmp78 + Ftmp106*Ftmp166 + Ftmp107*Ftmp163*M[25] + Ftmp108*Ftmp161*(Ftmp117*Ftmp122 - 10*Ftmp55 + 3) + Ftmp111*Ftmp29*Ftmp35 - Ftmp119*Ftmp159 - Ftmp119*Ftmp71 - Ftmp119*Ftmp74 - Ftmp131*Ftmp157*y*M[12] + Ftmp142*Ftmp15 - Ftmp142*Ftmp87 - Ftmp142*Ftmp88 - Ftmp142*Ftmp89 + Ftmp142*x*M[4] + Ftmp143*Ftmp144 + Ftmp143*Ftmp23 - Ftmp145*Ftmp147 - Ftmp145*Ftmp34 - Ftmp145*(Ftmp12 - Ftmp126*Ftmp13 + Ftmp21)*M[24] - Ftmp148*Ftmp37 + Ftmp149*Ftmp150 - Ftmp151*Ftmp18 - Ftmp152*M[17] + Ftmp153*Ftmp154 + Ftmp153*Ftmp41 + Ftmp155*M[27] + Ftmp156*M[10] + Ftmp157*Ftmp37*Ftmp99*M[30] + Ftmp158*Ftmp64 + Ftmp158*Ftmp66 + Ftmp158*Ftmp68 + Ftmp159*Ftmp164 + Ftmp16*Ftmp26 - Ftmp16*Ftmp54 - Ftmp16*Ftmp96 - Ftmp16*Ftmp97 - Ftmp160*Ftmp161*z*M[16] - Ftmp160*Ftmp19*M[15] - Ftmp162*Ftmp91 - Ftmp162*Ftmp92 - Ftmp162*Ftmp93 + Ftmp163*Ftmp95*M[26] + Ftmp164*Ftmp71 + Ftmp164*Ftmp74 + Ftmp165*Ftmp98 + Ftmp169*(Ftmp129 + Ftmp167 + Ftmp168)*M[22] + Ftmp169*(Ftmp136 + Ftmp168 + Ftmp170)*M[31] + Ftmp2*Ftmp88 - Ftmp36*M[14] - Ftmp38*Ftmp46*M[11] - Ftmp4*M[7] + Ftmp44*Ftmp50 - Ftmp5*M[4] + Ftmp52*M[23] - Ftmp6*M[0] + Ftmp62*M[17] - Ftmp65*M[1] - Ftmp7*y - Ftmp82*M[32] + Ftmp90*M[13] + M[1]);
#pragma omp atomic
F[2] += Ftmp0*(Ftmp1*Ftmp171*Ftmp174*M[8] - Ftmp10*Ftmp102 - Ftmp10*Ftmp105 - Ftmp10*Ftmp106 + Ftmp10*Ftmp21*M[7] + Ftmp101*Ftmp165 + Ftmp102*Ftmp178 + Ftmp103*Ftmp77*z + Ftmp104*Ftmp71 + Ftmp104*Ftmp74 + Ftmp105*Ftmp178 + Ftmp106*Ftmp178 + Ftmp109*Ftmp29*Ftmp37 + Ftmp11*x + Ftmp110*Ftmp179*(Ftmp117*Ftmp124 - 10*Ftmp59 + 3) + 210*Ftmp112*Ftmp174*M[28] - Ftmp115*Ftmp119 + Ftmp115*Ftmp177 + Ftmp119*Ftmp21*M[5] - Ftmp119*Ftmp79 - Ftmp119*Ftmp80 + Ftmp121*Ftmp78 + Ftmp123*Ftmp78 + Ftmp125*Ftmp78 - Ftmp131*Ftmp175*M[14] + Ftmp135*Ftmp78 + Ftmp139*Ftmp78 + Ftmp141*Ftmp78 + Ftmp144*Ftmp171 - Ftmp147*Ftmp172 - Ftmp148*Ftmp35 + Ftmp149*Ftmp173 + Ftmp150*Ftmp154 - Ftmp151*Ftmp38 - Ftmp152*M[16] + Ftmp155*M[26] + Ftmp156*M[11] - Ftmp158*M[7] - Ftmp162*Ftmp87 - Ftmp162*Ftmp88 - Ftmp162*Ftmp89 - Ftmp169*Ftmp175*M[17] + Ftmp171*Ftmp20 - Ftmp172*(Ftmp12 - Ftmp128*Ftmp13 + Ftmp18)*M[22] + Ftmp173*Ftmp27*M[30] - 30*Ftmp174*Ftmp22*M[18] + Ftmp174*Ftmp35*Ftmp99*M[32] + 210*Ftmp175*Ftmp51*Ftmp94 - Ftmp176*Ftmp91 - Ftmp176*Ftmp92 - Ftmp176*Ftmp93 + Ftmp177*Ftmp79 + Ftmp177*Ftmp80 + Ftmp181*(Ftmp127 + Ftmp170 + Ftmp180)*M[31] + Ftmp181*(Ftmp137 + Ftmp167 + Ftmp180)*M[24] + Ftmp2*Ftmp93 - Ftmp27*Ftmp31*z + Ftmp3*Ftmp90 - Ftmp32*Ftmp47*z + Ftmp35*Ftmp49*M[23] + Ftmp4*Ftmp64 + Ftmp4*Ftmp66 + Ftmp4*Ftmp68 - Ftmp4*x*M[0] - Ftmp4*y*M[1] + Ftmp45*M[25] - Ftmp48*M[11] - Ftmp5*M[5] - Ftmp54*Ftmp90 + Ftmp58*M[16] - Ftmp67*M[2] - Ftmp90*Ftmp96 - Ftmp90*Ftmp97 + M[2]);

}

void S2Mc_4(double x, double y, double z, double * S, double * M) {
double Mtmp0 = x*S[0];
double Mtmp1 = z*S[2];
double Mtmp2 = -Mtmp1;
double Mtmp3 = x*S[1];
double Mtmp4 = y*S[0];
double Mtmp5 = x*S[2];
double Mtmp6 = z*S[0];
double Mtmp7 = y*S[1];
double Mtmp8 = y*S[2];
double Mtmp9 = z*S[1];
double Mtmp10 = Mtmp1*x;
double Mtmp11 = pow(x, 2);
double Mtmp12 = pow(z, 2);
double Mtmp13 = (1.0/2.0)*S[0];
double Mtmp14 = Mtmp12*Mtmp13;
double Mtmp15 = Mtmp0*y;
double Mtmp16 = Mtmp1*y;
double Mtmp17 = (1.0/2.0)*S[1];
double Mtmp18 = Mtmp12*Mtmp17;
double Mtmp19 = (1.0/2.0)*S[2];
double Mtmp20 = -Mtmp12*Mtmp19;
double Mtmp21 = Mtmp3*y;
double Mtmp22 = pow(y, 2);
double Mtmp23 = pow(x, 3);
double Mtmp24 = pow(z, 3);
double Mtmp25 = Mtmp24*S[2];
double Mtmp26 = 3*Mtmp12;
double Mtmp27 = Mtmp0*Mtmp26;
double Mtmp28 = 3*Mtmp11;
double Mtmp29 = Mtmp1*Mtmp28;
double Mtmp30 = (1.0/2.0)*Mtmp12;
double Mtmp31 = Mtmp10*y + (1.0/2.0)*Mtmp12*Mtmp4 + Mtmp3*Mtmp30;
double Mtmp32 = Mtmp24*S[0];
double Mtmp33 = Mtmp26*Mtmp7;
double Mtmp34 = 3*Mtmp22;
double Mtmp35 = Mtmp1*Mtmp34;
double Mtmp36 = Mtmp24*S[1];
double Mtmp37 = (1.0/2.0)*Mtmp11;
double Mtmp38 = pow(y, 3);
double Mtmp39 = (1.0/2.0)*Mtmp22;
#pragma omp atomic
M[0] += S[0];
#pragma omp atomic
M[1] += S[1];
#pragma omp atomic
M[2] += S[2];
#pragma omp atomic
M[3] += Mtmp0 + Mtmp2;
#pragma omp atomic
M[4] += Mtmp3 + Mtmp4;
#pragma omp atomic
M[5] += Mtmp5 + Mtmp6;
#pragma omp atomic
M[6] += Mtmp2 + Mtmp7;
#pragma omp atomic
M[7] += Mtmp8 + Mtmp9;
#pragma omp atomic
M[8] += -Mtmp10 + (1.0/2.0)*Mtmp11*S[0] - Mtmp14;
#pragma omp atomic
M[9] += Mtmp11*Mtmp17 + Mtmp15 - Mtmp16 - Mtmp18;
#pragma omp atomic
M[10] += Mtmp0*z + Mtmp11*Mtmp19 + Mtmp20;
#pragma omp atomic
M[11] += -Mtmp10 + Mtmp13*Mtmp22 - Mtmp14 + Mtmp21;
#pragma omp atomic
M[12] += Mtmp3*z + Mtmp4*z + Mtmp5*y;
#pragma omp atomic
M[13] += -Mtmp16 - Mtmp18 + (1.0/2.0)*Mtmp22*S[1];
#pragma omp atomic
M[14] += Mtmp19*Mtmp22 + Mtmp20 + Mtmp7*z;
#pragma omp atomic
M[15] += (1.0/6.0)*Mtmp23*S[0] + (1.0/6.0)*Mtmp25 - 1.0/6.0*Mtmp27 - 1.0/6.0*Mtmp29;
#pragma omp atomic
M[16] += (1.0/2.0)*Mtmp11*y*S[0] + (1.0/6.0)*Mtmp23*S[1] - Mtmp31;
#pragma omp atomic
M[17] += (1.0/6.0)*Mtmp23*S[2] - 1.0/6.0*Mtmp26*Mtmp5 + (1.0/6.0)*Mtmp28*Mtmp6 - 1.0/6.0*Mtmp32;
#pragma omp atomic
M[18] += (1.0/2.0)*Mtmp11*y*S[1] + (1.0/2.0)*Mtmp22*x*S[0] + (1.0/3.0)*Mtmp24*S[2] - 1.0/6.0*Mtmp27 - 1.0/6.0*Mtmp29 - 1.0/6.0*Mtmp33 - 1.0/6.0*Mtmp35;
#pragma omp atomic
M[19] += Mtmp15*z - Mtmp30*Mtmp8 - 1.0/6.0*Mtmp36 + Mtmp37*Mtmp8 + Mtmp37*Mtmp9;
#pragma omp atomic
M[20] += (1.0/2.0)*Mtmp22*x*S[1] - Mtmp31 + (1.0/6.0)*Mtmp38*S[0];
#pragma omp atomic
M[21] += Mtmp21*z - Mtmp30*Mtmp5 - 1.0/6.0*Mtmp32 + Mtmp39*Mtmp5 + Mtmp39*Mtmp6;
#pragma omp atomic
M[22] += (1.0/6.0)*Mtmp25 - 1.0/6.0*Mtmp33 - 1.0/6.0*Mtmp35 + (1.0/6.0)*Mtmp38*S[1];
#pragma omp atomic
M[23] += -1.0/6.0*Mtmp26*Mtmp8 + (1.0/6.0)*Mtmp34*Mtmp9 - 1.0/6.0*Mtmp36 + (1.0/6.0)*Mtmp38*S[2];

}

void M2Mc_4(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = z*M[2];
double Mstmp2 = -Mstmp1;
double Mstmp3 = x*M[1];
double Mstmp4 = y*M[0];
double Mstmp5 = x*M[2];
double Mstmp6 = z*M[0];
double Mstmp7 = y*M[1];
double Mstmp8 = y*M[2];
double Mstmp9 = z*M[1];
double Mstmp10 = x*M[3];
double Mstmp11 = pow(x, 2);
double Mstmp12 = (1.0/2.0)*M[0];
double Mstmp13 = z*M[5];
double Mstmp14 = pow(z, 2);
double Mstmp15 = Mstmp1*x;
double Mstmp16 = -Mstmp12*Mstmp14 - Mstmp13 - Mstmp15;
double Mstmp17 = x*M[4];
double Mstmp18 = y*M[3];
double Mstmp19 = (1.0/2.0)*M[1];
double Mstmp20 = Mstmp0*y;
double Mstmp21 = z*M[7];
double Mstmp22 = -Mstmp1*y - Mstmp14*Mstmp19 - Mstmp21;
double Mstmp23 = x*M[5];
double Mstmp24 = (1.0/2.0)*M[2];
double Mstmp25 = -Mstmp14*Mstmp24;
double Mstmp26 = x*M[6];
double Mstmp27 = y*M[4];
double Mstmp28 = pow(y, 2);
double Mstmp29 = Mstmp3*y;
double Mstmp30 = x*M[7];
double Mstmp31 = y*M[6];
double Mstmp32 = z*M[10];
double Mstmp33 = Mstmp13*x;
double Mstmp34 = (1.0/2.0)*M[3];
double Mstmp35 = (1.0/6.0)*pow(x, 3);
double Mstmp36 = Mstmp14*Mstmp34;
double Mstmp37 = pow(z, 3);
double Mstmp38 = (1.0/6.0)*Mstmp37;
double Mstmp39 = Mstmp38*M[2];
double Mstmp40 = (1.0/2.0)*Mstmp14;
double Mstmp41 = Mstmp0*Mstmp40;
double Mstmp42 = (1.0/2.0)*Mstmp11;
double Mstmp43 = Mstmp1*Mstmp42;
double Mstmp44 = (1.0/2.0)*M[4];
double Mstmp45 = -Mstmp13*y - Mstmp14*Mstmp44 - Mstmp15*y - Mstmp21*x - Mstmp3*Mstmp40 - Mstmp4*Mstmp40 - z*M[12];
double Mstmp46 = (1.0/2.0)*M[5];
double Mstmp47 = -Mstmp14*Mstmp46 - Mstmp38*M[0] - Mstmp40*Mstmp5;
double Mstmp48 = z*M[14];
double Mstmp49 = Mstmp21*y;
double Mstmp50 = (1.0/2.0)*M[6];
double Mstmp51 = Mstmp14*Mstmp50;
double Mstmp52 = Mstmp40*Mstmp7;
double Mstmp53 = (1.0/2.0)*Mstmp28;
double Mstmp54 = Mstmp1*Mstmp53;
double Mstmp55 = (1.0/2.0)*M[7];
double Mstmp56 = -Mstmp14*Mstmp55 - Mstmp38*M[1] - Mstmp40*Mstmp8;
double Mstmp57 = (1.0/6.0)*pow(y, 3);
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += Mstmp0 + Mstmp2 + M[3];
#pragma omp atomic
Ms[4] += Mstmp3 + Mstmp4 + M[4];
#pragma omp atomic
Ms[5] += Mstmp5 + Mstmp6 + M[5];
#pragma omp atomic
Ms[6] += Mstmp2 + Mstmp7 + M[6];
#pragma omp atomic
Ms[7] += Mstmp8 + Mstmp9 + M[7];
#pragma omp atomic
Ms[8] += Mstmp10 + Mstmp11*Mstmp12 + Mstmp16 + M[8];
#pragma omp atomic
Ms[9] += Mstmp11*Mstmp19 + Mstmp17 + Mstmp18 + Mstmp20 + Mstmp22 + M[9];
#pragma omp atomic
Ms[10] += Mstmp0*z + Mstmp11*Mstmp24 + Mstmp23 + Mstmp25 + z*M[3] + M[10];
#pragma omp atomic
Ms[11] += Mstmp12*Mstmp28 + Mstmp16 + Mstmp26 + Mstmp27 + Mstmp29 + M[11];
#pragma omp atomic
Ms[12] += Mstmp3*z + Mstmp30 + Mstmp4*z + Mstmp5*y + y*M[5] + z*M[4] + M[12];
#pragma omp atomic
Ms[13] += Mstmp19*Mstmp28 + Mstmp22 + Mstmp31 + M[13];
#pragma omp atomic
Ms[14] += Mstmp24*Mstmp28 + Mstmp25 + Mstmp7*z + y*M[7] + z*M[6] + M[14];
#pragma omp atomic
Ms[15] += Mstmp11*Mstmp34 - Mstmp32 - Mstmp33 + Mstmp35*M[0] - Mstmp36 + Mstmp39 - Mstmp41 - Mstmp43 + x*M[8] + M[15];
#pragma omp atomic
Ms[16] += Mstmp10*y + Mstmp11*Mstmp44 + Mstmp35*M[1] + Mstmp4*Mstmp42 + Mstmp45 + x*M[9] + y*M[8] + M[16];
#pragma omp atomic
Ms[17] += Mstmp10*z + Mstmp11*Mstmp46 + Mstmp35*M[2] + Mstmp42*Mstmp6 + Mstmp47 + x*M[10] + z*M[8] + M[17];
#pragma omp atomic
Ms[18] += (1.0/2.0)*Mstmp11*y*M[1] + (1.0/2.0)*Mstmp11*M[6] + (1.0/2.0)*Mstmp28*x*M[0] + (1.0/2.0)*Mstmp28*M[3] - Mstmp32 - Mstmp33 - Mstmp36 + (1.0/3.0)*Mstmp37*M[2] - Mstmp41 - Mstmp43 - Mstmp48 - Mstmp49 - Mstmp51 - Mstmp52 - Mstmp54 + x*y*M[4] + x*M[11] + y*M[9] + M[18];
#pragma omp atomic
Ms[19] += Mstmp11*Mstmp55 + Mstmp17*z + Mstmp18*z + Mstmp20*z + Mstmp23*y + Mstmp42*Mstmp8 + Mstmp42*Mstmp9 + Mstmp56 + x*M[12] + y*M[10] + z*M[9] + M[19];
#pragma omp atomic
Ms[20] += Mstmp26*y + Mstmp28*Mstmp44 + Mstmp3*Mstmp53 + Mstmp45 + Mstmp57*M[0] + x*M[13] + y*M[11] + M[20];
#pragma omp atomic
Ms[21] += Mstmp26*z + Mstmp27*z + Mstmp28*Mstmp46 + Mstmp29*z + Mstmp30*y + Mstmp47 + Mstmp5*Mstmp53 + Mstmp53*Mstmp6 + x*M[14] + y*M[12] + z*M[11] + M[21];
#pragma omp atomic
Ms[22] += Mstmp28*Mstmp50 + Mstmp39 - Mstmp48 - Mstmp49 - Mstmp51 - Mstmp52 - Mstmp54 + Mstmp57*M[1] + y*M[13] + M[22];
#pragma omp atomic
Ms[23] += Mstmp28*Mstmp55 + Mstmp31*z + Mstmp53*Mstmp9 + Mstmp56 + Mstmp57*M[2] + y*M[14] + z*M[13] + M[23];

}

void L2Lc_4(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[13];
double Lstmp4 = Lstmp3*y;
double Lstmp5 = (1.0/2.0)*pow(x, 2);
double Lstmp6 = (1.0/2.0)*pow(y, 2);
double Lstmp7 = x*L[12];
double Lstmp8 = y*L[10];
double Lstmp9 = z*L[11];
double Lstmp10 = z*L[15];
double Lstmp11 = L[4] + L[7];
double Lstmp12 = (1.0/2.0)*pow(z, 2);
double Lstmp13 = L[11] + L[15];
double Lstmp14 = L[9] + L[12];
double Lstmp15 = Lstmp12*Lstmp14;
double Lstmp16 = L[10] + L[14];
double Lstmp17 = Lstmp12*Lstmp16;
double Lstmp18 = y*L[12];
double Lstmp19 = y*L[13];
double Lstmp20 = Lstmp14*z;
double Lstmp21 = Lstmp16*z;
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp10*Lstmp6 - Lstmp11*Lstmp12 - 1.0/6.0*Lstmp13*pow(z, 3) - Lstmp15*x - Lstmp17*y + Lstmp2*y + Lstmp4*x + Lstmp5*Lstmp8 + Lstmp5*Lstmp9 + Lstmp5*L[4] + Lstmp6*Lstmp7 + Lstmp6*L[7] + (1.0/6.0)*pow(x, 3)*L[9] + x*L[1] + (1.0/6.0)*pow(y, 3)*L[14] + y*L[2] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 - Lstmp15 + Lstmp4 + Lstmp5*L[9] + Lstmp6*L[12] + Lstmp8*x + Lstmp9*x + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp10*y - Lstmp17 + Lstmp18*x + Lstmp2 + Lstmp3*x + Lstmp5*L[10] + Lstmp6*L[14] + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += -Lstmp11*z - Lstmp12*Lstmp13 + Lstmp19*x - Lstmp20*x - Lstmp21*y + Lstmp5*L[11] + Lstmp6*L[15] + x*L[6] + y*L[8] + L[3];
#pragma omp atomic
Ls[4] += Lstmp8 + Lstmp9 + x*L[9] + L[4];
#pragma omp atomic
Ls[5] += Lstmp18 + Lstmp3 + x*L[10] + L[5];
#pragma omp atomic
Ls[6] += Lstmp19 - Lstmp20 + x*L[11] + L[6];
#pragma omp atomic
Ls[7] += Lstmp10 + Lstmp7 + y*L[14] + L[7];
#pragma omp atomic
Ls[8] += -Lstmp21 + x*L[13] + y*L[15] + L[8];
#pragma omp atomic
Ls[9] += L[9];
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

}

void L2Pc_4(double x, double y, double z, double * L, double * F) {
double Ftmp0 = x*y;
double Ftmp1 = x*z;
double Ftmp2 = y*z;
double Ftmp3 = (1.0/2.0)*pow(x, 2);
double Ftmp4 = (1.0/2.0)*pow(y, 2);
double Ftmp5 = pow(z, 2);
double Ftmp6 = L[9] + L[12];
double Ftmp7 = L[10] + L[14];
#pragma omp atomic
F[0] += -Ftmp0*L[10] - Ftmp1*L[11] - Ftmp2*L[13] - Ftmp3*L[9] - Ftmp4*L[12] + (1.0/2.0)*Ftmp5*Ftmp6 - x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -Ftmp0*L[12] - Ftmp1*L[13] - Ftmp2*L[15] - Ftmp3*L[10] - Ftmp4*L[14] + (1.0/2.0)*Ftmp5*Ftmp7 - x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -Ftmp0*L[13] - Ftmp3*L[11] - Ftmp4*L[15] + (1.0/2.0)*Ftmp5*(L[11] + L[15]) + Ftmp6*x*z + Ftmp7*y*z - x*L[6] - y*L[8] + z*(L[4] + L[7]) - L[3];

}

void M2Pc_4(double x, double y, double z, double * M, double * F) {
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(Rinv, 3);
double Ftmp1 = pow(Rinv, 2);
double Ftmp2 = 3*Ftmp1;
double Ftmp3 = y*M[4];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = Ftmp2*x;
double Ftmp6 = Ftmp5*y;
double Ftmp7 = Ftmp4*M[2];
double Ftmp8 = pow(Rinv, 4);
double Ftmp9 = 15*Ftmp8;
double Ftmp10 = Ftmp9*y;
double Ftmp11 = Ftmp10*M[12];
double Ftmp12 = pow(x, 2);
double Ftmp13 = Ftmp1*Ftmp12;
double Ftmp14 = 3*Ftmp13;
double Ftmp15 = Ftmp10*M[7];
double Ftmp16 = x*z;
double Ftmp17 = 6*x;
double Ftmp18 = pow(y, 2);
double Ftmp19 = Ftmp18*Ftmp8;
double Ftmp20 = Ftmp19*M[6];
double Ftmp21 = Ftmp12*Ftmp8;
double Ftmp22 = 15*Ftmp21;
double Ftmp23 = z*M[5];
double Ftmp24 = pow(y, 3);
double Ftmp25 = Ftmp24*x;
double Ftmp26 = pow(Rinv, 6);
double Ftmp27 = 30*Ftmp26;
double Ftmp28 = Ftmp27*M[13];
double Ftmp29 = Ftmp18*M[14];
double Ftmp30 = Ftmp27*z;
double Ftmp31 = Ftmp30*x;
double Ftmp32 = y*z;
double Ftmp33 = Ftmp12*Ftmp26;
double Ftmp34 = 105*Ftmp33;
double Ftmp35 = 210*pow(Rinv, 8);
double Ftmp36 = Ftmp35*z;
double Ftmp37 = Ftmp25*Ftmp36;
double Ftmp38 = Ftmp12*Ftmp18;
double Ftmp39 = Ftmp27*Ftmp38;
double Ftmp40 = Ftmp24*Ftmp35;
double Ftmp41 = Ftmp36*Ftmp38;
double Ftmp42 = 5*Ftmp13;
double Ftmp43 = (Ftmp42 - 3)*M[8];
double Ftmp44 = Ftmp1*Ftmp18;
double Ftmp45 = 5*Ftmp44;
double Ftmp46 = Ftmp45 - 1;
double Ftmp47 = Ftmp2*Ftmp46;
double Ftmp48 = Ftmp13 - 1;
double Ftmp49 = (Ftmp14 - 1)*M[3];
double Ftmp50 = 3*Ftmp44;
double Ftmp51 = (Ftmp50 - 1)*M[6];
double Ftmp52 = y*M[16];
double Ftmp53 = 7*Ftmp13;
double Ftmp54 = Ftmp53 - 3;
double Ftmp55 = Ftmp54*Ftmp9;
double Ftmp56 = 7*Ftmp44;
double Ftmp57 = Ftmp56 - 3;
double Ftmp58 = Ftmp57*M[20];
double Ftmp59 = Ftmp55*M[17];
double Ftmp60 = Ftmp9*z;
double Ftmp61 = (Ftmp56 - 1)*M[21];
double Ftmp62 = 30*Ftmp48;
double Ftmp63 = Ftmp8*x;
double Ftmp64 = Ftmp62*Ftmp63;
double Ftmp65 = Ftmp42 - 1;
double Ftmp66 = Ftmp65*M[9];
double Ftmp67 = Ftmp10*x;
double Ftmp68 = (Ftmp45 - 3)*M[13];
double Ftmp69 = Ftmp60*x;
double Ftmp70 = Ftmp65*M[10];
double Ftmp71 = Ftmp46*M[14];
double Ftmp72 = 210*Ftmp48;
double Ftmp73 = Ftmp32*x;
double Ftmp74 = Ftmp26*Ftmp73;
double Ftmp75 = Ftmp46*M[11];
double Ftmp76 = Ftmp33*Ftmp72;
double Ftmp77 = z*M[17];
double Ftmp78 = (Ftmp53 - 1)*M[19];
double Ftmp79 = 105*Ftmp26;
double Ftmp80 = Ftmp78*Ftmp79;
double Ftmp81 = Ftmp57*M[23];
double Ftmp82 = Ftmp79*Ftmp81;
double Ftmp83 = Ftmp18*Ftmp26;
double Ftmp84 = Ftmp83*x;
double Ftmp85 = 60*M[22];
double Ftmp86 = Ftmp57*Ftmp85;
double Ftmp87 = Ftmp34*Ftmp54;
double Ftmp88 = pow(x, 4);
double Ftmp89 = 7*Ftmp8;
double Ftmp90 = 60*M[15];
double Ftmp91 = Ftmp9*x;
double Ftmp92 = 35*Ftmp8;
double Ftmp93 = (-30*Ftmp13 + Ftmp88*Ftmp92 + 3)*M[15];
double Ftmp94 = pow(y, 4);
double Ftmp95 = (-30*Ftmp44 + Ftmp92*Ftmp94 + 3)*M[22];
double Ftmp96 = 14*Ftmp38*Ftmp8 + 1;
double Ftmp97 = 30*M[18];
double Ftmp98 = (Ftmp38*Ftmp92 - Ftmp42 - Ftmp45 + 1)*M[18];
double Ftmp99 = 6*y;
double Ftmp100 = Ftmp21*M[3];
double Ftmp101 = Ftmp27*y;
double Ftmp102 = pow(x, 3);
double Ftmp103 = Ftmp102*M[8];
double Ftmp104 = 105*Ftmp83;
double Ftmp105 = Ftmp12*M[10];
double Ftmp106 = Ftmp101*z;
double Ftmp107 = Ftmp102*Ftmp35;
double Ftmp108 = Ftmp107*M[17];
double Ftmp109 = Ftmp2*Ftmp65;
double Ftmp110 = Ftmp44 - 1;
double Ftmp111 = Ftmp2*y;
double Ftmp112 = x*M[16];
double Ftmp113 = 30*Ftmp110;
double Ftmp114 = Ftmp8*y;
double Ftmp115 = Ftmp10*z;
double Ftmp116 = 210*Ftmp110;
double Ftmp117 = Ftmp18*Ftmp9;
double Ftmp118 = Ftmp83*z;
double Ftmp119 = Ftmp54*M[17];
double Ftmp120 = Ftmp73*Ftmp79;
double Ftmp121 = 105*Ftmp58;
double Ftmp122 = Ftmp33*Ftmp54*Ftmp90;
double Ftmp123 = Ftmp104*z;
double Ftmp124 = pow(z, 2);
double Ftmp125 = 6*z;
double Ftmp126 = Ftmp124*x;
double Ftmp127 = Ftmp126*Ftmp79;
double Ftmp128 = Ftmp124*Ftmp27;
double Ftmp129 = Ftmp124*y;
double Ftmp130 = Ftmp124*Ftmp9;
#pragma omp atomic
F[0] += Ftmp0*(Ftmp1*Ftmp17*Ftmp48*M[3] - Ftmp10*Ftmp58 + Ftmp11*z + Ftmp12*Ftmp40*M[20] - Ftmp14*M[0] + Ftmp15*Ftmp16 + Ftmp17*Ftmp20 - Ftmp2*Ftmp3 + Ftmp2*Ftmp43 - Ftmp21*Ftmp62*M[8] + Ftmp22*Ftmp23 + Ftmp22*Ftmp3 - Ftmp22*Ftmp43 - Ftmp22*Ftmp75 - Ftmp25*Ftmp28 - Ftmp29*Ftmp31 - Ftmp32*Ftmp34*M[12] + Ftmp34*Ftmp58*y + Ftmp34*Ftmp61*z + Ftmp37*M[23] - Ftmp39*M[11] - Ftmp4*M[5] + Ftmp41*M[21] + Ftmp47*M[11] + Ftmp49*Ftmp5 + Ftmp5*Ftmp51 - Ftmp52*Ftmp55 + Ftmp52*Ftmp76 + Ftmp52*Ftmp87 - Ftmp59*z - Ftmp6*M[1] - Ftmp60*Ftmp61 + Ftmp63*Ftmp90*(-10*Ftmp13 + Ftmp88*Ftmp89 + 3) + Ftmp63*Ftmp97*(-Ftmp13 - 8*Ftmp44 + Ftmp96) - Ftmp64*y*M[9] - Ftmp64*z*M[10] - Ftmp66*Ftmp67 - Ftmp67*Ftmp68 - Ftmp69*Ftmp70 - Ftmp69*Ftmp71 - Ftmp7*x + Ftmp72*Ftmp74*M[19] + Ftmp73*Ftmp80 + Ftmp73*Ftmp82 + Ftmp76*Ftmp77 + Ftmp77*Ftmp87 + Ftmp84*Ftmp86 + Ftmp91*Ftmp93 + Ftmp91*Ftmp95 + Ftmp91*Ftmp98 + M[0]);
#pragma omp atomic
F[1] += Ftmp0*(Ftmp1*Ftmp110*Ftmp99*M[6] + Ftmp10*Ftmp93 + Ftmp10*Ftmp95 + Ftmp10*Ftmp98 + Ftmp100*Ftmp99 - Ftmp101*Ftmp103 + Ftmp104*Ftmp112*Ftmp54 - Ftmp104*Ftmp16*M[12] - Ftmp105*Ftmp106 + Ftmp107*Ftmp18*M[16] + Ftmp108*Ftmp32 + Ftmp109*M[9] + Ftmp111*Ftmp49 + Ftmp111*Ftmp51 - Ftmp112*Ftmp55 - Ftmp113*Ftmp114*z*M[14] - Ftmp113*Ftmp19*M[13] - Ftmp113*Ftmp63*y*M[11] + Ftmp114*Ftmp85*(-10*Ftmp44 + Ftmp89*Ftmp94 + 3) + Ftmp114*Ftmp97*(-8*Ftmp13 - Ftmp44 + Ftmp96) - Ftmp115*Ftmp70 - Ftmp115*Ftmp71 + Ftmp116*Ftmp118*M[23] + Ftmp116*Ftmp74*M[21] + Ftmp116*Ftmp84*M[20] - Ftmp117*Ftmp66 - Ftmp117*Ftmp68 + Ftmp119*Ftmp120 + Ftmp120*Ftmp61 + Ftmp121*Ftmp84 + Ftmp122*y + Ftmp123*Ftmp78 + Ftmp123*Ftmp81 + Ftmp18*Ftmp60*M[7] + Ftmp18*Ftmp91*M[4] + Ftmp2*Ftmp68 + Ftmp23*Ftmp67 - Ftmp39*M[9] - Ftmp4*M[7] + Ftmp41*M[19] - Ftmp43*Ftmp67 - Ftmp5*M[4] - Ftmp50*M[1] - Ftmp58*Ftmp91 - Ftmp6*M[0] - Ftmp60*Ftmp78 - Ftmp60*Ftmp81 - Ftmp67*Ftmp75 + Ftmp69*M[12] - Ftmp7*y + M[1]);
#pragma omp atomic
F[2] += Ftmp0*(-Ftmp10*Ftmp78 - Ftmp10*Ftmp81 + Ftmp100*Ftmp125 - Ftmp103*Ftmp30 - Ftmp105*Ftmp128 - Ftmp106*Ftmp12*M[9] + Ftmp107*Ftmp52*z + Ftmp108*Ftmp124 + Ftmp109*M[10] + Ftmp11*x - Ftmp111*M[7] - Ftmp115*Ftmp66 - Ftmp115*Ftmp68 + Ftmp118*Ftmp86 + Ftmp119*Ftmp127 + Ftmp12*Ftmp129*Ftmp35*M[19] + Ftmp121*Ftmp74 + Ftmp122*z + Ftmp124*Ftmp15 - Ftmp124*Ftmp2*M[2] + Ftmp124*Ftmp40*M[23] + Ftmp124*Ftmp91*M[5] + Ftmp125*Ftmp20 + Ftmp126*Ftmp18*Ftmp35*M[21] + Ftmp127*Ftmp61 - Ftmp127*y*M[12] - Ftmp128*Ftmp29 + Ftmp129*Ftmp80 + Ftmp129*Ftmp82 - Ftmp130*Ftmp70 - Ftmp130*Ftmp71 + Ftmp16*Ftmp52*Ftmp54*Ftmp79 - Ftmp18*Ftmp31*M[11] - Ftmp24*Ftmp28*z + Ftmp3*Ftmp69 - Ftmp30*(Ftmp12 - 14*Ftmp13*Ftmp18 + Ftmp18)*M[18] + Ftmp37*M[20] + Ftmp4*Ftmp49 + Ftmp4*Ftmp51 - Ftmp4*x*M[0] - Ftmp4*y*M[1] - Ftmp43*Ftmp69 + Ftmp47*M[14] - Ftmp5*M[5] - Ftmp59*x + Ftmp60*Ftmp93 + Ftmp60*Ftmp95 + Ftmp60*Ftmp98 - Ftmp61*Ftmp91 - Ftmp69*Ftmp75 + M[2]);

}

void M2Lc_4(double x, double y, double z, double * M, double * L) {
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double D[30];
double Dtmp0 = pow(Rinv, 3);
double Dtmp1 = pow(x, 2);
double Dtmp2 = pow(Rinv, 2);
double Dtmp3 = 3*Dtmp2;
double Dtmp4 = 3*pow(Rinv, 5);
double Dtmp5 = Dtmp4*x;
double Dtmp6 = pow(y, 2);
double Dtmp7 = Dtmp4*y;
double Dtmp8 = 5*Dtmp2;
double Dtmp9 = Dtmp1*Dtmp8;
double Dtmp10 = Dtmp9 - 1;
double Dtmp11 = Dtmp4*z;
double Dtmp12 = Dtmp6*Dtmp8;
double Dtmp13 = Dtmp12 - 1;
double Dtmp14 = pow(Rinv, 7);
double Dtmp15 = 15*Dtmp14*x*y;
double Dtmp16 = Dtmp1*Dtmp2;
double Dtmp17 = 35*pow(Rinv, 4);
double Dtmp18 = 7*Dtmp16;
double Dtmp19 = Dtmp18 - 3;
double Dtmp20 = 15*Dtmp14*z;
double Dtmp21 = Dtmp20*x;
double Dtmp22 = Dtmp20*y;
double Dtmp23 = Dtmp2*Dtmp6;
double Dtmp24 = 7*Dtmp23;
double Dtmp25 = Dtmp24 - 3;
D[0] = -Dtmp0*x;
D[1] = -Dtmp0*y;
D[2] = -Dtmp0*z;
D[3] = Dtmp0*(Dtmp1*Dtmp3 - 1);
D[4] = Dtmp5*y;
D[5] = Dtmp5*z;
D[6] = Dtmp0*(Dtmp3*Dtmp6 - 1);
D[7] = Dtmp7*z;
D[8] = -D[3] - D[6];
D[9] = -Dtmp5*(Dtmp9 - 3);
D[10] = -Dtmp10*Dtmp7;
D[11] = -Dtmp10*Dtmp11;
D[12] = -Dtmp13*Dtmp5;
D[13] = -Dtmp15*z;
D[14] = -D[9] - D[12];
D[15] = -Dtmp7*(Dtmp12 - 3);
D[16] = -Dtmp11*Dtmp13;
D[17] = -D[10] - D[15];
D[18] = Dtmp4*(-30*Dtmp16 + Dtmp17*pow(x, 4) + 3);
D[19] = Dtmp15*Dtmp19;
D[20] = Dtmp19*Dtmp21;
D[21] = Dtmp4*(Dtmp1*Dtmp17*Dtmp6 - Dtmp12 - Dtmp9 + 1);
D[22] = Dtmp22*(Dtmp18 - 1);
D[23] = -D[18] - D[21];
D[24] = Dtmp15*Dtmp25;
D[25] = Dtmp21*(Dtmp24 - 1);
D[26] = -D[19] - D[24];
D[27] = Dtmp4*(Dtmp17*pow(y, 4) - 30*Dtmp23 + 3);
D[28] = Dtmp22*Dtmp25;
D[29] = -D[21] - D[27];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[9]*M[8] + D[10]*M[9] + D[11]*M[10] + D[12]*M[11] + D[13]*M[12] + D[15]*M[13] + D[16]*M[14] + D[18]*M[15] + D[19]*M[16] + D[20]*M[17] + D[21]*M[18] + D[22]*M[19] + D[24]*M[20] + D[25]*M[21] + D[27]*M[22] + D[28]*M[23];
#pragma omp atomic
L[1] += D[3]*M[0] + D[4]*M[1] + D[5]*M[2] + D[9]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[18]*M[8] + D[19]*M[9] + D[20]*M[10] + D[21]*M[11] + D[22]*M[12] + D[24]*M[13] + D[25]*M[14];
#pragma omp atomic
L[2] += D[4]*M[0] + D[6]*M[1] + D[7]*M[2] + D[10]*M[3] + D[12]*M[4] + D[13]*M[5] + D[15]*M[6] + D[16]*M[7] + D[19]*M[8] + D[21]*M[9] + D[22]*M[10] + D[24]*M[11] + D[25]*M[12] + D[27]*M[13] + D[28]*M[14];
#pragma omp atomic
L[3] += D[5]*M[0] + D[7]*M[1] + D[8]*M[2] + D[11]*M[3] + D[13]*M[4] + D[14]*M[5] + D[16]*M[6] + D[17]*M[7] + D[20]*M[8] + D[22]*M[9] + D[23]*M[10] + D[25]*M[11] + D[26]*M[12] + D[28]*M[13] + D[29]*M[14];
#pragma omp atomic
L[4] += D[9]*M[0] + D[10]*M[1] + D[11]*M[2] + D[18]*M[3] + D[19]*M[4] + D[20]*M[5] + D[21]*M[6] + D[22]*M[7];
#pragma omp atomic
L[5] += D[10]*M[0] + D[12]*M[1] + D[13]*M[2] + D[19]*M[3] + D[21]*M[4] + D[22]*M[5] + D[24]*M[6] + D[25]*M[7];
#pragma omp atomic
L[6] += D[11]*M[0] + D[13]*M[1] + D[14]*M[2] + D[20]*M[3] + D[22]*M[4] + D[23]*M[5] + D[25]*M[6] + D[26]*M[7];
#pragma omp atomic
L[7] += D[12]*M[0] + D[15]*M[1] + D[16]*M[2] + D[21]*M[3] + D[24]*M[4] + D[25]*M[5] + D[27]*M[6] + D[28]*M[7];
#pragma omp atomic
L[8] += D[13]*M[0] + D[16]*M[1] + D[17]*M[2] + D[22]*M[3] + D[25]*M[4] + D[26]*M[5] + D[28]*M[6] + D[29]*M[7];
#pragma omp atomic
L[9] += D[18]*M[0] + D[19]*M[1] + D[20]*M[2];
#pragma omp atomic
L[10] += D[19]*M[0] + D[21]*M[1] + D[22]*M[2];
#pragma omp atomic
L[11] += D[20]*M[0] + D[22]*M[1] + D[23]*M[2];
#pragma omp atomic
L[12] += D[21]*M[0] + D[24]*M[1] + D[25]*M[2];
#pragma omp atomic
L[13] += D[22]*M[0] + D[25]*M[1] + D[26]*M[2];
#pragma omp atomic
L[14] += D[24]*M[0] + D[27]*M[1] + D[28]*M[2];
#pragma omp atomic
L[15] += D[25]*M[0] + D[28]*M[1] + D[29]*M[2];

}

void S2M_5(double x, double y, double z, double * S, double * M) {
double Mtmp0 = x*S[1];
double Mtmp1 = y*S[0];
double Mtmp2 = Mtmp0 + Mtmp1;
double Mtmp3 = x*S[2];
double Mtmp4 = z*S[0];
double Mtmp5 = Mtmp3 + Mtmp4;
double Mtmp6 = y*S[2];
double Mtmp7 = z*S[1];
double Mtmp8 = Mtmp6 + Mtmp7;
double Mtmp9 = pow(x, 2);
double Mtmp10 = (1.0/2.0)*Mtmp0;
double Mtmp11 = (1.0/2.0)*Mtmp3;
double Mtmp12 = (1.0/2.0)*Mtmp1;
double Mtmp13 = Mtmp1*z;
double Mtmp14 = Mtmp3*y;
double Mtmp15 = Mtmp0*z;
double Mtmp16 = Mtmp14 + Mtmp15;
double Mtmp17 = pow(y, 2);
double Mtmp18 = pow(z, 2);
double Mtmp19 = pow(x, 3);
double Mtmp20 = 3*Mtmp1;
double Mtmp21 = (1.0/6.0)*Mtmp9;
double Mtmp22 = 3*Mtmp4;
double Mtmp23 = (1.0/2.0)*x;
double Mtmp24 = Mtmp11*y;
double Mtmp25 = Mtmp10*z;
double Mtmp26 = 3*Mtmp0;
double Mtmp27 = (1.0/6.0)*Mtmp17;
double Mtmp28 = Mtmp12*z;
double Mtmp29 = 3*Mtmp3;
double Mtmp30 = (1.0/6.0)*Mtmp18;
double Mtmp31 = pow(y, 3);
double Mtmp32 = 3*Mtmp7;
double Mtmp33 = y*z;
double Mtmp34 = 3*Mtmp6;
double Mtmp35 = pow(z, 3);
double Mtmp36 = (1.0/24.0)*Mtmp19;
double Mtmp37 = (1.0/12.0)*Mtmp9;
double Mtmp38 = (1.0/12.0)*x;
double Mtmp39 = 2*Mtmp15;
double Mtmp40 = 2*Mtmp13;
double Mtmp41 = (1.0/4.0)*x;
double Mtmp42 = 2*Mtmp14;
double Mtmp43 = (1.0/24.0)*Mtmp31;
double Mtmp44 = (1.0/24.0)*Mtmp35;
#pragma omp atomic
M[0] += S[0];
#pragma omp atomic
M[1] += S[1];
#pragma omp atomic
M[2] += S[2];
#pragma omp atomic
M[3] += x*S[0];
#pragma omp atomic
M[4] += Mtmp2;
#pragma omp atomic
M[5] += Mtmp5;
#pragma omp atomic
M[6] += y*S[1];
#pragma omp atomic
M[7] += Mtmp8;
#pragma omp atomic
M[8] += z*S[2];
#pragma omp atomic
M[9] += (1.0/2.0)*Mtmp9*S[0];
#pragma omp atomic
M[10] += x*(Mtmp1 + Mtmp10);
#pragma omp atomic
M[11] += x*(Mtmp11 + Mtmp4);
#pragma omp atomic
M[12] += y*(Mtmp0 + Mtmp12);
#pragma omp atomic
M[13] += Mtmp13 + Mtmp16;
#pragma omp atomic
M[14] += z*(Mtmp3 + (1.0/2.0)*Mtmp4);
#pragma omp atomic
M[15] += (1.0/2.0)*Mtmp17*S[1];
#pragma omp atomic
M[16] += y*((1.0/2.0)*Mtmp6 + Mtmp7);
#pragma omp atomic
M[17] += z*(Mtmp6 + (1.0/2.0)*Mtmp7);
#pragma omp atomic
M[18] += (1.0/2.0)*Mtmp18*S[2];
#pragma omp atomic
M[19] += (1.0/6.0)*Mtmp19*S[0];
#pragma omp atomic
M[20] += Mtmp21*(Mtmp0 + Mtmp20);
#pragma omp atomic
M[21] += Mtmp21*(Mtmp22 + Mtmp3);
#pragma omp atomic
M[22] += Mtmp2*Mtmp23*y;
#pragma omp atomic
M[23] += x*(Mtmp13 + Mtmp24 + Mtmp25);
#pragma omp atomic
M[24] += Mtmp23*Mtmp5*z;
#pragma omp atomic
M[25] += Mtmp27*(Mtmp1 + Mtmp26);
#pragma omp atomic
M[26] += y*(Mtmp15 + Mtmp24 + Mtmp28);
#pragma omp atomic
M[27] += z*(Mtmp14 + Mtmp25 + Mtmp28);
#pragma omp atomic
M[28] += Mtmp30*(Mtmp29 + Mtmp4);
#pragma omp atomic
M[29] += (1.0/6.0)*Mtmp31*S[1];
#pragma omp atomic
M[30] += Mtmp27*(Mtmp32 + Mtmp6);
#pragma omp atomic
M[31] += (1.0/2.0)*Mtmp33*Mtmp8;
#pragma omp atomic
M[32] += Mtmp30*(Mtmp34 + Mtmp7);
#pragma omp atomic
M[33] += (1.0/6.0)*Mtmp35*S[2];
#pragma omp atomic
M[34] += (1.0/24.0)*pow(x, 4)*S[0];
#pragma omp atomic
M[35] += Mtmp36*(Mtmp0 + 4*Mtmp1);
#pragma omp atomic
M[36] += Mtmp36*(Mtmp3 + 4*Mtmp4);
#pragma omp atomic
M[37] += Mtmp37*y*(2*Mtmp0 + Mtmp20);
#pragma omp atomic
M[38] += Mtmp21*(3*Mtmp13 + Mtmp16);
#pragma omp atomic
M[39] += Mtmp37*z*(Mtmp22 + 2*Mtmp3);
#pragma omp atomic
M[40] += Mtmp17*Mtmp38*(2*Mtmp1 + Mtmp26);
#pragma omp atomic
M[41] += Mtmp41*y*(Mtmp14 + Mtmp39 + Mtmp40);
#pragma omp atomic
M[42] += Mtmp41*z*(Mtmp15 + Mtmp40 + Mtmp42);
#pragma omp atomic
M[43] += Mtmp18*Mtmp38*(Mtmp29 + 2*Mtmp4);
#pragma omp atomic
M[44] += Mtmp43*(4*Mtmp0 + Mtmp1);
#pragma omp atomic
M[45] += Mtmp27*(Mtmp13 + Mtmp14 + 3*Mtmp15);
#pragma omp atomic
M[46] += (1.0/4.0)*Mtmp33*(Mtmp13 + Mtmp39 + Mtmp42);
#pragma omp atomic
M[47] += Mtmp30*(Mtmp13 + 3*Mtmp14 + Mtmp15);
#pragma omp atomic
M[48] += Mtmp44*(4*Mtmp3 + Mtmp4);
#pragma omp atomic
M[49] += (1.0/24.0)*pow(y, 4)*S[1];
#pragma omp atomic
M[50] += Mtmp43*(Mtmp6 + 4*Mtmp7);
#pragma omp atomic
M[51] += (1.0/12.0)*Mtmp17*z*(Mtmp32 + 2*Mtmp6);
#pragma omp atomic
M[52] += (1.0/12.0)*Mtmp18*y*(Mtmp34 + 2*Mtmp7);
#pragma omp atomic
M[53] += Mtmp44*(4*Mtmp6 + Mtmp7);
#pragma omp atomic
M[54] += (1.0/24.0)*pow(z, 4)*S[2];

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
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double D[55];
double Dtmp0 = pow(Rinv, 3);
double Dtmp1 = pow(x, 2);
double Dtmp2 = pow(Rinv, 2);
double Dtmp3 = 3*Dtmp2;
double Dtmp4 = Dtmp1*Dtmp3 - 1;
double Dtmp5 = 3*pow(Rinv, 5);
double Dtmp6 = Dtmp5*x;
double Dtmp7 = pow(y, 2);
double Dtmp8 = Dtmp3*Dtmp7 - 1;
double Dtmp9 = Dtmp5*y;
double Dtmp10 = 5*Dtmp2;
double Dtmp11 = Dtmp1*Dtmp10;
double Dtmp12 = Dtmp11 - 1;
double Dtmp13 = Dtmp5*z;
double Dtmp14 = Dtmp10*Dtmp7;
double Dtmp15 = Dtmp14 - 1;
double Dtmp16 = pow(Rinv, 7);
double Dtmp17 = 15*Dtmp16;
double Dtmp18 = Dtmp17*x;
double Dtmp19 = Dtmp18*y;
double Dtmp20 = Dtmp1*Dtmp2;
double Dtmp21 = pow(x, 4);
double Dtmp22 = pow(Rinv, 4);
double Dtmp23 = 35*Dtmp22;
double Dtmp24 = 7*Dtmp20;
double Dtmp25 = Dtmp24 - 3;
double Dtmp26 = Dtmp18*z;
double Dtmp27 = Dtmp1*Dtmp7;
double Dtmp28 = Dtmp17*y;
double Dtmp29 = Dtmp28*z;
double Dtmp30 = Dtmp2*Dtmp7;
double Dtmp31 = 7*Dtmp30;
double Dtmp32 = Dtmp31 - 3;
double Dtmp33 = pow(y, 4);
double Dtmp34 = Dtmp21*Dtmp22;
double Dtmp35 = 45*Dtmp16;
double Dtmp36 = Dtmp35*(-14*Dtmp20 + 21*Dtmp34 + 1);
double Dtmp37 = -Dtmp24;
double Dtmp38 = 63*Dtmp22*Dtmp27;
double Dtmp39 = Dtmp38 + 3;
double Dtmp40 = 315*pow(Rinv, 9)*x*y*z;
double Dtmp41 = -Dtmp31;
double Dtmp42 = Dtmp22*Dtmp33;
double Dtmp43 = Dtmp35*(-14*Dtmp30 + 21*Dtmp42 + 1);
D[0] = -Dtmp0*x;
D[1] = -Dtmp0*y;
D[2] = -Dtmp0*z;
D[3] = Dtmp0*Dtmp4;
D[4] = Dtmp6*y;
D[5] = Dtmp6*z;
D[6] = Dtmp0*Dtmp8;
D[7] = Dtmp9*z;
D[8] = -D[3] - D[6];
D[9] = -Dtmp6*(Dtmp11 - 3);
D[10] = -Dtmp12*Dtmp9;
D[11] = -Dtmp12*Dtmp13;
D[12] = -Dtmp15*Dtmp6;
D[13] = -Dtmp19*z;
D[14] = -D[9] - D[12];
D[15] = -Dtmp9*(Dtmp14 - 3);
D[16] = -Dtmp13*Dtmp15;
D[17] = -D[10] - D[15];
D[18] = -D[11] - D[16];
D[19] = Dtmp5*(-30*Dtmp20 + Dtmp21*Dtmp23 + 3);
D[20] = Dtmp19*Dtmp25;
D[21] = Dtmp25*Dtmp26;
D[22] = Dtmp5*(-Dtmp11 - Dtmp14 + Dtmp23*Dtmp27 + 1);
D[23] = Dtmp29*(Dtmp24 - 1);
D[24] = -D[19] - D[22];
D[25] = Dtmp19*Dtmp32;
D[26] = Dtmp26*(Dtmp31 - 1);
D[27] = -D[20] - D[25];
D[28] = -D[21] - D[26];
D[29] = Dtmp5*(Dtmp23*Dtmp33 - 30*Dtmp30 + 3);
D[30] = Dtmp29*Dtmp32;
D[31] = -D[22] - D[29];
D[32] = -D[23] - D[30];
D[33] = -D[24] - D[31];
D[34] = -Dtmp18*(-70*Dtmp20 + 63*Dtmp34 + 15);
D[35] = -Dtmp36*y;
D[36] = -Dtmp36*z;
D[37] = -Dtmp18*(-21*Dtmp30 + Dtmp37 + Dtmp39);
D[38] = -Dtmp4*Dtmp40;
D[39] = -D[34] - D[37];
D[40] = -Dtmp28*(-21*Dtmp20 + Dtmp39 + Dtmp41);
D[41] = -Dtmp17*z*(Dtmp37 + Dtmp38 + Dtmp41 + 1);
D[42] = -D[35] - D[40];
D[43] = -D[36] - D[41];
D[44] = -Dtmp43*x;
D[45] = -Dtmp40*Dtmp8;
D[46] = -D[37] - D[44];
D[47] = -D[38] - D[45];
D[48] = -D[39] - D[46];
D[49] = -Dtmp28*(-70*Dtmp30 + 63*Dtmp42 + 15);
D[50] = -Dtmp43*z;
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
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(Rinv, 3);
double Ftmp1 = pow(Rinv, 2);
double Ftmp2 = 3*Ftmp1;
double Ftmp3 = Ftmp2*z;
double Ftmp4 = Ftmp2*x;
double Ftmp5 = Ftmp4*y;
double Ftmp6 = Ftmp3*M[2];
double Ftmp7 = pow(Rinv, 4);
double Ftmp8 = pow(x, 2);
double Ftmp9 = Ftmp1*Ftmp8;
double Ftmp10 = 3*Ftmp9;
double Ftmp11 = pow(y, 2);
double Ftmp12 = pow(z, 2);
double Ftmp13 = pow(Rinv, 6);
double Ftmp14 = 30*x;
double Ftmp15 = Ftmp13*Ftmp14;
double Ftmp16 = pow(y, 3);
double Ftmp17 = Ftmp16*M[15];
double Ftmp18 = pow(z, 3);
double Ftmp19 = Ftmp18*M[18];
double Ftmp20 = Ftmp15*y;
double Ftmp21 = Ftmp12*Ftmp20;
double Ftmp22 = z*M[16];
double Ftmp23 = Ftmp13*Ftmp8;
double Ftmp24 = 105*Ftmp23;
double Ftmp25 = z*M[13];
double Ftmp26 = pow(Rinv, 8);
double Ftmp27 = 30*Ftmp23;
double Ftmp28 = Ftmp11*M[12];
double Ftmp29 = Ftmp12*Ftmp27;
double Ftmp30 = Ftmp8*y;
double Ftmp31 = Ftmp30*M[47];
double Ftmp32 = 1890*pow(Rinv, 10);
double Ftmp33 = Ftmp18*Ftmp32;
double Ftmp34 = z*M[45];
double Ftmp35 = 5*Ftmp9;
double Ftmp36 = Ftmp35 - 3;
double Ftmp37 = Ftmp1*Ftmp11;
double Ftmp38 = 5*Ftmp37;
double Ftmp39 = Ftmp38 - 1;
double Ftmp40 = Ftmp1*Ftmp12;
double Ftmp41 = 5*Ftmp40;
double Ftmp42 = Ftmp41 - 1;
double Ftmp43 = Ftmp9 - 1;
double Ftmp44 = Ftmp10 - 1;
double Ftmp45 = 3*Ftmp11;
double Ftmp46 = Ftmp1*Ftmp45;
double Ftmp47 = Ftmp46 - 1;
double Ftmp48 = 3*Ftmp12;
double Ftmp49 = Ftmp1*Ftmp48;
double Ftmp50 = Ftmp49 - 1;
double Ftmp51 = 15*Ftmp7;
double Ftmp52 = Ftmp51*y;
double Ftmp53 = 7*Ftmp9;
double Ftmp54 = Ftmp53 - 3;
double Ftmp55 = Ftmp54*M[20];
double Ftmp56 = 7*Ftmp37;
double Ftmp57 = Ftmp56 - 3;
double Ftmp58 = Ftmp57*M[25];
double Ftmp59 = 7*Ftmp40;
double Ftmp60 = Ftmp59 - 1;
double Ftmp61 = Ftmp60*M[27];
double Ftmp62 = Ftmp51*z;
double Ftmp63 = Ftmp54*M[21];
double Ftmp64 = Ftmp56 - 1;
double Ftmp65 = Ftmp64*M[26];
double Ftmp66 = Ftmp59 - 3;
double Ftmp67 = Ftmp62*Ftmp66;
double Ftmp68 = Ftmp43*Ftmp7;
double Ftmp69 = Ftmp14*Ftmp68;
double Ftmp70 = y*M[10];
double Ftmp71 = z*M[11];
double Ftmp72 = 30*Ftmp8;
double Ftmp73 = Ftmp52*x;
double Ftmp74 = Ftmp35 - 1;
double Ftmp75 = Ftmp74*M[10];
double Ftmp76 = Ftmp38 - 3;
double Ftmp77 = Ftmp76*M[15];
double Ftmp78 = Ftmp42*M[17];
double Ftmp79 = Ftmp62*x;
double Ftmp80 = Ftmp39*M[16];
double Ftmp81 = Ftmp41 - 3;
double Ftmp82 = Ftmp81*M[18];
double Ftmp83 = Ftmp36*M[9];
double Ftmp84 = Ftmp51*Ftmp8;
double Ftmp85 = Ftmp39*M[12];
double Ftmp86 = Ftmp42*M[14];
double Ftmp87 = Ftmp53 - 1;
double Ftmp88 = Ftmp26*x;
double Ftmp89 = 420*M[49];
double Ftmp90 = Ftmp16*Ftmp89*(9*Ftmp37 - 5);
double Ftmp91 = 420*M[54];
double Ftmp92 = Ftmp18*Ftmp91*(9*Ftmp40 - 5);
double Ftmp93 = Ftmp26*z;
double Ftmp94 = Ftmp30*Ftmp93;
double Ftmp95 = 1260*Ftmp50;
double Ftmp96 = Ftmp95*M[53];
double Ftmp97 = Ftmp88*y;
double Ftmp98 = Ftmp12*Ftmp97;
double Ftmp99 = 1260*M[50];
double Ftmp100 = Ftmp47*Ftmp99;
double Ftmp101 = Ftmp11*Ftmp88;
double Ftmp102 = Ftmp101*z;
double Ftmp103 = Ftmp44*M[38];
double Ftmp104 = 2835*Ftmp30;
double Ftmp105 = Ftmp34*Ftmp47;
double Ftmp106 = Ftmp26*Ftmp8;
double Ftmp107 = Ftmp106*Ftmp11;
double Ftmp108 = 1260*M[44];
double Ftmp109 = Ftmp108*Ftmp47;
double Ftmp110 = Ftmp106*Ftmp12;
double Ftmp111 = Ftmp95*M[48];
double Ftmp112 = Ftmp7*pow(x, 4);
double Ftmp113 = 63*Ftmp112 - 70*Ftmp9 + 15;
double Ftmp114 = 14*Ftmp37;
double Ftmp115 = -Ftmp114;
double Ftmp116 = pow(y, 4);
double Ftmp117 = 21*Ftmp7;
double Ftmp118 = Ftmp115 + Ftmp116*Ftmp117 + 1;
double Ftmp119 = -14*Ftmp40;
double Ftmp120 = pow(z, 4);
double Ftmp121 = Ftmp117*Ftmp120 + Ftmp119 + 1;
double Ftmp122 = -10*Ftmp9;
double Ftmp123 = Ftmp122 + 3;
double Ftmp124 = 35*Ftmp112 - 30*Ftmp9 + 3;
double Ftmp125 = 35*Ftmp7;
double Ftmp126 = Ftmp116*Ftmp125 - 30*Ftmp37 + 3;
double Ftmp127 = Ftmp120*Ftmp125 - 30*Ftmp40 + 3;
double Ftmp128 = 3*Ftmp112 - 4*Ftmp9 + 1;
double Ftmp129 = 1260*M[35];
double Ftmp130 = Ftmp13*x;
double Ftmp131 = Ftmp130*y;
double Ftmp132 = 14*Ftmp9;
double Ftmp133 = -Ftmp132;
double Ftmp134 = 21*Ftmp112 + Ftmp133 + 1;
double Ftmp135 = Ftmp134*M[35];
double Ftmp136 = 315*Ftmp131;
double Ftmp137 = 63*Ftmp7;
double Ftmp138 = Ftmp116*Ftmp137 - 70*Ftmp37 + 15;
double Ftmp139 = Ftmp138*M[49];
double Ftmp140 = 105*Ftmp131;
double Ftmp141 = Ftmp121*M[53];
double Ftmp142 = Ftmp130*z;
double Ftmp143 = 1260*M[36];
double Ftmp144 = Ftmp134*M[36];
double Ftmp145 = 315*Ftmp142;
double Ftmp146 = Ftmp118*M[50];
double Ftmp147 = Ftmp120*Ftmp137 - 70*Ftmp40 + 15;
double Ftmp148 = Ftmp147*M[54];
double Ftmp149 = 105*Ftmp142;
double Ftmp150 = 420*M[34];
double Ftmp151 = Ftmp113*M[34];
double Ftmp152 = Ftmp118*M[44];
double Ftmp153 = 315*Ftmp23;
double Ftmp154 = Ftmp121*M[48];
double Ftmp155 = Ftmp12*Ftmp37;
double Ftmp156 = -18*Ftmp155;
double Ftmp157 = -21*Ftmp37;
double Ftmp158 = Ftmp137*Ftmp8;
double Ftmp159 = Ftmp11*Ftmp158;
double Ftmp160 = -Ftmp53;
double Ftmp161 = Ftmp160 + 3;
double Ftmp162 = Ftmp157 + Ftmp159 + Ftmp161;
double Ftmp163 = -21*Ftmp40;
double Ftmp164 = Ftmp12*Ftmp158;
double Ftmp165 = Ftmp161 + Ftmp163 + Ftmp164;
double Ftmp166 = 8*Ftmp37;
double Ftmp167 = -Ftmp166;
double Ftmp168 = Ftmp11*Ftmp7;
double Ftmp169 = 14*Ftmp8;
double Ftmp170 = Ftmp168*Ftmp169;
double Ftmp171 = 1 - Ftmp9;
double Ftmp172 = -Ftmp38;
double Ftmp173 = Ftmp125*Ftmp8;
double Ftmp174 = 1 - Ftmp35;
double Ftmp175 = Ftmp11*Ftmp173 + Ftmp172 + Ftmp174;
double Ftmp176 = -8*Ftmp40;
double Ftmp177 = Ftmp12*Ftmp7;
double Ftmp178 = Ftmp169*Ftmp177;
double Ftmp179 = -Ftmp41;
double Ftmp180 = Ftmp12*Ftmp173 + Ftmp174 + Ftmp179;
double Ftmp181 = Ftmp11*Ftmp12;
double Ftmp182 = Ftmp125*Ftmp181 + Ftmp172 + Ftmp179 + 1;
double Ftmp183 = 18*Ftmp8;
double Ftmp184 = Ftmp168*Ftmp183;
double Ftmp185 = -Ftmp10;
double Ftmp186 = -10*Ftmp37;
double Ftmp187 = Ftmp186 + 3;
double Ftmp188 = 210*Ftmp131;
double Ftmp189 = -21*Ftmp9;
double Ftmp190 = -Ftmp56;
double Ftmp191 = Ftmp190 + 3;
double Ftmp192 = Ftmp159 + Ftmp189 + Ftmp191;
double Ftmp193 = Ftmp192*M[40];
double Ftmp194 = -10*Ftmp40;
double Ftmp195 = Ftmp177*Ftmp183;
double Ftmp196 = Ftmp171 + Ftmp195;
double Ftmp197 = Ftmp160 + 1;
double Ftmp198 = -Ftmp59;
double Ftmp199 = Ftmp164 + Ftmp198;
double Ftmp200 = Ftmp197 + Ftmp199;
double Ftmp201 = Ftmp200*M[42];
double Ftmp202 = Ftmp137*Ftmp181;
double Ftmp203 = Ftmp163 + Ftmp191 + Ftmp202;
double Ftmp204 = Ftmp203*M[51];
double Ftmp205 = Ftmp171 + Ftmp184;
double Ftmp206 = 210*Ftmp142;
double Ftmp207 = Ftmp159 + Ftmp190 + Ftmp197;
double Ftmp208 = Ftmp207*M[41];
double Ftmp209 = Ftmp194 + 3;
double Ftmp210 = Ftmp189 + Ftmp199 + 3;
double Ftmp211 = Ftmp210*M[43];
double Ftmp212 = Ftmp157 + Ftmp198 + Ftmp202 + 3;
double Ftmp213 = Ftmp212*M[52];
double Ftmp214 = -12*Ftmp37;
double Ftmp215 = 210*Ftmp23;
double Ftmp216 = Ftmp162*M[37];
double Ftmp217 = -12*Ftmp40;
double Ftmp218 = Ftmp165*M[39];
double Ftmp219 = Ftmp12*Ftmp168;
double Ftmp220 = 4*Ftmp40;
double Ftmp221 = 2*Ftmp220 - 2;
double Ftmp222 = Ftmp221*Ftmp37;
double Ftmp223 = 20*Ftmp219 + Ftmp222 + Ftmp39*Ftmp60;
double Ftmp224 = Ftmp223*M[46];
double Ftmp225 = 5*Ftmp60;
double Ftmp226 = 7*Ftmp39;
double Ftmp227 = pow(x, 3);
double Ftmp228 = Ftmp227*y;
double Ftmp229 = 30*Ftmp13;
double Ftmp230 = Ftmp229*M[9];
double Ftmp231 = Ftmp229*y;
double Ftmp232 = Ftmp11*Ftmp13;
double Ftmp233 = 105*Ftmp232;
double Ftmp234 = Ftmp11*M[10];
double Ftmp235 = Ftmp181*Ftmp229;
double Ftmp236 = Ftmp227*z;
double Ftmp237 = Ftmp32*M[38];
double Ftmp238 = Ftmp37 - 1;
double Ftmp239 = Ftmp51*x;
double Ftmp240 = Ftmp87*M[23];
double Ftmp241 = Ftmp57*M[30];
double Ftmp242 = Ftmp14*Ftmp7;
double Ftmp243 = 30*Ftmp238;
double Ftmp244 = Ftmp52*z;
double Ftmp245 = Ftmp51*Ftmp74;
double Ftmp246 = Ftmp11*Ftmp51;
double Ftmp247 = Ftmp150*Ftmp26*(9*Ftmp9 - 5);
double Ftmp248 = 2835*Ftmp102;
double Ftmp249 = Ftmp50*M[47];
double Ftmp250 = Ftmp143*Ftmp44;
double Ftmp251 = Ftmp129*Ftmp44;
double Ftmp252 = Ftmp181*Ftmp26;
double Ftmp253 = Ftmp116*Ftmp7;
double Ftmp254 = Ftmp12 + Ftmp8;
double Ftmp255 = 3*Ftmp253 - 4*Ftmp37 + 1;
double Ftmp256 = y*z;
double Ftmp257 = 315*Ftmp13*Ftmp256;
double Ftmp258 = Ftmp13*Ftmp256;
double Ftmp259 = 105*Ftmp258;
double Ftmp260 = 315*Ftmp232;
double Ftmp261 = 18*Ftmp9;
double Ftmp262 = -Ftmp12*Ftmp261;
double Ftmp263 = 3*Ftmp8;
double Ftmp264 = -8*Ftmp9;
double Ftmp265 = 1 - Ftmp37;
double Ftmp266 = 14*Ftmp219;
double Ftmp267 = -Ftmp46;
double Ftmp268 = Ftmp184 + Ftmp265;
double Ftmp269 = 210*Ftmp258;
double Ftmp270 = 18*Ftmp219;
double Ftmp271 = -12*Ftmp9;
double Ftmp272 = 210*Ftmp232;
double Ftmp273 = Ftmp229*z;
double Ftmp274 = Ftmp12*Ftmp13;
double Ftmp275 = 105*Ftmp274;
double Ftmp276 = Ftmp15*z;
double Ftmp277 = Ftmp12*M[45];
double Ftmp278 = Ftmp40 - 1;
double Ftmp279 = Ftmp278*z;
double Ftmp280 = Ftmp12*Ftmp51;
double Ftmp281 = 2835*Ftmp98;
double Ftmp282 = Ftmp120*Ftmp7;
double Ftmp283 = Ftmp11 + Ftmp8;
double Ftmp284 = -1260*Ftmp220 + 3780*Ftmp282 + 1260;
double Ftmp285 = 315*Ftmp274;
double Ftmp286 = -Ftmp11*Ftmp261;
double Ftmp287 = 1 - Ftmp40;
double Ftmp288 = -Ftmp49;
double Ftmp289 = Ftmp195 + Ftmp287;
double Ftmp290 = 210*Ftmp274;
double Ftmp291 = -Ftmp278;
#pragma omp atomic
F[0] += Ftmp0*(3*Ftmp1*Ftmp36*M[9] + 3*Ftmp1*Ftmp39*M[12] + 3*Ftmp1*Ftmp42*M[14] + 6*Ftmp1*Ftmp43*x*M[3] + 3*Ftmp1*Ftmp44*x*M[3] + 3*Ftmp1*Ftmp47*x*M[6] + 3*Ftmp1*Ftmp50*x*M[8] - Ftmp10*M[0] - Ftmp100*Ftmp102 - Ftmp103*Ftmp104*Ftmp93 - Ftmp104*Ftmp105*Ftmp26 - Ftmp107*Ftmp109 + 60*Ftmp11*Ftmp13*Ftmp57*x*M[29] - Ftmp11*Ftmp15*Ftmp22 + 210*Ftmp11*Ftmp26*Ftmp8*z*M[26] + 6*Ftmp11*Ftmp7*x*M[6] - Ftmp110*Ftmp111 + 15*Ftmp113*Ftmp7*M[34] + 45*Ftmp118*Ftmp7*M[44] + 60*Ftmp12*Ftmp13*Ftmp66*x*M[33] + 210*Ftmp12*Ftmp26*Ftmp8*y*M[27] + 6*Ftmp12*Ftmp7*x*M[8] + 45*Ftmp121*Ftmp7*M[48] + 15*Ftmp124*Ftmp7*x*M[19] + 15*Ftmp126*Ftmp7*x*M[29] + 15*Ftmp127*Ftmp7*x*M[33] - Ftmp128*Ftmp129*Ftmp131 - Ftmp128*Ftmp142*Ftmp143 + 210*Ftmp13*Ftmp43*Ftmp8*y*M[20] + 210*Ftmp13*Ftmp43*Ftmp8*z*M[21] + 210*Ftmp13*Ftmp43*x*y*z*M[23] + 315*Ftmp13*Ftmp44*y*z*M[38] + 315*Ftmp13*Ftmp47*y*z*M[45] + 315*Ftmp13*Ftmp50*y*z*M[47] + 105*Ftmp13*Ftmp54*Ftmp8*y*M[20] + 105*Ftmp13*Ftmp54*Ftmp8*z*M[21] + 105*Ftmp13*Ftmp57*Ftmp8*y*M[25] + 105*Ftmp13*Ftmp57*x*y*z*M[30] + 105*Ftmp13*Ftmp60*Ftmp8*y*M[27] + 105*Ftmp13*Ftmp64*Ftmp8*z*M[26] + 105*Ftmp13*Ftmp66*Ftmp8*z*M[28] + 105*Ftmp13*Ftmp66*x*y*z*M[32] + 105*Ftmp13*Ftmp87*x*y*z*M[23] - Ftmp135*Ftmp136 - Ftmp136*Ftmp141 - Ftmp139*Ftmp140 - Ftmp140*Ftmp193 - Ftmp140*Ftmp201 - Ftmp140*Ftmp204 - Ftmp144*Ftmp145 - Ftmp145*Ftmp146 - Ftmp148*Ftmp149 - Ftmp149*Ftmp208 - Ftmp149*Ftmp211 - Ftmp149*Ftmp213 - Ftmp15*Ftmp17 - Ftmp15*Ftmp19 - Ftmp15*(Ftmp11 - Ftmp114*Ftmp12 + Ftmp12)*M[31] - Ftmp150*Ftmp23*(9*Ftmp112 + Ftmp133 + 5) - Ftmp151*Ftmp24 - Ftmp152*Ftmp153 - Ftmp153*Ftmp154 + 210*Ftmp16*Ftmp26*Ftmp8*M[25] + 210*Ftmp16*Ftmp26*x*z*M[30] - Ftmp16*Ftmp32*Ftmp34*Ftmp8 + 15*Ftmp162*Ftmp7*M[37] + 15*Ftmp165*Ftmp7*M[39] + 15*Ftmp175*Ftmp7*x*M[22] + 210*Ftmp18*Ftmp26*Ftmp8*M[28] + 210*Ftmp18*Ftmp26*x*y*M[32] + 15*Ftmp180*Ftmp7*x*M[24] + 15*Ftmp182*Ftmp7*x*M[31] - Ftmp188*(Ftmp194 + Ftmp196)*M[42] - Ftmp188*(Ftmp184 + Ftmp185 + Ftmp187)*M[40] - Ftmp2*y*M[4] - Ftmp206*(Ftmp186 + Ftmp205)*M[41] - Ftmp206*(Ftmp185 + Ftmp195 + Ftmp209)*M[43] - Ftmp21*M[17] - Ftmp215*(Ftmp196 + Ftmp217)*M[39] - Ftmp215*(Ftmp205 + Ftmp214)*M[37] - Ftmp216*Ftmp24 - Ftmp218*Ftmp24 + 15*Ftmp223*Ftmp7*M[46] - Ftmp224*Ftmp24 - Ftmp24*Ftmp25*y - Ftmp26*Ftmp72*(Ftmp11*Ftmp221 + Ftmp11*Ftmp225 + Ftmp12*Ftmp226 + 48*Ftmp155)*M[46] + 210*Ftmp26*x*y*(Ftmp11 + Ftmp156 + Ftmp48)*M[51] + 210*Ftmp26*x*z*(Ftmp12 + Ftmp156 + Ftmp45)*M[52] - Ftmp27*Ftmp28 - Ftmp29*M[14] - Ftmp3*M[5] - Ftmp31*Ftmp33 - 2835*Ftmp31*Ftmp50*Ftmp93 - 1890*Ftmp43*Ftmp94*M[38] - Ftmp5*M[1] - Ftmp52*Ftmp55 - Ftmp52*Ftmp58 - Ftmp52*Ftmp61 - Ftmp6*x - Ftmp62*Ftmp63 - Ftmp62*Ftmp65 - Ftmp67*M[28] - Ftmp68*Ftmp72*M[9] - Ftmp69*Ftmp70 - Ftmp69*Ftmp71 + 15*Ftmp7*Ftmp8*y*M[4] + 15*Ftmp7*Ftmp8*z*M[5] + 15*Ftmp7*x*y*z*M[7] + 60*Ftmp7*x*(7*Ftmp112 + Ftmp123)*M[19] + 30*Ftmp7*x*(Ftmp167 + Ftmp170 + Ftmp171)*M[22] + 30*Ftmp7*x*(Ftmp171 + Ftmp176 + Ftmp178)*M[24] + 15*Ftmp7*y*z*M[13] - Ftmp73*Ftmp75 - Ftmp73*Ftmp77 - Ftmp73*Ftmp78 - Ftmp74*Ftmp79*M[11] - Ftmp79*Ftmp80 - Ftmp79*Ftmp82 - Ftmp83*Ftmp84 - Ftmp84*Ftmp85 - Ftmp84*Ftmp86 - Ftmp88*Ftmp90 - Ftmp88*Ftmp92 - Ftmp96*Ftmp98 + M[0]);
#pragma omp atomic
F[1] += Ftmp0*(6*Ftmp1*Ftmp238*y*M[6] + 3*Ftmp1*Ftmp42*M[17] + 3*Ftmp1*Ftmp44*y*M[3] + 3*Ftmp1*Ftmp47*y*M[6] + 3*Ftmp1*Ftmp50*y*M[8] + 3*Ftmp1*Ftmp74*M[10] + 3*Ftmp1*Ftmp76*M[15] - 2835*Ftmp101*Ftmp105 - 1890*Ftmp101*Ftmp238*Ftmp34 - Ftmp103*Ftmp248 - Ftmp107*Ftmp251 - Ftmp108*Ftmp131*Ftmp255 + 210*Ftmp11*Ftmp12*Ftmp26*x*M[27] + 210*Ftmp11*Ftmp13*Ftmp238*x*M[25] + 210*Ftmp11*Ftmp13*Ftmp238*z*M[30] + 105*Ftmp11*Ftmp13*Ftmp54*x*M[20] + 105*Ftmp11*Ftmp13*Ftmp57*x*M[25] + 105*Ftmp11*Ftmp13*Ftmp57*z*M[30] + 105*Ftmp11*Ftmp13*Ftmp60*x*M[27] + 105*Ftmp11*Ftmp13*Ftmp66*z*M[32] + 105*Ftmp11*Ftmp13*Ftmp87*z*M[23] + 210*Ftmp11*Ftmp18*Ftmp26*M[32] + 210*Ftmp11*Ftmp227*Ftmp26*M[20] - Ftmp11*Ftmp236*Ftmp237 + 210*Ftmp11*Ftmp26*Ftmp8*z*M[23] + 210*Ftmp11*Ftmp26*(Ftmp254 + Ftmp262)*M[42] - Ftmp11*Ftmp33*x*M[47] + 15*Ftmp11*Ftmp7*x*M[4] + 15*Ftmp11*Ftmp7*z*M[7] - Ftmp111*Ftmp98 + 60*Ftmp12*Ftmp13*Ftmp66*y*M[33] + 6*Ftmp12*Ftmp7*y*M[8] + 45*Ftmp121*Ftmp7*M[53] + 15*Ftmp124*Ftmp7*y*M[19] + 15*Ftmp126*Ftmp7*y*M[29] + 15*Ftmp127*Ftmp7*y*M[33] + 210*Ftmp13*Ftmp238*x*y*z*M[26] + 315*Ftmp13*Ftmp44*x*z*M[38] + 315*Ftmp13*Ftmp47*x*z*M[45] + 315*Ftmp13*Ftmp50*x*z*M[47] + 60*Ftmp13*Ftmp54*Ftmp8*y*M[19] + 105*Ftmp13*Ftmp54*x*y*z*M[21] + 105*Ftmp13*Ftmp64*x*y*z*M[26] + 105*Ftmp13*Ftmp66*x*y*z*M[28] + 45*Ftmp134*Ftmp7*M[35] - Ftmp135*Ftmp260 - Ftmp136*Ftmp152 - Ftmp136*Ftmp154 + 15*Ftmp138*Ftmp7*M[49] - Ftmp139*Ftmp233 - Ftmp140*Ftmp151 - Ftmp140*Ftmp216 - Ftmp140*Ftmp218 - Ftmp140*Ftmp224 - Ftmp141*Ftmp260 - Ftmp144*Ftmp257 - Ftmp146*Ftmp257 - Ftmp148*Ftmp259 - Ftmp168*Ftmp243*M[15] + 15*Ftmp175*Ftmp7*y*M[22] + 210*Ftmp18*Ftmp26*x*y*M[28] + 15*Ftmp180*Ftmp7*y*M[24] + 15*Ftmp182*Ftmp7*y*M[31] - Ftmp188*(Ftmp123 + Ftmp184 + Ftmp267)*M[37] - Ftmp19*Ftmp231 + 15*Ftmp192*Ftmp7*M[40] - Ftmp193*Ftmp233 - Ftmp20*(48*Ftmp219 + Ftmp222 + Ftmp225*Ftmp238 + Ftmp39*Ftmp59 - 28*Ftmp40 + 2)*M[46] + 15*Ftmp200*Ftmp7*M[42] - Ftmp201*Ftmp233 + 15*Ftmp203*Ftmp7*M[51] - Ftmp204*Ftmp233 - Ftmp208*Ftmp259 - Ftmp21*M[14] - Ftmp211*Ftmp259 - Ftmp213*Ftmp259 - Ftmp22*Ftmp243*Ftmp7*y + 210*Ftmp227*Ftmp26*y*z*M[21] - Ftmp228*Ftmp230 - Ftmp228*Ftmp247 - Ftmp231*(-Ftmp12*Ftmp132 + Ftmp254)*M[24] - Ftmp232*Ftmp89*(Ftmp115 + 9*Ftmp253 + 5) - Ftmp233*Ftmp25*x - Ftmp234*Ftmp245 - Ftmp234*Ftmp27 - Ftmp235*M[17] - Ftmp238*Ftmp242*y*M[12] - Ftmp239*Ftmp55 - Ftmp239*Ftmp58 - Ftmp239*Ftmp61 - Ftmp240*Ftmp62 - Ftmp241*Ftmp62 - Ftmp244*Ftmp80 - Ftmp244*Ftmp82 - Ftmp246*Ftmp77 - Ftmp246*Ftmp78 - Ftmp248*Ftmp249 - Ftmp250*Ftmp94 - Ftmp252*Ftmp96 - Ftmp255*Ftmp258*Ftmp99 - Ftmp26*Ftmp92*y + 210*Ftmp26*x*y*(Ftmp262 + Ftmp48 + Ftmp8)*M[39] + 210*Ftmp26*y*z*(Ftmp12 + Ftmp262 + Ftmp263)*M[43] - Ftmp269*(Ftmp122 + Ftmp268)*M[41] - Ftmp269*(Ftmp209 + Ftmp267 + Ftmp270)*M[52] - Ftmp27*Ftmp71*y - Ftmp272*(Ftmp268 + Ftmp271)*M[40] - Ftmp272*(Ftmp217 + Ftmp265 + Ftmp270)*M[51] - Ftmp3*M[7] - Ftmp4*M[4] - Ftmp46*M[1] - Ftmp5*M[0] - Ftmp52*Ftmp71*Ftmp74 - Ftmp6*y - Ftmp67*M[32] + 6*Ftmp7*Ftmp8*y*M[3] + 15*Ftmp7*x*y*z*M[5] + 15*Ftmp7*x*z*M[13] + 60*Ftmp7*y*(Ftmp187 + 7*Ftmp253)*M[29] + 30*Ftmp7*y*(Ftmp170 + Ftmp264 + Ftmp265)*M[22] + 30*Ftmp7*y*(Ftmp176 + Ftmp265 + Ftmp266)*M[31] - Ftmp73*Ftmp83 - Ftmp73*Ftmp85 - Ftmp73*Ftmp86 + M[1]);
#pragma omp atomic
F[2] += Ftmp0*(6*Ftmp1*Ftmp278*z*M[8] + 3*Ftmp1*Ftmp39*M[16] + 3*Ftmp1*Ftmp44*z*M[3] + 3*Ftmp1*Ftmp47*z*M[6] + 3*Ftmp1*Ftmp50*z*M[8] + 3*Ftmp1*Ftmp74*M[11] + 3*Ftmp1*Ftmp81*M[18] - Ftmp100*Ftmp252 - Ftmp102*Ftmp109 - Ftmp103*Ftmp281 + 210*Ftmp11*Ftmp12*Ftmp26*x*M[26] + 60*Ftmp11*Ftmp13*Ftmp57*z*M[29] + 6*Ftmp11*Ftmp7*z*M[6] - Ftmp110*Ftmp250 + 45*Ftmp118*Ftmp7*M[50] + 210*Ftmp12*Ftmp13*Ftmp278*x*M[28] + 210*Ftmp12*Ftmp13*Ftmp278*y*M[32] + 105*Ftmp12*Ftmp13*Ftmp54*x*M[21] + 105*Ftmp12*Ftmp13*Ftmp57*y*M[30] + 105*Ftmp12*Ftmp13*Ftmp64*x*M[26] + 105*Ftmp12*Ftmp13*Ftmp66*x*M[28] + 105*Ftmp12*Ftmp13*Ftmp66*y*M[32] + 105*Ftmp12*Ftmp13*Ftmp87*y*M[23] + 210*Ftmp12*Ftmp16*Ftmp26*M[30] + 210*Ftmp12*Ftmp227*Ftmp26*M[21] - Ftmp12*Ftmp228*Ftmp237 - Ftmp12*Ftmp245*M[11] + 210*Ftmp12*Ftmp26*Ftmp8*y*M[23] + 210*Ftmp12*Ftmp26*(Ftmp283 + Ftmp286)*M[41] + 15*Ftmp12*Ftmp7*x*M[5] + 15*Ftmp12*Ftmp7*y*M[7] + 15*Ftmp124*Ftmp7*z*M[19] + 15*Ftmp126*Ftmp7*z*M[29] + 15*Ftmp127*Ftmp7*z*M[33] + 210*Ftmp13*Ftmp278*x*y*z*M[27] + 315*Ftmp13*Ftmp44*x*y*M[38] + 315*Ftmp13*Ftmp47*x*y*M[45] + 315*Ftmp13*Ftmp50*x*y*M[47] + 60*Ftmp13*Ftmp54*Ftmp8*z*M[19] + 105*Ftmp13*Ftmp54*x*y*z*M[20] + 105*Ftmp13*Ftmp57*x*y*z*M[25] + 105*Ftmp13*Ftmp60*x*y*z*M[27] + 45*Ftmp134*Ftmp7*M[36] - Ftmp135*Ftmp257 - Ftmp139*Ftmp259 - Ftmp141*Ftmp257 - Ftmp142*Ftmp284*M[48] - Ftmp144*Ftmp285 - Ftmp145*Ftmp152 - Ftmp145*Ftmp154 - Ftmp146*Ftmp285 + 15*Ftmp147*Ftmp7*M[54] - Ftmp148*Ftmp275 - Ftmp149*Ftmp151 - Ftmp149*Ftmp216 - Ftmp149*Ftmp218 - Ftmp149*Ftmp224 + 210*Ftmp16*Ftmp26*x*z*M[25] - Ftmp16*Ftmp277*Ftmp32*x - Ftmp17*Ftmp273 + 15*Ftmp175*Ftmp7*z*M[22] - 30*Ftmp177*Ftmp278*M[18] + 15*Ftmp180*Ftmp7*z*M[24] + 15*Ftmp182*Ftmp7*z*M[31] - Ftmp193*Ftmp259 - Ftmp2*y*M[7] - Ftmp201*Ftmp259 - Ftmp204*Ftmp259 - Ftmp206*(Ftmp123 + Ftmp195 + Ftmp288)*M[39] + 15*Ftmp207*Ftmp7*M[41] - Ftmp208*Ftmp275 + 15*Ftmp210*Ftmp7*M[43] - Ftmp211*Ftmp275 + 15*Ftmp212*Ftmp7*M[52] - Ftmp213*Ftmp275 + 210*Ftmp227*Ftmp26*y*z*M[20] - Ftmp230*Ftmp236 - Ftmp235*M[16] - Ftmp236*Ftmp247 - Ftmp239*Ftmp63 - Ftmp239*Ftmp65 - Ftmp239*Ftmp66*M[28] - Ftmp240*Ftmp52 - Ftmp241*Ftmp52 - Ftmp242*Ftmp279*M[14] - Ftmp244*Ftmp75 - Ftmp244*Ftmp77 - Ftmp244*Ftmp78 - Ftmp249*Ftmp281 - Ftmp251*Ftmp94 - Ftmp258*Ftmp284*M[53] + 210*Ftmp26*x*z*(Ftmp286 + Ftmp45 + Ftmp8)*M[37] + 210*Ftmp26*y*z*(Ftmp11 + Ftmp263 + Ftmp286)*M[40] - Ftmp269*(Ftmp122 + Ftmp289)*M[42] - Ftmp269*(Ftmp187 + Ftmp270 + Ftmp288)*M[51] - Ftmp27*Ftmp70*z - Ftmp273*(-Ftmp11*Ftmp132 + Ftmp283)*M[22] - Ftmp274*Ftmp91*(Ftmp119 + 9*Ftmp282 + 5) - Ftmp275*x*y*M[13] - Ftmp276*Ftmp28 - Ftmp276*(-Ftmp166*Ftmp291 + 40*Ftmp219 + Ftmp222 - Ftmp226*Ftmp291 - 20*Ftmp37 + Ftmp38*Ftmp60)*M[46] - 2835*Ftmp277*Ftmp47*Ftmp97 - 1890*Ftmp278*Ftmp98*M[47] - 30*Ftmp279*Ftmp7*y*M[17] - Ftmp280*Ftmp80 - Ftmp280*Ftmp82 - Ftmp29*M[11] - Ftmp290*(Ftmp271 + Ftmp289)*M[43] - Ftmp290*(Ftmp214 + Ftmp270 + Ftmp287)*M[52] - Ftmp3*x*M[0] - Ftmp3*y*M[1] - Ftmp4*M[5] - Ftmp49*M[2] - Ftmp52*Ftmp66*M[32] + 6*Ftmp7*Ftmp8*z*M[3] + 15*Ftmp7*x*y*z*M[4] + 15*Ftmp7*x*y*M[13] + 60*Ftmp7*z*(Ftmp209 + 7*Ftmp282)*M[33] + 30*Ftmp7*z*(Ftmp167 + Ftmp266 + Ftmp287)*M[31] + 30*Ftmp7*z*(Ftmp178 + Ftmp264 + Ftmp287)*M[24] - Ftmp79*Ftmp83 - Ftmp79*Ftmp85 - Ftmp79*Ftmp86 - Ftmp90*Ftmp93 + M[2]);

}

void S2Mc_5(double x, double y, double z, double * S, double * M) {
double Mtmp0 = x*S[0];
double Mtmp1 = z*S[2];
double Mtmp2 = -Mtmp1;
double Mtmp3 = x*S[1];
double Mtmp4 = y*S[0];
double Mtmp5 = x*S[2];
double Mtmp6 = z*S[0];
double Mtmp7 = y*S[1];
double Mtmp8 = y*S[2];
double Mtmp9 = z*S[1];
double Mtmp10 = Mtmp1*x;
double Mtmp11 = pow(x, 2);
double Mtmp12 = pow(z, 2);
double Mtmp13 = (1.0/2.0)*S[0];
double Mtmp14 = Mtmp12*Mtmp13;
double Mtmp15 = Mtmp0*y;
double Mtmp16 = Mtmp1*y;
double Mtmp17 = (1.0/2.0)*S[1];
double Mtmp18 = Mtmp12*Mtmp17;
double Mtmp19 = Mtmp0*z;
double Mtmp20 = (1.0/2.0)*S[2];
double Mtmp21 = -Mtmp12*Mtmp20;
double Mtmp22 = Mtmp3*y;
double Mtmp23 = pow(y, 2);
double Mtmp24 = Mtmp5*y;
double Mtmp25 = Mtmp3*z;
double Mtmp26 = Mtmp4*z;
double Mtmp27 = Mtmp7*z;
double Mtmp28 = pow(x, 3);
double Mtmp29 = pow(z, 3);
double Mtmp30 = Mtmp29*S[2];
double Mtmp31 = 3*Mtmp12;
double Mtmp32 = Mtmp0*Mtmp31;
double Mtmp33 = 3*Mtmp11;
double Mtmp34 = Mtmp1*Mtmp33;
double Mtmp35 = (1.0/2.0)*Mtmp12;
double Mtmp36 = Mtmp10*y + (1.0/2.0)*Mtmp12*Mtmp4 + Mtmp3*Mtmp35;
double Mtmp37 = Mtmp29*S[0];
double Mtmp38 = Mtmp31*Mtmp7;
double Mtmp39 = 3*Mtmp23;
double Mtmp40 = Mtmp1*Mtmp39;
double Mtmp41 = Mtmp29*S[1];
double Mtmp42 = (1.0/2.0)*Mtmp11;
double Mtmp43 = pow(y, 3);
double Mtmp44 = (1.0/2.0)*Mtmp23;
double Mtmp45 = pow(x, 4);
double Mtmp46 = Mtmp11*S[0];
double Mtmp47 = 6*Mtmp12;
double Mtmp48 = 4*Mtmp28;
double Mtmp49 = pow(z, 4);
double Mtmp50 = 4*Mtmp29;
double Mtmp51 = Mtmp49*S[0] + Mtmp5*Mtmp50;
double Mtmp52 = Mtmp11*Mtmp47;
double Mtmp53 = 12*Mtmp12;
double Mtmp54 = Mtmp11*Mtmp16;
double Mtmp55 = Mtmp49*S[1] + Mtmp50*Mtmp8;
double Mtmp56 = Mtmp49*S[2];
double Mtmp57 = 6*Mtmp23;
double Mtmp58 = Mtmp23*S[0];
double Mtmp59 = -Mtmp24*Mtmp31 - Mtmp29*Mtmp3 - Mtmp29*Mtmp4;
double Mtmp60 = Mtmp31*S[1];
double Mtmp61 = 2*Mtmp29;
double Mtmp62 = Mtmp31*S[2];
double Mtmp63 = pow(y, 4);
double Mtmp64 = 4*Mtmp43;
double Mtmp65 = Mtmp23*Mtmp47;
#pragma omp atomic
M[0] += S[0];
#pragma omp atomic
M[1] += S[1];
#pragma omp atomic
M[2] += S[2];
#pragma omp atomic
M[3] += Mtmp0 + Mtmp2;
#pragma omp atomic
M[4] += Mtmp3 + Mtmp4;
#pragma omp atomic
M[5] += Mtmp5 + Mtmp6;
#pragma omp atomic
M[6] += Mtmp2 + Mtmp7;
#pragma omp atomic
M[7] += Mtmp8 + Mtmp9;
#pragma omp atomic
M[8] += -Mtmp10 + (1.0/2.0)*Mtmp11*S[0] - Mtmp14;
#pragma omp atomic
M[9] += Mtmp11*Mtmp17 + Mtmp15 - Mtmp16 - Mtmp18;
#pragma omp atomic
M[10] += Mtmp11*Mtmp20 + Mtmp19 + Mtmp21;
#pragma omp atomic
M[11] += -Mtmp10 + Mtmp13*Mtmp23 - Mtmp14 + Mtmp22;
#pragma omp atomic
M[12] += Mtmp24 + Mtmp25 + Mtmp26;
#pragma omp atomic
M[13] += -Mtmp16 - Mtmp18 + (1.0/2.0)*Mtmp23*S[1];
#pragma omp atomic
M[14] += Mtmp20*Mtmp23 + Mtmp21 + Mtmp27;
#pragma omp atomic
M[15] += (1.0/6.0)*Mtmp28*S[0] + (1.0/6.0)*Mtmp30 - 1.0/6.0*Mtmp32 - 1.0/6.0*Mtmp34;
#pragma omp atomic
M[16] += (1.0/2.0)*Mtmp11*y*S[0] + (1.0/6.0)*Mtmp28*S[1] - Mtmp36;
#pragma omp atomic
M[17] += (1.0/6.0)*Mtmp28*S[2] - 1.0/6.0*Mtmp31*Mtmp5 + (1.0/6.0)*Mtmp33*Mtmp6 - 1.0/6.0*Mtmp37;
#pragma omp atomic
M[18] += (1.0/2.0)*Mtmp11*y*S[1] + (1.0/2.0)*Mtmp23*x*S[0] + (1.0/3.0)*Mtmp29*S[2] - 1.0/6.0*Mtmp32 - 1.0/6.0*Mtmp34 - 1.0/6.0*Mtmp38 - 1.0/6.0*Mtmp40;
#pragma omp atomic
M[19] += Mtmp15*z - Mtmp35*Mtmp8 - 1.0/6.0*Mtmp41 + Mtmp42*Mtmp8 + Mtmp42*Mtmp9;
#pragma omp atomic
M[20] += (1.0/2.0)*Mtmp23*x*S[1] - Mtmp36 + (1.0/6.0)*Mtmp43*S[0];
#pragma omp atomic
M[21] += Mtmp22*z - Mtmp35*Mtmp5 - 1.0/6.0*Mtmp37 + Mtmp44*Mtmp5 + Mtmp44*Mtmp6;
#pragma omp atomic
M[22] += (1.0/6.0)*Mtmp30 - 1.0/6.0*Mtmp38 - 1.0/6.0*Mtmp40 + (1.0/6.0)*Mtmp43*S[1];
#pragma omp atomic
M[23] += -1.0/6.0*Mtmp31*Mtmp8 + (1.0/6.0)*Mtmp39*Mtmp9 - 1.0/6.0*Mtmp41 + (1.0/6.0)*Mtmp43*S[2];
#pragma omp atomic
M[24] += -1.0/24.0*Mtmp1*Mtmp48 + (1.0/24.0)*Mtmp45*S[0] - 1.0/24.0*Mtmp46*Mtmp47 + (1.0/24.0)*Mtmp51;
#pragma omp atomic
M[25] += -1.0/24.0*Mtmp15*Mtmp53 + (1.0/24.0)*Mtmp4*Mtmp48 + (1.0/24.0)*Mtmp45*S[1] - 1.0/24.0*Mtmp52*S[1] - 1.0/2.0*Mtmp54 + (1.0/24.0)*Mtmp55;
#pragma omp atomic
M[26] += -1.0/24.0*Mtmp0*Mtmp50 + (1.0/24.0)*Mtmp45*S[2] + (1.0/24.0)*Mtmp48*Mtmp6 - 1.0/24.0*Mtmp52*S[2] + (1.0/24.0)*Mtmp56;
#pragma omp atomic
M[27] += -1.0/6.0*Mtmp1*Mtmp28 - 1.0/12.0*Mtmp10*Mtmp57 + (1.0/4.0)*Mtmp11*Mtmp23*S[0] - 1.0/12.0*Mtmp22*Mtmp47 + (1.0/6.0)*Mtmp28*y*S[1] + (1.0/3.0)*Mtmp29*x*S[2] - 1.0/12.0*Mtmp31*Mtmp46 - 1.0/12.0*Mtmp31*Mtmp58 + (1.0/12.0)*Mtmp49*S[0];
#pragma omp atomic
M[28] += (1.0/6.0)*Mtmp26*Mtmp33 + (1.0/6.0)*Mtmp28*Mtmp8 + (1.0/6.0)*Mtmp28*Mtmp9 + (1.0/6.0)*Mtmp59;
#pragma omp atomic
M[29] += -1.0/6.0*Mtmp1*Mtmp43 + (1.0/4.0)*Mtmp11*Mtmp23*S[1] - 1.0/12.0*Mtmp11*Mtmp60 - 1.0/12.0*Mtmp15*Mtmp47 - 1.0/12.0*Mtmp23*Mtmp60 + (1.0/3.0)*Mtmp29*y*S[2] + (1.0/6.0)*Mtmp43*x*S[0] + (1.0/12.0)*Mtmp49*S[1] - 1.0/2.0*Mtmp54;
#pragma omp atomic
M[30] += -1.0/12.0*Mtmp0*Mtmp61 + (1.0/2.0)*Mtmp11*Mtmp27 - 1.0/12.0*Mtmp11*Mtmp62 + (1.0/12.0)*Mtmp19*Mtmp57 + (1.0/12.0)*Mtmp23*Mtmp33*S[2] - 1.0/12.0*Mtmp23*Mtmp62 + (1.0/12.0)*Mtmp56 - 1.0/12.0*Mtmp61*Mtmp7;
#pragma omp atomic
M[31] += -1.0/2.0*Mtmp10*Mtmp23 - 1.0/24.0*Mtmp22*Mtmp53 + (1.0/24.0)*Mtmp3*Mtmp64 - 1.0/24.0*Mtmp47*Mtmp58 + (1.0/24.0)*Mtmp51 + (1.0/24.0)*Mtmp63*S[0];
#pragma omp atomic
M[32] += (1.0/6.0)*Mtmp25*Mtmp39 + (1.0/6.0)*Mtmp43*Mtmp5 + (1.0/6.0)*Mtmp43*Mtmp6 + (1.0/6.0)*Mtmp59;
#pragma omp atomic
M[33] += -1.0/24.0*Mtmp1*Mtmp64 + (1.0/24.0)*Mtmp55 + (1.0/24.0)*Mtmp63*S[1] - 1.0/24.0*Mtmp65*S[1];
#pragma omp atomic
M[34] += -1.0/24.0*Mtmp50*Mtmp7 + (1.0/24.0)*Mtmp56 + (1.0/24.0)*Mtmp63*S[2] + (1.0/24.0)*Mtmp64*Mtmp9 - 1.0/24.0*Mtmp65*S[2];

}

void M2Mc_5(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = z*M[2];
double Mstmp2 = -Mstmp1;
double Mstmp3 = x*M[1];
double Mstmp4 = y*M[0];
double Mstmp5 = x*M[2];
double Mstmp6 = z*M[0];
double Mstmp7 = y*M[1];
double Mstmp8 = y*M[2];
double Mstmp9 = z*M[1];
double Mstmp10 = x*M[3];
double Mstmp11 = pow(x, 2);
double Mstmp12 = (1.0/2.0)*M[0];
double Mstmp13 = z*M[5];
double Mstmp14 = pow(z, 2);
double Mstmp15 = Mstmp1*x;
double Mstmp16 = -Mstmp12*Mstmp14 - Mstmp13 - Mstmp15;
double Mstmp17 = x*M[4];
double Mstmp18 = y*M[3];
double Mstmp19 = (1.0/2.0)*M[1];
double Mstmp20 = Mstmp0*y;
double Mstmp21 = z*M[7];
double Mstmp22 = Mstmp1*y;
double Mstmp23 = -Mstmp14*Mstmp19 - Mstmp21 - Mstmp22;
double Mstmp24 = x*M[5];
double Mstmp25 = z*M[3];
double Mstmp26 = Mstmp0*z;
double Mstmp27 = (1.0/2.0)*M[2];
double Mstmp28 = -Mstmp14*Mstmp27;
double Mstmp29 = x*M[6];
double Mstmp30 = y*M[4];
double Mstmp31 = pow(y, 2);
double Mstmp32 = Mstmp3*y;
double Mstmp33 = x*M[7];
double Mstmp34 = y*M[5];
double Mstmp35 = z*M[4];
double Mstmp36 = Mstmp5*y;
double Mstmp37 = Mstmp3*z;
double Mstmp38 = Mstmp4*z;
double Mstmp39 = y*M[6];
double Mstmp40 = y*M[7];
double Mstmp41 = z*M[6];
double Mstmp42 = Mstmp7*z;
double Mstmp43 = x*M[8];
double Mstmp44 = z*M[10];
double Mstmp45 = Mstmp13*x;
double Mstmp46 = (1.0/2.0)*M[3];
double Mstmp47 = pow(x, 3);
double Mstmp48 = (1.0/6.0)*Mstmp47;
double Mstmp49 = Mstmp14*Mstmp46;
double Mstmp50 = pow(z, 3);
double Mstmp51 = (1.0/6.0)*Mstmp50;
double Mstmp52 = Mstmp51*M[2];
double Mstmp53 = (1.0/2.0)*Mstmp14;
double Mstmp54 = Mstmp0*Mstmp53;
double Mstmp55 = (1.0/2.0)*Mstmp11;
double Mstmp56 = Mstmp1*Mstmp55;
double Mstmp57 = x*M[9];
double Mstmp58 = y*M[8];
double Mstmp59 = (1.0/2.0)*M[4];
double Mstmp60 = Mstmp10*y;
double Mstmp61 = z*M[12];
double Mstmp62 = Mstmp21*x;
double Mstmp63 = -Mstmp13*y - Mstmp14*Mstmp59 - Mstmp15*y - Mstmp3*Mstmp53 - Mstmp4*Mstmp53 - Mstmp61 - Mstmp62;
double Mstmp64 = x*M[10];
double Mstmp65 = (1.0/2.0)*M[5];
double Mstmp66 = -Mstmp14*Mstmp65 - Mstmp5*Mstmp53 - Mstmp51*M[0];
double Mstmp67 = z*M[14];
double Mstmp68 = Mstmp21*y;
double Mstmp69 = (1.0/2.0)*M[6];
double Mstmp70 = Mstmp14*Mstmp69;
double Mstmp71 = Mstmp53*Mstmp7;
double Mstmp72 = (1.0/2.0)*Mstmp31;
double Mstmp73 = Mstmp1*Mstmp72;
double Mstmp74 = x*M[12];
double Mstmp75 = (1.0/2.0)*M[7];
double Mstmp76 = -Mstmp14*Mstmp75 - Mstmp51*M[1] - Mstmp53*Mstmp8;
double Mstmp77 = x*M[13];
double Mstmp78 = y*M[11];
double Mstmp79 = pow(y, 3);
double Mstmp80 = (1.0/6.0)*Mstmp79;
double Mstmp81 = Mstmp29*y;
double Mstmp82 = x*M[14];
double Mstmp83 = y*M[13];
double Mstmp84 = (1.0/2.0)*M[8];
double Mstmp85 = z*M[17];
double Mstmp86 = Mstmp14*Mstmp84;
double Mstmp87 = pow(x, 4);
double Mstmp88 = (1.0/24.0)*M[0];
double Mstmp89 = Mstmp44*x;
double Mstmp90 = Mstmp10*Mstmp53;
double Mstmp91 = Mstmp13*Mstmp55;
double Mstmp92 = (1.0/4.0)*Mstmp14;
double Mstmp93 = Mstmp11*Mstmp92;
double Mstmp94 = Mstmp93*M[0];
double Mstmp95 = Mstmp1*Mstmp48;
double Mstmp96 = pow(z, 4);
double Mstmp97 = Mstmp5*Mstmp51 + Mstmp51*M[5] + Mstmp88*Mstmp96;
double Mstmp98 = (1.0/2.0)*M[9];
double Mstmp99 = z*M[19];
double Mstmp100 = Mstmp14*Mstmp98;
double Mstmp101 = (1.0/24.0)*M[1];
double Mstmp102 = Mstmp61*x;
double Mstmp103 = Mstmp44*y;
double Mstmp104 = Mstmp17*Mstmp53;
double Mstmp105 = Mstmp18*Mstmp53;
double Mstmp106 = Mstmp21*Mstmp55;
double Mstmp107 = Mstmp93*M[1];
double Mstmp108 = Mstmp45*y;
double Mstmp109 = Mstmp20*Mstmp53;
double Mstmp110 = Mstmp22*Mstmp55;
double Mstmp111 = Mstmp101*Mstmp96 + Mstmp51*Mstmp8 + Mstmp51*M[7];
double Mstmp112 = (1.0/2.0)*M[10];
double Mstmp113 = (1.0/24.0)*M[2];
double Mstmp114 = Mstmp113*Mstmp96;
double Mstmp115 = -Mstmp0*Mstmp51 - Mstmp112*Mstmp14 - Mstmp24*Mstmp53 - Mstmp51*M[3] - Mstmp93*M[2];
double Mstmp116 = z*M[21];
double Mstmp117 = Mstmp67*x;
double Mstmp118 = Mstmp61*y;
double Mstmp119 = Mstmp62*y;
double Mstmp120 = (1.0/2.0)*M[11];
double Mstmp121 = Mstmp120*Mstmp14;
double Mstmp122 = Mstmp29*Mstmp53;
double Mstmp123 = Mstmp30*Mstmp53;
double Mstmp124 = Mstmp13*Mstmp72;
double Mstmp125 = Mstmp32*Mstmp53;
double Mstmp126 = Mstmp15*Mstmp72;
double Mstmp127 = Mstmp31*Mstmp92*M[0];
double Mstmp128 = (1.0/2.0)*M[12];
double Mstmp129 = -Mstmp128*Mstmp14 - Mstmp3*Mstmp51 - Mstmp33*Mstmp53 - Mstmp34*Mstmp53 - Mstmp36*Mstmp53 - Mstmp4*Mstmp51 - Mstmp51*M[4];
double Mstmp130 = z*M[23];
double Mstmp131 = Mstmp67*y;
double Mstmp132 = (1.0/2.0)*M[13];
double Mstmp133 = Mstmp132*Mstmp14;
double Mstmp134 = Mstmp39*Mstmp53;
double Mstmp135 = Mstmp21*Mstmp72;
double Mstmp136 = Mstmp1*Mstmp80;
double Mstmp137 = Mstmp31*Mstmp92*M[1];
double Mstmp138 = (1.0/2.0)*M[14];
double Mstmp139 = Mstmp31*M[2];
double Mstmp140 = -Mstmp138*Mstmp14 - Mstmp139*Mstmp92 - Mstmp40*Mstmp53 - Mstmp51*Mstmp7 - Mstmp51*M[6];
double Mstmp141 = pow(y, 4);
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += Mstmp0 + Mstmp2 + M[3];
#pragma omp atomic
Ms[4] += Mstmp3 + Mstmp4 + M[4];
#pragma omp atomic
Ms[5] += Mstmp5 + Mstmp6 + M[5];
#pragma omp atomic
Ms[6] += Mstmp2 + Mstmp7 + M[6];
#pragma omp atomic
Ms[7] += Mstmp8 + Mstmp9 + M[7];
#pragma omp atomic
Ms[8] += Mstmp10 + Mstmp11*Mstmp12 + Mstmp16 + M[8];
#pragma omp atomic
Ms[9] += Mstmp11*Mstmp19 + Mstmp17 + Mstmp18 + Mstmp20 + Mstmp23 + M[9];
#pragma omp atomic
Ms[10] += Mstmp11*Mstmp27 + Mstmp24 + Mstmp25 + Mstmp26 + Mstmp28 + M[10];
#pragma omp atomic
Ms[11] += Mstmp12*Mstmp31 + Mstmp16 + Mstmp29 + Mstmp30 + Mstmp32 + M[11];
#pragma omp atomic
Ms[12] += Mstmp33 + Mstmp34 + Mstmp35 + Mstmp36 + Mstmp37 + Mstmp38 + M[12];
#pragma omp atomic
Ms[13] += Mstmp19*Mstmp31 + Mstmp23 + Mstmp39 + M[13];
#pragma omp atomic
Ms[14] += Mstmp27*Mstmp31 + Mstmp28 + Mstmp40 + Mstmp41 + Mstmp42 + M[14];
#pragma omp atomic
Ms[15] += Mstmp11*Mstmp46 + Mstmp43 - Mstmp44 - Mstmp45 + Mstmp48*M[0] - Mstmp49 + Mstmp52 - Mstmp54 - Mstmp56 + M[15];
#pragma omp atomic
Ms[16] += Mstmp11*Mstmp59 + Mstmp4*Mstmp55 + Mstmp48*M[1] + Mstmp57 + Mstmp58 + Mstmp60 + Mstmp63 + M[16];
#pragma omp atomic
Ms[17] += Mstmp10*z + Mstmp11*Mstmp65 + Mstmp48*M[2] + Mstmp55*Mstmp6 + Mstmp64 + Mstmp66 + z*M[8] + M[17];
#pragma omp atomic
Ms[18] += (1.0/2.0)*Mstmp11*y*M[1] + (1.0/2.0)*Mstmp11*M[6] + (1.0/2.0)*Mstmp31*x*M[0] + (1.0/2.0)*Mstmp31*M[3] - Mstmp44 - Mstmp45 - Mstmp49 + (1.0/3.0)*Mstmp50*M[2] - Mstmp54 - Mstmp56 - Mstmp67 - Mstmp68 - Mstmp70 - Mstmp71 - Mstmp73 + x*y*M[4] + x*M[11] + y*M[9] + M[18];
#pragma omp atomic
Ms[19] += Mstmp11*Mstmp75 + Mstmp17*z + Mstmp18*z + Mstmp20*z + Mstmp24*y + Mstmp55*Mstmp8 + Mstmp55*Mstmp9 + Mstmp74 + Mstmp76 + y*M[10] + z*M[9] + M[19];
#pragma omp atomic
Ms[20] += Mstmp3*Mstmp72 + Mstmp31*Mstmp59 + Mstmp63 + Mstmp77 + Mstmp78 + Mstmp80*M[0] + Mstmp81 + M[20];
#pragma omp atomic
Ms[21] += Mstmp29*z + Mstmp30*z + Mstmp31*Mstmp65 + Mstmp32*z + Mstmp33*y + Mstmp5*Mstmp72 + Mstmp6*Mstmp72 + Mstmp66 + Mstmp82 + y*M[12] + z*M[11] + M[21];
#pragma omp atomic
Ms[22] += Mstmp31*Mstmp69 + Mstmp52 - Mstmp67 - Mstmp68 - Mstmp70 - Mstmp71 - Mstmp73 + Mstmp80*M[1] + Mstmp83 + M[22];
#pragma omp atomic
Ms[23] += Mstmp31*Mstmp75 + Mstmp39*z + Mstmp72*Mstmp9 + Mstmp76 + Mstmp80*M[2] + y*M[14] + z*M[13] + M[23];
#pragma omp atomic
Ms[24] += Mstmp11*Mstmp84 + Mstmp48*M[3] - Mstmp85 - Mstmp86 + Mstmp87*Mstmp88 - Mstmp89 - Mstmp90 - Mstmp91 - Mstmp94 - Mstmp95 + Mstmp97 + x*M[15] + M[24];
#pragma omp atomic
Ms[25] += -Mstmp100 + Mstmp101*Mstmp87 - Mstmp102 - Mstmp103 - Mstmp104 - Mstmp105 - Mstmp106 - Mstmp107 - Mstmp108 - Mstmp109 + Mstmp11*Mstmp98 - Mstmp110 + Mstmp111 + Mstmp18*Mstmp55 + Mstmp4*Mstmp48 + Mstmp43*y + Mstmp48*M[4] - Mstmp99 + x*M[16] + y*M[15] + M[25];
#pragma omp atomic
Ms[26] += Mstmp11*Mstmp112 + Mstmp113*Mstmp87 + Mstmp114 + Mstmp115 + Mstmp25*Mstmp55 + Mstmp43*z + Mstmp48*Mstmp6 + Mstmp48*M[5] + x*M[17] + z*M[15] + M[26];
#pragma omp atomic
Ms[27] += (1.0/4.0)*Mstmp11*Mstmp31*M[0] + (1.0/2.0)*Mstmp11*y*M[4] + (1.0/2.0)*Mstmp11*M[11] - Mstmp116 - Mstmp117 - Mstmp118 - Mstmp119 - Mstmp121 - Mstmp122 - Mstmp123 - Mstmp124 - Mstmp125 - Mstmp126 - Mstmp127 + (1.0/2.0)*Mstmp31*x*M[3] + (1.0/2.0)*Mstmp31*M[8] + (1.0/6.0)*Mstmp47*y*M[1] + (1.0/6.0)*Mstmp47*M[6] + (1.0/3.0)*Mstmp50*x*M[2] + (1.0/3.0)*Mstmp50*M[5] - Mstmp85 - Mstmp86 - Mstmp89 - Mstmp90 - Mstmp91 - Mstmp94 - Mstmp95 + (1.0/12.0)*Mstmp96*M[0] + x*y*M[9] + x*M[18] + y*M[16] + M[27];
#pragma omp atomic
Ms[28] += Mstmp11*Mstmp128 + Mstmp129 + Mstmp34*Mstmp55 + Mstmp35*Mstmp55 + Mstmp38*Mstmp55 + Mstmp48*Mstmp8 + Mstmp48*Mstmp9 + Mstmp48*M[7] + Mstmp57*z + Mstmp58*z + Mstmp60*z + Mstmp64*y + x*M[19] + y*M[17] + z*M[16] + M[28];
#pragma omp atomic
Ms[29] += -Mstmp100 - Mstmp102 - Mstmp103 - Mstmp104 - Mstmp105 - Mstmp106 - Mstmp107 - Mstmp108 - Mstmp109 + (1.0/4.0)*Mstmp11*Mstmp31*M[1] + (1.0/2.0)*Mstmp11*y*M[6] + (1.0/2.0)*Mstmp11*M[13] - Mstmp110 - Mstmp130 - Mstmp131 - Mstmp133 - Mstmp134 - Mstmp135 - Mstmp136 - Mstmp137 + (1.0/2.0)*Mstmp31*x*M[4] + (1.0/2.0)*Mstmp31*M[9] + (1.0/3.0)*Mstmp50*y*M[2] + (1.0/3.0)*Mstmp50*M[7] + (1.0/6.0)*Mstmp79*x*M[0] + (1.0/6.0)*Mstmp79*M[3] + (1.0/12.0)*Mstmp96*M[1] - Mstmp99 + x*y*M[11] + x*M[20] + y*M[18] + M[29];
#pragma omp atomic
Ms[30] += Mstmp11*Mstmp138 + (1.0/4.0)*Mstmp11*Mstmp139 + Mstmp112*Mstmp31 + Mstmp115 + Mstmp140 + Mstmp17*y*z + Mstmp24*Mstmp72 + Mstmp25*Mstmp72 + Mstmp26*Mstmp72 + Mstmp40*Mstmp55 + Mstmp41*Mstmp55 + Mstmp42*Mstmp55 + Mstmp74*y + (1.0/12.0)*Mstmp96*M[2] + x*z*M[11] + x*M[21] + y*z*M[9] + y*M[19] + z*M[18] + M[30];
#pragma omp atomic
Ms[31] += -Mstmp116 - Mstmp117 - Mstmp118 - Mstmp119 + Mstmp120*Mstmp31 - Mstmp121 - Mstmp122 - Mstmp123 - Mstmp124 - Mstmp125 - Mstmp126 - Mstmp127 + Mstmp141*Mstmp88 + Mstmp29*Mstmp72 + Mstmp3*Mstmp80 + Mstmp77*y + Mstmp80*M[4] + Mstmp97 + x*M[22] + y*M[20] + M[31];
#pragma omp atomic
Ms[32] += Mstmp128*Mstmp31 + Mstmp129 + Mstmp33*Mstmp72 + Mstmp35*Mstmp72 + Mstmp37*Mstmp72 + Mstmp5*Mstmp80 + Mstmp6*Mstmp80 + Mstmp77*z + Mstmp78*z + Mstmp80*M[5] + Mstmp81*z + Mstmp82*y + x*M[23] + y*M[21] + z*M[20] + M[32];
#pragma omp atomic
Ms[33] += Mstmp101*Mstmp141 + Mstmp111 - Mstmp130 - Mstmp131 + Mstmp132*Mstmp31 - Mstmp133 - Mstmp134 - Mstmp135 - Mstmp136 - Mstmp137 + Mstmp80*M[6] + y*M[22] + M[33];
#pragma omp atomic
Ms[34] += Mstmp113*Mstmp141 + Mstmp114 + Mstmp138*Mstmp31 + Mstmp140 + Mstmp41*Mstmp72 + Mstmp80*Mstmp9 + Mstmp80*M[7] + Mstmp83*z + y*M[23] + z*M[22] + M[34];

}

void L2Lc_5(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[13];
double Lstmp4 = Lstmp3*y;
double Lstmp5 = pow(x, 2);
double Lstmp6 = (1.0/2.0)*Lstmp5;
double Lstmp7 = (1.0/6.0)*pow(x, 3);
double Lstmp8 = pow(y, 2);
double Lstmp9 = (1.0/2.0)*Lstmp8;
double Lstmp10 = (1.0/6.0)*pow(y, 3);
double Lstmp11 = x*L[12];
double Lstmp12 = x*L[21];
double Lstmp13 = y*L[10];
double Lstmp14 = z*L[11];
double Lstmp15 = y*L[17];
double Lstmp16 = z*L[18];
double Lstmp17 = z*L[15];
double Lstmp18 = z*L[24];
double Lstmp19 = z*L[22];
double Lstmp20 = Lstmp19*x;
double Lstmp21 = z*L[20];
double Lstmp22 = Lstmp21*y;
double Lstmp23 = L[4] + L[7];
double Lstmp24 = pow(z, 2);
double Lstmp25 = (1.0/2.0)*Lstmp24;
double Lstmp26 = L[11] + L[15];
double Lstmp27 = (1.0/6.0)*pow(z, 3);
double Lstmp28 = L[9] + L[12];
double Lstmp29 = Lstmp25*Lstmp28;
double Lstmp30 = L[18] + L[22];
double Lstmp31 = Lstmp27*Lstmp30;
double Lstmp32 = L[10] + L[14];
double Lstmp33 = Lstmp25*Lstmp32;
double Lstmp34 = L[20] + L[24];
double Lstmp35 = Lstmp27*Lstmp34;
double Lstmp36 = L[17] + L[21];
double Lstmp37 = Lstmp25*Lstmp36;
double Lstmp38 = Lstmp37*y;
double Lstmp39 = L[16] + L[19];
double Lstmp40 = (1.0/4.0)*Lstmp24;
double Lstmp41 = L[19] + L[23];
double Lstmp42 = L[16] + 2*L[19] + L[23];
double Lstmp43 = x*L[19];
double Lstmp44 = Lstmp25*Lstmp39;
double Lstmp45 = y*L[12];
double Lstmp46 = Lstmp19*y;
double Lstmp47 = y*L[19];
double Lstmp48 = Lstmp25*Lstmp41;
double Lstmp49 = y*L[13];
double Lstmp50 = x*L[22];
double Lstmp51 = y*L[20];
double Lstmp52 = Lstmp28*z;
double Lstmp53 = Lstmp32*z;
double Lstmp54 = Lstmp36*z;
double Lstmp55 = Lstmp54*y;
double Lstmp56 = Lstmp25*Lstmp30;
double Lstmp57 = Lstmp39*z;
double Lstmp58 = Lstmp25*Lstmp34;
double Lstmp59 = Lstmp41*z;
double Lstmp60 = y*L[21];
double Lstmp61 = y*L[22];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp10*Lstmp12 + Lstmp10*Lstmp18 + Lstmp10*L[14] + Lstmp11*Lstmp9 + Lstmp13*Lstmp6 + Lstmp14*Lstmp6 + Lstmp15*Lstmp7 + Lstmp16*Lstmp7 + Lstmp17*Lstmp9 + Lstmp2*y + Lstmp20*Lstmp9 + Lstmp22*Lstmp6 - Lstmp23*Lstmp25 - Lstmp26*Lstmp27 - Lstmp29*x - Lstmp31*x - Lstmp33*y - Lstmp35*y - Lstmp38*x - Lstmp39*Lstmp40*Lstmp5 + Lstmp4*x - Lstmp40*Lstmp41*Lstmp8 + (1.0/24.0)*Lstmp42*pow(z, 4) + (1.0/4.0)*Lstmp5*Lstmp8*L[19] + Lstmp6*L[4] + Lstmp7*L[9] + Lstmp9*L[7] + (1.0/24.0)*pow(x, 4)*L[16] + x*L[1] + (1.0/24.0)*pow(y, 4)*L[23] + y*L[2] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + Lstmp10*L[21] + Lstmp13*x + Lstmp14*x + Lstmp15*Lstmp6 + Lstmp16*Lstmp6 + Lstmp19*Lstmp9 + Lstmp22*x - Lstmp29 - Lstmp31 - Lstmp38 + Lstmp4 + Lstmp43*Lstmp9 - Lstmp44*x + Lstmp6*L[9] + Lstmp7*L[16] + Lstmp9*L[12] + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp10*L[23] + Lstmp12*Lstmp9 + Lstmp17*y + Lstmp18*Lstmp9 + Lstmp2 + Lstmp21*Lstmp6 + Lstmp3*x - Lstmp33 - Lstmp35 - Lstmp37*x + Lstmp45*x + Lstmp46*x + Lstmp47*Lstmp6 - Lstmp48*y + Lstmp6*L[10] + Lstmp7*L[17] + Lstmp9*L[14] + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += Lstmp10*L[24] - Lstmp23*z - Lstmp25*Lstmp26 + Lstmp27*Lstmp42 + Lstmp49*x + Lstmp50*Lstmp9 + Lstmp51*Lstmp6 - Lstmp52*x - Lstmp53*y - Lstmp55*x - Lstmp56*x - Lstmp57*Lstmp6 - Lstmp58*y - Lstmp59*Lstmp9 + Lstmp6*L[11] + Lstmp7*L[18] + Lstmp9*L[15] + x*L[6] + y*L[8] + L[3];
#pragma omp atomic
Ls[4] += Lstmp13 + Lstmp14 + Lstmp15*x + Lstmp16*x + Lstmp22 - Lstmp44 + Lstmp6*L[16] + Lstmp9*L[19] + x*L[9] + L[4];
#pragma omp atomic
Ls[5] += Lstmp21*x + Lstmp3 - Lstmp37 + Lstmp45 + Lstmp46 + Lstmp47*x + Lstmp6*L[17] + Lstmp9*L[21] + x*L[10] + L[5];
#pragma omp atomic
Ls[6] += Lstmp49 + Lstmp51*x - Lstmp52 - Lstmp55 - Lstmp56 - Lstmp57*x + Lstmp6*L[18] + Lstmp9*L[22] + x*L[11] + L[6];
#pragma omp atomic
Ls[7] += Lstmp11 + Lstmp17 + Lstmp18*y + Lstmp20 - Lstmp48 + Lstmp6*L[19] + Lstmp60*x + Lstmp9*L[23] + y*L[14] + L[7];
#pragma omp atomic
Ls[8] += -Lstmp53 - Lstmp54*x - Lstmp58 - Lstmp59*y + Lstmp6*L[20] + Lstmp61*x + Lstmp9*L[24] + x*L[13] + y*L[15] + L[8];
#pragma omp atomic
Ls[9] += Lstmp15 + Lstmp16 + x*L[16] + L[9];
#pragma omp atomic
Ls[10] += Lstmp21 + Lstmp47 + x*L[17] + L[10];
#pragma omp atomic
Ls[11] += Lstmp51 - Lstmp57 + x*L[18] + L[11];
#pragma omp atomic
Ls[12] += Lstmp19 + Lstmp43 + Lstmp60 + L[12];
#pragma omp atomic
Ls[13] += -Lstmp54 + Lstmp61 + x*L[20] + L[13];
#pragma omp atomic
Ls[14] += Lstmp12 + Lstmp18 + y*L[23] + L[14];
#pragma omp atomic
Ls[15] += Lstmp50 - Lstmp59 + y*L[24] + L[15];
#pragma omp atomic
Ls[16] += L[16];
#pragma omp atomic
Ls[17] += L[17];
#pragma omp atomic
Ls[18] += L[18];
#pragma omp atomic
Ls[19] += L[19];
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

}

void L2Pc_5(double x, double y, double z, double * L, double * F) {
double Ftmp0 = x*y;
double Ftmp1 = x*z;
double Ftmp2 = y*z;
double Ftmp3 = Ftmp0*z;
double Ftmp4 = pow(x, 2);
double Ftmp5 = (1.0/2.0)*Ftmp4;
double Ftmp6 = (1.0/6.0)*pow(x, 3);
double Ftmp7 = pow(y, 2);
double Ftmp8 = (1.0/2.0)*Ftmp7;
double Ftmp9 = (1.0/6.0)*pow(y, 3);
double Ftmp10 = Ftmp8*x;
double Ftmp11 = Ftmp5*y;
double Ftmp12 = Ftmp5*z;
double Ftmp13 = Ftmp8*z;
double Ftmp14 = pow(z, 2);
double Ftmp15 = L[9] + L[12];
double Ftmp16 = pow(z, 3);
double Ftmp17 = L[18] + L[22];
double Ftmp18 = L[16] + L[19];
double Ftmp19 = L[17] + L[21];
double Ftmp20 = L[10] + L[14];
double Ftmp21 = L[20] + L[24];
double Ftmp22 = L[19] + L[23];
#pragma omp atomic
F[0] += -Ftmp0*L[10] - Ftmp1*L[11] - Ftmp10*L[19] - Ftmp11*L[17] - Ftmp12*L[18] - Ftmp13*L[22] + (1.0/2.0)*Ftmp14*Ftmp15 + (1.0/2.0)*Ftmp14*Ftmp18*x + (1.0/2.0)*Ftmp14*Ftmp19*y + (1.0/6.0)*Ftmp16*Ftmp17 - Ftmp2*L[13] - Ftmp3*L[20] - Ftmp5*L[9] - Ftmp6*L[16] - Ftmp8*L[12] - Ftmp9*L[21] - x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -Ftmp0*L[12] - Ftmp1*L[13] - Ftmp10*L[21] - Ftmp11*L[19] - Ftmp12*L[20] - Ftmp13*L[24] + (1.0/2.0)*Ftmp14*Ftmp19*x + (1.0/2.0)*Ftmp14*Ftmp20 + (1.0/2.0)*Ftmp14*Ftmp22*y + (1.0/6.0)*Ftmp16*Ftmp21 - Ftmp2*L[15] - Ftmp3*L[22] - Ftmp5*L[10] - Ftmp6*L[17] - Ftmp8*L[14] - Ftmp9*L[23] - x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -Ftmp0*L[13] - Ftmp10*L[22] - Ftmp11*L[20] + (1.0/2.0)*Ftmp14*Ftmp17*x + (1.0/2.0)*Ftmp14*Ftmp21*y + (1.0/2.0)*Ftmp14*(L[11] + L[15]) + Ftmp15*x*z - 1.0/6.0*Ftmp16*(L[16] + 2*L[19] + L[23]) + (1.0/2.0)*Ftmp18*Ftmp4*z + Ftmp19*x*y*z + Ftmp20*y*z + (1.0/2.0)*Ftmp22*Ftmp7*z - Ftmp5*L[11] - Ftmp6*L[18] - Ftmp8*L[15] - Ftmp9*L[24] - x*L[6] - y*L[8] + z*(L[4] + L[7]) - L[3];

}

void M2Pc_5(double x, double y, double z, double * M, double * F) {
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(Rinv, 3);
double Ftmp1 = pow(Rinv, 2);
double Ftmp2 = 3*Ftmp1;
double Ftmp3 = Ftmp2*z;
double Ftmp4 = Ftmp2*x;
double Ftmp5 = Ftmp4*y;
double Ftmp6 = Ftmp3*M[2];
double Ftmp7 = pow(Rinv, 4);
double Ftmp8 = pow(x, 2);
double Ftmp9 = Ftmp1*Ftmp8;
double Ftmp10 = 3*Ftmp9;
double Ftmp11 = pow(y, 2);
double Ftmp12 = pow(Rinv, 6);
double Ftmp13 = 30*x;
double Ftmp14 = Ftmp12*Ftmp13;
double Ftmp15 = pow(y, 3);
double Ftmp16 = Ftmp15*M[13];
double Ftmp17 = Ftmp11*z;
double Ftmp18 = Ftmp14*Ftmp17;
double Ftmp19 = Ftmp12*Ftmp8;
double Ftmp20 = 105*Ftmp19;
double Ftmp21 = z*M[12];
double Ftmp22 = pow(Rinv, 8);
double Ftmp23 = 30*Ftmp11;
double Ftmp24 = Ftmp19*Ftmp23;
double Ftmp25 = z*M[32];
double Ftmp26 = 1890*pow(Rinv, 10);
double Ftmp27 = 5*Ftmp9;
double Ftmp28 = Ftmp27 - 3;
double Ftmp29 = Ftmp1*Ftmp11;
double Ftmp30 = 5*Ftmp29;
double Ftmp31 = Ftmp30 - 1;
double Ftmp32 = Ftmp9 - 1;
double Ftmp33 = Ftmp10 - 1;
double Ftmp34 = 3*Ftmp29;
double Ftmp35 = Ftmp34 - 1;
double Ftmp36 = 15*Ftmp7;
double Ftmp37 = Ftmp36*y;
double Ftmp38 = 7*Ftmp9;
double Ftmp39 = Ftmp38 - 3;
double Ftmp40 = Ftmp39*M[16];
double Ftmp41 = 7*Ftmp29;
double Ftmp42 = Ftmp41 - 3;
double Ftmp43 = Ftmp42*M[20];
double Ftmp44 = Ftmp36*z;
double Ftmp45 = Ftmp39*M[17];
double Ftmp46 = Ftmp41 - 1;
double Ftmp47 = Ftmp46*M[21];
double Ftmp48 = Ftmp32*y;
double Ftmp49 = Ftmp13*Ftmp7;
double Ftmp50 = Ftmp7*Ftmp8;
double Ftmp51 = 30*M[8];
double Ftmp52 = Ftmp27 - 1;
double Ftmp53 = Ftmp52*M[9];
double Ftmp54 = Ftmp37*x;
double Ftmp55 = Ftmp30 - 3;
double Ftmp56 = Ftmp55*M[13];
double Ftmp57 = Ftmp44*x;
double Ftmp58 = Ftmp52*M[10];
double Ftmp59 = Ftmp31*M[14];
double Ftmp60 = Ftmp28*M[8];
double Ftmp61 = Ftmp36*Ftmp8;
double Ftmp62 = Ftmp31*M[11];
double Ftmp63 = Ftmp38 - 1;
double Ftmp64 = Ftmp15*x;
double Ftmp65 = 420*M[33];
double Ftmp66 = Ftmp22*Ftmp65*(9*Ftmp29 - 5);
double Ftmp67 = z*M[28];
double Ftmp68 = 1890*Ftmp22;
double Ftmp69 = Ftmp22*Ftmp35;
double Ftmp70 = 1260*M[34];
double Ftmp71 = Ftmp69*Ftmp70;
double Ftmp72 = Ftmp17*x;
double Ftmp73 = Ftmp22*y;
double Ftmp74 = 2835*Ftmp69;
double Ftmp75 = Ftmp25*Ftmp74;
double Ftmp76 = Ftmp11*Ftmp8;
double Ftmp77 = 1260*M[31];
double Ftmp78 = Ftmp69*Ftmp77;
double Ftmp79 = Ftmp7*pow(x, 4);
double Ftmp80 = 63*Ftmp79 - 70*Ftmp9 + 15;
double Ftmp81 = -14*Ftmp29;
double Ftmp82 = Ftmp7*pow(y, 4);
double Ftmp83 = Ftmp81 + 21*Ftmp82 + 1;
double Ftmp84 = -10*Ftmp9;
double Ftmp85 = Ftmp84 + 3;
double Ftmp86 = 35*Ftmp79 - 30*Ftmp9 + 3;
double Ftmp87 = -30*Ftmp29 + 35*Ftmp82 + 3;
double Ftmp88 = 3*Ftmp79 - 4*Ftmp9 + 1;
double Ftmp89 = Ftmp12*y;
double Ftmp90 = Ftmp89*x;
double Ftmp91 = 1260*M[25];
double Ftmp92 = 14*Ftmp9;
double Ftmp93 = -Ftmp92;
double Ftmp94 = 21*Ftmp79 + Ftmp93 + 1;
double Ftmp95 = 315*Ftmp94;
double Ftmp96 = Ftmp95*M[25];
double Ftmp97 = -70*Ftmp29 + 63*Ftmp82 + 15;
double Ftmp98 = Ftmp97*M[33];
double Ftmp99 = 105*Ftmp90;
double Ftmp100 = Ftmp12*z;
double Ftmp101 = Ftmp100*x;
double Ftmp102 = 1260*M[26];
double Ftmp103 = Ftmp95*M[26];
double Ftmp104 = 315*Ftmp83*M[34];
double Ftmp105 = 420*M[24];
double Ftmp106 = Ftmp80*M[24];
double Ftmp107 = 315*Ftmp83*M[31];
double Ftmp108 = -Ftmp38;
double Ftmp109 = Ftmp11*Ftmp50;
double Ftmp110 = 63*Ftmp109;
double Ftmp111 = Ftmp110 + 3;
double Ftmp112 = Ftmp108 + Ftmp111 - 21*Ftmp29;
double Ftmp113 = 14*Ftmp109;
double Ftmp114 = 1 - Ftmp9;
double Ftmp115 = 35*Ftmp109 - Ftmp27 - Ftmp30 + 1;
double Ftmp116 = 18*Ftmp109;
double Ftmp117 = -10*Ftmp29;
double Ftmp118 = Ftmp117 + 3;
double Ftmp119 = 210*M[29];
double Ftmp120 = -Ftmp41;
double Ftmp121 = Ftmp111 + Ftmp120 - 21*Ftmp9;
double Ftmp122 = Ftmp121*M[29];
double Ftmp123 = Ftmp114 + Ftmp116;
double Ftmp124 = 210*M[30];
double Ftmp125 = Ftmp108 + Ftmp110 + Ftmp120 + 1;
double Ftmp126 = 105*Ftmp125*M[30];
double Ftmp127 = Ftmp112*M[27];
double Ftmp128 = pow(x, 3);
double Ftmp129 = Ftmp128*Ftmp51;
double Ftmp130 = Ftmp11*Ftmp12;
double Ftmp131 = 105*Ftmp130;
double Ftmp132 = 30*y*z;
double Ftmp133 = Ftmp132*Ftmp19;
double Ftmp134 = Ftmp17*M[28];
double Ftmp135 = Ftmp128*Ftmp26;
double Ftmp136 = Ftmp29 - 1;
double Ftmp137 = Ftmp36*x;
double Ftmp138 = Ftmp63*M[19];
double Ftmp139 = Ftmp42*M[23];
double Ftmp140 = Ftmp136*Ftmp7;
double Ftmp141 = Ftmp37*z;
double Ftmp142 = Ftmp11*Ftmp36;
double Ftmp143 = Ftmp105*Ftmp128*(9*Ftmp9 - 5);
double Ftmp144 = Ftmp11*x;
double Ftmp145 = Ftmp22*Ftmp33;
double Ftmp146 = 2835*Ftmp145*x;
double Ftmp147 = y*z;
double Ftmp148 = Ftmp102*Ftmp145*Ftmp8;
double Ftmp149 = Ftmp145*Ftmp91;
double Ftmp150 = -4*Ftmp29 + 3*Ftmp82 + 1;
double Ftmp151 = Ftmp89*z;
double Ftmp152 = 1 - Ftmp29;
double Ftmp153 = Ftmp116 + Ftmp152;
double Ftmp154 = pow(z, 2);
double Ftmp155 = 30*Ftmp100;
double Ftmp156 = Ftmp154*M[10];
double Ftmp157 = Ftmp12*Ftmp154;
double Ftmp158 = Ftmp154*M[32];
double Ftmp159 = Ftmp154*y*M[28];
double Ftmp160 = Ftmp11 + Ftmp8;
double Ftmp161 = 105*Ftmp101;
double Ftmp162 = 105*Ftmp151;
double Ftmp163 = -18*Ftmp11*Ftmp9;
#pragma omp atomic
F[0] += Ftmp0*(3*Ftmp1*Ftmp28*M[8] + 3*Ftmp1*Ftmp31*M[11] + 6*Ftmp1*Ftmp32*x*M[3] + 3*Ftmp1*Ftmp33*x*M[3] + 3*Ftmp1*Ftmp35*x*M[6] - Ftmp10*M[0] - Ftmp101*Ftmp102*Ftmp88 - Ftmp101*Ftmp103 - Ftmp101*Ftmp104 - Ftmp101*Ftmp124*(Ftmp117 + Ftmp123) - Ftmp101*Ftmp126 - Ftmp105*Ftmp19*(9*Ftmp79 + Ftmp93 + 5) - Ftmp106*Ftmp20 - Ftmp107*Ftmp19 + 60*Ftmp11*Ftmp12*Ftmp42*x*M[22] + 210*Ftmp11*Ftmp22*Ftmp8*z*M[21] + 6*Ftmp11*Ftmp7*x*M[6] + 15*Ftmp112*Ftmp7*M[27] + 15*Ftmp115*Ftmp7*x*M[18] - Ftmp119*Ftmp90*(-Ftmp10 + Ftmp116 + Ftmp118) + 210*Ftmp12*Ftmp32*Ftmp8*y*M[16] + 210*Ftmp12*Ftmp32*Ftmp8*z*M[17] + 210*Ftmp12*Ftmp32*x*y*z*M[19] + 315*Ftmp12*Ftmp33*y*z*M[28] + 315*Ftmp12*Ftmp35*y*z*M[32] + 105*Ftmp12*Ftmp39*Ftmp8*y*M[16] + 105*Ftmp12*Ftmp39*Ftmp8*z*M[17] + 105*Ftmp12*Ftmp42*Ftmp8*y*M[20] + 105*Ftmp12*Ftmp42*x*y*z*M[23] + 105*Ftmp12*Ftmp46*Ftmp8*z*M[21] + 105*Ftmp12*Ftmp63*x*y*z*M[19] - Ftmp122*Ftmp99 - Ftmp127*Ftmp20 - Ftmp14*Ftmp16 + 210*Ftmp15*Ftmp22*Ftmp8*M[20] + 210*Ftmp15*Ftmp22*x*z*M[23] - Ftmp15*Ftmp25*Ftmp26*Ftmp8 - Ftmp18*M[14] - 210*Ftmp19*(Ftmp123 - 12*Ftmp29)*M[27] - Ftmp2*y*M[4] - Ftmp20*Ftmp21*y - Ftmp24*M[11] - Ftmp3*M[5] - Ftmp32*Ftmp49*z*M[10] - Ftmp32*Ftmp50*Ftmp51 - 2835*Ftmp33*Ftmp67*Ftmp73*Ftmp8 - Ftmp37*Ftmp40 - Ftmp37*Ftmp43 - Ftmp44*Ftmp45 - Ftmp44*Ftmp47 - Ftmp48*Ftmp49*M[9] - Ftmp48*Ftmp67*Ftmp68*Ftmp8 - Ftmp5*M[1] - Ftmp53*Ftmp54 - Ftmp54*Ftmp56 - Ftmp57*Ftmp58 - Ftmp57*Ftmp59 - Ftmp6*x - Ftmp60*Ftmp61 - Ftmp61*Ftmp62 - Ftmp64*Ftmp66 + 15*Ftmp7*Ftmp8*y*M[4] + 15*Ftmp7*Ftmp8*z*M[5] + 15*Ftmp7*Ftmp80*M[24] + 45*Ftmp7*Ftmp83*M[31] + 15*Ftmp7*Ftmp86*x*M[15] + 15*Ftmp7*Ftmp87*x*M[22] + 15*Ftmp7*x*y*z*M[7] + 60*Ftmp7*x*(7*Ftmp79 + Ftmp85)*M[15] + 30*Ftmp7*x*(Ftmp113 + Ftmp114 - 8*Ftmp29)*M[18] + 15*Ftmp7*y*z*M[12] - Ftmp71*Ftmp72 - Ftmp75*Ftmp8*y - Ftmp76*Ftmp78 - Ftmp88*Ftmp90*Ftmp91 - Ftmp90*Ftmp96 - Ftmp98*Ftmp99 + M[0]);
#pragma omp atomic
F[1] += Ftmp0*(6*Ftmp1*Ftmp136*y*M[6] + 3*Ftmp1*Ftmp33*y*M[3] + 3*Ftmp1*Ftmp35*y*M[6] + 3*Ftmp1*Ftmp52*M[9] + 3*Ftmp1*Ftmp55*M[13] - Ftmp103*Ftmp151 - Ftmp104*Ftmp151 - Ftmp106*Ftmp99 - Ftmp107*Ftmp90 + 210*Ftmp11*Ftmp12*Ftmp136*x*M[20] + 210*Ftmp11*Ftmp12*Ftmp136*z*M[23] + 105*Ftmp11*Ftmp12*Ftmp39*x*M[16] + 105*Ftmp11*Ftmp12*Ftmp42*x*M[20] + 105*Ftmp11*Ftmp12*Ftmp42*z*M[23] + 105*Ftmp11*Ftmp12*Ftmp63*z*M[19] + 210*Ftmp11*Ftmp128*Ftmp22*M[16] + 210*Ftmp11*Ftmp22*Ftmp8*z*M[19] + 15*Ftmp11*Ftmp7*x*M[4] + 15*Ftmp11*Ftmp7*z*M[7] + 15*Ftmp115*Ftmp7*y*M[18] - Ftmp119*Ftmp130*(Ftmp153 - 12*Ftmp9) + 210*Ftmp12*Ftmp136*x*y*z*M[21] + 315*Ftmp12*Ftmp33*x*z*M[28] + 315*Ftmp12*Ftmp35*x*z*M[32] + 60*Ftmp12*Ftmp39*Ftmp8*y*M[15] + 105*Ftmp12*Ftmp39*x*y*z*M[17] + 105*Ftmp12*Ftmp46*x*y*z*M[21] + 15*Ftmp121*Ftmp7*M[29] - Ftmp122*Ftmp131 - Ftmp124*Ftmp151*(Ftmp153 + Ftmp84) - Ftmp126*Ftmp151 - Ftmp127*Ftmp99 + 210*Ftmp128*Ftmp22*y*z*M[17] - Ftmp129*Ftmp89 - Ftmp13*Ftmp140*y*M[11] - Ftmp130*Ftmp65*(Ftmp81 + 9*Ftmp82 + 5) - Ftmp130*Ftmp96 - Ftmp131*Ftmp21*x - Ftmp131*Ftmp98 - Ftmp132*Ftmp140*M[14] - Ftmp133*M[10] - Ftmp134*Ftmp135 - Ftmp134*Ftmp146 - Ftmp136*Ftmp144*Ftmp25*Ftmp68 - Ftmp137*Ftmp40 - Ftmp137*Ftmp43 - Ftmp138*Ftmp44 - Ftmp139*Ftmp44 - Ftmp140*Ftmp23*M[13] - Ftmp141*Ftmp58 - Ftmp141*Ftmp59 - Ftmp142*Ftmp53 - Ftmp142*Ftmp56 - Ftmp143*Ftmp73 - Ftmp144*Ftmp75 - Ftmp147*Ftmp148 - Ftmp149*Ftmp76 - Ftmp150*Ftmp151*Ftmp70 - Ftmp150*Ftmp77*Ftmp90 - Ftmp24*M[9] - Ftmp3*M[7] - Ftmp34*M[1] - Ftmp4*M[4] - Ftmp5*M[0] - Ftmp54*Ftmp60 - Ftmp54*Ftmp62 - Ftmp6*y + 6*Ftmp7*Ftmp8*y*M[3] + 15*Ftmp7*Ftmp86*y*M[15] + 15*Ftmp7*Ftmp87*y*M[22] + 45*Ftmp7*Ftmp94*M[25] + 15*Ftmp7*Ftmp97*M[33] + 15*Ftmp7*x*y*z*M[5] + 15*Ftmp7*x*z*M[12] + 60*Ftmp7*y*(Ftmp118 + 7*Ftmp82)*M[22] + 30*Ftmp7*y*(Ftmp113 + Ftmp152 - 8*Ftmp9)*M[18] - 210*Ftmp90*(Ftmp116 - Ftmp34 + Ftmp85)*M[27] + M[1]);
#pragma omp atomic
F[2] += Ftmp0*(3*Ftmp1*Ftmp31*M[14] + 3*Ftmp1*Ftmp33*z*M[3] + 3*Ftmp1*Ftmp35*z*M[6] + 3*Ftmp1*Ftmp52*M[10] - Ftmp100*Ftmp129 - Ftmp101*Ftmp107 - Ftmp103*Ftmp157 - Ftmp104*Ftmp157 - Ftmp106*Ftmp161 + 60*Ftmp11*Ftmp12*Ftmp42*z*M[22] + 210*Ftmp11*Ftmp154*Ftmp22*x*M[21] - Ftmp11*Ftmp154*Ftmp71 + 6*Ftmp11*Ftmp7*z*M[6] + 15*Ftmp115*Ftmp7*z*M[18] + 105*Ftmp12*Ftmp154*Ftmp39*x*M[17] + 105*Ftmp12*Ftmp154*Ftmp42*y*M[23] + 105*Ftmp12*Ftmp154*Ftmp46*x*M[21] + 105*Ftmp12*Ftmp154*Ftmp63*y*M[19] + 315*Ftmp12*Ftmp33*x*y*M[28] + 315*Ftmp12*Ftmp35*x*y*M[32] + 60*Ftmp12*Ftmp39*Ftmp8*z*M[15] + 105*Ftmp12*Ftmp39*x*y*z*M[16] + 105*Ftmp12*Ftmp42*x*y*z*M[20] - Ftmp122*Ftmp162 + 15*Ftmp125*Ftmp7*M[30] - Ftmp126*Ftmp157 - Ftmp127*Ftmp161 + 210*Ftmp128*Ftmp154*Ftmp22*M[17] + 210*Ftmp128*Ftmp22*y*z*M[16] - Ftmp133*M[9] - Ftmp135*Ftmp159 - Ftmp137*Ftmp45 - Ftmp137*Ftmp47 - Ftmp138*Ftmp37 - Ftmp139*Ftmp37 - Ftmp141*Ftmp53 - Ftmp141*Ftmp56 - Ftmp143*Ftmp22*z - Ftmp146*Ftmp159 - Ftmp147*Ftmp149*Ftmp8 - Ftmp148*Ftmp154 + 210*Ftmp15*Ftmp154*Ftmp22*M[23] + 210*Ftmp15*Ftmp22*x*z*M[20] - Ftmp15*Ftmp66*z - Ftmp151*Ftmp96 - Ftmp154*Ftmp2*M[2] + 210*Ftmp154*Ftmp22*Ftmp8*y*M[19] + 210*Ftmp154*Ftmp22*(Ftmp160 + Ftmp163)*M[30] - Ftmp154*Ftmp36*Ftmp59 + 15*Ftmp154*Ftmp7*x*M[5] + 15*Ftmp154*Ftmp7*y*M[7] - Ftmp154*Ftmp99*M[12] - Ftmp155*Ftmp16 - Ftmp155*(-Ftmp11*Ftmp92 + Ftmp160)*M[18] - 30*Ftmp156*Ftmp19 - Ftmp156*Ftmp36*Ftmp52 - Ftmp157*Ftmp23*M[14] - Ftmp158*Ftmp26*Ftmp64 - Ftmp158*Ftmp74*x*y - Ftmp162*Ftmp98 - Ftmp18*M[11] - Ftmp2*y*M[7] + 210*Ftmp22*x*z*(3*Ftmp11 + Ftmp163 + Ftmp8)*M[27] + 210*Ftmp22*y*z*(Ftmp11 + Ftmp163 + 3*Ftmp8)*M[29] - Ftmp3*x*M[0] - Ftmp3*y*M[1] - Ftmp4*M[5] - Ftmp57*Ftmp60 - Ftmp57*Ftmp62 + 6*Ftmp7*Ftmp8*z*M[3] + 45*Ftmp7*Ftmp83*M[34] + 15*Ftmp7*Ftmp86*z*M[15] + 15*Ftmp7*Ftmp87*z*M[22] + 45*Ftmp7*Ftmp94*M[26] + 15*Ftmp7*x*y*z*M[4] + 15*Ftmp7*x*y*M[12] - Ftmp72*Ftmp78 + M[2]);

}

void M2Lc_5(double x, double y, double z, double * M, double * L) {
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double D[45];
double Dtmp0 = pow(Rinv, 3);
double Dtmp1 = pow(x, 2);
double Dtmp2 = pow(Rinv, 2);
double Dtmp3 = 3*Dtmp2;
double Dtmp4 = Dtmp1*Dtmp3 - 1;
double Dtmp5 = 3*pow(Rinv, 5);
double Dtmp6 = Dtmp5*x;
double Dtmp7 = pow(y, 2);
double Dtmp8 = Dtmp3*Dtmp7 - 1;
double Dtmp9 = Dtmp5*y;
double Dtmp10 = 5*Dtmp2;
double Dtmp11 = Dtmp1*Dtmp10;
double Dtmp12 = Dtmp11 - 1;
double Dtmp13 = Dtmp5*z;
double Dtmp14 = Dtmp10*Dtmp7;
double Dtmp15 = Dtmp14 - 1;
double Dtmp16 = pow(Rinv, 7);
double Dtmp17 = 15*Dtmp16;
double Dtmp18 = Dtmp17*x;
double Dtmp19 = Dtmp18*y;
double Dtmp20 = Dtmp1*Dtmp2;
double Dtmp21 = pow(x, 4);
double Dtmp22 = pow(Rinv, 4);
double Dtmp23 = 35*Dtmp22;
double Dtmp24 = 7*Dtmp20;
double Dtmp25 = Dtmp24 - 3;
double Dtmp26 = Dtmp18*z;
double Dtmp27 = Dtmp1*Dtmp7;
double Dtmp28 = Dtmp17*y;
double Dtmp29 = Dtmp28*z;
double Dtmp30 = Dtmp2*Dtmp7;
double Dtmp31 = 7*Dtmp30;
double Dtmp32 = Dtmp31 - 3;
double Dtmp33 = pow(y, 4);
double Dtmp34 = Dtmp21*Dtmp22;
double Dtmp35 = 45*Dtmp16;
double Dtmp36 = Dtmp35*(-14*Dtmp20 + 21*Dtmp34 + 1);
double Dtmp37 = -Dtmp24;
double Dtmp38 = 63*Dtmp22*Dtmp27;
double Dtmp39 = Dtmp38 + 3;
double Dtmp40 = 315*pow(Rinv, 9)*x*y*z;
double Dtmp41 = -Dtmp31;
double Dtmp42 = Dtmp22*Dtmp33;
double Dtmp43 = Dtmp35*(-14*Dtmp30 + 21*Dtmp42 + 1);
D[0] = -Dtmp0*x;
D[1] = -Dtmp0*y;
D[2] = -Dtmp0*z;
D[3] = Dtmp0*Dtmp4;
D[4] = Dtmp6*y;
D[5] = Dtmp6*z;
D[6] = Dtmp0*Dtmp8;
D[7] = Dtmp9*z;
D[8] = -D[3] - D[6];
D[9] = -Dtmp6*(Dtmp11 - 3);
D[10] = -Dtmp12*Dtmp9;
D[11] = -Dtmp12*Dtmp13;
D[12] = -Dtmp15*Dtmp6;
D[13] = -Dtmp19*z;
D[14] = -D[9] - D[12];
D[15] = -Dtmp9*(Dtmp14 - 3);
D[16] = -Dtmp13*Dtmp15;
D[17] = -D[10] - D[15];
D[18] = Dtmp5*(-30*Dtmp20 + Dtmp21*Dtmp23 + 3);
D[19] = Dtmp19*Dtmp25;
D[20] = Dtmp25*Dtmp26;
D[21] = Dtmp5*(-Dtmp11 - Dtmp14 + Dtmp23*Dtmp27 + 1);
D[22] = Dtmp29*(Dtmp24 - 1);
D[23] = -D[18] - D[21];
D[24] = Dtmp19*Dtmp32;
D[25] = Dtmp26*(Dtmp31 - 1);
D[26] = -D[19] - D[24];
D[27] = Dtmp5*(Dtmp23*Dtmp33 - 30*Dtmp30 + 3);
D[28] = Dtmp29*Dtmp32;
D[29] = -D[21] - D[27];
D[30] = -Dtmp18*(-70*Dtmp20 + 63*Dtmp34 + 15);
D[31] = -Dtmp36*y;
D[32] = -Dtmp36*z;
D[33] = -Dtmp18*(-21*Dtmp30 + Dtmp37 + Dtmp39);
D[34] = -Dtmp4*Dtmp40;
D[35] = -D[30] - D[33];
D[36] = -Dtmp28*(-21*Dtmp20 + Dtmp39 + Dtmp41);
D[37] = -Dtmp17*z*(Dtmp37 + Dtmp38 + Dtmp41 + 1);
D[38] = -D[31] - D[36];
D[39] = -Dtmp43*x;
D[40] = -Dtmp40*Dtmp8;
D[41] = -D[33] - D[39];
D[42] = -Dtmp28*(-70*Dtmp30 + 63*Dtmp42 + 15);
D[43] = -Dtmp43*z;
D[44] = -D[36] - D[42];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[9]*M[8] + D[10]*M[9] + D[11]*M[10] + D[12]*M[11] + D[13]*M[12] + D[15]*M[13] + D[16]*M[14] + D[18]*M[15] + D[19]*M[16] + D[20]*M[17] + D[21]*M[18] + D[22]*M[19] + D[24]*M[20] + D[25]*M[21] + D[27]*M[22] + D[28]*M[23] + D[30]*M[24] + D[31]*M[25] + D[32]*M[26] + D[33]*M[27] + D[34]*M[28] + D[36]*M[29] + D[37]*M[30] + D[39]*M[31] + D[40]*M[32] + D[42]*M[33] + D[43]*M[34];
#pragma omp atomic
L[1] += D[3]*M[0] + D[4]*M[1] + D[5]*M[2] + D[9]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[18]*M[8] + D[19]*M[9] + D[20]*M[10] + D[21]*M[11] + D[22]*M[12] + D[24]*M[13] + D[25]*M[14] + D[30]*M[15] + D[31]*M[16] + D[32]*M[17] + D[33]*M[18] + D[34]*M[19] + D[36]*M[20] + D[37]*M[21] + D[39]*M[22] + D[40]*M[23];
#pragma omp atomic
L[2] += D[4]*M[0] + D[6]*M[1] + D[7]*M[2] + D[10]*M[3] + D[12]*M[4] + D[13]*M[5] + D[15]*M[6] + D[16]*M[7] + D[19]*M[8] + D[21]*M[9] + D[22]*M[10] + D[24]*M[11] + D[25]*M[12] + D[27]*M[13] + D[28]*M[14] + D[31]*M[15] + D[33]*M[16] + D[34]*M[17] + D[36]*M[18] + D[37]*M[19] + D[39]*M[20] + D[40]*M[21] + D[42]*M[22] + D[43]*M[23];
#pragma omp atomic
L[3] += D[5]*M[0] + D[7]*M[1] + D[8]*M[2] + D[11]*M[3] + D[13]*M[4] + D[14]*M[5] + D[16]*M[6] + D[17]*M[7] + D[20]*M[8] + D[22]*M[9] + D[23]*M[10] + D[25]*M[11] + D[26]*M[12] + D[28]*M[13] + D[29]*M[14] + D[32]*M[15] + D[34]*M[16] + D[35]*M[17] + D[37]*M[18] + D[38]*M[19] + D[40]*M[20] + D[41]*M[21] + D[43]*M[22] + D[44]*M[23];
#pragma omp atomic
L[4] += D[9]*M[0] + D[10]*M[1] + D[11]*M[2] + D[18]*M[3] + D[19]*M[4] + D[20]*M[5] + D[21]*M[6] + D[22]*M[7] + D[30]*M[8] + D[31]*M[9] + D[32]*M[10] + D[33]*M[11] + D[34]*M[12] + D[36]*M[13] + D[37]*M[14];
#pragma omp atomic
L[5] += D[10]*M[0] + D[12]*M[1] + D[13]*M[2] + D[19]*M[3] + D[21]*M[4] + D[22]*M[5] + D[24]*M[6] + D[25]*M[7] + D[31]*M[8] + D[33]*M[9] + D[34]*M[10] + D[36]*M[11] + D[37]*M[12] + D[39]*M[13] + D[40]*M[14];
#pragma omp atomic
L[6] += D[11]*M[0] + D[13]*M[1] + D[14]*M[2] + D[20]*M[3] + D[22]*M[4] + D[23]*M[5] + D[25]*M[6] + D[26]*M[7] + D[32]*M[8] + D[34]*M[9] + D[35]*M[10] + D[37]*M[11] + D[38]*M[12] + D[40]*M[13] + D[41]*M[14];
#pragma omp atomic
L[7] += D[12]*M[0] + D[15]*M[1] + D[16]*M[2] + D[21]*M[3] + D[24]*M[4] + D[25]*M[5] + D[27]*M[6] + D[28]*M[7] + D[33]*M[8] + D[36]*M[9] + D[37]*M[10] + D[39]*M[11] + D[40]*M[12] + D[42]*M[13] + D[43]*M[14];
#pragma omp atomic
L[8] += D[13]*M[0] + D[16]*M[1] + D[17]*M[2] + D[22]*M[3] + D[25]*M[4] + D[26]*M[5] + D[28]*M[6] + D[29]*M[7] + D[34]*M[8] + D[37]*M[9] + D[38]*M[10] + D[40]*M[11] + D[41]*M[12] + D[43]*M[13] + D[44]*M[14];
#pragma omp atomic
L[9] += D[18]*M[0] + D[19]*M[1] + D[20]*M[2] + D[30]*M[3] + D[31]*M[4] + D[32]*M[5] + D[33]*M[6] + D[34]*M[7];
#pragma omp atomic
L[10] += D[19]*M[0] + D[21]*M[1] + D[22]*M[2] + D[31]*M[3] + D[33]*M[4] + D[34]*M[5] + D[36]*M[6] + D[37]*M[7];
#pragma omp atomic
L[11] += D[20]*M[0] + D[22]*M[1] + D[23]*M[2] + D[32]*M[3] + D[34]*M[4] + D[35]*M[5] + D[37]*M[6] + D[38]*M[7];
#pragma omp atomic
L[12] += D[21]*M[0] + D[24]*M[1] + D[25]*M[2] + D[33]*M[3] + D[36]*M[4] + D[37]*M[5] + D[39]*M[6] + D[40]*M[7];
#pragma omp atomic
L[13] += D[22]*M[0] + D[25]*M[1] + D[26]*M[2] + D[34]*M[3] + D[37]*M[4] + D[38]*M[5] + D[40]*M[6] + D[41]*M[7];
#pragma omp atomic
L[14] += D[24]*M[0] + D[27]*M[1] + D[28]*M[2] + D[36]*M[3] + D[39]*M[4] + D[40]*M[5] + D[42]*M[6] + D[43]*M[7];
#pragma omp atomic
L[15] += D[25]*M[0] + D[28]*M[1] + D[29]*M[2] + D[37]*M[3] + D[40]*M[4] + D[41]*M[5] + D[43]*M[6] + D[44]*M[7];
#pragma omp atomic
L[16] += D[30]*M[0] + D[31]*M[1] + D[32]*M[2];
#pragma omp atomic
L[17] += D[31]*M[0] + D[33]*M[1] + D[34]*M[2];
#pragma omp atomic
L[18] += D[32]*M[0] + D[34]*M[1] + D[35]*M[2];
#pragma omp atomic
L[19] += D[33]*M[0] + D[36]*M[1] + D[37]*M[2];
#pragma omp atomic
L[20] += D[34]*M[0] + D[37]*M[1] + D[38]*M[2];
#pragma omp atomic
L[21] += D[36]*M[0] + D[39]*M[1] + D[40]*M[2];
#pragma omp atomic
L[22] += D[37]*M[0] + D[40]*M[1] + D[41]*M[2];
#pragma omp atomic
L[23] += D[39]*M[0] + D[42]*M[1] + D[43]*M[2];
#pragma omp atomic
L[24] += D[40]*M[0] + D[43]*M[1] + D[44]*M[2];

}

void S2M_6(double x, double y, double z, double * S, double * M) {
double Mtmp0 = x*S[1];
double Mtmp1 = y*S[0];
double Mtmp2 = Mtmp0 + Mtmp1;
double Mtmp3 = x*S[2];
double Mtmp4 = z*S[0];
double Mtmp5 = Mtmp3 + Mtmp4;
double Mtmp6 = y*S[2];
double Mtmp7 = z*S[1];
double Mtmp8 = Mtmp6 + Mtmp7;
double Mtmp9 = pow(x, 2);
double Mtmp10 = (1.0/2.0)*Mtmp0;
double Mtmp11 = (1.0/2.0)*Mtmp3;
double Mtmp12 = (1.0/2.0)*Mtmp1;
double Mtmp13 = Mtmp1*z;
double Mtmp14 = Mtmp3*y;
double Mtmp15 = Mtmp0*z;
double Mtmp16 = Mtmp14 + Mtmp15;
double Mtmp17 = Mtmp13 + Mtmp16;
double Mtmp18 = pow(y, 2);
double Mtmp19 = pow(z, 2);
double Mtmp20 = pow(x, 3);
double Mtmp21 = 3*Mtmp1;
double Mtmp22 = (1.0/6.0)*Mtmp9;
double Mtmp23 = 3*Mtmp4;
double Mtmp24 = (1.0/2.0)*x;
double Mtmp25 = Mtmp11*y;
double Mtmp26 = Mtmp10*z;
double Mtmp27 = 3*Mtmp0;
double Mtmp28 = (1.0/6.0)*Mtmp18;
double Mtmp29 = Mtmp12*z;
double Mtmp30 = 3*Mtmp3;
double Mtmp31 = (1.0/6.0)*Mtmp19;
double Mtmp32 = pow(y, 3);
double Mtmp33 = 3*Mtmp7;
double Mtmp34 = y*z;
double Mtmp35 = 3*Mtmp6;
double Mtmp36 = pow(z, 3);
double Mtmp37 = pow(x, 4);
double Mtmp38 = (1.0/24.0)*Mtmp20;
double Mtmp39 = 2*Mtmp0;
double Mtmp40 = (1.0/12.0)*Mtmp9;
double Mtmp41 = Mtmp40*y;
double Mtmp42 = 3*Mtmp13;
double Mtmp43 = 2*Mtmp3;
double Mtmp44 = Mtmp40*z;
double Mtmp45 = 2*Mtmp1;
double Mtmp46 = (1.0/12.0)*x;
double Mtmp47 = Mtmp18*Mtmp46;
double Mtmp48 = 2*Mtmp13;
double Mtmp49 = 2*Mtmp15;
double Mtmp50 = Mtmp14 + Mtmp49;
double Mtmp51 = (1.0/4.0)*x;
double Mtmp52 = 2*Mtmp14;
double Mtmp53 = Mtmp15 + Mtmp52;
double Mtmp54 = 2*Mtmp4;
double Mtmp55 = Mtmp19*Mtmp46;
double Mtmp56 = (1.0/24.0)*Mtmp32;
double Mtmp57 = 3*Mtmp15;
double Mtmp58 = Mtmp14 + Mtmp57;
double Mtmp59 = Mtmp13 + Mtmp52;
double Mtmp60 = 3*Mtmp14;
double Mtmp61 = Mtmp15 + Mtmp60;
double Mtmp62 = (1.0/24.0)*Mtmp36;
double Mtmp63 = pow(y, 4);
double Mtmp64 = 2*Mtmp6;
double Mtmp65 = (1.0/12.0)*Mtmp18;
double Mtmp66 = Mtmp65*z;
double Mtmp67 = 2*Mtmp7;
double Mtmp68 = (1.0/12.0)*Mtmp19*y;
double Mtmp69 = pow(z, 4);
double Mtmp70 = (1.0/120.0)*Mtmp37;
double Mtmp71 = (1.0/120.0)*Mtmp63;
double Mtmp72 = (1.0/120.0)*Mtmp69;
#pragma omp atomic
M[0] += S[0];
#pragma omp atomic
M[1] += S[1];
#pragma omp atomic
M[2] += S[2];
#pragma omp atomic
M[3] += x*S[0];
#pragma omp atomic
M[4] += Mtmp2;
#pragma omp atomic
M[5] += Mtmp5;
#pragma omp atomic
M[6] += y*S[1];
#pragma omp atomic
M[7] += Mtmp8;
#pragma omp atomic
M[8] += z*S[2];
#pragma omp atomic
M[9] += (1.0/2.0)*Mtmp9*S[0];
#pragma omp atomic
M[10] += x*(Mtmp1 + Mtmp10);
#pragma omp atomic
M[11] += x*(Mtmp11 + Mtmp4);
#pragma omp atomic
M[12] += y*(Mtmp0 + Mtmp12);
#pragma omp atomic
M[13] += Mtmp17;
#pragma omp atomic
M[14] += z*(Mtmp3 + (1.0/2.0)*Mtmp4);
#pragma omp atomic
M[15] += (1.0/2.0)*Mtmp18*S[1];
#pragma omp atomic
M[16] += y*((1.0/2.0)*Mtmp6 + Mtmp7);
#pragma omp atomic
M[17] += z*(Mtmp6 + (1.0/2.0)*Mtmp7);
#pragma omp atomic
M[18] += (1.0/2.0)*Mtmp19*S[2];
#pragma omp atomic
M[19] += (1.0/6.0)*Mtmp20*S[0];
#pragma omp atomic
M[20] += Mtmp22*(Mtmp0 + Mtmp21);
#pragma omp atomic
M[21] += Mtmp22*(Mtmp23 + Mtmp3);
#pragma omp atomic
M[22] += Mtmp2*Mtmp24*y;
#pragma omp atomic
M[23] += x*(Mtmp13 + Mtmp25 + Mtmp26);
#pragma omp atomic
M[24] += Mtmp24*Mtmp5*z;
#pragma omp atomic
M[25] += Mtmp28*(Mtmp1 + Mtmp27);
#pragma omp atomic
M[26] += y*(Mtmp15 + Mtmp25 + Mtmp29);
#pragma omp atomic
M[27] += z*(Mtmp14 + Mtmp26 + Mtmp29);
#pragma omp atomic
M[28] += Mtmp31*(Mtmp30 + Mtmp4);
#pragma omp atomic
M[29] += (1.0/6.0)*Mtmp32*S[1];
#pragma omp atomic
M[30] += Mtmp28*(Mtmp33 + Mtmp6);
#pragma omp atomic
M[31] += (1.0/2.0)*Mtmp34*Mtmp8;
#pragma omp atomic
M[32] += Mtmp31*(Mtmp35 + Mtmp7);
#pragma omp atomic
M[33] += (1.0/6.0)*Mtmp36*S[2];
#pragma omp atomic
M[34] += (1.0/24.0)*Mtmp37*S[0];
#pragma omp atomic
M[35] += Mtmp38*(Mtmp0 + 4*Mtmp1);
#pragma omp atomic
M[36] += Mtmp38*(Mtmp3 + 4*Mtmp4);
#pragma omp atomic
M[37] += Mtmp41*(Mtmp21 + Mtmp39);
#pragma omp atomic
M[38] += Mtmp22*(Mtmp16 + Mtmp42);
#pragma omp atomic
M[39] += Mtmp44*(Mtmp23 + Mtmp43);
#pragma omp atomic
M[40] += Mtmp47*(Mtmp27 + Mtmp45);
#pragma omp atomic
M[41] += Mtmp51*y*(Mtmp48 + Mtmp50);
#pragma omp atomic
M[42] += Mtmp51*z*(Mtmp48 + Mtmp53);
#pragma omp atomic
M[43] += Mtmp55*(Mtmp30 + Mtmp54);
#pragma omp atomic
M[44] += Mtmp56*(4*Mtmp0 + Mtmp1);
#pragma omp atomic
M[45] += Mtmp28*(Mtmp13 + Mtmp58);
#pragma omp atomic
M[46] += (1.0/4.0)*Mtmp34*(Mtmp49 + Mtmp59);
#pragma omp atomic
M[47] += Mtmp31*(Mtmp13 + Mtmp61);
#pragma omp atomic
M[48] += Mtmp62*(4*Mtmp3 + Mtmp4);
#pragma omp atomic
M[49] += (1.0/24.0)*Mtmp63*S[1];
#pragma omp atomic
M[50] += Mtmp56*(Mtmp6 + 4*Mtmp7);
#pragma omp atomic
M[51] += Mtmp66*(Mtmp33 + Mtmp64);
#pragma omp atomic
M[52] += Mtmp68*(Mtmp35 + Mtmp67);
#pragma omp atomic
M[53] += Mtmp62*(4*Mtmp6 + Mtmp7);
#pragma omp atomic
M[54] += (1.0/24.0)*Mtmp69*S[2];
#pragma omp atomic
M[55] += (1.0/120.0)*pow(x, 5)*S[0];
#pragma omp atomic
M[56] += Mtmp70*(Mtmp0 + 5*Mtmp1);
#pragma omp atomic
M[57] += Mtmp70*(Mtmp3 + 5*Mtmp4);
#pragma omp atomic
M[58] += Mtmp38*y*(Mtmp0 + Mtmp45);
#pragma omp atomic
M[59] += Mtmp38*(4*Mtmp13 + Mtmp16);
#pragma omp atomic
M[60] += Mtmp38*z*(Mtmp3 + Mtmp54);
#pragma omp atomic
M[61] += Mtmp18*Mtmp2*Mtmp40;
#pragma omp atomic
M[62] += Mtmp41*(Mtmp42 + Mtmp50);
#pragma omp atomic
M[63] += Mtmp44*(Mtmp42 + Mtmp53);
#pragma omp atomic
M[64] += Mtmp19*Mtmp40*Mtmp5;
#pragma omp atomic
M[65] += Mtmp56*x*(Mtmp1 + Mtmp39);
#pragma omp atomic
M[66] += Mtmp47*(Mtmp48 + Mtmp58);
#pragma omp atomic
M[67] += Mtmp17*Mtmp34*Mtmp51;
#pragma omp atomic
M[68] += Mtmp55*(Mtmp48 + Mtmp61);
#pragma omp atomic
M[69] += Mtmp62*x*(Mtmp4 + Mtmp43);
#pragma omp atomic
M[70] += Mtmp71*(5*Mtmp0 + Mtmp1);
#pragma omp atomic
M[71] += Mtmp56*(Mtmp13 + Mtmp14 + 4*Mtmp15);
#pragma omp atomic
M[72] += Mtmp66*(Mtmp57 + Mtmp59);
#pragma omp atomic
M[73] += Mtmp68*(Mtmp13 + Mtmp49 + Mtmp60);
#pragma omp atomic
M[74] += Mtmp62*(Mtmp13 + 4*Mtmp14 + Mtmp15);
#pragma omp atomic
M[75] += Mtmp72*(5*Mtmp3 + Mtmp4);
#pragma omp atomic
M[76] += (1.0/120.0)*pow(y, 5)*S[1];
#pragma omp atomic
M[77] += Mtmp71*(Mtmp6 + 5*Mtmp7);
#pragma omp atomic
M[78] += Mtmp56*z*(Mtmp6 + Mtmp67);
#pragma omp atomic
M[79] += Mtmp19*Mtmp65*Mtmp8;
#pragma omp atomic
M[80] += Mtmp62*y*(Mtmp64 + Mtmp7);
#pragma omp atomic
M[81] += Mtmp72*(5*Mtmp6 + Mtmp7);
#pragma omp atomic
M[82] += (1.0/120.0)*pow(z, 5)*S[2];

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
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double D[83];
double Dtmp0 = pow(Rinv, 3);
double Dtmp1 = pow(x, 2);
double Dtmp2 = pow(Rinv, 2);
double Dtmp3 = 3*Dtmp2;
double Dtmp4 = Dtmp1*Dtmp3;
double Dtmp5 = Dtmp4 - 1;
double Dtmp6 = 3*pow(Rinv, 5);
double Dtmp7 = Dtmp6*x;
double Dtmp8 = pow(y, 2);
double Dtmp9 = Dtmp3*Dtmp8;
double Dtmp10 = Dtmp9 - 1;
double Dtmp11 = Dtmp6*y;
double Dtmp12 = 5*Dtmp2;
double Dtmp13 = Dtmp1*Dtmp12;
double Dtmp14 = Dtmp13 - 1;
double Dtmp15 = Dtmp6*z;
double Dtmp16 = Dtmp12*Dtmp8;
double Dtmp17 = Dtmp16 - 1;
double Dtmp18 = pow(Rinv, 7);
double Dtmp19 = 15*Dtmp18;
double Dtmp20 = Dtmp19*x;
double Dtmp21 = Dtmp20*y;
double Dtmp22 = Dtmp1*Dtmp2;
double Dtmp23 = -30*Dtmp22;
double Dtmp24 = pow(x, 4);
double Dtmp25 = pow(Rinv, 4);
double Dtmp26 = 35*Dtmp25;
double Dtmp27 = 7*Dtmp22;
double Dtmp28 = Dtmp27 - 3;
double Dtmp29 = Dtmp20*z;
double Dtmp30 = Dtmp1*Dtmp8;
double Dtmp31 = Dtmp27 - 1;
double Dtmp32 = Dtmp19*y;
double Dtmp33 = Dtmp32*z;
double Dtmp34 = Dtmp2*Dtmp8;
double Dtmp35 = 7*Dtmp34;
double Dtmp36 = Dtmp35 - 3;
double Dtmp37 = Dtmp35 - 1;
double Dtmp38 = -30*Dtmp34;
double Dtmp39 = pow(y, 4);
double Dtmp40 = Dtmp24*Dtmp25;
double Dtmp41 = 14*Dtmp22;
double Dtmp42 = 21*Dtmp40;
double Dtmp43 = 45*Dtmp18;
double Dtmp44 = Dtmp43*(-Dtmp41 + Dtmp42 + 1);
double Dtmp45 = -Dtmp27;
double Dtmp46 = Dtmp25*Dtmp30;
double Dtmp47 = 63*Dtmp46;
double Dtmp48 = Dtmp47 + 3;
double Dtmp49 = pow(Rinv, 9);
double Dtmp50 = 315*Dtmp49*y;
double Dtmp51 = Dtmp50*z;
double Dtmp52 = Dtmp51*x;
double Dtmp53 = -Dtmp35;
double Dtmp54 = 14*Dtmp34;
double Dtmp55 = Dtmp25*Dtmp39;
double Dtmp56 = 21*Dtmp55;
double Dtmp57 = -Dtmp54 + Dtmp56 + 1;
double Dtmp58 = Dtmp43*Dtmp57;
double Dtmp59 = 231*pow(Rinv, 6);
double Dtmp60 = 33*Dtmp40;
double Dtmp61 = x*(Dtmp23 + Dtmp60 + 5);
double Dtmp62 = 315*Dtmp49*z;
double Dtmp63 = -126*Dtmp46;
double Dtmp64 = -Dtmp9;
double Dtmp65 = 1 - Dtmp4;
double Dtmp66 = 33*Dtmp46;
double Dtmp67 = Dtmp62*x;
double Dtmp68 = Dtmp38 + 33*Dtmp55 + 5;
D[0] = -Dtmp0*x;
D[1] = -Dtmp0*y;
D[2] = -Dtmp0*z;
D[3] = Dtmp0*Dtmp5;
D[4] = Dtmp7*y;
D[5] = Dtmp7*z;
D[6] = Dtmp0*Dtmp10;
D[7] = Dtmp11*z;
D[8] = -D[3] - D[6];
D[9] = -Dtmp7*(Dtmp13 - 3);
D[10] = -Dtmp11*Dtmp14;
D[11] = -Dtmp14*Dtmp15;
D[12] = -Dtmp17*Dtmp7;
D[13] = -Dtmp21*z;
D[14] = -D[9] - D[12];
D[15] = -Dtmp11*(Dtmp16 - 3);
D[16] = -Dtmp15*Dtmp17;
D[17] = -D[10] - D[15];
D[18] = -D[11] - D[16];
D[19] = Dtmp6*(Dtmp23 + Dtmp24*Dtmp26 + 3);
D[20] = Dtmp21*Dtmp28;
D[21] = Dtmp28*Dtmp29;
D[22] = Dtmp6*(-Dtmp13 - Dtmp16 + Dtmp26*Dtmp30 + 1);
D[23] = Dtmp31*Dtmp33;
D[24] = -D[19] - D[22];
D[25] = Dtmp21*Dtmp36;
D[26] = Dtmp29*Dtmp37;
D[27] = -D[20] - D[25];
D[28] = -D[21] - D[26];
D[29] = Dtmp6*(Dtmp26*Dtmp39 + Dtmp38 + 3);
D[30] = Dtmp33*Dtmp36;
D[31] = -D[22] - D[29];
D[32] = -D[23] - D[30];
D[33] = -D[24] - D[31];
D[34] = -Dtmp20*(-70*Dtmp22 + 63*Dtmp40 + 15);
D[35] = -Dtmp44*y;
D[36] = -Dtmp44*z;
D[37] = -Dtmp20*(-21*Dtmp34 + Dtmp45 + Dtmp48);
D[38] = -Dtmp5*Dtmp52;
D[39] = -D[34] - D[37];
D[40] = -Dtmp32*(-21*Dtmp22 + Dtmp48 + Dtmp53);
D[41] = -Dtmp19*z*(Dtmp45 + Dtmp47 + Dtmp53 + 1);
D[42] = -D[35] - D[40];
D[43] = -D[36] - D[41];
D[44] = -Dtmp58*x;
D[45] = -Dtmp10*Dtmp52;
D[46] = -D[37] - D[44];
D[47] = -D[38] - D[45];
D[48] = -D[39] - D[46];
D[49] = -Dtmp32*(-70*Dtmp34 + 63*Dtmp55 + 15);
D[50] = -Dtmp58*z;
D[51] = -D[40] - D[49];
D[52] = -D[41] - D[50];
D[53] = -D[42] - D[51];
D[54] = -D[43] - D[52];
D[55] = Dtmp43*(105*Dtmp22 - 315*Dtmp40 + Dtmp59*pow(x, 6) - 5);
D[56] = Dtmp50*Dtmp61;
D[57] = Dtmp61*Dtmp62;
D[58] = Dtmp43*(Dtmp24*Dtmp59*Dtmp8 + Dtmp37 + Dtmp41 - Dtmp42 + Dtmp63);
D[59] = Dtmp51*(-18*Dtmp22 + Dtmp60 + 1);
D[60] = -D[55] - D[58];
D[61] = 945*Dtmp49*x*y*(11*Dtmp46 + Dtmp64 + Dtmp65);
D[62] = Dtmp67*(-9*Dtmp34 + Dtmp65 + Dtmp66);
D[63] = -D[56] - D[61];
D[64] = -D[57] - D[62];
D[65] = Dtmp43*(Dtmp1*Dtmp39*Dtmp59 + Dtmp31 + Dtmp54 - Dtmp56 + Dtmp63);
D[66] = Dtmp51*(-9*Dtmp22 + Dtmp64 + Dtmp66 + 1);
D[67] = -D[58] - D[65];
D[68] = -D[59] - D[66];
D[69] = -D[60] - D[67];
D[70] = Dtmp50*Dtmp68*x;
D[71] = Dtmp67*(4*Dtmp10*Dtmp34 + Dtmp57);
D[72] = -D[61] - D[70];
D[73] = -D[62] - D[71];
D[74] = -D[63] - D[72];
D[75] = -D[64] - D[73];
D[76] = Dtmp43*(105*Dtmp34 - 315*Dtmp55 + Dtmp59*pow(y, 6) - 5);
D[77] = Dtmp51*Dtmp68;
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
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(Rinv, 3);
double Ftmp1 = pow(Rinv, 2);
double Ftmp2 = 3*Ftmp1;
double Ftmp3 = y*M[4];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = Ftmp2*x;
double Ftmp6 = Ftmp5*y;
double Ftmp7 = Ftmp4*M[2];
double Ftmp8 = pow(Rinv, 4);
double Ftmp9 = 15*Ftmp8;
double Ftmp10 = Ftmp9*y;
double Ftmp11 = Ftmp10*M[13];
double Ftmp12 = pow(x, 2);
double Ftmp13 = Ftmp1*Ftmp12;
double Ftmp14 = 3*Ftmp13;
double Ftmp15 = Ftmp10*M[7];
double Ftmp16 = x*z;
double Ftmp17 = 6*x;
double Ftmp18 = pow(y, 2);
double Ftmp19 = Ftmp18*Ftmp8;
double Ftmp20 = Ftmp19*M[6];
double Ftmp21 = pow(z, 2);
double Ftmp22 = Ftmp21*Ftmp8;
double Ftmp23 = Ftmp22*M[8];
double Ftmp24 = Ftmp12*Ftmp9;
double Ftmp25 = Ftmp9*z;
double Ftmp26 = pow(Rinv, 6);
double Ftmp27 = 30*x;
double Ftmp28 = Ftmp26*Ftmp27;
double Ftmp29 = pow(y, 3);
double Ftmp30 = Ftmp29*M[15];
double Ftmp31 = pow(z, 3);
double Ftmp32 = Ftmp31*M[18];
double Ftmp33 = Ftmp28*y;
double Ftmp34 = Ftmp21*Ftmp33;
double Ftmp35 = Ftmp26*z;
double Ftmp36 = Ftmp27*Ftmp35;
double Ftmp37 = Ftmp12*Ftmp26;
double Ftmp38 = 105*Ftmp37;
double Ftmp39 = Ftmp38*y;
double Ftmp40 = pow(Rinv, 8);
double Ftmp41 = 210*Ftmp40;
double Ftmp42 = Ftmp31*Ftmp41;
double Ftmp43 = Ftmp42*M[32];
double Ftmp44 = x*y;
double Ftmp45 = Ftmp29*M[30];
double Ftmp46 = Ftmp16*Ftmp41;
double Ftmp47 = 30*Ftmp37;
double Ftmp48 = Ftmp18*M[12];
double Ftmp49 = Ftmp21*Ftmp47;
double Ftmp50 = Ftmp12*Ftmp29;
double Ftmp51 = Ftmp41*M[25];
double Ftmp52 = Ftmp42*M[28];
double Ftmp53 = Ftmp21*Ftmp41;
double Ftmp54 = Ftmp12*y;
double Ftmp55 = Ftmp53*Ftmp54;
double Ftmp56 = Ftmp18*Ftmp41;
double Ftmp57 = Ftmp12*Ftmp56*z;
double Ftmp58 = Ftmp54*M[47];
double Ftmp59 = pow(Rinv, 10);
double Ftmp60 = 1890*Ftmp59;
double Ftmp61 = Ftmp31*Ftmp60;
double Ftmp62 = z*M[45];
double Ftmp63 = 5*Ftmp13;
double Ftmp64 = (Ftmp63 - 3)*M[9];
double Ftmp65 = Ftmp1*Ftmp18;
double Ftmp66 = 5*Ftmp65;
double Ftmp67 = Ftmp66 - 1;
double Ftmp68 = Ftmp2*Ftmp67;
double Ftmp69 = Ftmp1*Ftmp21;
double Ftmp70 = 5*Ftmp69;
double Ftmp71 = Ftmp70 - 1;
double Ftmp72 = Ftmp2*Ftmp71;
double Ftmp73 = 6*Ftmp1;
double Ftmp74 = Ftmp13 - 1;
double Ftmp75 = Ftmp74*x;
double Ftmp76 = Ftmp14 - 1;
double Ftmp77 = Ftmp76*M[3];
double Ftmp78 = 3*Ftmp18;
double Ftmp79 = Ftmp1*Ftmp78;
double Ftmp80 = Ftmp79 - 1;
double Ftmp81 = Ftmp80*M[6];
double Ftmp82 = 3*Ftmp21;
double Ftmp83 = Ftmp1*Ftmp82;
double Ftmp84 = Ftmp83 - 1;
double Ftmp85 = Ftmp84*M[8];
double Ftmp86 = 7*Ftmp13;
double Ftmp87 = Ftmp86 - 3;
double Ftmp88 = Ftmp87*M[20];
double Ftmp89 = 7*Ftmp65;
double Ftmp90 = Ftmp89 - 3;
double Ftmp91 = Ftmp90*M[25];
double Ftmp92 = 7*Ftmp69;
double Ftmp93 = Ftmp92 - 1;
double Ftmp94 = Ftmp93*M[27];
double Ftmp95 = Ftmp87*M[21];
double Ftmp96 = Ftmp89 - 1;
double Ftmp97 = Ftmp96*M[26];
double Ftmp98 = Ftmp92 - 3;
double Ftmp99 = Ftmp25*Ftmp98;
double Ftmp100 = Ftmp74*Ftmp8;
double Ftmp101 = Ftmp100*Ftmp27;
double Ftmp102 = y*M[10];
double Ftmp103 = z*M[11];
double Ftmp104 = 30*Ftmp12;
double Ftmp105 = Ftmp10*x;
double Ftmp106 = Ftmp63 - 1;
double Ftmp107 = Ftmp106*M[10];
double Ftmp108 = (Ftmp66 - 3)*M[15];
double Ftmp109 = Ftmp71*M[17];
double Ftmp110 = Ftmp25*x;
double Ftmp111 = Ftmp67*M[16];
double Ftmp112 = (Ftmp70 - 3)*M[18];
double Ftmp113 = 315*Ftmp26;
double Ftmp114 = Ftmp113*y;
double Ftmp115 = Ftmp114*z;
double Ftmp116 = Ftmp76*M[38];
double Ftmp117 = Ftmp84*M[47];
double Ftmp118 = 210*Ftmp35;
double Ftmp119 = Ftmp118*y;
double Ftmp120 = Ftmp67*M[12];
double Ftmp121 = Ftmp71*M[14];
double Ftmp122 = 210*Ftmp37;
double Ftmp123 = Ftmp122*Ftmp74;
double Ftmp124 = z*M[21];
double Ftmp125 = Ftmp86 - 1;
double Ftmp126 = Ftmp125*M[23];
double Ftmp127 = 105*Ftmp26;
double Ftmp128 = Ftmp127*x;
double Ftmp129 = Ftmp128*y;
double Ftmp130 = Ftmp129*z;
double Ftmp131 = Ftmp90*M[30];
double Ftmp132 = Ftmp98*M[32];
double Ftmp133 = Ftmp26*x;
double Ftmp134 = 60*Ftmp133;
double Ftmp135 = Ftmp18*Ftmp90*M[29];
double Ftmp136 = Ftmp21*Ftmp98*M[33];
double Ftmp137 = Ftmp38*z;
double Ftmp138 = Ftmp98*M[28];
double Ftmp139 = Ftmp40*x;
double Ftmp140 = 420*Ftmp139;
double Ftmp141 = 9*Ftmp65;
double Ftmp142 = Ftmp29*(Ftmp141 - 5)*M[49];
double Ftmp143 = 9*Ftmp69;
double Ftmp144 = Ftmp31*(Ftmp143 - 5)*M[54];
double Ftmp145 = Ftmp40*z;
double Ftmp146 = 1890*Ftmp145;
double Ftmp147 = 1260*M[53];
double Ftmp148 = Ftmp21*Ftmp84;
double Ftmp149 = Ftmp139*y;
double Ftmp150 = Ftmp148*Ftmp149;
double Ftmp151 = Ftmp18*Ftmp80;
double Ftmp152 = 1260*M[50];
double Ftmp153 = Ftmp145*x;
double Ftmp154 = 2835*Ftmp145;
double Ftmp155 = Ftmp154*Ftmp54;
double Ftmp156 = Ftmp80*M[45];
double Ftmp157 = 11*Ftmp69;
double Ftmp158 = Ftmp31*(Ftmp157 - 5);
double Ftmp159 = 3780*Ftmp59;
double Ftmp160 = Ftmp158*Ftmp159*M[81];
double Ftmp161 = Ftmp16*Ftmp29;
double Ftmp162 = 11*Ftmp65;
double Ftmp163 = Ftmp59*(Ftmp162 - 5);
double Ftmp164 = 3780*M[77];
double Ftmp165 = Ftmp163*Ftmp164;
double Ftmp166 = Ftmp12*Ftmp40;
double Ftmp167 = 1260*M[44];
double Ftmp168 = Ftmp151*Ftmp167;
double Ftmp169 = 1260*Ftmp166;
double Ftmp170 = Ftmp148*M[48];
double Ftmp171 = 3780*M[70];
double Ftmp172 = Ftmp163*Ftmp171;
double Ftmp173 = Ftmp12*Ftmp59;
double Ftmp174 = 3780*M[75];
double Ftmp175 = Ftmp158*Ftmp174;
double Ftmp176 = Ftmp159*(Ftmp157 - 3)*M[74];
double Ftmp177 = Ftmp21*Ftmp54;
double Ftmp178 = Ftmp173*z;
double Ftmp179 = Ftmp178*Ftmp18;
double Ftmp180 = (Ftmp162 - 3)*M[71];
double Ftmp181 = pow(x, 4);
double Ftmp182 = Ftmp181*Ftmp8;
double Ftmp183 = 63*Ftmp182;
double Ftmp184 = (-70*Ftmp13 + Ftmp183 + 15)*M[34];
double Ftmp185 = 45*Ftmp8;
double Ftmp186 = pow(y, 4);
double Ftmp187 = 21*Ftmp8;
double Ftmp188 = Ftmp186*Ftmp187;
double Ftmp189 = 14*Ftmp65;
double Ftmp190 = -Ftmp189;
double Ftmp191 = Ftmp190 + 1;
double Ftmp192 = Ftmp188 + Ftmp191;
double Ftmp193 = Ftmp192*M[44];
double Ftmp194 = pow(z, 4);
double Ftmp195 = Ftmp187*Ftmp194;
double Ftmp196 = 14*Ftmp69;
double Ftmp197 = -Ftmp196;
double Ftmp198 = Ftmp197 + 1;
double Ftmp199 = Ftmp195 + Ftmp198;
double Ftmp200 = Ftmp185*Ftmp199;
double Ftmp201 = 10*Ftmp13;
double Ftmp202 = -Ftmp201;
double Ftmp203 = Ftmp202 + 3;
double Ftmp204 = 60*Ftmp8;
double Ftmp205 = Ftmp9*x;
double Ftmp206 = -30*Ftmp13;
double Ftmp207 = (35*Ftmp182 + Ftmp206 + 3)*M[19];
double Ftmp208 = -30*Ftmp65;
double Ftmp209 = 35*Ftmp8;
double Ftmp210 = (Ftmp186*Ftmp209 + Ftmp208 + 3)*M[29];
double Ftmp211 = -30*Ftmp69;
double Ftmp212 = (Ftmp194*Ftmp209 + Ftmp211 + 3)*M[33];
double Ftmp213 = 33*Ftmp182;
double Ftmp214 = Ftmp206 + Ftmp213 + 5;
double Ftmp215 = Ftmp214*M[56];
double Ftmp216 = 33*Ftmp8;
double Ftmp217 = Ftmp186*Ftmp216 + Ftmp208 + 5;
double Ftmp218 = Ftmp217*M[70];
double Ftmp219 = Ftmp194*Ftmp216;
double Ftmp220 = (Ftmp219 - 18*Ftmp69 + 1)*M[74];
double Ftmp221 = Ftmp113*z;
double Ftmp222 = Ftmp214*M[57];
double Ftmp223 = Ftmp211 + Ftmp219 + 5;
double Ftmp224 = Ftmp221*Ftmp223;
double Ftmp225 = -5040*Ftmp13 + 3780*Ftmp182 + 1260;
double Ftmp226 = Ftmp133*y;
double Ftmp227 = 21*Ftmp182;
double Ftmp228 = 14*Ftmp13;
double Ftmp229 = -Ftmp228;
double Ftmp230 = Ftmp229 + 1;
double Ftmp231 = Ftmp227 + Ftmp230;
double Ftmp232 = Ftmp231*M[35];
double Ftmp233 = Ftmp114*x;
double Ftmp234 = 63*Ftmp8;
double Ftmp235 = Ftmp186*Ftmp234;
double Ftmp236 = (Ftmp235 - 70*Ftmp65 + 15)*M[49];
double Ftmp237 = Ftmp199*M[53];
double Ftmp238 = Ftmp35*x;
double Ftmp239 = Ftmp231*M[36];
double Ftmp240 = Ftmp221*x;
double Ftmp241 = Ftmp192*M[50];
double Ftmp242 = Ftmp194*Ftmp234;
double Ftmp243 = (Ftmp242 - 70*Ftmp69 + 15)*M[54];
double Ftmp244 = Ftmp128*z;
double Ftmp245 = 420*M[34];
double Ftmp246 = 315*Ftmp37;
double Ftmp247 = Ftmp199*M[48];
double Ftmp248 = 11*Ftmp182;
double Ftmp249 = Ftmp229 + 3;
double Ftmp250 = Ftmp145*Ftmp44;
double Ftmp251 = 3780*M[59];
double Ftmp252 = 18*Ftmp13;
double Ftmp253 = (Ftmp213 - Ftmp252 + 1)*M[59];
double Ftmp254 = Ftmp154*Ftmp44;
double Ftmp255 = Ftmp217*M[77];
double Ftmp256 = Ftmp223*M[81];
double Ftmp257 = Ftmp21*Ftmp65;
double Ftmp258 = -18*Ftmp257;
double Ftmp259 = Ftmp41*Ftmp44;
double Ftmp260 = Ftmp139*Ftmp18;
double Ftmp261 = 1890*Ftmp260;
double Ftmp262 = Ftmp217*M[76];
double Ftmp263 = Ftmp139*Ftmp21;
double Ftmp264 = 1890*M[82];
double Ftmp265 = Ftmp223*Ftmp264;
double Ftmp266 = Ftmp40*Ftmp54;
double Ftmp267 = -60480*Ftmp13 + 3780*Ftmp248 + 18900;
double Ftmp268 = 2835*Ftmp266;
double Ftmp269 = Ftmp12*Ftmp145;
double Ftmp270 = 2835*Ftmp269;
double Ftmp271 = Ftmp223*M[75];
double Ftmp272 = Ftmp60*z;
double Ftmp273 = Ftmp272*Ftmp44;
double Ftmp274 = pow(x, 6);
double Ftmp275 = 33*Ftmp26;
double Ftmp276 = 1890*Ftmp26;
double Ftmp277 = Ftmp113*x;
double Ftmp278 = 231*Ftmp26;
double Ftmp279 = (105*Ftmp13 - 315*Ftmp182 + Ftmp274*Ftmp278 - 5)*M[55];
double Ftmp280 = 315*Ftmp8;
double Ftmp281 = pow(y, 6);
double Ftmp282 = (-Ftmp186*Ftmp280 + Ftmp278*Ftmp281 + 105*Ftmp65 - 5)*M[76];
double Ftmp283 = pow(z, 6);
double Ftmp284 = (-Ftmp194*Ftmp280 + Ftmp278*Ftmp283 + 105*Ftmp69 - 5)*M[82];
double Ftmp285 = -21*Ftmp65;
double Ftmp286 = Ftmp12*Ftmp234;
double Ftmp287 = Ftmp18*Ftmp286;
double Ftmp288 = -Ftmp86;
double Ftmp289 = Ftmp288 + 3;
double Ftmp290 = (Ftmp285 + Ftmp287 + Ftmp289)*M[37];
double Ftmp291 = 21*Ftmp69;
double Ftmp292 = -Ftmp291;
double Ftmp293 = Ftmp21*Ftmp286;
double Ftmp294 = (Ftmp289 + Ftmp292 + Ftmp293)*M[39];
double Ftmp295 = 8*Ftmp65;
double Ftmp296 = -Ftmp295;
double Ftmp297 = 14*Ftmp12;
double Ftmp298 = Ftmp19*Ftmp297;
double Ftmp299 = -Ftmp13;
double Ftmp300 = Ftmp299 + 1;
double Ftmp301 = Ftmp27*Ftmp8;
double Ftmp302 = -Ftmp66;
double Ftmp303 = Ftmp12*Ftmp209;
double Ftmp304 = 1 - Ftmp63;
double Ftmp305 = (Ftmp18*Ftmp303 + Ftmp302 + Ftmp304)*M[22];
double Ftmp306 = 8*Ftmp69;
double Ftmp307 = -Ftmp306;
double Ftmp308 = Ftmp22*Ftmp297;
double Ftmp309 = -Ftmp70;
double Ftmp310 = (Ftmp21*Ftmp303 + Ftmp304 + Ftmp309)*M[24];
double Ftmp311 = Ftmp18*Ftmp21;
double Ftmp312 = (Ftmp209*Ftmp311 + Ftmp302 + Ftmp309 + 1)*M[31];
double Ftmp313 = 945*Ftmp26;
double Ftmp314 = -Ftmp79;
double Ftmp315 = 11*Ftmp12;
double Ftmp316 = -Ftmp14;
double Ftmp317 = Ftmp316 + 1;
double Ftmp318 = Ftmp19*Ftmp315 + Ftmp314 + Ftmp317;
double Ftmp319 = y*M[61];
double Ftmp320 = Ftmp318*Ftmp319;
double Ftmp321 = Ftmp12*Ftmp216;
double Ftmp322 = Ftmp21*Ftmp321;
double Ftmp323 = (-Ftmp143 + Ftmp317 + Ftmp322)*M[63];
double Ftmp324 = Ftmp18*Ftmp321;
double Ftmp325 = (-Ftmp141 + Ftmp317 + Ftmp324)*M[62];
double Ftmp326 = Ftmp313*z;
double Ftmp327 = -Ftmp83;
double Ftmp328 = (Ftmp22*Ftmp315 + Ftmp317 + Ftmp327)*M[64];
double Ftmp329 = 18*Ftmp12;
double Ftmp330 = Ftmp19*Ftmp329;
double Ftmp331 = 10*Ftmp65;
double Ftmp332 = -Ftmp331;
double Ftmp333 = Ftmp332 + 3;
double Ftmp334 = 210*Ftmp26;
double Ftmp335 = Ftmp334*Ftmp44;
double Ftmp336 = -21*Ftmp13;
double Ftmp337 = -Ftmp89;
double Ftmp338 = Ftmp337 + 3;
double Ftmp339 = (Ftmp287 + Ftmp336 + Ftmp338)*M[40];
double Ftmp340 = 10*Ftmp69;
double Ftmp341 = -Ftmp340;
double Ftmp342 = Ftmp22*Ftmp329;
double Ftmp343 = Ftmp300 + Ftmp342;
double Ftmp344 = Ftmp288 + 1;
double Ftmp345 = -Ftmp92;
double Ftmp346 = Ftmp293 + Ftmp345;
double Ftmp347 = (Ftmp344 + Ftmp346)*M[42];
double Ftmp348 = Ftmp234*Ftmp311;
double Ftmp349 = (Ftmp292 + Ftmp338 + Ftmp348)*M[51];
double Ftmp350 = Ftmp300 + Ftmp330;
double Ftmp351 = Ftmp118*x;
double Ftmp352 = (Ftmp287 + Ftmp337 + Ftmp344)*M[41];
double Ftmp353 = Ftmp341 + 3;
double Ftmp354 = (Ftmp336 + Ftmp346 + 3)*M[43];
double Ftmp355 = (Ftmp285 + Ftmp345 + Ftmp348 + 3)*M[52];
double Ftmp356 = -12*Ftmp65;
double Ftmp357 = -12*Ftmp69;
double Ftmp358 = 22*Ftmp12;
double Ftmp359 = Ftmp19*Ftmp358;
double Ftmp360 = Ftmp316 + 3;
double Ftmp361 = Ftmp359 + Ftmp360;
double Ftmp362 = Ftmp146*Ftmp44;
double Ftmp363 = 9*Ftmp13;
double Ftmp364 = 1 - Ftmp363;
double Ftmp365 = (Ftmp314 + Ftmp324 + Ftmp364)*M[66];
double Ftmp366 = Ftmp22*Ftmp358;
double Ftmp367 = (Ftmp322 + Ftmp327 + Ftmp364)*M[68];
double Ftmp368 = Ftmp19*Ftmp21;
double Ftmp369 = (Ftmp314 + Ftmp327 + 11*Ftmp368 + 1)*M[79];
double Ftmp370 = Ftmp369*y;
double Ftmp371 = 8505*Ftmp370;
double Ftmp372 = 8505*Ftmp320;
double Ftmp373 = 1890*M[63];
double Ftmp374 = 1890*Ftmp269;
double Ftmp375 = 8505*Ftmp328;
double Ftmp376 = Ftmp197 + 3;
double Ftmp377 = 4*Ftmp65;
double Ftmp378 = (Ftmp192 + Ftmp377*Ftmp80)*M[71];
double Ftmp379 = -36*Ftmp257;
double Ftmp380 = Ftmp194*Ftmp8;
double Ftmp381 = 99*Ftmp380;
double Ftmp382 = 2*Ftmp21;
double Ftmp383 = -Ftmp194*Ftmp73 + Ftmp382;
double Ftmp384 = 630*x;
double Ftmp385 = Ftmp384*Ftmp40;
double Ftmp386 = Ftmp186*Ftmp8;
double Ftmp387 = 99*Ftmp21;
double Ftmp388 = 2*Ftmp18;
double Ftmp389 = -Ftmp186*Ftmp73 + Ftmp388;
double Ftmp390 = 36*Ftmp12;
double Ftmp391 = -Ftmp19*Ftmp390;
double Ftmp392 = 99*Ftmp26;
double Ftmp393 = Ftmp12*Ftmp392;
double Ftmp394 = Ftmp186*Ftmp393;
double Ftmp395 = 39*Ftmp8;
double Ftmp396 = 20*Ftmp65;
double Ftmp397 = -Ftmp186*Ftmp395 + Ftmp396;
double Ftmp398 = Ftmp26*Ftmp384;
double Ftmp399 = -Ftmp22*Ftmp390;
double Ftmp400 = Ftmp194*Ftmp393;
double Ftmp401 = -Ftmp194*Ftmp395 + 20*Ftmp69;
double Ftmp402 = 126*Ftmp12;
double Ftmp403 = -Ftmp19*Ftmp402;
double Ftmp404 = Ftmp186*Ftmp278;
double Ftmp405 = -Ftmp188 + Ftmp189;
double Ftmp406 = (Ftmp12*Ftmp404 + Ftmp125 + Ftmp403 + Ftmp405)*M[65];
double Ftmp407 = -Ftmp22*Ftmp402;
double Ftmp408 = Ftmp194*Ftmp278;
double Ftmp409 = -Ftmp195 + Ftmp196;
double Ftmp410 = (Ftmp12*Ftmp408 + Ftmp125 + Ftmp407 + Ftmp409)*M[69];
double Ftmp411 = 19*Ftmp65;
double Ftmp412 = Ftmp181*Ftmp392;
double Ftmp413 = Ftmp18*Ftmp412;
double Ftmp414 = 102*Ftmp12;
double Ftmp415 = -Ftmp19*Ftmp414 - 2;
double Ftmp416 = 8*Ftmp13;
double Ftmp417 = -6*Ftmp182 + Ftmp416;
double Ftmp418 = Ftmp181*Ftmp278;
double Ftmp419 = -Ftmp227 + Ftmp228;
double Ftmp420 = (Ftmp18*Ftmp418 + Ftmp403 + Ftmp419 + Ftmp96)*M[58];
double Ftmp421 = -Ftmp22*Ftmp414;
double Ftmp422 = Ftmp21*Ftmp412;
double Ftmp423 = 19*Ftmp69 - 2;
double Ftmp424 = (Ftmp21*Ftmp418 + Ftmp407 + Ftmp419 + Ftmp93)*M[60];
double Ftmp425 = 126*Ftmp368;
double Ftmp426 = -Ftmp425;
double Ftmp427 = (Ftmp18*Ftmp408 + Ftmp409 + Ftmp426 + Ftmp96)*M[80];
double Ftmp428 = (Ftmp21*Ftmp404 + Ftmp405 + Ftmp426 + Ftmp93)*M[78];
double Ftmp429 = 4*Ftmp69;
double Ftmp430 = Ftmp429 - 1;
double Ftmp431 = Ftmp388*Ftmp430;
double Ftmp432 = Ftmp1*Ftmp431;
double Ftmp433 = (20*Ftmp368 + Ftmp432 + Ftmp67*Ftmp93)*M[46];
double Ftmp434 = Ftmp143 - 1;
double Ftmp435 = 28*Ftmp368 + Ftmp432 + Ftmp434*Ftmp90;
double Ftmp436 = y*M[72];
double Ftmp437 = Ftmp435*Ftmp436;
double Ftmp438 = 945*Ftmp437;
double Ftmp439 = Ftmp18*Ftmp93;
double Ftmp440 = 7*Ftmp67;
double Ftmp441 = 7*Ftmp434;
double Ftmp442 = 9*Ftmp90;
double Ftmp443 = Ftmp311*Ftmp37;
double Ftmp444 = 297*Ftmp443;
double Ftmp445 = -Ftmp330 + Ftmp340 + Ftmp444;
double Ftmp446 = Ftmp331 - Ftmp342;
double Ftmp447 = x*M[67];
double Ftmp448 = Ftmp127*(-Ftmp287 - Ftmp293 - Ftmp348 + 693*Ftmp443 + Ftmp86 + Ftmp89 + Ftmp93);
double Ftmp449 = Ftmp1*Ftmp382 - 1;
double Ftmp450 = Ftmp295*Ftmp449 + Ftmp331*Ftmp430 + Ftmp331*Ftmp93;
double Ftmp451 = (Ftmp440*Ftmp84 + Ftmp450)*M[73];
double Ftmp452 = 45*Ftmp451;
double Ftmp453 = 405*Ftmp451;
double Ftmp454 = 35*Ftmp84;
double Ftmp455 = 21*Ftmp67;
double Ftmp456 = 90*M[73];
double Ftmp457 = 6*y;
double Ftmp458 = Ftmp12*Ftmp8*M[3];
double Ftmp459 = 30*Ftmp26;
double Ftmp460 = pow(x, 3);
double Ftmp461 = Ftmp460*y;
double Ftmp462 = Ftmp459*y;
double Ftmp463 = Ftmp41*Ftmp461;
double Ftmp464 = Ftmp18*M[10];
double Ftmp465 = Ftmp311*Ftmp459;
double Ftmp466 = Ftmp311*x;
double Ftmp467 = Ftmp41*Ftmp466;
double Ftmp468 = Ftmp18*x;
double Ftmp469 = Ftmp18*Ftmp460;
double Ftmp470 = Ftmp60*M[38];
double Ftmp471 = Ftmp106*Ftmp2;
double Ftmp472 = Ftmp65 - 1;
double Ftmp473 = Ftmp472*y;
double Ftmp474 = Ftmp2*y;
double Ftmp475 = 30*Ftmp472;
double Ftmp476 = Ftmp8*y;
double Ftmp477 = Ftmp476*z;
double Ftmp478 = Ftmp10*z;
double Ftmp479 = Ftmp106*Ftmp9;
double Ftmp480 = Ftmp18*Ftmp9;
double Ftmp481 = Ftmp18*Ftmp334;
double Ftmp482 = Ftmp128*Ftmp18;
double Ftmp483 = 60*y;
double Ftmp484 = Ftmp37*Ftmp87*M[19];
double Ftmp485 = Ftmp127*Ftmp18;
double Ftmp486 = Ftmp485*z;
double Ftmp487 = Ftmp245*(Ftmp363 - 5);
double Ftmp488 = Ftmp40*y;
double Ftmp489 = Ftmp154*Ftmp468;
double Ftmp490 = x*M[45];
double Ftmp491 = 1260*Ftmp145*Ftmp54*Ftmp76;
double Ftmp492 = 11*Ftmp13;
double Ftmp493 = Ftmp159*(Ftmp492 - 5);
double Ftmp494 = Ftmp493*M[57];
double Ftmp495 = Ftmp461*z;
double Ftmp496 = Ftmp169*Ftmp76;
double Ftmp497 = Ftmp311*Ftmp40;
double Ftmp498 = Ftmp493*M[56];
double Ftmp499 = Ftmp492 - 3;
double Ftmp500 = Ftmp185*Ftmp231;
double Ftmp501 = Ftmp12 + Ftmp21;
double Ftmp502 = -Ftmp377 + 3*Ftmp386 + 1;
double Ftmp503 = Ftmp35*y;
double Ftmp504 = y*z;
double Ftmp505 = Ftmp127*Ftmp504;
double Ftmp506 = Ftmp113*Ftmp18;
double Ftmp507 = 420*Ftmp26;
double Ftmp508 = -Ftmp21*Ftmp252;
double Ftmp509 = Ftmp12 + Ftmp82;
double Ftmp510 = 2835*Ftmp260;
double Ftmp511 = 11*Ftmp386 - 16*Ftmp65 + 5;
double Ftmp512 = Ftmp214*M[55];
double Ftmp513 = 3*Ftmp12;
double Ftmp514 = Ftmp21 + Ftmp513;
double Ftmp515 = Ftmp41*Ftmp504;
double Ftmp516 = Ftmp21*Ftmp488;
double Ftmp517 = Ftmp145*Ftmp18;
double Ftmp518 = 2835*Ftmp517;
double Ftmp519 = Ftmp13*Ftmp21;
double Ftmp520 = -22*Ftmp519;
double Ftmp521 = -Ftmp416;
double Ftmp522 = -Ftmp65;
double Ftmp523 = Ftmp522 + 1;
double Ftmp524 = 30*M[22];
double Ftmp525 = 14*Ftmp368;
double Ftmp526 = 30*M[31];
double Ftmp527 = Ftmp313*x;
double Ftmp528 = Ftmp318*M[61];
double Ftmp529 = Ftmp330 + Ftmp523;
double Ftmp530 = 18*Ftmp368;
double Ftmp531 = -12*Ftmp13;
double Ftmp532 = Ftmp314 + Ftmp359;
double Ftmp533 = Ftmp531 + 3;
double Ftmp534 = 1890*Ftmp517;
double Ftmp535 = 22*Ftmp368;
double Ftmp536 = -36*Ftmp519;
double Ftmp537 = 630*y;
double Ftmp538 = Ftmp40*Ftmp537;
double Ftmp539 = 2*Ftmp12 - Ftmp181*Ftmp73;
double Ftmp540 = 20*Ftmp13 - 39*Ftmp182;
double Ftmp541 = Ftmp26*Ftmp537;
double Ftmp542 = -36*Ftmp368;
double Ftmp543 = Ftmp18*Ftmp194*Ftmp392;
double Ftmp544 = 19*Ftmp13;
double Ftmp545 = Ftmp295 - 6*Ftmp386;
double Ftmp546 = -102*Ftmp368;
double Ftmp547 = Ftmp186*Ftmp21*Ftmp392;
double Ftmp548 = Ftmp435*M[72];
double Ftmp549 = -Ftmp472;
double Ftmp550 = 117*Ftmp12;
double Ftmp551 = Ftmp201 - Ftmp530;
double Ftmp552 = y*M[67];
double Ftmp553 = Ftmp432 + 2;
double Ftmp554 = 6*z;
double Ftmp555 = 30*Ftmp35;
double Ftmp556 = Ftmp21*Ftmp60;
double Ftmp557 = Ftmp69 - 1;
double Ftmp558 = Ftmp557*z;
double Ftmp559 = 30*Ftmp557;
double Ftmp560 = Ftmp21*Ftmp9;
double Ftmp561 = Ftmp21*Ftmp334;
double Ftmp562 = Ftmp557*Ftmp561;
double Ftmp563 = Ftmp128*Ftmp21;
double Ftmp564 = Ftmp127*Ftmp21;
double Ftmp565 = Ftmp564*y;
double Ftmp566 = 1890*Ftmp263;
double Ftmp567 = 2835*Ftmp263;
double Ftmp568 = Ftmp567*y;
double Ftmp569 = Ftmp12 + Ftmp18;
double Ftmp570 = 3*Ftmp380 - Ftmp429 + 1;
double Ftmp571 = Ftmp113*Ftmp21;
double Ftmp572 = -Ftmp18*Ftmp252;
double Ftmp573 = 11*Ftmp380;
double Ftmp574 = Ftmp12 + Ftmp78;
double Ftmp575 = Ftmp573 - 16*Ftmp69 + 5;
double Ftmp576 = Ftmp18 + Ftmp513;
double Ftmp577 = 2835*Ftmp516;
double Ftmp578 = Ftmp13*Ftmp18;
double Ftmp579 = -22*Ftmp578;
double Ftmp580 = -Ftmp69;
double Ftmp581 = Ftmp580 + 1;
double Ftmp582 = Ftmp8*z;
double Ftmp583 = Ftmp342 + Ftmp581;
double Ftmp584 = Ftmp327 + Ftmp366;
double Ftmp585 = 1890*Ftmp516;
double Ftmp586 = -36*Ftmp578;
double Ftmp587 = 630*Ftmp145;
double Ftmp588 = 630*Ftmp35;
double Ftmp589 = Ftmp306 - 6*Ftmp380 - 2;
double Ftmp590 = -Ftmp557;
double Ftmp591 = -Ftmp295*Ftmp590 + Ftmp432;
#pragma omp atomic
F[0] += Ftmp0*(-Ftmp10*Ftmp88 - Ftmp10*Ftmp91 - Ftmp10*Ftmp94 - Ftmp100*Ftmp104*M[9] - Ftmp101*Ftmp102 - Ftmp101*Ftmp103 - Ftmp104*Ftmp40*(Ftmp21*Ftmp440 + 48*Ftmp257 + Ftmp431 + 5*Ftmp439)*M[46] - Ftmp105*Ftmp107 - Ftmp105*Ftmp108 - Ftmp105*Ftmp109 - Ftmp106*Ftmp110*M[11] + Ftmp11*z - Ftmp110*Ftmp111 - Ftmp110*Ftmp112 - Ftmp114*Ftmp215 - Ftmp114*Ftmp218 - Ftmp114*Ftmp220 - Ftmp114*Ftmp323 + Ftmp114*Ftmp62*Ftmp80 + Ftmp115*Ftmp116 + Ftmp115*Ftmp117 - Ftmp116*Ftmp155 + Ftmp119*Ftmp75*M[23] + Ftmp12*Ftmp25*M[5] + Ftmp12*Ftmp52 - Ftmp120*Ftmp24 - Ftmp121*Ftmp24 - Ftmp122*(Ftmp343 + Ftmp357)*M[39] - Ftmp122*(Ftmp350 + Ftmp356)*M[37] + Ftmp123*Ftmp124 + Ftmp123*y*M[20] + Ftmp126*Ftmp130 - Ftmp127*Ftmp437 - Ftmp129*Ftmp236 - Ftmp129*Ftmp339 - Ftmp129*Ftmp347 - Ftmp129*Ftmp349 + Ftmp130*Ftmp131 + Ftmp130*Ftmp132 + Ftmp134*Ftmp135 + Ftmp134*Ftmp136 + Ftmp137*Ftmp138 + Ftmp137*Ftmp95 + Ftmp137*Ftmp97 - Ftmp14*M[0] - Ftmp140*Ftmp142 - Ftmp140*Ftmp144 - Ftmp146*Ftmp54*Ftmp74*M[38] - Ftmp147*Ftmp150 + Ftmp15*Ftmp16 - Ftmp151*Ftmp152*Ftmp153 + Ftmp153*Ftmp371 - Ftmp154*Ftmp58*Ftmp84 - Ftmp155*Ftmp156 + Ftmp160*Ftmp44 + Ftmp161*Ftmp165 - Ftmp166*Ftmp168 + 1890*Ftmp166*Ftmp319*(Ftmp190 + Ftmp361) + Ftmp166*Ftmp372 + Ftmp166*Ftmp438 - Ftmp169*Ftmp170 + Ftmp17*Ftmp20 + Ftmp17*Ftmp23 + Ftmp172*Ftmp50 + Ftmp173*Ftmp175 + 210*Ftmp173*Ftmp436*(Ftmp18*Ftmp441 + Ftmp21*Ftmp442 + 64*Ftmp257 + Ftmp431) + Ftmp176*Ftmp177 + Ftmp178*Ftmp456*(10*Ftmp18*Ftmp430 + 8*Ftmp18*Ftmp449 + Ftmp18*Ftmp454 + Ftmp21*Ftmp455 + 126*Ftmp257 + 10*Ftmp439) + 3780*Ftmp179*Ftmp180 - Ftmp18*Ftmp36*M[16] - Ftmp184*Ftmp38 + Ftmp184*Ftmp9 + Ftmp185*Ftmp193 - Ftmp193*Ftmp246 - Ftmp2*Ftmp3 + Ftmp2*Ftmp64 + Ftmp200*M[48] + Ftmp204*x*(7*Ftmp182 + Ftmp203)*M[19] + Ftmp205*Ftmp207 + Ftmp205*Ftmp210 + Ftmp205*Ftmp212 + Ftmp205*Ftmp305 + Ftmp205*Ftmp310 + Ftmp205*Ftmp312 + Ftmp215*Ftmp268 + Ftmp218*Ftmp268 + Ftmp220*Ftmp268 - Ftmp221*Ftmp222 - Ftmp221*Ftmp325 - Ftmp221*Ftmp378 + Ftmp222*Ftmp270 - Ftmp224*M[75] - Ftmp225*Ftmp226*M[35] - Ftmp225*Ftmp238*M[36] - Ftmp232*Ftmp233 - Ftmp233*Ftmp237 - Ftmp239*Ftmp240 + Ftmp24*Ftmp3 - Ftmp24*Ftmp64 - Ftmp240*Ftmp241 - Ftmp243*Ftmp244 - Ftmp244*Ftmp352 - Ftmp244*Ftmp354 - Ftmp244*Ftmp355 - Ftmp245*Ftmp37*(9*Ftmp182 + Ftmp229 + 5) - Ftmp246*Ftmp247 - Ftmp25*Ftmp95 - Ftmp25*Ftmp97 + Ftmp250*Ftmp251*(Ftmp248 + Ftmp249) + Ftmp253*Ftmp254 + Ftmp254*Ftmp255 + Ftmp254*Ftmp256 + Ftmp254*Ftmp365 + Ftmp254*Ftmp367 + Ftmp259*(Ftmp18 + Ftmp258 + Ftmp82)*M[51] + Ftmp261*Ftmp262 + Ftmp263*Ftmp265 + Ftmp266*Ftmp267*M[56] + Ftmp266*Ftmp373*(Ftmp198 + Ftmp299 + Ftmp366) + Ftmp267*Ftmp269*M[57] + Ftmp268*Ftmp323 + Ftmp269*Ftmp375 + Ftmp269*Ftmp453 + Ftmp270*Ftmp271 + Ftmp270*Ftmp325 + Ftmp270*Ftmp378 - Ftmp273*(-22*Ftmp257 + Ftmp78 + Ftmp82)*M[79] + Ftmp276*x*(35*Ftmp13 - Ftmp183 + Ftmp274*Ftmp275 - 5)*M[55] + Ftmp277*Ftmp279 + Ftmp277*Ftmp282 + Ftmp277*Ftmp284 + Ftmp277*Ftmp406 + Ftmp277*Ftmp410 + Ftmp277*Ftmp420 + Ftmp277*Ftmp424 + Ftmp277*Ftmp427 + Ftmp277*Ftmp428 - Ftmp28*Ftmp30 - Ftmp28*Ftmp32 - Ftmp28*(Ftmp18 - Ftmp189*Ftmp21 + Ftmp21)*M[31] - Ftmp290*Ftmp38 + Ftmp290*Ftmp9 - Ftmp294*Ftmp38 + Ftmp294*Ftmp9 + Ftmp301*(Ftmp296 + Ftmp298 + Ftmp300)*M[22] + Ftmp301*(Ftmp300 + Ftmp307 + Ftmp308)*M[24] - Ftmp313*Ftmp320 - Ftmp326*Ftmp328 + Ftmp334*Ftmp447*(-117*Ftmp368 + Ftmp445 + Ftmp446 + Ftmp74) - Ftmp335*(Ftmp341 + Ftmp343)*M[42] - Ftmp335*(Ftmp316 + Ftmp330 + Ftmp333)*M[40] - Ftmp34*M[17] - Ftmp35*Ftmp452 - Ftmp351*(Ftmp332 + Ftmp350)*M[41] - Ftmp351*(Ftmp316 + Ftmp342 + Ftmp353)*M[43] + Ftmp362*(Ftmp356 + Ftmp361)*M[66] + Ftmp362*(Ftmp357 + Ftmp360 + Ftmp366)*M[68] + Ftmp374*(Ftmp191 + Ftmp299 + Ftmp359)*M[62] + Ftmp374*(Ftmp316 + Ftmp366 + Ftmp376)*M[64] - Ftmp38*Ftmp433 + Ftmp385*(Ftmp21 + Ftmp379 + Ftmp386*Ftmp387 + Ftmp389)*M[78] + Ftmp385*(Ftmp18*Ftmp381 + Ftmp18 + Ftmp379 + Ftmp383)*M[80] + Ftmp39*Ftmp88 + Ftmp39*Ftmp91 + Ftmp39*Ftmp94 - Ftmp39*z*M[13] + Ftmp398*(Ftmp391 + Ftmp394 + Ftmp397 + Ftmp74)*M[65] + Ftmp398*(Ftmp399 + Ftmp400 + Ftmp401 + Ftmp74)*M[69] + Ftmp398*(Ftmp411 + Ftmp413 + Ftmp415 + Ftmp417)*M[58] + Ftmp398*(Ftmp417 + Ftmp421 + Ftmp422 + Ftmp423)*M[60] - Ftmp4*M[5] + Ftmp43*Ftmp44 + Ftmp433*Ftmp9 + Ftmp447*Ftmp448 + Ftmp45*Ftmp46 + Ftmp46*(Ftmp21 + Ftmp258 + Ftmp78)*M[52] - Ftmp47*Ftmp48 - Ftmp49*M[14] + Ftmp5*Ftmp77 + Ftmp5*Ftmp81 + Ftmp5*Ftmp85 + Ftmp50*Ftmp51 - Ftmp50*Ftmp60*Ftmp62 + Ftmp55*M[27] + Ftmp57*M[26] - Ftmp58*Ftmp61 - Ftmp6*M[1] + Ftmp68*M[12] - Ftmp7*x + Ftmp72*M[14] + Ftmp73*Ftmp75*M[3] - Ftmp99*M[28] + M[0]);
#pragma omp atomic
F[1] += Ftmp0*(-Ftmp10*Ftmp103*Ftmp106 + Ftmp10*Ftmp207 + Ftmp10*Ftmp210 + Ftmp10*Ftmp212 + Ftmp10*Ftmp305 + Ftmp10*Ftmp310 + Ftmp10*Ftmp312 - Ftmp103*Ftmp47*y - Ftmp105*Ftmp120 - Ftmp105*Ftmp121 - Ftmp105*Ftmp64 + Ftmp105*z*M[5] + Ftmp108*Ftmp2 - Ftmp108*Ftmp480 - Ftmp109*Ftmp480 + Ftmp110*M[13] - Ftmp111*Ftmp478 - Ftmp112*Ftmp478 + Ftmp114*Ftmp279 + Ftmp114*Ftmp282 + Ftmp114*Ftmp284 + Ftmp114*Ftmp406 + Ftmp114*Ftmp410 + Ftmp114*Ftmp420 + Ftmp114*Ftmp424 + Ftmp114*Ftmp427 + Ftmp114*Ftmp428 - Ftmp115*Ftmp239 - Ftmp115*Ftmp241 + Ftmp116*Ftmp240 - Ftmp116*Ftmp489 + Ftmp117*Ftmp240 - Ftmp117*Ftmp489 + Ftmp118*Ftmp18*Ftmp472*M[30] - Ftmp119*(Ftmp202 + Ftmp529)*M[41] - Ftmp119*(Ftmp314 + Ftmp353 + Ftmp530)*M[52] + Ftmp124*Ftmp463 - Ftmp126*Ftmp25 + Ftmp126*Ftmp486 + Ftmp128*Ftmp439*M[27] - Ftmp128*Ftmp548 - Ftmp129*Ftmp184 - Ftmp129*Ftmp290 - Ftmp129*Ftmp294 - Ftmp129*Ftmp433 + Ftmp130*Ftmp138 + Ftmp130*Ftmp95 + Ftmp130*Ftmp97 - Ftmp131*Ftmp25 + Ftmp131*Ftmp486 + Ftmp132*Ftmp486 + Ftmp136*Ftmp26*Ftmp483 - 420*Ftmp144*Ftmp488 - Ftmp146*Ftmp468*Ftmp472*M[45] - Ftmp147*Ftmp497*Ftmp84 - 1260*Ftmp149*Ftmp170 - Ftmp151*Ftmp154*Ftmp490 - Ftmp152*Ftmp502*Ftmp503 + Ftmp156*Ftmp240 + Ftmp160*Ftmp18 + Ftmp164*Ftmp511*Ftmp517 - Ftmp167*Ftmp226*Ftmp502 + Ftmp171*Ftmp260*Ftmp511 + Ftmp175*Ftmp44*Ftmp59 + Ftmp176*Ftmp466 + Ftmp179*Ftmp251*Ftmp499 + Ftmp18*Ftmp205*M[4] - Ftmp18*Ftmp244*M[13] + Ftmp18*Ftmp25*M[7] - Ftmp18*Ftmp272*(Ftmp514 + Ftmp520)*M[68] + Ftmp18*Ftmp43 - Ftmp18*Ftmp496*M[35] - Ftmp18*Ftmp507*(Ftmp190 + 9*Ftmp386 + 5)*M[49] - Ftmp19*Ftmp475*M[15] - Ftmp193*Ftmp233 + Ftmp200*M[53] + Ftmp204*y*(Ftmp333 + 7*Ftmp386)*M[29] - Ftmp205*Ftmp88 - Ftmp205*Ftmp91 - Ftmp205*Ftmp94 - Ftmp215*Ftmp277 + Ftmp215*Ftmp510 - Ftmp218*Ftmp277 + Ftmp218*Ftmp510 - Ftmp220*Ftmp277 + Ftmp220*Ftmp510 - Ftmp221*Ftmp253 - Ftmp221*Ftmp255 - Ftmp221*Ftmp365 - Ftmp221*Ftmp367 + Ftmp222*Ftmp254 - Ftmp224*M[81] + Ftmp23*Ftmp457 - Ftmp232*Ftmp506 - Ftmp233*Ftmp247 - Ftmp236*Ftmp485 + Ftmp236*Ftmp9 - Ftmp237*Ftmp506 - Ftmp243*Ftmp505 + Ftmp250*Ftmp375 + Ftmp250*Ftmp453 + Ftmp250*Ftmp456*(Ftmp291*Ftmp67 + Ftmp425 + Ftmp450 - Ftmp454*Ftmp549 - 126*Ftmp69 + 28) + 1260*Ftmp250*(Ftmp1*Ftmp388*Ftmp80 + Ftmp188 - 6*Ftmp549*Ftmp65 - 34*Ftmp65 + 9)*M[71] + Ftmp253*Ftmp518 + Ftmp254*Ftmp271 + Ftmp254*Ftmp325 + Ftmp254*Ftmp378 + Ftmp255*Ftmp518 + Ftmp256*Ftmp518 + Ftmp259*(Ftmp508 + Ftmp509)*M[39] + 8505*Ftmp260*Ftmp528 + 945*Ftmp260*Ftmp548 + Ftmp261*(Ftmp249 + Ftmp532)*M[61] + Ftmp265*Ftmp516 + 1890*Ftmp266*Ftmp512 - Ftmp273*(Ftmp513 + Ftmp520 + Ftmp82)*M[64] + Ftmp276*y*(-Ftmp235 + Ftmp275*Ftmp281 + 35*Ftmp65 - 5)*M[76] - Ftmp277*Ftmp323 - Ftmp301*Ftmp473*M[12] - Ftmp32*Ftmp462 + Ftmp323*Ftmp510 - Ftmp326*Ftmp369 - Ftmp33*(48*Ftmp368 - 5*Ftmp549*Ftmp93 + Ftmp553 + Ftmp67*Ftmp92 - 28*Ftmp69)*M[46] + Ftmp334*Ftmp552*(-Ftmp22*Ftmp550 + Ftmp445 + Ftmp472 + Ftmp551) - Ftmp335*(Ftmp203 + Ftmp314 + Ftmp330)*M[37] - Ftmp339*Ftmp485 + Ftmp339*Ftmp9 - Ftmp34*M[14] - Ftmp347*Ftmp485 + Ftmp347*Ftmp9 - Ftmp349*Ftmp485 + Ftmp349*Ftmp9 + Ftmp351*Ftmp473*M[26] - Ftmp352*Ftmp505 - Ftmp354*Ftmp505 - Ftmp355*Ftmp505 + Ftmp362*(Ftmp532 + Ftmp533)*M[62] + Ftmp365*Ftmp518 + Ftmp367*Ftmp518 + 8505*Ftmp369*Ftmp517 - Ftmp4*M[7] - Ftmp40*Ftmp461*Ftmp487 + Ftmp44*Ftmp52 + Ftmp448*Ftmp552 + Ftmp457*Ftmp458 - Ftmp459*Ftmp461*M[9] + Ftmp460*Ftmp56*M[20] - Ftmp462*(-Ftmp21*Ftmp228 + Ftmp501)*M[24] - Ftmp464*Ftmp47 - Ftmp464*Ftmp479 - Ftmp465*M[17] + Ftmp467*M[27] - Ftmp468*Ftmp60*(Ftmp509 + Ftmp520)*M[63] - Ftmp468*Ftmp61*M[47] - Ftmp469*Ftmp470*z + Ftmp469*Ftmp498 + Ftmp471*M[10] + Ftmp472*Ftmp481*x*M[25] + Ftmp473*Ftmp73*M[6] + Ftmp474*Ftmp77 + Ftmp474*Ftmp81 + Ftmp474*Ftmp85 - Ftmp475*Ftmp477*M[16] + Ftmp476*Ftmp524*(Ftmp298 + Ftmp521 + Ftmp523) + Ftmp476*Ftmp526*(Ftmp307 + Ftmp523 + Ftmp525) - Ftmp481*(Ftmp529 + Ftmp531)*M[40] - Ftmp481*(Ftmp357 + Ftmp523 + Ftmp530)*M[51] + Ftmp482*Ftmp88 + Ftmp482*Ftmp91 + Ftmp483*Ftmp484 - Ftmp491*M[36] + Ftmp494*Ftmp495 - Ftmp5*M[4] + Ftmp500*M[35] + Ftmp515*(Ftmp508 + Ftmp514)*M[43] - Ftmp527*Ftmp528 + Ftmp534*(Ftmp230 + Ftmp359 + Ftmp522)*M[66] + Ftmp534*(Ftmp314 + Ftmp376 + Ftmp535)*M[79] + Ftmp538*(Ftmp12*Ftmp381 + Ftmp12 + Ftmp383 + Ftmp536)*M[69] + Ftmp538*(Ftmp182*Ftmp387 + Ftmp21 + Ftmp536 + Ftmp539)*M[60] + Ftmp541*(Ftmp391 + Ftmp413 + Ftmp472 + Ftmp540)*M[58] + Ftmp541*(Ftmp394 + Ftmp415 + Ftmp544 + Ftmp545)*M[65] + Ftmp541*(Ftmp401 + Ftmp472 + Ftmp542 + Ftmp543)*M[80] + Ftmp541*(Ftmp423 + Ftmp545 + Ftmp546 + Ftmp547)*M[78] + Ftmp56*x*(Ftmp143*Ftmp90 + 64*Ftmp368 - Ftmp441*Ftmp549 + Ftmp553 - 36*Ftmp69)*M[72] + Ftmp56*(Ftmp501 + Ftmp508)*M[42] + Ftmp57*M[23] - Ftmp6*M[0] - Ftmp7*y + Ftmp72*M[17] - Ftmp79*M[1] - Ftmp99*M[32] + M[1]);
#pragma omp atomic
F[2] += Ftmp0*(-Ftmp10*Ftmp126 - Ftmp10*Ftmp131 - Ftmp10*Ftmp132 - Ftmp102*Ftmp47*z - Ftmp107*Ftmp478 - Ftmp108*Ftmp478 - Ftmp109*Ftmp478 + Ftmp11*x - Ftmp110*Ftmp120 - Ftmp110*Ftmp121 + Ftmp110*Ftmp3 - Ftmp110*Ftmp64 - Ftmp111*Ftmp560 + Ftmp112*Ftmp2 - Ftmp112*Ftmp560 - Ftmp114*Ftmp253 - Ftmp114*Ftmp255 - Ftmp114*Ftmp256 - Ftmp114*Ftmp365 - Ftmp114*Ftmp367 - Ftmp115*Ftmp232 - Ftmp115*Ftmp237 + Ftmp116*Ftmp233 - Ftmp116*Ftmp568 + Ftmp117*Ftmp233 + Ftmp118*Ftmp44*Ftmp557*M[27] + Ftmp118*(-Ftmp19*Ftmp550 + Ftmp444 + Ftmp446 + Ftmp551 + Ftmp557)*M[67] - Ftmp119*(Ftmp202 + Ftmp583)*M[42] - Ftmp119*(Ftmp327 + Ftmp333 + Ftmp530)*M[51] + Ftmp126*Ftmp565 - Ftmp129*Ftmp21*M[13] + Ftmp130*Ftmp88 + Ftmp130*Ftmp91 + Ftmp130*Ftmp94 + Ftmp131*Ftmp565 + Ftmp132*Ftmp565 - Ftmp133*Ftmp452 + 60*Ftmp135*Ftmp35 - Ftmp138*Ftmp205 + Ftmp138*Ftmp563 - 420*Ftmp142*Ftmp145 - Ftmp145*Ftmp460*Ftmp487 - Ftmp147*Ftmp503*Ftmp570 + Ftmp15*Ftmp21 - 2835*Ftmp150*M[47] - Ftmp152*Ftmp497*Ftmp80 - Ftmp153*Ftmp168 + Ftmp153*Ftmp372 + Ftmp153*Ftmp438 + Ftmp156*Ftmp233 - Ftmp156*Ftmp568 + Ftmp159*Ftmp177*Ftmp499*M[59] + Ftmp159*Ftmp180*Ftmp466 - Ftmp16*Ftmp319*Ftmp60*(Ftmp513 + Ftmp579 + Ftmp78) + Ftmp161*Ftmp172 + Ftmp161*Ftmp51 + Ftmp165*Ftmp21*Ftmp29 + Ftmp174*Ftmp263*Ftmp575 - Ftmp184*Ftmp244 + Ftmp185*Ftmp241 - Ftmp193*Ftmp240 + Ftmp20*Ftmp554 + Ftmp204*z*(Ftmp353 + 7*Ftmp380)*M[33] + Ftmp205*Ftmp21*M[5] - Ftmp205*Ftmp95 - Ftmp205*Ftmp97 + Ftmp207*Ftmp25 + Ftmp21*Ftmp371*Ftmp40 + Ftmp21*Ftmp460*Ftmp494 - Ftmp21*Ftmp461*Ftmp470 - Ftmp21*Ftmp479*M[11] - Ftmp21*Ftmp496*M[36] - Ftmp21*Ftmp507*(Ftmp197 + 9*Ftmp380 + 5)*M[54] + Ftmp210*Ftmp25 + Ftmp212*Ftmp25 + Ftmp215*Ftmp254 + Ftmp218*Ftmp254 - Ftmp22*Ftmp559*M[18] + Ftmp220*Ftmp254 + Ftmp221*Ftmp279 + Ftmp221*Ftmp282 + Ftmp221*Ftmp284 + Ftmp221*Ftmp406 + Ftmp221*Ftmp410 + Ftmp221*Ftmp420 + Ftmp221*Ftmp424 + Ftmp221*Ftmp427 + Ftmp221*Ftmp428 - Ftmp222*Ftmp277 + Ftmp222*Ftmp567 - Ftmp236*Ftmp505 - 1260*Ftmp238*Ftmp570*M[48] - Ftmp239*Ftmp571 - Ftmp240*Ftmp247 - Ftmp241*Ftmp571 - Ftmp243*Ftmp564 + Ftmp243*Ftmp9 - Ftmp244*Ftmp290 - Ftmp244*Ftmp294 - Ftmp244*Ftmp433 + Ftmp25*Ftmp305 + Ftmp25*Ftmp310 + Ftmp25*Ftmp312 + Ftmp250*Ftmp373*(Ftmp533 + Ftmp584) + 3780*Ftmp250*(Ftmp376 + Ftmp573)*M[74] + Ftmp253*Ftmp577 + Ftmp254*Ftmp323 + Ftmp255*Ftmp577 + Ftmp256*Ftmp577 + Ftmp262*Ftmp534 + Ftmp263*Ftmp375 + Ftmp263*Ftmp453 + Ftmp263*Ftmp456*(Ftmp450 + Ftmp454*Ftmp65 - Ftmp455*Ftmp590 - 126*Ftmp590*Ftmp65) + Ftmp264*Ftmp35*(-Ftmp242 + Ftmp275*Ftmp283 + 35*Ftmp69 - 5) - Ftmp271*Ftmp277 + Ftmp271*Ftmp567 - Ftmp277*Ftmp325 - Ftmp277*Ftmp378 - Ftmp29*Ftmp490*Ftmp556 - Ftmp30*Ftmp555 - Ftmp301*Ftmp558*M[14] - Ftmp313*Ftmp370 + Ftmp325*Ftmp567 - Ftmp328*Ftmp527 - Ftmp339*Ftmp505 - Ftmp347*Ftmp505 - Ftmp349*Ftmp505 - Ftmp35*Ftmp524*(-Ftmp18*Ftmp228 + Ftmp569) - Ftmp351*(Ftmp203 + Ftmp327 + Ftmp342)*M[39] - Ftmp352*Ftmp564 + Ftmp352*Ftmp9 - Ftmp354*Ftmp564 + Ftmp354*Ftmp9 - Ftmp355*Ftmp564 + Ftmp355*Ftmp9 - Ftmp36*Ftmp48 - Ftmp36*(40*Ftmp368 - Ftmp396 - Ftmp440*Ftmp590 + Ftmp591 + Ftmp66*Ftmp93)*M[46] + Ftmp365*Ftmp577 + Ftmp367*Ftmp577 + Ftmp374*Ftmp512 + Ftmp378*Ftmp567 + Ftmp4*Ftmp77 + Ftmp4*Ftmp81 + Ftmp4*Ftmp85 - Ftmp4*x*M[0] - Ftmp4*y*M[1] + Ftmp436*Ftmp46*(56*Ftmp368 + Ftmp434*Ftmp89 - Ftmp442*Ftmp590 + Ftmp591 - 28*Ftmp65) + Ftmp448*z*M[67] + Ftmp45*Ftmp53 + Ftmp458*Ftmp554 + Ftmp46*(Ftmp572 + Ftmp574)*M[37] + Ftmp460*Ftmp53*M[21] - Ftmp460*Ftmp555*M[9] + Ftmp463*z*M[20] - Ftmp465*M[16] + Ftmp467*M[26] + Ftmp471*M[11] - Ftmp474*M[7] - Ftmp477*Ftmp559*M[17] + 60*Ftmp484*z - Ftmp49*M[11] - Ftmp491*M[35] + Ftmp495*Ftmp498 - Ftmp5*M[5] + Ftmp500*M[36] + Ftmp515*(Ftmp572 + Ftmp576)*M[40] + 3780*Ftmp516*Ftmp575*M[81] + Ftmp526*Ftmp582*(Ftmp296 + Ftmp525 + Ftmp581) + Ftmp53*(Ftmp569 + Ftmp572)*M[41] + Ftmp55*M[23] - Ftmp556*x*(Ftmp574 + Ftmp579)*M[62] - Ftmp556*y*(Ftmp576 + Ftmp579)*M[66] - Ftmp557*Ftmp566*y*M[47] + Ftmp558*Ftmp73*M[8] - Ftmp561*(Ftmp531 + Ftmp583)*M[43] - Ftmp561*(Ftmp356 + Ftmp530 + Ftmp581)*M[52] + Ftmp562*x*M[28] + Ftmp562*y*M[32] + Ftmp563*Ftmp95 + Ftmp563*Ftmp97 + Ftmp566*(Ftmp249 + Ftmp584)*M[64] + 30*Ftmp582*(Ftmp308 + Ftmp521 + Ftmp581)*M[24] + Ftmp585*(Ftmp230 + Ftmp366 + Ftmp580)*M[68] + Ftmp585*(Ftmp190 + Ftmp327 + Ftmp535 + 3)*M[79] + Ftmp587*(99*Ftmp12*Ftmp386 + Ftmp12 + Ftmp389 + Ftmp586)*M[65] + Ftmp587*(99*Ftmp18*Ftmp182 + Ftmp18 + Ftmp539 + Ftmp586)*M[58] + Ftmp588*(Ftmp397 + Ftmp542 + Ftmp547 + Ftmp557)*M[78] + Ftmp588*(Ftmp399 + Ftmp422 + Ftmp540 + Ftmp557)*M[60] + Ftmp588*(Ftmp400 + Ftmp421 + Ftmp544 + Ftmp589)*M[69] + Ftmp588*(Ftmp411 + Ftmp543 + Ftmp546 + Ftmp589)*M[80] + Ftmp68*M[16] - Ftmp83*M[2] + M[2]);

}

void S2Mc_6(double x, double y, double z, double * S, double * M) {
double Mtmp0 = x*S[0];
double Mtmp1 = z*S[2];
double Mtmp2 = -Mtmp1;
double Mtmp3 = x*S[1];
double Mtmp4 = y*S[0];
double Mtmp5 = x*S[2];
double Mtmp6 = z*S[0];
double Mtmp7 = y*S[1];
double Mtmp8 = y*S[2];
double Mtmp9 = z*S[1];
double Mtmp10 = Mtmp1*x;
double Mtmp11 = pow(x, 2);
double Mtmp12 = pow(z, 2);
double Mtmp13 = (1.0/2.0)*S[0];
double Mtmp14 = Mtmp12*Mtmp13;
double Mtmp15 = Mtmp0*y;
double Mtmp16 = Mtmp1*y;
double Mtmp17 = (1.0/2.0)*S[1];
double Mtmp18 = Mtmp12*Mtmp17;
double Mtmp19 = Mtmp0*z;
double Mtmp20 = (1.0/2.0)*S[2];
double Mtmp21 = -Mtmp12*Mtmp20;
double Mtmp22 = Mtmp3*y;
double Mtmp23 = pow(y, 2);
double Mtmp24 = Mtmp5*y;
double Mtmp25 = Mtmp3*z;
double Mtmp26 = Mtmp4*z;
double Mtmp27 = Mtmp7*z;
double Mtmp28 = pow(x, 3);
double Mtmp29 = Mtmp28*S[0];
double Mtmp30 = pow(z, 3);
double Mtmp31 = Mtmp30*S[2];
double Mtmp32 = 3*Mtmp12;
double Mtmp33 = Mtmp0*Mtmp32;
double Mtmp34 = 3*Mtmp11;
double Mtmp35 = Mtmp1*Mtmp34;
double Mtmp36 = (1.0/2.0)*Mtmp12;
double Mtmp37 = Mtmp10*y + (1.0/2.0)*Mtmp12*Mtmp4 + Mtmp3*Mtmp36;
double Mtmp38 = Mtmp28*S[2];
double Mtmp39 = Mtmp30*S[0];
double Mtmp40 = Mtmp32*Mtmp7;
double Mtmp41 = 3*Mtmp23;
double Mtmp42 = Mtmp1*Mtmp41;
double Mtmp43 = Mtmp30*S[1];
double Mtmp44 = (1.0/2.0)*Mtmp11;
double Mtmp45 = pow(y, 3);
double Mtmp46 = (1.0/2.0)*Mtmp23;
double Mtmp47 = Mtmp45*S[1];
double Mtmp48 = Mtmp45*S[2];
double Mtmp49 = pow(x, 4);
double Mtmp50 = Mtmp11*S[0];
double Mtmp51 = 6*Mtmp12;
double Mtmp52 = 4*Mtmp28;
double Mtmp53 = pow(z, 4);
double Mtmp54 = 4*Mtmp30;
double Mtmp55 = Mtmp5*Mtmp54 + Mtmp53*S[0];
double Mtmp56 = Mtmp11*Mtmp51;
double Mtmp57 = 12*Mtmp12;
double Mtmp58 = Mtmp11*Mtmp16;
double Mtmp59 = Mtmp53*S[1] + Mtmp54*Mtmp8;
double Mtmp60 = Mtmp53*S[2];
double Mtmp61 = 2*Mtmp28;
double Mtmp62 = 6*Mtmp23;
double Mtmp63 = Mtmp23*S[0];
double Mtmp64 = -Mtmp24*Mtmp32 - Mtmp3*Mtmp30 - Mtmp30*Mtmp4;
double Mtmp65 = 2*Mtmp45;
double Mtmp66 = Mtmp32*S[1];
double Mtmp67 = 2*Mtmp30;
double Mtmp68 = Mtmp32*S[2];
double Mtmp69 = pow(y, 4);
double Mtmp70 = 4*Mtmp45;
double Mtmp71 = Mtmp23*Mtmp51;
double Mtmp72 = pow(x, 5);
double Mtmp73 = pow(z, 5);
double Mtmp74 = Mtmp73*S[2];
double Mtmp75 = -Mtmp74;
double Mtmp76 = 10*Mtmp12;
double Mtmp77 = 5*Mtmp49;
double Mtmp78 = -Mtmp1*Mtmp77 - Mtmp29*Mtmp76;
double Mtmp79 = 5*Mtmp53;
double Mtmp80 = 10*Mtmp11;
double Mtmp81 = Mtmp0*Mtmp79 + Mtmp31*Mtmp80;
double Mtmp82 = Mtmp12*Mtmp28*S[1];
double Mtmp83 = 30*Mtmp12;
double Mtmp84 = Mtmp11*Mtmp4;
double Mtmp85 = 20*Mtmp28;
double Mtmp86 = 20*Mtmp30;
double Mtmp87 = Mtmp24*Mtmp86 + 5*Mtmp3*Mtmp53 + 5*Mtmp4*Mtmp53;
double Mtmp88 = Mtmp5*Mtmp79 + Mtmp73*S[0];
double Mtmp89 = 10*Mtmp53;
double Mtmp90 = 10*Mtmp23;
double Mtmp91 = 20*Mtmp31;
double Mtmp92 = Mtmp31*Mtmp90 + Mtmp7*Mtmp79;
double Mtmp93 = Mtmp23*Mtmp83;
double Mtmp94 = Mtmp11*Mtmp83;
double Mtmp95 = Mtmp11*Mtmp23;
double Mtmp96 = -Mtmp0*Mtmp93 - 30*Mtmp1*Mtmp95 - Mtmp7*Mtmp94 - 3*Mtmp74;
double Mtmp97 = Mtmp73*S[1] + Mtmp79*Mtmp8;
double Mtmp98 = Mtmp12*Mtmp45*S[0];
double Mtmp99 = 5*Mtmp39;
double Mtmp100 = 5*Mtmp38;
double Mtmp101 = 15*Mtmp12;
double Mtmp102 = 10*Mtmp30;
double Mtmp103 = 15*Mtmp95;
double Mtmp104 = 5*Mtmp69;
double Mtmp105 = -Mtmp1*Mtmp104 - Mtmp47*Mtmp76;
double Mtmp106 = 5*Mtmp43;
double Mtmp107 = 5*Mtmp48;
double Mtmp108 = pow(y, 5);
double Mtmp109 = 20*Mtmp45;
#pragma omp atomic
M[0] += S[0];
#pragma omp atomic
M[1] += S[1];
#pragma omp atomic
M[2] += S[2];
#pragma omp atomic
M[3] += Mtmp0 + Mtmp2;
#pragma omp atomic
M[4] += Mtmp3 + Mtmp4;
#pragma omp atomic
M[5] += Mtmp5 + Mtmp6;
#pragma omp atomic
M[6] += Mtmp2 + Mtmp7;
#pragma omp atomic
M[7] += Mtmp8 + Mtmp9;
#pragma omp atomic
M[8] += -Mtmp10 + (1.0/2.0)*Mtmp11*S[0] - Mtmp14;
#pragma omp atomic
M[9] += Mtmp11*Mtmp17 + Mtmp15 - Mtmp16 - Mtmp18;
#pragma omp atomic
M[10] += Mtmp11*Mtmp20 + Mtmp19 + Mtmp21;
#pragma omp atomic
M[11] += -Mtmp10 + Mtmp13*Mtmp23 - Mtmp14 + Mtmp22;
#pragma omp atomic
M[12] += Mtmp24 + Mtmp25 + Mtmp26;
#pragma omp atomic
M[13] += -Mtmp16 - Mtmp18 + (1.0/2.0)*Mtmp23*S[1];
#pragma omp atomic
M[14] += Mtmp20*Mtmp23 + Mtmp21 + Mtmp27;
#pragma omp atomic
M[15] += (1.0/6.0)*Mtmp29 + (1.0/6.0)*Mtmp31 - 1.0/6.0*Mtmp33 - 1.0/6.0*Mtmp35;
#pragma omp atomic
M[16] += (1.0/2.0)*Mtmp11*y*S[0] + (1.0/6.0)*Mtmp28*S[1] - Mtmp37;
#pragma omp atomic
M[17] += -1.0/6.0*Mtmp32*Mtmp5 + (1.0/6.0)*Mtmp34*Mtmp6 + (1.0/6.0)*Mtmp38 - 1.0/6.0*Mtmp39;
#pragma omp atomic
M[18] += (1.0/2.0)*Mtmp11*y*S[1] + (1.0/2.0)*Mtmp23*x*S[0] + (1.0/3.0)*Mtmp30*S[2] - 1.0/6.0*Mtmp33 - 1.0/6.0*Mtmp35 - 1.0/6.0*Mtmp40 - 1.0/6.0*Mtmp42;
#pragma omp atomic
M[19] += Mtmp15*z - Mtmp36*Mtmp8 - 1.0/6.0*Mtmp43 + Mtmp44*Mtmp8 + Mtmp44*Mtmp9;
#pragma omp atomic
M[20] += (1.0/2.0)*Mtmp23*x*S[1] - Mtmp37 + (1.0/6.0)*Mtmp45*S[0];
#pragma omp atomic
M[21] += Mtmp22*z - Mtmp36*Mtmp5 - 1.0/6.0*Mtmp39 + Mtmp46*Mtmp5 + Mtmp46*Mtmp6;
#pragma omp atomic
M[22] += (1.0/6.0)*Mtmp31 - 1.0/6.0*Mtmp40 - 1.0/6.0*Mtmp42 + (1.0/6.0)*Mtmp47;
#pragma omp atomic
M[23] += -1.0/6.0*Mtmp32*Mtmp8 + (1.0/6.0)*Mtmp41*Mtmp9 - 1.0/6.0*Mtmp43 + (1.0/6.0)*Mtmp48;
#pragma omp atomic
M[24] += -1.0/24.0*Mtmp1*Mtmp52 + (1.0/24.0)*Mtmp49*S[0] - 1.0/24.0*Mtmp50*Mtmp51 + (1.0/24.0)*Mtmp55;
#pragma omp atomic
M[25] += -1.0/24.0*Mtmp15*Mtmp57 + (1.0/24.0)*Mtmp4*Mtmp52 + (1.0/24.0)*Mtmp49*S[1] - 1.0/24.0*Mtmp56*S[1] - 1.0/2.0*Mtmp58 + (1.0/24.0)*Mtmp59;
#pragma omp atomic
M[26] += -1.0/24.0*Mtmp0*Mtmp54 + (1.0/24.0)*Mtmp49*S[2] + (1.0/24.0)*Mtmp52*Mtmp6 - 1.0/24.0*Mtmp56*S[2] + (1.0/24.0)*Mtmp60;
#pragma omp atomic
M[27] += -1.0/12.0*Mtmp1*Mtmp61 - 1.0/12.0*Mtmp10*Mtmp62 + (1.0/4.0)*Mtmp11*Mtmp23*S[0] - 1.0/12.0*Mtmp22*Mtmp51 + (1.0/6.0)*Mtmp28*y*S[1] + (1.0/3.0)*Mtmp30*x*S[2] - 1.0/12.0*Mtmp32*Mtmp50 - 1.0/12.0*Mtmp32*Mtmp63 + (1.0/12.0)*Mtmp53*S[0];
#pragma omp atomic
M[28] += (1.0/6.0)*Mtmp26*Mtmp34 + (1.0/6.0)*Mtmp28*Mtmp8 + (1.0/6.0)*Mtmp28*Mtmp9 + (1.0/6.0)*Mtmp64;
#pragma omp atomic
M[29] += -1.0/12.0*Mtmp1*Mtmp65 + (1.0/4.0)*Mtmp11*Mtmp23*S[1] - 1.0/12.0*Mtmp11*Mtmp66 - 1.0/12.0*Mtmp15*Mtmp51 - 1.0/12.0*Mtmp23*Mtmp66 + (1.0/3.0)*Mtmp30*y*S[2] + (1.0/6.0)*Mtmp45*x*S[0] + (1.0/12.0)*Mtmp53*S[1] - 1.0/2.0*Mtmp58;
#pragma omp atomic
M[30] += -1.0/12.0*Mtmp0*Mtmp67 + (1.0/2.0)*Mtmp11*Mtmp27 - 1.0/12.0*Mtmp11*Mtmp68 + (1.0/12.0)*Mtmp19*Mtmp62 + (1.0/12.0)*Mtmp23*Mtmp34*S[2] - 1.0/12.0*Mtmp23*Mtmp68 + (1.0/12.0)*Mtmp60 - 1.0/12.0*Mtmp67*Mtmp7;
#pragma omp atomic
M[31] += -1.0/2.0*Mtmp10*Mtmp23 - 1.0/24.0*Mtmp22*Mtmp57 + (1.0/24.0)*Mtmp3*Mtmp70 - 1.0/24.0*Mtmp51*Mtmp63 + (1.0/24.0)*Mtmp55 + (1.0/24.0)*Mtmp69*S[0];
#pragma omp atomic
M[32] += (1.0/6.0)*Mtmp25*Mtmp41 + (1.0/6.0)*Mtmp45*Mtmp5 + (1.0/6.0)*Mtmp45*Mtmp6 + (1.0/6.0)*Mtmp64;
#pragma omp atomic
M[33] += -1.0/24.0*Mtmp1*Mtmp70 + (1.0/24.0)*Mtmp59 + (1.0/24.0)*Mtmp69*S[1] - 1.0/24.0*Mtmp71*S[1];
#pragma omp atomic
M[34] += -1.0/24.0*Mtmp54*Mtmp7 + (1.0/24.0)*Mtmp60 + (1.0/24.0)*Mtmp69*S[2] + (1.0/24.0)*Mtmp70*Mtmp9 - 1.0/24.0*Mtmp71*S[2];
#pragma omp atomic
M[35] += (1.0/120.0)*Mtmp72*S[0] + (1.0/120.0)*Mtmp75 + (1.0/120.0)*Mtmp78 + (1.0/120.0)*Mtmp81;
#pragma omp atomic
M[36] += -1.0/120.0*Mtmp16*Mtmp85 + (1.0/120.0)*Mtmp4*Mtmp77 + (1.0/120.0)*Mtmp72*S[1] - 1.0/12.0*Mtmp82 - 1.0/120.0*Mtmp83*Mtmp84 + (1.0/120.0)*Mtmp87;
#pragma omp atomic
M[37] += -1.0/120.0*Mtmp38*Mtmp76 - 1.0/120.0*Mtmp39*Mtmp80 + (1.0/120.0)*Mtmp6*Mtmp77 + (1.0/120.0)*Mtmp72*S[2] + (1.0/120.0)*Mtmp88;
#pragma omp atomic
M[38] += (1.0/120.0)*Mtmp0*Mtmp89 + (1.0/120.0)*Mtmp11*Mtmp91 + (1.0/120.0)*Mtmp29*Mtmp90 + (1.0/120.0)*Mtmp7*Mtmp77 + (1.0/120.0)*Mtmp78 + (1.0/120.0)*Mtmp92 + (1.0/120.0)*Mtmp96;
#pragma omp atomic
M[39] += -1.0/120.0*Mtmp15*Mtmp86 + (1.0/120.0)*Mtmp26*Mtmp85 - 1.0/120.0*Mtmp43*Mtmp80 + (1.0/120.0)*Mtmp77*Mtmp8 + (1.0/120.0)*Mtmp77*Mtmp9 - 1.0/120.0*Mtmp8*Mtmp94 + (1.0/120.0)*Mtmp97;
#pragma omp atomic
M[40] += -1.0/12.0*Mtmp10*Mtmp65 + (1.0/12.0)*Mtmp11*Mtmp45*S[0] - 1.0/12.0*Mtmp16*Mtmp61 + (1.0/12.0)*Mtmp23*Mtmp28*S[1] - 1.0/12.0*Mtmp23*Mtmp3*Mtmp32 + (1.0/3.0)*Mtmp30*x*y*S[2] - 1.0/12.0*Mtmp32*Mtmp84 + (1.0/12.0)*Mtmp53*x*S[1] + (1.0/12.0)*Mtmp53*y*S[0] - 1.0/12.0*Mtmp82 - 1.0/12.0*Mtmp98;
#pragma omp atomic
M[41] += -1.0/60.0*Mtmp100*Mtmp12 + (1.0/60.0)*Mtmp100*Mtmp23 - 1.0/60.0*Mtmp101*Mtmp23*Mtmp5 - 1.0/60.0*Mtmp102*Mtmp22 + (1.0/60.0)*Mtmp103*Mtmp6 - 1.0/60.0*Mtmp11*Mtmp99 - 1.0/60.0*Mtmp23*Mtmp99 + (1.0/6.0)*Mtmp27*Mtmp28 + (1.0/60.0)*Mtmp88;
#pragma omp atomic
M[42] += (1.0/120.0)*Mtmp0*Mtmp104 + (1.0/120.0)*Mtmp105 + (1.0/120.0)*Mtmp23*Mtmp91 + (1.0/120.0)*Mtmp47*Mtmp80 + (1.0/120.0)*Mtmp7*Mtmp89 + (1.0/120.0)*Mtmp81 + (1.0/120.0)*Mtmp96;
#pragma omp atomic
M[43] += -1.0/60.0*Mtmp101*Mtmp11*Mtmp8 - 1.0/60.0*Mtmp102*Mtmp15 + (1.0/60.0)*Mtmp103*Mtmp9 - 1.0/60.0*Mtmp106*Mtmp11 - 1.0/60.0*Mtmp106*Mtmp23 + (1.0/60.0)*Mtmp107*Mtmp11 - 1.0/60.0*Mtmp107*Mtmp12 + (1.0/6.0)*Mtmp19*Mtmp45 + (1.0/60.0)*Mtmp97;
#pragma omp atomic
M[44] += -1.0/120.0*Mtmp10*Mtmp109 + (1.0/120.0)*Mtmp104*Mtmp3 + (1.0/120.0)*Mtmp108*S[0] - 1.0/120.0*Mtmp3*Mtmp93 + (1.0/120.0)*Mtmp87 - 1.0/12.0*Mtmp98;
#pragma omp atomic
M[45] += (1.0/120.0)*Mtmp104*Mtmp5 + (1.0/120.0)*Mtmp104*Mtmp6 + (1.0/120.0)*Mtmp109*Mtmp25 - 1.0/120.0*Mtmp22*Mtmp86 - 1.0/120.0*Mtmp39*Mtmp90 - 1.0/120.0*Mtmp5*Mtmp93 + (1.0/120.0)*Mtmp88;
#pragma omp atomic
M[46] += (1.0/120.0)*Mtmp105 + (1.0/120.0)*Mtmp108*S[1] + (1.0/120.0)*Mtmp75 + (1.0/120.0)*Mtmp92;
#pragma omp atomic
M[47] += (1.0/120.0)*Mtmp104*Mtmp9 + (1.0/120.0)*Mtmp108*S[2] - 1.0/120.0*Mtmp43*Mtmp90 - 1.0/120.0*Mtmp48*Mtmp76 + (1.0/120.0)*Mtmp97;

}

void M2Mc_6(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = z*M[2];
double Mstmp2 = -Mstmp1;
double Mstmp3 = x*M[1];
double Mstmp4 = y*M[0];
double Mstmp5 = x*M[2];
double Mstmp6 = z*M[0];
double Mstmp7 = y*M[1];
double Mstmp8 = y*M[2];
double Mstmp9 = z*M[1];
double Mstmp10 = x*M[3];
double Mstmp11 = pow(x, 2);
double Mstmp12 = (1.0/2.0)*M[0];
double Mstmp13 = z*M[5];
double Mstmp14 = pow(z, 2);
double Mstmp15 = Mstmp1*x;
double Mstmp16 = -Mstmp12*Mstmp14 - Mstmp13 - Mstmp15;
double Mstmp17 = x*M[4];
double Mstmp18 = y*M[3];
double Mstmp19 = (1.0/2.0)*M[1];
double Mstmp20 = Mstmp0*y;
double Mstmp21 = z*M[7];
double Mstmp22 = Mstmp1*y;
double Mstmp23 = -Mstmp14*Mstmp19 - Mstmp21 - Mstmp22;
double Mstmp24 = x*M[5];
double Mstmp25 = z*M[3];
double Mstmp26 = Mstmp0*z;
double Mstmp27 = (1.0/2.0)*M[2];
double Mstmp28 = -Mstmp14*Mstmp27;
double Mstmp29 = x*M[6];
double Mstmp30 = y*M[4];
double Mstmp31 = pow(y, 2);
double Mstmp32 = Mstmp3*y;
double Mstmp33 = x*M[7];
double Mstmp34 = y*M[5];
double Mstmp35 = z*M[4];
double Mstmp36 = Mstmp5*y;
double Mstmp37 = Mstmp3*z;
double Mstmp38 = Mstmp4*z;
double Mstmp39 = y*M[6];
double Mstmp40 = y*M[7];
double Mstmp41 = z*M[6];
double Mstmp42 = Mstmp7*z;
double Mstmp43 = x*M[8];
double Mstmp44 = z*M[10];
double Mstmp45 = Mstmp13*x;
double Mstmp46 = (1.0/2.0)*M[3];
double Mstmp47 = pow(x, 3);
double Mstmp48 = (1.0/6.0)*Mstmp47;
double Mstmp49 = Mstmp14*Mstmp46;
double Mstmp50 = pow(z, 3);
double Mstmp51 = (1.0/6.0)*Mstmp50;
double Mstmp52 = Mstmp51*M[2];
double Mstmp53 = (1.0/2.0)*Mstmp14;
double Mstmp54 = Mstmp0*Mstmp53;
double Mstmp55 = (1.0/2.0)*Mstmp11;
double Mstmp56 = Mstmp1*Mstmp55;
double Mstmp57 = x*M[9];
double Mstmp58 = y*M[8];
double Mstmp59 = (1.0/2.0)*M[4];
double Mstmp60 = Mstmp10*y;
double Mstmp61 = z*M[12];
double Mstmp62 = Mstmp21*x;
double Mstmp63 = Mstmp13*y;
double Mstmp64 = -Mstmp14*Mstmp59 - Mstmp15*y - Mstmp3*Mstmp53 - Mstmp4*Mstmp53 - Mstmp61 - Mstmp62 - Mstmp63;
double Mstmp65 = x*M[10];
double Mstmp66 = z*M[8];
double Mstmp67 = (1.0/2.0)*M[5];
double Mstmp68 = Mstmp10*z;
double Mstmp69 = -Mstmp14*Mstmp67 - Mstmp5*Mstmp53 - Mstmp51*M[0];
double Mstmp70 = z*M[14];
double Mstmp71 = Mstmp21*y;
double Mstmp72 = (1.0/2.0)*M[6];
double Mstmp73 = Mstmp14*Mstmp72;
double Mstmp74 = Mstmp53*Mstmp7;
double Mstmp75 = (1.0/2.0)*Mstmp31;
double Mstmp76 = Mstmp1*Mstmp75;
double Mstmp77 = x*M[12];
double Mstmp78 = y*M[10];
double Mstmp79 = z*M[9];
double Mstmp80 = (1.0/2.0)*M[7];
double Mstmp81 = Mstmp24*y;
double Mstmp82 = Mstmp17*z;
double Mstmp83 = Mstmp18*z;
double Mstmp84 = -Mstmp14*Mstmp80 - Mstmp51*M[1] - Mstmp53*Mstmp8;
double Mstmp85 = x*M[13];
double Mstmp86 = y*M[11];
double Mstmp87 = pow(y, 3);
double Mstmp88 = (1.0/6.0)*Mstmp87;
double Mstmp89 = Mstmp29*y;
double Mstmp90 = x*M[14];
double Mstmp91 = y*M[12];
double Mstmp92 = z*M[11];
double Mstmp93 = Mstmp33*y;
double Mstmp94 = Mstmp29*z;
double Mstmp95 = Mstmp30*z;
double Mstmp96 = y*M[13];
double Mstmp97 = y*M[14];
double Mstmp98 = z*M[13];
double Mstmp99 = Mstmp39*z;
double Mstmp100 = x*M[15];
double Mstmp101 = (1.0/2.0)*M[8];
double Mstmp102 = z*M[17];
double Mstmp103 = Mstmp101*Mstmp14;
double Mstmp104 = pow(x, 4);
double Mstmp105 = (1.0/24.0)*M[0];
double Mstmp106 = Mstmp44*x;
double Mstmp107 = Mstmp10*Mstmp53;
double Mstmp108 = Mstmp13*Mstmp55;
double Mstmp109 = (1.0/4.0)*Mstmp14;
double Mstmp110 = Mstmp11*M[0];
double Mstmp111 = Mstmp109*Mstmp110;
double Mstmp112 = Mstmp1*Mstmp48;
double Mstmp113 = pow(z, 4);
double Mstmp114 = Mstmp105*Mstmp113 + Mstmp5*Mstmp51 + Mstmp51*M[5];
double Mstmp115 = x*M[16];
double Mstmp116 = y*M[15];
double Mstmp117 = (1.0/2.0)*M[9];
double Mstmp118 = z*M[19];
double Mstmp119 = Mstmp117*Mstmp14;
double Mstmp120 = (1.0/24.0)*M[1];
double Mstmp121 = Mstmp43*y;
double Mstmp122 = Mstmp61*x;
double Mstmp123 = Mstmp44*y;
double Mstmp124 = Mstmp17*Mstmp53;
double Mstmp125 = Mstmp18*Mstmp53;
double Mstmp126 = Mstmp21*Mstmp55;
double Mstmp127 = Mstmp109*Mstmp11;
double Mstmp128 = Mstmp127*M[1];
double Mstmp129 = Mstmp45*y;
double Mstmp130 = Mstmp20*Mstmp53;
double Mstmp131 = Mstmp22*Mstmp55;
double Mstmp132 = Mstmp113*Mstmp120 + Mstmp51*Mstmp8 + Mstmp51*M[7];
double Mstmp133 = x*M[17];
double Mstmp134 = (1.0/2.0)*M[10];
double Mstmp135 = (1.0/24.0)*M[2];
double Mstmp136 = Mstmp113*Mstmp135;
double Mstmp137 = -Mstmp0*Mstmp51 - Mstmp127*M[2] - Mstmp134*Mstmp14 - Mstmp24*Mstmp53 - Mstmp51*M[3];
double Mstmp138 = z*M[21];
double Mstmp139 = Mstmp70*x;
double Mstmp140 = Mstmp61*y;
double Mstmp141 = Mstmp62*y;
double Mstmp142 = (1.0/2.0)*M[11];
double Mstmp143 = Mstmp14*Mstmp142;
double Mstmp144 = Mstmp29*Mstmp53;
double Mstmp145 = Mstmp30*Mstmp53;
double Mstmp146 = Mstmp13*Mstmp75;
double Mstmp147 = Mstmp32*Mstmp53;
double Mstmp148 = Mstmp15*Mstmp75;
double Mstmp149 = Mstmp109*Mstmp31;
double Mstmp150 = Mstmp149*M[0];
double Mstmp151 = x*M[19];
double Mstmp152 = (1.0/2.0)*M[12];
double Mstmp153 = -Mstmp14*Mstmp152 - Mstmp3*Mstmp51 - Mstmp33*Mstmp53 - Mstmp34*Mstmp53 - Mstmp36*Mstmp53 - Mstmp4*Mstmp51 - Mstmp51*M[4];
double Mstmp154 = z*M[23];
double Mstmp155 = Mstmp70*y;
double Mstmp156 = (1.0/2.0)*M[13];
double Mstmp157 = Mstmp14*Mstmp156;
double Mstmp158 = Mstmp39*Mstmp53;
double Mstmp159 = Mstmp21*Mstmp75;
double Mstmp160 = Mstmp1*Mstmp88;
double Mstmp161 = Mstmp149*M[1];
double Mstmp162 = x*M[21];
double Mstmp163 = (1.0/2.0)*M[14];
double Mstmp164 = (1.0/12.0)*Mstmp113;
double Mstmp165 = x*M[11];
double Mstmp166 = y*M[9];
double Mstmp167 = (1.0/4.0)*Mstmp11*Mstmp31;
double Mstmp168 = Mstmp17*y;
double Mstmp169 = -Mstmp14*Mstmp163 - Mstmp149*M[2] - Mstmp40*Mstmp53 - Mstmp51*Mstmp7 - Mstmp51*M[6];
double Mstmp170 = x*M[22];
double Mstmp171 = y*M[20];
double Mstmp172 = pow(y, 4);
double Mstmp173 = Mstmp85*y;
double Mstmp174 = x*M[23];
double Mstmp175 = y*M[22];
double Mstmp176 = z*M[26];
double Mstmp177 = Mstmp102*x;
double Mstmp178 = (1.0/2.0)*M[15];
double Mstmp179 = (1.0/24.0)*M[3];
double Mstmp180 = pow(x, 5);
double Mstmp181 = (1.0/120.0)*Mstmp180;
double Mstmp182 = Mstmp14*Mstmp178;
double Mstmp183 = pow(z, 5);
double Mstmp184 = (1.0/120.0)*M[2];
double Mstmp185 = -Mstmp183*Mstmp184;
double Mstmp186 = Mstmp43*Mstmp53;
double Mstmp187 = (1.0/24.0)*Mstmp113;
double Mstmp188 = Mstmp44*Mstmp55;
double Mstmp189 = Mstmp13*Mstmp48;
double Mstmp190 = (1.0/24.0)*Mstmp104;
double Mstmp191 = Mstmp1*Mstmp190;
double Mstmp192 = Mstmp127*M[3];
double Mstmp193 = (1.0/12.0)*Mstmp50;
double Mstmp194 = Mstmp193*M[2];
double Mstmp195 = (1.0/12.0)*Mstmp47;
double Mstmp196 = Mstmp14*Mstmp195;
double Mstmp197 = Mstmp196*M[0];
double Mstmp198 = (1.0/2.0)*M[16];
double Mstmp199 = z*M[28];
double Mstmp200 = Mstmp14*Mstmp198;
double Mstmp201 = (1.0/24.0)*M[4];
double Mstmp202 = Mstmp118*x;
double Mstmp203 = Mstmp102*y;
double Mstmp204 = Mstmp53*Mstmp57;
double Mstmp205 = Mstmp53*Mstmp58;
double Mstmp206 = Mstmp55*Mstmp61;
double Mstmp207 = Mstmp127*M[4];
double Mstmp208 = Mstmp21*Mstmp48;
double Mstmp209 = Mstmp196*M[1];
double Mstmp210 = Mstmp106*y;
double Mstmp211 = Mstmp53*Mstmp60;
double Mstmp212 = Mstmp55*Mstmp63;
double Mstmp213 = Mstmp127*Mstmp4;
double Mstmp214 = Mstmp22*Mstmp48;
double Mstmp215 = Mstmp113*Mstmp201 + Mstmp187*Mstmp3 + Mstmp187*Mstmp4 + Mstmp33*Mstmp51 + Mstmp34*Mstmp51 + Mstmp36*Mstmp51 + Mstmp51*M[12];
double Mstmp216 = (1.0/2.0)*M[17];
double Mstmp217 = (1.0/24.0)*M[5];
double Mstmp218 = (1.0/120.0)*Mstmp183;
double Mstmp219 = Mstmp113*Mstmp217 + Mstmp187*Mstmp5 + Mstmp218*M[0];
double Mstmp220 = -Mstmp10*Mstmp51 - Mstmp110*Mstmp193 - Mstmp127*M[5] - Mstmp14*Mstmp216 - Mstmp196*M[2] - Mstmp51*M[8] - Mstmp53*Mstmp65;
double Mstmp221 = Mstmp0*Mstmp149 + Mstmp1*Mstmp167 + Mstmp118*y + Mstmp122*y + Mstmp127*Mstmp7 + Mstmp127*M[6] + Mstmp138*x + (1.0/2.0)*Mstmp14*M[18] + Mstmp149*M[3] + Mstmp165*Mstmp53 + Mstmp166*Mstmp53 + Mstmp168*Mstmp53 + (1.0/40.0)*Mstmp183*M[2] + Mstmp44*Mstmp75 + Mstmp45*Mstmp75 + Mstmp55*Mstmp70 + Mstmp55*Mstmp71 + z*M[30];
double Mstmp222 = (1.0/2.0)*M[19];
double Mstmp223 = (1.0/24.0)*M[7];
double Mstmp224 = Mstmp113*Mstmp223 + Mstmp187*Mstmp8 + Mstmp218*M[1];
double Mstmp225 = -Mstmp11*Mstmp193*M[1] - Mstmp127*Mstmp8 - Mstmp127*M[7] - Mstmp14*Mstmp222 - Mstmp17*Mstmp51 - Mstmp18*Mstmp51 - Mstmp20*Mstmp51 - Mstmp51*M[9] - Mstmp53*Mstmp77 - Mstmp53*Mstmp78 - Mstmp53*Mstmp81;
double Mstmp226 = z*M[32];
double Mstmp227 = Mstmp154*x;
double Mstmp228 = Mstmp138*y;
double Mstmp229 = Mstmp139*y;
double Mstmp230 = (1.0/2.0)*M[20];
double Mstmp231 = Mstmp14*Mstmp230;
double Mstmp232 = Mstmp53*Mstmp85;
double Mstmp233 = Mstmp53*Mstmp86;
double Mstmp234 = Mstmp61*Mstmp75;
double Mstmp235 = Mstmp13*Mstmp88;
double Mstmp236 = Mstmp53*Mstmp89;
double Mstmp237 = Mstmp62*Mstmp75;
double Mstmp238 = Mstmp15*Mstmp88;
double Mstmp239 = Mstmp149*M[4];
double Mstmp240 = (1.0/12.0)*Mstmp87;
double Mstmp241 = Mstmp14*Mstmp240;
double Mstmp242 = Mstmp241*M[0];
double Mstmp243 = Mstmp149*Mstmp3;
double Mstmp244 = (1.0/2.0)*M[21];
double Mstmp245 = (1.0/60.0)*Mstmp183;
double Mstmp246 = Mstmp193*Mstmp31;
double Mstmp247 = -Mstmp14*Mstmp244 - Mstmp149*Mstmp5 - Mstmp149*M[5] - Mstmp246*M[0] - Mstmp29*Mstmp51 - Mstmp30*Mstmp51 - Mstmp32*Mstmp51 - Mstmp51*M[11] - Mstmp53*Mstmp90 - Mstmp53*Mstmp91 - Mstmp53*Mstmp93;
double Mstmp248 = z*M[34];
double Mstmp249 = (1.0/2.0)*M[22];
double Mstmp250 = Mstmp14*Mstmp249;
double Mstmp251 = Mstmp154*y;
double Mstmp252 = Mstmp53*Mstmp96;
double Mstmp253 = Mstmp70*Mstmp75;
double Mstmp254 = Mstmp149*M[6];
double Mstmp255 = Mstmp21*Mstmp88;
double Mstmp256 = Mstmp241*M[1];
double Mstmp257 = (1.0/24.0)*Mstmp172;
double Mstmp258 = Mstmp1*Mstmp257;
double Mstmp259 = (1.0/2.0)*M[23];
double Mstmp260 = -Mstmp14*Mstmp259 - Mstmp149*M[7] - Mstmp241*M[2] - Mstmp246*M[1] - Mstmp39*Mstmp51 - Mstmp51*M[13] - Mstmp53*Mstmp97;
double Mstmp261 = pow(y, 5);
double Mstmp262 = (1.0/120.0)*Mstmp261;
double Mstmp263 = (1.0/24.0)*M[6];
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += Mstmp0 + Mstmp2 + M[3];
#pragma omp atomic
Ms[4] += Mstmp3 + Mstmp4 + M[4];
#pragma omp atomic
Ms[5] += Mstmp5 + Mstmp6 + M[5];
#pragma omp atomic
Ms[6] += Mstmp2 + Mstmp7 + M[6];
#pragma omp atomic
Ms[7] += Mstmp8 + Mstmp9 + M[7];
#pragma omp atomic
Ms[8] += Mstmp10 + Mstmp11*Mstmp12 + Mstmp16 + M[8];
#pragma omp atomic
Ms[9] += Mstmp11*Mstmp19 + Mstmp17 + Mstmp18 + Mstmp20 + Mstmp23 + M[9];
#pragma omp atomic
Ms[10] += Mstmp11*Mstmp27 + Mstmp24 + Mstmp25 + Mstmp26 + Mstmp28 + M[10];
#pragma omp atomic
Ms[11] += Mstmp12*Mstmp31 + Mstmp16 + Mstmp29 + Mstmp30 + Mstmp32 + M[11];
#pragma omp atomic
Ms[12] += Mstmp33 + Mstmp34 + Mstmp35 + Mstmp36 + Mstmp37 + Mstmp38 + M[12];
#pragma omp atomic
Ms[13] += Mstmp19*Mstmp31 + Mstmp23 + Mstmp39 + M[13];
#pragma omp atomic
Ms[14] += Mstmp27*Mstmp31 + Mstmp28 + Mstmp40 + Mstmp41 + Mstmp42 + M[14];
#pragma omp atomic
Ms[15] += Mstmp11*Mstmp46 + Mstmp43 - Mstmp44 - Mstmp45 + Mstmp48*M[0] - Mstmp49 + Mstmp52 - Mstmp54 - Mstmp56 + M[15];
#pragma omp atomic
Ms[16] += Mstmp11*Mstmp59 + Mstmp4*Mstmp55 + Mstmp48*M[1] + Mstmp57 + Mstmp58 + Mstmp60 + Mstmp64 + M[16];
#pragma omp atomic
Ms[17] += Mstmp11*Mstmp67 + Mstmp48*M[2] + Mstmp55*Mstmp6 + Mstmp65 + Mstmp66 + Mstmp68 + Mstmp69 + M[17];
#pragma omp atomic
Ms[18] += (1.0/2.0)*Mstmp11*y*M[1] + (1.0/2.0)*Mstmp11*M[6] + (1.0/2.0)*Mstmp31*x*M[0] + (1.0/2.0)*Mstmp31*M[3] - Mstmp44 - Mstmp45 - Mstmp49 + (1.0/3.0)*Mstmp50*M[2] - Mstmp54 - Mstmp56 - Mstmp70 - Mstmp71 - Mstmp73 - Mstmp74 - Mstmp76 + x*y*M[4] + x*M[11] + y*M[9] + M[18];
#pragma omp atomic
Ms[19] += Mstmp11*Mstmp80 + Mstmp20*z + Mstmp55*Mstmp8 + Mstmp55*Mstmp9 + Mstmp77 + Mstmp78 + Mstmp79 + Mstmp81 + Mstmp82 + Mstmp83 + Mstmp84 + M[19];
#pragma omp atomic
Ms[20] += Mstmp3*Mstmp75 + Mstmp31*Mstmp59 + Mstmp64 + Mstmp85 + Mstmp86 + Mstmp88*M[0] + Mstmp89 + M[20];
#pragma omp atomic
Ms[21] += Mstmp31*Mstmp67 + Mstmp32*z + Mstmp5*Mstmp75 + Mstmp6*Mstmp75 + Mstmp69 + Mstmp90 + Mstmp91 + Mstmp92 + Mstmp93 + Mstmp94 + Mstmp95 + M[21];
#pragma omp atomic
Ms[22] += Mstmp31*Mstmp72 + Mstmp52 - Mstmp70 - Mstmp71 - Mstmp73 - Mstmp74 - Mstmp76 + Mstmp88*M[1] + Mstmp96 + M[22];
#pragma omp atomic
Ms[23] += Mstmp31*Mstmp80 + Mstmp75*Mstmp9 + Mstmp84 + Mstmp88*M[2] + Mstmp97 + Mstmp98 + Mstmp99 + M[23];
#pragma omp atomic
Ms[24] += Mstmp100 + Mstmp101*Mstmp11 - Mstmp102 - Mstmp103 + Mstmp104*Mstmp105 - Mstmp106 - Mstmp107 - Mstmp108 - Mstmp111 - Mstmp112 + Mstmp114 + Mstmp48*M[3] + M[24];
#pragma omp atomic
Ms[25] += Mstmp104*Mstmp120 + Mstmp11*Mstmp117 + Mstmp115 + Mstmp116 - Mstmp118 - Mstmp119 + Mstmp121 - Mstmp122 - Mstmp123 - Mstmp124 - Mstmp125 - Mstmp126 - Mstmp128 - Mstmp129 - Mstmp130 - Mstmp131 + Mstmp132 + Mstmp18*Mstmp55 + Mstmp4*Mstmp48 + Mstmp48*M[4] + M[25];
#pragma omp atomic
Ms[26] += Mstmp104*Mstmp135 + Mstmp11*Mstmp134 + Mstmp133 + Mstmp136 + Mstmp137 + Mstmp25*Mstmp55 + Mstmp43*z + Mstmp48*Mstmp6 + Mstmp48*M[5] + z*M[15] + M[26];
#pragma omp atomic
Ms[27] += -Mstmp102 - Mstmp103 - Mstmp106 - Mstmp107 - Mstmp108 + (1.0/4.0)*Mstmp11*Mstmp31*M[0] + (1.0/2.0)*Mstmp11*y*M[4] + (1.0/2.0)*Mstmp11*M[11] - Mstmp111 - Mstmp112 + (1.0/12.0)*Mstmp113*M[0] - Mstmp138 - Mstmp139 - Mstmp140 - Mstmp141 - Mstmp143 - Mstmp144 - Mstmp145 - Mstmp146 - Mstmp147 - Mstmp148 - Mstmp150 + (1.0/2.0)*Mstmp31*x*M[3] + (1.0/2.0)*Mstmp31*M[8] + (1.0/6.0)*Mstmp47*y*M[1] + (1.0/6.0)*Mstmp47*M[6] + (1.0/3.0)*Mstmp50*x*M[2] + (1.0/3.0)*Mstmp50*M[5] + x*y*M[9] + x*M[18] + y*M[16] + M[27];
#pragma omp atomic
Ms[28] += Mstmp11*Mstmp152 + Mstmp151 + Mstmp153 + Mstmp34*Mstmp55 + Mstmp35*Mstmp55 + Mstmp38*Mstmp55 + Mstmp48*Mstmp8 + Mstmp48*Mstmp9 + Mstmp48*M[7] + Mstmp57*z + Mstmp58*z + Mstmp60*z + Mstmp65*y + y*M[17] + z*M[16] + M[28];
#pragma omp atomic
Ms[29] += (1.0/4.0)*Mstmp11*Mstmp31*M[1] + (1.0/2.0)*Mstmp11*y*M[6] + (1.0/2.0)*Mstmp11*M[13] + (1.0/12.0)*Mstmp113*M[1] - Mstmp118 - Mstmp119 - Mstmp122 - Mstmp123 - Mstmp124 - Mstmp125 - Mstmp126 - Mstmp128 - Mstmp129 - Mstmp130 - Mstmp131 - Mstmp154 - Mstmp155 - Mstmp157 - Mstmp158 - Mstmp159 - Mstmp160 - Mstmp161 + (1.0/2.0)*Mstmp31*x*M[4] + (1.0/2.0)*Mstmp31*M[9] + (1.0/3.0)*Mstmp50*y*M[2] + (1.0/3.0)*Mstmp50*M[7] + (1.0/6.0)*Mstmp87*x*M[0] + (1.0/6.0)*Mstmp87*M[3] + x*y*M[11] + x*M[20] + y*M[18] + M[29];
#pragma omp atomic
Ms[30] += Mstmp11*Mstmp163 + Mstmp134*Mstmp31 + Mstmp137 + Mstmp162 + Mstmp164*M[2] + Mstmp165*z + Mstmp166*z + Mstmp167*M[2] + Mstmp168*z + Mstmp169 + Mstmp24*Mstmp75 + Mstmp25*Mstmp75 + Mstmp26*Mstmp75 + Mstmp40*Mstmp55 + Mstmp41*Mstmp55 + Mstmp42*Mstmp55 + Mstmp77*y + y*M[19] + z*M[18] + M[30];
#pragma omp atomic
Ms[31] += Mstmp105*Mstmp172 + Mstmp114 - Mstmp138 - Mstmp139 - Mstmp140 - Mstmp141 + Mstmp142*Mstmp31 - Mstmp143 - Mstmp144 - Mstmp145 - Mstmp146 - Mstmp147 - Mstmp148 - Mstmp150 + Mstmp170 + Mstmp171 + Mstmp173 + Mstmp29*Mstmp75 + Mstmp3*Mstmp88 + Mstmp88*M[4] + M[31];
#pragma omp atomic
Ms[32] += Mstmp152*Mstmp31 + Mstmp153 + Mstmp174 + Mstmp33*Mstmp75 + Mstmp35*Mstmp75 + Mstmp37*Mstmp75 + Mstmp5*Mstmp88 + Mstmp6*Mstmp88 + Mstmp85*z + Mstmp86*z + Mstmp88*M[5] + Mstmp89*z + Mstmp90*y + y*M[21] + z*M[20] + M[32];
#pragma omp atomic
Ms[33] += Mstmp120*Mstmp172 + Mstmp132 - Mstmp154 - Mstmp155 + Mstmp156*Mstmp31 - Mstmp157 - Mstmp158 - Mstmp159 - Mstmp160 - Mstmp161 + Mstmp175 + Mstmp88*M[6] + M[33];
#pragma omp atomic
Ms[34] += Mstmp135*Mstmp172 + Mstmp136 + Mstmp163*Mstmp31 + Mstmp169 + Mstmp41*Mstmp75 + Mstmp88*Mstmp9 + Mstmp88*M[7] + Mstmp96*z + y*M[23] + z*M[22] + M[34];
#pragma omp atomic
Ms[35] += Mstmp0*Mstmp187 + Mstmp104*Mstmp179 + Mstmp11*Mstmp178 + Mstmp11*Mstmp194 + Mstmp113*Mstmp179 - Mstmp176 - Mstmp177 + Mstmp181*M[0] - Mstmp182 + Mstmp185 - Mstmp186 - Mstmp188 - Mstmp189 - Mstmp191 - Mstmp192 - Mstmp197 + Mstmp24*Mstmp51 + Mstmp48*M[8] + Mstmp51*M[10] + x*M[24] + M[35];
#pragma omp atomic
Ms[36] += Mstmp100*y + Mstmp104*Mstmp201 + Mstmp11*Mstmp198 + Mstmp18*Mstmp48 + Mstmp181*M[1] + Mstmp190*Mstmp4 - Mstmp199 - Mstmp200 - Mstmp202 - Mstmp203 - Mstmp204 - Mstmp205 - Mstmp206 - Mstmp207 - Mstmp208 - Mstmp209 - Mstmp210 - Mstmp211 - Mstmp212 - Mstmp213 - Mstmp214 + Mstmp215 + Mstmp48*M[9] + Mstmp55*Mstmp58 + x*M[25] + y*M[24] + M[36];
#pragma omp atomic
Ms[37] += Mstmp100*z + Mstmp104*Mstmp217 + Mstmp11*Mstmp216 + Mstmp180*Mstmp184 + Mstmp190*Mstmp6 + Mstmp219 + Mstmp220 + Mstmp25*Mstmp48 + Mstmp48*M[10] + Mstmp55*Mstmp66 + x*M[26] + z*M[24] + M[37];
#pragma omp atomic
Ms[38] += (1.0/24.0)*Mstmp104*y*M[1] + (1.0/24.0)*Mstmp104*M[6] + (1.0/4.0)*Mstmp11*Mstmp31*M[3] + (1.0/6.0)*Mstmp11*Mstmp50*M[2] + (1.0/2.0)*Mstmp11*y*M[9] + (1.0/2.0)*Mstmp11*M[18] + (1.0/12.0)*Mstmp113*x*M[0] + (1.0/24.0)*Mstmp113*y*M[1] + (1.0/12.0)*Mstmp113*M[3] + (1.0/24.0)*Mstmp113*M[6] - Mstmp176 - Mstmp177 - Mstmp182 - Mstmp186 - Mstmp188 - Mstmp189 - Mstmp191 - Mstmp192 - Mstmp197 - Mstmp221 + (1.0/12.0)*Mstmp31*Mstmp47*M[0] + (1.0/12.0)*Mstmp31*Mstmp50*M[2] + (1.0/2.0)*Mstmp31*x*M[8] + (1.0/2.0)*Mstmp31*M[15] + (1.0/6.0)*Mstmp47*y*M[4] + (1.0/6.0)*Mstmp47*M[11] + (1.0/3.0)*Mstmp50*x*M[5] + (1.0/6.0)*Mstmp50*y*M[7] + (1.0/3.0)*Mstmp50*M[10] + (1.0/6.0)*Mstmp50*M[14] + x*y*M[16] + x*M[27] + y*M[25] + M[38];
#pragma omp atomic
Ms[39] += Mstmp104*Mstmp223 + Mstmp11*Mstmp222 + Mstmp115*z + Mstmp116*z + Mstmp121*z + Mstmp133*y + Mstmp190*Mstmp8 + Mstmp190*Mstmp9 + Mstmp224 + Mstmp225 + Mstmp34*Mstmp48 + Mstmp35*Mstmp48 + Mstmp38*Mstmp48 + Mstmp48*M[12] + Mstmp55*Mstmp78 + Mstmp55*Mstmp79 + Mstmp55*Mstmp83 + x*M[28] + y*M[26] + z*M[25] + M[39];
#pragma omp atomic
Ms[40] += (1.0/4.0)*Mstmp11*Mstmp31*M[4] + (1.0/12.0)*Mstmp11*Mstmp87*M[0] + (1.0/2.0)*Mstmp11*y*M[11] + (1.0/2.0)*Mstmp11*M[20] + (1.0/12.0)*Mstmp113*x*M[1] + (1.0/12.0)*Mstmp113*y*M[0] + (1.0/12.0)*Mstmp113*M[4] - Mstmp199 - Mstmp200 - Mstmp202 - Mstmp203 - Mstmp204 - Mstmp205 - Mstmp206 - Mstmp207 - Mstmp208 - Mstmp209 - Mstmp210 - Mstmp211 - Mstmp212 - Mstmp213 - Mstmp214 - Mstmp226 - Mstmp227 - Mstmp228 - Mstmp229 - Mstmp231 - Mstmp232 - Mstmp233 - Mstmp234 - Mstmp235 - Mstmp236 - Mstmp237 - Mstmp238 - Mstmp239 - Mstmp242 - Mstmp243 + (1.0/12.0)*Mstmp31*Mstmp47*M[1] + (1.0/2.0)*Mstmp31*x*M[9] + (1.0/2.0)*Mstmp31*M[16] + (1.0/6.0)*Mstmp47*y*M[6] + (1.0/6.0)*Mstmp47*M[13] + (1.0/3.0)*Mstmp50*x*y*M[2] + (1.0/3.0)*Mstmp50*x*M[7] + (1.0/3.0)*Mstmp50*y*M[5] + (1.0/3.0)*Mstmp50*M[12] + (1.0/6.0)*Mstmp87*x*M[3] + (1.0/6.0)*Mstmp87*M[8] + x*y*M[18] + x*M[29] + y*M[27] + M[40];
#pragma omp atomic
Ms[41] += Mstmp11*Mstmp244 + Mstmp151*y + Mstmp164*Mstmp5 + Mstmp164*M[5] + Mstmp167*Mstmp6 + Mstmp167*M[5] + Mstmp195*Mstmp31*M[2] + Mstmp216*Mstmp31 + Mstmp220 + Mstmp245*M[0] + Mstmp247 + Mstmp40*Mstmp48 + Mstmp41*Mstmp48 + Mstmp42*Mstmp48 + Mstmp48*M[14] + Mstmp55*Mstmp91 + Mstmp55*Mstmp92 + Mstmp55*Mstmp95 + Mstmp57*y*z + Mstmp65*Mstmp75 + Mstmp66*Mstmp75 + Mstmp68*Mstmp75 + x*z*M[18] + x*M[30] + y*z*M[16] + y*M[28] + z*M[27] + M[41];
#pragma omp atomic
Ms[42] += (1.0/4.0)*Mstmp11*Mstmp31*M[6] + (1.0/12.0)*Mstmp11*Mstmp50*M[2] + (1.0/12.0)*Mstmp11*Mstmp87*M[1] + (1.0/2.0)*Mstmp11*y*M[13] + (1.0/2.0)*Mstmp11*M[22] + (1.0/24.0)*Mstmp113*x*M[0] + (1.0/12.0)*Mstmp113*y*M[1] + (1.0/24.0)*Mstmp113*M[3] + (1.0/12.0)*Mstmp113*M[6] + (1.0/24.0)*Mstmp172*x*M[0] + (1.0/24.0)*Mstmp172*M[3] - Mstmp221 - Mstmp248 - Mstmp250 - Mstmp251 - Mstmp252 - Mstmp253 - Mstmp254 - Mstmp255 - Mstmp256 - Mstmp258 + (1.0/6.0)*Mstmp31*Mstmp50*M[2] + (1.0/2.0)*Mstmp31*x*M[11] + (1.0/2.0)*Mstmp31*M[18] + (1.0/6.0)*Mstmp50*x*M[5] + (1.0/3.0)*Mstmp50*y*M[7] + (1.0/6.0)*Mstmp50*M[10] + (1.0/3.0)*Mstmp50*M[14] + (1.0/6.0)*Mstmp87*x*M[4] + (1.0/6.0)*Mstmp87*M[9] + x*y*M[20] + x*M[31] + y*M[29] + M[42];
#pragma omp atomic
Ms[43] += Mstmp11*Mstmp240*M[2] + Mstmp11*Mstmp259 + Mstmp162*y + Mstmp164*Mstmp8 + Mstmp164*M[7] + Mstmp165*y*z + Mstmp167*Mstmp9 + Mstmp167*M[7] + Mstmp222*Mstmp31 + Mstmp225 + Mstmp24*Mstmp88 + Mstmp245*M[1] + Mstmp25*Mstmp88 + Mstmp26*Mstmp88 + Mstmp260 + Mstmp55*Mstmp97 + Mstmp55*Mstmp98 + Mstmp55*Mstmp99 + Mstmp75*Mstmp77 + Mstmp75*Mstmp79 + Mstmp75*Mstmp82 + Mstmp88*M[10] + x*z*M[20] + x*M[32] + y*z*M[18] + y*M[30] + z*M[29] + M[43];
#pragma omp atomic
Ms[44] += Mstmp170*y + Mstmp172*Mstmp201 + Mstmp215 - Mstmp226 - Mstmp227 - Mstmp228 - Mstmp229 + Mstmp230*Mstmp31 - Mstmp231 - Mstmp232 - Mstmp233 - Mstmp234 - Mstmp235 - Mstmp236 - Mstmp237 - Mstmp238 - Mstmp239 - Mstmp242 - Mstmp243 + Mstmp257*Mstmp3 + Mstmp262*M[0] + Mstmp29*Mstmp88 + Mstmp75*Mstmp85 + Mstmp88*M[11] + x*M[33] + y*M[31] + M[44];
#pragma omp atomic
Ms[45] += Mstmp170*z + Mstmp171*z + Mstmp172*Mstmp217 + Mstmp173*z + Mstmp174*y + Mstmp219 + Mstmp244*Mstmp31 + Mstmp247 + Mstmp257*Mstmp5 + Mstmp257*Mstmp6 + Mstmp33*Mstmp88 + Mstmp35*Mstmp88 + Mstmp37*Mstmp88 + Mstmp75*Mstmp90 + Mstmp75*Mstmp92 + Mstmp75*Mstmp94 + Mstmp88*M[12] + x*M[34] + y*M[32] + z*M[31] + M[45];
#pragma omp atomic
Ms[46] += Mstmp113*Mstmp263 + Mstmp172*Mstmp263 + Mstmp185 + Mstmp187*Mstmp7 + Mstmp194*Mstmp31 - Mstmp248 + Mstmp249*Mstmp31 - Mstmp250 - Mstmp251 - Mstmp252 - Mstmp253 - Mstmp254 - Mstmp255 - Mstmp256 - Mstmp258 + Mstmp262*M[1] + Mstmp40*Mstmp51 + Mstmp51*M[14] + Mstmp88*M[13] + y*M[33] + M[46];
#pragma omp atomic
Ms[47] += Mstmp172*Mstmp223 + Mstmp175*z + Mstmp184*Mstmp261 + Mstmp224 + Mstmp257*Mstmp9 + Mstmp259*Mstmp31 + Mstmp260 + Mstmp41*Mstmp88 + Mstmp75*Mstmp98 + Mstmp88*M[14] + y*M[34] + z*M[33] + M[47];

}

void L2Lc_6(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[13];
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
double Lstmp15 = x*L[12];
double Lstmp16 = x*L[21];
double Lstmp17 = x*L[32];
double Lstmp18 = y*L[10];
double Lstmp19 = z*L[11];
double Lstmp20 = y*L[17];
double Lstmp21 = z*L[18];
double Lstmp22 = y*L[26];
double Lstmp23 = z*L[27];
double Lstmp24 = z*L[15];
double Lstmp25 = z*L[24];
double Lstmp26 = z*L[35];
double Lstmp27 = z*L[22];
double Lstmp28 = Lstmp27*x;
double Lstmp29 = z*L[33];
double Lstmp30 = Lstmp29*x;
double Lstmp31 = z*L[20];
double Lstmp32 = Lstmp31*y;
double Lstmp33 = z*L[29];
double Lstmp34 = Lstmp33*y;
double Lstmp35 = (1.0/4.0)*Lstmp10*Lstmp5;
double Lstmp36 = z*L[31];
double Lstmp37 = L[4] + L[7];
double Lstmp38 = pow(z, 2);
double Lstmp39 = (1.0/2.0)*Lstmp38;
double Lstmp40 = L[11] + L[15];
double Lstmp41 = pow(z, 3);
double Lstmp42 = (1.0/6.0)*Lstmp41;
double Lstmp43 = L[9] + L[12];
double Lstmp44 = Lstmp39*Lstmp43;
double Lstmp45 = L[18] + L[22];
double Lstmp46 = Lstmp42*Lstmp45;
double Lstmp47 = L[10] + L[14];
double Lstmp48 = Lstmp39*Lstmp47;
double Lstmp49 = L[20] + L[24];
double Lstmp50 = Lstmp42*Lstmp49;
double Lstmp51 = L[17] + L[21];
double Lstmp52 = Lstmp39*Lstmp51;
double Lstmp53 = Lstmp52*y;
double Lstmp54 = L[29] + L[33];
double Lstmp55 = Lstmp42*Lstmp54;
double Lstmp56 = Lstmp55*y;
double Lstmp57 = L[16] + L[19];
double Lstmp58 = (1.0/4.0)*Lstmp38;
double Lstmp59 = Lstmp5*Lstmp58;
double Lstmp60 = L[27] + L[31];
double Lstmp61 = (1.0/12.0)*Lstmp41;
double Lstmp62 = L[25] + L[28];
double Lstmp63 = (1.0/12.0)*Lstmp38;
double Lstmp64 = L[19] + L[23];
double Lstmp65 = Lstmp10*Lstmp58;
double Lstmp66 = L[31] + L[35];
double Lstmp67 = L[30] + L[34];
double Lstmp68 = L[28] + L[32];
double Lstmp69 = Lstmp65*Lstmp68;
double Lstmp70 = L[26] + L[30];
double Lstmp71 = Lstmp59*Lstmp70;
double Lstmp72 = L[16] + 2*L[19] + L[23];
double Lstmp73 = (1.0/24.0)*pow(z, 4);
double Lstmp74 = L[27] + 2*L[31] + L[35];
double Lstmp75 = L[25] + 2*L[28] + L[32];
double Lstmp76 = Lstmp73*Lstmp75;
double Lstmp77 = L[26] + 2*L[30] + L[34];
double Lstmp78 = Lstmp73*Lstmp77;
double Lstmp79 = x*L[19];
double Lstmp80 = x*L[30];
double Lstmp81 = Lstmp36*x;
double Lstmp82 = Lstmp39*Lstmp57;
double Lstmp83 = Lstmp42*Lstmp60;
double Lstmp84 = Lstmp39*Lstmp70;
double Lstmp85 = Lstmp84*y;
double Lstmp86 = y*L[12];
double Lstmp87 = Lstmp27*y;
double Lstmp88 = y*L[19];
double Lstmp89 = y*L[28];
double Lstmp90 = Lstmp36*y;
double Lstmp91 = Lstmp39*Lstmp64;
double Lstmp92 = Lstmp42*Lstmp66;
double Lstmp93 = Lstmp39*Lstmp68;
double Lstmp94 = Lstmp93*y;
double Lstmp95 = y*L[13];
double Lstmp96 = x*L[22];
double Lstmp97 = x*L[33];
double Lstmp98 = y*L[20];
double Lstmp99 = y*L[29];
double Lstmp100 = Lstmp43*z;
double Lstmp101 = Lstmp47*z;
double Lstmp102 = Lstmp51*z;
double Lstmp103 = Lstmp102*y;
double Lstmp104 = Lstmp39*Lstmp45;
double Lstmp105 = Lstmp57*z;
double Lstmp106 = Lstmp62*z;
double Lstmp107 = Lstmp39*Lstmp49;
double Lstmp108 = Lstmp64*z;
double Lstmp109 = Lstmp67*z;
double Lstmp110 = Lstmp39*Lstmp54;
double Lstmp111 = Lstmp110*y;
double Lstmp112 = Lstmp68*z;
double Lstmp113 = Lstmp112*x;
double Lstmp114 = Lstmp70*z;
double Lstmp115 = Lstmp114*y;
double Lstmp116 = Lstmp42*Lstmp75;
double Lstmp117 = Lstmp42*Lstmp77;
double Lstmp118 = x*L[28];
double Lstmp119 = Lstmp39*Lstmp62;
double Lstmp120 = x*L[31];
double Lstmp121 = Lstmp39*Lstmp60;
double Lstmp122 = y*L[21];
double Lstmp123 = Lstmp29*y;
double Lstmp124 = y*L[30];
double Lstmp125 = Lstmp39*Lstmp67;
double Lstmp126 = y*L[22];
double Lstmp127 = y*L[31];
double Lstmp128 = Lstmp112*y;
double Lstmp129 = Lstmp39*Lstmp66;
double Lstmp130 = y*L[32];
double Lstmp131 = y*L[33];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x - Lstmp10*Lstmp61*Lstmp66 + (1.0/12.0)*Lstmp10*Lstmp7*L[28] + Lstmp11*Lstmp15 + Lstmp11*Lstmp24 + Lstmp11*Lstmp28 + Lstmp11*L[7] + (1.0/12.0)*Lstmp12*Lstmp5*L[30] - Lstmp12*Lstmp63*Lstmp67 + Lstmp13*Lstmp16 + Lstmp13*Lstmp25 + Lstmp13*Lstmp30 + Lstmp13*L[14] + Lstmp14*Lstmp17 + Lstmp14*Lstmp26 + Lstmp14*L[23] + Lstmp18*Lstmp6 + Lstmp19*Lstmp6 + Lstmp2*y + Lstmp20*Lstmp8 + Lstmp21*Lstmp8 + Lstmp22*Lstmp9 + Lstmp23*Lstmp9 + Lstmp32*Lstmp6 + Lstmp34*Lstmp8 + Lstmp35*Lstmp36 + Lstmp35*L[19] - Lstmp37*Lstmp39 + Lstmp4*x - Lstmp40*Lstmp42 - Lstmp44*x - Lstmp46*x - Lstmp48*y - Lstmp5*Lstmp60*Lstmp61 - Lstmp50*y - Lstmp53*x - Lstmp56*x - Lstmp57*Lstmp59 + Lstmp6*L[4] - Lstmp62*Lstmp63*Lstmp7 - Lstmp64*Lstmp65 - Lstmp69*x - Lstmp71*y + Lstmp72*Lstmp73 + (1.0/120.0)*Lstmp74*pow(z, 5) + Lstmp76*x + Lstmp78*y + Lstmp8*L[9] + Lstmp9*L[16] + (1.0/120.0)*pow(x, 5)*L[25] + x*L[1] + (1.0/120.0)*pow(y, 5)*L[34] + y*L[2] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + Lstmp11*Lstmp27 + Lstmp11*Lstmp79 + Lstmp11*Lstmp81 + Lstmp11*L[12] + Lstmp13*Lstmp29 + Lstmp13*Lstmp80 + Lstmp13*L[21] + Lstmp14*L[32] + Lstmp18*x + Lstmp19*x + Lstmp20*Lstmp6 + Lstmp21*Lstmp6 + Lstmp22*Lstmp8 + Lstmp23*Lstmp8 + Lstmp32*x + Lstmp34*Lstmp6 + Lstmp35*L[28] + Lstmp4 - Lstmp44 - Lstmp46 - Lstmp53 - Lstmp56 - Lstmp59*Lstmp62 + Lstmp6*L[9] - Lstmp69 + Lstmp76 + Lstmp8*L[16] - Lstmp82*x - Lstmp83*x - Lstmp85*x + Lstmp9*L[25] + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp11*Lstmp16 + Lstmp11*Lstmp25 + Lstmp11*Lstmp30 + Lstmp11*L[14] + Lstmp13*Lstmp17 + Lstmp13*Lstmp26 + Lstmp13*L[23] + Lstmp14*L[34] + Lstmp2 + Lstmp24*y + Lstmp3*x + Lstmp31*Lstmp6 + Lstmp33*Lstmp8 + Lstmp35*L[30] - Lstmp48 - Lstmp50 - Lstmp52*x - Lstmp55*x + Lstmp6*Lstmp88 + Lstmp6*Lstmp90 + Lstmp6*L[10] - Lstmp65*Lstmp67 - Lstmp71 + Lstmp78 + Lstmp8*Lstmp89 + Lstmp8*L[17] + Lstmp86*x + Lstmp87*x + Lstmp9*L[26] - Lstmp91*y - Lstmp92*y - Lstmp94*x + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += -Lstmp100*x - Lstmp101*y - Lstmp103*x - Lstmp104*x - Lstmp105*Lstmp6 - Lstmp106*Lstmp8 - Lstmp107*y - Lstmp108*Lstmp11 - Lstmp109*Lstmp13 - Lstmp11*Lstmp113 + Lstmp11*Lstmp96 + Lstmp11*L[15] - Lstmp111*x - Lstmp115*Lstmp6 + Lstmp116*x + Lstmp117*y + Lstmp13*Lstmp97 + Lstmp13*L[24] + Lstmp14*L[35] + Lstmp35*L[31] - Lstmp37*z - Lstmp39*Lstmp40 + Lstmp42*Lstmp72 - Lstmp59*Lstmp60 + Lstmp6*Lstmp98 + Lstmp6*L[11] - Lstmp65*Lstmp66 + Lstmp73*Lstmp74 + Lstmp8*Lstmp99 + Lstmp8*L[18] + Lstmp9*L[27] + Lstmp95*x + x*L[6] + y*L[8] + L[3];
#pragma omp atomic
Ls[4] += Lstmp11*Lstmp118 + Lstmp11*Lstmp36 + Lstmp11*L[19] - Lstmp119*x + Lstmp13*L[30] + Lstmp18 + Lstmp19 + Lstmp20*x + Lstmp21*x + Lstmp22*Lstmp6 + Lstmp23*Lstmp6 + Lstmp32 + Lstmp34*x + Lstmp6*L[16] + Lstmp8*L[25] - Lstmp82 - Lstmp83 - Lstmp85 + x*L[9] + L[4];
#pragma omp atomic
Ls[5] += Lstmp11*Lstmp29 + Lstmp11*Lstmp80 + Lstmp11*L[21] + Lstmp13*L[32] + Lstmp3 + Lstmp31*x + Lstmp33*Lstmp6 - Lstmp52 - Lstmp55 + Lstmp6*Lstmp89 + Lstmp6*L[17] + Lstmp8*L[26] - Lstmp84*x + Lstmp86 + Lstmp87 + Lstmp88*x + Lstmp90*x - Lstmp94 + x*L[10] + L[5];
#pragma omp atomic
Ls[6] += -Lstmp100 - Lstmp103 - Lstmp104 - Lstmp105*x - Lstmp106*Lstmp6 - Lstmp11*Lstmp112 + Lstmp11*Lstmp120 + Lstmp11*L[22] - Lstmp111 - Lstmp115*x + Lstmp116 - Lstmp121*x + Lstmp13*L[33] + Lstmp6*Lstmp99 + Lstmp6*L[18] + Lstmp8*L[27] + Lstmp95 + Lstmp98*x + x*L[11] + L[6];
#pragma omp atomic
Ls[7] += Lstmp11*Lstmp17 + Lstmp11*Lstmp26 + Lstmp11*L[23] + Lstmp122*x + Lstmp123*x + Lstmp124*Lstmp6 - Lstmp125*y + Lstmp13*L[34] + Lstmp15 + Lstmp24 + Lstmp25*y + Lstmp28 + Lstmp36*Lstmp6 + Lstmp6*L[19] + Lstmp8*L[28] - Lstmp91 - Lstmp92 - Lstmp93*x + y*L[14] + L[7];
#pragma omp atomic
Ls[8] += -Lstmp101 - Lstmp102*x - Lstmp107 - Lstmp108*y - Lstmp109*Lstmp11 + Lstmp11*Lstmp97 + Lstmp11*L[24] - Lstmp110*x - Lstmp114*Lstmp6 + Lstmp117 + Lstmp126*x + Lstmp127*Lstmp6 - Lstmp128*x - Lstmp129*y + Lstmp13*L[35] + Lstmp6*L[20] + Lstmp8*L[29] + x*L[13] + y*L[15] + L[8];
#pragma omp atomic
Ls[9] += Lstmp11*L[28] - Lstmp119 + Lstmp20 + Lstmp21 + Lstmp22*x + Lstmp23*x + Lstmp34 + Lstmp6*L[25] + x*L[16] + L[9];
#pragma omp atomic
Ls[10] += Lstmp11*L[30] + Lstmp31 + Lstmp33*x + Lstmp6*L[26] - Lstmp84 + Lstmp88 + Lstmp89*x + Lstmp90 + x*L[17] + L[10];
#pragma omp atomic
Ls[11] += -Lstmp105 - Lstmp106*x + Lstmp11*L[31] - Lstmp115 - Lstmp121 + Lstmp6*L[27] + Lstmp98 + Lstmp99*x + x*L[18] + L[11];
#pragma omp atomic
Ls[12] += Lstmp11*L[32] + Lstmp122 + Lstmp123 + Lstmp124*x + Lstmp27 + Lstmp6*L[28] + Lstmp79 + Lstmp81 - Lstmp93 + L[12];
#pragma omp atomic
Ls[13] += -Lstmp102 + Lstmp11*L[33] - Lstmp110 - Lstmp114*x + Lstmp126 + Lstmp127*x - Lstmp128 + Lstmp6*L[29] + x*L[20] + L[13];
#pragma omp atomic
Ls[14] += Lstmp11*L[34] - Lstmp125 + Lstmp130*x + Lstmp16 + Lstmp25 + Lstmp26*y + Lstmp30 + Lstmp6*L[30] + y*L[23] + L[14];
#pragma omp atomic
Ls[15] += -Lstmp108 - Lstmp109*y + Lstmp11*L[35] - Lstmp113 - Lstmp129 + Lstmp131*x + Lstmp6*L[31] + Lstmp96 + y*L[24] + L[15];
#pragma omp atomic
Ls[16] += Lstmp22 + Lstmp23 + x*L[25] + L[16];
#pragma omp atomic
Ls[17] += Lstmp33 + Lstmp89 + x*L[26] + L[17];
#pragma omp atomic
Ls[18] += -Lstmp106 + Lstmp99 + x*L[27] + L[18];
#pragma omp atomic
Ls[19] += Lstmp118 + Lstmp124 + Lstmp36 + L[19];
#pragma omp atomic
Ls[20] += -Lstmp114 + Lstmp127 + x*L[29] + L[20];
#pragma omp atomic
Ls[21] += Lstmp130 + Lstmp29 + Lstmp80 + L[21];
#pragma omp atomic
Ls[22] += -Lstmp112 + Lstmp120 + Lstmp131 + L[22];
#pragma omp atomic
Ls[23] += Lstmp17 + Lstmp26 + y*L[34] + L[23];
#pragma omp atomic
Ls[24] += -Lstmp109 + Lstmp97 + y*L[35] + L[24];
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
#pragma omp atomic
Ls[35] += L[35];

}

void L2Pc_6(double x, double y, double z, double * L, double * F) {
double Ftmp0 = x*y;
double Ftmp1 = x*z;
double Ftmp2 = y*z;
double Ftmp3 = Ftmp0*z;
double Ftmp4 = pow(x, 2);
double Ftmp5 = (1.0/2.0)*Ftmp4;
double Ftmp6 = pow(x, 3);
double Ftmp7 = (1.0/6.0)*Ftmp6;
double Ftmp8 = (1.0/24.0)*pow(x, 4);
double Ftmp9 = pow(y, 2);
double Ftmp10 = (1.0/2.0)*Ftmp9;
double Ftmp11 = pow(y, 3);
double Ftmp12 = (1.0/6.0)*Ftmp11;
double Ftmp13 = (1.0/24.0)*pow(y, 4);
double Ftmp14 = Ftmp10*x;
double Ftmp15 = Ftmp12*x;
double Ftmp16 = Ftmp5*y;
double Ftmp17 = Ftmp5*z;
double Ftmp18 = Ftmp7*y;
double Ftmp19 = Ftmp7*z;
double Ftmp20 = Ftmp10*z;
double Ftmp21 = Ftmp12*z;
double Ftmp22 = Ftmp1*Ftmp10;
double Ftmp23 = Ftmp2*Ftmp5;
double Ftmp24 = (1.0/4.0)*Ftmp4*Ftmp9;
double Ftmp25 = pow(z, 2);
double Ftmp26 = L[9] + L[12];
double Ftmp27 = pow(z, 3);
double Ftmp28 = L[18] + L[22];
double Ftmp29 = L[16] + L[19];
double Ftmp30 = L[27] + L[31];
double Ftmp31 = L[17] + L[21];
double Ftmp32 = L[29] + L[33];
double Ftmp33 = L[26] + L[30];
double Ftmp34 = L[25] + L[28];
double Ftmp35 = L[28] + L[32];
double Ftmp36 = L[25] + 2*L[28] + L[32];
double Ftmp37 = (1.0/24.0)*pow(z, 4);
double Ftmp38 = L[10] + L[14];
double Ftmp39 = L[20] + L[24];
double Ftmp40 = L[19] + L[23];
double Ftmp41 = L[31] + L[35];
double Ftmp42 = L[30] + L[34];
double Ftmp43 = L[26] + 2*L[30] + L[34];
double Ftmp44 = (1.0/6.0)*Ftmp27;
#pragma omp atomic
F[0] += -Ftmp0*L[10] - Ftmp1*L[11] - Ftmp10*L[12] - Ftmp12*L[21] - Ftmp13*L[32] - Ftmp14*L[19] - Ftmp15*L[30] - Ftmp16*L[17] - Ftmp17*L[18] - Ftmp18*L[26] - Ftmp19*L[27] - Ftmp2*L[13] - Ftmp20*L[22] - Ftmp21*L[33] - Ftmp22*L[31] - Ftmp23*L[29] - Ftmp24*L[28] + (1.0/2.0)*Ftmp25*Ftmp26 + (1.0/2.0)*Ftmp25*Ftmp29*x + (1.0/2.0)*Ftmp25*Ftmp31*y + (1.0/2.0)*Ftmp25*Ftmp33*x*y + (1.0/4.0)*Ftmp25*Ftmp34*Ftmp4 + (1.0/4.0)*Ftmp25*Ftmp35*Ftmp9 + (1.0/6.0)*Ftmp27*Ftmp28 + (1.0/6.0)*Ftmp27*Ftmp30*x + (1.0/6.0)*Ftmp27*Ftmp32*y - Ftmp3*L[20] - Ftmp36*Ftmp37 - Ftmp5*L[9] - Ftmp7*L[16] - Ftmp8*L[25] - x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -Ftmp0*L[12] - Ftmp1*L[13] - Ftmp10*L[14] - Ftmp12*L[23] - Ftmp13*L[34] - Ftmp14*L[21] - Ftmp15*L[32] - Ftmp16*L[19] - Ftmp17*L[20] - Ftmp18*L[28] - Ftmp19*L[29] - Ftmp2*L[15] - Ftmp20*L[24] - Ftmp21*L[35] - Ftmp22*L[33] - Ftmp23*L[31] - Ftmp24*L[30] + (1.0/2.0)*Ftmp25*Ftmp31*x + (1.0/4.0)*Ftmp25*Ftmp33*Ftmp4 + (1.0/2.0)*Ftmp25*Ftmp35*x*y + (1.0/2.0)*Ftmp25*Ftmp38 + (1.0/2.0)*Ftmp25*Ftmp40*y + (1.0/4.0)*Ftmp25*Ftmp42*Ftmp9 + (1.0/6.0)*Ftmp27*Ftmp32*x + (1.0/6.0)*Ftmp27*Ftmp39 + (1.0/6.0)*Ftmp27*Ftmp41*y - Ftmp3*L[22] - Ftmp37*Ftmp43 - Ftmp5*L[10] - Ftmp7*L[17] - Ftmp8*L[26] - x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -Ftmp0*L[13] - Ftmp10*L[15] + (1.0/6.0)*Ftmp11*Ftmp42*z - Ftmp12*L[24] - Ftmp13*L[35] - Ftmp14*L[22] - Ftmp15*L[33] - Ftmp16*L[20] - Ftmp18*L[29] - Ftmp24*L[31] + (1.0/2.0)*Ftmp25*Ftmp28*x + (1.0/4.0)*Ftmp25*Ftmp30*Ftmp4 + (1.0/2.0)*Ftmp25*Ftmp32*x*y + (1.0/2.0)*Ftmp25*Ftmp39*y + (1.0/4.0)*Ftmp25*Ftmp41*Ftmp9 + (1.0/2.0)*Ftmp25*(L[11] + L[15]) + Ftmp26*x*z + (1.0/2.0)*Ftmp29*Ftmp4*z + Ftmp31*x*y*z + (1.0/2.0)*Ftmp33*Ftmp4*y*z + (1.0/6.0)*Ftmp34*Ftmp6*z + (1.0/2.0)*Ftmp35*Ftmp9*x*z - Ftmp36*Ftmp44*x - Ftmp37*(L[27] + 2*L[31] + L[35]) + Ftmp38*y*z + (1.0/2.0)*Ftmp40*Ftmp9*z - Ftmp43*Ftmp44*y - Ftmp44*(L[16] + 2*L[19] + L[23]) - Ftmp5*L[11] - Ftmp7*L[18] - Ftmp8*L[27] - x*L[6] - y*L[8] + z*(L[4] + L[7]) - L[3];

}

void M2Pc_6(double x, double y, double z, double * M, double * F) {
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(Rinv, 3);
double Ftmp1 = pow(Rinv, 2);
double Ftmp2 = 3*Ftmp1;
double Ftmp3 = y*M[4];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = Ftmp2*x;
double Ftmp6 = Ftmp5*y;
double Ftmp7 = Ftmp4*M[2];
double Ftmp8 = pow(Rinv, 4);
double Ftmp9 = 15*Ftmp8;
double Ftmp10 = Ftmp9*y;
double Ftmp11 = Ftmp10*M[12];
double Ftmp12 = pow(x, 2);
double Ftmp13 = Ftmp1*Ftmp12;
double Ftmp14 = 3*Ftmp13;
double Ftmp15 = Ftmp10*M[7];
double Ftmp16 = x*z;
double Ftmp17 = Ftmp8*x;
double Ftmp18 = pow(y, 2);
double Ftmp19 = 6*Ftmp18*M[6];
double Ftmp20 = Ftmp12*Ftmp9;
double Ftmp21 = Ftmp9*z;
double Ftmp22 = pow(Rinv, 6);
double Ftmp23 = 30*x;
double Ftmp24 = Ftmp22*Ftmp23;
double Ftmp25 = pow(y, 3);
double Ftmp26 = Ftmp25*M[13];
double Ftmp27 = Ftmp18*z;
double Ftmp28 = Ftmp24*Ftmp27;
double Ftmp29 = 105*M[12];
double Ftmp30 = Ftmp22*y;
double Ftmp31 = Ftmp12*Ftmp30;
double Ftmp32 = Ftmp31*z;
double Ftmp33 = Ftmp25*M[23];
double Ftmp34 = pow(Rinv, 8);
double Ftmp35 = 210*Ftmp34;
double Ftmp36 = Ftmp16*Ftmp35;
double Ftmp37 = Ftmp12*Ftmp22;
double Ftmp38 = 30*Ftmp18;
double Ftmp39 = Ftmp37*Ftmp38;
double Ftmp40 = Ftmp12*Ftmp25;
double Ftmp41 = Ftmp35*M[20];
double Ftmp42 = Ftmp12*Ftmp18;
double Ftmp43 = Ftmp35*z;
double Ftmp44 = Ftmp42*Ftmp43;
double Ftmp45 = z*M[32];
double Ftmp46 = pow(Rinv, 10);
double Ftmp47 = Ftmp40*Ftmp46;
double Ftmp48 = 5*Ftmp13;
double Ftmp49 = (Ftmp48 - 3)*M[8];
double Ftmp50 = Ftmp1*Ftmp18;
double Ftmp51 = 5*Ftmp50;
double Ftmp52 = Ftmp51 - 1;
double Ftmp53 = Ftmp2*Ftmp52;
double Ftmp54 = 6*Ftmp1;
double Ftmp55 = Ftmp13 - 1;
double Ftmp56 = Ftmp55*x;
double Ftmp57 = Ftmp14 - 1;
double Ftmp58 = Ftmp57*M[3];
double Ftmp59 = 3*Ftmp50;
double Ftmp60 = Ftmp59 - 1;
double Ftmp61 = Ftmp60*M[6];
double Ftmp62 = 7*Ftmp13;
double Ftmp63 = Ftmp62 - 3;
double Ftmp64 = Ftmp63*M[16];
double Ftmp65 = 7*Ftmp50;
double Ftmp66 = Ftmp65 - 3;
double Ftmp67 = Ftmp66*M[20];
double Ftmp68 = Ftmp63*M[17];
double Ftmp69 = Ftmp65 - 1;
double Ftmp70 = Ftmp69*M[21];
double Ftmp71 = Ftmp8*y;
double Ftmp72 = Ftmp23*Ftmp55;
double Ftmp73 = Ftmp8*z;
double Ftmp74 = Ftmp12*Ftmp8;
double Ftmp75 = 30*M[8];
double Ftmp76 = Ftmp48 - 1;
double Ftmp77 = Ftmp76*M[9];
double Ftmp78 = Ftmp10*x;
double Ftmp79 = (Ftmp51 - 3)*M[13];
double Ftmp80 = Ftmp21*x;
double Ftmp81 = Ftmp76*M[10];
double Ftmp82 = Ftmp52*M[14];
double Ftmp83 = 315*Ftmp30;
double Ftmp84 = Ftmp83*z;
double Ftmp85 = Ftmp57*M[28];
double Ftmp86 = Ftmp60*Ftmp83;
double Ftmp87 = Ftmp30*z;
double Ftmp88 = Ftmp52*M[11];
double Ftmp89 = 210*Ftmp37;
double Ftmp90 = Ftmp62 - 1;
double Ftmp91 = Ftmp90*M[19];
double Ftmp92 = Ftmp30*x;
double Ftmp93 = 105*Ftmp92;
double Ftmp94 = Ftmp93*z;
double Ftmp95 = Ftmp66*M[23];
double Ftmp96 = Ftmp22*x;
double Ftmp97 = Ftmp18*Ftmp96;
double Ftmp98 = 60*M[22];
double Ftmp99 = Ftmp66*Ftmp98;
double Ftmp100 = 105*Ftmp31;
double Ftmp101 = 105*Ftmp37;
double Ftmp102 = Ftmp101*z;
double Ftmp103 = Ftmp34*x;
double Ftmp104 = 9*Ftmp50;
double Ftmp105 = 420*M[33];
double Ftmp106 = Ftmp105*Ftmp25*(Ftmp104 - 5);
double Ftmp107 = Ftmp34*z;
double Ftmp108 = 1890*Ftmp12;
double Ftmp109 = Ftmp107*Ftmp108;
double Ftmp110 = 1260*M[34];
double Ftmp111 = Ftmp107*Ftmp18;
double Ftmp112 = Ftmp111*Ftmp60*x;
double Ftmp113 = 2835*Ftmp12;
double Ftmp114 = Ftmp107*Ftmp113;
double Ftmp115 = Ftmp114*y;
double Ftmp116 = Ftmp60*M[32];
double Ftmp117 = 11*Ftmp50;
double Ftmp118 = Ftmp117 - 5;
double Ftmp119 = 3780*M[47];
double Ftmp120 = Ftmp118*Ftmp119;
double Ftmp121 = Ftmp16*Ftmp25;
double Ftmp122 = Ftmp121*Ftmp46;
double Ftmp123 = 1260*M[31];
double Ftmp124 = Ftmp34*Ftmp60;
double Ftmp125 = 3780*M[44];
double Ftmp126 = Ftmp118*Ftmp125;
double Ftmp127 = Ftmp42*Ftmp46*z;
double Ftmp128 = 3780*(Ftmp117 - 3)*M[45];
double Ftmp129 = pow(x, 4);
double Ftmp130 = Ftmp129*Ftmp8;
double Ftmp131 = 63*Ftmp130;
double Ftmp132 = (-70*Ftmp13 + Ftmp131 + 15)*M[24];
double Ftmp133 = 45*Ftmp8;
double Ftmp134 = pow(y, 4);
double Ftmp135 = Ftmp134*Ftmp8;
double Ftmp136 = 21*Ftmp135;
double Ftmp137 = 14*Ftmp50;
double Ftmp138 = -Ftmp137;
double Ftmp139 = Ftmp138 + 1;
double Ftmp140 = Ftmp136 + Ftmp139;
double Ftmp141 = Ftmp140*M[31];
double Ftmp142 = -10*Ftmp13;
double Ftmp143 = Ftmp142 + 3;
double Ftmp144 = 60*M[15];
double Ftmp145 = Ftmp9*x;
double Ftmp146 = -30*Ftmp13;
double Ftmp147 = (35*Ftmp130 + Ftmp146 + 3)*M[15];
double Ftmp148 = -30*Ftmp50;
double Ftmp149 = (35*Ftmp135 + Ftmp148 + 3)*M[22];
double Ftmp150 = 33*Ftmp130;
double Ftmp151 = Ftmp146 + Ftmp150 + 5;
double Ftmp152 = Ftmp151*M[36];
double Ftmp153 = 33*Ftmp135 + Ftmp148 + 5;
double Ftmp154 = Ftmp153*M[44];
double Ftmp155 = 315*Ftmp22;
double Ftmp156 = z*M[37];
double Ftmp157 = Ftmp151*Ftmp156;
double Ftmp158 = -5040*Ftmp13 + 3780*Ftmp130 + 1260;
double Ftmp159 = x*M[25];
double Ftmp160 = 21*Ftmp130;
double Ftmp161 = 14*Ftmp13;
double Ftmp162 = -Ftmp161;
double Ftmp163 = Ftmp162 + 1;
double Ftmp164 = Ftmp160 + Ftmp163;
double Ftmp165 = Ftmp164*Ftmp83;
double Ftmp166 = 63*Ftmp135;
double Ftmp167 = (Ftmp166 - 70*Ftmp50 + 15)*M[33];
double Ftmp168 = Ftmp96*z;
double Ftmp169 = Ftmp155*z;
double Ftmp170 = Ftmp169*x;
double Ftmp171 = Ftmp140*M[34];
double Ftmp172 = 420*M[24];
double Ftmp173 = 11*Ftmp130;
double Ftmp174 = Ftmp162 + 3;
double Ftmp175 = Ftmp107*x;
double Ftmp176 = Ftmp175*y;
double Ftmp177 = 3780*M[39];
double Ftmp178 = 18*Ftmp13;
double Ftmp179 = (Ftmp150 - Ftmp178 + 1)*M[39];
double Ftmp180 = 2835*Ftmp176;
double Ftmp181 = Ftmp153*M[47];
double Ftmp182 = Ftmp103*Ftmp18;
double Ftmp183 = 1890*M[46];
double Ftmp184 = Ftmp153*Ftmp183;
double Ftmp185 = -16*Ftmp13 + Ftmp173 + 5;
double Ftmp186 = Ftmp34*y;
double Ftmp187 = Ftmp12*Ftmp186;
double Ftmp188 = 3780*M[36];
double Ftmp189 = Ftmp113*Ftmp186;
double Ftmp190 = Ftmp12*Ftmp34;
double Ftmp191 = 3780*Ftmp156;
double Ftmp192 = pow(x, 6);
double Ftmp193 = 33*Ftmp22;
double Ftmp194 = 1890*M[35];
double Ftmp195 = Ftmp155*x;
double Ftmp196 = 231*Ftmp22;
double Ftmp197 = (105*Ftmp13 - 315*Ftmp130 + Ftmp192*Ftmp196 - 5)*M[35];
double Ftmp198 = pow(y, 6);
double Ftmp199 = (-315*Ftmp135 + Ftmp196*Ftmp198 + 105*Ftmp50 - 5)*M[46];
double Ftmp200 = -Ftmp62;
double Ftmp201 = Ftmp18*Ftmp74;
double Ftmp202 = 63*Ftmp201;
double Ftmp203 = Ftmp202 + 3;
double Ftmp204 = (Ftmp200 + Ftmp203 - 21*Ftmp50)*M[27];
double Ftmp205 = 8*Ftmp50;
double Ftmp206 = 14*Ftmp201;
double Ftmp207 = -Ftmp13;
double Ftmp208 = Ftmp207 + 1;
double Ftmp209 = Ftmp8*M[18];
double Ftmp210 = (35*Ftmp201 - Ftmp48 - Ftmp51 + 1)*M[18];
double Ftmp211 = -Ftmp59;
double Ftmp212 = -Ftmp14;
double Ftmp213 = Ftmp212 + 1;
double Ftmp214 = (11*Ftmp201 + Ftmp211 + Ftmp213)*M[40];
double Ftmp215 = 945*Ftmp214;
double Ftmp216 = 33*Ftmp201;
double Ftmp217 = (-Ftmp104 + Ftmp213 + Ftmp216)*M[41];
double Ftmp218 = 18*Ftmp201;
double Ftmp219 = -10*Ftmp50;
double Ftmp220 = Ftmp219 + 3;
double Ftmp221 = 210*M[29];
double Ftmp222 = -Ftmp65;
double Ftmp223 = (-21*Ftmp13 + Ftmp203 + Ftmp222)*M[29];
double Ftmp224 = Ftmp208 + Ftmp218;
double Ftmp225 = 210*M[30];
double Ftmp226 = (Ftmp200 + Ftmp202 + Ftmp222 + 1)*M[30];
double Ftmp227 = 105*Ftmp226;
double Ftmp228 = -12*Ftmp50;
double Ftmp229 = 22*Ftmp201;
double Ftmp230 = Ftmp229 + 3;
double Ftmp231 = Ftmp212 + Ftmp230;
double Ftmp232 = 1890*M[43];
double Ftmp233 = Ftmp232*y;
double Ftmp234 = 9*Ftmp13;
double Ftmp235 = (Ftmp211 + Ftmp216 - Ftmp234 + 1)*M[43];
double Ftmp236 = 8505*Ftmp214;
double Ftmp237 = 4*Ftmp50;
double Ftmp238 = (Ftmp140 + Ftmp237*Ftmp60)*M[45];
double Ftmp239 = -36*Ftmp201;
double Ftmp240 = Ftmp12*Ftmp134;
double Ftmp241 = 99*Ftmp22;
double Ftmp242 = Ftmp240*Ftmp241;
double Ftmp243 = 630*M[42];
double Ftmp244 = -126*Ftmp201;
double Ftmp245 = (-Ftmp136 + Ftmp137 + Ftmp196*Ftmp240 + Ftmp244 + Ftmp90)*M[42];
double Ftmp246 = 8*Ftmp13;
double Ftmp247 = Ftmp129*Ftmp18;
double Ftmp248 = Ftmp241*Ftmp247;
double Ftmp249 = -102*Ftmp201 - 2;
double Ftmp250 = 630*M[38];
double Ftmp251 = (-Ftmp160 + Ftmp161 + Ftmp196*Ftmp247 + Ftmp244 + Ftmp69)*M[38];
double Ftmp252 = 6*Ftmp74*M[3];
double Ftmp253 = pow(x, 3);
double Ftmp254 = Ftmp253*Ftmp75;
double Ftmp255 = 30*Ftmp32;
double Ftmp256 = Ftmp253*Ftmp35;
double Ftmp257 = Ftmp256*M[17];
double Ftmp258 = y*z;
double Ftmp259 = Ftmp256*M[16];
double Ftmp260 = Ftmp253*Ftmp46;
double Ftmp261 = 1890*Ftmp260*M[28];
double Ftmp262 = Ftmp2*Ftmp76;
double Ftmp263 = Ftmp50 - 1;
double Ftmp264 = Ftmp2*y;
double Ftmp265 = Ftmp263*Ftmp8;
double Ftmp266 = 30*y;
double Ftmp267 = Ftmp10*z;
double Ftmp268 = 210*Ftmp92;
double Ftmp269 = Ftmp18*Ftmp9;
double Ftmp270 = 210*Ftmp263;
double Ftmp271 = Ftmp18*Ftmp22;
double Ftmp272 = Ftmp271*z;
double Ftmp273 = 105*Ftmp97;
double Ftmp274 = Ftmp144*Ftmp63;
double Ftmp275 = 105*Ftmp271;
double Ftmp276 = Ftmp275*z;
double Ftmp277 = Ftmp172*Ftmp253*(Ftmp234 - 5);
double Ftmp278 = x*M[32];
double Ftmp279 = 1890*Ftmp278;
double Ftmp280 = 2835*Ftmp111;
double Ftmp281 = Ftmp280*x;
double Ftmp282 = 1260*Ftmp57;
double Ftmp283 = Ftmp12*y;
double Ftmp284 = Ftmp107*Ftmp282*Ftmp283;
double Ftmp285 = 11*Ftmp13;
double Ftmp286 = Ftmp260*(Ftmp285 - 5);
double Ftmp287 = Ftmp18*M[25];
double Ftmp288 = Ftmp190*Ftmp282;
double Ftmp289 = Ftmp188*Ftmp286;
double Ftmp290 = Ftmp177*(Ftmp285 - 3);
double Ftmp291 = Ftmp133*Ftmp164;
double Ftmp292 = 3*Ftmp135 - Ftmp237 + 1;
double Ftmp293 = Ftmp83*x;
double Ftmp294 = Ftmp165*z;
double Ftmp295 = Ftmp155*Ftmp164;
double Ftmp296 = 2835*Ftmp103;
double Ftmp297 = 2835*Ftmp182;
double Ftmp298 = 11*Ftmp135 - 16*Ftmp50 + 5;
double Ftmp299 = Ftmp151*Ftmp194;
double Ftmp300 = -Ftmp50;
double Ftmp301 = Ftmp300 + 1;
double Ftmp302 = Ftmp218 + Ftmp301;
double Ftmp303 = -12*Ftmp13;
double Ftmp304 = 1890*M[41];
double Ftmp305 = 1890*M[40];
double Ftmp306 = pow(z, 2);
double Ftmp307 = Ftmp22*z;
double Ftmp308 = 30*Ftmp307;
double Ftmp309 = Ftmp306*M[10];
double Ftmp310 = Ftmp22*Ftmp306;
double Ftmp311 = Ftmp306*Ftmp35;
double Ftmp312 = Ftmp306*Ftmp46;
double Ftmp313 = Ftmp25*Ftmp312;
double Ftmp314 = Ftmp306*y;
double Ftmp315 = 105*Ftmp306;
double Ftmp316 = Ftmp315*Ftmp96;
double Ftmp317 = Ftmp30*Ftmp315;
double Ftmp318 = Ftmp306*M[26];
double Ftmp319 = Ftmp312*x;
double Ftmp320 = Ftmp151*M[37];
double Ftmp321 = Ftmp12 + Ftmp18;
double Ftmp322 = 105*Ftmp168;
double Ftmp323 = 105*Ftmp87;
double Ftmp324 = -Ftmp178*Ftmp18;
double Ftmp325 = 3*Ftmp18;
double Ftmp326 = Ftmp12 + Ftmp325;
double Ftmp327 = 2835*Ftmp306;
double Ftmp328 = Ftmp103*Ftmp327;
double Ftmp329 = 3*Ftmp12;
double Ftmp330 = Ftmp18 + Ftmp329;
double Ftmp331 = Ftmp186*Ftmp327;
double Ftmp332 = Ftmp13*Ftmp18;
double Ftmp333 = -22*Ftmp332;
double Ftmp334 = -36*Ftmp332;
#pragma omp atomic
F[0] += Ftmp0*(-Ftmp10*Ftmp64 - Ftmp10*Ftmp67 + Ftmp100*Ftmp64 + Ftmp100*Ftmp67 - Ftmp101*Ftmp132 - Ftmp101*Ftmp204 + Ftmp102*Ftmp68 + Ftmp102*Ftmp70 - Ftmp103*Ftmp106 + Ftmp108*Ftmp186*(Ftmp138 + Ftmp231)*M[40] - Ftmp109*Ftmp55*y*M[28] + Ftmp109*(Ftmp139 + Ftmp207 + Ftmp229)*M[41] + Ftmp11*z - Ftmp110*Ftmp112 + Ftmp113*Ftmp157*Ftmp34 + Ftmp114*Ftmp217 + Ftmp114*Ftmp238 - Ftmp115*Ftmp116 - Ftmp115*Ftmp85 + Ftmp12*Ftmp21*M[5] + Ftmp120*Ftmp122 - Ftmp123*Ftmp124*Ftmp42 + Ftmp126*Ftmp47 + Ftmp127*Ftmp128 + Ftmp132*Ftmp9 + Ftmp133*Ftmp141 - Ftmp14*M[0] - 315*Ftmp141*Ftmp37 + Ftmp144*Ftmp17*(7*Ftmp130 + Ftmp143) + Ftmp145*Ftmp147 + Ftmp145*Ftmp149 + Ftmp145*Ftmp210 + Ftmp15*Ftmp16 + Ftmp152*Ftmp189 - Ftmp152*Ftmp83 + Ftmp154*Ftmp189 - Ftmp154*Ftmp83 - Ftmp155*Ftmp157 - Ftmp158*Ftmp159*Ftmp30 - Ftmp158*Ftmp168*M[26] - Ftmp159*Ftmp165 - Ftmp164*Ftmp170*M[26] - Ftmp167*Ftmp93 - Ftmp168*Ftmp225*(Ftmp219 + Ftmp224) - Ftmp168*Ftmp227 - Ftmp169*Ftmp217 - Ftmp169*Ftmp238 + Ftmp17*Ftmp19 - Ftmp170*Ftmp171 - Ftmp172*Ftmp37*(9*Ftmp130 + Ftmp162 + 5) + Ftmp175*Ftmp233*(Ftmp228 + Ftmp231) + Ftmp176*Ftmp177*(Ftmp173 + Ftmp174) + Ftmp179*Ftmp180 + Ftmp180*Ftmp181 + Ftmp180*Ftmp235 + Ftmp182*Ftmp184 + Ftmp185*Ftmp187*Ftmp188 + Ftmp185*Ftmp190*Ftmp191 + Ftmp187*Ftmp236 + Ftmp194*Ftmp96*(35*Ftmp13 - Ftmp131 + Ftmp192*Ftmp193 - 5) + Ftmp195*Ftmp197 + Ftmp195*Ftmp199 + Ftmp195*Ftmp245 + Ftmp195*Ftmp251 - Ftmp2*Ftmp3 + Ftmp2*Ftmp49 + Ftmp20*Ftmp3 - Ftmp20*Ftmp49 - Ftmp20*Ftmp88 + Ftmp204*Ftmp9 + Ftmp209*Ftmp23*(-Ftmp205 + Ftmp206 + Ftmp208) - Ftmp21*Ftmp68 - Ftmp21*Ftmp70 - Ftmp215*Ftmp30 - Ftmp221*Ftmp92*(Ftmp212 + Ftmp218 + Ftmp220) - Ftmp223*Ftmp93 - Ftmp24*Ftmp26 + Ftmp243*Ftmp96*(-39*Ftmp135 + Ftmp239 + Ftmp242 + 20*Ftmp50 + Ftmp55) + Ftmp250*Ftmp96*(-6*Ftmp130 + Ftmp246 + Ftmp248 + Ftmp249 + 19*Ftmp50) - Ftmp28*M[14] - Ftmp29*Ftmp32 + 210*Ftmp31*Ftmp55*M[16] + Ftmp33*Ftmp36 - Ftmp39*M[11] - Ftmp4*M[5] + Ftmp40*Ftmp41 + Ftmp44*M[21] - 1890*Ftmp45*Ftmp47 + Ftmp45*Ftmp86 + Ftmp5*Ftmp58 + Ftmp5*Ftmp61 + Ftmp53*M[11] + Ftmp54*Ftmp56*M[3] - Ftmp55*Ftmp74*Ftmp75 + Ftmp55*Ftmp89*z*M[17] + 210*Ftmp56*Ftmp87*M[19] - Ftmp6*M[1] - Ftmp7*x - Ftmp71*Ftmp72*M[9] - Ftmp72*Ftmp73*M[10] - Ftmp77*Ftmp78 - Ftmp78*Ftmp79 - Ftmp80*Ftmp81 - Ftmp80*Ftmp82 + Ftmp84*Ftmp85 - Ftmp89*(Ftmp224 + Ftmp228)*M[27] + Ftmp91*Ftmp94 + Ftmp94*Ftmp95 + Ftmp97*Ftmp99 + M[0]);
#pragma omp atomic
F[1] += Ftmp0*(Ftmp10*Ftmp147 + Ftmp10*Ftmp149 + Ftmp10*Ftmp210 - Ftmp105*Ftmp271*(9*Ftmp135 + Ftmp138 + 5) - Ftmp110*Ftmp292*Ftmp87 + Ftmp111*Ftmp119*Ftmp298 + Ftmp111*Ftmp232*(Ftmp163 + Ftmp229 + Ftmp300) - Ftmp111*Ftmp263*Ftmp279 + Ftmp116*Ftmp170 - Ftmp116*Ftmp281 - Ftmp123*Ftmp292*Ftmp92 + Ftmp125*Ftmp182*Ftmp298 + Ftmp127*Ftmp290 - Ftmp132*Ftmp93 - Ftmp141*Ftmp293 + Ftmp145*Ftmp18*M[4] - Ftmp145*Ftmp64 - Ftmp145*Ftmp67 - Ftmp152*Ftmp195 + Ftmp152*Ftmp297 - Ftmp154*Ftmp195 + Ftmp154*Ftmp297 + Ftmp157*Ftmp296*y - Ftmp167*Ftmp275 + Ftmp167*Ftmp9 - Ftmp168*Ftmp18*Ftmp29 - Ftmp169*Ftmp179 - Ftmp169*Ftmp181 - Ftmp169*Ftmp235 + Ftmp170*Ftmp85 - Ftmp171*Ftmp84 + Ftmp176*Ftmp304*(Ftmp211 + Ftmp230 + Ftmp303) + 1260*Ftmp176*(Ftmp136 + 6*Ftmp263*Ftmp50 + 2*Ftmp50*Ftmp60 - 34*Ftmp50 + 9)*M[45] + Ftmp179*Ftmp280 + Ftmp18*Ftmp21*M[7] + Ftmp18*Ftmp259 + Ftmp18*Ftmp289 + Ftmp180*Ftmp217 + Ftmp180*Ftmp238 + Ftmp181*Ftmp280 + Ftmp182*Ftmp236 + Ftmp182*Ftmp305*(Ftmp174 + Ftmp211 + Ftmp229) + Ftmp183*Ftmp30*(-Ftmp166 + Ftmp193*Ftmp198 + 35*Ftmp50 - 5) - Ftmp186*Ftmp277 + Ftmp187*Ftmp299 + Ftmp191*Ftmp286*y + Ftmp197*Ftmp83 + Ftmp199*Ftmp83 + Ftmp2*Ftmp79 - Ftmp204*Ftmp93 + Ftmp209*Ftmp266*(Ftmp206 - Ftmp246 + Ftmp301) - Ftmp21*Ftmp91 - Ftmp21*Ftmp95 - Ftmp215*Ftmp96 - Ftmp221*Ftmp271*(Ftmp302 + Ftmp303) - Ftmp223*Ftmp275 + Ftmp223*Ftmp9 - Ftmp225*Ftmp87*(Ftmp142 + Ftmp302) - Ftmp227*Ftmp87 - Ftmp23*Ftmp265*y*M[11] + Ftmp235*Ftmp280 + Ftmp243*Ftmp30*(19*Ftmp13 - 6*Ftmp135 + Ftmp205 + Ftmp242 + Ftmp249) + Ftmp245*Ftmp83 + Ftmp250*Ftmp30*(20*Ftmp13 - 39*Ftmp130 + Ftmp239 + Ftmp248 + Ftmp263) + Ftmp251*Ftmp83 + Ftmp252*y - Ftmp254*Ftmp30 - Ftmp255*M[10] + Ftmp257*Ftmp258 - Ftmp261*Ftmp27 + Ftmp262*M[9] + Ftmp263*Ftmp268*z*M[21] + Ftmp263*Ftmp54*y*M[6] + Ftmp264*Ftmp58 + Ftmp264*Ftmp61 - Ftmp265*Ftmp266*z*M[14] - Ftmp265*Ftmp38*M[13] - Ftmp267*Ftmp81 - Ftmp267*Ftmp82 - Ftmp268*(Ftmp143 + Ftmp211 + Ftmp218)*M[27] - Ftmp269*Ftmp77 - Ftmp269*Ftmp79 + Ftmp270*Ftmp272*M[23] + Ftmp270*Ftmp97*M[20] + Ftmp273*Ftmp64 + Ftmp273*Ftmp67 + Ftmp274*Ftmp31 + Ftmp276*Ftmp91 + Ftmp276*Ftmp95 - Ftmp281*Ftmp85 - Ftmp284*M[26] - Ftmp287*Ftmp288 - Ftmp287*Ftmp295 + Ftmp291*M[25] - Ftmp294*M[26] - Ftmp39*M[9] - Ftmp4*M[7] + Ftmp44*M[19] - Ftmp49*Ftmp78 - Ftmp5*M[4] - Ftmp59*M[1] - Ftmp6*M[0] + Ftmp68*Ftmp94 - Ftmp7*y + Ftmp70*Ftmp94 + Ftmp71*Ftmp98*(7*Ftmp135 + Ftmp220) - Ftmp78*Ftmp88 + Ftmp78*z*M[5] + Ftmp80*M[12] + M[1]);
#pragma omp atomic
F[2] += Ftmp0*(-Ftmp10*Ftmp91 - Ftmp10*Ftmp95 - Ftmp106*Ftmp107 + Ftmp107*Ftmp12*Ftmp299 + Ftmp107*Ftmp243*(99*Ftmp12*Ftmp135 + Ftmp12 - Ftmp134*Ftmp54 + 2*Ftmp18 + Ftmp334) + Ftmp107*Ftmp250*(2*Ftmp12 - Ftmp129*Ftmp54 + 99*Ftmp130*Ftmp18 + Ftmp18 + Ftmp334) - Ftmp107*Ftmp277 + Ftmp11*x - Ftmp110*Ftmp124*Ftmp18*Ftmp306 + Ftmp111*Ftmp184 - Ftmp112*Ftmp123 + Ftmp120*Ftmp313 + Ftmp121*Ftmp41 + Ftmp122*Ftmp126 - 2835*Ftmp124*Ftmp278*Ftmp314 + Ftmp128*Ftmp18*Ftmp319 - Ftmp132*Ftmp322 + Ftmp133*Ftmp171 - Ftmp141*Ftmp170 + Ftmp145*Ftmp306*M[5] - Ftmp145*Ftmp68 - Ftmp145*Ftmp70 + Ftmp147*Ftmp21 + Ftmp149*Ftmp21 + Ftmp15*Ftmp306 + Ftmp152*Ftmp180 + Ftmp154*Ftmp180 - Ftmp155*Ftmp171*Ftmp306 - Ftmp16*Ftmp305*Ftmp46*y*(Ftmp325 + Ftmp329 + Ftmp333) - Ftmp167*Ftmp323 + Ftmp169*Ftmp197 + Ftmp169*Ftmp199 + Ftmp169*Ftmp245 + Ftmp169*Ftmp251 + Ftmp176*Ftmp236 + Ftmp179*Ftmp331 - Ftmp179*Ftmp83 + Ftmp18*Ftmp311*x*M[21] + Ftmp181*Ftmp331 - Ftmp181*Ftmp83 + Ftmp19*Ftmp73 - Ftmp195*Ftmp217 - Ftmp195*Ftmp238 - Ftmp195*Ftmp320 - Ftmp2*Ftmp306*M[2] - Ftmp204*Ftmp322 + Ftmp21*Ftmp210 + Ftmp217*Ftmp328 - Ftmp223*Ftmp323 + Ftmp226*Ftmp9 - Ftmp227*Ftmp310 - Ftmp233*Ftmp312*(Ftmp330 + Ftmp333) + Ftmp235*Ftmp331 - Ftmp235*Ftmp83 + Ftmp238*Ftmp328 + Ftmp252*z - Ftmp254*Ftmp307 - Ftmp255*M[9] + Ftmp257*Ftmp306 + Ftmp258*Ftmp259 + Ftmp258*Ftmp289 - Ftmp26*Ftmp308 - Ftmp261*Ftmp314 + Ftmp262*M[10] - Ftmp264*M[7] - Ftmp267*Ftmp77 - Ftmp267*Ftmp79 + Ftmp272*Ftmp99 + Ftmp274*Ftmp37*z + Ftmp278*Ftmp86 - Ftmp279*Ftmp313 - Ftmp28*M[11] + Ftmp283*Ftmp290*Ftmp312 + Ftmp283*Ftmp311*M[19] - Ftmp284*M[25] + 3780*Ftmp286*Ftmp306*M[37] - Ftmp288*Ftmp318 + Ftmp291*M[26] + Ftmp293*Ftmp85 - Ftmp294*M[25] - Ftmp295*Ftmp318 - Ftmp296*Ftmp314*Ftmp85 + Ftmp3*Ftmp80 - Ftmp304*Ftmp319*(Ftmp326 + Ftmp333) - Ftmp306*Ftmp82*Ftmp9 - Ftmp306*Ftmp93*M[12] - Ftmp308*(-Ftmp161*Ftmp18 + Ftmp321)*M[18] - 30*Ftmp309*Ftmp37 - Ftmp309*Ftmp76*Ftmp9 - Ftmp310*Ftmp38*M[14] + Ftmp311*Ftmp33 + Ftmp311*(Ftmp321 + Ftmp324)*M[30] + Ftmp316*Ftmp68 + Ftmp316*Ftmp70 + Ftmp317*Ftmp91 + Ftmp317*Ftmp95 + Ftmp320*Ftmp328 + Ftmp36*(Ftmp324 + Ftmp326)*M[27] + Ftmp4*Ftmp58 + Ftmp4*Ftmp61 - Ftmp4*x*M[0] - Ftmp4*y*M[1] + Ftmp43*y*(Ftmp324 + Ftmp330)*M[29] - Ftmp49*Ftmp80 - Ftmp5*M[5] + Ftmp53*M[14] + Ftmp64*Ftmp94 + Ftmp67*Ftmp94 - Ftmp80*Ftmp88 + M[2]);

}

void M2Lc_6(double x, double y, double z, double * M, double * L) {
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double D[63];
double Dtmp0 = pow(Rinv, 3);
double Dtmp1 = pow(x, 2);
double Dtmp2 = pow(Rinv, 2);
double Dtmp3 = 3*Dtmp2;
double Dtmp4 = Dtmp1*Dtmp3;
double Dtmp5 = Dtmp4 - 1;
double Dtmp6 = 3*pow(Rinv, 5);
double Dtmp7 = Dtmp6*x;
double Dtmp8 = pow(y, 2);
double Dtmp9 = Dtmp3*Dtmp8;
double Dtmp10 = Dtmp9 - 1;
double Dtmp11 = Dtmp6*y;
double Dtmp12 = 5*Dtmp2;
double Dtmp13 = Dtmp1*Dtmp12;
double Dtmp14 = Dtmp13 - 1;
double Dtmp15 = Dtmp6*z;
double Dtmp16 = Dtmp12*Dtmp8;
double Dtmp17 = Dtmp16 - 1;
double Dtmp18 = pow(Rinv, 7);
double Dtmp19 = 15*Dtmp18;
double Dtmp20 = Dtmp19*x;
double Dtmp21 = Dtmp20*y;
double Dtmp22 = Dtmp1*Dtmp2;
double Dtmp23 = -30*Dtmp22;
double Dtmp24 = pow(x, 4);
double Dtmp25 = pow(Rinv, 4);
double Dtmp26 = 35*Dtmp25;
double Dtmp27 = 7*Dtmp22;
double Dtmp28 = Dtmp27 - 3;
double Dtmp29 = Dtmp20*z;
double Dtmp30 = Dtmp1*Dtmp8;
double Dtmp31 = Dtmp27 - 1;
double Dtmp32 = Dtmp19*y;
double Dtmp33 = Dtmp32*z;
double Dtmp34 = Dtmp2*Dtmp8;
double Dtmp35 = 7*Dtmp34;
double Dtmp36 = Dtmp35 - 3;
double Dtmp37 = Dtmp35 - 1;
double Dtmp38 = -30*Dtmp34;
double Dtmp39 = pow(y, 4);
double Dtmp40 = Dtmp24*Dtmp25;
double Dtmp41 = 14*Dtmp22;
double Dtmp42 = 21*Dtmp40;
double Dtmp43 = 45*Dtmp18;
double Dtmp44 = Dtmp43*(-Dtmp41 + Dtmp42 + 1);
double Dtmp45 = -Dtmp27;
double Dtmp46 = Dtmp25*Dtmp30;
double Dtmp47 = 63*Dtmp46;
double Dtmp48 = Dtmp47 + 3;
double Dtmp49 = pow(Rinv, 9);
double Dtmp50 = 315*Dtmp49*y;
double Dtmp51 = Dtmp50*z;
double Dtmp52 = Dtmp51*x;
double Dtmp53 = -Dtmp35;
double Dtmp54 = 14*Dtmp34;
double Dtmp55 = Dtmp25*Dtmp39;
double Dtmp56 = 21*Dtmp55;
double Dtmp57 = -Dtmp54 + Dtmp56 + 1;
double Dtmp58 = Dtmp43*Dtmp57;
double Dtmp59 = 231*pow(Rinv, 6);
double Dtmp60 = 33*Dtmp40;
double Dtmp61 = x*(Dtmp23 + Dtmp60 + 5);
double Dtmp62 = 315*Dtmp49*z;
double Dtmp63 = -126*Dtmp46;
double Dtmp64 = -Dtmp9;
double Dtmp65 = 1 - Dtmp4;
double Dtmp66 = 33*Dtmp46;
double Dtmp67 = Dtmp62*x;
double Dtmp68 = Dtmp38 + 33*Dtmp55 + 5;
D[0] = -Dtmp0*x;
D[1] = -Dtmp0*y;
D[2] = -Dtmp0*z;
D[3] = Dtmp0*Dtmp5;
D[4] = Dtmp7*y;
D[5] = Dtmp7*z;
D[6] = Dtmp0*Dtmp10;
D[7] = Dtmp11*z;
D[8] = -D[3] - D[6];
D[9] = -Dtmp7*(Dtmp13 - 3);
D[10] = -Dtmp11*Dtmp14;
D[11] = -Dtmp14*Dtmp15;
D[12] = -Dtmp17*Dtmp7;
D[13] = -Dtmp21*z;
D[14] = -D[9] - D[12];
D[15] = -Dtmp11*(Dtmp16 - 3);
D[16] = -Dtmp15*Dtmp17;
D[17] = -D[10] - D[15];
D[18] = Dtmp6*(Dtmp23 + Dtmp24*Dtmp26 + 3);
D[19] = Dtmp21*Dtmp28;
D[20] = Dtmp28*Dtmp29;
D[21] = Dtmp6*(-Dtmp13 - Dtmp16 + Dtmp26*Dtmp30 + 1);
D[22] = Dtmp31*Dtmp33;
D[23] = -D[18] - D[21];
D[24] = Dtmp21*Dtmp36;
D[25] = Dtmp29*Dtmp37;
D[26] = -D[19] - D[24];
D[27] = Dtmp6*(Dtmp26*Dtmp39 + Dtmp38 + 3);
D[28] = Dtmp33*Dtmp36;
D[29] = -D[21] - D[27];
D[30] = -Dtmp20*(-70*Dtmp22 + 63*Dtmp40 + 15);
D[31] = -Dtmp44*y;
D[32] = -Dtmp44*z;
D[33] = -Dtmp20*(-21*Dtmp34 + Dtmp45 + Dtmp48);
D[34] = -Dtmp5*Dtmp52;
D[35] = -D[30] - D[33];
D[36] = -Dtmp32*(-21*Dtmp22 + Dtmp48 + Dtmp53);
D[37] = -Dtmp19*z*(Dtmp45 + Dtmp47 + Dtmp53 + 1);
D[38] = -D[31] - D[36];
D[39] = -Dtmp58*x;
D[40] = -Dtmp10*Dtmp52;
D[41] = -D[33] - D[39];
D[42] = -Dtmp32*(-70*Dtmp34 + 63*Dtmp55 + 15);
D[43] = -Dtmp58*z;
D[44] = -D[36] - D[42];
D[45] = Dtmp43*(105*Dtmp22 - 315*Dtmp40 + Dtmp59*pow(x, 6) - 5);
D[46] = Dtmp50*Dtmp61;
D[47] = Dtmp61*Dtmp62;
D[48] = Dtmp43*(Dtmp24*Dtmp59*Dtmp8 + Dtmp37 + Dtmp41 - Dtmp42 + Dtmp63);
D[49] = Dtmp51*(-18*Dtmp22 + Dtmp60 + 1);
D[50] = -D[45] - D[48];
D[51] = 945*Dtmp49*x*y*(11*Dtmp46 + Dtmp64 + Dtmp65);
D[52] = Dtmp67*(-9*Dtmp34 + Dtmp65 + Dtmp66);
D[53] = -D[46] - D[51];
D[54] = Dtmp43*(Dtmp1*Dtmp39*Dtmp59 + Dtmp31 + Dtmp54 - Dtmp56 + Dtmp63);
D[55] = Dtmp51*(-9*Dtmp22 + Dtmp64 + Dtmp66 + 1);
D[56] = -D[48] - D[54];
D[57] = Dtmp50*Dtmp68*x;
D[58] = Dtmp67*(4*Dtmp10*Dtmp34 + Dtmp57);
D[59] = -D[51] - D[57];
D[60] = Dtmp43*(105*Dtmp34 - 315*Dtmp55 + Dtmp59*pow(y, 6) - 5);
D[61] = Dtmp51*Dtmp68;
D[62] = -D[54] - D[60];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[9]*M[8] + D[10]*M[9] + D[11]*M[10] + D[12]*M[11] + D[13]*M[12] + D[15]*M[13] + D[16]*M[14] + D[18]*M[15] + D[19]*M[16] + D[20]*M[17] + D[21]*M[18] + D[22]*M[19] + D[24]*M[20] + D[25]*M[21] + D[27]*M[22] + D[28]*M[23] + D[30]*M[24] + D[31]*M[25] + D[32]*M[26] + D[33]*M[27] + D[34]*M[28] + D[36]*M[29] + D[37]*M[30] + D[39]*M[31] + D[40]*M[32] + D[42]*M[33] + D[43]*M[34] + D[45]*M[35] + D[46]*M[36] + D[47]*M[37] + D[48]*M[38] + D[49]*M[39] + D[51]*M[40] + D[52]*M[41] + D[54]*M[42] + D[55]*M[43] + D[57]*M[44] + D[58]*M[45] + D[60]*M[46] + D[61]*M[47];
#pragma omp atomic
L[1] += D[3]*M[0] + D[4]*M[1] + D[5]*M[2] + D[9]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[18]*M[8] + D[19]*M[9] + D[20]*M[10] + D[21]*M[11] + D[22]*M[12] + D[24]*M[13] + D[25]*M[14] + D[30]*M[15] + D[31]*M[16] + D[32]*M[17] + D[33]*M[18] + D[34]*M[19] + D[36]*M[20] + D[37]*M[21] + D[39]*M[22] + D[40]*M[23] + D[45]*M[24] + D[46]*M[25] + D[47]*M[26] + D[48]*M[27] + D[49]*M[28] + D[51]*M[29] + D[52]*M[30] + D[54]*M[31] + D[55]*M[32] + D[57]*M[33] + D[58]*M[34];
#pragma omp atomic
L[2] += D[4]*M[0] + D[6]*M[1] + D[7]*M[2] + D[10]*M[3] + D[12]*M[4] + D[13]*M[5] + D[15]*M[6] + D[16]*M[7] + D[19]*M[8] + D[21]*M[9] + D[22]*M[10] + D[24]*M[11] + D[25]*M[12] + D[27]*M[13] + D[28]*M[14] + D[31]*M[15] + D[33]*M[16] + D[34]*M[17] + D[36]*M[18] + D[37]*M[19] + D[39]*M[20] + D[40]*M[21] + D[42]*M[22] + D[43]*M[23] + D[46]*M[24] + D[48]*M[25] + D[49]*M[26] + D[51]*M[27] + D[52]*M[28] + D[54]*M[29] + D[55]*M[30] + D[57]*M[31] + D[58]*M[32] + D[60]*M[33] + D[61]*M[34];
#pragma omp atomic
L[3] += D[5]*M[0] + D[7]*M[1] + D[8]*M[2] + D[11]*M[3] + D[13]*M[4] + D[14]*M[5] + D[16]*M[6] + D[17]*M[7] + D[20]*M[8] + D[22]*M[9] + D[23]*M[10] + D[25]*M[11] + D[26]*M[12] + D[28]*M[13] + D[29]*M[14] + D[32]*M[15] + D[34]*M[16] + D[35]*M[17] + D[37]*M[18] + D[38]*M[19] + D[40]*M[20] + D[41]*M[21] + D[43]*M[22] + D[44]*M[23] + D[47]*M[24] + D[49]*M[25] + D[50]*M[26] + D[52]*M[27] + D[53]*M[28] + D[55]*M[29] + D[56]*M[30] + D[58]*M[31] + D[59]*M[32] + D[61]*M[33] + D[62]*M[34];
#pragma omp atomic
L[4] += D[9]*M[0] + D[10]*M[1] + D[11]*M[2] + D[18]*M[3] + D[19]*M[4] + D[20]*M[5] + D[21]*M[6] + D[22]*M[7] + D[30]*M[8] + D[31]*M[9] + D[32]*M[10] + D[33]*M[11] + D[34]*M[12] + D[36]*M[13] + D[37]*M[14] + D[45]*M[15] + D[46]*M[16] + D[47]*M[17] + D[48]*M[18] + D[49]*M[19] + D[51]*M[20] + D[52]*M[21] + D[54]*M[22] + D[55]*M[23];
#pragma omp atomic
L[5] += D[10]*M[0] + D[12]*M[1] + D[13]*M[2] + D[19]*M[3] + D[21]*M[4] + D[22]*M[5] + D[24]*M[6] + D[25]*M[7] + D[31]*M[8] + D[33]*M[9] + D[34]*M[10] + D[36]*M[11] + D[37]*M[12] + D[39]*M[13] + D[40]*M[14] + D[46]*M[15] + D[48]*M[16] + D[49]*M[17] + D[51]*M[18] + D[52]*M[19] + D[54]*M[20] + D[55]*M[21] + D[57]*M[22] + D[58]*M[23];
#pragma omp atomic
L[6] += D[11]*M[0] + D[13]*M[1] + D[14]*M[2] + D[20]*M[3] + D[22]*M[4] + D[23]*M[5] + D[25]*M[6] + D[26]*M[7] + D[32]*M[8] + D[34]*M[9] + D[35]*M[10] + D[37]*M[11] + D[38]*M[12] + D[40]*M[13] + D[41]*M[14] + D[47]*M[15] + D[49]*M[16] + D[50]*M[17] + D[52]*M[18] + D[53]*M[19] + D[55]*M[20] + D[56]*M[21] + D[58]*M[22] + D[59]*M[23];
#pragma omp atomic
L[7] += D[12]*M[0] + D[15]*M[1] + D[16]*M[2] + D[21]*M[3] + D[24]*M[4] + D[25]*M[5] + D[27]*M[6] + D[28]*M[7] + D[33]*M[8] + D[36]*M[9] + D[37]*M[10] + D[39]*M[11] + D[40]*M[12] + D[42]*M[13] + D[43]*M[14] + D[48]*M[15] + D[51]*M[16] + D[52]*M[17] + D[54]*M[18] + D[55]*M[19] + D[57]*M[20] + D[58]*M[21] + D[60]*M[22] + D[61]*M[23];
#pragma omp atomic
L[8] += D[13]*M[0] + D[16]*M[1] + D[17]*M[2] + D[22]*M[3] + D[25]*M[4] + D[26]*M[5] + D[28]*M[6] + D[29]*M[7] + D[34]*M[8] + D[37]*M[9] + D[38]*M[10] + D[40]*M[11] + D[41]*M[12] + D[43]*M[13] + D[44]*M[14] + D[49]*M[15] + D[52]*M[16] + D[53]*M[17] + D[55]*M[18] + D[56]*M[19] + D[58]*M[20] + D[59]*M[21] + D[61]*M[22] + D[62]*M[23];
#pragma omp atomic
L[9] += D[18]*M[0] + D[19]*M[1] + D[20]*M[2] + D[30]*M[3] + D[31]*M[4] + D[32]*M[5] + D[33]*M[6] + D[34]*M[7] + D[45]*M[8] + D[46]*M[9] + D[47]*M[10] + D[48]*M[11] + D[49]*M[12] + D[51]*M[13] + D[52]*M[14];
#pragma omp atomic
L[10] += D[19]*M[0] + D[21]*M[1] + D[22]*M[2] + D[31]*M[3] + D[33]*M[4] + D[34]*M[5] + D[36]*M[6] + D[37]*M[7] + D[46]*M[8] + D[48]*M[9] + D[49]*M[10] + D[51]*M[11] + D[52]*M[12] + D[54]*M[13] + D[55]*M[14];
#pragma omp atomic
L[11] += D[20]*M[0] + D[22]*M[1] + D[23]*M[2] + D[32]*M[3] + D[34]*M[4] + D[35]*M[5] + D[37]*M[6] + D[38]*M[7] + D[47]*M[8] + D[49]*M[9] + D[50]*M[10] + D[52]*M[11] + D[53]*M[12] + D[55]*M[13] + D[56]*M[14];
#pragma omp atomic
L[12] += D[21]*M[0] + D[24]*M[1] + D[25]*M[2] + D[33]*M[3] + D[36]*M[4] + D[37]*M[5] + D[39]*M[6] + D[40]*M[7] + D[48]*M[8] + D[51]*M[9] + D[52]*M[10] + D[54]*M[11] + D[55]*M[12] + D[57]*M[13] + D[58]*M[14];
#pragma omp atomic
L[13] += D[22]*M[0] + D[25]*M[1] + D[26]*M[2] + D[34]*M[3] + D[37]*M[4] + D[38]*M[5] + D[40]*M[6] + D[41]*M[7] + D[49]*M[8] + D[52]*M[9] + D[53]*M[10] + D[55]*M[11] + D[56]*M[12] + D[58]*M[13] + D[59]*M[14];
#pragma omp atomic
L[14] += D[24]*M[0] + D[27]*M[1] + D[28]*M[2] + D[36]*M[3] + D[39]*M[4] + D[40]*M[5] + D[42]*M[6] + D[43]*M[7] + D[51]*M[8] + D[54]*M[9] + D[55]*M[10] + D[57]*M[11] + D[58]*M[12] + D[60]*M[13] + D[61]*M[14];
#pragma omp atomic
L[15] += D[25]*M[0] + D[28]*M[1] + D[29]*M[2] + D[37]*M[3] + D[40]*M[4] + D[41]*M[5] + D[43]*M[6] + D[44]*M[7] + D[52]*M[8] + D[55]*M[9] + D[56]*M[10] + D[58]*M[11] + D[59]*M[12] + D[61]*M[13] + D[62]*M[14];
#pragma omp atomic
L[16] += D[30]*M[0] + D[31]*M[1] + D[32]*M[2] + D[45]*M[3] + D[46]*M[4] + D[47]*M[5] + D[48]*M[6] + D[49]*M[7];
#pragma omp atomic
L[17] += D[31]*M[0] + D[33]*M[1] + D[34]*M[2] + D[46]*M[3] + D[48]*M[4] + D[49]*M[5] + D[51]*M[6] + D[52]*M[7];
#pragma omp atomic
L[18] += D[32]*M[0] + D[34]*M[1] + D[35]*M[2] + D[47]*M[3] + D[49]*M[4] + D[50]*M[5] + D[52]*M[6] + D[53]*M[7];
#pragma omp atomic
L[19] += D[33]*M[0] + D[36]*M[1] + D[37]*M[2] + D[48]*M[3] + D[51]*M[4] + D[52]*M[5] + D[54]*M[6] + D[55]*M[7];
#pragma omp atomic
L[20] += D[34]*M[0] + D[37]*M[1] + D[38]*M[2] + D[49]*M[3] + D[52]*M[4] + D[53]*M[5] + D[55]*M[6] + D[56]*M[7];
#pragma omp atomic
L[21] += D[36]*M[0] + D[39]*M[1] + D[40]*M[2] + D[51]*M[3] + D[54]*M[4] + D[55]*M[5] + D[57]*M[6] + D[58]*M[7];
#pragma omp atomic
L[22] += D[37]*M[0] + D[40]*M[1] + D[41]*M[2] + D[52]*M[3] + D[55]*M[4] + D[56]*M[5] + D[58]*M[6] + D[59]*M[7];
#pragma omp atomic
L[23] += D[39]*M[0] + D[42]*M[1] + D[43]*M[2] + D[54]*M[3] + D[57]*M[4] + D[58]*M[5] + D[60]*M[6] + D[61]*M[7];
#pragma omp atomic
L[24] += D[40]*M[0] + D[43]*M[1] + D[44]*M[2] + D[55]*M[3] + D[58]*M[4] + D[59]*M[5] + D[61]*M[6] + D[62]*M[7];
#pragma omp atomic
L[25] += D[45]*M[0] + D[46]*M[1] + D[47]*M[2];
#pragma omp atomic
L[26] += D[46]*M[0] + D[48]*M[1] + D[49]*M[2];
#pragma omp atomic
L[27] += D[47]*M[0] + D[49]*M[1] + D[50]*M[2];
#pragma omp atomic
L[28] += D[48]*M[0] + D[51]*M[1] + D[52]*M[2];
#pragma omp atomic
L[29] += D[49]*M[0] + D[52]*M[1] + D[53]*M[2];
#pragma omp atomic
L[30] += D[51]*M[0] + D[54]*M[1] + D[55]*M[2];
#pragma omp atomic
L[31] += D[52]*M[0] + D[55]*M[1] + D[56]*M[2];
#pragma omp atomic
L[32] += D[54]*M[0] + D[57]*M[1] + D[58]*M[2];
#pragma omp atomic
L[33] += D[55]*M[0] + D[58]*M[1] + D[59]*M[2];
#pragma omp atomic
L[34] += D[57]*M[0] + D[60]*M[1] + D[61]*M[2];
#pragma omp atomic
L[35] += D[58]*M[0] + D[61]*M[1] + D[62]*M[2];

}

void S2M_7(double x, double y, double z, double * S, double * M) {
double Mtmp0 = x*S[1];
double Mtmp1 = y*S[0];
double Mtmp2 = Mtmp0 + Mtmp1;
double Mtmp3 = x*S[2];
double Mtmp4 = z*S[0];
double Mtmp5 = Mtmp3 + Mtmp4;
double Mtmp6 = y*S[2];
double Mtmp7 = z*S[1];
double Mtmp8 = Mtmp6 + Mtmp7;
double Mtmp9 = pow(x, 2);
double Mtmp10 = (1.0/2.0)*Mtmp0;
double Mtmp11 = (1.0/2.0)*Mtmp3;
double Mtmp12 = (1.0/2.0)*Mtmp1;
double Mtmp13 = Mtmp1*z;
double Mtmp14 = Mtmp3*y;
double Mtmp15 = Mtmp0*z;
double Mtmp16 = Mtmp14 + Mtmp15;
double Mtmp17 = Mtmp13 + Mtmp16;
double Mtmp18 = pow(y, 2);
double Mtmp19 = pow(z, 2);
double Mtmp20 = pow(x, 3);
double Mtmp21 = 3*Mtmp1;
double Mtmp22 = (1.0/6.0)*Mtmp9;
double Mtmp23 = 3*Mtmp4;
double Mtmp24 = (1.0/2.0)*x;
double Mtmp25 = Mtmp11*y;
double Mtmp26 = Mtmp10*z;
double Mtmp27 = 3*Mtmp0;
double Mtmp28 = (1.0/6.0)*Mtmp18;
double Mtmp29 = Mtmp12*z;
double Mtmp30 = 3*Mtmp3;
double Mtmp31 = (1.0/6.0)*Mtmp19;
double Mtmp32 = pow(y, 3);
double Mtmp33 = 3*Mtmp7;
double Mtmp34 = y*z;
double Mtmp35 = 3*Mtmp6;
double Mtmp36 = pow(z, 3);
double Mtmp37 = pow(x, 4);
double Mtmp38 = 4*Mtmp1;
double Mtmp39 = (1.0/24.0)*Mtmp20;
double Mtmp40 = 4*Mtmp4;
double Mtmp41 = 2*Mtmp0;
double Mtmp42 = (1.0/12.0)*Mtmp9;
double Mtmp43 = Mtmp42*y;
double Mtmp44 = 3*Mtmp13;
double Mtmp45 = 2*Mtmp3;
double Mtmp46 = Mtmp42*z;
double Mtmp47 = 2*Mtmp1;
double Mtmp48 = (1.0/12.0)*x;
double Mtmp49 = Mtmp18*Mtmp48;
double Mtmp50 = 2*Mtmp13;
double Mtmp51 = 2*Mtmp15;
double Mtmp52 = Mtmp14 + Mtmp51;
double Mtmp53 = (1.0/4.0)*x;
double Mtmp54 = 2*Mtmp14;
double Mtmp55 = Mtmp15 + Mtmp54;
double Mtmp56 = 2*Mtmp4;
double Mtmp57 = Mtmp19*Mtmp48;
double Mtmp58 = 4*Mtmp0;
double Mtmp59 = (1.0/24.0)*Mtmp32;
double Mtmp60 = 3*Mtmp15;
double Mtmp61 = Mtmp14 + Mtmp60;
double Mtmp62 = Mtmp13 + Mtmp54;
double Mtmp63 = 3*Mtmp14;
double Mtmp64 = Mtmp15 + Mtmp63;
double Mtmp65 = 4*Mtmp3;
double Mtmp66 = (1.0/24.0)*Mtmp36;
double Mtmp67 = pow(y, 4);
double Mtmp68 = 4*Mtmp7;
double Mtmp69 = 2*Mtmp6;
double Mtmp70 = (1.0/12.0)*Mtmp18;
double Mtmp71 = Mtmp70*z;
double Mtmp72 = 2*Mtmp7;
double Mtmp73 = Mtmp19*y;
double Mtmp74 = (1.0/12.0)*Mtmp73;
double Mtmp75 = 4*Mtmp6;
double Mtmp76 = pow(z, 4);
double Mtmp77 = pow(x, 5);
double Mtmp78 = 5*Mtmp1;
double Mtmp79 = (1.0/120.0)*Mtmp37;
double Mtmp80 = 5*Mtmp4;
double Mtmp81 = 4*Mtmp13;
double Mtmp82 = 5*Mtmp0;
double Mtmp83 = (1.0/120.0)*Mtmp67;
double Mtmp84 = 4*Mtmp15;
double Mtmp85 = Mtmp14 + Mtmp84;
double Mtmp86 = Mtmp51 + Mtmp63;
double Mtmp87 = 4*Mtmp14;
double Mtmp88 = Mtmp15 + Mtmp87;
double Mtmp89 = 5*Mtmp3;
double Mtmp90 = (1.0/120.0)*Mtmp76;
double Mtmp91 = pow(y, 5);
double Mtmp92 = 5*Mtmp7;
double Mtmp93 = 5*Mtmp6;
double Mtmp94 = pow(z, 5);
double Mtmp95 = (1.0/720.0)*Mtmp77;
double Mtmp96 = (1.0/240.0)*Mtmp37;
double Mtmp97 = (1.0/144.0)*Mtmp20;
double Mtmp98 = (1.0/48.0)*Mtmp20;
double Mtmp99 = (1.0/144.0)*Mtmp9;
double Mtmp100 = (1.0/36.0)*Mtmp9;
double Mtmp101 = (1.0/240.0)*x;
double Mtmp102 = (1.0/48.0)*x;
double Mtmp103 = (1.0/24.0)*x;
double Mtmp104 = (1.0/720.0)*Mtmp91;
double Mtmp105 = (1.0/720.0)*Mtmp94;
#pragma omp atomic
M[0] += S[0];
#pragma omp atomic
M[1] += S[1];
#pragma omp atomic
M[2] += S[2];
#pragma omp atomic
M[3] += x*S[0];
#pragma omp atomic
M[4] += Mtmp2;
#pragma omp atomic
M[5] += Mtmp5;
#pragma omp atomic
M[6] += y*S[1];
#pragma omp atomic
M[7] += Mtmp8;
#pragma omp atomic
M[8] += z*S[2];
#pragma omp atomic
M[9] += (1.0/2.0)*Mtmp9*S[0];
#pragma omp atomic
M[10] += x*(Mtmp1 + Mtmp10);
#pragma omp atomic
M[11] += x*(Mtmp11 + Mtmp4);
#pragma omp atomic
M[12] += y*(Mtmp0 + Mtmp12);
#pragma omp atomic
M[13] += Mtmp17;
#pragma omp atomic
M[14] += z*(Mtmp3 + (1.0/2.0)*Mtmp4);
#pragma omp atomic
M[15] += (1.0/2.0)*Mtmp18*S[1];
#pragma omp atomic
M[16] += y*((1.0/2.0)*Mtmp6 + Mtmp7);
#pragma omp atomic
M[17] += z*(Mtmp6 + (1.0/2.0)*Mtmp7);
#pragma omp atomic
M[18] += (1.0/2.0)*Mtmp19*S[2];
#pragma omp atomic
M[19] += (1.0/6.0)*Mtmp20*S[0];
#pragma omp atomic
M[20] += Mtmp22*(Mtmp0 + Mtmp21);
#pragma omp atomic
M[21] += Mtmp22*(Mtmp23 + Mtmp3);
#pragma omp atomic
M[22] += Mtmp2*Mtmp24*y;
#pragma omp atomic
M[23] += x*(Mtmp13 + Mtmp25 + Mtmp26);
#pragma omp atomic
M[24] += Mtmp24*Mtmp5*z;
#pragma omp atomic
M[25] += Mtmp28*(Mtmp1 + Mtmp27);
#pragma omp atomic
M[26] += y*(Mtmp15 + Mtmp25 + Mtmp29);
#pragma omp atomic
M[27] += z*(Mtmp14 + Mtmp26 + Mtmp29);
#pragma omp atomic
M[28] += Mtmp31*(Mtmp30 + Mtmp4);
#pragma omp atomic
M[29] += (1.0/6.0)*Mtmp32*S[1];
#pragma omp atomic
M[30] += Mtmp28*(Mtmp33 + Mtmp6);
#pragma omp atomic
M[31] += (1.0/2.0)*Mtmp34*Mtmp8;
#pragma omp atomic
M[32] += Mtmp31*(Mtmp35 + Mtmp7);
#pragma omp atomic
M[33] += (1.0/6.0)*Mtmp36*S[2];
#pragma omp atomic
M[34] += (1.0/24.0)*Mtmp37*S[0];
#pragma omp atomic
M[35] += Mtmp39*(Mtmp0 + Mtmp38);
#pragma omp atomic
M[36] += Mtmp39*(Mtmp3 + Mtmp40);
#pragma omp atomic
M[37] += Mtmp43*(Mtmp21 + Mtmp41);
#pragma omp atomic
M[38] += Mtmp22*(Mtmp16 + Mtmp44);
#pragma omp atomic
M[39] += Mtmp46*(Mtmp23 + Mtmp45);
#pragma omp atomic
M[40] += Mtmp49*(Mtmp27 + Mtmp47);
#pragma omp atomic
M[41] += Mtmp53*y*(Mtmp50 + Mtmp52);
#pragma omp atomic
M[42] += Mtmp53*z*(Mtmp50 + Mtmp55);
#pragma omp atomic
M[43] += Mtmp57*(Mtmp30 + Mtmp56);
#pragma omp atomic
M[44] += Mtmp59*(Mtmp1 + Mtmp58);
#pragma omp atomic
M[45] += Mtmp28*(Mtmp13 + Mtmp61);
#pragma omp atomic
M[46] += (1.0/4.0)*Mtmp34*(Mtmp51 + Mtmp62);
#pragma omp atomic
M[47] += Mtmp31*(Mtmp13 + Mtmp64);
#pragma omp atomic
M[48] += Mtmp66*(Mtmp4 + Mtmp65);
#pragma omp atomic
M[49] += (1.0/24.0)*Mtmp67*S[1];
#pragma omp atomic
M[50] += Mtmp59*(Mtmp6 + Mtmp68);
#pragma omp atomic
M[51] += Mtmp71*(Mtmp33 + Mtmp69);
#pragma omp atomic
M[52] += Mtmp74*(Mtmp35 + Mtmp72);
#pragma omp atomic
M[53] += Mtmp66*(Mtmp7 + Mtmp75);
#pragma omp atomic
M[54] += (1.0/24.0)*Mtmp76*S[2];
#pragma omp atomic
M[55] += (1.0/120.0)*Mtmp77*S[0];
#pragma omp atomic
M[56] += Mtmp79*(Mtmp0 + Mtmp78);
#pragma omp atomic
M[57] += Mtmp79*(Mtmp3 + Mtmp80);
#pragma omp atomic
M[58] += Mtmp39*y*(Mtmp0 + Mtmp47);
#pragma omp atomic
M[59] += Mtmp39*(Mtmp16 + Mtmp81);
#pragma omp atomic
M[60] += Mtmp39*z*(Mtmp3 + Mtmp56);
#pragma omp atomic
M[61] += Mtmp18*Mtmp2*Mtmp42;
#pragma omp atomic
M[62] += Mtmp43*(Mtmp44 + Mtmp52);
#pragma omp atomic
M[63] += Mtmp46*(Mtmp44 + Mtmp55);
#pragma omp atomic
M[64] += Mtmp19*Mtmp42*Mtmp5;
#pragma omp atomic
M[65] += Mtmp59*x*(Mtmp1 + Mtmp41);
#pragma omp atomic
M[66] += Mtmp49*(Mtmp50 + Mtmp61);
#pragma omp atomic
M[67] += Mtmp17*Mtmp34*Mtmp53;
#pragma omp atomic
M[68] += Mtmp57*(Mtmp50 + Mtmp64);
#pragma omp atomic
M[69] += Mtmp66*x*(Mtmp4 + Mtmp45);
#pragma omp atomic
M[70] += Mtmp83*(Mtmp1 + Mtmp82);
#pragma omp atomic
M[71] += Mtmp59*(Mtmp13 + Mtmp85);
#pragma omp atomic
M[72] += Mtmp71*(Mtmp60 + Mtmp62);
#pragma omp atomic
M[73] += Mtmp74*(Mtmp13 + Mtmp86);
#pragma omp atomic
M[74] += Mtmp66*(Mtmp13 + Mtmp88);
#pragma omp atomic
M[75] += Mtmp90*(Mtmp4 + Mtmp89);
#pragma omp atomic
M[76] += (1.0/120.0)*Mtmp91*S[1];
#pragma omp atomic
M[77] += Mtmp83*(Mtmp6 + Mtmp92);
#pragma omp atomic
M[78] += Mtmp59*z*(Mtmp6 + Mtmp72);
#pragma omp atomic
M[79] += Mtmp19*Mtmp70*Mtmp8;
#pragma omp atomic
M[80] += Mtmp66*y*(Mtmp69 + Mtmp7);
#pragma omp atomic
M[81] += Mtmp90*(Mtmp7 + Mtmp93);
#pragma omp atomic
M[82] += (1.0/120.0)*Mtmp94*S[2];
#pragma omp atomic
M[83] += (1.0/720.0)*pow(x, 6)*S[0];
#pragma omp atomic
M[84] += Mtmp95*(Mtmp0 + 6*Mtmp1);
#pragma omp atomic
M[85] += Mtmp95*(Mtmp3 + 6*Mtmp4);
#pragma omp atomic
M[86] += Mtmp96*y*(Mtmp41 + Mtmp78);
#pragma omp atomic
M[87] += Mtmp79*(5*Mtmp13 + Mtmp16);
#pragma omp atomic
M[88] += Mtmp96*z*(Mtmp45 + Mtmp80);
#pragma omp atomic
M[89] += Mtmp18*Mtmp97*(Mtmp27 + Mtmp38);
#pragma omp atomic
M[90] += Mtmp98*y*(Mtmp52 + Mtmp81);
#pragma omp atomic
M[91] += Mtmp98*z*(Mtmp55 + Mtmp81);
#pragma omp atomic
M[92] += Mtmp19*Mtmp97*(Mtmp30 + Mtmp40);
#pragma omp atomic
M[93] += Mtmp32*Mtmp99*(Mtmp21 + Mtmp58);
#pragma omp atomic
M[94] += Mtmp100*Mtmp18*(Mtmp44 + Mtmp61);
#pragma omp atomic
M[95] += (1.0/24.0)*Mtmp34*Mtmp9*(Mtmp44 + Mtmp51 + Mtmp54);
#pragma omp atomic
M[96] += Mtmp100*Mtmp19*(Mtmp44 + Mtmp64);
#pragma omp atomic
M[97] += Mtmp36*Mtmp99*(Mtmp23 + Mtmp65);
#pragma omp atomic
M[98] += Mtmp101*Mtmp67*(Mtmp47 + Mtmp82);
#pragma omp atomic
M[99] += Mtmp102*Mtmp32*(Mtmp50 + Mtmp85);
#pragma omp atomic
M[100] += Mtmp103*Mtmp18*z*(Mtmp50 + Mtmp54 + Mtmp60);
#pragma omp atomic
M[101] += Mtmp103*Mtmp73*(Mtmp50 + Mtmp86);
#pragma omp atomic
M[102] += Mtmp102*Mtmp36*(Mtmp50 + Mtmp88);
#pragma omp atomic
M[103] += Mtmp101*Mtmp76*(Mtmp56 + Mtmp89);
#pragma omp atomic
M[104] += Mtmp104*(6*Mtmp0 + Mtmp1);
#pragma omp atomic
M[105] += Mtmp83*(Mtmp13 + Mtmp14 + 5*Mtmp15);
#pragma omp atomic
M[106] += (1.0/48.0)*Mtmp32*z*(Mtmp62 + Mtmp84);
#pragma omp atomic
M[107] += (1.0/36.0)*Mtmp18*Mtmp19*(Mtmp13 + Mtmp60 + Mtmp63);
#pragma omp atomic
M[108] += (1.0/48.0)*Mtmp36*y*(Mtmp13 + Mtmp51 + Mtmp87);
#pragma omp atomic
M[109] += Mtmp90*(Mtmp13 + 5*Mtmp14 + Mtmp15);
#pragma omp atomic
M[110] += Mtmp105*(6*Mtmp3 + Mtmp4);
#pragma omp atomic
M[111] += (1.0/720.0)*pow(y, 6)*S[1];
#pragma omp atomic
M[112] += Mtmp104*(Mtmp6 + 6*Mtmp7);
#pragma omp atomic
M[113] += (1.0/240.0)*Mtmp67*z*(Mtmp69 + Mtmp92);
#pragma omp atomic
M[114] += (1.0/144.0)*Mtmp19*Mtmp32*(Mtmp35 + Mtmp68);
#pragma omp atomic
M[115] += (1.0/144.0)*Mtmp18*Mtmp36*(Mtmp33 + Mtmp75);
#pragma omp atomic
M[116] += (1.0/240.0)*Mtmp76*y*(Mtmp72 + Mtmp93);
#pragma omp atomic
M[117] += Mtmp105*(6*Mtmp6 + Mtmp7);
#pragma omp atomic
M[118] += (1.0/720.0)*pow(z, 6)*S[2];

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
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double D[119];
double Dtmp0 = pow(Rinv, 3);
double Dtmp1 = pow(x, 2);
double Dtmp2 = pow(Rinv, 2);
double Dtmp3 = 3*Dtmp2;
double Dtmp4 = Dtmp1*Dtmp3;
double Dtmp5 = Dtmp4 - 1;
double Dtmp6 = 3*pow(Rinv, 5);
double Dtmp7 = Dtmp6*x;
double Dtmp8 = pow(y, 2);
double Dtmp9 = Dtmp3*Dtmp8;
double Dtmp10 = Dtmp9 - 1;
double Dtmp11 = Dtmp6*y;
double Dtmp12 = 5*Dtmp2;
double Dtmp13 = Dtmp1*Dtmp12;
double Dtmp14 = Dtmp13 - 1;
double Dtmp15 = Dtmp6*z;
double Dtmp16 = Dtmp12*Dtmp8;
double Dtmp17 = Dtmp16 - 1;
double Dtmp18 = pow(Rinv, 7);
double Dtmp19 = 15*Dtmp18;
double Dtmp20 = Dtmp19*x;
double Dtmp21 = Dtmp20*y;
double Dtmp22 = Dtmp1*Dtmp2;
double Dtmp23 = 30*Dtmp22;
double Dtmp24 = -Dtmp23;
double Dtmp25 = pow(x, 4);
double Dtmp26 = pow(Rinv, 4);
double Dtmp27 = 35*Dtmp26;
double Dtmp28 = 7*Dtmp22;
double Dtmp29 = Dtmp28 - 3;
double Dtmp30 = Dtmp20*z;
double Dtmp31 = Dtmp1*Dtmp8;
double Dtmp32 = Dtmp28 - 1;
double Dtmp33 = Dtmp19*y;
double Dtmp34 = Dtmp33*z;
double Dtmp35 = Dtmp2*Dtmp8;
double Dtmp36 = 7*Dtmp35;
double Dtmp37 = Dtmp36 - 3;
double Dtmp38 = Dtmp36 - 1;
double Dtmp39 = 30*Dtmp35;
double Dtmp40 = -Dtmp39;
double Dtmp41 = pow(y, 4);
double Dtmp42 = Dtmp25*Dtmp26;
double Dtmp43 = 14*Dtmp22;
double Dtmp44 = 21*Dtmp42;
double Dtmp45 = 45*Dtmp18;
double Dtmp46 = Dtmp45*(-Dtmp43 + Dtmp44 + 1);
double Dtmp47 = -Dtmp28;
double Dtmp48 = Dtmp26*Dtmp31;
double Dtmp49 = 63*Dtmp48;
double Dtmp50 = Dtmp49 + 3;
double Dtmp51 = pow(Rinv, 9);
double Dtmp52 = 315*Dtmp51;
double Dtmp53 = Dtmp52*x;
double Dtmp54 = Dtmp53*z;
double Dtmp55 = Dtmp54*y;
double Dtmp56 = -Dtmp36;
double Dtmp57 = 14*Dtmp35;
double Dtmp58 = Dtmp26*Dtmp41;
double Dtmp59 = 21*Dtmp58;
double Dtmp60 = -Dtmp57 + Dtmp59 + 1;
double Dtmp61 = Dtmp45*Dtmp60;
double Dtmp62 = pow(x, 6);
double Dtmp63 = pow(Rinv, 6);
double Dtmp64 = 231*Dtmp63;
double Dtmp65 = 33*Dtmp42;
double Dtmp66 = Dtmp53*(Dtmp24 + Dtmp65 + 5);
double Dtmp67 = -126*Dtmp48;
double Dtmp68 = Dtmp25*Dtmp8;
double Dtmp69 = 18*Dtmp22;
double Dtmp70 = Dtmp52*z;
double Dtmp71 = Dtmp70*y;
double Dtmp72 = -Dtmp9;
double Dtmp73 = 1 - Dtmp4;
double Dtmp74 = 945*Dtmp51;
double Dtmp75 = Dtmp74*y;
double Dtmp76 = 9*Dtmp35;
double Dtmp77 = 33*Dtmp48;
double Dtmp78 = Dtmp1*Dtmp41;
double Dtmp79 = 9*Dtmp22;
double Dtmp80 = 33*Dtmp58;
double Dtmp81 = Dtmp40 + Dtmp80 + 5;
double Dtmp82 = 4*Dtmp35;
double Dtmp83 = pow(y, 6);
double Dtmp84 = 429*Dtmp63;
double Dtmp85 = Dtmp62*Dtmp84;
double Dtmp86 = Dtmp52*(135*Dtmp22 - 495*Dtmp42 + Dtmp85 - 5);
double Dtmp87 = -Dtmp65;
double Dtmp88 = Dtmp68*Dtmp84 + Dtmp87;
double Dtmp89 = -330*Dtmp48 - 5;
double Dtmp90 = 945*pow(Rinv, 11)*x*y*z;
double Dtmp91 = -66*Dtmp48;
double Dtmp92 = 143*Dtmp63;
double Dtmp93 = -198*Dtmp48 - 1;
double Dtmp94 = -Dtmp80;
double Dtmp95 = 18*Dtmp35 + Dtmp94;
double Dtmp96 = Dtmp78*Dtmp84;
double Dtmp97 = Dtmp52*y;
double Dtmp98 = Dtmp83*Dtmp84;
double Dtmp99 = 135*Dtmp35 - 495*Dtmp58 + Dtmp98 - 5;
D[0] = -Dtmp0*x;
D[1] = -Dtmp0*y;
D[2] = -Dtmp0*z;
D[3] = Dtmp0*Dtmp5;
D[4] = Dtmp7*y;
D[5] = Dtmp7*z;
D[6] = Dtmp0*Dtmp10;
D[7] = Dtmp11*z;
D[8] = -D[3] - D[6];
D[9] = -Dtmp7*(Dtmp13 - 3);
D[10] = -Dtmp11*Dtmp14;
D[11] = -Dtmp14*Dtmp15;
D[12] = -Dtmp17*Dtmp7;
D[13] = -Dtmp21*z;
D[14] = -D[9] - D[12];
D[15] = -Dtmp11*(Dtmp16 - 3);
D[16] = -Dtmp15*Dtmp17;
D[17] = -D[10] - D[15];
D[18] = -D[11] - D[16];
D[19] = Dtmp6*(Dtmp24 + Dtmp25*Dtmp27 + 3);
D[20] = Dtmp21*Dtmp29;
D[21] = Dtmp29*Dtmp30;
D[22] = Dtmp6*(-Dtmp13 - Dtmp16 + Dtmp27*Dtmp31 + 1);
D[23] = Dtmp32*Dtmp34;
D[24] = -D[19] - D[22];
D[25] = Dtmp21*Dtmp37;
D[26] = Dtmp30*Dtmp38;
D[27] = -D[20] - D[25];
D[28] = -D[21] - D[26];
D[29] = Dtmp6*(Dtmp27*Dtmp41 + Dtmp40 + 3);
D[30] = Dtmp34*Dtmp37;
D[31] = -D[22] - D[29];
D[32] = -D[23] - D[30];
D[33] = -D[24] - D[31];
D[34] = -Dtmp20*(-70*Dtmp22 + 63*Dtmp42 + 15);
D[35] = -Dtmp46*y;
D[36] = -Dtmp46*z;
D[37] = -Dtmp20*(-21*Dtmp35 + Dtmp47 + Dtmp50);
D[38] = -Dtmp5*Dtmp55;
D[39] = -D[34] - D[37];
D[40] = -Dtmp33*(-21*Dtmp22 + Dtmp50 + Dtmp56);
D[41] = -Dtmp19*z*(Dtmp47 + Dtmp49 + Dtmp56 + 1);
D[42] = -D[35] - D[40];
D[43] = -D[36] - D[41];
D[44] = -Dtmp61*x;
D[45] = -Dtmp10*Dtmp55;
D[46] = -D[37] - D[44];
D[47] = -D[38] - D[45];
D[48] = -D[39] - D[46];
D[49] = -Dtmp33*(-70*Dtmp35 + 63*Dtmp58 + 15);
D[50] = -Dtmp61*z;
D[51] = -D[40] - D[49];
D[52] = -D[41] - D[50];
D[53] = -D[42] - D[51];
D[54] = -D[43] - D[52];
D[55] = Dtmp45*(105*Dtmp22 - 315*Dtmp42 + Dtmp62*Dtmp64 - 5);
D[56] = Dtmp66*y;
D[57] = Dtmp66*z;
D[58] = Dtmp45*(Dtmp38 + Dtmp43 - Dtmp44 + Dtmp64*Dtmp68 + Dtmp67);
D[59] = Dtmp71*(Dtmp65 - Dtmp69 + 1);
D[60] = -D[55] - D[58];
D[61] = Dtmp75*x*(11*Dtmp48 + Dtmp72 + Dtmp73);
D[62] = Dtmp54*(Dtmp73 - Dtmp76 + Dtmp77);
D[63] = -D[56] - D[61];
D[64] = -D[57] - D[62];
D[65] = Dtmp45*(Dtmp32 + Dtmp57 - Dtmp59 + Dtmp64*Dtmp78 + Dtmp67);
D[66] = Dtmp71*(Dtmp72 + Dtmp77 - Dtmp79 + 1);
D[67] = -D[58] - D[65];
D[68] = -D[59] - D[66];
D[69] = -D[60] - D[67];
D[70] = Dtmp53*Dtmp81*y;
D[71] = Dtmp54*(Dtmp10*Dtmp82 + Dtmp60);
D[72] = -D[61] - D[70];
D[73] = -D[62] - D[71];
D[74] = -D[63] - D[72];
D[75] = -D[64] - D[73];
D[76] = Dtmp45*(105*Dtmp35 - 315*Dtmp58 + Dtmp64*Dtmp83 - 5);
D[77] = Dtmp71*Dtmp81;
D[78] = -D[65] - D[76];
D[79] = -D[66] - D[77];
D[80] = -D[67] - D[78];
D[81] = -D[68] - D[79];
D[82] = -D[69] - D[80];
D[83] = -Dtmp53*(315*Dtmp22 - 693*Dtmp42 + Dtmp85 - 35);
D[84] = -Dtmp86*y;
D[85] = -Dtmp86*z;
D[86] = -Dtmp53*(Dtmp23 + 45*Dtmp35 + Dtmp88 + Dtmp89);
D[87] = -Dtmp90*(-110*Dtmp22 + 143*Dtmp42 + 15);
D[88] = -D[83] - D[86];
D[89] = -Dtmp75*(Dtmp10 + Dtmp68*Dtmp92 + Dtmp69 + Dtmp87 + Dtmp91);
D[90] = -Dtmp70*(Dtmp69 + Dtmp76 + Dtmp88 + Dtmp93);
D[91] = -D[84] - D[89];
D[92] = -D[85] - D[90];
D[93] = -Dtmp74*x*(Dtmp5 + Dtmp78*Dtmp92 + Dtmp91 + Dtmp95);
D[94] = -Dtmp90*(-33*Dtmp22 - 33*Dtmp35 + 143*Dtmp48 + 9);
D[95] = -D[86] - D[93];
D[96] = -D[87] - D[94];
D[97] = -D[88] - D[95];
D[98] = -Dtmp97*(45*Dtmp22 + Dtmp39 + Dtmp89 + Dtmp94 + Dtmp96);
D[99] = -Dtmp70*(Dtmp79 + Dtmp93 + Dtmp95 + Dtmp96);
D[100] = -D[89] - D[98];
D[101] = -D[90] - D[99];
D[102] = -D[91] - D[100];
D[103] = -D[92] - D[101];
D[104] = -Dtmp53*Dtmp99;
D[105] = -Dtmp90*(-90*Dtmp35 + 99*Dtmp58 + Dtmp82*(11*Dtmp35 - 5) + 15);
D[106] = -D[93] - D[104];
D[107] = -D[94] - D[105];
D[108] = -D[95] - D[106];
D[109] = -D[96] - D[107];
D[110] = -D[97] - D[108];
D[111] = -Dtmp97*(315*Dtmp35 - 693*Dtmp58 + Dtmp98 - 35);
D[112] = -Dtmp70*Dtmp99;
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
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(Rinv, 3);
double Ftmp1 = pow(Rinv, 2);
double Ftmp2 = 3*Ftmp1;
double Ftmp3 = Ftmp2*z;
double Ftmp4 = Ftmp2*x;
double Ftmp5 = Ftmp4*y;
double Ftmp6 = Ftmp3*M[2];
double Ftmp7 = pow(Rinv, 4);
double Ftmp8 = pow(x, 2);
double Ftmp9 = Ftmp1*Ftmp8;
double Ftmp10 = 3*Ftmp9;
double Ftmp11 = pow(y, 2);
double Ftmp12 = pow(z, 2);
double Ftmp13 = pow(y, 3);
double Ftmp14 = Ftmp13*x;
double Ftmp15 = pow(Rinv, 6);
double Ftmp16 = 30*Ftmp15;
double Ftmp17 = Ftmp16*M[15];
double Ftmp18 = Ftmp16*x;
double Ftmp19 = pow(z, 3);
double Ftmp20 = Ftmp19*M[18];
double Ftmp21 = Ftmp12*M[17];
double Ftmp22 = Ftmp18*y;
double Ftmp23 = z*M[16];
double Ftmp24 = Ftmp11*Ftmp18;
double Ftmp25 = Ftmp15*Ftmp8;
double Ftmp26 = 105*Ftmp25;
double Ftmp27 = pow(Rinv, 8);
double Ftmp28 = Ftmp11*Ftmp8;
double Ftmp29 = Ftmp16*Ftmp28;
double Ftmp30 = Ftmp12*M[14];
double Ftmp31 = 30*Ftmp25;
double Ftmp32 = pow(Rinv, 10);
double Ftmp33 = Ftmp32*Ftmp8;
double Ftmp34 = 1890*Ftmp33;
double Ftmp35 = Ftmp19*y;
double Ftmp36 = Ftmp13*z;
double Ftmp37 = 5*Ftmp9;
double Ftmp38 = Ftmp37 - 3;
double Ftmp39 = Ftmp1*Ftmp11;
double Ftmp40 = 5*Ftmp39;
double Ftmp41 = Ftmp40 - 1;
double Ftmp42 = Ftmp1*Ftmp12;
double Ftmp43 = 5*Ftmp42;
double Ftmp44 = Ftmp43 - 1;
double Ftmp45 = Ftmp9 - 1;
double Ftmp46 = Ftmp10 - 1;
double Ftmp47 = 3*Ftmp11;
double Ftmp48 = Ftmp1*Ftmp47;
double Ftmp49 = Ftmp48 - 1;
double Ftmp50 = 3*Ftmp12;
double Ftmp51 = Ftmp1*Ftmp50;
double Ftmp52 = Ftmp51 - 1;
double Ftmp53 = 15*Ftmp7;
double Ftmp54 = Ftmp53*y;
double Ftmp55 = 7*Ftmp9;
double Ftmp56 = Ftmp55 - 3;
double Ftmp57 = Ftmp56*M[20];
double Ftmp58 = 7*Ftmp39;
double Ftmp59 = Ftmp58 - 3;
double Ftmp60 = Ftmp59*M[25];
double Ftmp61 = 7*Ftmp42;
double Ftmp62 = Ftmp61 - 1;
double Ftmp63 = Ftmp62*M[27];
double Ftmp64 = Ftmp53*z;
double Ftmp65 = Ftmp56*M[21];
double Ftmp66 = Ftmp58 - 1;
double Ftmp67 = Ftmp66*M[26];
double Ftmp68 = Ftmp61 - 3;
double Ftmp69 = Ftmp64*Ftmp68;
double Ftmp70 = 30*Ftmp45;
double Ftmp71 = Ftmp7*x;
double Ftmp72 = Ftmp70*Ftmp71;
double Ftmp73 = y*M[10];
double Ftmp74 = Ftmp7*Ftmp8;
double Ftmp75 = Ftmp37 - 1;
double Ftmp76 = Ftmp75*M[10];
double Ftmp77 = Ftmp54*x;
double Ftmp78 = Ftmp40 - 3;
double Ftmp79 = Ftmp78*M[15];
double Ftmp80 = Ftmp44*M[17];
double Ftmp81 = Ftmp64*x;
double Ftmp82 = Ftmp75*M[11];
double Ftmp83 = Ftmp41*M[16];
double Ftmp84 = Ftmp43 - 3;
double Ftmp85 = Ftmp84*M[18];
double Ftmp86 = Ftmp38*M[9];
double Ftmp87 = 15*Ftmp8;
double Ftmp88 = Ftmp7*Ftmp87;
double Ftmp89 = Ftmp41*M[12];
double Ftmp90 = Ftmp44*M[14];
double Ftmp91 = Ftmp55 - 1;
double Ftmp92 = 9*Ftmp39;
double Ftmp93 = 420*Ftmp27;
double Ftmp94 = Ftmp93*(Ftmp92 - 5)*M[49];
double Ftmp95 = Ftmp19*x;
double Ftmp96 = 9*Ftmp42;
double Ftmp97 = Ftmp93*(Ftmp96 - 5)*M[54];
double Ftmp98 = Ftmp27*Ftmp8;
double Ftmp99 = 1890*Ftmp98;
double Ftmp100 = y*z;
double Ftmp101 = Ftmp12*Ftmp27;
double Ftmp102 = x*y;
double Ftmp103 = 1260*Ftmp52;
double Ftmp104 = Ftmp101*Ftmp102*Ftmp103;
double Ftmp105 = x*z;
double Ftmp106 = Ftmp11*Ftmp27;
double Ftmp107 = 1260*Ftmp49;
double Ftmp108 = Ftmp106*Ftmp107*M[50];
double Ftmp109 = Ftmp46*M[38];
double Ftmp110 = 2835*Ftmp98;
double Ftmp111 = Ftmp100*Ftmp110;
double Ftmp112 = Ftmp49*M[45];
double Ftmp113 = Ftmp52*M[47];
double Ftmp114 = 11*Ftmp42;
double Ftmp115 = Ftmp114 - 5;
double Ftmp116 = 11*Ftmp39;
double Ftmp117 = Ftmp116 - 5;
double Ftmp118 = Ftmp27*Ftmp28;
double Ftmp119 = Ftmp107*M[44];
double Ftmp120 = Ftmp103*Ftmp12;
double Ftmp121 = Ftmp114 - 3;
double Ftmp122 = Ftmp116 - 3;
double Ftmp123 = pow(Rinv, 12);
double Ftmp124 = Ftmp123*Ftmp8;
double Ftmp125 = 41580*M[109];
double Ftmp126 = Ftmp125*(13*Ftmp42 - 5);
double Ftmp127 = 41580*(13*Ftmp39 - 5)*M[105];
double Ftmp128 = pow(x, 4);
double Ftmp129 = Ftmp128*Ftmp7;
double Ftmp130 = 63*Ftmp129;
double Ftmp131 = Ftmp130 - 70*Ftmp9 + 15;
double Ftmp132 = pow(y, 4);
double Ftmp133 = 21*Ftmp7;
double Ftmp134 = Ftmp132*Ftmp133;
double Ftmp135 = 14*Ftmp39;
double Ftmp136 = -Ftmp135;
double Ftmp137 = Ftmp136 + 1;
double Ftmp138 = Ftmp134 + Ftmp137;
double Ftmp139 = pow(z, 4);
double Ftmp140 = Ftmp133*Ftmp139;
double Ftmp141 = 14*Ftmp42;
double Ftmp142 = -Ftmp141;
double Ftmp143 = Ftmp142 + 1;
double Ftmp144 = Ftmp140 + Ftmp143;
double Ftmp145 = 10*Ftmp9;
double Ftmp146 = -Ftmp145;
double Ftmp147 = Ftmp146 + 3;
double Ftmp148 = 30*Ftmp9;
double Ftmp149 = -Ftmp148;
double Ftmp150 = 35*Ftmp129 + Ftmp149 + 3;
double Ftmp151 = 30*Ftmp39;
double Ftmp152 = -Ftmp151;
double Ftmp153 = Ftmp132*Ftmp7;
double Ftmp154 = Ftmp152 + 35*Ftmp153 + 3;
double Ftmp155 = 30*Ftmp42;
double Ftmp156 = -Ftmp155;
double Ftmp157 = Ftmp139*Ftmp7;
double Ftmp158 = Ftmp156 + 35*Ftmp157 + 3;
double Ftmp159 = 315*Ftmp15;
double Ftmp160 = Ftmp159*y;
double Ftmp161 = 33*Ftmp129;
double Ftmp162 = Ftmp149 + Ftmp161 + 5;
double Ftmp163 = Ftmp162*M[56];
double Ftmp164 = 33*Ftmp153;
double Ftmp165 = Ftmp152 + Ftmp164 + 5;
double Ftmp166 = Ftmp165*M[70];
double Ftmp167 = 18*Ftmp42;
double Ftmp168 = -Ftmp167;
double Ftmp169 = 33*Ftmp157;
double Ftmp170 = Ftmp168 + Ftmp169 + 1;
double Ftmp171 = Ftmp170*M[74];
double Ftmp172 = Ftmp159*z;
double Ftmp173 = Ftmp162*M[57];
double Ftmp174 = Ftmp156 + Ftmp169 + 5;
double Ftmp175 = Ftmp172*Ftmp174;
double Ftmp176 = 1260*Ftmp15;
double Ftmp177 = Ftmp176*(3*Ftmp129 - 4*Ftmp9 + 1);
double Ftmp178 = 21*Ftmp129;
double Ftmp179 = 14*Ftmp9;
double Ftmp180 = -Ftmp179;
double Ftmp181 = Ftmp180 + 1;
double Ftmp182 = Ftmp178 + Ftmp181;
double Ftmp183 = Ftmp182*M[35];
double Ftmp184 = Ftmp160*x;
double Ftmp185 = 63*Ftmp153;
double Ftmp186 = Ftmp185 - 70*Ftmp39 + 15;
double Ftmp187 = Ftmp186*M[49];
double Ftmp188 = 105*Ftmp15;
double Ftmp189 = Ftmp188*x;
double Ftmp190 = Ftmp189*y;
double Ftmp191 = Ftmp144*M[53];
double Ftmp192 = Ftmp182*M[36];
double Ftmp193 = Ftmp172*x;
double Ftmp194 = Ftmp138*M[50];
double Ftmp195 = 63*Ftmp157;
double Ftmp196 = Ftmp195 - 70*Ftmp42 + 15;
double Ftmp197 = Ftmp196*M[54];
double Ftmp198 = Ftmp189*z;
double Ftmp199 = 143*Ftmp129;
double Ftmp200 = Ftmp199 - 110*Ftmp9 + 15;
double Ftmp201 = 143*Ftmp157;
double Ftmp202 = Ftmp201 - 110*Ftmp42 + 15;
double Ftmp203 = 420*M[34];
double Ftmp204 = Ftmp131*M[34];
double Ftmp205 = Ftmp138*M[44];
double Ftmp206 = Ftmp159*Ftmp8;
double Ftmp207 = Ftmp144*M[48];
double Ftmp208 = 11*Ftmp129;
double Ftmp209 = Ftmp180 + 3;
double Ftmp210 = 18*Ftmp9;
double Ftmp211 = -Ftmp210;
double Ftmp212 = Ftmp161 + Ftmp211 + 1;
double Ftmp213 = 18*Ftmp39;
double Ftmp214 = -Ftmp12*Ftmp213;
double Ftmp215 = -16*Ftmp9;
double Ftmp216 = Ftmp208 + Ftmp215 + 5;
double Ftmp217 = 143*Ftmp153;
double Ftmp218 = 5670*M[111];
double Ftmp219 = Ftmp218*(Ftmp217 - 154*Ftmp39 + 35);
double Ftmp220 = 5670*M[118];
double Ftmp221 = Ftmp220*Ftmp32*(Ftmp201 - 154*Ftmp42 + 35);
double Ftmp222 = 5670*Ftmp102;
double Ftmp223 = Ftmp222*Ftmp32;
double Ftmp224 = Ftmp12*Ftmp202;
double Ftmp225 = Ftmp224*M[117];
double Ftmp226 = Ftmp217 - 110*Ftmp39 + 15;
double Ftmp227 = Ftmp226*M[112];
double Ftmp228 = 5670*Ftmp105;
double Ftmp229 = Ftmp228*Ftmp32;
double Ftmp230 = Ftmp11*Ftmp229;
double Ftmp231 = Ftmp100*Ftmp33;
double Ftmp232 = 41580*M[87];
double Ftmp233 = Ftmp200*M[87];
double Ftmp234 = 10395*Ftmp231;
double Ftmp235 = Ftmp202*M[109];
double Ftmp236 = Ftmp226*M[104];
double Ftmp237 = 5670*Ftmp32;
double Ftmp238 = Ftmp237*Ftmp28;
double Ftmp239 = 5670*Ftmp33;
double Ftmp240 = Ftmp224*M[110];
double Ftmp241 = Ftmp12*Ftmp39;
double Ftmp242 = 1890*Ftmp32;
double Ftmp243 = Ftmp102*Ftmp242;
double Ftmp244 = Ftmp243*z;
double Ftmp245 = pow(x, 6);
double Ftmp246 = 429*Ftmp15;
double Ftmp247 = Ftmp245*Ftmp246;
double Ftmp248 = -693*Ftmp129 + Ftmp247 + 315*Ftmp9 - 35;
double Ftmp249 = pow(y, 6);
double Ftmp250 = Ftmp246*Ftmp249;
double Ftmp251 = -495*Ftmp153 + Ftmp250 + 135*Ftmp39 - 5;
double Ftmp252 = pow(z, 6);
double Ftmp253 = Ftmp246*Ftmp252;
double Ftmp254 = -495*Ftmp157 + Ftmp253 + 135*Ftmp42 - 5;
double Ftmp255 = Ftmp15*Ftmp245;
double Ftmp256 = -315*Ftmp129 + 231*Ftmp255 + 105*Ftmp9 - 5;
double Ftmp257 = 231*Ftmp15;
double Ftmp258 = -315*Ftmp153 + Ftmp249*Ftmp257 + 105*Ftmp39 - 5;
double Ftmp259 = -315*Ftmp157 + Ftmp252*Ftmp257 + 105*Ftmp42 - 5;
double Ftmp260 = 125*Ftmp9;
double Ftmp261 = 143*Ftmp255;
double Ftmp262 = -253*Ftmp129 + Ftmp260 + Ftmp261 - 15;
double Ftmp263 = Ftmp222*Ftmp27;
double Ftmp264 = -495*Ftmp129 + Ftmp247 + 135*Ftmp9 - 5;
double Ftmp265 = Ftmp264*M[84];
double Ftmp266 = 2835*Ftmp27;
double Ftmp267 = Ftmp102*Ftmp266;
double Ftmp268 = -693*Ftmp153 + Ftmp250 + 315*Ftmp39 - 35;
double Ftmp269 = Ftmp268*M[111];
double Ftmp270 = Ftmp254*M[117];
double Ftmp271 = Ftmp228*Ftmp27;
double Ftmp272 = Ftmp264*M[85];
double Ftmp273 = Ftmp105*Ftmp266;
double Ftmp274 = Ftmp251*M[112];
double Ftmp275 = -693*Ftmp157 + Ftmp253 + 315*Ftmp42 - 35;
double Ftmp276 = Ftmp275*M[118];
double Ftmp277 = 21*Ftmp39;
double Ftmp278 = -Ftmp277;
double Ftmp279 = Ftmp11*Ftmp7;
double Ftmp280 = 63*Ftmp8;
double Ftmp281 = Ftmp279*Ftmp280;
double Ftmp282 = -Ftmp55;
double Ftmp283 = Ftmp282 + 3;
double Ftmp284 = Ftmp278 + Ftmp281 + Ftmp283;
double Ftmp285 = 21*Ftmp42;
double Ftmp286 = -Ftmp285;
double Ftmp287 = Ftmp12*Ftmp7;
double Ftmp288 = Ftmp280*Ftmp287;
double Ftmp289 = Ftmp283 + Ftmp286 + Ftmp288;
double Ftmp290 = 5670*M[83];
double Ftmp291 = Ftmp248*M[83];
double Ftmp292 = Ftmp251*M[104];
double Ftmp293 = Ftmp254*M[110];
double Ftmp294 = 8*Ftmp39;
double Ftmp295 = -Ftmp294;
double Ftmp296 = Ftmp279*Ftmp8;
double Ftmp297 = 14*Ftmp296;
double Ftmp298 = -Ftmp9;
double Ftmp299 = Ftmp298 + 1;
double Ftmp300 = -Ftmp40;
double Ftmp301 = 35*Ftmp11;
double Ftmp302 = 1 - Ftmp37;
double Ftmp303 = Ftmp300 + Ftmp301*Ftmp74 + Ftmp302;
double Ftmp304 = 8*Ftmp42;
double Ftmp305 = -Ftmp304;
double Ftmp306 = Ftmp287*Ftmp8;
double Ftmp307 = 14*Ftmp306;
double Ftmp308 = -Ftmp43;
double Ftmp309 = Ftmp302 + 35*Ftmp306 + Ftmp308;
double Ftmp310 = Ftmp287*Ftmp301 + Ftmp300 + Ftmp308 + 1;
double Ftmp311 = 945*Ftmp15;
double Ftmp312 = Ftmp311*y;
double Ftmp313 = -Ftmp48;
double Ftmp314 = 11*Ftmp8;
double Ftmp315 = -Ftmp10;
double Ftmp316 = Ftmp315 + 1;
double Ftmp317 = Ftmp279*Ftmp314 + Ftmp313 + Ftmp316;
double Ftmp318 = Ftmp317*M[61];
double Ftmp319 = 33*Ftmp8;
double Ftmp320 = Ftmp287*Ftmp319;
double Ftmp321 = Ftmp316 + Ftmp320 - Ftmp96;
double Ftmp322 = Ftmp321*M[63];
double Ftmp323 = Ftmp279*Ftmp319;
double Ftmp324 = Ftmp316 + Ftmp323 - Ftmp92;
double Ftmp325 = Ftmp324*M[62];
double Ftmp326 = Ftmp311*z;
double Ftmp327 = -Ftmp51;
double Ftmp328 = Ftmp287*Ftmp314 + Ftmp316 + Ftmp327;
double Ftmp329 = Ftmp328*M[64];
double Ftmp330 = 18*Ftmp8;
double Ftmp331 = Ftmp279*Ftmp330;
double Ftmp332 = 10*Ftmp11;
double Ftmp333 = Ftmp1*Ftmp332;
double Ftmp334 = -Ftmp333;
double Ftmp335 = Ftmp334 + 3;
double Ftmp336 = 210*Ftmp15;
double Ftmp337 = Ftmp102*Ftmp336;
double Ftmp338 = -21*Ftmp9;
double Ftmp339 = -Ftmp58;
double Ftmp340 = Ftmp339 + 3;
double Ftmp341 = Ftmp281 + Ftmp338 + Ftmp340;
double Ftmp342 = Ftmp341*M[40];
double Ftmp343 = 10*Ftmp12;
double Ftmp344 = Ftmp1*Ftmp343;
double Ftmp345 = -Ftmp344;
double Ftmp346 = Ftmp287*Ftmp330;
double Ftmp347 = Ftmp282 + 1;
double Ftmp348 = -Ftmp61;
double Ftmp349 = Ftmp288 + Ftmp348;
double Ftmp350 = Ftmp347 + Ftmp349;
double Ftmp351 = Ftmp350*M[42];
double Ftmp352 = Ftmp12*Ftmp279;
double Ftmp353 = 63*Ftmp352;
double Ftmp354 = Ftmp286 + Ftmp340 + Ftmp353;
double Ftmp355 = Ftmp354*M[51];
double Ftmp356 = Ftmp299 + Ftmp331;
double Ftmp357 = Ftmp105*Ftmp336;
double Ftmp358 = Ftmp281 + Ftmp339 + Ftmp347;
double Ftmp359 = Ftmp358*M[41];
double Ftmp360 = Ftmp345 + 3;
double Ftmp361 = Ftmp338 + Ftmp349 + 3;
double Ftmp362 = Ftmp361*M[43];
double Ftmp363 = Ftmp278 + Ftmp348 + Ftmp353 + 3;
double Ftmp364 = Ftmp363*M[52];
double Ftmp365 = 9 - 33*Ftmp9;
double Ftmp366 = 143*Ftmp296 + Ftmp365 - 33*Ftmp39;
double Ftmp367 = 33*Ftmp42;
double Ftmp368 = 143*Ftmp306 + Ftmp365 - Ftmp367;
double Ftmp369 = 12*Ftmp39;
double Ftmp370 = -Ftmp369;
double Ftmp371 = 210*Ftmp25;
double Ftmp372 = Ftmp284*M[37];
double Ftmp373 = 12*Ftmp42;
double Ftmp374 = -Ftmp373;
double Ftmp375 = Ftmp374 + 1;
double Ftmp376 = Ftmp289*M[39];
double Ftmp377 = 22*Ftmp8;
double Ftmp378 = Ftmp279*Ftmp377;
double Ftmp379 = Ftmp315 + 3;
double Ftmp380 = Ftmp378 + Ftmp379;
double Ftmp381 = 9*Ftmp9;
double Ftmp382 = 1 - Ftmp381;
double Ftmp383 = Ftmp313 + Ftmp323 + Ftmp382;
double Ftmp384 = Ftmp287*Ftmp377;
double Ftmp385 = Ftmp320 + Ftmp327 + Ftmp382;
double Ftmp386 = Ftmp313 + Ftmp327 + 11*Ftmp352 + 1;
double Ftmp387 = Ftmp142 + 3;
double Ftmp388 = -16*Ftmp39;
double Ftmp389 = 26*Ftmp296;
double Ftmp390 = 20790*Ftmp231;
double Ftmp391 = Ftmp366*M[94];
double Ftmp392 = 16*Ftmp42;
double Ftmp393 = -Ftmp392;
double Ftmp394 = 26*Ftmp306;
double Ftmp395 = Ftmp368*M[96];
double Ftmp396 = 4*Ftmp39;
double Ftmp397 = Ftmp138 + Ftmp396*Ftmp49;
double Ftmp398 = Ftmp397*M[71];
double Ftmp399 = 99*Ftmp153;
double Ftmp400 = Ftmp117*Ftmp396 - 90*Ftmp39 + Ftmp399 + 15;
double Ftmp401 = 36*Ftmp39;
double Ftmp402 = -Ftmp12*Ftmp401;
double Ftmp403 = Ftmp11*Ftmp157;
double Ftmp404 = 2*Ftmp12;
double Ftmp405 = 6*Ftmp1;
double Ftmp406 = -Ftmp139*Ftmp405 + Ftmp404;
double Ftmp407 = 2*Ftmp11;
double Ftmp408 = -Ftmp132*Ftmp405 + Ftmp407;
double Ftmp409 = -44*Ftmp241;
double Ftmp410 = 22*Ftmp1;
double Ftmp411 = -Ftmp139*Ftmp410;
double Ftmp412 = 6*Ftmp12 + Ftmp411;
double Ftmp413 = 6*Ftmp11;
double Ftmp414 = -Ftmp132*Ftmp410;
double Ftmp415 = Ftmp413 + Ftmp414;
double Ftmp416 = Ftmp400*M[105];
double Ftmp417 = 15*Ftmp12;
double Ftmp418 = -220*Ftmp241;
double Ftmp419 = Ftmp12*Ftmp153;
double Ftmp420 = Ftmp332 + Ftmp414;
double Ftmp421 = 15*Ftmp11;
double Ftmp422 = Ftmp343 + Ftmp411;
double Ftmp423 = Ftmp105*Ftmp242;
double Ftmp424 = 66*Ftmp8;
double Ftmp425 = -Ftmp279*Ftmp424;
double Ftmp426 = Ftmp132*Ftmp8;
double Ftmp427 = 143*Ftmp15;
double Ftmp428 = Ftmp426*Ftmp427;
double Ftmp429 = -Ftmp164;
double Ftmp430 = Ftmp213 + Ftmp429;
double Ftmp431 = Ftmp425 + Ftmp428 + Ftmp430 + Ftmp46;
double Ftmp432 = -Ftmp287*Ftmp424;
double Ftmp433 = Ftmp139*Ftmp8;
double Ftmp434 = Ftmp427*Ftmp433;
double Ftmp435 = -Ftmp169;
double Ftmp436 = Ftmp167 + Ftmp435;
double Ftmp437 = Ftmp432 + Ftmp434 + Ftmp436 + Ftmp46;
double Ftmp438 = 45*Ftmp39;
double Ftmp439 = 330*Ftmp8;
double Ftmp440 = -Ftmp279*Ftmp439 - 5;
double Ftmp441 = -Ftmp161;
double Ftmp442 = Ftmp128*Ftmp246;
double Ftmp443 = Ftmp11*Ftmp442;
double Ftmp444 = Ftmp441 + Ftmp443;
double Ftmp445 = Ftmp148 + Ftmp438 + Ftmp440 + Ftmp444;
double Ftmp446 = -Ftmp287*Ftmp439;
double Ftmp447 = Ftmp12*Ftmp442;
double Ftmp448 = Ftmp441 + Ftmp447;
double Ftmp449 = 45*Ftmp42 - 5;
double Ftmp450 = Ftmp148 + Ftmp446 + Ftmp448 + Ftmp449;
double Ftmp451 = -36*Ftmp296;
double Ftmp452 = 99*Ftmp15;
double Ftmp453 = Ftmp426*Ftmp452;
double Ftmp454 = 20*Ftmp39;
double Ftmp455 = -39*Ftmp153 + Ftmp454;
double Ftmp456 = -36*Ftmp306;
double Ftmp457 = Ftmp433*Ftmp452;
double Ftmp458 = -39*Ftmp157 + 20*Ftmp42;
double Ftmp459 = -126*Ftmp296;
double Ftmp460 = -Ftmp134 + Ftmp135;
double Ftmp461 = Ftmp257*Ftmp426 + Ftmp459 + Ftmp460 + Ftmp91;
double Ftmp462 = -126*Ftmp306;
double Ftmp463 = -Ftmp140 + Ftmp141;
double Ftmp464 = Ftmp257*Ftmp433 + Ftmp462 + Ftmp463 + Ftmp91;
double Ftmp465 = 19*Ftmp39;
double Ftmp466 = Ftmp128*Ftmp452;
double Ftmp467 = Ftmp11*Ftmp466;
double Ftmp468 = -102*Ftmp296 - 2;
double Ftmp469 = 8*Ftmp9;
double Ftmp470 = -6*Ftmp129 + Ftmp469;
double Ftmp471 = Ftmp128*Ftmp257;
double Ftmp472 = -Ftmp178 + Ftmp179;
double Ftmp473 = Ftmp11*Ftmp471 + Ftmp459 + Ftmp472 + Ftmp66;
double Ftmp474 = -102*Ftmp306;
double Ftmp475 = Ftmp12*Ftmp466;
double Ftmp476 = 19*Ftmp42 - 2;
double Ftmp477 = Ftmp12*Ftmp471 + Ftmp462 + Ftmp472 + Ftmp62;
double Ftmp478 = 126*Ftmp352;
double Ftmp479 = -Ftmp478;
double Ftmp480 = Ftmp11*Ftmp139;
double Ftmp481 = Ftmp257*Ftmp480 + Ftmp463 + Ftmp479 + Ftmp66;
double Ftmp482 = Ftmp12*Ftmp132;
double Ftmp483 = Ftmp257*Ftmp482 + Ftmp460 + Ftmp479 + Ftmp62;
double Ftmp484 = -44*Ftmp306;
double Ftmp485 = -44*Ftmp296;
double Ftmp486 = Ftmp246*Ftmp426;
double Ftmp487 = 45*Ftmp9;
double Ftmp488 = Ftmp151 + Ftmp429;
double Ftmp489 = Ftmp440 + Ftmp486 + Ftmp487 + Ftmp488;
double Ftmp490 = Ftmp489*M[98];
double Ftmp491 = -220*Ftmp296;
double Ftmp492 = 15*Ftmp9 - 15;
double Ftmp493 = -165*Ftmp153 + 120*Ftmp39;
double Ftmp494 = 1890*Ftmp27;
double Ftmp495 = Ftmp102*Ftmp494;
double Ftmp496 = Ftmp246*Ftmp433;
double Ftmp497 = -198*Ftmp306;
double Ftmp498 = Ftmp381 - 1;
double Ftmp499 = Ftmp436 + Ftmp496 + Ftmp497 + Ftmp498;
double Ftmp500 = Ftmp499*M[102];
double Ftmp501 = Ftmp128*Ftmp427;
double Ftmp502 = Ftmp11*Ftmp501;
double Ftmp503 = Ftmp210 + Ftmp441;
double Ftmp504 = Ftmp425 + Ftmp49 + Ftmp502 + Ftmp503;
double Ftmp505 = Ftmp504*M[89];
double Ftmp506 = 8505*Ftmp27;
double Ftmp507 = Ftmp102*Ftmp506;
double Ftmp508 = 69*Ftmp39;
double Ftmp509 = -418*Ftmp296;
double Ftmp510 = Ftmp508 + Ftmp509;
double Ftmp511 = 84*Ftmp9;
double Ftmp512 = Ftmp443 + Ftmp511;
double Ftmp513 = -66*Ftmp129 - 18;
double Ftmp514 = Ftmp96 - 1;
double Ftmp515 = Ftmp210 + Ftmp448 + Ftmp497 + Ftmp514;
double Ftmp516 = Ftmp515*M[91];
double Ftmp517 = -22*Ftmp129;
double Ftmp518 = Ftmp447 + Ftmp517;
double Ftmp519 = 69*Ftmp42;
double Ftmp520 = -418*Ftmp306;
double Ftmp521 = Ftmp519 + Ftmp520;
double Ftmp522 = 28*Ftmp9 - 6;
double Ftmp523 = -66*Ftmp352;
double Ftmp524 = Ftmp427*Ftmp480 + Ftmp436 + Ftmp49 + Ftmp523;
double Ftmp525 = Ftmp524*M[115];
double Ftmp526 = -330*Ftmp352;
double Ftmp527 = Ftmp246*Ftmp482;
double Ftmp528 = Ftmp449 + Ftmp488 + Ftmp526 + Ftmp527;
double Ftmp529 = Ftmp528*M[113];
double Ftmp530 = -198*Ftmp296;
double Ftmp531 = Ftmp430 + Ftmp486 + Ftmp498 + Ftmp530;
double Ftmp532 = Ftmp531*M[99];
double Ftmp533 = Ftmp155 + Ftmp435 - 5;
double Ftmp534 = Ftmp446 + Ftmp487 + Ftmp496 + Ftmp533;
double Ftmp535 = Ftmp534*M[103];
double Ftmp536 = -220*Ftmp306;
double Ftmp537 = -165*Ftmp157 + 120*Ftmp42;
double Ftmp538 = Ftmp105*Ftmp494;
double Ftmp539 = Ftmp210 + Ftmp444 + Ftmp530 + Ftmp92 - 1;
double Ftmp540 = Ftmp539*M[90];
double Ftmp541 = Ftmp443 + Ftmp517;
double Ftmp542 = Ftmp12*Ftmp501;
double Ftmp543 = Ftmp432 + Ftmp503 + Ftmp52 + Ftmp542;
double Ftmp544 = Ftmp543*M[92];
double Ftmp545 = Ftmp105*Ftmp506;
double Ftmp546 = Ftmp447 + Ftmp511;
double Ftmp547 = Ftmp246*Ftmp480;
double Ftmp548 = Ftmp438 + Ftmp526 + Ftmp533 + Ftmp547;
double Ftmp549 = Ftmp548*M[116];
double Ftmp550 = Ftmp427*Ftmp482 + Ftmp430 + Ftmp52 + Ftmp523;
double Ftmp551 = Ftmp550*M[114];
double Ftmp552 = 4*Ftmp12;
double Ftmp553 = Ftmp1*Ftmp552;
double Ftmp554 = Ftmp553 - 1;
double Ftmp555 = Ftmp407*Ftmp554;
double Ftmp556 = Ftmp1*Ftmp555;
double Ftmp557 = 20*Ftmp352 + Ftmp41*Ftmp62 + Ftmp556;
double Ftmp558 = Ftmp431*M[93];
double Ftmp559 = 8505*Ftmp98;
double Ftmp560 = -209*Ftmp153;
double Ftmp561 = 132*Ftmp8;
double Ftmp562 = -Ftmp279*Ftmp561;
double Ftmp563 = Ftmp10 - 3;
double Ftmp564 = 84*Ftmp39;
double Ftmp565 = Ftmp486 + Ftmp564;
double Ftmp566 = Ftmp437*M[97];
double Ftmp567 = -Ftmp287*Ftmp561;
double Ftmp568 = 84*Ftmp42;
double Ftmp569 = -209*Ftmp157 + Ftmp568;
double Ftmp570 = Ftmp445*M[86];
double Ftmp571 = 125*Ftmp39;
double Ftmp572 = 32*Ftmp9;
double Ftmp573 = 506*Ftmp8;
double Ftmp574 = -Ftmp279*Ftmp573 - 10;
double Ftmp575 = Ftmp450*M[88];
double Ftmp576 = -Ftmp287*Ftmp573;
double Ftmp577 = 125*Ftmp42;
double Ftmp578 = Ftmp577 - 10;
double Ftmp579 = 28*Ftmp352 + Ftmp514*Ftmp59 + Ftmp556;
double Ftmp580 = Ftmp557*M[46];
double Ftmp581 = 5*Ftmp62;
double Ftmp582 = 7*Ftmp41;
double Ftmp583 = Ftmp11*Ftmp514;
double Ftmp584 = 9*Ftmp12;
double Ftmp585 = Ftmp12*Ftmp28;
double Ftmp586 = Ftmp246*Ftmp585;
double Ftmp587 = -Ftmp323 + Ftmp514 + Ftmp586;
double Ftmp588 = -Ftmp320 + Ftmp92;
double Ftmp589 = Ftmp10 - 99*Ftmp352 + Ftmp587 + Ftmp588;
double Ftmp590 = Ftmp15*Ftmp585;
double Ftmp591 = 297*Ftmp590;
double Ftmp592 = -Ftmp331 + Ftmp344 + Ftmp591;
double Ftmp593 = Ftmp333 - Ftmp346;
double Ftmp594 = -Ftmp281 - Ftmp288 - Ftmp353 + Ftmp55 + Ftmp58 + 693*Ftmp590 + Ftmp62;
double Ftmp595 = -Ftmp378 + Ftmp586;
double Ftmp596 = Ftmp141 + Ftmp595;
double Ftmp597 = -Ftmp384;
double Ftmp598 = Ftmp135 + Ftmp597;
double Ftmp599 = -33*Ftmp352 + Ftmp381;
double Ftmp600 = -99*Ftmp306 + Ftmp48 + Ftmp587 + Ftmp599;
double Ftmp601 = Ftmp600*M[100];
double Ftmp602 = -165*Ftmp352 + Ftmp563;
double Ftmp603 = 36*Ftmp42;
double Ftmp604 = Ftmp595 + Ftmp603;
double Ftmp605 = Ftmp369 + Ftmp432;
double Ftmp606 = -99*Ftmp296 + Ftmp52 + Ftmp586 + Ftmp588 + Ftmp599;
double Ftmp607 = Ftmp606*M[101];
double Ftmp608 = Ftmp373 + Ftmp425 + Ftmp586;
double Ftmp609 = Ftmp401 + Ftmp597;
double Ftmp610 = Ftmp589*M[95];
double Ftmp611 = Ftmp1*Ftmp404 - 1;
double Ftmp612 = Ftmp294*Ftmp611;
double Ftmp613 = Ftmp333*Ftmp554 + Ftmp333*Ftmp62 + Ftmp612;
double Ftmp614 = Ftmp52*Ftmp582 + Ftmp613;
double Ftmp615 = 45*Ftmp15*Ftmp614;
double Ftmp616 = Ftmp135*Ftmp514;
double Ftmp617 = Ftmp135*Ftmp554 + Ftmp612 + Ftmp616;
double Ftmp618 = 3*Ftmp121*Ftmp59 + Ftmp617;
double Ftmp619 = 3465*Ftmp618*M[107];
double Ftmp620 = -Ftmp553;
double Ftmp621 = Ftmp620 + 1;
double Ftmp622 = 18*Ftmp352;
double Ftmp623 = Ftmp313 + Ftmp622;
double Ftmp624 = 4*Ftmp621 + 4*Ftmp623;
double Ftmp625 = Ftmp39*Ftmp624;
double Ftmp626 = Ftmp138*Ftmp514 + 56*Ftmp352*Ftmp49 + Ftmp625;
double Ftmp627 = Ftmp626*M[106];
double Ftmp628 = 8*Ftmp11;
double Ftmp629 = Ftmp611*Ftmp628;
double Ftmp630 = 21*Ftmp41;
double Ftmp631 = Ftmp11*Ftmp554;
double Ftmp632 = 14*Ftmp583;
double Ftmp633 = 21*Ftmp121;
double Ftmp634 = 33*Ftmp59;
double Ftmp635 = z*M[107];
double Ftmp636 = 160*Ftmp611;
double Ftmp637 = 280*Ftmp52;
double Ftmp638 = 16*Ftmp157 + Ftmp375;
double Ftmp639 = Ftmp294*Ftmp638 + Ftmp454*Ftmp554*Ftmp62;
double Ftmp640 = Ftmp170*Ftmp582 + Ftmp352*Ftmp636 + Ftmp352*Ftmp637 + Ftmp639;
double Ftmp641 = Ftmp640*M[108];
double Ftmp642 = 112*Ftmp49;
double Ftmp643 = 3 - Ftmp304;
double Ftmp644 = 32*Ftmp39;
double Ftmp645 = 320*Ftmp611;
double Ftmp646 = 560*Ftmp52;
double Ftmp647 = 140*Ftmp554;
double Ftmp648 = 80*Ftmp62;
double Ftmp649 = 20*Ftmp62;
double Ftmp650 = 42*Ftmp41;
double Ftmp651 = -Ftmp121*Ftmp650;
double Ftmp652 = Ftmp16*y;
double Ftmp653 = pow(x, 3);
double Ftmp654 = Ftmp653*M[9];
double Ftmp655 = Ftmp31*M[11];
double Ftmp656 = Ftmp11*Ftmp16;
double Ftmp657 = Ftmp11*Ftmp242;
double Ftmp658 = Ftmp657*z;
double Ftmp659 = Ftmp39 - 1;
double Ftmp660 = Ftmp53*x;
double Ftmp661 = Ftmp91*M[23];
double Ftmp662 = Ftmp59*M[30];
double Ftmp663 = 30*Ftmp71;
double Ftmp664 = 30*Ftmp659;
double Ftmp665 = Ftmp7*y;
double Ftmp666 = Ftmp54*z;
double Ftmp667 = Ftmp421*Ftmp7;
double Ftmp668 = Ftmp653*y;
double Ftmp669 = Ftmp203*Ftmp27*(Ftmp381 - 5);
double Ftmp670 = 1890*Ftmp106;
double Ftmp671 = 2835*Ftmp106;
double Ftmp672 = Ftmp105*Ftmp671;
double Ftmp673 = 1260*Ftmp46;
double Ftmp674 = Ftmp673*Ftmp98*M[36];
double Ftmp675 = 11*Ftmp9;
double Ftmp676 = Ftmp675 - 5;
double Ftmp677 = Ftmp673*M[35];
double Ftmp678 = Ftmp675 - 3;
double Ftmp679 = Ftmp11*Ftmp123;
double Ftmp680 = Ftmp653*z;
double Ftmp681 = Ftmp232*(13*Ftmp9 - 5);
double Ftmp682 = Ftmp159*x;
double Ftmp683 = Ftmp12 + Ftmp8;
double Ftmp684 = Ftmp212*M[59];
double Ftmp685 = Ftmp165*M[77];
double Ftmp686 = 3*Ftmp153 - Ftmp396 + 1;
double Ftmp687 = Ftmp176*Ftmp686;
double Ftmp688 = Ftmp160*z;
double Ftmp689 = Ftmp100*Ftmp188;
double Ftmp690 = Ftmp11*Ftmp159;
double Ftmp691 = 420*Ftmp15;
double Ftmp692 = Ftmp11*Ftmp188;
double Ftmp693 = -Ftmp12*Ftmp210;
double Ftmp694 = Ftmp50 + Ftmp8;
double Ftmp695 = 11*Ftmp153 + Ftmp388 + 5;
double Ftmp696 = 3*Ftmp8;
double Ftmp697 = Ftmp12 + Ftmp696;
double Ftmp698 = Ftmp290*Ftmp32*(Ftmp199 - 154*Ftmp9 + 35);
double Ftmp699 = Ftmp11*Ftmp32;
double Ftmp700 = Ftmp105*Ftmp699;
double Ftmp701 = 10395*Ftmp700;
double Ftmp702 = Ftmp200*Ftmp239*M[85];
double Ftmp703 = Ftmp12*Ftmp9;
double Ftmp704 = -22*Ftmp703;
double Ftmp705 = Ftmp200*M[84];
double Ftmp706 = Ftmp11*Ftmp237;
double Ftmp707 = Ftmp50 + Ftmp696;
double Ftmp708 = 35*Ftmp39;
double Ftmp709 = Ftmp15*Ftmp249;
double Ftmp710 = 143*Ftmp709;
double Ftmp711 = -253*Ftmp153 + Ftmp571 + Ftmp710 - 15;
double Ftmp712 = Ftmp100*Ftmp266;
double Ftmp713 = 5670*Ftmp100*Ftmp27;
double Ftmp714 = -Ftmp469;
double Ftmp715 = -Ftmp39;
double Ftmp716 = Ftmp715 + 1;
double Ftmp717 = 14*Ftmp352;
double Ftmp718 = Ftmp311*x;
double Ftmp719 = Ftmp383*M[66];
double Ftmp720 = Ftmp385*M[68];
double Ftmp721 = Ftmp386*M[79];
double Ftmp722 = Ftmp331 + Ftmp716;
double Ftmp723 = Ftmp100*Ftmp336;
double Ftmp724 = 12*Ftmp9;
double Ftmp725 = -Ftmp724;
double Ftmp726 = Ftmp11*Ftmp336;
double Ftmp727 = Ftmp313 + 3;
double Ftmp728 = 22*Ftmp352;
double Ftmp729 = 36*Ftmp9;
double Ftmp730 = -Ftmp12*Ftmp729;
double Ftmp731 = Ftmp157*Ftmp8;
double Ftmp732 = Ftmp12*Ftmp129;
double Ftmp733 = -Ftmp128*Ftmp405 + 2*Ftmp8;
double Ftmp734 = -44*Ftmp703;
double Ftmp735 = -Ftmp128*Ftmp410;
double Ftmp736 = Ftmp735 + 6*Ftmp8;
double Ftmp737 = Ftmp100*Ftmp237;
double Ftmp738 = 429*Ftmp732;
double Ftmp739 = -220*Ftmp703;
double Ftmp740 = Ftmp735 + 10*Ftmp8;
double Ftmp741 = 429*Ftmp731;
double Ftmp742 = Ftmp100*Ftmp242;
double Ftmp743 = -132*Ftmp703;
double Ftmp744 = -39*Ftmp129 + 20*Ftmp9;
double Ftmp745 = 36*Ftmp352;
double Ftmp746 = -Ftmp745;
double Ftmp747 = Ftmp452*Ftmp480;
double Ftmp748 = 19*Ftmp9;
double Ftmp749 = -6*Ftmp153 + Ftmp294;
double Ftmp750 = -102*Ftmp352;
double Ftmp751 = Ftmp452*Ftmp482;
double Ftmp752 = -55*Ftmp129 + 24*Ftmp9;
double Ftmp753 = -66*Ftmp153 - 18;
double Ftmp754 = 69*Ftmp9;
double Ftmp755 = Ftmp509 + Ftmp754;
double Ftmp756 = Ftmp1*Ftmp421 - 15;
double Ftmp757 = -165*Ftmp129 + 120*Ftmp9;
double Ftmp758 = -22*Ftmp153;
double Ftmp759 = 28*Ftmp39;
double Ftmp760 = Ftmp100*Ftmp494;
double Ftmp761 = Ftmp100*Ftmp506;
double Ftmp762 = -220*Ftmp352;
double Ftmp763 = -418*Ftmp352;
double Ftmp764 = Ftmp527 + Ftmp564;
double Ftmp765 = Ftmp644 + Ftmp758;
double Ftmp766 = 8505*Ftmp106;
double Ftmp767 = -209*Ftmp129;
double Ftmp768 = Ftmp48 - 3;
double Ftmp769 = -132*Ftmp352;
double Ftmp770 = -506*Ftmp352;
double Ftmp771 = -Ftmp659;
double Ftmp772 = Ftmp1*Ftmp407;
double Ftmp773 = Ftmp145 - Ftmp622;
double Ftmp774 = -Ftmp728;
double Ftmp775 = Ftmp179 + Ftmp774;
double Ftmp776 = -165*Ftmp306 + Ftmp768;
double Ftmp777 = Ftmp523 + Ftmp724;
double Ftmp778 = Ftmp729 + Ftmp774;
double Ftmp779 = 28*Ftmp42;
double Ftmp780 = Ftmp556 + 2;
double Ftmp781 = 35*Ftmp771;
double Ftmp782 = 630*x;
double Ftmp783 = Ftmp102*Ftmp27;
double Ftmp784 = 405*Ftmp641;
double Ftmp785 = Ftmp352*Ftmp642 + Ftmp625;
double Ftmp786 = Ftmp352*Ftmp645 + Ftmp352*Ftmp646 + Ftmp639;
double Ftmp787 = Ftmp16*z;
double Ftmp788 = Ftmp12*Ftmp242;
double Ftmp789 = Ftmp42 - 1;
double Ftmp790 = Ftmp789*z;
double Ftmp791 = Ftmp417*Ftmp7;
double Ftmp792 = 1890*Ftmp101;
double Ftmp793 = 2835*Ftmp101;
double Ftmp794 = Ftmp102*Ftmp793;
double Ftmp795 = Ftmp12*Ftmp123;
double Ftmp796 = Ftmp11 + Ftmp8;
double Ftmp797 = Ftmp176*(3*Ftmp157 + Ftmp621);
double Ftmp798 = Ftmp12*Ftmp159;
double Ftmp799 = Ftmp12*Ftmp188;
double Ftmp800 = -Ftmp11*Ftmp210;
double Ftmp801 = 11*Ftmp157;
double Ftmp802 = Ftmp387 + Ftmp801;
double Ftmp803 = Ftmp47 + Ftmp8;
double Ftmp804 = Ftmp393 + Ftmp801 + 5;
double Ftmp805 = Ftmp11 + Ftmp696;
double Ftmp806 = Ftmp102*Ftmp32;
double Ftmp807 = Ftmp12*Ftmp806;
double Ftmp808 = 10395*Ftmp807;
double Ftmp809 = Ftmp11*Ftmp9;
double Ftmp810 = -22*Ftmp809;
double Ftmp811 = Ftmp47 + Ftmp696;
double Ftmp812 = Ftmp15*Ftmp252;
double Ftmp813 = 143*Ftmp812;
double Ftmp814 = -253*Ftmp157 + Ftmp577 + Ftmp813 - 15;
double Ftmp815 = -Ftmp42;
double Ftmp816 = Ftmp815 + 1;
double Ftmp817 = Ftmp346 + Ftmp816;
double Ftmp818 = Ftmp12*Ftmp336;
double Ftmp819 = Ftmp327 + Ftmp384;
double Ftmp820 = Ftmp327 + 3;
double Ftmp821 = -Ftmp11*Ftmp729;
double Ftmp822 = Ftmp11*Ftmp129;
double Ftmp823 = -44*Ftmp809;
double Ftmp824 = 429*Ftmp822;
double Ftmp825 = -220*Ftmp809;
double Ftmp826 = 429*Ftmp153*Ftmp8;
double Ftmp827 = -132*Ftmp809;
double Ftmp828 = -6*Ftmp157 + Ftmp304 - 2;
double Ftmp829 = Ftmp520 + Ftmp754;
double Ftmp830 = -66*Ftmp157 + Ftmp568 - 18;
double Ftmp831 = Ftmp1*Ftmp417 - 15;
double Ftmp832 = -22*Ftmp157;
double Ftmp833 = Ftmp496 + Ftmp832;
double Ftmp834 = 32*Ftmp42 - 10;
double Ftmp835 = 8505*Ftmp101;
double Ftmp836 = Ftmp51 - 3;
double Ftmp837 = -165*Ftmp296 + Ftmp586 + Ftmp836;
double Ftmp838 = -Ftmp789;
double Ftmp839 = -Ftmp294*Ftmp838 + Ftmp556;
double Ftmp840 = Ftmp39*Ftmp838;
double Ftmp841 = Ftmp39*Ftmp789;
#pragma omp atomic
F[0] += Ftmp0*(3*Ftmp1*Ftmp38*M[9] + 3*Ftmp1*Ftmp41*M[12] + 3*Ftmp1*Ftmp44*M[14] + 6*Ftmp1*Ftmp45*x*M[3] + 3*Ftmp1*Ftmp46*x*M[3] + 3*Ftmp1*Ftmp49*x*M[6] + 3*Ftmp1*Ftmp52*x*M[8] - Ftmp10*M[0] - Ftmp100*Ftmp45*Ftmp99*M[38] - Ftmp102*Ftmp177*M[35] - Ftmp104*M[53] - Ftmp105*Ftmp108 - Ftmp105*Ftmp177*M[36] - Ftmp109*Ftmp111 + 3780*Ftmp11*Ftmp122*Ftmp32*Ftmp8*z*M[71] + 60*Ftmp11*Ftmp15*Ftmp59*x*M[29] + 1890*Ftmp11*Ftmp165*Ftmp27*x*M[76] + 210*Ftmp11*Ftmp27*Ftmp8*z*M[26] + 6*Ftmp11*Ftmp7*x*M[6] - Ftmp110*Ftmp291 - Ftmp110*Ftmp292 - Ftmp110*Ftmp293 - Ftmp110*Ftmp570 - Ftmp110*Ftmp575 - Ftmp110*Ftmp610 - Ftmp110*Ftmp627 - Ftmp111*Ftmp112 - Ftmp111*Ftmp113 + 3780*Ftmp115*Ftmp19*Ftmp32*Ftmp8*M[75] + 3780*Ftmp115*Ftmp19*Ftmp32*x*y*M[81] + 3780*Ftmp117*Ftmp13*Ftmp32*Ftmp8*M[70] + 3780*Ftmp117*Ftmp13*Ftmp32*x*z*M[77] - Ftmp118*Ftmp119 + 3780*Ftmp12*Ftmp121*Ftmp32*Ftmp8*y*M[74] + 60*Ftmp12*Ftmp15*Ftmp68*x*M[33] + 1890*Ftmp12*Ftmp174*Ftmp27*x*M[82] + 210*Ftmp12*Ftmp27*Ftmp8*y*M[27] + 6*Ftmp12*Ftmp7*x*M[8] - Ftmp120*Ftmp98*M[48] - Ftmp124*Ftmp126*Ftmp35 - Ftmp124*Ftmp127*Ftmp36 - 630*Ftmp124*Ftmp635*y*(Ftmp11*Ftmp633 + Ftmp12*Ftmp634 + 198*Ftmp241 + Ftmp629 + 14*Ftmp631 + Ftmp632) + 210*Ftmp13*Ftmp27*Ftmp8*M[25] + 210*Ftmp13*Ftmp27*x*z*M[30] + 15*Ftmp131*Ftmp7*M[34] + 45*Ftmp138*Ftmp7*M[44] - Ftmp14*Ftmp17 - Ftmp14*Ftmp219*Ftmp32 - Ftmp14*Ftmp94 + 45*Ftmp144*Ftmp7*M[48] + 315*Ftmp15*Ftmp248*M[83] + 315*Ftmp15*Ftmp251*M[104] + 315*Ftmp15*Ftmp254*M[110] + 315*Ftmp15*Ftmp256*x*M[55] + 315*Ftmp15*Ftmp258*x*M[76] + 315*Ftmp15*Ftmp259*x*M[82] + 945*Ftmp15*Ftmp431*M[93] + 945*Ftmp15*Ftmp437*M[97] + 315*Ftmp15*Ftmp445*M[86] + 210*Ftmp15*Ftmp45*Ftmp8*y*M[20] + 210*Ftmp15*Ftmp45*Ftmp8*z*M[21] + 210*Ftmp15*Ftmp45*x*y*z*M[23] + 315*Ftmp15*Ftmp450*M[88] + 315*Ftmp15*Ftmp46*y*z*M[38] + 315*Ftmp15*Ftmp461*x*M[65] + 315*Ftmp15*Ftmp464*x*M[69] + 315*Ftmp15*Ftmp473*x*M[58] + 315*Ftmp15*Ftmp477*x*M[60] + 315*Ftmp15*Ftmp481*x*M[80] + 315*Ftmp15*Ftmp483*x*M[78] + 315*Ftmp15*Ftmp49*y*z*M[45] + 315*Ftmp15*Ftmp52*y*z*M[47] + 105*Ftmp15*Ftmp56*Ftmp8*y*M[20] + 105*Ftmp15*Ftmp56*Ftmp8*z*M[21] + 315*Ftmp15*Ftmp589*M[95] + 105*Ftmp15*Ftmp59*Ftmp8*y*M[25] + 105*Ftmp15*Ftmp59*x*y*z*M[30] + 105*Ftmp15*Ftmp594*x*M[67] + 105*Ftmp15*Ftmp62*Ftmp8*y*M[27] + 315*Ftmp15*Ftmp626*M[106] + 45*Ftmp15*Ftmp640*M[108] + 105*Ftmp15*Ftmp66*Ftmp8*z*M[26] + 105*Ftmp15*Ftmp68*Ftmp8*z*M[28] + 105*Ftmp15*Ftmp68*x*y*z*M[32] + 105*Ftmp15*Ftmp91*x*y*z*M[23] + 1890*Ftmp15*x*(-Ftmp130 + 33*Ftmp255 + 35*Ftmp9 - 5)*M[55] + 210*Ftmp15*x*(-117*Ftmp352 + Ftmp45 + Ftmp592 + Ftmp593)*M[67] + 630*Ftmp15*x*(Ftmp45 + Ftmp451 + Ftmp453 + Ftmp455)*M[65] + 630*Ftmp15*x*(Ftmp45 + Ftmp456 + Ftmp457 + Ftmp458)*M[69] + 630*Ftmp15*x*(Ftmp465 + Ftmp467 + Ftmp468 + Ftmp470)*M[58] + 630*Ftmp15*x*(Ftmp470 + Ftmp474 + Ftmp475 + Ftmp476)*M[60] + 15*Ftmp150*Ftmp7*x*M[19] + 15*Ftmp154*Ftmp7*x*M[29] + 15*Ftmp158*Ftmp7*x*M[33] - Ftmp160*Ftmp163 - Ftmp160*Ftmp166 - Ftmp160*Ftmp171 - Ftmp160*Ftmp322 + 2835*Ftmp162*Ftmp27*Ftmp8*y*M[56] + 2835*Ftmp162*Ftmp27*Ftmp8*z*M[57] + 2835*Ftmp165*Ftmp27*Ftmp8*y*M[70] + 2835*Ftmp165*Ftmp27*x*y*z*M[77] + 2835*Ftmp170*Ftmp27*Ftmp8*y*M[74] - Ftmp172*Ftmp173 - Ftmp172*Ftmp325 - Ftmp172*Ftmp398 + 2835*Ftmp174*Ftmp27*Ftmp8*z*M[75] + 2835*Ftmp174*Ftmp27*x*y*z*M[81] - Ftmp175*M[75] - Ftmp18*Ftmp20 - Ftmp18*(Ftmp11 - Ftmp12*Ftmp135 + Ftmp12)*M[31] - Ftmp183*Ftmp184 - Ftmp184*Ftmp191 - Ftmp187*Ftmp190 - Ftmp188*Ftmp579*y*M[72] + 210*Ftmp19*Ftmp27*Ftmp8*M[28] + 210*Ftmp19*Ftmp27*x*y*M[32] - Ftmp190*Ftmp342 - Ftmp190*Ftmp351 - Ftmp190*Ftmp355 - Ftmp192*Ftmp193 - Ftmp193*Ftmp194 - Ftmp197*Ftmp198 - Ftmp198*Ftmp359 - Ftmp198*Ftmp362 - Ftmp198*Ftmp364 - Ftmp2*y*M[4] + 945*Ftmp200*Ftmp27*y*z*M[87] + 945*Ftmp202*Ftmp27*y*z*M[109] - Ftmp203*Ftmp25*(9*Ftmp129 + Ftmp180 + 5) - Ftmp204*Ftmp26 - Ftmp205*Ftmp206 - Ftmp206*Ftmp207 - Ftmp21*Ftmp22 + 2835*Ftmp212*Ftmp27*x*y*z*M[59] + 3780*Ftmp216*Ftmp27*Ftmp8*y*M[56] + 3780*Ftmp216*Ftmp27*Ftmp8*z*M[57] - Ftmp221*Ftmp95 - Ftmp223*Ftmp225 - Ftmp223*(Ftmp11*Ftmp201 + Ftmp11 + Ftmp409 + Ftmp412)*M[115] - Ftmp227*Ftmp230 - Ftmp229*(Ftmp12*Ftmp217 + Ftmp12 + Ftmp409 + Ftmp415)*M[114] - Ftmp23*Ftmp24 - Ftmp231*Ftmp232*(13*Ftmp129 + Ftmp211 + 5) - Ftmp231*Ftmp619 - Ftmp233*Ftmp234 - Ftmp234*Ftmp235 - Ftmp234*Ftmp391 - Ftmp234*Ftmp395 - Ftmp234*Ftmp416 - Ftmp236*Ftmp238 - Ftmp239*Ftmp240 - Ftmp243*(Ftmp417 + Ftmp418 + 429*Ftmp419 + Ftmp420)*M[113] - Ftmp244*(-22*Ftmp241 + Ftmp47 + Ftmp50)*M[79] - Ftmp26*Ftmp372 - Ftmp26*Ftmp376 - Ftmp26*Ftmp580 - Ftmp26*y*z*M[13] - Ftmp262*Ftmp263*M[84] - Ftmp262*Ftmp271*M[85] - Ftmp263*(-55*Ftmp157 + 24*Ftmp42 + Ftmp434 + Ftmp45 + Ftmp484)*M[102] - Ftmp265*Ftmp267 - Ftmp267*Ftmp269 - Ftmp267*Ftmp270 - Ftmp267*Ftmp490 - Ftmp267*Ftmp500 - Ftmp267*Ftmp516 - Ftmp267*Ftmp529 - Ftmp267*Ftmp601 + 8505*Ftmp27*Ftmp317*Ftmp8*y*M[61] + 2835*Ftmp27*Ftmp321*Ftmp8*y*M[63] + 2835*Ftmp27*Ftmp324*Ftmp8*z*M[62] + 8505*Ftmp27*Ftmp328*Ftmp8*z*M[64] + 945*Ftmp27*Ftmp366*y*z*M[94] + 945*Ftmp27*Ftmp368*y*z*M[96] + 2835*Ftmp27*Ftmp383*x*y*z*M[66] + 2835*Ftmp27*Ftmp385*x*y*z*M[68] + 8505*Ftmp27*Ftmp386*x*y*z*M[79] + 2835*Ftmp27*Ftmp397*Ftmp8*z*M[71] + 945*Ftmp27*Ftmp400*y*z*M[105] + 945*Ftmp27*Ftmp579*Ftmp8*y*M[72] + 405*Ftmp27*Ftmp614*Ftmp8*z*M[73] + 315*Ftmp27*Ftmp618*y*z*M[107] + 1890*Ftmp27*Ftmp8*y*(Ftmp136 + Ftmp380)*M[61] + 1890*Ftmp27*Ftmp8*y*(Ftmp143 + Ftmp298 + Ftmp384)*M[63] + 1890*Ftmp27*Ftmp8*z*(Ftmp137 + Ftmp298 + Ftmp378)*M[62] + 1890*Ftmp27*Ftmp8*z*(Ftmp315 + Ftmp384 + Ftmp387)*M[64] + 3780*Ftmp27*x*y*z*(Ftmp208 + Ftmp209)*M[59] + 1890*Ftmp27*x*y*z*(Ftmp370 + Ftmp380)*M[66] + 1890*Ftmp27*x*y*z*(Ftmp374 + Ftmp379 + Ftmp384)*M[68] + 210*Ftmp27*x*y*(Ftmp11 + Ftmp214 + Ftmp50)*M[51] + 210*Ftmp27*x*z*(Ftmp12 + Ftmp214 + Ftmp47)*M[52] + 630*Ftmp27*x*(Ftmp11 + Ftmp402 + 99*Ftmp403 + Ftmp406)*M[80] + 630*Ftmp27*x*(Ftmp12*Ftmp399 + Ftmp12 + Ftmp402 + Ftmp408)*M[78] - Ftmp271*(-55*Ftmp153 + 24*Ftmp39 + Ftmp428 + Ftmp45 + Ftmp485)*M[99] - Ftmp272*Ftmp273 - Ftmp273*Ftmp274 - Ftmp273*Ftmp276 - Ftmp273*Ftmp532 - Ftmp273*Ftmp535 - Ftmp273*Ftmp540 - Ftmp273*Ftmp549 - Ftmp273*Ftmp607 + 15*Ftmp284*Ftmp7*M[37] + 15*Ftmp289*Ftmp7*M[39] - Ftmp29*M[12] - Ftmp290*Ftmp98*(-297*Ftmp129 + Ftmp261 + 189*Ftmp9 - 35) - Ftmp3*M[5] - Ftmp30*Ftmp31 + 15*Ftmp303*Ftmp7*x*M[22] + 15*Ftmp309*Ftmp7*x*M[24] + 15*Ftmp310*Ftmp7*x*M[31] - Ftmp312*Ftmp318 + 210*Ftmp32*Ftmp8*y*(64*Ftmp241 + Ftmp555 + 7*Ftmp583 + Ftmp584*Ftmp59)*M[72] + 90*Ftmp32*Ftmp8*z*(Ftmp12*Ftmp630 + 126*Ftmp241 + Ftmp301*Ftmp52 + Ftmp332*Ftmp554 + Ftmp332*Ftmp62 + Ftmp629)*M[73] - Ftmp326*Ftmp329 - 630*Ftmp33*(Ftmp11*Ftmp624 + Ftmp138*Ftmp584 + Ftmp241*Ftmp642 - Ftmp396*(Ftmp402 + Ftmp47 + Ftmp552) + 168*Ftmp419 + Ftmp49*Ftmp632)*M[106] - 90*Ftmp33*(-Ftmp12*Ftmp643*Ftmp644 - Ftmp12*Ftmp651 + Ftmp170*Ftmp301 + Ftmp241*Ftmp645 + Ftmp241*Ftmp646 + Ftmp241*Ftmp647 + Ftmp241*Ftmp648 + 1160*Ftmp403 + Ftmp628*Ftmp638 + Ftmp631*Ftmp649)*M[108] - Ftmp337*(Ftmp299 + Ftmp345 + Ftmp346)*M[42] - Ftmp337*(Ftmp315 + Ftmp331 + Ftmp335)*M[40] - Ftmp34*Ftmp35*M[47] - Ftmp34*Ftmp36*M[45] - Ftmp357*(Ftmp334 + Ftmp356)*M[41] - Ftmp357*(Ftmp315 + Ftmp346 + Ftmp360)*M[43] - Ftmp371*(Ftmp356 + Ftmp370)*M[37] - Ftmp371*(Ftmp298 + Ftmp346 + Ftmp375)*M[39] - Ftmp390*(Ftmp379 + Ftmp388 + Ftmp389)*M[94] - Ftmp390*(Ftmp379 + Ftmp393 + Ftmp394)*M[96] - Ftmp423*(429*Ftmp403 + Ftmp418 + Ftmp421 + Ftmp422)*M[116] - Ftmp495*(Ftmp510 + Ftmp512 + Ftmp513)*M[89] - Ftmp495*(Ftmp518 + Ftmp521 + Ftmp522)*M[91] - Ftmp495*(Ftmp602 + Ftmp604 + Ftmp605)*M[100] - Ftmp495*(Ftmp486 + Ftmp491 + Ftmp492 + Ftmp493)*M[98] - Ftmp5*M[1] - Ftmp505*Ftmp507 - Ftmp507*Ftmp525 - Ftmp538*(Ftmp510 + Ftmp522 + Ftmp541)*M[90] - Ftmp538*(Ftmp513 + Ftmp521 + Ftmp546)*M[92] - Ftmp538*(Ftmp602 + Ftmp608 + Ftmp609)*M[101] - Ftmp538*(Ftmp492 + Ftmp496 + Ftmp536 + Ftmp537)*M[103] - Ftmp54*Ftmp57 - Ftmp54*Ftmp60 - Ftmp54*Ftmp63 - Ftmp544*Ftmp545 - Ftmp545*Ftmp551 + 15*Ftmp557*Ftmp7*M[46] - Ftmp558*Ftmp559 - Ftmp559*Ftmp566 - Ftmp6*x - Ftmp615*z*M[73] - Ftmp64*Ftmp65 - Ftmp64*Ftmp67 - 405*Ftmp641*Ftmp98 - Ftmp69*M[28] + 15*Ftmp7*Ftmp8*y*M[4] + 15*Ftmp7*Ftmp8*z*M[5] + 15*Ftmp7*x*y*z*M[7] + 60*Ftmp7*x*(7*Ftmp129 + Ftmp147)*M[19] + 30*Ftmp7*x*(Ftmp295 + Ftmp297 + Ftmp299)*M[22] + 30*Ftmp7*x*(Ftmp299 + Ftmp305 + Ftmp307)*M[24] + 15*Ftmp7*y*z*M[13] - Ftmp70*Ftmp74*M[9] - Ftmp72*Ftmp73 - Ftmp72*z*M[11] - Ftmp76*Ftmp77 - Ftmp77*Ftmp79 - Ftmp77*Ftmp80 - Ftmp81*Ftmp82 - Ftmp81*Ftmp83 - Ftmp81*Ftmp85 - Ftmp86*Ftmp88 - Ftmp88*Ftmp89 - Ftmp88*Ftmp90 - Ftmp95*Ftmp97 - 30*Ftmp98*(Ftmp11*Ftmp581 + Ftmp12*Ftmp582 + 48*Ftmp241 + Ftmp555)*M[46] - Ftmp99*(-209*Ftmp352 + Ftmp45 + Ftmp596 + Ftmp598)*M[95] - Ftmp99*(Ftmp496 + Ftmp563 + Ftmp567 + Ftmp569)*M[97] - Ftmp99*(Ftmp518 + Ftmp572 + Ftmp576 + Ftmp578)*M[88] - Ftmp99*(Ftmp541 + Ftmp571 + Ftmp572 + Ftmp574)*M[86] - Ftmp99*(Ftmp560 + Ftmp562 + Ftmp563 + Ftmp565)*M[93] + M[0]);
#pragma omp atomic
F[1] += Ftmp0*(3*Ftmp1*Ftmp44*M[17] + 3*Ftmp1*Ftmp46*y*M[3] + 3*Ftmp1*Ftmp49*y*M[6] + 3*Ftmp1*Ftmp52*y*M[8] + 6*Ftmp1*Ftmp659*y*M[6] + 3*Ftmp1*Ftmp75*M[10] + 3*Ftmp1*Ftmp78*M[15] - Ftmp100*Ftmp655 - Ftmp100*Ftmp674 - Ftmp100*Ftmp687*M[50] - Ftmp100*Ftmp702 - Ftmp102*Ftmp687*M[44] - Ftmp104*M[48] - Ftmp105*Ftmp659*Ftmp670*M[45] - Ftmp106*Ftmp120*M[53] - Ftmp106*Ftmp218*(-297*Ftmp153 + 189*Ftmp39 + Ftmp710 - 35) - Ftmp109*Ftmp672 + 3780*Ftmp11*Ftmp115*Ftmp19*Ftmp32*M[81] + 3780*Ftmp11*Ftmp12*Ftmp121*Ftmp32*x*M[74] + 210*Ftmp11*Ftmp12*Ftmp27*x*M[27] + 20790*Ftmp11*Ftmp123*x*z*(-26*Ftmp703 + Ftmp707)*M[96] + 105*Ftmp11*Ftmp15*Ftmp56*x*M[20] + 105*Ftmp11*Ftmp15*Ftmp59*x*M[25] + 105*Ftmp11*Ftmp15*Ftmp59*z*M[30] + 105*Ftmp11*Ftmp15*Ftmp62*x*M[27] + 210*Ftmp11*Ftmp15*Ftmp659*x*M[25] + 210*Ftmp11*Ftmp15*Ftmp659*z*M[30] + 105*Ftmp11*Ftmp15*Ftmp68*z*M[32] + 105*Ftmp11*Ftmp15*Ftmp91*z*M[23] + 2835*Ftmp11*Ftmp162*Ftmp27*x*M[56] + 2835*Ftmp11*Ftmp165*Ftmp27*x*M[70] + 2835*Ftmp11*Ftmp165*Ftmp27*z*M[77] + 2835*Ftmp11*Ftmp170*Ftmp27*x*M[74] + 2835*Ftmp11*Ftmp174*Ftmp27*z*M[81] + 210*Ftmp11*Ftmp19*Ftmp27*M[32] - Ftmp11*Ftmp198*M[13] + 2835*Ftmp11*Ftmp212*Ftmp27*z*M[59] + 8505*Ftmp11*Ftmp27*Ftmp317*x*M[61] + 2835*Ftmp11*Ftmp27*Ftmp321*x*M[63] + 2835*Ftmp11*Ftmp27*Ftmp383*z*M[66] + 2835*Ftmp11*Ftmp27*Ftmp385*z*M[68] + 8505*Ftmp11*Ftmp27*Ftmp386*z*M[79] + 945*Ftmp11*Ftmp27*Ftmp579*x*M[72] + 210*Ftmp11*Ftmp27*Ftmp653*M[20] + 3780*Ftmp11*Ftmp27*Ftmp695*x*M[70] + 3780*Ftmp11*Ftmp27*Ftmp695*z*M[77] + 210*Ftmp11*Ftmp27*Ftmp8*z*M[23] + 1890*Ftmp11*Ftmp27*x*(Ftmp209 + Ftmp313 + Ftmp378)*M[61] + 210*Ftmp11*Ftmp27*x*(64*Ftmp352 - 7*Ftmp514*Ftmp771 + Ftmp59*Ftmp96 - Ftmp603 + Ftmp780)*M[72] + 1890*Ftmp11*Ftmp27*z*(Ftmp181 + Ftmp378 + Ftmp715)*M[66] + 1890*Ftmp11*Ftmp27*z*(Ftmp313 + Ftmp387 + Ftmp728)*M[79] + 210*Ftmp11*Ftmp27*(Ftmp683 + Ftmp693)*M[42] + 3780*Ftmp11*Ftmp32*Ftmp653*Ftmp676*M[56] + 3780*Ftmp11*Ftmp32*Ftmp678*Ftmp8*z*M[59] - Ftmp11*Ftmp691*(Ftmp136 + 9*Ftmp153 + 5)*M[49] + 15*Ftmp11*Ftmp7*x*M[4] + 15*Ftmp11*Ftmp7*z*M[7] - Ftmp112*Ftmp672 - Ftmp113*Ftmp672 + 3780*Ftmp115*Ftmp19*Ftmp32*x*y*M[75] - Ftmp118*Ftmp677 + 60*Ftmp12*Ftmp15*Ftmp68*y*M[33] + 1890*Ftmp12*Ftmp174*Ftmp27*y*M[82] + 6*Ftmp12*Ftmp7*y*M[8] - Ftmp126*Ftmp679*Ftmp95 + 45*Ftmp144*Ftmp7*M[53] + 315*Ftmp15*Ftmp254*M[117] + 315*Ftmp15*Ftmp256*y*M[55] + 315*Ftmp15*Ftmp258*y*M[76] + 315*Ftmp15*Ftmp259*y*M[82] + 315*Ftmp15*Ftmp264*M[84] + 315*Ftmp15*Ftmp268*M[111] + 315*Ftmp15*Ftmp46*x*z*M[38] + 315*Ftmp15*Ftmp461*y*M[65] + 315*Ftmp15*Ftmp464*y*M[69] + 315*Ftmp15*Ftmp473*y*M[58] + 315*Ftmp15*Ftmp477*y*M[60] + 315*Ftmp15*Ftmp481*y*M[80] + 315*Ftmp15*Ftmp483*y*M[78] + 315*Ftmp15*Ftmp489*M[98] + 315*Ftmp15*Ftmp49*x*z*M[45] + 315*Ftmp15*Ftmp499*M[102] + 945*Ftmp15*Ftmp504*M[89] + 315*Ftmp15*Ftmp515*M[91] + 315*Ftmp15*Ftmp52*x*z*M[47] + 945*Ftmp15*Ftmp524*M[115] + 315*Ftmp15*Ftmp528*M[113] + 60*Ftmp15*Ftmp56*Ftmp8*y*M[19] + 105*Ftmp15*Ftmp56*x*y*z*M[21] + 105*Ftmp15*Ftmp594*y*M[67] + 315*Ftmp15*Ftmp600*M[100] + 210*Ftmp15*Ftmp659*x*y*z*M[26] + 105*Ftmp15*Ftmp66*x*y*z*M[26] + 105*Ftmp15*Ftmp68*x*y*z*M[28] + 1890*Ftmp15*y*(-Ftmp185 + Ftmp708 + 33*Ftmp709 - 5)*M[76] + 210*Ftmp15*y*(-117*Ftmp306 + Ftmp592 + Ftmp659 + Ftmp773)*M[67] + 630*Ftmp15*y*(Ftmp451 + Ftmp467 + Ftmp659 + Ftmp744)*M[58] + 630*Ftmp15*y*(Ftmp453 + Ftmp468 + Ftmp748 + Ftmp749)*M[65] + 630*Ftmp15*y*(Ftmp458 + Ftmp659 + Ftmp746 + Ftmp747)*M[80] + 630*Ftmp15*y*(Ftmp476 + Ftmp749 + Ftmp750 + Ftmp751)*M[78] + 15*Ftmp150*Ftmp7*y*M[19] + 15*Ftmp154*Ftmp7*y*M[29] + 15*Ftmp158*Ftmp7*y*M[33] + 1890*Ftmp162*Ftmp27*Ftmp8*y*M[55] + 2835*Ftmp162*Ftmp27*x*y*z*M[57] - Ftmp163*Ftmp682 - Ftmp166*Ftmp682 - Ftmp171*Ftmp682 - Ftmp172*Ftmp684 - Ftmp172*Ftmp685 - Ftmp172*Ftmp719 - Ftmp172*Ftmp720 + 2835*Ftmp174*Ftmp27*x*y*z*M[75] - Ftmp175*M[81] + 45*Ftmp182*Ftmp7*M[35] - Ftmp183*Ftmp690 - Ftmp184*Ftmp205 - Ftmp184*Ftmp207 + 15*Ftmp186*Ftmp7*M[49] - Ftmp187*Ftmp692 - Ftmp189*Ftmp579*M[72] + 210*Ftmp19*Ftmp27*x*y*M[28] - Ftmp190*Ftmp204 - Ftmp190*Ftmp372 - Ftmp190*Ftmp376 - Ftmp190*Ftmp580 - Ftmp191*Ftmp690 - Ftmp192*Ftmp688 - Ftmp194*Ftmp688 - Ftmp197*Ftmp689 - Ftmp20*Ftmp652 + 945*Ftmp200*Ftmp27*x*z*M[87] + 945*Ftmp202*Ftmp27*x*z*M[109] - Ftmp21*Ftmp656 - Ftmp22*Ftmp30 - Ftmp22*(48*Ftmp352 + Ftmp41*Ftmp61 - Ftmp581*Ftmp771 - Ftmp779 + Ftmp780)*M[46] - Ftmp221*Ftmp35 - Ftmp223*Ftmp240 - Ftmp223*(Ftmp201*Ftmp8 + Ftmp412 + Ftmp734 + Ftmp8)*M[97] - Ftmp225*Ftmp706 - Ftmp23*Ftmp664*Ftmp665 - Ftmp233*Ftmp701 - Ftmp235*Ftmp701 - Ftmp238*Ftmp705 - Ftmp243*(Ftmp417 + Ftmp738 + Ftmp739 + Ftmp740)*M[88] - Ftmp244*(Ftmp704 + Ftmp707)*M[64] - Ftmp263*Ftmp711*M[104] - Ftmp265*Ftmp671 - Ftmp267*Ftmp291 - Ftmp267*Ftmp292 - Ftmp267*Ftmp293 - Ftmp267*Ftmp570 - Ftmp267*Ftmp575 - Ftmp267*Ftmp610 - Ftmp267*Ftmp627 - Ftmp269*Ftmp671 + 2835*Ftmp27*Ftmp324*x*y*z*M[62] + 8505*Ftmp27*Ftmp328*x*y*z*M[64] + 945*Ftmp27*Ftmp366*x*z*M[94] + 945*Ftmp27*Ftmp368*x*z*M[96] + 2835*Ftmp27*Ftmp397*x*y*z*M[71] + 945*Ftmp27*Ftmp400*x*z*M[105] + 405*Ftmp27*Ftmp614*x*y*z*M[73] + 315*Ftmp27*Ftmp618*x*z*M[107] + 210*Ftmp27*Ftmp653*y*z*M[21] + 1890*Ftmp27*x*y*z*(Ftmp378 + Ftmp725 + Ftmp727)*M[62] + 1260*Ftmp27*x*y*z*(-Ftmp1*Ftmp413*Ftmp771 + Ftmp134 - 34*Ftmp39 + Ftmp49*Ftmp772 + 9)*M[71] + 90*Ftmp27*x*y*z*(Ftmp285*Ftmp41 - 126*Ftmp42 + Ftmp478 - Ftmp52*Ftmp781 + Ftmp613 + 28)*M[73] + 210*Ftmp27*x*y*(Ftmp693 + Ftmp694)*M[39] + 630*Ftmp27*x*y*(56*Ftmp1*Ftmp12*Ftmp49 + 72*Ftmp11*Ftmp12*Ftmp7 - Ftmp138*Ftmp96 - 168*Ftmp352*Ftmp659 - Ftmp369 - Ftmp392 - Ftmp396*(-22*Ftmp42 + Ftmp727 + Ftmp745) - 14*Ftmp514*Ftmp686 - Ftmp785 + 4)*M[106] + 210*Ftmp27*y*z*(Ftmp693 + Ftmp697)*M[43] + 630*Ftmp27*y*(Ftmp12 + Ftmp730 + 99*Ftmp732 + Ftmp733)*M[60] + 630*Ftmp27*y*(Ftmp406 + Ftmp730 + 99*Ftmp731 + Ftmp8)*M[69] - Ftmp270*Ftmp671 - Ftmp272*Ftmp712 - Ftmp274*Ftmp712 - Ftmp276*Ftmp712 - Ftmp279*Ftmp664*M[15] - Ftmp29*M[10] - Ftmp3*M[7] + 15*Ftmp303*Ftmp7*y*M[22] + 15*Ftmp309*Ftmp7*y*M[24] + 15*Ftmp310*Ftmp7*y*M[31] - Ftmp318*Ftmp718 + 3780*Ftmp32*Ftmp653*Ftmp676*y*z*M[57] - Ftmp322*Ftmp682 - Ftmp326*Ftmp721 - Ftmp337*(Ftmp147 + Ftmp313 + Ftmp331)*M[37] + 15*Ftmp341*Ftmp7*M[40] - Ftmp342*Ftmp692 - Ftmp35*Ftmp97 + 15*Ftmp350*Ftmp7*M[42] - Ftmp351*Ftmp692 + 15*Ftmp354*Ftmp7*M[51] - Ftmp355*Ftmp692 - Ftmp359*Ftmp689 - Ftmp362*Ftmp689 - Ftmp364*Ftmp689 - Ftmp391*Ftmp701 - Ftmp395*Ftmp701 - Ftmp4*M[4] - Ftmp416*Ftmp701 - Ftmp48*M[1] - Ftmp490*Ftmp671 - Ftmp495*(Ftmp565 + Ftmp753 + Ftmp755)*M[93] - Ftmp495*(Ftmp604 + Ftmp776 + Ftmp777)*M[95] - Ftmp495*(Ftmp443 + Ftmp491 + Ftmp756 + Ftmp757)*M[86] - Ftmp5*M[0] - Ftmp500*Ftmp671 - Ftmp505*Ftmp766 - Ftmp507*Ftmp558 - Ftmp507*Ftmp566 - Ftmp516*Ftmp671 - Ftmp525*Ftmp766 - Ftmp529*Ftmp671 - Ftmp532*Ftmp712 - Ftmp535*Ftmp712 - Ftmp540*Ftmp712 - Ftmp544*Ftmp761 - Ftmp549*Ftmp712 - Ftmp551*Ftmp761 - Ftmp57*Ftmp660 - Ftmp6*y - Ftmp60*Ftmp660 - Ftmp601*Ftmp671 - Ftmp607*Ftmp712 - Ftmp619*Ftmp700 - Ftmp63*Ftmp660 - Ftmp635*Ftmp699*Ftmp782*(198*Ftmp352 + Ftmp367*Ftmp59 - 198*Ftmp42 + Ftmp617 - Ftmp633*Ftmp771 + 36) - Ftmp64*Ftmp661 - Ftmp64*Ftmp662 - Ftmp652*Ftmp654 - Ftmp652*(-Ftmp12*Ftmp179 + Ftmp683)*M[24] - Ftmp653*Ftmp658*M[38] - Ftmp657*Ftmp95*M[47] - Ftmp657*x*(Ftmp694 + Ftmp704)*M[63] - Ftmp657*(Ftmp412 + Ftmp696 + Ftmp741 + Ftmp743)*M[102] - Ftmp657*(Ftmp50 + Ftmp736 + Ftmp738 + Ftmp743)*M[91] - Ftmp658*(Ftmp697 + Ftmp704)*M[68] - Ftmp659*Ftmp663*y*M[12] - Ftmp666*Ftmp82 - Ftmp666*Ftmp83 - Ftmp666*Ftmp85 - Ftmp667*Ftmp76 - Ftmp667*Ftmp79 - Ftmp667*Ftmp80 - Ftmp668*Ftmp669 - Ftmp668*Ftmp698 - Ftmp670*(Ftmp260 + Ftmp486 + Ftmp574 + Ftmp765)*M[98] - Ftmp670*(-209*Ftmp306 + Ftmp596 + Ftmp659 + Ftmp775)*M[100] - Ftmp670*(Ftmp512 + Ftmp562 + Ftmp767 + Ftmp768)*M[89] - Ftmp670*(Ftmp527 + Ftmp578 + Ftmp765 + Ftmp770)*M[113] - Ftmp670*(Ftmp547 + Ftmp569 + Ftmp768 + Ftmp769)*M[115] - Ftmp679*Ftmp680*Ftmp681 - Ftmp69*M[32] + 6*Ftmp7*Ftmp8*y*M[3] + 15*Ftmp7*x*y*z*M[5] + 15*Ftmp7*x*z*M[13] + 60*Ftmp7*y*(7*Ftmp153 + Ftmp335)*M[29] + 30*Ftmp7*y*(Ftmp297 + Ftmp714 + Ftmp716)*M[22] + 30*Ftmp7*y*(Ftmp305 + Ftmp716 + Ftmp717)*M[31] - 20790*Ftmp700*(Ftmp215 + Ftmp389 + Ftmp727)*M[94] - 3780*Ftmp700*(Ftmp117*Ftmp772 - 22*Ftmp39*Ftmp771 - 166*Ftmp39 + Ftmp399 + 55)*M[105] - Ftmp711*Ftmp713*M[112] - Ftmp713*(Ftmp485 + Ftmp502 + Ftmp659 + Ftmp752)*M[90] - Ftmp723*(Ftmp146 + Ftmp722)*M[41] - Ftmp723*(Ftmp360 + Ftmp623)*M[52] - Ftmp726*(Ftmp722 + Ftmp725)*M[40] - Ftmp726*(Ftmp375 + Ftmp622 + Ftmp715)*M[51] - Ftmp737*(Ftmp12*Ftmp199 + Ftmp12 + Ftmp734 + Ftmp736)*M[92] - Ftmp742*(Ftmp422 + Ftmp739 + Ftmp741 + Ftmp87)*M[103] - Ftmp760*(Ftmp608 + Ftmp776 + Ftmp778)*M[101] - Ftmp760*(Ftmp519 + Ftmp753 + Ftmp763 + Ftmp764)*M[114] - Ftmp760*(Ftmp537 + Ftmp547 + Ftmp756 + Ftmp762)*M[116] - Ftmp760*(Ftmp486 + Ftmp755 + Ftmp758 + Ftmp759 - 6)*M[99] - Ftmp77*Ftmp86 - Ftmp77*Ftmp89 - Ftmp77*Ftmp90 - Ftmp783*Ftmp784 - 90*Ftmp783*(1160*Ftmp15*Ftmp480 - 128*Ftmp157 - Ftmp170*Ftmp781 + 140*Ftmp287*Ftmp631 - 32*Ftmp352*Ftmp643 + Ftmp352*Ftmp648 - Ftmp42*Ftmp636 - Ftmp42*Ftmp637 - Ftmp42*Ftmp651 + 96*Ftmp42 - Ftmp554*Ftmp649 + Ftmp786 - 8)*M[108] + M[1]);
#pragma omp atomic
F[2] += Ftmp0*(3*Ftmp1*Ftmp41*M[16] + 3*Ftmp1*Ftmp46*z*M[3] + 3*Ftmp1*Ftmp49*z*M[6] + 3*Ftmp1*Ftmp52*z*M[8] + 3*Ftmp1*Ftmp75*M[11] + 6*Ftmp1*Ftmp789*z*M[8] + 3*Ftmp1*Ftmp84*M[18] - Ftmp100*Ftmp239*Ftmp705 - Ftmp100*Ftmp677*Ftmp98 - Ftmp100*Ftmp797*M[53] - Ftmp101*Ftmp220*(-297*Ftmp157 + 189*Ftmp42 + Ftmp813 - 35) - Ftmp102*Ftmp789*Ftmp792*M[47] - Ftmp105*Ftmp106*Ftmp119 - Ftmp105*Ftmp27*Ftmp784 - Ftmp105*Ftmp797*M[48] - Ftmp108*Ftmp12 - Ftmp109*Ftmp794 + 3780*Ftmp11*Ftmp12*Ftmp122*Ftmp32*x*M[71] + 210*Ftmp11*Ftmp12*Ftmp27*x*M[26] + 60*Ftmp11*Ftmp15*Ftmp59*z*M[29] + 1890*Ftmp11*Ftmp165*Ftmp27*z*M[76] + 6*Ftmp11*Ftmp7*z*M[6] - Ftmp112*Ftmp794 - Ftmp113*Ftmp794 + 3780*Ftmp117*Ftmp12*Ftmp13*Ftmp32*M[77] + 3780*Ftmp117*Ftmp13*Ftmp32*x*z*M[70] + 20790*Ftmp12*Ftmp123*x*y*(-26*Ftmp809 + Ftmp811)*M[94] + 210*Ftmp12*Ftmp13*Ftmp27*M[30] + 105*Ftmp12*Ftmp15*Ftmp56*x*M[21] + 105*Ftmp12*Ftmp15*Ftmp59*y*M[30] + 105*Ftmp12*Ftmp15*Ftmp66*x*M[26] + 105*Ftmp12*Ftmp15*Ftmp68*x*M[28] + 105*Ftmp12*Ftmp15*Ftmp68*y*M[32] + 210*Ftmp12*Ftmp15*Ftmp789*x*M[28] + 210*Ftmp12*Ftmp15*Ftmp789*y*M[32] + 105*Ftmp12*Ftmp15*Ftmp91*y*M[23] + 2835*Ftmp12*Ftmp162*Ftmp27*x*M[57] + 2835*Ftmp12*Ftmp165*Ftmp27*y*M[77] + 2835*Ftmp12*Ftmp174*Ftmp27*x*M[75] + 2835*Ftmp12*Ftmp174*Ftmp27*y*M[81] - Ftmp12*Ftmp190*M[13] + 2835*Ftmp12*Ftmp212*Ftmp27*y*M[59] - Ftmp12*Ftmp227*Ftmp706 + 2835*Ftmp12*Ftmp27*Ftmp324*x*M[62] + 8505*Ftmp12*Ftmp27*Ftmp328*x*M[64] + 2835*Ftmp12*Ftmp27*Ftmp383*y*M[66] + 2835*Ftmp12*Ftmp27*Ftmp385*y*M[68] + 8505*Ftmp12*Ftmp27*Ftmp386*y*M[79] + 2835*Ftmp12*Ftmp27*Ftmp397*x*M[71] + 405*Ftmp12*Ftmp27*Ftmp614*x*M[73] + 210*Ftmp12*Ftmp27*Ftmp653*M[21] + 210*Ftmp12*Ftmp27*Ftmp8*y*M[23] + 3780*Ftmp12*Ftmp27*Ftmp804*x*M[75] + 3780*Ftmp12*Ftmp27*Ftmp804*y*M[81] + 1890*Ftmp12*Ftmp27*x*(Ftmp209 + Ftmp819)*M[64] + 90*Ftmp12*Ftmp27*x*(Ftmp52*Ftmp708 + Ftmp613 - Ftmp630*Ftmp838 - 126*Ftmp840)*M[73] + 1890*Ftmp12*Ftmp27*y*(Ftmp136 + Ftmp728 + Ftmp820)*M[79] + 1890*Ftmp12*Ftmp27*y*(Ftmp181 + Ftmp384 + Ftmp815)*M[68] + 210*Ftmp12*Ftmp27*(Ftmp796 + Ftmp800)*M[41] + 3780*Ftmp12*Ftmp32*Ftmp653*Ftmp676*M[57] + 3780*Ftmp12*Ftmp32*Ftmp678*Ftmp8*y*M[59] - Ftmp12*Ftmp32*Ftmp782*y*(Ftmp121*Ftmp277 + Ftmp617 - Ftmp634*Ftmp838 - 198*Ftmp840)*M[107] - Ftmp12*Ftmp655 - Ftmp12*Ftmp656*M[16] - Ftmp12*Ftmp674 - Ftmp12*Ftmp691*(Ftmp142 + 9*Ftmp157 + 5)*M[54] + 15*Ftmp12*Ftmp7*x*M[5] + 15*Ftmp12*Ftmp7*y*M[7] - Ftmp12*Ftmp702 - Ftmp125*Ftmp807*(13*Ftmp157 + Ftmp168 + 5) - Ftmp127*Ftmp14*Ftmp795 + 210*Ftmp13*Ftmp27*x*z*M[25] + 45*Ftmp138*Ftmp7*M[50] - Ftmp14*Ftmp788*M[45] + 315*Ftmp15*Ftmp251*M[112] + 315*Ftmp15*Ftmp256*z*M[55] + 315*Ftmp15*Ftmp258*z*M[76] + 315*Ftmp15*Ftmp259*z*M[82] + 315*Ftmp15*Ftmp264*M[85] + 315*Ftmp15*Ftmp275*M[118] + 315*Ftmp15*Ftmp46*x*y*M[38] + 315*Ftmp15*Ftmp461*z*M[65] + 315*Ftmp15*Ftmp464*z*M[69] + 315*Ftmp15*Ftmp473*z*M[58] + 315*Ftmp15*Ftmp477*z*M[60] + 315*Ftmp15*Ftmp481*z*M[80] + 315*Ftmp15*Ftmp483*z*M[78] + 315*Ftmp15*Ftmp49*x*y*M[45] + 315*Ftmp15*Ftmp52*x*y*M[47] + 315*Ftmp15*Ftmp531*M[99] + 315*Ftmp15*Ftmp534*M[103] + 315*Ftmp15*Ftmp539*M[90] + 945*Ftmp15*Ftmp543*M[92] + 315*Ftmp15*Ftmp548*M[116] + 945*Ftmp15*Ftmp550*M[114] + 60*Ftmp15*Ftmp56*Ftmp8*z*M[19] + 105*Ftmp15*Ftmp56*x*y*z*M[20] + 105*Ftmp15*Ftmp59*x*y*z*M[25] + 105*Ftmp15*Ftmp594*z*M[67] + 315*Ftmp15*Ftmp606*M[101] + 105*Ftmp15*Ftmp62*x*y*z*M[27] + 210*Ftmp15*Ftmp789*x*y*z*M[27] + 1890*Ftmp15*z*(-Ftmp195 + 35*Ftmp42 + 33*Ftmp812 - 5)*M[82] + 630*Ftmp15*z*(Ftmp455 + Ftmp746 + Ftmp751 + Ftmp789)*M[78] + 630*Ftmp15*z*(Ftmp456 + Ftmp475 + Ftmp744 + Ftmp789)*M[60] + 630*Ftmp15*z*(Ftmp457 + Ftmp474 + Ftmp748 + Ftmp828)*M[69] + 630*Ftmp15*z*(Ftmp465 + Ftmp747 + Ftmp750 + Ftmp828)*M[80] + 210*Ftmp15*z*(-117*Ftmp296 + Ftmp591 + Ftmp593 + Ftmp773 + Ftmp789)*M[67] + 15*Ftmp150*Ftmp7*z*M[19] + 15*Ftmp154*Ftmp7*z*M[29] + 15*Ftmp158*Ftmp7*z*M[33] - Ftmp160*Ftmp174*M[81] - Ftmp160*Ftmp684 - Ftmp160*Ftmp685 - Ftmp160*Ftmp719 - Ftmp160*Ftmp720 + 1890*Ftmp162*Ftmp27*Ftmp8*z*M[55] + 2835*Ftmp162*Ftmp27*x*y*z*M[56] + 2835*Ftmp165*Ftmp27*x*y*z*M[70] - Ftmp17*Ftmp36 + 2835*Ftmp170*Ftmp27*x*y*z*M[74] - Ftmp173*Ftmp682 - Ftmp174*Ftmp682*M[75] - Ftmp18*z*(40*Ftmp352 + Ftmp40*Ftmp62 - Ftmp454 - Ftmp582*Ftmp838 + Ftmp839)*M[46] + 45*Ftmp182*Ftmp7*M[36] - Ftmp183*Ftmp688 - Ftmp187*Ftmp689 - Ftmp191*Ftmp688 - Ftmp192*Ftmp798 - Ftmp193*Ftmp205 - Ftmp193*Ftmp207 - Ftmp194*Ftmp798 + 15*Ftmp196*Ftmp7*M[54] - Ftmp197*Ftmp799 - Ftmp198*Ftmp204 - Ftmp198*Ftmp372 - Ftmp198*Ftmp376 - Ftmp198*Ftmp580 - Ftmp2*y*M[7] + 945*Ftmp200*Ftmp27*x*y*M[87] + 945*Ftmp202*Ftmp27*x*y*M[109] - Ftmp219*Ftmp32*Ftmp36 - 10395*Ftmp224*Ftmp806*M[109] - Ftmp229*(Ftmp217*Ftmp8 + Ftmp415 + Ftmp8 + Ftmp823)*M[93] - Ftmp230*Ftmp236 - Ftmp233*Ftmp808 - Ftmp24*z*M[12] - Ftmp244*(Ftmp810 + Ftmp811)*M[61] - Ftmp265*Ftmp712 - Ftmp269*Ftmp712 + 8505*Ftmp27*Ftmp317*x*y*z*M[61] + 2835*Ftmp27*Ftmp321*x*y*z*M[63] + 945*Ftmp27*Ftmp366*x*y*M[94] + 945*Ftmp27*Ftmp368*x*y*M[96] + 945*Ftmp27*Ftmp400*x*y*M[105] + 945*Ftmp27*Ftmp579*x*y*z*M[72] + 315*Ftmp27*Ftmp618*x*y*M[107] + 210*Ftmp27*Ftmp653*y*z*M[20] + 3780*Ftmp27*Ftmp802*x*y*z*M[74] + 1890*Ftmp27*x*y*z*(Ftmp725 + Ftmp819 + 3)*M[63] + 210*Ftmp27*x*y*z*(56*Ftmp352 + Ftmp514*Ftmp58 - 9*Ftmp59*Ftmp838 - Ftmp759 + Ftmp839)*M[72] + 210*Ftmp27*x*z*(Ftmp800 + Ftmp803)*M[37] + 630*Ftmp27*x*z*(56*Ftmp1*Ftmp11*Ftmp49 - 9*Ftmp138*Ftmp789 - 168*Ftmp15*Ftmp482 - Ftmp396*(Ftmp278 + Ftmp620 + Ftmp745 + 4) - Ftmp49*Ftmp616 - Ftmp785)*M[106] + 90*Ftmp27*x*z*(280*Ftmp1*Ftmp11*Ftmp52 + 160*Ftmp1*Ftmp11*Ftmp611 - Ftmp170*Ftmp708 - 1160*Ftmp352*Ftmp789 - Ftmp644*(-Ftmp114 + 8*Ftmp157 + 3) - Ftmp647*Ftmp841 - Ftmp648*Ftmp841 - Ftmp650*Ftmp802 - Ftmp786)*M[108] + 210*Ftmp27*y*z*(Ftmp800 + Ftmp805)*M[40] + 630*Ftmp27*z*(Ftmp11 + Ftmp733 + Ftmp821 + 99*Ftmp822)*M[58] + 630*Ftmp27*z*(Ftmp399*Ftmp8 + Ftmp408 + Ftmp8 + Ftmp821)*M[65] - Ftmp270*Ftmp712 - Ftmp271*Ftmp814*M[110] - Ftmp272*Ftmp793 - Ftmp273*Ftmp291 - Ftmp273*Ftmp292 - Ftmp273*Ftmp293 - Ftmp273*Ftmp570 - Ftmp273*Ftmp575 - Ftmp273*Ftmp610 - Ftmp273*Ftmp627 - Ftmp274*Ftmp793 - Ftmp276*Ftmp793 - 30*Ftmp287*Ftmp789*M[18] - Ftmp3*x*M[0] - Ftmp3*y*M[1] + 15*Ftmp303*Ftmp7*z*M[22] + 15*Ftmp309*Ftmp7*z*M[24] - Ftmp31*Ftmp73*z + 15*Ftmp310*Ftmp7*z*M[31] - Ftmp312*Ftmp721 + 3780*Ftmp32*Ftmp653*Ftmp676*y*z*M[56] - Ftmp325*Ftmp682 - Ftmp329*Ftmp718 - Ftmp342*Ftmp689 - Ftmp351*Ftmp689 - Ftmp355*Ftmp689 - Ftmp357*(Ftmp147 + Ftmp327 + Ftmp346)*M[39] + 15*Ftmp358*Ftmp7*M[41] - Ftmp359*Ftmp799 - Ftmp36*Ftmp94 + 15*Ftmp361*Ftmp7*M[43] - Ftmp362*Ftmp799 + 15*Ftmp363*Ftmp7*M[52] - Ftmp364*Ftmp799 - Ftmp391*Ftmp808 - Ftmp395*Ftmp808 - Ftmp398*Ftmp682 - Ftmp4*M[5] - Ftmp416*Ftmp808 - Ftmp423*(Ftmp421 + Ftmp740 + Ftmp824 + Ftmp825)*M[86] - Ftmp490*Ftmp712 - Ftmp500*Ftmp712 - Ftmp505*Ftmp761 - Ftmp51*M[2] - Ftmp516*Ftmp712 - Ftmp525*Ftmp761 - Ftmp529*Ftmp712 - Ftmp532*Ftmp793 - Ftmp535*Ftmp793 - Ftmp538*(Ftmp496 + Ftmp829 + Ftmp830)*M[97] - Ftmp538*(Ftmp609 + Ftmp777 + Ftmp837)*M[95] - Ftmp538*(Ftmp447 + Ftmp536 + Ftmp757 + Ftmp831)*M[88] - Ftmp54*Ftmp661 - Ftmp54*Ftmp662 - Ftmp54*Ftmp68*M[32] - Ftmp540*Ftmp793 - Ftmp544*Ftmp835 - Ftmp545*Ftmp558 - Ftmp545*Ftmp566 - Ftmp549*Ftmp793 - Ftmp551*Ftmp835 - Ftmp601*Ftmp712 - Ftmp607*Ftmp793 - Ftmp615*x*M[73] - Ftmp619*Ftmp807 - Ftmp65*Ftmp660 - Ftmp654*Ftmp787 - Ftmp660*Ftmp67 - Ftmp660*Ftmp68*M[28] - Ftmp663*Ftmp790*M[14] - 30*Ftmp665*Ftmp790*M[17] - Ftmp666*Ftmp76 - Ftmp666*Ftmp79 - Ftmp666*Ftmp80 - Ftmp668*Ftmp681*Ftmp795 - Ftmp668*Ftmp788*M[38] - Ftmp669*Ftmp680 - Ftmp680*Ftmp698 + 6*Ftmp7*Ftmp8*z*M[3] + 15*Ftmp7*x*y*z*M[4] + 15*Ftmp7*x*y*M[13] + 60*Ftmp7*z*(7*Ftmp157 + Ftmp360)*M[33] + 30*Ftmp7*z*(Ftmp295 + Ftmp717 + Ftmp816)*M[31] + 30*Ftmp7*z*(Ftmp307 + Ftmp714 + Ftmp816)*M[24] - Ftmp713*Ftmp814*M[117] - Ftmp713*(Ftmp484 + Ftmp542 + Ftmp752 + Ftmp789)*M[91] - Ftmp723*(Ftmp146 + Ftmp817)*M[42] - Ftmp723*(Ftmp327 + Ftmp335 + Ftmp622)*M[51] - Ftmp737*(Ftmp11*Ftmp199 + Ftmp11 + Ftmp736 + Ftmp823)*M[89] - Ftmp742*(Ftmp420 + Ftmp825 + Ftmp826 + Ftmp87)*M[98] - Ftmp760*(Ftmp605 + Ftmp778 + Ftmp837)*M[100] - Ftmp760*(Ftmp493 + Ftmp527 + Ftmp762 + Ftmp831)*M[113] - Ftmp760*(Ftmp508 + Ftmp547 + Ftmp763 + Ftmp830)*M[115] - Ftmp760*(Ftmp779 + Ftmp829 + Ftmp833 - 6)*M[102] - Ftmp787*(-Ftmp11*Ftmp179 + Ftmp796)*M[22] - Ftmp788*x*(Ftmp803 + Ftmp810)*M[62] - Ftmp788*y*(Ftmp805 + Ftmp810)*M[66] - Ftmp788*(Ftmp415 + Ftmp696 + Ftmp826 + Ftmp827)*M[99] - Ftmp788*(Ftmp47 + Ftmp736 + Ftmp824 + Ftmp827)*M[90] - Ftmp791*Ftmp82 - Ftmp791*Ftmp83 - Ftmp791*Ftmp85 - Ftmp792*(Ftmp260 + Ftmp576 + Ftmp833 + Ftmp834)*M[103] - Ftmp792*(Ftmp546 + Ftmp567 + Ftmp767 + Ftmp836)*M[92] - Ftmp792*(Ftmp560 + Ftmp764 + Ftmp769 + Ftmp836)*M[114] - Ftmp792*(-209*Ftmp296 + Ftmp586 + Ftmp598 + Ftmp775 + Ftmp789)*M[101] - Ftmp792*(Ftmp547 + Ftmp571 + Ftmp770 + Ftmp832 + Ftmp834)*M[116] - 20790*Ftmp807*(Ftmp215 + Ftmp394 + Ftmp820)*M[96] - Ftmp81*Ftmp86 - Ftmp81*Ftmp89 - Ftmp81*Ftmp90 - Ftmp818*(Ftmp725 + Ftmp817)*M[43] - Ftmp818*(Ftmp370 + Ftmp622 + Ftmp816)*M[52] + M[2]);

}

void S2Mc_7(double x, double y, double z, double * S, double * M) {
double Mtmp0 = x*S[0];
double Mtmp1 = z*S[2];
double Mtmp2 = -Mtmp1;
double Mtmp3 = x*S[1];
double Mtmp4 = y*S[0];
double Mtmp5 = x*S[2];
double Mtmp6 = z*S[0];
double Mtmp7 = y*S[1];
double Mtmp8 = y*S[2];
double Mtmp9 = z*S[1];
double Mtmp10 = Mtmp1*x;
double Mtmp11 = pow(x, 2);
double Mtmp12 = pow(z, 2);
double Mtmp13 = (1.0/2.0)*S[0];
double Mtmp14 = Mtmp12*Mtmp13;
double Mtmp15 = Mtmp0*y;
double Mtmp16 = Mtmp1*y;
double Mtmp17 = (1.0/2.0)*S[1];
double Mtmp18 = Mtmp12*Mtmp17;
double Mtmp19 = Mtmp0*z;
double Mtmp20 = (1.0/2.0)*S[2];
double Mtmp21 = -Mtmp12*Mtmp20;
double Mtmp22 = Mtmp3*y;
double Mtmp23 = pow(y, 2);
double Mtmp24 = Mtmp5*y;
double Mtmp25 = Mtmp3*z;
double Mtmp26 = Mtmp4*z;
double Mtmp27 = Mtmp7*z;
double Mtmp28 = pow(x, 3);
double Mtmp29 = Mtmp28*S[0];
double Mtmp30 = pow(z, 3);
double Mtmp31 = Mtmp30*S[2];
double Mtmp32 = 3*Mtmp12;
double Mtmp33 = Mtmp0*Mtmp32;
double Mtmp34 = 3*Mtmp11;
double Mtmp35 = Mtmp1*Mtmp34;
double Mtmp36 = (1.0/2.0)*Mtmp12;
double Mtmp37 = Mtmp10*y + (1.0/2.0)*Mtmp12*Mtmp4 + Mtmp3*Mtmp36;
double Mtmp38 = Mtmp28*S[2];
double Mtmp39 = Mtmp30*S[0];
double Mtmp40 = Mtmp32*Mtmp7;
double Mtmp41 = 3*Mtmp23;
double Mtmp42 = Mtmp1*Mtmp41;
double Mtmp43 = Mtmp30*S[1];
double Mtmp44 = (1.0/2.0)*Mtmp11;
double Mtmp45 = pow(y, 3);
double Mtmp46 = (1.0/2.0)*Mtmp23;
double Mtmp47 = Mtmp45*S[1];
double Mtmp48 = Mtmp45*S[2];
double Mtmp49 = pow(x, 4);
double Mtmp50 = Mtmp49*S[0];
double Mtmp51 = Mtmp11*S[0];
double Mtmp52 = 6*Mtmp12;
double Mtmp53 = 4*Mtmp28;
double Mtmp54 = pow(z, 4);
double Mtmp55 = Mtmp54*S[0];
double Mtmp56 = 4*Mtmp30;
double Mtmp57 = Mtmp5*Mtmp56 + Mtmp55;
double Mtmp58 = Mtmp49*S[1];
double Mtmp59 = Mtmp11*Mtmp52;
double Mtmp60 = 12*Mtmp12;
double Mtmp61 = Mtmp11*Mtmp16;
double Mtmp62 = Mtmp54*S[1];
double Mtmp63 = Mtmp56*Mtmp8 + Mtmp62;
double Mtmp64 = Mtmp49*S[2];
double Mtmp65 = Mtmp54*S[2];
double Mtmp66 = 2*Mtmp28;
double Mtmp67 = 6*Mtmp23;
double Mtmp68 = Mtmp23*S[0];
double Mtmp69 = Mtmp28*Mtmp8;
double Mtmp70 = Mtmp28*Mtmp9;
double Mtmp71 = Mtmp3*Mtmp30;
double Mtmp72 = Mtmp30*Mtmp4;
double Mtmp73 = -Mtmp24*Mtmp32 - Mtmp71 - Mtmp72;
double Mtmp74 = 2*Mtmp45;
double Mtmp75 = Mtmp32*S[1];
double Mtmp76 = 2*Mtmp30;
double Mtmp77 = Mtmp32*S[2];
double Mtmp78 = pow(y, 4);
double Mtmp79 = Mtmp78*S[0];
double Mtmp80 = 4*Mtmp45;
double Mtmp81 = Mtmp45*Mtmp5;
double Mtmp82 = Mtmp45*Mtmp6;
double Mtmp83 = Mtmp78*S[1];
double Mtmp84 = Mtmp23*Mtmp52;
double Mtmp85 = Mtmp78*S[2];
double Mtmp86 = pow(x, 5);
double Mtmp87 = pow(z, 5);
double Mtmp88 = Mtmp87*S[2];
double Mtmp89 = -Mtmp88;
double Mtmp90 = 10*Mtmp12;
double Mtmp91 = 5*Mtmp49;
double Mtmp92 = -Mtmp1*Mtmp91 - Mtmp29*Mtmp90;
double Mtmp93 = 5*Mtmp54;
double Mtmp94 = 10*Mtmp11;
double Mtmp95 = Mtmp0*Mtmp93 + Mtmp31*Mtmp94;
double Mtmp96 = Mtmp12*Mtmp28*S[1];
double Mtmp97 = 30*Mtmp11;
double Mtmp98 = Mtmp12*Mtmp97;
double Mtmp99 = 20*Mtmp28;
double Mtmp100 = 20*Mtmp30;
double Mtmp101 = Mtmp100*Mtmp24 + 5*Mtmp3*Mtmp54 + 5*Mtmp4*Mtmp54;
double Mtmp102 = Mtmp5*Mtmp93 + Mtmp87*S[0];
double Mtmp103 = 10*Mtmp54;
double Mtmp104 = 10*Mtmp23;
double Mtmp105 = 20*Mtmp31;
double Mtmp106 = Mtmp104*Mtmp31 + Mtmp7*Mtmp93;
double Mtmp107 = 30*Mtmp23;
double Mtmp108 = Mtmp107*Mtmp12;
double Mtmp109 = Mtmp1*Mtmp23;
double Mtmp110 = -Mtmp0*Mtmp108 - Mtmp109*Mtmp97 - Mtmp7*Mtmp98 - 3*Mtmp88;
double Mtmp111 = Mtmp8*Mtmp93 + Mtmp87*S[1];
double Mtmp112 = Mtmp12*Mtmp45*S[0];
double Mtmp113 = 5*Mtmp39;
double Mtmp114 = 5*Mtmp38;
double Mtmp115 = 15*Mtmp12;
double Mtmp116 = Mtmp23*Mtmp5;
double Mtmp117 = 10*Mtmp30;
double Mtmp118 = 10*Mtmp28;
double Mtmp119 = 15*Mtmp11;
double Mtmp120 = Mtmp119*Mtmp23;
double Mtmp121 = 5*Mtmp78;
double Mtmp122 = -Mtmp1*Mtmp121 - Mtmp47*Mtmp90;
double Mtmp123 = 5*Mtmp43;
double Mtmp124 = 5*Mtmp48;
double Mtmp125 = Mtmp11*Mtmp8;
double Mtmp126 = 10*Mtmp45;
double Mtmp127 = pow(y, 5);
double Mtmp128 = 20*Mtmp45;
double Mtmp129 = pow(z, 6);
double Mtmp130 = Mtmp129*S[0];
double Mtmp131 = pow(x, 6);
double Mtmp132 = 6*Mtmp87;
double Mtmp133 = Mtmp132*Mtmp5;
double Mtmp134 = 6*Mtmp86;
double Mtmp135 = Mtmp1*Mtmp134;
double Mtmp136 = Mtmp115*Mtmp50;
double Mtmp137 = -15*Mtmp11*Mtmp54*S[0] - 20*Mtmp28*Mtmp30*S[2];
double Mtmp138 = Mtmp129*S[1];
double Mtmp139 = Mtmp115*Mtmp58;
double Mtmp140 = Mtmp132*Mtmp8;
double Mtmp141 = 60*Mtmp28;
double Mtmp142 = Mtmp12*Mtmp141;
double Mtmp143 = Mtmp142*Mtmp4;
double Mtmp144 = 30*Mtmp49;
double Mtmp145 = Mtmp144*Mtmp16;
double Mtmp146 = 30*Mtmp54;
double Mtmp147 = 60*Mtmp30;
double Mtmp148 = Mtmp119*Mtmp62 + Mtmp125*Mtmp147 + Mtmp146*Mtmp15;
double Mtmp149 = Mtmp129*S[2];
double Mtmp150 = -Mtmp149;
double Mtmp151 = -Mtmp115*Mtmp64 - 20*Mtmp29*Mtmp30;
double Mtmp152 = Mtmp0*Mtmp132 + Mtmp119*Mtmp65;
double Mtmp153 = 3*Mtmp130;
double Mtmp154 = 18*Mtmp87;
double Mtmp155 = Mtmp154*Mtmp5;
double Mtmp156 = 15*Mtmp23;
double Mtmp157 = 40*Mtmp31;
double Mtmp158 = 90*Mtmp12*Mtmp23;
double Mtmp159 = Mtmp158*Mtmp51;
double Mtmp160 = Mtmp142*Mtmp7;
double Mtmp161 = Mtmp109*Mtmp141;
double Mtmp162 = Mtmp116*Mtmp147 + Mtmp146*Mtmp22 + Mtmp156*Mtmp55;
double Mtmp163 = 3*Mtmp86;
double Mtmp164 = 30*Mtmp12;
double Mtmp165 = 3*Mtmp87;
double Mtmp166 = Mtmp165*Mtmp3 + Mtmp165*Mtmp4 + 15*Mtmp24*Mtmp54;
double Mtmp167 = 3*Mtmp138;
double Mtmp168 = Mtmp154*Mtmp8;
double Mtmp169 = 60*Mtmp45;
double Mtmp170 = Mtmp12*Mtmp169;
double Mtmp171 = Mtmp0*Mtmp170;
double Mtmp172 = Mtmp11*Mtmp169;
double Mtmp173 = Mtmp1*Mtmp172;
double Mtmp174 = Mtmp11*Mtmp158;
double Mtmp175 = Mtmp174*S[1];
double Mtmp176 = -15*Mtmp23*Mtmp54*S[1] - 20*Mtmp30*Mtmp45*S[2];
double Mtmp177 = 12*Mtmp87;
double Mtmp178 = Mtmp132*Mtmp7 + Mtmp156*Mtmp65;
double Mtmp179 = -Mtmp0*Mtmp147*Mtmp23 - Mtmp11*Mtmp147*Mtmp7 - 3*Mtmp149 - Mtmp174*S[2];
double Mtmp180 = Mtmp115*Mtmp79;
double Mtmp181 = 30*Mtmp78;
double Mtmp182 = Mtmp10*Mtmp181;
double Mtmp183 = Mtmp170*Mtmp3;
double Mtmp184 = Mtmp115*Mtmp83;
double Mtmp185 = 6*Mtmp127;
double Mtmp186 = Mtmp1*Mtmp185;
double Mtmp187 = -Mtmp115*Mtmp85 - 20*Mtmp30*Mtmp47;
double Mtmp188 = pow(y, 6);
double Mtmp189 = 3*Mtmp127;
#pragma omp atomic
M[0] += S[0];
#pragma omp atomic
M[1] += S[1];
#pragma omp atomic
M[2] += S[2];
#pragma omp atomic
M[3] += Mtmp0 + Mtmp2;
#pragma omp atomic
M[4] += Mtmp3 + Mtmp4;
#pragma omp atomic
M[5] += Mtmp5 + Mtmp6;
#pragma omp atomic
M[6] += Mtmp2 + Mtmp7;
#pragma omp atomic
M[7] += Mtmp8 + Mtmp9;
#pragma omp atomic
M[8] += -Mtmp10 + (1.0/2.0)*Mtmp11*S[0] - Mtmp14;
#pragma omp atomic
M[9] += Mtmp11*Mtmp17 + Mtmp15 - Mtmp16 - Mtmp18;
#pragma omp atomic
M[10] += Mtmp11*Mtmp20 + Mtmp19 + Mtmp21;
#pragma omp atomic
M[11] += -Mtmp10 + Mtmp13*Mtmp23 - Mtmp14 + Mtmp22;
#pragma omp atomic
M[12] += Mtmp24 + Mtmp25 + Mtmp26;
#pragma omp atomic
M[13] += -Mtmp16 - Mtmp18 + (1.0/2.0)*Mtmp23*S[1];
#pragma omp atomic
M[14] += Mtmp20*Mtmp23 + Mtmp21 + Mtmp27;
#pragma omp atomic
M[15] += (1.0/6.0)*Mtmp29 + (1.0/6.0)*Mtmp31 - 1.0/6.0*Mtmp33 - 1.0/6.0*Mtmp35;
#pragma omp atomic
M[16] += (1.0/2.0)*Mtmp11*y*S[0] + (1.0/6.0)*Mtmp28*S[1] - Mtmp37;
#pragma omp atomic
M[17] += -1.0/6.0*Mtmp32*Mtmp5 + (1.0/6.0)*Mtmp34*Mtmp6 + (1.0/6.0)*Mtmp38 - 1.0/6.0*Mtmp39;
#pragma omp atomic
M[18] += (1.0/2.0)*Mtmp11*y*S[1] + (1.0/2.0)*Mtmp23*x*S[0] + (1.0/3.0)*Mtmp30*S[2] - 1.0/6.0*Mtmp33 - 1.0/6.0*Mtmp35 - 1.0/6.0*Mtmp40 - 1.0/6.0*Mtmp42;
#pragma omp atomic
M[19] += Mtmp15*z - Mtmp36*Mtmp8 - 1.0/6.0*Mtmp43 + Mtmp44*Mtmp8 + Mtmp44*Mtmp9;
#pragma omp atomic
M[20] += (1.0/2.0)*Mtmp23*x*S[1] - Mtmp37 + (1.0/6.0)*Mtmp45*S[0];
#pragma omp atomic
M[21] += Mtmp22*z - Mtmp36*Mtmp5 - 1.0/6.0*Mtmp39 + Mtmp46*Mtmp5 + Mtmp46*Mtmp6;
#pragma omp atomic
M[22] += (1.0/6.0)*Mtmp31 - 1.0/6.0*Mtmp40 - 1.0/6.0*Mtmp42 + (1.0/6.0)*Mtmp47;
#pragma omp atomic
M[23] += -1.0/6.0*Mtmp32*Mtmp8 + (1.0/6.0)*Mtmp41*Mtmp9 - 1.0/6.0*Mtmp43 + (1.0/6.0)*Mtmp48;
#pragma omp atomic
M[24] += -1.0/24.0*Mtmp1*Mtmp53 + (1.0/24.0)*Mtmp50 - 1.0/24.0*Mtmp51*Mtmp52 + (1.0/24.0)*Mtmp57;
#pragma omp atomic
M[25] += -1.0/24.0*Mtmp15*Mtmp60 + (1.0/24.0)*Mtmp4*Mtmp53 + (1.0/24.0)*Mtmp58 - 1.0/24.0*Mtmp59*S[1] - 1.0/2.0*Mtmp61 + (1.0/24.0)*Mtmp63;
#pragma omp atomic
M[26] += -1.0/24.0*Mtmp0*Mtmp56 + (1.0/24.0)*Mtmp53*Mtmp6 - 1.0/24.0*Mtmp59*S[2] + (1.0/24.0)*Mtmp64 + (1.0/24.0)*Mtmp65;
#pragma omp atomic
M[27] += -1.0/12.0*Mtmp1*Mtmp66 - 1.0/12.0*Mtmp10*Mtmp67 + (1.0/4.0)*Mtmp11*Mtmp23*S[0] - 1.0/12.0*Mtmp22*Mtmp52 + (1.0/6.0)*Mtmp28*y*S[1] + (1.0/3.0)*Mtmp30*x*S[2] - 1.0/12.0*Mtmp32*Mtmp51 - 1.0/12.0*Mtmp32*Mtmp68 + (1.0/12.0)*Mtmp54*S[0];
#pragma omp atomic
M[28] += (1.0/6.0)*Mtmp26*Mtmp34 + (1.0/6.0)*Mtmp69 + (1.0/6.0)*Mtmp70 + (1.0/6.0)*Mtmp73;
#pragma omp atomic
M[29] += -1.0/12.0*Mtmp1*Mtmp74 + (1.0/4.0)*Mtmp11*Mtmp23*S[1] - 1.0/12.0*Mtmp11*Mtmp75 - 1.0/12.0*Mtmp15*Mtmp52 - 1.0/12.0*Mtmp23*Mtmp75 + (1.0/3.0)*Mtmp30*y*S[2] + (1.0/6.0)*Mtmp45*x*S[0] + (1.0/12.0)*Mtmp54*S[1] - 1.0/2.0*Mtmp61;
#pragma omp atomic
M[30] += -1.0/12.0*Mtmp0*Mtmp76 + (1.0/2.0)*Mtmp11*Mtmp27 - 1.0/12.0*Mtmp11*Mtmp77 + (1.0/12.0)*Mtmp19*Mtmp67 + (1.0/12.0)*Mtmp23*Mtmp34*S[2] - 1.0/12.0*Mtmp23*Mtmp77 + (1.0/12.0)*Mtmp65 - 1.0/12.0*Mtmp7*Mtmp76;
#pragma omp atomic
M[31] += -1.0/2.0*Mtmp10*Mtmp23 - 1.0/24.0*Mtmp22*Mtmp60 + (1.0/24.0)*Mtmp3*Mtmp80 - 1.0/24.0*Mtmp52*Mtmp68 + (1.0/24.0)*Mtmp57 + (1.0/24.0)*Mtmp79;
#pragma omp atomic
M[32] += (1.0/6.0)*Mtmp25*Mtmp41 + (1.0/6.0)*Mtmp73 + (1.0/6.0)*Mtmp81 + (1.0/6.0)*Mtmp82;
#pragma omp atomic
M[33] += -1.0/24.0*Mtmp1*Mtmp80 + (1.0/24.0)*Mtmp63 + (1.0/24.0)*Mtmp83 - 1.0/24.0*Mtmp84*S[1];
#pragma omp atomic
M[34] += -1.0/24.0*Mtmp56*Mtmp7 + (1.0/24.0)*Mtmp65 + (1.0/24.0)*Mtmp80*Mtmp9 - 1.0/24.0*Mtmp84*S[2] + (1.0/24.0)*Mtmp85;
#pragma omp atomic
M[35] += (1.0/120.0)*Mtmp86*S[0] + (1.0/120.0)*Mtmp89 + (1.0/120.0)*Mtmp92 + (1.0/120.0)*Mtmp95;
#pragma omp atomic
M[36] += (1.0/120.0)*Mtmp101 - 1.0/120.0*Mtmp16*Mtmp99 + (1.0/120.0)*Mtmp4*Mtmp91 - 1.0/120.0*Mtmp4*Mtmp98 + (1.0/120.0)*Mtmp86*S[1] - 1.0/12.0*Mtmp96;
#pragma omp atomic
M[37] += (1.0/120.0)*Mtmp102 - 1.0/120.0*Mtmp38*Mtmp90 - 1.0/120.0*Mtmp39*Mtmp94 + (1.0/120.0)*Mtmp6*Mtmp91 + (1.0/120.0)*Mtmp86*S[2];
#pragma omp atomic
M[38] += (1.0/120.0)*Mtmp0*Mtmp103 + (1.0/120.0)*Mtmp104*Mtmp29 + (1.0/120.0)*Mtmp105*Mtmp11 + (1.0/120.0)*Mtmp106 + (1.0/120.0)*Mtmp110 + (1.0/120.0)*Mtmp7*Mtmp91 + (1.0/120.0)*Mtmp92;
#pragma omp atomic
M[39] += -1.0/120.0*Mtmp100*Mtmp15 + (1.0/120.0)*Mtmp111 + (1.0/120.0)*Mtmp26*Mtmp99 - 1.0/120.0*Mtmp43*Mtmp94 + (1.0/120.0)*Mtmp8*Mtmp91 - 1.0/120.0*Mtmp8*Mtmp98 + (1.0/120.0)*Mtmp9*Mtmp91;
#pragma omp atomic
M[40] += -1.0/12.0*Mtmp10*Mtmp74 - 1.0/12.0*Mtmp11*Mtmp32*Mtmp4 + (1.0/12.0)*Mtmp11*Mtmp45*S[0] - 1.0/12.0*Mtmp112 - 1.0/12.0*Mtmp16*Mtmp66 + (1.0/12.0)*Mtmp23*Mtmp28*S[1] - 1.0/12.0*Mtmp23*Mtmp3*Mtmp32 + (1.0/3.0)*Mtmp30*x*y*S[2] + (1.0/12.0)*Mtmp54*x*S[1] + (1.0/12.0)*Mtmp54*y*S[0] - 1.0/12.0*Mtmp96;
#pragma omp atomic
M[41] += (1.0/60.0)*Mtmp102 - 1.0/60.0*Mtmp11*Mtmp113 - 1.0/60.0*Mtmp113*Mtmp23 - 1.0/60.0*Mtmp114*Mtmp12 + (1.0/60.0)*Mtmp114*Mtmp23 - 1.0/60.0*Mtmp115*Mtmp116 - 1.0/60.0*Mtmp117*Mtmp22 + (1.0/60.0)*Mtmp118*Mtmp27 + (1.0/60.0)*Mtmp120*Mtmp6;
#pragma omp atomic
M[42] += (1.0/120.0)*Mtmp0*Mtmp121 + (1.0/120.0)*Mtmp103*Mtmp7 + (1.0/120.0)*Mtmp105*Mtmp23 + (1.0/120.0)*Mtmp110 + (1.0/120.0)*Mtmp122 + (1.0/120.0)*Mtmp47*Mtmp94 + (1.0/120.0)*Mtmp95;
#pragma omp atomic
M[43] += -1.0/60.0*Mtmp11*Mtmp123 + (1.0/60.0)*Mtmp11*Mtmp124 + (1.0/60.0)*Mtmp111 - 1.0/60.0*Mtmp115*Mtmp125 - 1.0/60.0*Mtmp117*Mtmp15 - 1.0/60.0*Mtmp12*Mtmp124 + (1.0/60.0)*Mtmp120*Mtmp9 - 1.0/60.0*Mtmp123*Mtmp23 + (1.0/60.0)*Mtmp126*Mtmp19;
#pragma omp atomic
M[44] += -1.0/120.0*Mtmp10*Mtmp128 + (1.0/120.0)*Mtmp101 - 1.0/120.0*Mtmp108*Mtmp3 - 1.0/12.0*Mtmp112 + (1.0/120.0)*Mtmp121*Mtmp3 + (1.0/120.0)*Mtmp127*S[0];
#pragma omp atomic
M[45] += -1.0/120.0*Mtmp100*Mtmp22 + (1.0/120.0)*Mtmp102 - 1.0/120.0*Mtmp104*Mtmp39 - 1.0/120.0*Mtmp108*Mtmp5 + (1.0/120.0)*Mtmp121*Mtmp5 + (1.0/120.0)*Mtmp121*Mtmp6 + (1.0/120.0)*Mtmp128*Mtmp25;
#pragma omp atomic
M[46] += (1.0/120.0)*Mtmp106 + (1.0/120.0)*Mtmp122 + (1.0/120.0)*Mtmp127*S[1] + (1.0/120.0)*Mtmp89;
#pragma omp atomic
M[47] += -1.0/120.0*Mtmp104*Mtmp43 + (1.0/120.0)*Mtmp111 + (1.0/120.0)*Mtmp121*Mtmp9 + (1.0/120.0)*Mtmp127*S[2] - 1.0/120.0*Mtmp48*Mtmp90;
#pragma omp atomic
M[48] += -1.0/720.0*Mtmp130 + (1.0/720.0)*Mtmp131*S[0] - 1.0/720.0*Mtmp133 - 1.0/720.0*Mtmp135 - 1.0/720.0*Mtmp136 - 1.0/720.0*Mtmp137;
#pragma omp atomic
M[49] += (1.0/720.0)*Mtmp131*S[1] + (1.0/720.0)*Mtmp134*Mtmp4 - 1.0/720.0*Mtmp138 - 1.0/720.0*Mtmp139 - 1.0/720.0*Mtmp140 - 1.0/720.0*Mtmp143 - 1.0/720.0*Mtmp145 + (1.0/720.0)*Mtmp148;
#pragma omp atomic
M[50] += (1.0/720.0)*Mtmp131*S[2] + (1.0/720.0)*Mtmp134*Mtmp6 + (1.0/720.0)*Mtmp150 + (1.0/720.0)*Mtmp151 + (1.0/720.0)*Mtmp152;
#pragma omp atomic
M[51] += (1.0/720.0)*Mtmp134*Mtmp7 - 1.0/720.0*Mtmp135 - 1.0/720.0*Mtmp136 - 1.0/720.0*Mtmp153 - 1.0/720.0*Mtmp155 + (1.0/720.0)*Mtmp156*Mtmp50 + (1.0/720.0)*Mtmp157*Mtmp28 - 1.0/720.0*Mtmp159 - 1.0/720.0*Mtmp160 - 1.0/720.0*Mtmp161 + (1.0/720.0)*Mtmp162 + (1.0/720.0)*Mtmp55*Mtmp97;
#pragma omp atomic
M[52] += -1.0/360.0*Mtmp118*Mtmp43 + (1.0/360.0)*Mtmp163*Mtmp8 + (1.0/360.0)*Mtmp163*Mtmp9 - 1.0/360.0*Mtmp164*Mtmp69 + (1.0/360.0)*Mtmp166 + (1.0/24.0)*Mtmp26*Mtmp49 - 1.0/360.0*Mtmp72*Mtmp97;
#pragma omp atomic
M[53] += (1.0/6.0)*Mtmp11*Mtmp30*y*S[2] + (1.0/24.0)*Mtmp11*Mtmp54*S[1] - 1.0/720.0*Mtmp139 - 1.0/720.0*Mtmp143 - 1.0/720.0*Mtmp145 - 1.0/720.0*Mtmp167 - 1.0/720.0*Mtmp168 - 1.0/720.0*Mtmp171 - 1.0/720.0*Mtmp173 - 1.0/720.0*Mtmp175 - 1.0/720.0*Mtmp176 + (1.0/48.0)*Mtmp23*Mtmp49*S[1] + (1.0/36.0)*Mtmp28*Mtmp45*S[0] + (1.0/12.0)*Mtmp54*x*y*S[0];
#pragma omp atomic
M[54] += (1.0/720.0)*Mtmp0*Mtmp177 + (1.0/720.0)*Mtmp141*Mtmp23*Mtmp6 + (1.0/720.0)*Mtmp144*Mtmp27 + (1.0/720.0)*Mtmp151 + (1.0/720.0)*Mtmp156*Mtmp64 + (1.0/720.0)*Mtmp178 + (1.0/720.0)*Mtmp179 + (1.0/720.0)*Mtmp65*Mtmp97;
#pragma omp atomic
M[55] += (1.0/48.0)*Mtmp11*Mtmp78*S[0] - 1.0/720.0*Mtmp137 - 1.0/720.0*Mtmp153 - 1.0/720.0*Mtmp155 - 1.0/720.0*Mtmp159 - 1.0/720.0*Mtmp160 - 1.0/720.0*Mtmp161 - 1.0/720.0*Mtmp180 - 1.0/720.0*Mtmp182 - 1.0/720.0*Mtmp183 + (1.0/6.0)*Mtmp23*Mtmp30*x*S[2] + (1.0/24.0)*Mtmp23*Mtmp54*S[0] + (1.0/36.0)*Mtmp28*Mtmp45*S[1] + (1.0/12.0)*Mtmp54*x*y*S[1];
#pragma omp atomic
M[56] += -1.0/180.0*Mtmp113*Mtmp45 + (1.0/180.0)*Mtmp114*Mtmp45 - 1.0/180.0*Mtmp115*Mtmp69 - 1.0/180.0*Mtmp115*Mtmp81 - 1.0/180.0*Mtmp119*Mtmp72 + (1.0/180.0)*Mtmp119*Mtmp82 - 1.0/180.0*Mtmp123*Mtmp28 + (1.0/180.0)*Mtmp156*Mtmp70 - 1.0/180.0*Mtmp156*Mtmp71 + (1.0/180.0)*Mtmp166;
#pragma omp atomic
M[57] += (1.0/720.0)*Mtmp0*Mtmp185 + (1.0/720.0)*Mtmp107*Mtmp62 + (1.0/720.0)*Mtmp119*Mtmp83 + (1.0/720.0)*Mtmp148 + (1.0/720.0)*Mtmp157*Mtmp45 - 1.0/720.0*Mtmp167 - 1.0/720.0*Mtmp168 - 1.0/720.0*Mtmp171 - 1.0/720.0*Mtmp173 - 1.0/720.0*Mtmp175 - 1.0/720.0*Mtmp184 - 1.0/720.0*Mtmp186;
#pragma omp atomic
M[58] += (1.0/720.0)*Mtmp107*Mtmp65 + (1.0/720.0)*Mtmp119*Mtmp85 + (1.0/720.0)*Mtmp152 + (1.0/720.0)*Mtmp172*Mtmp9 + (1.0/720.0)*Mtmp177*Mtmp7 + (1.0/720.0)*Mtmp179 + (1.0/720.0)*Mtmp181*Mtmp19 + (1.0/720.0)*Mtmp187;
#pragma omp atomic
M[59] += -1.0/720.0*Mtmp130 - 1.0/720.0*Mtmp133 + (1.0/720.0)*Mtmp162 - 1.0/720.0*Mtmp180 - 1.0/720.0*Mtmp182 - 1.0/720.0*Mtmp183 + (1.0/720.0)*Mtmp185*Mtmp3 + (1.0/720.0)*Mtmp188*S[0];
#pragma omp atomic
M[60] += -1.0/360.0*Mtmp107*Mtmp71 - 1.0/360.0*Mtmp126*Mtmp39 - 1.0/360.0*Mtmp164*Mtmp81 + (1.0/360.0)*Mtmp166 + (1.0/360.0)*Mtmp189*Mtmp5 + (1.0/360.0)*Mtmp189*Mtmp6 + (1.0/24.0)*Mtmp25*Mtmp78;
#pragma omp atomic
M[61] += -1.0/720.0*Mtmp138 - 1.0/720.0*Mtmp140 - 1.0/720.0*Mtmp176 - 1.0/720.0*Mtmp184 - 1.0/720.0*Mtmp186 + (1.0/720.0)*Mtmp188*S[1];
#pragma omp atomic
M[62] += (1.0/720.0)*Mtmp150 + (1.0/720.0)*Mtmp178 + (1.0/720.0)*Mtmp185*Mtmp9 + (1.0/720.0)*Mtmp187 + (1.0/720.0)*Mtmp188*S[2];

}

void M2Mc_7(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = z*M[2];
double Mstmp2 = -Mstmp1;
double Mstmp3 = x*M[1];
double Mstmp4 = y*M[0];
double Mstmp5 = x*M[2];
double Mstmp6 = z*M[0];
double Mstmp7 = y*M[1];
double Mstmp8 = y*M[2];
double Mstmp9 = z*M[1];
double Mstmp10 = x*M[3];
double Mstmp11 = pow(x, 2);
double Mstmp12 = (1.0/2.0)*M[0];
double Mstmp13 = z*M[5];
double Mstmp14 = pow(z, 2);
double Mstmp15 = Mstmp1*x;
double Mstmp16 = -Mstmp12*Mstmp14 - Mstmp13 - Mstmp15;
double Mstmp17 = x*M[4];
double Mstmp18 = y*M[3];
double Mstmp19 = (1.0/2.0)*M[1];
double Mstmp20 = Mstmp0*y;
double Mstmp21 = z*M[7];
double Mstmp22 = Mstmp1*y;
double Mstmp23 = -Mstmp14*Mstmp19 - Mstmp21 - Mstmp22;
double Mstmp24 = x*M[5];
double Mstmp25 = z*M[3];
double Mstmp26 = Mstmp0*z;
double Mstmp27 = (1.0/2.0)*M[2];
double Mstmp28 = -Mstmp14*Mstmp27;
double Mstmp29 = x*M[6];
double Mstmp30 = y*M[4];
double Mstmp31 = pow(y, 2);
double Mstmp32 = Mstmp3*y;
double Mstmp33 = x*M[7];
double Mstmp34 = y*M[5];
double Mstmp35 = z*M[4];
double Mstmp36 = Mstmp5*y;
double Mstmp37 = Mstmp3*z;
double Mstmp38 = Mstmp4*z;
double Mstmp39 = y*M[6];
double Mstmp40 = y*M[7];
double Mstmp41 = z*M[6];
double Mstmp42 = Mstmp7*z;
double Mstmp43 = x*M[8];
double Mstmp44 = z*M[10];
double Mstmp45 = Mstmp13*x;
double Mstmp46 = (1.0/2.0)*M[3];
double Mstmp47 = pow(x, 3);
double Mstmp48 = (1.0/6.0)*Mstmp47;
double Mstmp49 = Mstmp14*Mstmp46;
double Mstmp50 = pow(z, 3);
double Mstmp51 = (1.0/6.0)*Mstmp50;
double Mstmp52 = Mstmp51*M[2];
double Mstmp53 = (1.0/2.0)*Mstmp14;
double Mstmp54 = Mstmp0*Mstmp53;
double Mstmp55 = (1.0/2.0)*Mstmp11;
double Mstmp56 = Mstmp1*Mstmp55;
double Mstmp57 = x*M[9];
double Mstmp58 = y*M[8];
double Mstmp59 = (1.0/2.0)*M[4];
double Mstmp60 = Mstmp10*y;
double Mstmp61 = z*M[12];
double Mstmp62 = Mstmp21*x;
double Mstmp63 = Mstmp13*y;
double Mstmp64 = -Mstmp14*Mstmp59 - Mstmp15*y - Mstmp3*Mstmp53 - Mstmp4*Mstmp53 - Mstmp61 - Mstmp62 - Mstmp63;
double Mstmp65 = x*M[10];
double Mstmp66 = z*M[8];
double Mstmp67 = (1.0/2.0)*M[5];
double Mstmp68 = Mstmp10*z;
double Mstmp69 = -Mstmp14*Mstmp67 - Mstmp5*Mstmp53 - Mstmp51*M[0];
double Mstmp70 = z*M[14];
double Mstmp71 = Mstmp21*y;
double Mstmp72 = (1.0/2.0)*M[6];
double Mstmp73 = Mstmp14*Mstmp72;
double Mstmp74 = Mstmp53*Mstmp7;
double Mstmp75 = (1.0/2.0)*Mstmp31;
double Mstmp76 = Mstmp1*Mstmp75;
double Mstmp77 = x*M[12];
double Mstmp78 = y*M[10];
double Mstmp79 = z*M[9];
double Mstmp80 = (1.0/2.0)*M[7];
double Mstmp81 = Mstmp24*y;
double Mstmp82 = Mstmp17*z;
double Mstmp83 = Mstmp18*z;
double Mstmp84 = -Mstmp14*Mstmp80 - Mstmp51*M[1] - Mstmp53*Mstmp8;
double Mstmp85 = x*M[13];
double Mstmp86 = y*M[11];
double Mstmp87 = pow(y, 3);
double Mstmp88 = (1.0/6.0)*Mstmp87;
double Mstmp89 = Mstmp29*y;
double Mstmp90 = x*M[14];
double Mstmp91 = y*M[12];
double Mstmp92 = z*M[11];
double Mstmp93 = Mstmp33*y;
double Mstmp94 = Mstmp29*z;
double Mstmp95 = Mstmp30*z;
double Mstmp96 = y*M[13];
double Mstmp97 = y*M[14];
double Mstmp98 = z*M[13];
double Mstmp99 = Mstmp39*z;
double Mstmp100 = x*M[15];
double Mstmp101 = (1.0/2.0)*M[8];
double Mstmp102 = z*M[17];
double Mstmp103 = Mstmp101*Mstmp14;
double Mstmp104 = pow(x, 4);
double Mstmp105 = (1.0/24.0)*M[0];
double Mstmp106 = Mstmp44*x;
double Mstmp107 = Mstmp10*Mstmp53;
double Mstmp108 = Mstmp13*Mstmp55;
double Mstmp109 = (1.0/4.0)*Mstmp14;
double Mstmp110 = Mstmp11*M[0];
double Mstmp111 = Mstmp109*Mstmp110;
double Mstmp112 = Mstmp1*Mstmp48;
double Mstmp113 = pow(z, 4);
double Mstmp114 = Mstmp105*Mstmp113 + Mstmp5*Mstmp51 + Mstmp51*M[5];
double Mstmp115 = x*M[16];
double Mstmp116 = y*M[15];
double Mstmp117 = (1.0/2.0)*M[9];
double Mstmp118 = z*M[19];
double Mstmp119 = Mstmp117*Mstmp14;
double Mstmp120 = (1.0/24.0)*M[1];
double Mstmp121 = Mstmp43*y;
double Mstmp122 = Mstmp61*x;
double Mstmp123 = Mstmp44*y;
double Mstmp124 = Mstmp17*Mstmp53;
double Mstmp125 = Mstmp18*Mstmp53;
double Mstmp126 = Mstmp21*Mstmp55;
double Mstmp127 = Mstmp109*Mstmp11;
double Mstmp128 = Mstmp127*M[1];
double Mstmp129 = Mstmp45*y;
double Mstmp130 = Mstmp20*Mstmp53;
double Mstmp131 = Mstmp22*Mstmp55;
double Mstmp132 = Mstmp113*Mstmp120 + Mstmp51*Mstmp8 + Mstmp51*M[7];
double Mstmp133 = x*M[17];
double Mstmp134 = z*M[15];
double Mstmp135 = (1.0/2.0)*M[10];
double Mstmp136 = (1.0/24.0)*M[2];
double Mstmp137 = Mstmp113*Mstmp136;
double Mstmp138 = Mstmp43*z;
double Mstmp139 = -Mstmp0*Mstmp51 - Mstmp127*M[2] - Mstmp135*Mstmp14 - Mstmp24*Mstmp53 - Mstmp51*M[3];
double Mstmp140 = z*M[21];
double Mstmp141 = Mstmp70*x;
double Mstmp142 = Mstmp61*y;
double Mstmp143 = Mstmp62*y;
double Mstmp144 = (1.0/2.0)*M[11];
double Mstmp145 = Mstmp14*Mstmp144;
double Mstmp146 = Mstmp29*Mstmp53;
double Mstmp147 = Mstmp30*Mstmp53;
double Mstmp148 = Mstmp13*Mstmp75;
double Mstmp149 = Mstmp32*Mstmp53;
double Mstmp150 = Mstmp15*Mstmp75;
double Mstmp151 = Mstmp109*Mstmp31;
double Mstmp152 = Mstmp151*M[0];
double Mstmp153 = x*M[19];
double Mstmp154 = y*M[17];
double Mstmp155 = z*M[16];
double Mstmp156 = (1.0/2.0)*M[12];
double Mstmp157 = Mstmp65*y;
double Mstmp158 = Mstmp57*z;
double Mstmp159 = Mstmp58*z;
double Mstmp160 = -Mstmp14*Mstmp156 - Mstmp3*Mstmp51 - Mstmp33*Mstmp53 - Mstmp34*Mstmp53 - Mstmp36*Mstmp53 - Mstmp4*Mstmp51 - Mstmp51*M[4];
double Mstmp161 = z*M[23];
double Mstmp162 = Mstmp70*y;
double Mstmp163 = (1.0/2.0)*M[13];
double Mstmp164 = Mstmp14*Mstmp163;
double Mstmp165 = Mstmp39*Mstmp53;
double Mstmp166 = Mstmp21*Mstmp75;
double Mstmp167 = Mstmp1*Mstmp88;
double Mstmp168 = Mstmp151*M[1];
double Mstmp169 = x*M[21];
double Mstmp170 = y*M[19];
double Mstmp171 = z*M[18];
double Mstmp172 = (1.0/2.0)*M[14];
double Mstmp173 = (1.0/12.0)*Mstmp113;
double Mstmp174 = Mstmp77*y;
double Mstmp175 = x*M[11];
double Mstmp176 = Mstmp175*z;
double Mstmp177 = y*M[9];
double Mstmp178 = Mstmp177*z;
double Mstmp179 = (1.0/4.0)*Mstmp11*Mstmp31;
double Mstmp180 = Mstmp17*y;
double Mstmp181 = -Mstmp14*Mstmp172 - Mstmp151*M[2] - Mstmp40*Mstmp53 - Mstmp51*Mstmp7 - Mstmp51*M[6];
double Mstmp182 = x*M[22];
double Mstmp183 = y*M[20];
double Mstmp184 = pow(y, 4);
double Mstmp185 = Mstmp85*y;
double Mstmp186 = x*M[23];
double Mstmp187 = y*M[21];
double Mstmp188 = z*M[20];
double Mstmp189 = Mstmp90*y;
double Mstmp190 = Mstmp85*z;
double Mstmp191 = Mstmp86*z;
double Mstmp192 = y*M[22];
double Mstmp193 = y*M[23];
double Mstmp194 = z*M[22];
double Mstmp195 = Mstmp96*z;
double Mstmp196 = x*M[24];
double Mstmp197 = z*M[26];
double Mstmp198 = Mstmp102*x;
double Mstmp199 = (1.0/2.0)*M[15];
double Mstmp200 = (1.0/24.0)*M[3];
double Mstmp201 = pow(x, 5);
double Mstmp202 = (1.0/120.0)*Mstmp201;
double Mstmp203 = Mstmp14*Mstmp199;
double Mstmp204 = pow(z, 5);
double Mstmp205 = (1.0/120.0)*M[2];
double Mstmp206 = -Mstmp204*Mstmp205;
double Mstmp207 = Mstmp43*Mstmp53;
double Mstmp208 = (1.0/24.0)*Mstmp113;
double Mstmp209 = Mstmp44*Mstmp55;
double Mstmp210 = Mstmp13*Mstmp48;
double Mstmp211 = (1.0/24.0)*Mstmp104;
double Mstmp212 = Mstmp1*Mstmp211;
double Mstmp213 = Mstmp127*M[3];
double Mstmp214 = (1.0/12.0)*Mstmp50;
double Mstmp215 = Mstmp11*M[2];
double Mstmp216 = (1.0/12.0)*Mstmp47;
double Mstmp217 = Mstmp14*M[0];
double Mstmp218 = Mstmp216*Mstmp217;
double Mstmp219 = x*M[25];
double Mstmp220 = y*M[24];
double Mstmp221 = (1.0/2.0)*M[16];
double Mstmp222 = z*M[28];
double Mstmp223 = Mstmp14*Mstmp221;
double Mstmp224 = (1.0/24.0)*M[4];
double Mstmp225 = Mstmp100*y;
double Mstmp226 = Mstmp118*x;
double Mstmp227 = Mstmp102*y;
double Mstmp228 = Mstmp53*Mstmp57;
double Mstmp229 = Mstmp53*Mstmp58;
double Mstmp230 = Mstmp55*Mstmp61;
double Mstmp231 = Mstmp127*M[4];
double Mstmp232 = Mstmp21*Mstmp48;
double Mstmp233 = Mstmp14*Mstmp216;
double Mstmp234 = Mstmp233*M[1];
double Mstmp235 = Mstmp106*y;
double Mstmp236 = Mstmp53*Mstmp60;
double Mstmp237 = Mstmp55*Mstmp63;
double Mstmp238 = Mstmp127*Mstmp4;
double Mstmp239 = Mstmp22*Mstmp48;
double Mstmp240 = Mstmp113*Mstmp224 + Mstmp208*Mstmp3 + Mstmp208*Mstmp4 + Mstmp33*Mstmp51 + Mstmp34*Mstmp51 + Mstmp36*Mstmp51 + Mstmp51*M[12];
double Mstmp241 = x*M[26];
double Mstmp242 = (1.0/2.0)*M[17];
double Mstmp243 = (1.0/24.0)*M[5];
double Mstmp244 = (1.0/120.0)*Mstmp204;
double Mstmp245 = Mstmp113*Mstmp243 + Mstmp208*Mstmp5 + Mstmp244*M[0];
double Mstmp246 = -Mstmp10*Mstmp51 - Mstmp110*Mstmp214 - Mstmp127*M[5] - Mstmp14*Mstmp242 - Mstmp233*M[2] - Mstmp51*M[8] - Mstmp53*Mstmp65;
double Mstmp247 = z*M[30];
double Mstmp248 = (1.0/40.0)*Mstmp204;
double Mstmp249 = Mstmp140*x;
double Mstmp250 = Mstmp0*Mstmp151 + Mstmp1*Mstmp179 + Mstmp118*y + Mstmp122*y + Mstmp127*Mstmp7 + Mstmp127*M[6] + (1.0/2.0)*Mstmp14*M[18] + Mstmp151*M[3] + Mstmp175*Mstmp53 + Mstmp177*Mstmp53 + Mstmp180*Mstmp53 + Mstmp247 + Mstmp248*M[2] + Mstmp249 + Mstmp44*Mstmp75 + Mstmp45*Mstmp75 + Mstmp55*Mstmp70 + Mstmp55*Mstmp71;
double Mstmp251 = x*M[28];
double Mstmp252 = (1.0/2.0)*M[19];
double Mstmp253 = (1.0/24.0)*M[7];
double Mstmp254 = Mstmp113*Mstmp253 + Mstmp208*Mstmp8 + Mstmp244*M[1];
double Mstmp255 = Mstmp11*M[1];
double Mstmp256 = -Mstmp127*Mstmp8 - Mstmp127*M[7] - Mstmp14*Mstmp252 - Mstmp17*Mstmp51 - Mstmp18*Mstmp51 - Mstmp20*Mstmp51 - Mstmp214*Mstmp255 - Mstmp51*M[9] - Mstmp53*Mstmp77 - Mstmp53*Mstmp78 - Mstmp53*Mstmp81;
double Mstmp257 = z*M[32];
double Mstmp258 = Mstmp161*x;
double Mstmp259 = Mstmp140*y;
double Mstmp260 = Mstmp141*y;
double Mstmp261 = (1.0/2.0)*M[20];
double Mstmp262 = Mstmp14*Mstmp261;
double Mstmp263 = Mstmp53*Mstmp85;
double Mstmp264 = Mstmp53*Mstmp86;
double Mstmp265 = Mstmp61*Mstmp75;
double Mstmp266 = Mstmp13*Mstmp88;
double Mstmp267 = Mstmp53*Mstmp89;
double Mstmp268 = Mstmp62*Mstmp75;
double Mstmp269 = Mstmp15*Mstmp88;
double Mstmp270 = Mstmp151*M[4];
double Mstmp271 = (1.0/12.0)*Mstmp87;
double Mstmp272 = Mstmp217*Mstmp271;
double Mstmp273 = Mstmp151*Mstmp3;
double Mstmp274 = x*M[30];
double Mstmp275 = (1.0/2.0)*M[21];
double Mstmp276 = (1.0/60.0)*Mstmp204;
double Mstmp277 = x*M[18];
double Mstmp278 = y*M[16];
double Mstmp279 = Mstmp216*Mstmp31;
double Mstmp280 = Mstmp57*y;
double Mstmp281 = Mstmp214*Mstmp31;
double Mstmp282 = -Mstmp14*Mstmp275 - Mstmp151*Mstmp5 - Mstmp151*M[5] - Mstmp281*M[0] - Mstmp29*Mstmp51 - Mstmp30*Mstmp51 - Mstmp32*Mstmp51 - Mstmp51*M[11] - Mstmp53*Mstmp90 - Mstmp53*Mstmp91 - Mstmp53*Mstmp93;
double Mstmp283 = z*M[34];
double Mstmp284 = (1.0/2.0)*M[22];
double Mstmp285 = Mstmp14*Mstmp284;
double Mstmp286 = Mstmp161*y;
double Mstmp287 = Mstmp53*Mstmp96;
double Mstmp288 = Mstmp70*Mstmp75;
double Mstmp289 = Mstmp151*M[6];
double Mstmp290 = Mstmp21*Mstmp88;
double Mstmp291 = Mstmp14*Mstmp271;
double Mstmp292 = Mstmp291*M[1];
double Mstmp293 = (1.0/24.0)*Mstmp184;
double Mstmp294 = Mstmp1*Mstmp293;
double Mstmp295 = x*M[32];
double Mstmp296 = (1.0/2.0)*M[23];
double Mstmp297 = x*M[20];
double Mstmp298 = y*M[18];
double Mstmp299 = Mstmp175*y;
double Mstmp300 = -Mstmp14*Mstmp296 - Mstmp151*M[7] - Mstmp281*M[1] - Mstmp291*M[2] - Mstmp39*Mstmp51 - Mstmp51*M[13] - Mstmp53*Mstmp97;
double Mstmp301 = x*M[33];
double Mstmp302 = y*M[31];
double Mstmp303 = pow(y, 5);
double Mstmp304 = (1.0/120.0)*Mstmp303;
double Mstmp305 = Mstmp182*y;
double Mstmp306 = x*M[34];
double Mstmp307 = y*M[33];
double Mstmp308 = (1.0/24.0)*M[6];
double Mstmp309 = (1.0/2.0)*M[24];
double Mstmp310 = z*M[37];
double Mstmp311 = Mstmp14*Mstmp309;
double Mstmp312 = (1.0/24.0)*M[8];
double Mstmp313 = pow(x, 6);
double Mstmp314 = (1.0/720.0)*M[0];
double Mstmp315 = Mstmp197*x;
double Mstmp316 = Mstmp100*Mstmp53;
double Mstmp317 = Mstmp102*Mstmp55;
double Mstmp318 = Mstmp127*M[8];
double Mstmp319 = Mstmp44*Mstmp48;
double Mstmp320 = Mstmp233*M[3];
double Mstmp321 = Mstmp13*Mstmp211;
double Mstmp322 = (1.0/48.0)*Mstmp104;
double Mstmp323 = Mstmp217*Mstmp322;
double Mstmp324 = Mstmp1*Mstmp202;
double Mstmp325 = Mstmp11*Mstmp214;
double Mstmp326 = (1.0/36.0)*Mstmp47;
double Mstmp327 = Mstmp50*M[2];
double Mstmp328 = (1.0/48.0)*Mstmp113;
double Mstmp329 = pow(z, 6);
double Mstmp330 = -Mstmp244*Mstmp5 - Mstmp244*M[5] - Mstmp314*Mstmp329;
double Mstmp331 = (1.0/2.0)*M[25];
double Mstmp332 = z*M[39];
double Mstmp333 = Mstmp14*Mstmp331;
double Mstmp334 = (1.0/24.0)*M[9];
double Mstmp335 = (1.0/720.0)*M[1];
double Mstmp336 = Mstmp222*x;
double Mstmp337 = Mstmp197*y;
double Mstmp338 = Mstmp115*Mstmp53;
double Mstmp339 = Mstmp116*Mstmp53;
double Mstmp340 = Mstmp118*Mstmp55;
double Mstmp341 = Mstmp127*M[9];
double Mstmp342 = Mstmp48*Mstmp61;
double Mstmp343 = Mstmp233*M[4];
double Mstmp344 = Mstmp21*Mstmp211;
double Mstmp345 = Mstmp14*Mstmp322;
double Mstmp346 = Mstmp345*M[1];
double Mstmp347 = Mstmp198*y;
double Mstmp348 = Mstmp121*Mstmp53;
double Mstmp349 = Mstmp123*Mstmp55;
double Mstmp350 = Mstmp127*Mstmp18;
double Mstmp351 = Mstmp48*Mstmp63;
double Mstmp352 = Mstmp233*Mstmp4;
double Mstmp353 = Mstmp211*Mstmp22;
double Mstmp354 = -Mstmp244*Mstmp8 - Mstmp244*M[7] - Mstmp329*Mstmp335;
double Mstmp355 = (1.0/2.0)*M[26];
double Mstmp356 = (1.0/720.0)*M[2];
double Mstmp357 = -Mstmp329*Mstmp356;
double Mstmp358 = (1.0/24.0)*M[10];
double Mstmp359 = Mstmp0*Mstmp244 + Mstmp113*Mstmp358 + Mstmp208*Mstmp24 + Mstmp215*Mstmp328 + Mstmp244*M[3];
double Mstmp360 = Mstmp326*Mstmp50;
double Mstmp361 = -Mstmp127*M[10] - Mstmp133*Mstmp53 - Mstmp14*Mstmp355 - Mstmp233*M[5] - Mstmp325*M[3] - Mstmp345*M[2] - Mstmp360*M[0] - Mstmp43*Mstmp51 - Mstmp51*M[15];
double Mstmp362 = (1.0/240.0)*Mstmp329;
double Mstmp363 = (1.0/8.0)*Mstmp14*Mstmp31;
double Mstmp364 = Mstmp1*Mstmp279 + Mstmp10*Mstmp151 + Mstmp102*Mstmp75 + Mstmp106*Mstmp75 + Mstmp110*Mstmp363 + Mstmp127*Mstmp30 + Mstmp127*M[11] + Mstmp13*Mstmp179 + (1.0/2.0)*Mstmp14*M[27] + Mstmp140*Mstmp55 + Mstmp142*Mstmp55 + Mstmp151*M[8] + Mstmp222*y + Mstmp226*y + Mstmp233*Mstmp7 + Mstmp233*M[6] + Mstmp247*x + Mstmp248*Mstmp5 + Mstmp248*M[5] + Mstmp277*Mstmp53 + Mstmp278*Mstmp53 + Mstmp280*Mstmp53 + Mstmp362*M[0] + Mstmp48*Mstmp70 + Mstmp48*Mstmp71 + z*M[41];
double Mstmp365 = (1.0/2.0)*M[28];
double Mstmp366 = (1.0/24.0)*M[12];
double Mstmp367 = Mstmp113*Mstmp366 + Mstmp208*Mstmp33 + Mstmp208*Mstmp34 + Mstmp208*Mstmp36 + Mstmp244*Mstmp3 + Mstmp244*Mstmp4 + Mstmp244*M[4];
double Mstmp368 = -Mstmp127*Mstmp34 - Mstmp127*M[12] - Mstmp14*Mstmp365 - Mstmp153*Mstmp53 - Mstmp154*Mstmp53 - Mstmp157*Mstmp53 - Mstmp233*Mstmp8 - Mstmp233*M[7] - Mstmp325*Mstmp4 - Mstmp325*M[4] - Mstmp360*M[1] - Mstmp51*Mstmp57 - Mstmp51*Mstmp58 - Mstmp51*Mstmp60 - Mstmp51*M[16];
double Mstmp369 = Mstmp11*Mstmp271;
double Mstmp370 = Mstmp0*Mstmp291 + Mstmp1*Mstmp369 + Mstmp118*Mstmp75 + Mstmp122*Mstmp75 + Mstmp127*Mstmp39 + Mstmp127*M[13] + (1.0/2.0)*Mstmp14*M[29] + Mstmp151*Mstmp17 + Mstmp151*M[9] + Mstmp161*Mstmp55 + Mstmp162*Mstmp55 + Mstmp179*Mstmp21 + Mstmp247*y + Mstmp248*Mstmp8 + Mstmp248*M[7] + Mstmp249*y + Mstmp255*Mstmp363 + Mstmp257*x + Mstmp291*M[3] + Mstmp297*Mstmp53 + Mstmp298*Mstmp53 + Mstmp299*Mstmp53 + Mstmp362*M[1] + Mstmp44*Mstmp88 + Mstmp45*Mstmp88 + z*M[43];
double Mstmp371 = (1.0/2.0)*M[30];
double Mstmp372 = (1.0/24.0)*M[14];
double Mstmp373 = Mstmp31*M[2];
double Mstmp374 = Mstmp113*Mstmp372 + Mstmp208*Mstmp40 + Mstmp244*Mstmp7 + Mstmp244*M[6] + Mstmp328*Mstmp373;
double Mstmp375 = -Mstmp0*Mstmp281 - Mstmp127*Mstmp40 - Mstmp127*M[14] - Mstmp14*Mstmp371 - Mstmp151*Mstmp24 - Mstmp151*M[10] - Mstmp169*Mstmp53 - Mstmp170*Mstmp53 - Mstmp174*Mstmp53 - Mstmp175*Mstmp51 - Mstmp177*Mstmp51 - Mstmp180*Mstmp51 - Mstmp215*Mstmp363 - Mstmp281*M[3] - Mstmp325*Mstmp7 - Mstmp325*M[6] - Mstmp362*M[2] - Mstmp51*M[18];
double Mstmp376 = z*M[45];
double Mstmp377 = (1.0/2.0)*M[31];
double Mstmp378 = Mstmp14*Mstmp377;
double Mstmp379 = Mstmp283*x;
double Mstmp380 = Mstmp257*y;
double Mstmp381 = Mstmp182*Mstmp53;
double Mstmp382 = Mstmp183*Mstmp53;
double Mstmp383 = Mstmp140*Mstmp75;
double Mstmp384 = Mstmp151*M[11];
double Mstmp385 = Mstmp61*Mstmp88;
double Mstmp386 = Mstmp291*M[4];
double Mstmp387 = Mstmp13*Mstmp293;
double Mstmp388 = (1.0/48.0)*Mstmp184;
double Mstmp389 = Mstmp217*Mstmp388;
double Mstmp390 = Mstmp258*y;
double Mstmp391 = Mstmp185*Mstmp53;
double Mstmp392 = Mstmp141*Mstmp75;
double Mstmp393 = Mstmp151*Mstmp29;
double Mstmp394 = Mstmp62*Mstmp88;
double Mstmp395 = Mstmp291*Mstmp3;
double Mstmp396 = Mstmp15*Mstmp293;
double Mstmp397 = (1.0/2.0)*M[32];
double Mstmp398 = (1.0/36.0)*Mstmp87;
double Mstmp399 = Mstmp398*Mstmp50;
double Mstmp400 = -Mstmp14*Mstmp397 - Mstmp151*Mstmp33 - Mstmp151*M[12] - Mstmp186*Mstmp53 - Mstmp187*Mstmp53 - Mstmp189*Mstmp53 - Mstmp281*Mstmp3 - Mstmp281*M[4] - Mstmp291*Mstmp5 - Mstmp291*M[5] - Mstmp399*M[0] - Mstmp51*Mstmp85 - Mstmp51*Mstmp86 - Mstmp51*Mstmp89 - Mstmp51*M[20];
double Mstmp401 = z*M[47];
double Mstmp402 = (1.0/2.0)*M[33];
double Mstmp403 = Mstmp14*Mstmp402;
double Mstmp404 = Mstmp283*y;
double Mstmp405 = Mstmp192*Mstmp53;
double Mstmp406 = Mstmp161*Mstmp75;
double Mstmp407 = Mstmp151*M[13];
double Mstmp408 = Mstmp70*Mstmp88;
double Mstmp409 = Mstmp291*M[6];
double Mstmp410 = Mstmp21*Mstmp293;
double Mstmp411 = Mstmp14*Mstmp388;
double Mstmp412 = Mstmp411*M[1];
double Mstmp413 = Mstmp1*Mstmp304;
double Mstmp414 = (1.0/2.0)*M[34];
double Mstmp415 = -Mstmp14*Mstmp414 - Mstmp151*M[14] - Mstmp193*Mstmp53 - Mstmp281*M[6] - Mstmp291*M[7] - Mstmp399*M[1] - Mstmp411*M[2] - Mstmp51*Mstmp96 - Mstmp51*M[22];
double Mstmp416 = (1.0/24.0)*M[11];
double Mstmp417 = pow(y, 6);
double Mstmp418 = (1.0/24.0)*M[13];
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += Mstmp0 + Mstmp2 + M[3];
#pragma omp atomic
Ms[4] += Mstmp3 + Mstmp4 + M[4];
#pragma omp atomic
Ms[5] += Mstmp5 + Mstmp6 + M[5];
#pragma omp atomic
Ms[6] += Mstmp2 + Mstmp7 + M[6];
#pragma omp atomic
Ms[7] += Mstmp8 + Mstmp9 + M[7];
#pragma omp atomic
Ms[8] += Mstmp10 + Mstmp11*Mstmp12 + Mstmp16 + M[8];
#pragma omp atomic
Ms[9] += Mstmp11*Mstmp19 + Mstmp17 + Mstmp18 + Mstmp20 + Mstmp23 + M[9];
#pragma omp atomic
Ms[10] += Mstmp11*Mstmp27 + Mstmp24 + Mstmp25 + Mstmp26 + Mstmp28 + M[10];
#pragma omp atomic
Ms[11] += Mstmp12*Mstmp31 + Mstmp16 + Mstmp29 + Mstmp30 + Mstmp32 + M[11];
#pragma omp atomic
Ms[12] += Mstmp33 + Mstmp34 + Mstmp35 + Mstmp36 + Mstmp37 + Mstmp38 + M[12];
#pragma omp atomic
Ms[13] += Mstmp19*Mstmp31 + Mstmp23 + Mstmp39 + M[13];
#pragma omp atomic
Ms[14] += Mstmp27*Mstmp31 + Mstmp28 + Mstmp40 + Mstmp41 + Mstmp42 + M[14];
#pragma omp atomic
Ms[15] += Mstmp11*Mstmp46 + Mstmp43 - Mstmp44 - Mstmp45 + Mstmp48*M[0] - Mstmp49 + Mstmp52 - Mstmp54 - Mstmp56 + M[15];
#pragma omp atomic
Ms[16] += Mstmp11*Mstmp59 + Mstmp4*Mstmp55 + Mstmp48*M[1] + Mstmp57 + Mstmp58 + Mstmp60 + Mstmp64 + M[16];
#pragma omp atomic
Ms[17] += Mstmp11*Mstmp67 + Mstmp48*M[2] + Mstmp55*Mstmp6 + Mstmp65 + Mstmp66 + Mstmp68 + Mstmp69 + M[17];
#pragma omp atomic
Ms[18] += (1.0/2.0)*Mstmp11*y*M[1] + (1.0/2.0)*Mstmp11*M[6] + (1.0/2.0)*Mstmp31*x*M[0] + (1.0/2.0)*Mstmp31*M[3] - Mstmp44 - Mstmp45 - Mstmp49 + (1.0/3.0)*Mstmp50*M[2] - Mstmp54 - Mstmp56 - Mstmp70 - Mstmp71 - Mstmp73 - Mstmp74 - Mstmp76 + x*y*M[4] + x*M[11] + y*M[9] + M[18];
#pragma omp atomic
Ms[19] += Mstmp11*Mstmp80 + Mstmp20*z + Mstmp55*Mstmp8 + Mstmp55*Mstmp9 + Mstmp77 + Mstmp78 + Mstmp79 + Mstmp81 + Mstmp82 + Mstmp83 + Mstmp84 + M[19];
#pragma omp atomic
Ms[20] += Mstmp3*Mstmp75 + Mstmp31*Mstmp59 + Mstmp64 + Mstmp85 + Mstmp86 + Mstmp88*M[0] + Mstmp89 + M[20];
#pragma omp atomic
Ms[21] += Mstmp31*Mstmp67 + Mstmp32*z + Mstmp5*Mstmp75 + Mstmp6*Mstmp75 + Mstmp69 + Mstmp90 + Mstmp91 + Mstmp92 + Mstmp93 + Mstmp94 + Mstmp95 + M[21];
#pragma omp atomic
Ms[22] += Mstmp31*Mstmp72 + Mstmp52 - Mstmp70 - Mstmp71 - Mstmp73 - Mstmp74 - Mstmp76 + Mstmp88*M[1] + Mstmp96 + M[22];
#pragma omp atomic
Ms[23] += Mstmp31*Mstmp80 + Mstmp75*Mstmp9 + Mstmp84 + Mstmp88*M[2] + Mstmp97 + Mstmp98 + Mstmp99 + M[23];
#pragma omp atomic
Ms[24] += Mstmp100 + Mstmp101*Mstmp11 - Mstmp102 - Mstmp103 + Mstmp104*Mstmp105 - Mstmp106 - Mstmp107 - Mstmp108 - Mstmp111 - Mstmp112 + Mstmp114 + Mstmp48*M[3] + M[24];
#pragma omp atomic
Ms[25] += Mstmp104*Mstmp120 + Mstmp11*Mstmp117 + Mstmp115 + Mstmp116 - Mstmp118 - Mstmp119 + Mstmp121 - Mstmp122 - Mstmp123 - Mstmp124 - Mstmp125 - Mstmp126 - Mstmp128 - Mstmp129 - Mstmp130 - Mstmp131 + Mstmp132 + Mstmp18*Mstmp55 + Mstmp4*Mstmp48 + Mstmp48*M[4] + M[25];
#pragma omp atomic
Ms[26] += Mstmp104*Mstmp136 + Mstmp11*Mstmp135 + Mstmp133 + Mstmp134 + Mstmp137 + Mstmp138 + Mstmp139 + Mstmp25*Mstmp55 + Mstmp48*Mstmp6 + Mstmp48*M[5] + M[26];
#pragma omp atomic
Ms[27] += -Mstmp102 - Mstmp103 - Mstmp106 - Mstmp107 - Mstmp108 + (1.0/4.0)*Mstmp11*Mstmp31*M[0] + (1.0/2.0)*Mstmp11*y*M[4] + (1.0/2.0)*Mstmp11*M[11] - Mstmp111 - Mstmp112 + (1.0/12.0)*Mstmp113*M[0] - Mstmp140 - Mstmp141 - Mstmp142 - Mstmp143 - Mstmp145 - Mstmp146 - Mstmp147 - Mstmp148 - Mstmp149 - Mstmp150 - Mstmp152 + (1.0/2.0)*Mstmp31*x*M[3] + (1.0/2.0)*Mstmp31*M[8] + (1.0/6.0)*Mstmp47*y*M[1] + (1.0/6.0)*Mstmp47*M[6] + (1.0/3.0)*Mstmp50*x*M[2] + (1.0/3.0)*Mstmp50*M[5] + x*y*M[9] + x*M[18] + y*M[16] + M[27];
#pragma omp atomic
Ms[28] += Mstmp11*Mstmp156 + Mstmp153 + Mstmp154 + Mstmp155 + Mstmp157 + Mstmp158 + Mstmp159 + Mstmp160 + Mstmp34*Mstmp55 + Mstmp35*Mstmp55 + Mstmp38*Mstmp55 + Mstmp48*Mstmp8 + Mstmp48*Mstmp9 + Mstmp48*M[7] + Mstmp60*z + M[28];
#pragma omp atomic
Ms[29] += (1.0/4.0)*Mstmp11*Mstmp31*M[1] + (1.0/2.0)*Mstmp11*y*M[6] + (1.0/2.0)*Mstmp11*M[13] + (1.0/12.0)*Mstmp113*M[1] - Mstmp118 - Mstmp119 - Mstmp122 - Mstmp123 - Mstmp124 - Mstmp125 - Mstmp126 - Mstmp128 - Mstmp129 - Mstmp130 - Mstmp131 - Mstmp161 - Mstmp162 - Mstmp164 - Mstmp165 - Mstmp166 - Mstmp167 - Mstmp168 + (1.0/2.0)*Mstmp31*x*M[4] + (1.0/2.0)*Mstmp31*M[9] + (1.0/3.0)*Mstmp50*y*M[2] + (1.0/3.0)*Mstmp50*M[7] + (1.0/6.0)*Mstmp87*x*M[0] + (1.0/6.0)*Mstmp87*M[3] + x*y*M[11] + x*M[20] + y*M[18] + M[29];
#pragma omp atomic
Ms[30] += Mstmp11*Mstmp172 + Mstmp135*Mstmp31 + Mstmp139 + Mstmp169 + Mstmp170 + Mstmp171 + Mstmp173*M[2] + Mstmp174 + Mstmp176 + Mstmp178 + Mstmp179*M[2] + Mstmp180*z + Mstmp181 + Mstmp24*Mstmp75 + Mstmp25*Mstmp75 + Mstmp26*Mstmp75 + Mstmp40*Mstmp55 + Mstmp41*Mstmp55 + Mstmp42*Mstmp55 + M[30];
#pragma omp atomic
Ms[31] += Mstmp105*Mstmp184 + Mstmp114 - Mstmp140 - Mstmp141 - Mstmp142 - Mstmp143 + Mstmp144*Mstmp31 - Mstmp145 - Mstmp146 - Mstmp147 - Mstmp148 - Mstmp149 - Mstmp150 - Mstmp152 + Mstmp182 + Mstmp183 + Mstmp185 + Mstmp29*Mstmp75 + Mstmp3*Mstmp88 + Mstmp88*M[4] + M[31];
#pragma omp atomic
Ms[32] += Mstmp156*Mstmp31 + Mstmp160 + Mstmp186 + Mstmp187 + Mstmp188 + Mstmp189 + Mstmp190 + Mstmp191 + Mstmp33*Mstmp75 + Mstmp35*Mstmp75 + Mstmp37*Mstmp75 + Mstmp5*Mstmp88 + Mstmp6*Mstmp88 + Mstmp88*M[5] + Mstmp89*z + M[32];
#pragma omp atomic
Ms[33] += Mstmp120*Mstmp184 + Mstmp132 - Mstmp161 - Mstmp162 + Mstmp163*Mstmp31 - Mstmp164 - Mstmp165 - Mstmp166 - Mstmp167 - Mstmp168 + Mstmp192 + Mstmp88*M[6] + M[33];
#pragma omp atomic
Ms[34] += Mstmp136*Mstmp184 + Mstmp137 + Mstmp172*Mstmp31 + Mstmp181 + Mstmp193 + Mstmp194 + Mstmp195 + Mstmp41*Mstmp75 + Mstmp88*Mstmp9 + Mstmp88*M[7] + M[34];
#pragma omp atomic
Ms[35] += Mstmp0*Mstmp208 + Mstmp104*Mstmp200 + Mstmp11*Mstmp199 + Mstmp113*Mstmp200 + Mstmp196 - Mstmp197 - Mstmp198 + Mstmp202*M[0] - Mstmp203 + Mstmp206 - Mstmp207 - Mstmp209 - Mstmp210 - Mstmp212 - Mstmp213 + Mstmp214*Mstmp215 - Mstmp218 + Mstmp24*Mstmp51 + Mstmp48*M[8] + Mstmp51*M[10] + M[35];
#pragma omp atomic
Ms[36] += Mstmp104*Mstmp224 + Mstmp11*Mstmp221 + Mstmp18*Mstmp48 + Mstmp202*M[1] + Mstmp211*Mstmp4 + Mstmp219 + Mstmp220 - Mstmp222 - Mstmp223 + Mstmp225 - Mstmp226 - Mstmp227 - Mstmp228 - Mstmp229 - Mstmp230 - Mstmp231 - Mstmp232 - Mstmp234 - Mstmp235 - Mstmp236 - Mstmp237 - Mstmp238 - Mstmp239 + Mstmp240 + Mstmp48*M[9] + Mstmp55*Mstmp58 + M[36];
#pragma omp atomic
Ms[37] += Mstmp100*z + Mstmp104*Mstmp243 + Mstmp11*Mstmp242 + Mstmp201*Mstmp205 + Mstmp211*Mstmp6 + Mstmp241 + Mstmp245 + Mstmp246 + Mstmp25*Mstmp48 + Mstmp48*M[10] + Mstmp55*Mstmp66 + z*M[24] + M[37];
#pragma omp atomic
Ms[38] += (1.0/24.0)*Mstmp104*y*M[1] + (1.0/24.0)*Mstmp104*M[6] + (1.0/4.0)*Mstmp11*Mstmp31*M[3] + (1.0/6.0)*Mstmp11*Mstmp50*M[2] + (1.0/2.0)*Mstmp11*y*M[9] + (1.0/2.0)*Mstmp11*M[18] + (1.0/12.0)*Mstmp113*x*M[0] + (1.0/24.0)*Mstmp113*y*M[1] + (1.0/12.0)*Mstmp113*M[3] + (1.0/24.0)*Mstmp113*M[6] - Mstmp197 - Mstmp198 - Mstmp203 - Mstmp207 - Mstmp209 - Mstmp210 - Mstmp212 - Mstmp213 - Mstmp218 - Mstmp250 + (1.0/12.0)*Mstmp31*Mstmp47*M[0] + (1.0/12.0)*Mstmp31*Mstmp50*M[2] + (1.0/2.0)*Mstmp31*x*M[8] + (1.0/2.0)*Mstmp31*M[15] + (1.0/6.0)*Mstmp47*y*M[4] + (1.0/6.0)*Mstmp47*M[11] + (1.0/3.0)*Mstmp50*x*M[5] + (1.0/6.0)*Mstmp50*y*M[7] + (1.0/3.0)*Mstmp50*M[10] + (1.0/6.0)*Mstmp50*M[14] + x*y*M[16] + x*M[27] + y*M[25] + M[38];
#pragma omp atomic
Ms[39] += Mstmp104*Mstmp253 + Mstmp11*Mstmp252 + Mstmp115*z + Mstmp116*z + Mstmp121*z + Mstmp133*y + Mstmp211*Mstmp8 + Mstmp211*Mstmp9 + Mstmp251 + Mstmp254 + Mstmp256 + Mstmp34*Mstmp48 + Mstmp35*Mstmp48 + Mstmp38*Mstmp48 + Mstmp48*M[12] + Mstmp55*Mstmp78 + Mstmp55*Mstmp79 + Mstmp55*Mstmp83 + y*M[26] + z*M[25] + M[39];
#pragma omp atomic
Ms[40] += (1.0/4.0)*Mstmp11*Mstmp31*M[4] + (1.0/12.0)*Mstmp11*Mstmp87*M[0] + (1.0/2.0)*Mstmp11*y*M[11] + (1.0/2.0)*Mstmp11*M[20] + (1.0/12.0)*Mstmp113*x*M[1] + (1.0/12.0)*Mstmp113*y*M[0] + (1.0/12.0)*Mstmp113*M[4] - Mstmp222 - Mstmp223 - Mstmp226 - Mstmp227 - Mstmp228 - Mstmp229 - Mstmp230 - Mstmp231 - Mstmp232 - Mstmp234 - Mstmp235 - Mstmp236 - Mstmp237 - Mstmp238 - Mstmp239 - Mstmp257 - Mstmp258 - Mstmp259 - Mstmp260 - Mstmp262 - Mstmp263 - Mstmp264 - Mstmp265 - Mstmp266 - Mstmp267 - Mstmp268 - Mstmp269 - Mstmp270 - Mstmp272 - Mstmp273 + (1.0/12.0)*Mstmp31*Mstmp47*M[1] + (1.0/2.0)*Mstmp31*x*M[9] + (1.0/2.0)*Mstmp31*M[16] + (1.0/6.0)*Mstmp47*y*M[6] + (1.0/6.0)*Mstmp47*M[13] + (1.0/3.0)*Mstmp50*x*y*M[2] + (1.0/3.0)*Mstmp50*x*M[7] + (1.0/3.0)*Mstmp50*y*M[5] + (1.0/3.0)*Mstmp50*M[12] + (1.0/6.0)*Mstmp87*x*M[3] + (1.0/6.0)*Mstmp87*M[8] + x*y*M[18] + x*M[29] + y*M[27] + M[40];
#pragma omp atomic
Ms[41] += Mstmp11*Mstmp275 + Mstmp153*y + Mstmp173*Mstmp5 + Mstmp173*M[5] + Mstmp179*Mstmp6 + Mstmp179*M[5] + Mstmp242*Mstmp31 + Mstmp246 + Mstmp274 + Mstmp276*M[0] + Mstmp277*z + Mstmp278*z + Mstmp279*M[2] + Mstmp280*z + Mstmp282 + Mstmp40*Mstmp48 + Mstmp41*Mstmp48 + Mstmp42*Mstmp48 + Mstmp48*M[14] + Mstmp55*Mstmp91 + Mstmp55*Mstmp92 + Mstmp55*Mstmp95 + Mstmp65*Mstmp75 + Mstmp66*Mstmp75 + Mstmp68*Mstmp75 + y*M[28] + z*M[27] + M[41];
#pragma omp atomic
Ms[42] += (1.0/4.0)*Mstmp11*Mstmp31*M[6] + (1.0/12.0)*Mstmp11*Mstmp50*M[2] + (1.0/12.0)*Mstmp11*Mstmp87*M[1] + (1.0/2.0)*Mstmp11*y*M[13] + (1.0/2.0)*Mstmp11*M[22] + (1.0/24.0)*Mstmp113*x*M[0] + (1.0/12.0)*Mstmp113*y*M[1] + (1.0/24.0)*Mstmp113*M[3] + (1.0/12.0)*Mstmp113*M[6] + (1.0/24.0)*Mstmp184*x*M[0] + (1.0/24.0)*Mstmp184*M[3] - Mstmp250 - Mstmp283 - Mstmp285 - Mstmp286 - Mstmp287 - Mstmp288 - Mstmp289 - Mstmp290 - Mstmp292 - Mstmp294 + (1.0/6.0)*Mstmp31*Mstmp50*M[2] + (1.0/2.0)*Mstmp31*x*M[11] + (1.0/2.0)*Mstmp31*M[18] + (1.0/6.0)*Mstmp50*x*M[5] + (1.0/3.0)*Mstmp50*y*M[7] + (1.0/6.0)*Mstmp50*M[10] + (1.0/3.0)*Mstmp50*M[14] + (1.0/6.0)*Mstmp87*x*M[4] + (1.0/6.0)*Mstmp87*M[9] + x*y*M[20] + x*M[31] + y*M[29] + M[42];
#pragma omp atomic
Ms[43] += Mstmp11*Mstmp296 + Mstmp169*y + Mstmp173*Mstmp8 + Mstmp173*M[7] + Mstmp179*Mstmp9 + Mstmp179*M[7] + Mstmp215*Mstmp271 + Mstmp24*Mstmp88 + Mstmp25*Mstmp88 + Mstmp252*Mstmp31 + Mstmp256 + Mstmp26*Mstmp88 + Mstmp276*M[1] + Mstmp295 + Mstmp297*z + Mstmp298*z + Mstmp299*z + Mstmp300 + Mstmp55*Mstmp97 + Mstmp55*Mstmp98 + Mstmp55*Mstmp99 + Mstmp75*Mstmp77 + Mstmp75*Mstmp79 + Mstmp75*Mstmp82 + Mstmp88*M[10] + y*M[30] + z*M[29] + M[43];
#pragma omp atomic
Ms[44] += Mstmp184*Mstmp224 + Mstmp240 - Mstmp257 - Mstmp258 - Mstmp259 - Mstmp260 + Mstmp261*Mstmp31 - Mstmp262 - Mstmp263 - Mstmp264 - Mstmp265 - Mstmp266 - Mstmp267 - Mstmp268 - Mstmp269 - Mstmp270 - Mstmp272 - Mstmp273 + Mstmp29*Mstmp88 + Mstmp293*Mstmp3 + Mstmp301 + Mstmp302 + Mstmp304*M[0] + Mstmp305 + Mstmp75*Mstmp85 + Mstmp88*M[11] + M[44];
#pragma omp atomic
Ms[45] += Mstmp182*z + Mstmp183*z + Mstmp184*Mstmp243 + Mstmp185*z + Mstmp186*y + Mstmp245 + Mstmp275*Mstmp31 + Mstmp282 + Mstmp293*Mstmp5 + Mstmp293*Mstmp6 + Mstmp306 + Mstmp33*Mstmp88 + Mstmp35*Mstmp88 + Mstmp37*Mstmp88 + Mstmp75*Mstmp90 + Mstmp75*Mstmp92 + Mstmp75*Mstmp94 + Mstmp88*M[12] + y*M[32] + z*M[31] + M[45];
#pragma omp atomic
Ms[46] += Mstmp113*Mstmp308 + Mstmp184*Mstmp308 + Mstmp206 + Mstmp208*Mstmp7 + Mstmp281*M[2] - Mstmp283 + Mstmp284*Mstmp31 - Mstmp285 - Mstmp286 - Mstmp287 - Mstmp288 - Mstmp289 - Mstmp290 - Mstmp292 - Mstmp294 + Mstmp304*M[1] + Mstmp307 + Mstmp40*Mstmp51 + Mstmp51*M[14] + Mstmp88*M[13] + M[46];
#pragma omp atomic
Ms[47] += Mstmp184*Mstmp253 + Mstmp192*z + Mstmp205*Mstmp303 + Mstmp254 + Mstmp293*Mstmp9 + Mstmp296*Mstmp31 + Mstmp300 + Mstmp41*Mstmp88 + Mstmp75*Mstmp98 + Mstmp88*M[14] + y*M[34] + z*M[33] + M[47];
#pragma omp atomic
Ms[48] += Mstmp10*Mstmp208 + Mstmp104*Mstmp312 + Mstmp11*Mstmp309 + Mstmp110*Mstmp328 + Mstmp113*Mstmp312 + Mstmp202*M[3] - Mstmp310 - Mstmp311 + Mstmp313*Mstmp314 - Mstmp315 - Mstmp316 - Mstmp317 - Mstmp318 - Mstmp319 - Mstmp320 - Mstmp321 - Mstmp323 - Mstmp324 + Mstmp325*M[5] + Mstmp326*Mstmp327 + Mstmp330 + Mstmp48*M[15] + Mstmp51*Mstmp65 + Mstmp51*M[17] + x*M[35] + M[48];
#pragma omp atomic
Ms[49] += Mstmp104*Mstmp334 + Mstmp11*Mstmp331 + Mstmp113*Mstmp334 + Mstmp116*Mstmp55 + Mstmp17*Mstmp208 + Mstmp18*Mstmp208 + Mstmp18*Mstmp211 + Mstmp196*y + Mstmp20*Mstmp208 + Mstmp202*Mstmp4 + Mstmp202*M[4] + Mstmp255*Mstmp328 + Mstmp313*Mstmp335 + Mstmp325*Mstmp8 + Mstmp325*M[7] - Mstmp332 - Mstmp333 - Mstmp336 - Mstmp337 - Mstmp338 - Mstmp339 - Mstmp340 - Mstmp341 - Mstmp342 - Mstmp343 - Mstmp344 - Mstmp346 - Mstmp347 - Mstmp348 - Mstmp349 - Mstmp350 - Mstmp351 - Mstmp352 - Mstmp353 + Mstmp354 + Mstmp48*Mstmp58 + Mstmp48*M[16] + Mstmp51*Mstmp77 + Mstmp51*Mstmp78 + Mstmp51*Mstmp81 + Mstmp51*M[19] + x*M[36] + y*M[35] + M[49];
#pragma omp atomic
Ms[50] += Mstmp104*Mstmp358 + Mstmp11*Mstmp355 + Mstmp134*Mstmp55 + Mstmp196*z + Mstmp202*Mstmp6 + Mstmp202*M[5] + Mstmp211*Mstmp25 + Mstmp313*Mstmp356 + Mstmp357 + Mstmp359 + Mstmp361 + Mstmp48*Mstmp66 + Mstmp48*M[17] + x*M[37] + z*M[35] + M[50];
#pragma omp atomic
Ms[51] += (1.0/48.0)*Mstmp104*Mstmp31*M[0] + (1.0/24.0)*Mstmp104*y*M[4] + (1.0/24.0)*Mstmp104*M[11] + (1.0/24.0)*Mstmp11*Mstmp113*M[0] + (1.0/4.0)*Mstmp11*Mstmp31*M[8] + (1.0/6.0)*Mstmp11*Mstmp50*M[5] + (1.0/2.0)*Mstmp11*y*M[16] + (1.0/2.0)*Mstmp11*M[27] + (1.0/48.0)*Mstmp113*Mstmp31*M[0] + (1.0/24.0)*Mstmp113*x*y*M[1] + (1.0/12.0)*Mstmp113*x*M[3] + (1.0/24.0)*Mstmp113*x*M[6] + (1.0/24.0)*Mstmp113*y*M[4] + (1.0/12.0)*Mstmp113*M[8] + (1.0/24.0)*Mstmp113*M[11] + (1.0/120.0)*Mstmp201*y*M[1] + (1.0/120.0)*Mstmp201*M[6] + (1.0/12.0)*Mstmp31*Mstmp47*M[3] + (1.0/12.0)*Mstmp31*Mstmp50*x*M[2] + (1.0/12.0)*Mstmp31*Mstmp50*M[5] + (1.0/2.0)*Mstmp31*x*M[15] + (1.0/2.0)*Mstmp31*M[24] - Mstmp310 - Mstmp311 - Mstmp315 - Mstmp316 - Mstmp317 - Mstmp318 - Mstmp319 - Mstmp320 - Mstmp321 - Mstmp323 - Mstmp324 - Mstmp364 + (1.0/18.0)*Mstmp47*Mstmp50*M[2] + (1.0/6.0)*Mstmp47*y*M[9] + (1.0/6.0)*Mstmp47*M[18] + (1.0/6.0)*Mstmp50*x*y*M[7] + (1.0/3.0)*Mstmp50*x*M[10] + (1.0/6.0)*Mstmp50*x*M[14] + (1.0/6.0)*Mstmp50*y*M[12] + (1.0/3.0)*Mstmp50*M[17] + (1.0/6.0)*Mstmp50*M[21] + x*y*M[25] + x*M[38] + y*M[36] + M[51];
#pragma omp atomic
Ms[52] += Mstmp104*Mstmp366 + Mstmp11*Mstmp365 + Mstmp154*Mstmp55 + Mstmp155*Mstmp55 + Mstmp159*Mstmp55 + Mstmp202*Mstmp8 + Mstmp202*Mstmp9 + Mstmp202*M[7] + Mstmp211*Mstmp34 + Mstmp211*Mstmp35 + Mstmp211*Mstmp38 + Mstmp219*z + Mstmp220*z + Mstmp225*z + Mstmp241*y + Mstmp367 + Mstmp368 + Mstmp48*Mstmp78 + Mstmp48*Mstmp79 + Mstmp48*Mstmp83 + Mstmp48*M[19] + x*M[39] + y*M[37] + z*M[36] + M[52];
#pragma omp atomic
Ms[53] += (1.0/48.0)*Mstmp104*Mstmp31*M[1] + (1.0/24.0)*Mstmp104*y*M[6] + (1.0/24.0)*Mstmp104*M[13] + (1.0/24.0)*Mstmp11*Mstmp113*M[1] + (1.0/4.0)*Mstmp11*Mstmp31*M[9] + (1.0/6.0)*Mstmp11*Mstmp50*y*M[2] + (1.0/6.0)*Mstmp11*Mstmp50*M[7] + (1.0/12.0)*Mstmp11*Mstmp87*M[3] + (1.0/2.0)*Mstmp11*y*M[18] + (1.0/2.0)*Mstmp11*M[29] + (1.0/48.0)*Mstmp113*Mstmp31*M[1] + (1.0/12.0)*Mstmp113*x*y*M[0] + (1.0/12.0)*Mstmp113*x*M[4] + (1.0/12.0)*Mstmp113*y*M[3] + (1.0/24.0)*Mstmp113*y*M[6] + (1.0/12.0)*Mstmp113*M[9] + (1.0/24.0)*Mstmp113*M[13] + (1.0/12.0)*Mstmp31*Mstmp47*M[4] + (1.0/12.0)*Mstmp31*Mstmp50*M[7] + (1.0/2.0)*Mstmp31*x*M[16] + (1.0/2.0)*Mstmp31*M[25] - Mstmp332 - Mstmp333 - Mstmp336 - Mstmp337 - Mstmp338 - Mstmp339 - Mstmp340 - Mstmp341 - Mstmp342 - Mstmp343 - Mstmp344 - Mstmp346 - Mstmp347 - Mstmp348 - Mstmp349 - Mstmp350 - Mstmp351 - Mstmp352 - Mstmp353 - Mstmp370 + (1.0/36.0)*Mstmp47*Mstmp87*M[0] + (1.0/6.0)*Mstmp47*y*M[11] + (1.0/6.0)*Mstmp47*M[20] + (1.0/36.0)*Mstmp50*Mstmp87*M[2] + (1.0/3.0)*Mstmp50*x*y*M[5] + (1.0/3.0)*Mstmp50*x*M[12] + (1.0/3.0)*Mstmp50*y*M[10] + (1.0/6.0)*Mstmp50*y*M[14] + (1.0/3.0)*Mstmp50*M[19] + (1.0/6.0)*Mstmp50*M[23] + (1.0/6.0)*Mstmp87*x*M[8] + (1.0/6.0)*Mstmp87*M[15] + x*y*M[27] + x*M[40] + y*M[38] + M[53];
#pragma omp atomic
Ms[54] += Mstmp0*Mstmp276 + Mstmp104*Mstmp372 + Mstmp11*Mstmp137 + Mstmp11*Mstmp371 + Mstmp115*y*z + Mstmp133*Mstmp75 + Mstmp134*Mstmp75 + Mstmp138*Mstmp75 + Mstmp170*Mstmp55 + Mstmp171*Mstmp55 + Mstmp173*Mstmp24 + Mstmp173*M[10] + Mstmp178*Mstmp55 + Mstmp179*Mstmp25 + Mstmp179*M[10] + Mstmp211*Mstmp40 + Mstmp211*Mstmp41 + Mstmp211*Mstmp42 + Mstmp251*y + Mstmp276*M[3] + Mstmp279*Mstmp6 + Mstmp279*M[5] + Mstmp31*Mstmp355 + Mstmp322*Mstmp373 + Mstmp361 + Mstmp374 + Mstmp375 + Mstmp48*Mstmp91 + Mstmp48*Mstmp92 + Mstmp48*Mstmp95 + Mstmp48*M[21] + x*z*M[27] + x*M[41] + y*z*M[25] + y*M[39] + z*M[38] + M[54];
#pragma omp atomic
Ms[55] += (1.0/48.0)*Mstmp11*Mstmp113*M[0] + (1.0/48.0)*Mstmp11*Mstmp184*M[0] + (1.0/4.0)*Mstmp11*Mstmp31*M[11] + (1.0/12.0)*Mstmp11*Mstmp50*M[5] + (1.0/12.0)*Mstmp11*Mstmp87*M[4] + (1.0/2.0)*Mstmp11*y*M[20] + (1.0/2.0)*Mstmp11*M[31] + (1.0/24.0)*Mstmp113*Mstmp31*M[0] + (1.0/12.0)*Mstmp113*x*y*M[1] + (1.0/24.0)*Mstmp113*x*M[3] + (1.0/12.0)*Mstmp113*x*M[6] + (1.0/12.0)*Mstmp113*y*M[4] + (1.0/24.0)*Mstmp113*M[8] + (1.0/12.0)*Mstmp113*M[11] + (1.0/24.0)*Mstmp184*x*M[3] + (1.0/24.0)*Mstmp184*M[8] + (1.0/12.0)*Mstmp31*Mstmp47*M[6] + (1.0/6.0)*Mstmp31*Mstmp50*x*M[2] + (1.0/6.0)*Mstmp31*Mstmp50*M[5] + (1.0/2.0)*Mstmp31*x*M[18] + (1.0/2.0)*Mstmp31*M[27] - Mstmp364 - Mstmp376 - Mstmp378 - Mstmp379 - Mstmp380 - Mstmp381 - Mstmp382 - Mstmp383 - Mstmp384 - Mstmp385 - Mstmp386 - Mstmp387 - Mstmp389 - Mstmp390 - Mstmp391 - Mstmp392 - Mstmp393 - Mstmp394 - Mstmp395 - Mstmp396 + (1.0/36.0)*Mstmp47*Mstmp50*M[2] + (1.0/36.0)*Mstmp47*Mstmp87*M[1] + (1.0/6.0)*Mstmp47*y*M[13] + (1.0/6.0)*Mstmp47*M[22] + (1.0/3.0)*Mstmp50*x*y*M[7] + (1.0/6.0)*Mstmp50*x*M[10] + (1.0/3.0)*Mstmp50*x*M[14] + (1.0/3.0)*Mstmp50*y*M[12] + (1.0/6.0)*Mstmp50*M[17] + (1.0/3.0)*Mstmp50*M[21] + (1.0/6.0)*Mstmp87*x*M[9] + (1.0/6.0)*Mstmp87*M[16] + x*y*M[29] + x*M[42] + y*M[40] + M[55];
#pragma omp atomic
Ms[56] += Mstmp11*Mstmp397 + Mstmp153*Mstmp75 + Mstmp155*Mstmp75 + Mstmp158*Mstmp75 + Mstmp173*Mstmp33 + Mstmp173*Mstmp34 + Mstmp173*Mstmp36 + Mstmp173*M[12] + Mstmp179*Mstmp35 + Mstmp179*M[12] + Mstmp187*Mstmp55 + Mstmp188*Mstmp55 + Mstmp191*Mstmp55 + Mstmp274*y + Mstmp276*Mstmp3 + Mstmp276*Mstmp4 + Mstmp276*M[4] + Mstmp277*y*z + Mstmp279*Mstmp9 + Mstmp279*M[7] + Mstmp31*Mstmp365 + Mstmp326*Mstmp87*M[2] + Mstmp368 + Mstmp369*Mstmp6 + Mstmp369*M[5] + Mstmp400 + Mstmp48*Mstmp97 + Mstmp48*Mstmp98 + Mstmp48*Mstmp99 + Mstmp48*M[23] + Mstmp65*Mstmp88 + Mstmp66*Mstmp88 + Mstmp68*Mstmp88 + Mstmp88*M[17] + x*z*M[29] + x*M[43] + y*z*M[27] + y*M[41] + z*M[40] + M[56];
#pragma omp atomic
Ms[57] += (1.0/48.0)*Mstmp11*Mstmp113*M[1] + (1.0/48.0)*Mstmp11*Mstmp184*M[1] + (1.0/4.0)*Mstmp11*Mstmp31*M[13] + (1.0/12.0)*Mstmp11*Mstmp50*y*M[2] + (1.0/12.0)*Mstmp11*Mstmp50*M[7] + (1.0/12.0)*Mstmp11*Mstmp87*M[6] + (1.0/2.0)*Mstmp11*y*M[22] + (1.0/2.0)*Mstmp11*M[33] + (1.0/24.0)*Mstmp113*Mstmp31*M[1] + (1.0/24.0)*Mstmp113*x*y*M[0] + (1.0/24.0)*Mstmp113*x*M[4] + (1.0/24.0)*Mstmp113*y*M[3] + (1.0/12.0)*Mstmp113*y*M[6] + (1.0/24.0)*Mstmp113*M[9] + (1.0/12.0)*Mstmp113*M[13] + (1.0/24.0)*Mstmp184*x*M[4] + (1.0/24.0)*Mstmp184*M[9] + (1.0/120.0)*Mstmp303*x*M[0] + (1.0/120.0)*Mstmp303*M[3] + (1.0/6.0)*Mstmp31*Mstmp50*M[7] + (1.0/2.0)*Mstmp31*x*M[20] + (1.0/2.0)*Mstmp31*M[29] - Mstmp370 - Mstmp401 - Mstmp403 - Mstmp404 - Mstmp405 - Mstmp406 - Mstmp407 - Mstmp408 - Mstmp409 - Mstmp410 - Mstmp412 - Mstmp413 + (1.0/18.0)*Mstmp50*Mstmp87*M[2] + (1.0/6.0)*Mstmp50*x*y*M[5] + (1.0/6.0)*Mstmp50*x*M[12] + (1.0/6.0)*Mstmp50*y*M[10] + (1.0/3.0)*Mstmp50*y*M[14] + (1.0/6.0)*Mstmp50*M[19] + (1.0/3.0)*Mstmp50*M[23] + (1.0/6.0)*Mstmp87*x*M[11] + (1.0/6.0)*Mstmp87*M[18] + x*y*M[31] + x*M[44] + y*M[42] + M[57];
#pragma omp atomic
Ms[58] += Mstmp11*Mstmp414 + Mstmp137*Mstmp31 + Mstmp169*Mstmp75 + Mstmp171*Mstmp75 + Mstmp173*Mstmp40 + Mstmp173*M[14] + Mstmp176*Mstmp75 + Mstmp179*Mstmp41 + Mstmp179*M[14] + Mstmp184*Mstmp358 + Mstmp193*Mstmp55 + Mstmp194*Mstmp55 + Mstmp195*Mstmp55 + Mstmp215*Mstmp388 + Mstmp24*Mstmp293 + Mstmp25*Mstmp293 + Mstmp26*Mstmp293 + Mstmp276*Mstmp7 + Mstmp276*M[6] + Mstmp295*y + Mstmp297*y*z + Mstmp31*Mstmp371 + Mstmp359 + Mstmp369*Mstmp9 + Mstmp369*M[7] + Mstmp375 + Mstmp415 + Mstmp77*Mstmp88 + Mstmp79*Mstmp88 + Mstmp82*Mstmp88 + Mstmp88*M[19] + x*z*M[31] + x*M[45] + y*z*M[29] + y*M[43] + z*M[42] + M[58];
#pragma omp atomic
Ms[59] += Mstmp113*Mstmp416 + Mstmp182*Mstmp75 + Mstmp184*Mstmp416 + Mstmp208*Mstmp29 + Mstmp208*Mstmp30 + Mstmp208*Mstmp32 + Mstmp281*Mstmp5 + Mstmp281*M[5] + Mstmp29*Mstmp293 + Mstmp3*Mstmp304 + Mstmp301*y + Mstmp304*M[4] + Mstmp31*Mstmp328*M[0] + Mstmp31*Mstmp377 + Mstmp314*Mstmp417 + Mstmp330 - Mstmp376 - Mstmp378 - Mstmp379 - Mstmp380 - Mstmp381 - Mstmp382 - Mstmp383 - Mstmp384 - Mstmp385 - Mstmp386 - Mstmp387 - Mstmp389 - Mstmp390 - Mstmp391 - Mstmp392 - Mstmp393 - Mstmp394 - Mstmp395 - Mstmp396 + Mstmp51*Mstmp90 + Mstmp51*Mstmp91 + Mstmp51*Mstmp93 + Mstmp51*M[21] + Mstmp85*Mstmp88 + Mstmp88*M[20] + x*M[46] + y*M[44] + M[59];
#pragma omp atomic
Ms[60] += Mstmp184*Mstmp366 + Mstmp186*Mstmp75 + Mstmp188*Mstmp75 + Mstmp190*Mstmp75 + Mstmp293*Mstmp33 + Mstmp293*Mstmp35 + Mstmp293*Mstmp37 + Mstmp301*z + Mstmp302*z + Mstmp304*Mstmp5 + Mstmp304*Mstmp6 + Mstmp304*M[5] + Mstmp305*z + Mstmp306*y + Mstmp31*Mstmp397 + Mstmp367 + Mstmp400 + Mstmp88*Mstmp90 + Mstmp88*Mstmp92 + Mstmp88*Mstmp94 + Mstmp88*M[21] + x*M[47] + y*M[45] + z*M[44] + M[60];
#pragma omp atomic
Ms[61] += Mstmp113*Mstmp418 + Mstmp184*Mstmp418 + Mstmp208*Mstmp39 + Mstmp281*M[7] + Mstmp304*M[6] + Mstmp31*Mstmp328*M[1] + Mstmp31*Mstmp402 + Mstmp327*Mstmp398 + Mstmp335*Mstmp417 + Mstmp354 - Mstmp401 - Mstmp403 - Mstmp404 - Mstmp405 - Mstmp406 - Mstmp407 - Mstmp408 - Mstmp409 - Mstmp410 - Mstmp412 - Mstmp413 + Mstmp51*Mstmp97 + Mstmp51*M[23] + Mstmp88*M[22] + y*M[46] + M[61];
#pragma omp atomic
Ms[62] += Mstmp184*Mstmp372 + Mstmp194*Mstmp75 + Mstmp293*Mstmp41 + Mstmp304*Mstmp9 + Mstmp304*M[7] + Mstmp307*z + Mstmp31*Mstmp414 + Mstmp356*Mstmp417 + Mstmp357 + Mstmp374 + Mstmp415 + Mstmp88*Mstmp98 + Mstmp88*M[23] + y*M[47] + z*M[46] + M[62];

}

void L2Lc_7(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[13];
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
double Lstmp19 = x*L[12];
double Lstmp20 = x*L[21];
double Lstmp21 = x*L[32];
double Lstmp22 = x*L[45];
double Lstmp23 = y*L[10];
double Lstmp24 = z*L[11];
double Lstmp25 = y*L[17];
double Lstmp26 = z*L[18];
double Lstmp27 = y*L[26];
double Lstmp28 = z*L[27];
double Lstmp29 = y*L[37];
double Lstmp30 = z*L[38];
double Lstmp31 = z*L[15];
double Lstmp32 = z*L[24];
double Lstmp33 = z*L[35];
double Lstmp34 = z*L[48];
double Lstmp35 = z*L[22];
double Lstmp36 = Lstmp35*x;
double Lstmp37 = z*L[33];
double Lstmp38 = Lstmp37*x;
double Lstmp39 = z*L[46];
double Lstmp40 = Lstmp39*x;
double Lstmp41 = z*L[20];
double Lstmp42 = Lstmp41*y;
double Lstmp43 = z*L[29];
double Lstmp44 = Lstmp43*y;
double Lstmp45 = z*L[40];
double Lstmp46 = Lstmp45*y;
double Lstmp47 = Lstmp12*Lstmp5;
double Lstmp48 = (1.0/4.0)*Lstmp47;
double Lstmp49 = (1.0/12.0)*Lstmp14*Lstmp5;
double Lstmp50 = (1.0/48.0)*Lstmp5;
double Lstmp51 = (1.0/12.0)*Lstmp12*Lstmp7;
double Lstmp52 = (1.0/48.0)*Lstmp12;
double Lstmp53 = z*L[31];
double Lstmp54 = z*L[44];
double Lstmp55 = z*L[42];
double Lstmp56 = L[4] + L[7];
double Lstmp57 = pow(z, 2);
double Lstmp58 = (1.0/2.0)*Lstmp57;
double Lstmp59 = L[11] + L[15];
double Lstmp60 = pow(z, 3);
double Lstmp61 = (1.0/6.0)*Lstmp60;
double Lstmp62 = L[9] + L[12];
double Lstmp63 = Lstmp58*Lstmp62;
double Lstmp64 = L[18] + L[22];
double Lstmp65 = Lstmp61*Lstmp64;
double Lstmp66 = L[10] + L[14];
double Lstmp67 = Lstmp58*Lstmp66;
double Lstmp68 = L[20] + L[24];
double Lstmp69 = Lstmp61*Lstmp68;
double Lstmp70 = L[17] + L[21];
double Lstmp71 = Lstmp58*Lstmp70;
double Lstmp72 = Lstmp71*y;
double Lstmp73 = L[29] + L[33];
double Lstmp74 = Lstmp61*Lstmp73;
double Lstmp75 = Lstmp74*y;
double Lstmp76 = L[16] + L[19];
double Lstmp77 = (1.0/4.0)*Lstmp57;
double Lstmp78 = Lstmp5*Lstmp77;
double Lstmp79 = L[27] + L[31];
double Lstmp80 = (1.0/12.0)*Lstmp60;
double Lstmp81 = Lstmp5*Lstmp80;
double Lstmp82 = L[25] + L[28];
double Lstmp83 = (1.0/12.0)*Lstmp57;
double Lstmp84 = Lstmp7*Lstmp83;
double Lstmp85 = L[38] + L[42];
double Lstmp86 = (1.0/36.0)*Lstmp60;
double Lstmp87 = L[36] + L[39];
double Lstmp88 = (1.0/48.0)*Lstmp57;
double Lstmp89 = L[19] + L[23];
double Lstmp90 = Lstmp12*Lstmp77;
double Lstmp91 = L[31] + L[35];
double Lstmp92 = Lstmp12*Lstmp80;
double Lstmp93 = L[30] + L[34];
double Lstmp94 = Lstmp14*Lstmp83;
double Lstmp95 = L[44] + L[48];
double Lstmp96 = L[43] + L[47];
double Lstmp97 = L[28] + L[32];
double Lstmp98 = Lstmp90*Lstmp97;
double Lstmp99 = L[42] + L[46];
double Lstmp100 = Lstmp92*Lstmp99;
double Lstmp101 = L[41] + L[45];
double Lstmp102 = Lstmp101*Lstmp94;
double Lstmp103 = L[26] + L[30];
double Lstmp104 = Lstmp103*Lstmp78;
double Lstmp105 = L[40] + L[44];
double Lstmp106 = Lstmp105*Lstmp81;
double Lstmp107 = L[37] + L[41];
double Lstmp108 = Lstmp107*Lstmp84;
double Lstmp109 = L[39] + L[43];
double Lstmp110 = L[16] + 2*L[19] + L[23];
double Lstmp111 = pow(z, 4);
double Lstmp112 = (1.0/24.0)*Lstmp111;
double Lstmp113 = L[27] + 2*L[31] + L[35];
double Lstmp114 = (1.0/120.0)*pow(z, 5);
double Lstmp115 = L[25] + 2*L[28] + L[32];
double Lstmp116 = Lstmp112*Lstmp115;
double Lstmp117 = L[38] + 2*L[42] + L[46];
double Lstmp118 = Lstmp114*Lstmp117;
double Lstmp119 = L[26] + 2*L[30] + L[34];
double Lstmp120 = Lstmp112*Lstmp119;
double Lstmp121 = L[40] + 2*L[44] + L[48];
double Lstmp122 = Lstmp114*Lstmp121;
double Lstmp123 = L[37] + 2*L[41] + L[45];
double Lstmp124 = Lstmp112*Lstmp123;
double Lstmp125 = Lstmp124*y;
double Lstmp126 = L[36] + 2*L[39] + L[43];
double Lstmp127 = L[39] + 2*L[43] + L[47];
double Lstmp128 = L[36] + 3*L[39] + 3*L[43] + L[47];
double Lstmp129 = x*L[19];
double Lstmp130 = x*L[30];
double Lstmp131 = x*L[43];
double Lstmp132 = Lstmp53*x;
double Lstmp133 = Lstmp54*x;
double Lstmp134 = Lstmp58*Lstmp76;
double Lstmp135 = Lstmp61*Lstmp79;
double Lstmp136 = Lstmp103*Lstmp58;
double Lstmp137 = Lstmp136*y;
double Lstmp138 = Lstmp105*Lstmp61;
double Lstmp139 = Lstmp138*y;
double Lstmp140 = Lstmp109*Lstmp90;
double Lstmp141 = Lstmp107*Lstmp78;
double Lstmp142 = Lstmp112*Lstmp126;
double Lstmp143 = y*L[12];
double Lstmp144 = Lstmp35*y;
double Lstmp145 = y*L[19];
double Lstmp146 = y*L[28];
double Lstmp147 = y*L[39];
double Lstmp148 = Lstmp53*y;
double Lstmp149 = Lstmp55*y;
double Lstmp150 = Lstmp58*Lstmp89;
double Lstmp151 = Lstmp61*Lstmp91;
double Lstmp152 = Lstmp58*Lstmp97;
double Lstmp153 = Lstmp152*y;
double Lstmp154 = Lstmp61*Lstmp99;
double Lstmp155 = Lstmp154*y;
double Lstmp156 = Lstmp101*Lstmp90;
double Lstmp157 = Lstmp109*Lstmp78;
double Lstmp158 = Lstmp112*Lstmp127;
double Lstmp159 = y*L[13];
double Lstmp160 = x*L[22];
double Lstmp161 = x*L[33];
double Lstmp162 = x*L[46];
double Lstmp163 = y*L[20];
double Lstmp164 = y*L[29];
double Lstmp165 = y*L[40];
double Lstmp166 = Lstmp62*z;
double Lstmp167 = Lstmp66*z;
double Lstmp168 = Lstmp70*z;
double Lstmp169 = Lstmp168*y;
double Lstmp170 = Lstmp58*Lstmp64;
double Lstmp171 = Lstmp76*z;
double Lstmp172 = Lstmp82*z;
double Lstmp173 = Lstmp87*z;
double Lstmp174 = Lstmp58*Lstmp68;
double Lstmp175 = Lstmp89*z;
double Lstmp176 = Lstmp93*z;
double Lstmp177 = Lstmp96*z;
double Lstmp178 = Lstmp58*Lstmp73;
double Lstmp179 = Lstmp178*y;
double Lstmp180 = Lstmp97*z;
double Lstmp181 = Lstmp180*x;
double Lstmp182 = Lstmp101*z;
double Lstmp183 = Lstmp182*x;
double Lstmp184 = Lstmp103*z;
double Lstmp185 = Lstmp184*y;
double Lstmp186 = Lstmp107*z;
double Lstmp187 = Lstmp186*y;
double Lstmp188 = Lstmp90*Lstmp99;
double Lstmp189 = Lstmp105*Lstmp78;
double Lstmp190 = Lstmp109*z;
double Lstmp191 = Lstmp115*Lstmp61;
double Lstmp192 = Lstmp112*Lstmp117;
double Lstmp193 = Lstmp119*Lstmp61;
double Lstmp194 = Lstmp112*Lstmp121;
double Lstmp195 = Lstmp123*Lstmp61;
double Lstmp196 = Lstmp195*y;
double Lstmp197 = x*L[28];
double Lstmp198 = x*L[41];
double Lstmp199 = Lstmp55*x;
double Lstmp200 = Lstmp58*Lstmp82;
double Lstmp201 = Lstmp61*Lstmp85;
double Lstmp202 = Lstmp107*Lstmp58;
double Lstmp203 = Lstmp202*y;
double Lstmp204 = Lstmp109*Lstmp58;
double Lstmp205 = Lstmp204*y;
double Lstmp206 = x*L[31];
double Lstmp207 = x*L[44];
double Lstmp208 = Lstmp58*Lstmp79;
double Lstmp209 = Lstmp105*Lstmp58;
double Lstmp210 = Lstmp209*y;
double Lstmp211 = Lstmp190*x;
double Lstmp212 = Lstmp126*Lstmp61;
double Lstmp213 = y*L[21];
double Lstmp214 = Lstmp37*y;
double Lstmp215 = y*L[30];
double Lstmp216 = y*L[41];
double Lstmp217 = Lstmp54*y;
double Lstmp218 = Lstmp58*Lstmp93;
double Lstmp219 = Lstmp61*Lstmp95;
double Lstmp220 = Lstmp101*Lstmp58;
double Lstmp221 = Lstmp220*y;
double Lstmp222 = y*L[22];
double Lstmp223 = y*L[31];
double Lstmp224 = y*L[42];
double Lstmp225 = Lstmp180*y;
double Lstmp226 = Lstmp58*Lstmp91;
double Lstmp227 = Lstmp58*Lstmp99;
double Lstmp228 = Lstmp227*y;
double Lstmp229 = Lstmp190*y;
double Lstmp230 = Lstmp127*Lstmp61;
double Lstmp231 = x*L[39];
double Lstmp232 = Lstmp58*Lstmp87;
double Lstmp233 = x*L[42];
double Lstmp234 = Lstmp58*Lstmp85;
double Lstmp235 = y*L[32];
double Lstmp236 = Lstmp39*y;
double Lstmp237 = y*L[43];
double Lstmp238 = Lstmp58*Lstmp96;
double Lstmp239 = y*L[33];
double Lstmp240 = y*L[44];
double Lstmp241 = Lstmp182*y;
double Lstmp242 = Lstmp58*Lstmp95;
double Lstmp243 = y*L[45];
double Lstmp244 = y*L[46];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp10*Lstmp27 + Lstmp10*Lstmp28 + Lstmp10*Lstmp46 + Lstmp10*L[16] - Lstmp100*x - Lstmp102*x - Lstmp104*y - Lstmp106*y - Lstmp108*y - 1.0/8.0*Lstmp109*Lstmp47*Lstmp57 + Lstmp11*Lstmp29 + Lstmp11*Lstmp30 + Lstmp11*L[25] + Lstmp110*Lstmp112 + Lstmp111*Lstmp126*Lstmp50 + Lstmp111*Lstmp127*Lstmp52 + Lstmp113*Lstmp114 + Lstmp116*x + Lstmp118*x + Lstmp120*y + Lstmp122*y + Lstmp125*x - 1.0/720.0*Lstmp128*pow(z, 6) + Lstmp13*Lstmp19 + Lstmp13*Lstmp31 + Lstmp13*Lstmp36 + Lstmp13*L[7] + (1.0/36.0)*Lstmp14*Lstmp7*L[41] - Lstmp14*Lstmp86*Lstmp95 + Lstmp15*Lstmp20 + Lstmp15*Lstmp32 + Lstmp15*Lstmp38 + Lstmp15*L[14] + Lstmp16*Lstmp50*L[43] - Lstmp16*Lstmp88*Lstmp96 + Lstmp17*Lstmp21 + Lstmp17*Lstmp33 + Lstmp17*Lstmp40 + Lstmp17*L[23] + Lstmp18*Lstmp22 + Lstmp18*Lstmp34 + Lstmp18*L[34] + Lstmp2*y + Lstmp23*Lstmp6 + Lstmp24*Lstmp6 + Lstmp25*Lstmp8 + Lstmp26*Lstmp8 + Lstmp4*x + Lstmp42*Lstmp6 + Lstmp44*Lstmp8 + Lstmp48*Lstmp53 + Lstmp48*L[19] + Lstmp49*Lstmp54 + Lstmp49*L[30] + Lstmp51*Lstmp55 + Lstmp51*L[28] + Lstmp52*Lstmp9*L[39] - Lstmp56*Lstmp58 - Lstmp59*Lstmp61 + Lstmp6*L[4] - Lstmp63*x - Lstmp65*x - Lstmp67*y - Lstmp69*y - Lstmp7*Lstmp85*Lstmp86 - Lstmp72*x - Lstmp75*x - Lstmp76*Lstmp78 - Lstmp79*Lstmp81 + Lstmp8*L[9] - Lstmp82*Lstmp84 - Lstmp87*Lstmp88*Lstmp9 - Lstmp89*Lstmp90 - Lstmp91*Lstmp92 - Lstmp93*Lstmp94 - Lstmp98*x + (1.0/720.0)*pow(x, 6)*L[36] + x*L[1] + (1.0/720.0)*pow(y, 6)*L[47] + y*L[2] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + Lstmp10*Lstmp29 + Lstmp10*Lstmp30 + Lstmp10*L[25] - Lstmp100 - Lstmp102 + Lstmp11*L[36] + Lstmp116 + Lstmp118 + Lstmp125 + Lstmp129*Lstmp13 + Lstmp13*Lstmp132 + Lstmp13*Lstmp35 + Lstmp13*L[12] + Lstmp130*Lstmp15 + Lstmp131*Lstmp17 + Lstmp133*Lstmp15 - Lstmp134*x - Lstmp135*x - Lstmp137*x - Lstmp139*x - Lstmp140*x - Lstmp141*y + Lstmp142*x + Lstmp15*Lstmp37 + Lstmp15*L[21] + Lstmp17*Lstmp39 + Lstmp17*L[32] + Lstmp18*L[45] + Lstmp23*x + Lstmp24*x + Lstmp25*Lstmp6 + Lstmp26*Lstmp6 + Lstmp27*Lstmp8 + Lstmp28*Lstmp8 + Lstmp4 + Lstmp42*x + Lstmp44*Lstmp6 + Lstmp46*Lstmp8 + Lstmp48*Lstmp55 + Lstmp48*L[28] + Lstmp49*L[41] + Lstmp51*L[39] + Lstmp6*L[9] - Lstmp63 - Lstmp65 - Lstmp72 - Lstmp75 - Lstmp78*Lstmp82 + Lstmp8*L[16] - Lstmp81*Lstmp85 - Lstmp84*Lstmp87 - Lstmp98 + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp10*Lstmp147 + Lstmp10*Lstmp45 + Lstmp10*L[26] - Lstmp104 - Lstmp106 - Lstmp108 + Lstmp11*L[37] + Lstmp120 + Lstmp122 + Lstmp124*x + Lstmp13*Lstmp20 + Lstmp13*Lstmp32 + Lstmp13*Lstmp38 + Lstmp13*L[14] + Lstmp143*x + Lstmp144*x + Lstmp145*Lstmp6 + Lstmp146*Lstmp8 + Lstmp148*Lstmp6 + Lstmp149*Lstmp8 + Lstmp15*Lstmp21 + Lstmp15*Lstmp33 + Lstmp15*Lstmp40 + Lstmp15*L[23] - Lstmp150*y - Lstmp151*y - Lstmp153*x - Lstmp155*x - Lstmp156*x - Lstmp157*y + Lstmp158*y + Lstmp17*Lstmp22 + Lstmp17*Lstmp34 + Lstmp17*L[34] + Lstmp18*L[47] + Lstmp2 + Lstmp3*x + Lstmp31*y + Lstmp41*Lstmp6 + Lstmp43*Lstmp8 + Lstmp48*Lstmp54 + Lstmp48*L[30] + Lstmp49*L[43] + Lstmp51*L[41] + Lstmp6*L[10] - Lstmp67 - Lstmp69 - Lstmp71*x - Lstmp74*x + Lstmp8*L[17] - Lstmp90*Lstmp93 - Lstmp92*Lstmp95 - Lstmp94*Lstmp96 + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += Lstmp10*Lstmp165 - Lstmp10*Lstmp173 + Lstmp10*L[27] + Lstmp11*L[38] + Lstmp110*Lstmp61 + Lstmp112*Lstmp113 - Lstmp114*Lstmp128 + Lstmp126*Lstmp81 + Lstmp127*Lstmp92 + Lstmp13*Lstmp160 - Lstmp13*Lstmp175 - Lstmp13*Lstmp181 + Lstmp13*L[15] + Lstmp15*Lstmp161 - Lstmp15*Lstmp176 - Lstmp15*Lstmp183 + Lstmp15*L[24] + Lstmp159*x + Lstmp162*Lstmp17 + Lstmp163*Lstmp6 + Lstmp164*Lstmp8 - Lstmp166*x - Lstmp167*y - Lstmp169*x - Lstmp17*Lstmp177 + Lstmp17*L[35] - Lstmp170*x - Lstmp171*Lstmp6 - Lstmp172*Lstmp8 - Lstmp174*y - Lstmp179*x + Lstmp18*L[48] - Lstmp185*Lstmp6 - Lstmp187*Lstmp8 - Lstmp188*x - Lstmp189*y - Lstmp190*Lstmp48 + Lstmp191*x + Lstmp192*x + Lstmp193*y + Lstmp194*y + Lstmp196*x + Lstmp48*L[31] + Lstmp49*L[44] + Lstmp51*L[42] - Lstmp56*z - Lstmp58*Lstmp59 + Lstmp6*L[11] - Lstmp78*Lstmp79 + Lstmp8*L[18] - Lstmp84*Lstmp85 - Lstmp90*Lstmp91 - Lstmp94*Lstmp95 + x*L[6] + y*L[8] + L[3];
#pragma omp atomic
Ls[4] += Lstmp10*L[36] + Lstmp13*Lstmp197 + Lstmp13*Lstmp199 + Lstmp13*Lstmp53 + Lstmp13*L[19] - Lstmp134 - Lstmp135 - Lstmp137 - Lstmp139 - Lstmp140 + Lstmp142 + Lstmp15*Lstmp198 + Lstmp15*Lstmp54 + Lstmp15*L[30] + Lstmp17*L[43] - Lstmp200*x - Lstmp201*x - Lstmp203*x + Lstmp23 + Lstmp24 + Lstmp25*x + Lstmp26*x + Lstmp27*Lstmp6 + Lstmp28*Lstmp6 + Lstmp29*Lstmp8 + Lstmp30*Lstmp8 + Lstmp42 + Lstmp44*x + Lstmp46*Lstmp6 + Lstmp48*L[39] + Lstmp6*L[16] - Lstmp78*Lstmp87 + Lstmp8*L[25] + x*L[9] + L[4];
#pragma omp atomic
Ls[5] += Lstmp10*L[37] + Lstmp124 + Lstmp13*Lstmp130 + Lstmp13*Lstmp133 + Lstmp13*Lstmp37 + Lstmp13*L[21] + Lstmp131*Lstmp15 - Lstmp136*x - Lstmp138*x - Lstmp141 + Lstmp143 + Lstmp144 + Lstmp145*x + Lstmp146*Lstmp6 + Lstmp147*Lstmp8 + Lstmp148*x + Lstmp149*Lstmp6 + Lstmp15*Lstmp39 + Lstmp15*L[32] - Lstmp153 - Lstmp155 - Lstmp156 + Lstmp17*L[45] - Lstmp205*x + Lstmp3 + Lstmp41*x + Lstmp43*Lstmp6 + Lstmp45*Lstmp8 + Lstmp48*L[41] + Lstmp6*L[17] - Lstmp71 - Lstmp74 + Lstmp8*L[26] + x*L[10] + L[5];
#pragma omp atomic
Ls[6] += Lstmp10*L[38] - Lstmp13*Lstmp180 + Lstmp13*Lstmp206 - Lstmp13*Lstmp211 + Lstmp13*L[22] - Lstmp15*Lstmp182 + Lstmp15*Lstmp207 + Lstmp15*L[33] + Lstmp159 + Lstmp163*x + Lstmp164*Lstmp6 + Lstmp165*Lstmp8 - Lstmp166 - Lstmp169 + Lstmp17*L[46] - Lstmp170 - Lstmp171*x - Lstmp172*Lstmp6 - Lstmp173*Lstmp8 - Lstmp179 - Lstmp185*x - Lstmp187*Lstmp6 - Lstmp188 + Lstmp191 + Lstmp192 + Lstmp196 - Lstmp208*x - Lstmp210*x + Lstmp212*x + Lstmp48*L[42] + Lstmp6*L[18] - Lstmp78*Lstmp85 + Lstmp8*L[27] + x*L[11] + L[6];
#pragma omp atomic
Ls[7] += Lstmp10*L[39] + Lstmp13*Lstmp21 + Lstmp13*Lstmp33 + Lstmp13*Lstmp40 + Lstmp13*L[23] + Lstmp15*Lstmp22 + Lstmp15*Lstmp34 + Lstmp15*L[34] - Lstmp150 - Lstmp151 - Lstmp152*x - Lstmp154*x - Lstmp157 + Lstmp158 + Lstmp17*L[47] + Lstmp19 + Lstmp213*x + Lstmp214*x + Lstmp215*Lstmp6 + Lstmp216*Lstmp8 + Lstmp217*Lstmp6 - Lstmp218*y - Lstmp219*y - Lstmp221*x + Lstmp31 + Lstmp32*y + Lstmp36 + Lstmp48*L[43] + Lstmp53*Lstmp6 + Lstmp55*Lstmp8 + Lstmp6*L[19] + Lstmp8*L[28] - Lstmp90*Lstmp96 + y*L[14] + L[7];
#pragma omp atomic
Ls[8] += Lstmp10*L[40] + Lstmp13*Lstmp161 - Lstmp13*Lstmp176 - Lstmp13*Lstmp183 + Lstmp13*L[24] + Lstmp15*Lstmp162 - Lstmp15*Lstmp177 + Lstmp15*L[35] - Lstmp167 - Lstmp168*x + Lstmp17*L[48] - Lstmp174 - Lstmp175*y - Lstmp178*x - Lstmp184*Lstmp6 - Lstmp186*Lstmp8 - Lstmp189 + Lstmp193 + Lstmp194 + Lstmp195*x + Lstmp222*x + Lstmp223*Lstmp6 + Lstmp224*Lstmp8 - Lstmp225*x - Lstmp226*y - Lstmp228*x - Lstmp229*Lstmp6 + Lstmp230*y + Lstmp48*L[44] + Lstmp6*L[20] + Lstmp8*L[29] - Lstmp90*Lstmp95 + x*L[13] + y*L[15] + L[8];
#pragma omp atomic
Ls[9] += Lstmp13*Lstmp231 + Lstmp13*Lstmp55 + Lstmp13*L[28] + Lstmp15*L[41] - Lstmp200 - Lstmp201 - Lstmp203 - Lstmp232*x + Lstmp25 + Lstmp26 + Lstmp27*x + Lstmp28*x + Lstmp29*Lstmp6 + Lstmp30*Lstmp6 + Lstmp44 + Lstmp46*x + Lstmp6*L[25] + Lstmp8*L[36] + x*L[16] + L[9];
#pragma omp atomic
Ls[10] += Lstmp13*Lstmp198 + Lstmp13*Lstmp54 + Lstmp13*L[30] - Lstmp136 - Lstmp138 + Lstmp145 + Lstmp146*x + Lstmp147*Lstmp6 + Lstmp148 + Lstmp149*x + Lstmp15*L[43] - Lstmp202*x - Lstmp205 + Lstmp41 + Lstmp43*x + Lstmp45*Lstmp6 + Lstmp6*L[26] + Lstmp8*L[37] + x*L[17] + L[10];
#pragma omp atomic
Ls[11] += -Lstmp13*Lstmp190 + Lstmp13*Lstmp233 + Lstmp13*L[31] + Lstmp15*L[44] + Lstmp163 + Lstmp164*x + Lstmp165*Lstmp6 - Lstmp171 - Lstmp172*x - Lstmp173*Lstmp6 - Lstmp185 - Lstmp187*x - Lstmp208 - Lstmp210 + Lstmp212 - Lstmp234*x + Lstmp6*L[27] + Lstmp8*L[38] + x*L[18] + L[11];
#pragma omp atomic
Ls[12] += Lstmp129 + Lstmp13*Lstmp131 + Lstmp13*Lstmp39 + Lstmp13*L[32] + Lstmp132 + Lstmp15*L[45] - Lstmp152 - Lstmp154 - Lstmp204*x + Lstmp213 + Lstmp214 + Lstmp215*x + Lstmp216*Lstmp6 + Lstmp217*x - Lstmp221 + Lstmp35 + Lstmp55*Lstmp6 + Lstmp6*L[28] + Lstmp8*L[39] + L[12];
#pragma omp atomic
Ls[13] += -Lstmp13*Lstmp182 + Lstmp13*Lstmp207 + Lstmp13*L[33] + Lstmp15*L[46] - Lstmp168 - Lstmp178 - Lstmp184*x - Lstmp186*Lstmp6 + Lstmp195 - Lstmp209*x + Lstmp222 + Lstmp223*x + Lstmp224*Lstmp6 - Lstmp225 - Lstmp228 - Lstmp229*x + Lstmp6*L[29] + Lstmp8*L[40] + x*L[20] + L[13];
#pragma omp atomic
Ls[14] += Lstmp13*Lstmp22 + Lstmp13*Lstmp34 + Lstmp13*L[34] + Lstmp15*L[47] + Lstmp20 - Lstmp218 - Lstmp219 - Lstmp220*x + Lstmp235*x + Lstmp236*x + Lstmp237*Lstmp6 - Lstmp238*y + Lstmp32 + Lstmp33*y + Lstmp38 + Lstmp54*Lstmp6 + Lstmp6*L[30] + Lstmp8*L[41] + y*L[23] + L[14];
#pragma omp atomic
Ls[15] += Lstmp13*Lstmp162 - Lstmp13*Lstmp177 + Lstmp13*L[35] + Lstmp15*L[48] + Lstmp160 - Lstmp175 - Lstmp176*y - Lstmp181 - Lstmp190*Lstmp6 - Lstmp226 - Lstmp227*x + Lstmp230 + Lstmp239*x + Lstmp240*Lstmp6 - Lstmp241*x - Lstmp242*y + Lstmp6*L[31] + Lstmp8*L[42] + y*L[24] + L[15];
#pragma omp atomic
Ls[16] += Lstmp13*L[39] - Lstmp232 + Lstmp27 + Lstmp28 + Lstmp29*x + Lstmp30*x + Lstmp46 + Lstmp6*L[36] + x*L[25] + L[16];
#pragma omp atomic
Ls[17] += Lstmp13*L[41] + Lstmp146 + Lstmp147*x + Lstmp149 - Lstmp202 + Lstmp43 + Lstmp45*x + Lstmp6*L[37] + x*L[26] + L[17];
#pragma omp atomic
Ls[18] += Lstmp13*L[42] + Lstmp164 + Lstmp165*x - Lstmp172 - Lstmp173*x - Lstmp187 - Lstmp234 + Lstmp6*L[38] + x*L[27] + L[18];
#pragma omp atomic
Ls[19] += Lstmp13*L[43] + Lstmp197 + Lstmp199 - Lstmp204 + Lstmp215 + Lstmp216*x + Lstmp217 + Lstmp53 + Lstmp6*L[39] + L[19];
#pragma omp atomic
Ls[20] += Lstmp13*L[44] - Lstmp184 - Lstmp186*x - Lstmp209 + Lstmp223 + Lstmp224*x - Lstmp229 + Lstmp6*L[40] + x*L[29] + L[20];
#pragma omp atomic
Ls[21] += Lstmp13*L[45] + Lstmp130 + Lstmp133 - Lstmp220 + Lstmp235 + Lstmp236 + Lstmp237*x + Lstmp37 + Lstmp6*L[41] + L[21];
#pragma omp atomic
Ls[22] += Lstmp13*L[46] - Lstmp180 + Lstmp206 - Lstmp211 - Lstmp227 + Lstmp239 + Lstmp240*x - Lstmp241 + Lstmp6*L[42] + L[22];
#pragma omp atomic
Ls[23] += Lstmp13*L[47] + Lstmp21 - Lstmp238 + Lstmp243*x + Lstmp33 + Lstmp34*y + Lstmp40 + Lstmp6*L[43] + y*L[34] + L[23];
#pragma omp atomic
Ls[24] += Lstmp13*L[48] + Lstmp161 - Lstmp176 - Lstmp177*y - Lstmp183 - Lstmp242 + Lstmp244*x + Lstmp6*L[44] + y*L[35] + L[24];
#pragma omp atomic
Ls[25] += Lstmp29 + Lstmp30 + x*L[36] + L[25];
#pragma omp atomic
Ls[26] += Lstmp147 + Lstmp45 + x*L[37] + L[26];
#pragma omp atomic
Ls[27] += Lstmp165 - Lstmp173 + x*L[38] + L[27];
#pragma omp atomic
Ls[28] += Lstmp216 + Lstmp231 + Lstmp55 + L[28];
#pragma omp atomic
Ls[29] += -Lstmp186 + Lstmp224 + x*L[40] + L[29];
#pragma omp atomic
Ls[30] += Lstmp198 + Lstmp237 + Lstmp54 + L[30];
#pragma omp atomic
Ls[31] += -Lstmp190 + Lstmp233 + Lstmp240 + L[31];
#pragma omp atomic
Ls[32] += Lstmp131 + Lstmp243 + Lstmp39 + L[32];
#pragma omp atomic
Ls[33] += -Lstmp182 + Lstmp207 + Lstmp244 + L[33];
#pragma omp atomic
Ls[34] += Lstmp22 + Lstmp34 + y*L[47] + L[34];
#pragma omp atomic
Ls[35] += Lstmp162 - Lstmp177 + y*L[48] + L[35];
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

}

void L2Pc_7(double x, double y, double z, double * L, double * F) {
double Ftmp0 = x*y;
double Ftmp1 = x*z;
double Ftmp2 = y*z;
double Ftmp3 = Ftmp0*z;
double Ftmp4 = pow(x, 2);
double Ftmp5 = (1.0/2.0)*Ftmp4;
double Ftmp6 = pow(x, 3);
double Ftmp7 = (1.0/6.0)*Ftmp6;
double Ftmp8 = pow(x, 4);
double Ftmp9 = (1.0/24.0)*Ftmp8;
double Ftmp10 = (1.0/120.0)*pow(x, 5);
double Ftmp11 = pow(y, 2);
double Ftmp12 = (1.0/2.0)*Ftmp11;
double Ftmp13 = pow(y, 3);
double Ftmp14 = (1.0/6.0)*Ftmp13;
double Ftmp15 = pow(y, 4);
double Ftmp16 = (1.0/24.0)*Ftmp15;
double Ftmp17 = (1.0/120.0)*pow(y, 5);
double Ftmp18 = Ftmp12*x;
double Ftmp19 = Ftmp14*x;
double Ftmp20 = Ftmp16*x;
double Ftmp21 = Ftmp5*y;
double Ftmp22 = Ftmp5*z;
double Ftmp23 = Ftmp7*y;
double Ftmp24 = Ftmp7*z;
double Ftmp25 = Ftmp9*y;
double Ftmp26 = Ftmp9*z;
double Ftmp27 = Ftmp12*z;
double Ftmp28 = Ftmp14*z;
double Ftmp29 = Ftmp16*z;
double Ftmp30 = Ftmp1*Ftmp12;
double Ftmp31 = Ftmp1*Ftmp14;
double Ftmp32 = Ftmp2*Ftmp5;
double Ftmp33 = Ftmp2*Ftmp7;
double Ftmp34 = (1.0/4.0)*Ftmp11*Ftmp4;
double Ftmp35 = (1.0/12.0)*Ftmp4;
double Ftmp36 = Ftmp13*Ftmp35;
double Ftmp37 = (1.0/12.0)*Ftmp11*Ftmp6;
double Ftmp38 = Ftmp34*z;
double Ftmp39 = pow(z, 2);
double Ftmp40 = L[9] + L[12];
double Ftmp41 = pow(z, 3);
double Ftmp42 = L[18] + L[22];
double Ftmp43 = L[16] + L[19];
double Ftmp44 = L[27] + L[31];
double Ftmp45 = L[17] + L[21];
double Ftmp46 = L[29] + L[33];
double Ftmp47 = L[26] + L[30];
double Ftmp48 = L[40] + L[44];
double Ftmp49 = L[25] + L[28];
double Ftmp50 = L[38] + L[42];
double Ftmp51 = L[36] + L[39];
double Ftmp52 = L[28] + L[32];
double Ftmp53 = L[42] + L[46];
double Ftmp54 = L[41] + L[45];
double Ftmp55 = L[39] + L[43];
double Ftmp56 = L[37] + L[41];
double Ftmp57 = L[25] + 2*L[28] + L[32];
double Ftmp58 = (1.0/24.0)*pow(z, 4);
double Ftmp59 = L[38] + 2*L[42] + L[46];
double Ftmp60 = pow(z, 5);
double Ftmp61 = (1.0/120.0)*Ftmp60;
double Ftmp62 = L[36] + 2*L[39] + L[43];
double Ftmp63 = Ftmp58*x;
double Ftmp64 = L[37] + 2*L[41] + L[45];
double Ftmp65 = Ftmp58*y;
double Ftmp66 = L[10] + L[14];
double Ftmp67 = L[20] + L[24];
double Ftmp68 = L[19] + L[23];
double Ftmp69 = L[31] + L[35];
double Ftmp70 = L[30] + L[34];
double Ftmp71 = L[44] + L[48];
double Ftmp72 = L[43] + L[47];
double Ftmp73 = L[26] + 2*L[30] + L[34];
double Ftmp74 = L[40] + 2*L[44] + L[48];
double Ftmp75 = L[39] + 2*L[43] + L[47];
double Ftmp76 = (1.0/6.0)*Ftmp41;
#pragma omp atomic
F[0] += -Ftmp0*L[10] - Ftmp1*L[11] - Ftmp10*L[36] + (1.0/4.0)*Ftmp11*Ftmp39*Ftmp52 + (1.0/4.0)*Ftmp11*Ftmp39*Ftmp55*x + (1.0/12.0)*Ftmp11*Ftmp41*Ftmp53 - Ftmp12*L[12] + (1.0/12.0)*Ftmp13*Ftmp39*Ftmp54 - Ftmp14*L[21] - Ftmp16*L[32] - Ftmp17*L[45] - Ftmp18*L[19] - Ftmp19*L[30] - Ftmp2*L[13] - Ftmp20*L[43] - Ftmp21*L[17] - Ftmp22*L[18] - Ftmp23*L[26] - Ftmp24*L[27] - Ftmp25*L[37] - Ftmp26*L[38] - Ftmp27*L[22] - Ftmp28*L[33] - Ftmp29*L[46] - Ftmp3*L[20] - Ftmp30*L[31] - Ftmp31*L[44] - Ftmp32*L[29] - Ftmp33*L[40] - Ftmp34*L[28] - Ftmp36*L[41] - Ftmp37*L[39] - Ftmp38*L[42] + (1.0/4.0)*Ftmp39*Ftmp4*Ftmp49 + (1.0/4.0)*Ftmp39*Ftmp4*Ftmp56*y + (1.0/2.0)*Ftmp39*Ftmp40 + (1.0/2.0)*Ftmp39*Ftmp43*x + (1.0/2.0)*Ftmp39*Ftmp45*y + (1.0/2.0)*Ftmp39*Ftmp47*x*y + (1.0/12.0)*Ftmp39*Ftmp51*Ftmp6 + (1.0/12.0)*Ftmp4*Ftmp41*Ftmp50 + (1.0/6.0)*Ftmp41*Ftmp42 + (1.0/6.0)*Ftmp41*Ftmp44*x + (1.0/6.0)*Ftmp41*Ftmp46*y + (1.0/6.0)*Ftmp41*Ftmp48*x*y - Ftmp5*L[9] - Ftmp57*Ftmp58 - Ftmp59*Ftmp61 - Ftmp62*Ftmp63 - Ftmp64*Ftmp65 - Ftmp7*L[16] - Ftmp9*L[25] - x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -Ftmp0*L[12] - Ftmp1*L[13] - Ftmp10*L[37] + (1.0/4.0)*Ftmp11*Ftmp39*Ftmp54*x + (1.0/4.0)*Ftmp11*Ftmp39*Ftmp70 + (1.0/12.0)*Ftmp11*Ftmp41*Ftmp71 - Ftmp12*L[14] + (1.0/12.0)*Ftmp13*Ftmp39*Ftmp72 - Ftmp14*L[23] - Ftmp16*L[34] - Ftmp17*L[47] - Ftmp18*L[21] - Ftmp19*L[32] - Ftmp2*L[15] - Ftmp20*L[45] - Ftmp21*L[19] - Ftmp22*L[20] - Ftmp23*L[28] - Ftmp24*L[29] - Ftmp25*L[39] - Ftmp26*L[40] - Ftmp27*L[24] - Ftmp28*L[35] - Ftmp29*L[48] - Ftmp3*L[22] - Ftmp30*L[33] - Ftmp31*L[46] - Ftmp32*L[31] - Ftmp33*L[42] - Ftmp34*L[30] - Ftmp36*L[43] - Ftmp37*L[41] - Ftmp38*L[44] + (1.0/4.0)*Ftmp39*Ftmp4*Ftmp47 + (1.0/4.0)*Ftmp39*Ftmp4*Ftmp55*y + (1.0/2.0)*Ftmp39*Ftmp45*x + (1.0/2.0)*Ftmp39*Ftmp52*x*y + (1.0/12.0)*Ftmp39*Ftmp56*Ftmp6 + (1.0/2.0)*Ftmp39*Ftmp66 + (1.0/2.0)*Ftmp39*Ftmp68*y + (1.0/12.0)*Ftmp4*Ftmp41*Ftmp48 + (1.0/6.0)*Ftmp41*Ftmp46*x + (1.0/6.0)*Ftmp41*Ftmp53*x*y + (1.0/6.0)*Ftmp41*Ftmp67 + (1.0/6.0)*Ftmp41*Ftmp69*y - Ftmp5*L[10] - Ftmp58*Ftmp73 - Ftmp61*Ftmp74 - Ftmp63*Ftmp64 - Ftmp65*Ftmp75 - Ftmp7*L[17] - Ftmp9*L[26] - x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -Ftmp0*Ftmp64*Ftmp76 - Ftmp0*L[13] - Ftmp10*L[38] + (1.0/4.0)*Ftmp11*Ftmp39*Ftmp53*x + (1.0/4.0)*Ftmp11*Ftmp39*Ftmp69 + (1.0/4.0)*Ftmp11*Ftmp4*Ftmp55*z - 1.0/12.0*Ftmp11*Ftmp41*Ftmp75 + (1.0/2.0)*Ftmp11*Ftmp52*x*z + (1.0/2.0)*Ftmp11*Ftmp68*z - Ftmp12*L[15] + (1.0/12.0)*Ftmp13*Ftmp39*Ftmp71 + (1.0/6.0)*Ftmp13*Ftmp54*x*z + (1.0/6.0)*Ftmp13*Ftmp70*z - Ftmp14*L[24] + (1.0/24.0)*Ftmp15*Ftmp72*z - Ftmp16*L[35] - Ftmp17*L[48] - Ftmp18*L[22] - Ftmp19*L[33] - Ftmp20*L[46] - Ftmp21*L[20] - Ftmp23*L[29] - Ftmp25*L[40] - Ftmp34*L[31] - Ftmp35*Ftmp41*Ftmp62 - Ftmp36*L[44] - Ftmp37*L[42] + (1.0/4.0)*Ftmp39*Ftmp4*Ftmp44 + (1.0/4.0)*Ftmp39*Ftmp4*Ftmp48*y + (1.0/2.0)*Ftmp39*Ftmp42*x + (1.0/2.0)*Ftmp39*Ftmp46*x*y + (1.0/12.0)*Ftmp39*Ftmp50*Ftmp6 + (1.0/2.0)*Ftmp39*Ftmp67*y + (1.0/2.0)*Ftmp39*(L[11] + L[15]) + (1.0/2.0)*Ftmp4*Ftmp43*z + (1.0/2.0)*Ftmp4*Ftmp47*y*z + Ftmp40*x*z + Ftmp45*x*y*z + (1.0/6.0)*Ftmp49*Ftmp6*z - Ftmp5*L[11] + (1.0/24.0)*Ftmp51*Ftmp8*z + (1.0/6.0)*Ftmp56*Ftmp6*y*z - Ftmp57*Ftmp76*x - Ftmp58*(L[27] + 2*L[31] + L[35]) - Ftmp59*Ftmp63 + (1.0/120.0)*Ftmp60*(L[36] + 3*L[39] + 3*L[43] + L[47]) - Ftmp65*Ftmp74 + Ftmp66*y*z - Ftmp7*L[18] - Ftmp73*Ftmp76*y - Ftmp76*(L[16] + 2*L[19] + L[23]) - Ftmp9*L[27] - x*L[6] - y*L[8] + z*(L[4] + L[7]) - L[3];

}

void M2Pc_7(double x, double y, double z, double * M, double * F) {
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double Ftmp0 = pow(Rinv, 3);
double Ftmp1 = pow(Rinv, 2);
double Ftmp2 = 3*Ftmp1;
double Ftmp3 = Ftmp2*z;
double Ftmp4 = Ftmp2*x;
double Ftmp5 = Ftmp4*y;
double Ftmp6 = Ftmp3*M[2];
double Ftmp7 = pow(Rinv, 4);
double Ftmp8 = pow(x, 2);
double Ftmp9 = Ftmp1*Ftmp8;
double Ftmp10 = 3*Ftmp9;
double Ftmp11 = pow(y, 2);
double Ftmp12 = pow(Rinv, 6);
double Ftmp13 = 30*x;
double Ftmp14 = Ftmp12*Ftmp13;
double Ftmp15 = pow(y, 3);
double Ftmp16 = Ftmp15*M[13];
double Ftmp17 = Ftmp11*z;
double Ftmp18 = Ftmp14*Ftmp17;
double Ftmp19 = Ftmp12*Ftmp8;
double Ftmp20 = 105*Ftmp19;
double Ftmp21 = pow(Rinv, 8);
double Ftmp22 = 30*Ftmp11;
double Ftmp23 = Ftmp19*Ftmp22;
double Ftmp24 = pow(Rinv, 10);
double Ftmp25 = 1890*Ftmp24;
double Ftmp26 = z*M[32];
double Ftmp27 = Ftmp15*Ftmp8;
double Ftmp28 = 5*Ftmp9;
double Ftmp29 = Ftmp28 - 3;
double Ftmp30 = Ftmp1*Ftmp11;
double Ftmp31 = 5*Ftmp30;
double Ftmp32 = Ftmp31 - 1;
double Ftmp33 = Ftmp9 - 1;
double Ftmp34 = Ftmp10 - 1;
double Ftmp35 = 3*Ftmp30;
double Ftmp36 = Ftmp35 - 1;
double Ftmp37 = 15*Ftmp7;
double Ftmp38 = Ftmp37*y;
double Ftmp39 = 7*Ftmp9;
double Ftmp40 = Ftmp39 - 3;
double Ftmp41 = Ftmp40*M[16];
double Ftmp42 = 7*Ftmp30;
double Ftmp43 = Ftmp42 - 3;
double Ftmp44 = Ftmp43*M[20];
double Ftmp45 = Ftmp37*z;
double Ftmp46 = Ftmp40*M[17];
double Ftmp47 = Ftmp42 - 1;
double Ftmp48 = Ftmp47*M[21];
double Ftmp49 = Ftmp13*Ftmp33;
double Ftmp50 = Ftmp7*Ftmp8;
double Ftmp51 = 30*M[8];
double Ftmp52 = Ftmp28 - 1;
double Ftmp53 = Ftmp52*M[9];
double Ftmp54 = Ftmp38*x;
double Ftmp55 = Ftmp31 - 3;
double Ftmp56 = Ftmp55*M[13];
double Ftmp57 = Ftmp45*x;
double Ftmp58 = Ftmp52*M[10];
double Ftmp59 = Ftmp32*M[14];
double Ftmp60 = Ftmp29*M[8];
double Ftmp61 = 15*Ftmp8;
double Ftmp62 = Ftmp61*Ftmp7;
double Ftmp63 = Ftmp32*M[11];
double Ftmp64 = Ftmp39 - 1;
double Ftmp65 = Ftmp21*x;
double Ftmp66 = 9*Ftmp30;
double Ftmp67 = 420*M[33];
double Ftmp68 = Ftmp15*Ftmp67*(Ftmp66 - 5);
double Ftmp69 = Ftmp21*Ftmp8;
double Ftmp70 = 1890*Ftmp69;
double Ftmp71 = 1260*M[34];
double Ftmp72 = Ftmp36*Ftmp71;
double Ftmp73 = Ftmp11*Ftmp21;
double Ftmp74 = x*z;
double Ftmp75 = Ftmp73*Ftmp74;
double Ftmp76 = 2835*Ftmp69;
double Ftmp77 = Ftmp76*y;
double Ftmp78 = Ftmp34*M[28];
double Ftmp79 = Ftmp78*z;
double Ftmp80 = Ftmp26*Ftmp36;
double Ftmp81 = 11*Ftmp30;
double Ftmp82 = Ftmp81 - 5;
double Ftmp83 = 1260*M[31];
double Ftmp84 = Ftmp36*Ftmp83;
double Ftmp85 = Ftmp81 - 3;
double Ftmp86 = pow(Rinv, 12);
double Ftmp87 = 41580*Ftmp86*(13*Ftmp30 - 5)*M[60];
double Ftmp88 = pow(x, 4);
double Ftmp89 = Ftmp7*Ftmp88;
double Ftmp90 = 63*Ftmp89;
double Ftmp91 = -70*Ftmp9 + Ftmp90 + 15;
double Ftmp92 = pow(y, 4);
double Ftmp93 = Ftmp7*Ftmp92;
double Ftmp94 = 21*Ftmp93;
double Ftmp95 = 14*Ftmp30;
double Ftmp96 = -Ftmp95;
double Ftmp97 = Ftmp96 + 1;
double Ftmp98 = Ftmp94 + Ftmp97;
double Ftmp99 = -10*Ftmp9;
double Ftmp100 = Ftmp99 + 3;
double Ftmp101 = 30*Ftmp9;
double Ftmp102 = -Ftmp101;
double Ftmp103 = Ftmp102 + 35*Ftmp89 + 3;
double Ftmp104 = 30*Ftmp30;
double Ftmp105 = -Ftmp104;
double Ftmp106 = Ftmp105 + 35*Ftmp93 + 3;
double Ftmp107 = 315*Ftmp12;
double Ftmp108 = Ftmp107*y;
double Ftmp109 = 33*Ftmp89;
double Ftmp110 = Ftmp102 + Ftmp109 + 5;
double Ftmp111 = Ftmp110*M[36];
double Ftmp112 = 33*Ftmp93;
double Ftmp113 = Ftmp105 + Ftmp112 + 5;
double Ftmp114 = Ftmp113*M[44];
double Ftmp115 = Ftmp107*z;
double Ftmp116 = Ftmp110*M[37];
double Ftmp117 = Ftmp12*x;
double Ftmp118 = Ftmp117*y;
double Ftmp119 = 3780*Ftmp89 - 5040*Ftmp9 + 1260;
double Ftmp120 = 21*Ftmp89;
double Ftmp121 = 14*Ftmp9;
double Ftmp122 = -Ftmp121;
double Ftmp123 = Ftmp122 + 1;
double Ftmp124 = Ftmp120 + Ftmp123;
double Ftmp125 = Ftmp124*x;
double Ftmp126 = Ftmp108*M[25];
double Ftmp127 = 63*Ftmp93;
double Ftmp128 = Ftmp127 - 70*Ftmp30 + 15;
double Ftmp129 = Ftmp128*M[33];
double Ftmp130 = 105*Ftmp118;
double Ftmp131 = Ftmp117*z;
double Ftmp132 = Ftmp98*M[34];
double Ftmp133 = Ftmp115*x;
double Ftmp134 = 143*Ftmp89;
double Ftmp135 = Ftmp134 - 110*Ftmp9 + 15;
double Ftmp136 = 420*M[24];
double Ftmp137 = Ftmp91*M[24];
double Ftmp138 = Ftmp98*M[31];
double Ftmp139 = 11*Ftmp89;
double Ftmp140 = Ftmp122 + 3;
double Ftmp141 = 18*Ftmp9;
double Ftmp142 = -Ftmp141;
double Ftmp143 = Ftmp109 + Ftmp142 + 1;
double Ftmp144 = -16*Ftmp9;
double Ftmp145 = Ftmp139 + Ftmp144 + 5;
double Ftmp146 = Ftmp15*Ftmp24;
double Ftmp147 = 143*Ftmp93;
double Ftmp148 = 5670*M[61];
double Ftmp149 = Ftmp148*(Ftmp147 - 154*Ftmp30 + 35);
double Ftmp150 = Ftmp147 - 110*Ftmp30 + 15;
double Ftmp151 = 5670*Ftmp150*M[62];
double Ftmp152 = Ftmp11*Ftmp24;
double Ftmp153 = Ftmp152*Ftmp74;
double Ftmp154 = 41580*M[52];
double Ftmp155 = Ftmp24*z;
double Ftmp156 = Ftmp155*y;
double Ftmp157 = Ftmp156*Ftmp8;
double Ftmp158 = Ftmp135*M[52];
double Ftmp159 = 10395*Ftmp157;
double Ftmp160 = Ftmp152*Ftmp8;
double Ftmp161 = 5670*M[59];
double Ftmp162 = Ftmp150*Ftmp161;
double Ftmp163 = pow(x, 6);
double Ftmp164 = 429*Ftmp12;
double Ftmp165 = Ftmp163*Ftmp164;
double Ftmp166 = Ftmp165 - 693*Ftmp89 + 315*Ftmp9 - 35;
double Ftmp167 = pow(y, 6);
double Ftmp168 = Ftmp164*Ftmp167;
double Ftmp169 = Ftmp168 + 135*Ftmp30 - 495*Ftmp93 - 5;
double Ftmp170 = Ftmp12*Ftmp163;
double Ftmp171 = 231*Ftmp170 - 315*Ftmp89 + 105*Ftmp9 - 5;
double Ftmp172 = Ftmp12*Ftmp167;
double Ftmp173 = 231*Ftmp172 + 105*Ftmp30 - 315*Ftmp93 - 5;
double Ftmp174 = 125*Ftmp9;
double Ftmp175 = 143*Ftmp170;
double Ftmp176 = 5670*Ftmp174 + 5670*Ftmp175 - 1434510*Ftmp89 - 85050;
double Ftmp177 = Ftmp65*y;
double Ftmp178 = Ftmp165 - 495*Ftmp89 + 135*Ftmp9 - 5;
double Ftmp179 = Ftmp178*M[49];
double Ftmp180 = 2835*Ftmp177;
double Ftmp181 = Ftmp168 + 315*Ftmp30 - 693*Ftmp93 - 35;
double Ftmp182 = Ftmp181*M[61];
double Ftmp183 = Ftmp21*z;
double Ftmp184 = Ftmp183*x;
double Ftmp185 = Ftmp178*M[50];
double Ftmp186 = 2835*Ftmp184;
double Ftmp187 = Ftmp169*M[62];
double Ftmp188 = -Ftmp39;
double Ftmp189 = Ftmp11*Ftmp50;
double Ftmp190 = 63*Ftmp189;
double Ftmp191 = Ftmp190 + 3;
double Ftmp192 = Ftmp188 + Ftmp191 - 21*Ftmp30;
double Ftmp193 = 5670*M[48];
double Ftmp194 = Ftmp166*M[48];
double Ftmp195 = Ftmp169*M[59];
double Ftmp196 = 8*Ftmp30;
double Ftmp197 = 14*Ftmp189;
double Ftmp198 = -Ftmp9;
double Ftmp199 = Ftmp198 + 1;
double Ftmp200 = 35*Ftmp189 - Ftmp28 - Ftmp31 + 1;
double Ftmp201 = y*M[40];
double Ftmp202 = -Ftmp35;
double Ftmp203 = -Ftmp10;
double Ftmp204 = Ftmp203 + 1;
double Ftmp205 = 11*Ftmp189 + Ftmp202 + Ftmp204;
double Ftmp206 = 945*Ftmp12*Ftmp205;
double Ftmp207 = 33*Ftmp189;
double Ftmp208 = Ftmp204 + Ftmp207 - Ftmp66;
double Ftmp209 = Ftmp208*M[41];
double Ftmp210 = 18*Ftmp189;
double Ftmp211 = -10*Ftmp30;
double Ftmp212 = Ftmp211 + 3;
double Ftmp213 = 210*M[29];
double Ftmp214 = -Ftmp42;
double Ftmp215 = Ftmp191 + Ftmp214 - 21*Ftmp9;
double Ftmp216 = Ftmp215*M[29];
double Ftmp217 = Ftmp199 + Ftmp210;
double Ftmp218 = 210*M[30];
double Ftmp219 = Ftmp188 + Ftmp190 + Ftmp214 + 1;
double Ftmp220 = 105*Ftmp219*M[30];
double Ftmp221 = 143*Ftmp189 - 33*Ftmp30 - 33*Ftmp9 + 9;
double Ftmp222 = -12*Ftmp30;
double Ftmp223 = Ftmp192*M[27];
double Ftmp224 = 22*Ftmp189;
double Ftmp225 = Ftmp203 + 3;
double Ftmp226 = Ftmp224 + Ftmp225;
double Ftmp227 = 9*Ftmp9;
double Ftmp228 = Ftmp202 + Ftmp207 - Ftmp227 + 1;
double Ftmp229 = -16*Ftmp30;
double Ftmp230 = 26*Ftmp189;
double Ftmp231 = 20790*M[56];
double Ftmp232 = Ftmp221*M[56];
double Ftmp233 = 4*Ftmp30;
double Ftmp234 = Ftmp233*Ftmp36 + Ftmp98;
double Ftmp235 = Ftmp234*M[45];
double Ftmp236 = 99*Ftmp93;
double Ftmp237 = Ftmp233*Ftmp82 + Ftmp236 - 90*Ftmp30 + 15;
double Ftmp238 = Ftmp237*M[60];
double Ftmp239 = -66*Ftmp189;
double Ftmp240 = Ftmp8*Ftmp92;
double Ftmp241 = 143*Ftmp12;
double Ftmp242 = Ftmp240*Ftmp241;
double Ftmp243 = -Ftmp112;
double Ftmp244 = Ftmp243 + 18*Ftmp30;
double Ftmp245 = Ftmp239 + Ftmp242 + Ftmp244 + Ftmp34;
double Ftmp246 = -330*Ftmp189 - 5;
double Ftmp247 = -Ftmp109;
double Ftmp248 = Ftmp11*Ftmp88;
double Ftmp249 = Ftmp164*Ftmp248;
double Ftmp250 = Ftmp247 + Ftmp249;
double Ftmp251 = Ftmp101 + Ftmp246 + Ftmp250 + 45*Ftmp30;
double Ftmp252 = -36*Ftmp189;
double Ftmp253 = 99*Ftmp12;
double Ftmp254 = Ftmp240*Ftmp253;
double Ftmp255 = -126*Ftmp189;
double Ftmp256 = 231*Ftmp12;
double Ftmp257 = Ftmp240*Ftmp256 + Ftmp255 + Ftmp64 - Ftmp94 + Ftmp95;
double Ftmp258 = 8*Ftmp9;
double Ftmp259 = Ftmp248*Ftmp253;
double Ftmp260 = -102*Ftmp189 - 2;
double Ftmp261 = -Ftmp120 + Ftmp121 + Ftmp248*Ftmp256 + Ftmp255 + Ftmp47;
double Ftmp262 = -44*Ftmp189;
double Ftmp263 = Ftmp164*Ftmp240;
double Ftmp264 = Ftmp104 + Ftmp243 + Ftmp246 + Ftmp263 + 45*Ftmp9;
double Ftmp265 = Ftmp264*M[57];
double Ftmp266 = -220*Ftmp189 - 15;
double Ftmp267 = 1890*Ftmp177;
double Ftmp268 = Ftmp241*Ftmp248;
double Ftmp269 = Ftmp141 + Ftmp239 + Ftmp247 + Ftmp268 + Ftmp36;
double Ftmp270 = 8505*Ftmp269*M[53];
double Ftmp271 = Ftmp270*y;
double Ftmp272 = -418*Ftmp189;
double Ftmp273 = Ftmp272 + 69*Ftmp30;
double Ftmp274 = Ftmp249 + 84*Ftmp9;
double Ftmp275 = -198*Ftmp189 - 1;
double Ftmp276 = Ftmp227 + Ftmp244 + Ftmp263 + Ftmp275;
double Ftmp277 = Ftmp276*M[58];
double Ftmp278 = Ftmp141 + Ftmp250 + Ftmp275 + Ftmp66;
double Ftmp279 = Ftmp278*M[54];
double Ftmp280 = Ftmp249 - 22*Ftmp89;
double Ftmp281 = Ftmp245*M[55];
double Ftmp282 = -132*Ftmp189 - 3;
double Ftmp283 = Ftmp263 + 84*Ftmp30;
double Ftmp284 = Ftmp251*M[51];
double Ftmp285 = 125*Ftmp30;
double Ftmp286 = -506*Ftmp189 - 10;
double Ftmp287 = Ftmp12*y;
double Ftmp288 = pow(x, 3);
double Ftmp289 = Ftmp288*Ftmp51;
double Ftmp290 = 105*Ftmp131;
double Ftmp291 = 30*y*z;
double Ftmp292 = Ftmp19*Ftmp291;
double Ftmp293 = Ftmp17*Ftmp288;
double Ftmp294 = Ftmp30 - 1;
double Ftmp295 = Ftmp37*x;
double Ftmp296 = Ftmp64*M[19];
double Ftmp297 = Ftmp43*M[23];
double Ftmp298 = Ftmp294*Ftmp7;
double Ftmp299 = Ftmp38*z;
double Ftmp300 = 15*Ftmp11;
double Ftmp301 = Ftmp300*Ftmp7;
double Ftmp302 = Ftmp288*y;
double Ftmp303 = Ftmp136*(Ftmp227 - 5);
double Ftmp304 = 1890*Ftmp73;
double Ftmp305 = 2835*Ftmp73;
double Ftmp306 = Ftmp305*x;
double Ftmp307 = 1260*Ftmp34;
double Ftmp308 = Ftmp307*Ftmp69*y*z;
double Ftmp309 = 11*Ftmp9;
double Ftmp310 = Ftmp309 - 5;
double Ftmp311 = Ftmp11*M[25];
double Ftmp312 = Ftmp307*Ftmp69;
double Ftmp313 = Ftmp309 - 3;
double Ftmp314 = Ftmp154*Ftmp86*(13*Ftmp9 - 5);
double Ftmp315 = Ftmp107*x;
double Ftmp316 = Ftmp143*M[39];
double Ftmp317 = Ftmp113*M[47];
double Ftmp318 = -Ftmp233 + 3*Ftmp93 + 1;
double Ftmp319 = Ftmp108*z;
double Ftmp320 = Ftmp287*z;
double Ftmp321 = Ftmp107*Ftmp124;
double Ftmp322 = Ftmp11*Ftmp12;
double Ftmp323 = 105*Ftmp322;
double Ftmp324 = Ftmp229 + 11*Ftmp93 + 5;
double Ftmp325 = Ftmp193*(Ftmp134 - 154*Ftmp9 + 35);
double Ftmp326 = 10395*Ftmp153;
double Ftmp327 = 5670*Ftmp135;
double Ftmp328 = Ftmp327*Ftmp8*M[50];
double Ftmp329 = Ftmp327*M[49];
double Ftmp330 = 143*Ftmp172;
double Ftmp331 = Ftmp285 + Ftmp330 - 253*Ftmp93 - 15;
double Ftmp332 = Ftmp183*y;
double Ftmp333 = 2835*Ftmp332;
double Ftmp334 = 5670*Ftmp332;
double Ftmp335 = -Ftmp30;
double Ftmp336 = Ftmp335 + 1;
double Ftmp337 = Ftmp228*M[43];
double Ftmp338 = Ftmp210 + Ftmp336;
double Ftmp339 = -12*Ftmp9;
double Ftmp340 = Ftmp202 + Ftmp224;
double Ftmp341 = 8505*Ftmp281;
double Ftmp342 = Ftmp272 + 69*Ftmp9;
double Ftmp343 = Ftmp263 - 22*Ftmp93;
double Ftmp344 = -Ftmp294*Ftmp30;
double Ftmp345 = 2*Ftmp30;
double Ftmp346 = pow(z, 2);
double Ftmp347 = Ftmp12*z;
double Ftmp348 = 30*Ftmp347;
double Ftmp349 = Ftmp346*M[10];
double Ftmp350 = Ftmp12*Ftmp346;
double Ftmp351 = Ftmp25*Ftmp346;
double Ftmp352 = Ftmp351*x;
double Ftmp353 = 2835*Ftmp21*Ftmp346;
double Ftmp354 = Ftmp353*x*y;
double Ftmp355 = Ftmp346*M[26];
double Ftmp356 = Ftmp11 + Ftmp8;
double Ftmp357 = 105*Ftmp320;
double Ftmp358 = -Ftmp11*Ftmp141;
double Ftmp359 = 3*Ftmp11;
double Ftmp360 = Ftmp359 + Ftmp8;
double Ftmp361 = 3*Ftmp8;
double Ftmp362 = Ftmp11 + Ftmp361;
double Ftmp363 = Ftmp24*Ftmp346;
double Ftmp364 = 10395*Ftmp363*x*y;
double Ftmp365 = Ftmp11*Ftmp9;
double Ftmp366 = -22*Ftmp365;
double Ftmp367 = Ftmp359 + Ftmp361;
double Ftmp368 = Ftmp25*z;
double Ftmp369 = Ftmp368*x;
double Ftmp370 = 6*Ftmp1;
double Ftmp371 = -36*Ftmp365;
double Ftmp372 = Ftmp11*Ftmp89;
double Ftmp373 = -44*Ftmp365;
double Ftmp374 = 22*Ftmp1;
double Ftmp375 = -Ftmp374*Ftmp92;
double Ftmp376 = 6*Ftmp11 + Ftmp375;
double Ftmp377 = 5670*Ftmp155;
double Ftmp378 = -Ftmp374*Ftmp88;
double Ftmp379 = Ftmp378 + 6*Ftmp8;
double Ftmp380 = -220*Ftmp365;
double Ftmp381 = 429*Ftmp372;
double Ftmp382 = 429*Ftmp8*Ftmp93;
double Ftmp383 = -132*Ftmp365;
#pragma omp atomic
F[0] += Ftmp0*(3*Ftmp1*Ftmp29*M[8] + 3*Ftmp1*Ftmp32*M[11] + 6*Ftmp1*Ftmp33*x*M[3] + 3*Ftmp1*Ftmp34*x*M[3] + 3*Ftmp1*Ftmp36*x*M[6] - Ftmp10*M[0] + 15*Ftmp103*Ftmp7*x*M[15] + 15*Ftmp106*Ftmp7*x*M[22] - Ftmp107*Ftmp138*Ftmp8 - Ftmp108*Ftmp111 - Ftmp108*Ftmp114 + 1890*Ftmp11*Ftmp113*Ftmp21*x*M[46] + 60*Ftmp11*Ftmp12*Ftmp43*x*M[22] + 210*Ftmp11*Ftmp21*Ftmp8*z*M[21] + 3780*Ftmp11*Ftmp24*Ftmp8*Ftmp85*z*M[45] - Ftmp11*Ftmp69*Ftmp84 + 6*Ftmp11*Ftmp7*x*M[6] + 2835*Ftmp110*Ftmp21*Ftmp8*y*M[36] + 2835*Ftmp110*Ftmp21*Ftmp8*z*M[37] + 2835*Ftmp113*Ftmp21*Ftmp8*y*M[44] + 2835*Ftmp113*Ftmp21*x*y*z*M[47] - Ftmp115*Ftmp116 - Ftmp115*Ftmp125*M[26] - Ftmp115*Ftmp209 - Ftmp115*Ftmp235 - Ftmp118*Ftmp119*M[25] - Ftmp118*Ftmp213*(Ftmp203 + Ftmp210 + Ftmp212) - Ftmp119*Ftmp131*M[26] + 315*Ftmp12*Ftmp166*M[48] + 315*Ftmp12*Ftmp169*M[59] + 315*Ftmp12*Ftmp171*x*M[35] + 315*Ftmp12*Ftmp173*x*M[46] + 945*Ftmp12*Ftmp245*M[55] + 315*Ftmp12*Ftmp251*M[51] + 315*Ftmp12*Ftmp257*x*M[42] + 315*Ftmp12*Ftmp261*x*M[38] + 210*Ftmp12*Ftmp33*Ftmp8*y*M[16] + 210*Ftmp12*Ftmp33*Ftmp8*z*M[17] + 210*Ftmp12*Ftmp33*x*y*z*M[19] + 315*Ftmp12*Ftmp34*y*z*M[28] + 315*Ftmp12*Ftmp36*y*z*M[32] + 105*Ftmp12*Ftmp40*Ftmp8*y*M[16] + 105*Ftmp12*Ftmp40*Ftmp8*z*M[17] + 105*Ftmp12*Ftmp43*Ftmp8*y*M[20] + 105*Ftmp12*Ftmp43*x*y*z*M[23] + 105*Ftmp12*Ftmp47*Ftmp8*z*M[21] + 105*Ftmp12*Ftmp64*x*y*z*M[19] + 1890*Ftmp12*x*(33*Ftmp170 + 35*Ftmp9 - Ftmp90 - 5)*M[35] + 630*Ftmp12*x*(Ftmp252 + Ftmp254 + 20*Ftmp30 + Ftmp33 - 39*Ftmp93)*M[42] + 630*Ftmp12*x*(Ftmp258 + Ftmp259 + Ftmp260 + 19*Ftmp30 - 6*Ftmp89)*M[38] - Ftmp125*Ftmp126 - Ftmp129*Ftmp130 - Ftmp130*Ftmp216 - Ftmp131*Ftmp218*(Ftmp211 + Ftmp217) - Ftmp131*Ftmp220 - Ftmp132*Ftmp133 + 945*Ftmp135*Ftmp21*y*z*M[52] - Ftmp136*Ftmp19*(Ftmp122 + 9*Ftmp89 + 5) - Ftmp137*Ftmp20 - Ftmp14*Ftmp16 + 2835*Ftmp143*Ftmp21*x*y*z*M[39] + 3780*Ftmp145*Ftmp21*Ftmp8*y*M[36] + 3780*Ftmp145*Ftmp21*Ftmp8*z*M[37] - Ftmp146*Ftmp149*x + 210*Ftmp15*Ftmp21*Ftmp8*M[20] + 210*Ftmp15*Ftmp21*x*z*M[23] + 3780*Ftmp15*Ftmp24*Ftmp8*Ftmp82*M[44] + 3780*Ftmp15*Ftmp24*Ftmp82*x*z*M[47] - Ftmp151*Ftmp153 - Ftmp154*Ftmp157*(Ftmp142 + 13*Ftmp89 + 5) - Ftmp157*Ftmp231*(Ftmp225 + Ftmp229 + Ftmp230) - Ftmp158*Ftmp159 - Ftmp159*Ftmp232 - Ftmp159*Ftmp238 - Ftmp160*Ftmp162 - Ftmp176*Ftmp177*M[49] - Ftmp176*Ftmp184*M[50] - Ftmp179*Ftmp180 - Ftmp18*M[14] - Ftmp180*Ftmp182 - Ftmp180*Ftmp265 - 1890*Ftmp184*(Ftmp273 + Ftmp280 + 28*Ftmp9 - 6)*M[54] - 5670*Ftmp184*(Ftmp242 + Ftmp262 + 24*Ftmp30 + Ftmp33 - 55*Ftmp93)*M[58] - Ftmp185*Ftmp186 - Ftmp186*Ftmp187 - Ftmp186*Ftmp277 - Ftmp186*Ftmp279 - 210*Ftmp19*(Ftmp217 + Ftmp222)*M[27] + 15*Ftmp192*Ftmp7*M[27] - Ftmp193*Ftmp69*(Ftmp175 - 297*Ftmp89 + 189*Ftmp9 - 35) - Ftmp194*Ftmp76 - Ftmp195*Ftmp76 - Ftmp2*y*M[4] - Ftmp20*Ftmp223 - Ftmp20*y*z*M[12] + 15*Ftmp200*Ftmp7*x*M[18] - Ftmp201*Ftmp206 + 8505*Ftmp205*Ftmp21*Ftmp8*y*M[40] + 2835*Ftmp208*Ftmp21*Ftmp8*z*M[41] + 945*Ftmp21*Ftmp221*y*z*M[56] + 2835*Ftmp21*Ftmp228*x*y*z*M[43] + 2835*Ftmp21*Ftmp234*Ftmp8*z*M[45] + 945*Ftmp21*Ftmp237*y*z*M[60] + 1890*Ftmp21*Ftmp8*y*(Ftmp226 + Ftmp96)*M[40] + 1890*Ftmp21*Ftmp8*z*(Ftmp198 + Ftmp224 + Ftmp97)*M[41] + 3780*Ftmp21*x*y*z*(Ftmp139 + Ftmp140)*M[39] + 1890*Ftmp21*x*y*z*(Ftmp222 + Ftmp226)*M[43] - Ftmp23*M[11] - Ftmp25*Ftmp26*Ftmp27 - Ftmp267*(Ftmp273 + Ftmp274 - 66*Ftmp89 - 18)*M[53] - Ftmp267*(Ftmp263 + Ftmp266 + 120*Ftmp30 + 15*Ftmp9 - 165*Ftmp93)*M[57] - Ftmp27*Ftmp87*z - Ftmp271*Ftmp65 - 8505*Ftmp281*Ftmp69 - Ftmp284*Ftmp76 - Ftmp3*M[5] - Ftmp33*Ftmp50*Ftmp51 - Ftmp33*Ftmp70*y*z*M[28] - Ftmp38*Ftmp41 - Ftmp38*Ftmp44 - Ftmp45*Ftmp46 - Ftmp45*Ftmp48 - Ftmp49*Ftmp7*y*M[9] - Ftmp49*Ftmp7*z*M[10] - Ftmp5*M[1] - Ftmp53*Ftmp54 - Ftmp54*Ftmp56 - Ftmp57*Ftmp58 - Ftmp57*Ftmp59 - Ftmp6*x - Ftmp60*Ftmp62 - Ftmp62*Ftmp63 - Ftmp65*Ftmp68 + 15*Ftmp7*Ftmp8*y*M[4] + 15*Ftmp7*Ftmp8*z*M[5] + 15*Ftmp7*Ftmp91*M[24] + 45*Ftmp7*Ftmp98*M[31] + 15*Ftmp7*x*y*z*M[7] + 60*Ftmp7*x*(Ftmp100 + 7*Ftmp89)*M[15] + 30*Ftmp7*x*(-Ftmp196 + Ftmp197 + Ftmp199)*M[18] + 15*Ftmp7*y*z*M[12] - Ftmp70*(Ftmp10 + Ftmp282 + Ftmp283 - 209*Ftmp93)*M[55] - Ftmp70*(Ftmp280 + Ftmp285 + Ftmp286 + 32*Ftmp9)*M[51] - Ftmp72*Ftmp75 - Ftmp77*Ftmp79 - Ftmp77*Ftmp80 + M[0]);
#pragma omp atomic
F[1] += Ftmp0*(6*Ftmp1*Ftmp294*y*M[6] + 3*Ftmp1*Ftmp34*y*M[3] + 3*Ftmp1*Ftmp36*y*M[6] + 3*Ftmp1*Ftmp52*M[9] + 3*Ftmp1*Ftmp55*M[13] + 15*Ftmp103*Ftmp7*y*M[15] + 15*Ftmp106*Ftmp7*y*M[22] - Ftmp108*Ftmp138*x + 2835*Ftmp11*Ftmp110*Ftmp21*x*M[36] + 2835*Ftmp11*Ftmp113*Ftmp21*x*M[44] + 2835*Ftmp11*Ftmp113*Ftmp21*z*M[47] + 210*Ftmp11*Ftmp12*Ftmp294*x*M[20] + 210*Ftmp11*Ftmp12*Ftmp294*z*M[23] + 105*Ftmp11*Ftmp12*Ftmp40*x*M[16] + 105*Ftmp11*Ftmp12*Ftmp43*x*M[20] + 105*Ftmp11*Ftmp12*Ftmp43*z*M[23] + 105*Ftmp11*Ftmp12*Ftmp64*z*M[19] + 2835*Ftmp11*Ftmp143*Ftmp21*z*M[39] + 8505*Ftmp11*Ftmp205*Ftmp21*x*M[40] + 2835*Ftmp11*Ftmp21*Ftmp228*z*M[43] + 210*Ftmp11*Ftmp21*Ftmp288*M[16] + 3780*Ftmp11*Ftmp21*Ftmp324*x*M[44] + 3780*Ftmp11*Ftmp21*Ftmp324*z*M[47] + 210*Ftmp11*Ftmp21*Ftmp8*z*M[19] + 1890*Ftmp11*Ftmp21*x*(Ftmp140 + Ftmp340)*M[40] + 1890*Ftmp11*Ftmp21*z*(Ftmp123 + Ftmp224 + Ftmp335)*M[43] + 3780*Ftmp11*Ftmp24*Ftmp288*Ftmp310*M[36] + 3780*Ftmp11*Ftmp24*Ftmp313*Ftmp8*z*M[39] - Ftmp11*Ftmp290*M[12] + 15*Ftmp11*Ftmp7*x*M[4] + 15*Ftmp11*Ftmp7*z*M[7] + 1890*Ftmp110*Ftmp21*Ftmp8*y*M[35] + 2835*Ftmp110*Ftmp21*x*y*z*M[37] - Ftmp111*Ftmp315 - Ftmp114*Ftmp315 - Ftmp115*Ftmp316 - Ftmp115*Ftmp317 - Ftmp115*Ftmp337 - Ftmp118*Ftmp318*Ftmp83 - 210*Ftmp118*(Ftmp100 + Ftmp202 + Ftmp210)*M[27] + 315*Ftmp12*Ftmp171*y*M[35] + 315*Ftmp12*Ftmp173*y*M[46] + 315*Ftmp12*Ftmp178*M[49] + 315*Ftmp12*Ftmp181*M[61] + 315*Ftmp12*Ftmp257*y*M[42] + 315*Ftmp12*Ftmp261*y*M[38] + 315*Ftmp12*Ftmp264*M[57] + 945*Ftmp12*Ftmp269*M[53] + 210*Ftmp12*Ftmp294*x*y*z*M[21] + 315*Ftmp12*Ftmp34*x*z*M[28] + 315*Ftmp12*Ftmp36*x*z*M[32] + 60*Ftmp12*Ftmp40*Ftmp8*y*M[15] + 105*Ftmp12*Ftmp40*x*y*z*M[17] + 105*Ftmp12*Ftmp47*x*y*z*M[21] + 1890*Ftmp12*y*(-Ftmp127 + 33*Ftmp172 + 35*Ftmp30 - 5)*M[46] + 630*Ftmp12*y*(Ftmp196 + Ftmp254 + Ftmp260 + 19*Ftmp9 - 6*Ftmp93)*M[42] + 630*Ftmp12*y*(Ftmp252 + Ftmp259 + Ftmp294 - 39*Ftmp89 + 20*Ftmp9)*M[38] - Ftmp124*Ftmp319*M[26] + 45*Ftmp124*Ftmp7*M[25] + 15*Ftmp128*Ftmp7*M[33] - Ftmp129*Ftmp323 - Ftmp13*Ftmp298*y*M[11] - Ftmp130*Ftmp137 - Ftmp130*Ftmp223 - Ftmp132*Ftmp319 + 945*Ftmp135*Ftmp21*x*z*M[52] - Ftmp148*Ftmp73*(189*Ftmp30 + Ftmp330 - 297*Ftmp93 - 35) - 3780*Ftmp152*x*z*(Ftmp236 - 166*Ftmp30 - 22*Ftmp344 + Ftmp345*Ftmp82 + 55)*M[60] - Ftmp153*Ftmp231*(Ftmp144 + Ftmp202 + Ftmp230 + 3) - Ftmp156*Ftmp328 - Ftmp158*Ftmp326 - Ftmp160*Ftmp329 - Ftmp161*Ftmp177*Ftmp331 - Ftmp177*Ftmp341 - Ftmp179*Ftmp305 - Ftmp180*Ftmp194 - Ftmp180*Ftmp195 - Ftmp180*Ftmp284 - Ftmp182*Ftmp305 - Ftmp185*Ftmp333 - Ftmp187*Ftmp333 + 15*Ftmp200*Ftmp7*y*M[18] - Ftmp206*x*M[40] + 2835*Ftmp208*Ftmp21*x*y*z*M[41] + 945*Ftmp21*Ftmp221*x*z*M[56] + 2835*Ftmp21*Ftmp234*x*y*z*M[45] + 945*Ftmp21*Ftmp237*x*z*M[60] + 210*Ftmp21*Ftmp288*y*z*M[17] - Ftmp21*Ftmp302*Ftmp303 + 1890*Ftmp21*x*y*z*(Ftmp339 + Ftmp340 + 3)*M[41] + 1260*Ftmp21*x*y*z*(-34*Ftmp30 - 6*Ftmp344 + Ftmp345*Ftmp36 + Ftmp94 + 9)*M[45] - Ftmp213*Ftmp322*(Ftmp338 + Ftmp339) + 15*Ftmp215*Ftmp7*M[29] - Ftmp216*Ftmp323 - Ftmp218*Ftmp320*(Ftmp338 + Ftmp99) - Ftmp22*Ftmp298*M[13] - Ftmp220*Ftmp320 - Ftmp23*M[9] - Ftmp232*Ftmp326 - Ftmp238*Ftmp326 + 3780*Ftmp24*Ftmp288*Ftmp310*y*z*M[37] - Ftmp24*Ftmp302*Ftmp325 - Ftmp25*Ftmp293*M[28] - Ftmp26*Ftmp294*Ftmp304*x - Ftmp265*Ftmp305 - Ftmp267*(Ftmp283 + Ftmp342 - 66*Ftmp93 - 18)*M[55] - Ftmp267*(Ftmp249 + Ftmp266 + 15*Ftmp30 - 165*Ftmp89 + 120*Ftmp9)*M[51] - Ftmp270*Ftmp73 - Ftmp277*Ftmp333 - Ftmp279*Ftmp333 - Ftmp287*Ftmp289 - Ftmp291*Ftmp298*M[14] - Ftmp292*M[10] - Ftmp293*Ftmp314 - Ftmp295*Ftmp41 - Ftmp295*Ftmp44 - Ftmp296*Ftmp45 - Ftmp297*Ftmp45 - Ftmp299*Ftmp58 - Ftmp299*Ftmp59 - Ftmp3*M[7] - Ftmp301*Ftmp53 - Ftmp301*Ftmp56 - Ftmp304*(Ftmp174 + Ftmp286 + 32*Ftmp30 + Ftmp343)*M[57] - Ftmp304*(Ftmp274 + Ftmp282 + Ftmp35 - 209*Ftmp89)*M[53] - Ftmp306*Ftmp79 - Ftmp306*Ftmp80 - Ftmp308*M[26] - Ftmp311*Ftmp312 - Ftmp311*Ftmp321 - Ftmp318*Ftmp320*Ftmp71 - Ftmp322*Ftmp67*(9*Ftmp93 + Ftmp96 + 5) - Ftmp331*Ftmp334*M[62] - 1890*Ftmp332*(28*Ftmp30 + Ftmp342 + Ftmp343 - 6)*M[58] - Ftmp334*(Ftmp262 + Ftmp268 + Ftmp294 - 55*Ftmp89 + 24*Ftmp9)*M[54] - Ftmp35*M[1] - Ftmp4*M[4] - Ftmp5*M[0] - Ftmp54*Ftmp60 - Ftmp54*Ftmp63 - Ftmp6*y + 6*Ftmp7*Ftmp8*y*M[3] + 15*Ftmp7*x*y*z*M[5] + 15*Ftmp7*x*z*M[12] + 60*Ftmp7*y*(Ftmp212 + 7*Ftmp93)*M[22] + 30*Ftmp7*y*(Ftmp197 - Ftmp258 + Ftmp336)*M[18] + M[1]);
#pragma omp atomic
F[2] += Ftmp0*(3*Ftmp1*Ftmp32*M[14] + 3*Ftmp1*Ftmp34*z*M[3] + 3*Ftmp1*Ftmp36*z*M[6] + 3*Ftmp1*Ftmp52*M[10] + 15*Ftmp103*Ftmp7*z*M[15] + 15*Ftmp106*Ftmp7*z*M[22] - Ftmp107*Ftmp132*Ftmp346 - Ftmp108*Ftmp316 - Ftmp108*Ftmp317 - Ftmp108*Ftmp337 + 1890*Ftmp11*Ftmp113*Ftmp21*z*M[46] + 60*Ftmp11*Ftmp12*Ftmp43*z*M[22] + 210*Ftmp11*Ftmp21*Ftmp346*x*M[21] + 3780*Ftmp11*Ftmp24*Ftmp346*Ftmp85*x*M[45] + 6*Ftmp11*Ftmp7*z*M[6] + 2835*Ftmp110*Ftmp21*Ftmp346*x*M[37] + 1890*Ftmp110*Ftmp21*Ftmp8*z*M[35] + 2835*Ftmp110*Ftmp21*x*y*z*M[36] + 2835*Ftmp113*Ftmp21*Ftmp346*y*M[47] + 2835*Ftmp113*Ftmp21*x*y*z*M[44] - Ftmp116*Ftmp315 + 315*Ftmp12*Ftmp169*M[62] + 315*Ftmp12*Ftmp171*z*M[35] + 315*Ftmp12*Ftmp173*z*M[46] + 315*Ftmp12*Ftmp178*M[50] + 315*Ftmp12*Ftmp257*z*M[42] + 315*Ftmp12*Ftmp261*z*M[38] + 315*Ftmp12*Ftmp276*M[58] + 315*Ftmp12*Ftmp278*M[54] + 315*Ftmp12*Ftmp34*x*y*M[28] + 105*Ftmp12*Ftmp346*Ftmp40*x*M[17] + 105*Ftmp12*Ftmp346*Ftmp43*y*M[23] + 105*Ftmp12*Ftmp346*Ftmp47*x*M[21] + 105*Ftmp12*Ftmp346*Ftmp64*y*M[19] + 315*Ftmp12*Ftmp36*x*y*M[32] + 60*Ftmp12*Ftmp40*Ftmp8*z*M[15] + 105*Ftmp12*Ftmp40*x*y*z*M[16] + 105*Ftmp12*Ftmp43*x*y*z*M[20] - Ftmp124*Ftmp126*z + 45*Ftmp124*Ftmp7*M[26] - Ftmp129*Ftmp357 - Ftmp130*Ftmp346*M[12] - Ftmp133*Ftmp138 + 945*Ftmp135*Ftmp21*x*y*M[52] - Ftmp137*Ftmp290 + 2835*Ftmp143*Ftmp21*Ftmp346*y*M[39] - Ftmp146*Ftmp149*z + 210*Ftmp15*Ftmp21*Ftmp346*M[23] + 210*Ftmp15*Ftmp21*x*z*M[20] + 3780*Ftmp15*Ftmp24*Ftmp346*Ftmp82*M[47] + 3780*Ftmp15*Ftmp24*Ftmp82*x*z*M[44] - Ftmp15*Ftmp346*Ftmp87*x - Ftmp15*Ftmp352*M[32] - Ftmp151*Ftmp152*Ftmp346 - Ftmp153*Ftmp162 - Ftmp155*Ftmp288*Ftmp325 - Ftmp157*Ftmp329 - Ftmp158*Ftmp364 - Ftmp16*Ftmp348 - Ftmp179*Ftmp333 - Ftmp18*M[11] - Ftmp182*Ftmp333 - Ftmp183*Ftmp271 - Ftmp183*Ftmp288*Ftmp303 - Ftmp183*Ftmp68 - Ftmp184*Ftmp341 - Ftmp185*Ftmp353 - Ftmp186*Ftmp194 - Ftmp186*Ftmp195 - Ftmp186*Ftmp284 - Ftmp187*Ftmp353 - 30*Ftmp19*Ftmp349 - Ftmp2*Ftmp346*M[2] - Ftmp2*y*M[7] + 15*Ftmp200*Ftmp7*z*M[18] - Ftmp201*Ftmp369*(Ftmp366 + Ftmp367) + 8505*Ftmp205*Ftmp21*x*y*z*M[40] + 2835*Ftmp208*Ftmp21*Ftmp346*x*M[41] - Ftmp209*Ftmp315 + 945*Ftmp21*Ftmp221*x*y*M[56] + 2835*Ftmp21*Ftmp228*Ftmp346*y*M[43] + 2835*Ftmp21*Ftmp234*Ftmp346*x*M[45] + 945*Ftmp21*Ftmp237*x*y*M[60] + 210*Ftmp21*Ftmp288*Ftmp346*M[17] + 210*Ftmp21*Ftmp288*y*z*M[16] + 210*Ftmp21*Ftmp346*Ftmp8*y*M[19] + 210*Ftmp21*Ftmp346*(Ftmp356 + Ftmp358)*M[30] + 210*Ftmp21*x*z*(Ftmp358 + Ftmp360)*M[27] + 210*Ftmp21*y*z*(Ftmp358 + Ftmp362)*M[29] + 630*Ftmp21*z*(Ftmp11 - Ftmp370*Ftmp88 + Ftmp371 + 99*Ftmp372 + 2*Ftmp8)*M[38] + 630*Ftmp21*z*(2*Ftmp11 + Ftmp236*Ftmp8 - Ftmp370*Ftmp92 + Ftmp371 + Ftmp8)*M[42] - Ftmp216*Ftmp357 + 15*Ftmp219*Ftmp7*M[30] - Ftmp22*Ftmp350*M[14] - Ftmp220*Ftmp350 - Ftmp223*Ftmp290 - Ftmp232*Ftmp364 - Ftmp235*Ftmp315 - Ftmp238*Ftmp364 + 3780*Ftmp24*Ftmp288*Ftmp310*Ftmp346*M[37] + 3780*Ftmp24*Ftmp288*Ftmp310*y*z*M[36] + 3780*Ftmp24*Ftmp313*Ftmp346*Ftmp8*y*M[39] - Ftmp265*Ftmp333 - Ftmp277*Ftmp353 - Ftmp279*Ftmp353 - Ftmp289*Ftmp347 - Ftmp292*M[9] - Ftmp295*Ftmp46 - Ftmp295*Ftmp48 - Ftmp296*Ftmp38 - Ftmp297*Ftmp38 - Ftmp299*Ftmp53 - Ftmp299*Ftmp56 - Ftmp3*x*M[0] - Ftmp3*y*M[1] - Ftmp302*Ftmp314*Ftmp346 - Ftmp302*Ftmp351*M[28] - Ftmp308*M[25] - Ftmp312*Ftmp355 - Ftmp321*Ftmp355 - Ftmp328*Ftmp363 - Ftmp346*Ftmp37*Ftmp59 + 15*Ftmp346*Ftmp7*x*M[5] + 15*Ftmp346*Ftmp7*y*M[7] - Ftmp346*Ftmp72*Ftmp73 + 20790*Ftmp346*Ftmp86*x*y*(-26*Ftmp365 + Ftmp367)*M[56] - Ftmp348*(-Ftmp11*Ftmp121 + Ftmp356)*M[18] - Ftmp349*Ftmp37*Ftmp52 - Ftmp351*y*(Ftmp362 + Ftmp366)*M[43] - Ftmp351*(Ftmp359 + Ftmp379 + Ftmp381 + Ftmp383)*M[54] - Ftmp351*(Ftmp361 + Ftmp376 + Ftmp382 + Ftmp383)*M[58] - Ftmp352*(Ftmp360 + Ftmp366)*M[41] - Ftmp354*Ftmp36*M[32] - Ftmp354*Ftmp78 - Ftmp368*y*(10*Ftmp11 + Ftmp375 + Ftmp380 + Ftmp382 + Ftmp61)*M[57] - Ftmp369*(Ftmp300 + Ftmp378 + Ftmp380 + Ftmp381 + 10*Ftmp8)*M[51] - Ftmp377*x*(Ftmp147*Ftmp8 + Ftmp373 + Ftmp376 + Ftmp8)*M[55] - Ftmp377*y*(Ftmp11*Ftmp134 + Ftmp11 + Ftmp373 + Ftmp379)*M[53] - Ftmp4*M[5] - Ftmp57*Ftmp60 - Ftmp57*Ftmp63 + 6*Ftmp7*Ftmp8*z*M[3] + 45*Ftmp7*Ftmp98*M[34] + 15*Ftmp7*x*y*z*M[4] + 15*Ftmp7*x*y*M[12] - Ftmp75*Ftmp84 + M[2]);

}

void M2Lc_7(double x, double y, double z, double * M, double * L) {
double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
double D[84];
double Dtmp0 = pow(Rinv, 3);
double Dtmp1 = pow(x, 2);
double Dtmp2 = pow(Rinv, 2);
double Dtmp3 = 3*Dtmp2;
double Dtmp4 = Dtmp1*Dtmp3;
double Dtmp5 = Dtmp4 - 1;
double Dtmp6 = 3*pow(Rinv, 5);
double Dtmp7 = Dtmp6*x;
double Dtmp8 = pow(y, 2);
double Dtmp9 = Dtmp3*Dtmp8;
double Dtmp10 = Dtmp9 - 1;
double Dtmp11 = Dtmp6*y;
double Dtmp12 = 5*Dtmp2;
double Dtmp13 = Dtmp1*Dtmp12;
double Dtmp14 = Dtmp13 - 1;
double Dtmp15 = Dtmp6*z;
double Dtmp16 = Dtmp12*Dtmp8;
double Dtmp17 = Dtmp16 - 1;
double Dtmp18 = pow(Rinv, 7);
double Dtmp19 = 15*Dtmp18;
double Dtmp20 = Dtmp19*x;
double Dtmp21 = Dtmp20*y;
double Dtmp22 = Dtmp1*Dtmp2;
double Dtmp23 = 30*Dtmp22;
double Dtmp24 = -Dtmp23;
double Dtmp25 = pow(x, 4);
double Dtmp26 = pow(Rinv, 4);
double Dtmp27 = 35*Dtmp26;
double Dtmp28 = 7*Dtmp22;
double Dtmp29 = Dtmp28 - 3;
double Dtmp30 = Dtmp20*z;
double Dtmp31 = Dtmp1*Dtmp8;
double Dtmp32 = Dtmp28 - 1;
double Dtmp33 = Dtmp19*y;
double Dtmp34 = Dtmp33*z;
double Dtmp35 = Dtmp2*Dtmp8;
double Dtmp36 = 7*Dtmp35;
double Dtmp37 = Dtmp36 - 3;
double Dtmp38 = Dtmp36 - 1;
double Dtmp39 = 30*Dtmp35;
double Dtmp40 = -Dtmp39;
double Dtmp41 = pow(y, 4);
double Dtmp42 = Dtmp25*Dtmp26;
double Dtmp43 = 14*Dtmp22;
double Dtmp44 = 21*Dtmp42;
double Dtmp45 = 45*Dtmp18;
double Dtmp46 = Dtmp45*(-Dtmp43 + Dtmp44 + 1);
double Dtmp47 = -Dtmp28;
double Dtmp48 = Dtmp26*Dtmp31;
double Dtmp49 = 63*Dtmp48;
double Dtmp50 = Dtmp49 + 3;
double Dtmp51 = pow(Rinv, 9);
double Dtmp52 = 315*Dtmp51;
double Dtmp53 = Dtmp52*x;
double Dtmp54 = Dtmp53*z;
double Dtmp55 = Dtmp54*y;
double Dtmp56 = -Dtmp36;
double Dtmp57 = 14*Dtmp35;
double Dtmp58 = Dtmp26*Dtmp41;
double Dtmp59 = 21*Dtmp58;
double Dtmp60 = -Dtmp57 + Dtmp59 + 1;
double Dtmp61 = Dtmp45*Dtmp60;
double Dtmp62 = pow(x, 6);
double Dtmp63 = pow(Rinv, 6);
double Dtmp64 = 231*Dtmp63;
double Dtmp65 = 33*Dtmp42;
double Dtmp66 = Dtmp53*(Dtmp24 + Dtmp65 + 5);
double Dtmp67 = -126*Dtmp48;
double Dtmp68 = Dtmp25*Dtmp8;
double Dtmp69 = 18*Dtmp22;
double Dtmp70 = Dtmp52*z;
double Dtmp71 = Dtmp70*y;
double Dtmp72 = -Dtmp9;
double Dtmp73 = 1 - Dtmp4;
double Dtmp74 = 945*Dtmp51;
double Dtmp75 = Dtmp74*y;
double Dtmp76 = 9*Dtmp35;
double Dtmp77 = 33*Dtmp48;
double Dtmp78 = Dtmp1*Dtmp41;
double Dtmp79 = 9*Dtmp22;
double Dtmp80 = 33*Dtmp58;
double Dtmp81 = Dtmp40 + Dtmp80 + 5;
double Dtmp82 = 4*Dtmp35;
double Dtmp83 = pow(y, 6);
double Dtmp84 = 429*Dtmp63;
double Dtmp85 = Dtmp62*Dtmp84;
double Dtmp86 = Dtmp52*(135*Dtmp22 - 495*Dtmp42 + Dtmp85 - 5);
double Dtmp87 = -Dtmp65;
double Dtmp88 = Dtmp68*Dtmp84 + Dtmp87;
double Dtmp89 = -330*Dtmp48 - 5;
double Dtmp90 = 945*pow(Rinv, 11)*x*y*z;
double Dtmp91 = -66*Dtmp48;
double Dtmp92 = 143*Dtmp63;
double Dtmp93 = -198*Dtmp48 - 1;
double Dtmp94 = -Dtmp80;
double Dtmp95 = 18*Dtmp35 + Dtmp94;
double Dtmp96 = Dtmp78*Dtmp84;
double Dtmp97 = Dtmp52*y;
double Dtmp98 = Dtmp83*Dtmp84;
double Dtmp99 = 135*Dtmp35 - 495*Dtmp58 + Dtmp98 - 5;
D[0] = -Dtmp0*x;
D[1] = -Dtmp0*y;
D[2] = -Dtmp0*z;
D[3] = Dtmp0*Dtmp5;
D[4] = Dtmp7*y;
D[5] = Dtmp7*z;
D[6] = Dtmp0*Dtmp10;
D[7] = Dtmp11*z;
D[8] = -D[3] - D[6];
D[9] = -Dtmp7*(Dtmp13 - 3);
D[10] = -Dtmp11*Dtmp14;
D[11] = -Dtmp14*Dtmp15;
D[12] = -Dtmp17*Dtmp7;
D[13] = -Dtmp21*z;
D[14] = -D[9] - D[12];
D[15] = -Dtmp11*(Dtmp16 - 3);
D[16] = -Dtmp15*Dtmp17;
D[17] = -D[10] - D[15];
D[18] = Dtmp6*(Dtmp24 + Dtmp25*Dtmp27 + 3);
D[19] = Dtmp21*Dtmp29;
D[20] = Dtmp29*Dtmp30;
D[21] = Dtmp6*(-Dtmp13 - Dtmp16 + Dtmp27*Dtmp31 + 1);
D[22] = Dtmp32*Dtmp34;
D[23] = -D[18] - D[21];
D[24] = Dtmp21*Dtmp37;
D[25] = Dtmp30*Dtmp38;
D[26] = -D[19] - D[24];
D[27] = Dtmp6*(Dtmp27*Dtmp41 + Dtmp40 + 3);
D[28] = Dtmp34*Dtmp37;
D[29] = -D[21] - D[27];
D[30] = -Dtmp20*(-70*Dtmp22 + 63*Dtmp42 + 15);
D[31] = -Dtmp46*y;
D[32] = -Dtmp46*z;
D[33] = -Dtmp20*(-21*Dtmp35 + Dtmp47 + Dtmp50);
D[34] = -Dtmp5*Dtmp55;
D[35] = -D[30] - D[33];
D[36] = -Dtmp33*(-21*Dtmp22 + Dtmp50 + Dtmp56);
D[37] = -Dtmp19*z*(Dtmp47 + Dtmp49 + Dtmp56 + 1);
D[38] = -D[31] - D[36];
D[39] = -Dtmp61*x;
D[40] = -Dtmp10*Dtmp55;
D[41] = -D[33] - D[39];
D[42] = -Dtmp33*(-70*Dtmp35 + 63*Dtmp58 + 15);
D[43] = -Dtmp61*z;
D[44] = -D[36] - D[42];
D[45] = Dtmp45*(105*Dtmp22 - 315*Dtmp42 + Dtmp62*Dtmp64 - 5);
D[46] = Dtmp66*y;
D[47] = Dtmp66*z;
D[48] = Dtmp45*(Dtmp38 + Dtmp43 - Dtmp44 + Dtmp64*Dtmp68 + Dtmp67);
D[49] = Dtmp71*(Dtmp65 - Dtmp69 + 1);
D[50] = -D[45] - D[48];
D[51] = Dtmp75*x*(11*Dtmp48 + Dtmp72 + Dtmp73);
D[52] = Dtmp54*(Dtmp73 - Dtmp76 + Dtmp77);
D[53] = -D[46] - D[51];
D[54] = Dtmp45*(Dtmp32 + Dtmp57 - Dtmp59 + Dtmp64*Dtmp78 + Dtmp67);
D[55] = Dtmp71*(Dtmp72 + Dtmp77 - Dtmp79 + 1);
D[56] = -D[48] - D[54];
D[57] = Dtmp53*Dtmp81*y;
D[58] = Dtmp54*(Dtmp10*Dtmp82 + Dtmp60);
D[59] = -D[51] - D[57];
D[60] = Dtmp45*(105*Dtmp35 - 315*Dtmp58 + Dtmp64*Dtmp83 - 5);
D[61] = Dtmp71*Dtmp81;
D[62] = -D[54] - D[60];
D[63] = -Dtmp53*(315*Dtmp22 - 693*Dtmp42 + Dtmp85 - 35);
D[64] = -Dtmp86*y;
D[65] = -Dtmp86*z;
D[66] = -Dtmp53*(Dtmp23 + 45*Dtmp35 + Dtmp88 + Dtmp89);
D[67] = -Dtmp90*(-110*Dtmp22 + 143*Dtmp42 + 15);
D[68] = -D[63] - D[66];
D[69] = -Dtmp75*(Dtmp10 + Dtmp68*Dtmp92 + Dtmp69 + Dtmp87 + Dtmp91);
D[70] = -Dtmp70*(Dtmp69 + Dtmp76 + Dtmp88 + Dtmp93);
D[71] = -D[64] - D[69];
D[72] = -Dtmp74*x*(Dtmp5 + Dtmp78*Dtmp92 + Dtmp91 + Dtmp95);
D[73] = -Dtmp90*(-33*Dtmp22 - 33*Dtmp35 + 143*Dtmp48 + 9);
D[74] = -D[66] - D[72];
D[75] = -Dtmp97*(45*Dtmp22 + Dtmp39 + Dtmp89 + Dtmp94 + Dtmp96);
D[76] = -Dtmp70*(Dtmp79 + Dtmp93 + Dtmp95 + Dtmp96);
D[77] = -D[69] - D[75];
D[78] = -Dtmp53*Dtmp99;
D[79] = -Dtmp90*(-90*Dtmp35 + 99*Dtmp58 + Dtmp82*(11*Dtmp35 - 5) + 15);
D[80] = -D[72] - D[78];
D[81] = -Dtmp97*(315*Dtmp35 - 693*Dtmp58 + Dtmp98 - 35);
D[82] = -Dtmp70*Dtmp99;
D[83] = -D[75] - D[81];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[9]*M[8] + D[10]*M[9] + D[11]*M[10] + D[12]*M[11] + D[13]*M[12] + D[15]*M[13] + D[16]*M[14] + D[18]*M[15] + D[19]*M[16] + D[20]*M[17] + D[21]*M[18] + D[22]*M[19] + D[24]*M[20] + D[25]*M[21] + D[27]*M[22] + D[28]*M[23] + D[30]*M[24] + D[31]*M[25] + D[32]*M[26] + D[33]*M[27] + D[34]*M[28] + D[36]*M[29] + D[37]*M[30] + D[39]*M[31] + D[40]*M[32] + D[42]*M[33] + D[43]*M[34] + D[45]*M[35] + D[46]*M[36] + D[47]*M[37] + D[48]*M[38] + D[49]*M[39] + D[51]*M[40] + D[52]*M[41] + D[54]*M[42] + D[55]*M[43] + D[57]*M[44] + D[58]*M[45] + D[60]*M[46] + D[61]*M[47] + D[63]*M[48] + D[64]*M[49] + D[65]*M[50] + D[66]*M[51] + D[67]*M[52] + D[69]*M[53] + D[70]*M[54] + D[72]*M[55] + D[73]*M[56] + D[75]*M[57] + D[76]*M[58] + D[78]*M[59] + D[79]*M[60] + D[81]*M[61] + D[82]*M[62];
#pragma omp atomic
L[1] += D[3]*M[0] + D[4]*M[1] + D[5]*M[2] + D[9]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[18]*M[8] + D[19]*M[9] + D[20]*M[10] + D[21]*M[11] + D[22]*M[12] + D[24]*M[13] + D[25]*M[14] + D[30]*M[15] + D[31]*M[16] + D[32]*M[17] + D[33]*M[18] + D[34]*M[19] + D[36]*M[20] + D[37]*M[21] + D[39]*M[22] + D[40]*M[23] + D[45]*M[24] + D[46]*M[25] + D[47]*M[26] + D[48]*M[27] + D[49]*M[28] + D[51]*M[29] + D[52]*M[30] + D[54]*M[31] + D[55]*M[32] + D[57]*M[33] + D[58]*M[34] + D[63]*M[35] + D[64]*M[36] + D[65]*M[37] + D[66]*M[38] + D[67]*M[39] + D[69]*M[40] + D[70]*M[41] + D[72]*M[42] + D[73]*M[43] + D[75]*M[44] + D[76]*M[45] + D[78]*M[46] + D[79]*M[47];
#pragma omp atomic
L[2] += D[4]*M[0] + D[6]*M[1] + D[7]*M[2] + D[10]*M[3] + D[12]*M[4] + D[13]*M[5] + D[15]*M[6] + D[16]*M[7] + D[19]*M[8] + D[21]*M[9] + D[22]*M[10] + D[24]*M[11] + D[25]*M[12] + D[27]*M[13] + D[28]*M[14] + D[31]*M[15] + D[33]*M[16] + D[34]*M[17] + D[36]*M[18] + D[37]*M[19] + D[39]*M[20] + D[40]*M[21] + D[42]*M[22] + D[43]*M[23] + D[46]*M[24] + D[48]*M[25] + D[49]*M[26] + D[51]*M[27] + D[52]*M[28] + D[54]*M[29] + D[55]*M[30] + D[57]*M[31] + D[58]*M[32] + D[60]*M[33] + D[61]*M[34] + D[64]*M[35] + D[66]*M[36] + D[67]*M[37] + D[69]*M[38] + D[70]*M[39] + D[72]*M[40] + D[73]*M[41] + D[75]*M[42] + D[76]*M[43] + D[78]*M[44] + D[79]*M[45] + D[81]*M[46] + D[82]*M[47];
#pragma omp atomic
L[3] += D[5]*M[0] + D[7]*M[1] + D[8]*M[2] + D[11]*M[3] + D[13]*M[4] + D[14]*M[5] + D[16]*M[6] + D[17]*M[7] + D[20]*M[8] + D[22]*M[9] + D[23]*M[10] + D[25]*M[11] + D[26]*M[12] + D[28]*M[13] + D[29]*M[14] + D[32]*M[15] + D[34]*M[16] + D[35]*M[17] + D[37]*M[18] + D[38]*M[19] + D[40]*M[20] + D[41]*M[21] + D[43]*M[22] + D[44]*M[23] + D[47]*M[24] + D[49]*M[25] + D[50]*M[26] + D[52]*M[27] + D[53]*M[28] + D[55]*M[29] + D[56]*M[30] + D[58]*M[31] + D[59]*M[32] + D[61]*M[33] + D[62]*M[34] + D[65]*M[35] + D[67]*M[36] + D[68]*M[37] + D[70]*M[38] + D[71]*M[39] + D[73]*M[40] + D[74]*M[41] + D[76]*M[42] + D[77]*M[43] + D[79]*M[44] + D[80]*M[45] + D[82]*M[46] + D[83]*M[47];
#pragma omp atomic
L[4] += D[9]*M[0] + D[10]*M[1] + D[11]*M[2] + D[18]*M[3] + D[19]*M[4] + D[20]*M[5] + D[21]*M[6] + D[22]*M[7] + D[30]*M[8] + D[31]*M[9] + D[32]*M[10] + D[33]*M[11] + D[34]*M[12] + D[36]*M[13] + D[37]*M[14] + D[45]*M[15] + D[46]*M[16] + D[47]*M[17] + D[48]*M[18] + D[49]*M[19] + D[51]*M[20] + D[52]*M[21] + D[54]*M[22] + D[55]*M[23] + D[63]*M[24] + D[64]*M[25] + D[65]*M[26] + D[66]*M[27] + D[67]*M[28] + D[69]*M[29] + D[70]*M[30] + D[72]*M[31] + D[73]*M[32] + D[75]*M[33] + D[76]*M[34];
#pragma omp atomic
L[5] += D[10]*M[0] + D[12]*M[1] + D[13]*M[2] + D[19]*M[3] + D[21]*M[4] + D[22]*M[5] + D[24]*M[6] + D[25]*M[7] + D[31]*M[8] + D[33]*M[9] + D[34]*M[10] + D[36]*M[11] + D[37]*M[12] + D[39]*M[13] + D[40]*M[14] + D[46]*M[15] + D[48]*M[16] + D[49]*M[17] + D[51]*M[18] + D[52]*M[19] + D[54]*M[20] + D[55]*M[21] + D[57]*M[22] + D[58]*M[23] + D[64]*M[24] + D[66]*M[25] + D[67]*M[26] + D[69]*M[27] + D[70]*M[28] + D[72]*M[29] + D[73]*M[30] + D[75]*M[31] + D[76]*M[32] + D[78]*M[33] + D[79]*M[34];
#pragma omp atomic
L[6] += D[11]*M[0] + D[13]*M[1] + D[14]*M[2] + D[20]*M[3] + D[22]*M[4] + D[23]*M[5] + D[25]*M[6] + D[26]*M[7] + D[32]*M[8] + D[34]*M[9] + D[35]*M[10] + D[37]*M[11] + D[38]*M[12] + D[40]*M[13] + D[41]*M[14] + D[47]*M[15] + D[49]*M[16] + D[50]*M[17] + D[52]*M[18] + D[53]*M[19] + D[55]*M[20] + D[56]*M[21] + D[58]*M[22] + D[59]*M[23] + D[65]*M[24] + D[67]*M[25] + D[68]*M[26] + D[70]*M[27] + D[71]*M[28] + D[73]*M[29] + D[74]*M[30] + D[76]*M[31] + D[77]*M[32] + D[79]*M[33] + D[80]*M[34];
#pragma omp atomic
L[7] += D[12]*M[0] + D[15]*M[1] + D[16]*M[2] + D[21]*M[3] + D[24]*M[4] + D[25]*M[5] + D[27]*M[6] + D[28]*M[7] + D[33]*M[8] + D[36]*M[9] + D[37]*M[10] + D[39]*M[11] + D[40]*M[12] + D[42]*M[13] + D[43]*M[14] + D[48]*M[15] + D[51]*M[16] + D[52]*M[17] + D[54]*M[18] + D[55]*M[19] + D[57]*M[20] + D[58]*M[21] + D[60]*M[22] + D[61]*M[23] + D[66]*M[24] + D[69]*M[25] + D[70]*M[26] + D[72]*M[27] + D[73]*M[28] + D[75]*M[29] + D[76]*M[30] + D[78]*M[31] + D[79]*M[32] + D[81]*M[33] + D[82]*M[34];
#pragma omp atomic
L[8] += D[13]*M[0] + D[16]*M[1] + D[17]*M[2] + D[22]*M[3] + D[25]*M[4] + D[26]*M[5] + D[28]*M[6] + D[29]*M[7] + D[34]*M[8] + D[37]*M[9] + D[38]*M[10] + D[40]*M[11] + D[41]*M[12] + D[43]*M[13] + D[44]*M[14] + D[49]*M[15] + D[52]*M[16] + D[53]*M[17] + D[55]*M[18] + D[56]*M[19] + D[58]*M[20] + D[59]*M[21] + D[61]*M[22] + D[62]*M[23] + D[67]*M[24] + D[70]*M[25] + D[71]*M[26] + D[73]*M[27] + D[74]*M[28] + D[76]*M[29] + D[77]*M[30] + D[79]*M[31] + D[80]*M[32] + D[82]*M[33] + D[83]*M[34];
#pragma omp atomic
L[9] += D[18]*M[0] + D[19]*M[1] + D[20]*M[2] + D[30]*M[3] + D[31]*M[4] + D[32]*M[5] + D[33]*M[6] + D[34]*M[7] + D[45]*M[8] + D[46]*M[9] + D[47]*M[10] + D[48]*M[11] + D[49]*M[12] + D[51]*M[13] + D[52]*M[14] + D[63]*M[15] + D[64]*M[16] + D[65]*M[17] + D[66]*M[18] + D[67]*M[19] + D[69]*M[20] + D[70]*M[21] + D[72]*M[22] + D[73]*M[23];
#pragma omp atomic
L[10] += D[19]*M[0] + D[21]*M[1] + D[22]*M[2] + D[31]*M[3] + D[33]*M[4] + D[34]*M[5] + D[36]*M[6] + D[37]*M[7] + D[46]*M[8] + D[48]*M[9] + D[49]*M[10] + D[51]*M[11] + D[52]*M[12] + D[54]*M[13] + D[55]*M[14] + D[64]*M[15] + D[66]*M[16] + D[67]*M[17] + D[69]*M[18] + D[70]*M[19] + D[72]*M[20] + D[73]*M[21] + D[75]*M[22] + D[76]*M[23];
#pragma omp atomic
L[11] += D[20]*M[0] + D[22]*M[1] + D[23]*M[2] + D[32]*M[3] + D[34]*M[4] + D[35]*M[5] + D[37]*M[6] + D[38]*M[7] + D[47]*M[8] + D[49]*M[9] + D[50]*M[10] + D[52]*M[11] + D[53]*M[12] + D[55]*M[13] + D[56]*M[14] + D[65]*M[15] + D[67]*M[16] + D[68]*M[17] + D[70]*M[18] + D[71]*M[19] + D[73]*M[20] + D[74]*M[21] + D[76]*M[22] + D[77]*M[23];
#pragma omp atomic
L[12] += D[21]*M[0] + D[24]*M[1] + D[25]*M[2] + D[33]*M[3] + D[36]*M[4] + D[37]*M[5] + D[39]*M[6] + D[40]*M[7] + D[48]*M[8] + D[51]*M[9] + D[52]*M[10] + D[54]*M[11] + D[55]*M[12] + D[57]*M[13] + D[58]*M[14] + D[66]*M[15] + D[69]*M[16] + D[70]*M[17] + D[72]*M[18] + D[73]*M[19] + D[75]*M[20] + D[76]*M[21] + D[78]*M[22] + D[79]*M[23];
#pragma omp atomic
L[13] += D[22]*M[0] + D[25]*M[1] + D[26]*M[2] + D[34]*M[3] + D[37]*M[4] + D[38]*M[5] + D[40]*M[6] + D[41]*M[7] + D[49]*M[8] + D[52]*M[9] + D[53]*M[10] + D[55]*M[11] + D[56]*M[12] + D[58]*M[13] + D[59]*M[14] + D[67]*M[15] + D[70]*M[16] + D[71]*M[17] + D[73]*M[18] + D[74]*M[19] + D[76]*M[20] + D[77]*M[21] + D[79]*M[22] + D[80]*M[23];
#pragma omp atomic
L[14] += D[24]*M[0] + D[27]*M[1] + D[28]*M[2] + D[36]*M[3] + D[39]*M[4] + D[40]*M[5] + D[42]*M[6] + D[43]*M[7] + D[51]*M[8] + D[54]*M[9] + D[55]*M[10] + D[57]*M[11] + D[58]*M[12] + D[60]*M[13] + D[61]*M[14] + D[69]*M[15] + D[72]*M[16] + D[73]*M[17] + D[75]*M[18] + D[76]*M[19] + D[78]*M[20] + D[79]*M[21] + D[81]*M[22] + D[82]*M[23];
#pragma omp atomic
L[15] += D[25]*M[0] + D[28]*M[1] + D[29]*M[2] + D[37]*M[3] + D[40]*M[4] + D[41]*M[5] + D[43]*M[6] + D[44]*M[7] + D[52]*M[8] + D[55]*M[9] + D[56]*M[10] + D[58]*M[11] + D[59]*M[12] + D[61]*M[13] + D[62]*M[14] + D[70]*M[15] + D[73]*M[16] + D[74]*M[17] + D[76]*M[18] + D[77]*M[19] + D[79]*M[20] + D[80]*M[21] + D[82]*M[22] + D[83]*M[23];
#pragma omp atomic
L[16] += D[30]*M[0] + D[31]*M[1] + D[32]*M[2] + D[45]*M[3] + D[46]*M[4] + D[47]*M[5] + D[48]*M[6] + D[49]*M[7] + D[63]*M[8] + D[64]*M[9] + D[65]*M[10] + D[66]*M[11] + D[67]*M[12] + D[69]*M[13] + D[70]*M[14];
#pragma omp atomic
L[17] += D[31]*M[0] + D[33]*M[1] + D[34]*M[2] + D[46]*M[3] + D[48]*M[4] + D[49]*M[5] + D[51]*M[6] + D[52]*M[7] + D[64]*M[8] + D[66]*M[9] + D[67]*M[10] + D[69]*M[11] + D[70]*M[12] + D[72]*M[13] + D[73]*M[14];
#pragma omp atomic
L[18] += D[32]*M[0] + D[34]*M[1] + D[35]*M[2] + D[47]*M[3] + D[49]*M[4] + D[50]*M[5] + D[52]*M[6] + D[53]*M[7] + D[65]*M[8] + D[67]*M[9] + D[68]*M[10] + D[70]*M[11] + D[71]*M[12] + D[73]*M[13] + D[74]*M[14];
#pragma omp atomic
L[19] += D[33]*M[0] + D[36]*M[1] + D[37]*M[2] + D[48]*M[3] + D[51]*M[4] + D[52]*M[5] + D[54]*M[6] + D[55]*M[7] + D[66]*M[8] + D[69]*M[9] + D[70]*M[10] + D[72]*M[11] + D[73]*M[12] + D[75]*M[13] + D[76]*M[14];
#pragma omp atomic
L[20] += D[34]*M[0] + D[37]*M[1] + D[38]*M[2] + D[49]*M[3] + D[52]*M[4] + D[53]*M[5] + D[55]*M[6] + D[56]*M[7] + D[67]*M[8] + D[70]*M[9] + D[71]*M[10] + D[73]*M[11] + D[74]*M[12] + D[76]*M[13] + D[77]*M[14];
#pragma omp atomic
L[21] += D[36]*M[0] + D[39]*M[1] + D[40]*M[2] + D[51]*M[3] + D[54]*M[4] + D[55]*M[5] + D[57]*M[6] + D[58]*M[7] + D[69]*M[8] + D[72]*M[9] + D[73]*M[10] + D[75]*M[11] + D[76]*M[12] + D[78]*M[13] + D[79]*M[14];
#pragma omp atomic
L[22] += D[37]*M[0] + D[40]*M[1] + D[41]*M[2] + D[52]*M[3] + D[55]*M[4] + D[56]*M[5] + D[58]*M[6] + D[59]*M[7] + D[70]*M[8] + D[73]*M[9] + D[74]*M[10] + D[76]*M[11] + D[77]*M[12] + D[79]*M[13] + D[80]*M[14];
#pragma omp atomic
L[23] += D[39]*M[0] + D[42]*M[1] + D[43]*M[2] + D[54]*M[3] + D[57]*M[4] + D[58]*M[5] + D[60]*M[6] + D[61]*M[7] + D[72]*M[8] + D[75]*M[9] + D[76]*M[10] + D[78]*M[11] + D[79]*M[12] + D[81]*M[13] + D[82]*M[14];
#pragma omp atomic
L[24] += D[40]*M[0] + D[43]*M[1] + D[44]*M[2] + D[55]*M[3] + D[58]*M[4] + D[59]*M[5] + D[61]*M[6] + D[62]*M[7] + D[73]*M[8] + D[76]*M[9] + D[77]*M[10] + D[79]*M[11] + D[80]*M[12] + D[82]*M[13] + D[83]*M[14];
#pragma omp atomic
L[25] += D[45]*M[0] + D[46]*M[1] + D[47]*M[2] + D[63]*M[3] + D[64]*M[4] + D[65]*M[5] + D[66]*M[6] + D[67]*M[7];
#pragma omp atomic
L[26] += D[46]*M[0] + D[48]*M[1] + D[49]*M[2] + D[64]*M[3] + D[66]*M[4] + D[67]*M[5] + D[69]*M[6] + D[70]*M[7];
#pragma omp atomic
L[27] += D[47]*M[0] + D[49]*M[1] + D[50]*M[2] + D[65]*M[3] + D[67]*M[4] + D[68]*M[5] + D[70]*M[6] + D[71]*M[7];
#pragma omp atomic
L[28] += D[48]*M[0] + D[51]*M[1] + D[52]*M[2] + D[66]*M[3] + D[69]*M[4] + D[70]*M[5] + D[72]*M[6] + D[73]*M[7];
#pragma omp atomic
L[29] += D[49]*M[0] + D[52]*M[1] + D[53]*M[2] + D[67]*M[3] + D[70]*M[4] + D[71]*M[5] + D[73]*M[6] + D[74]*M[7];
#pragma omp atomic
L[30] += D[51]*M[0] + D[54]*M[1] + D[55]*M[2] + D[69]*M[3] + D[72]*M[4] + D[73]*M[5] + D[75]*M[6] + D[76]*M[7];
#pragma omp atomic
L[31] += D[52]*M[0] + D[55]*M[1] + D[56]*M[2] + D[70]*M[3] + D[73]*M[4] + D[74]*M[5] + D[76]*M[6] + D[77]*M[7];
#pragma omp atomic
L[32] += D[54]*M[0] + D[57]*M[1] + D[58]*M[2] + D[72]*M[3] + D[75]*M[4] + D[76]*M[5] + D[78]*M[6] + D[79]*M[7];
#pragma omp atomic
L[33] += D[55]*M[0] + D[58]*M[1] + D[59]*M[2] + D[73]*M[3] + D[76]*M[4] + D[77]*M[5] + D[79]*M[6] + D[80]*M[7];
#pragma omp atomic
L[34] += D[57]*M[0] + D[60]*M[1] + D[61]*M[2] + D[75]*M[3] + D[78]*M[4] + D[79]*M[5] + D[81]*M[6] + D[82]*M[7];
#pragma omp atomic
L[35] += D[58]*M[0] + D[61]*M[1] + D[62]*M[2] + D[76]*M[3] + D[79]*M[4] + D[80]*M[5] + D[82]*M[6] + D[83]*M[7];
#pragma omp atomic
L[36] += D[63]*M[0] + D[64]*M[1] + D[65]*M[2];
#pragma omp atomic
L[37] += D[64]*M[0] + D[66]*M[1] + D[67]*M[2];
#pragma omp atomic
L[38] += D[65]*M[0] + D[67]*M[1] + D[68]*M[2];
#pragma omp atomic
L[39] += D[66]*M[0] + D[69]*M[1] + D[70]*M[2];
#pragma omp atomic
L[40] += D[67]*M[0] + D[70]*M[1] + D[71]*M[2];
#pragma omp atomic
L[41] += D[69]*M[0] + D[72]*M[1] + D[73]*M[2];
#pragma omp atomic
L[42] += D[70]*M[0] + D[73]*M[1] + D[74]*M[2];
#pragma omp atomic
L[43] += D[72]*M[0] + D[75]*M[1] + D[76]*M[2];
#pragma omp atomic
L[44] += D[73]*M[0] + D[76]*M[1] + D[77]*M[2];
#pragma omp atomic
L[45] += D[75]*M[0] + D[78]*M[1] + D[79]*M[2];
#pragma omp atomic
L[46] += D[76]*M[0] + D[79]*M[1] + D[80]*M[2];
#pragma omp atomic
L[47] += D[78]*M[0] + D[81]*M[1] + D[82]*M[2];
#pragma omp atomic
L[48] += D[79]*M[0] + D[82]*M[1] + D[83]*M[2];

}

void S2M(double x, double y, double z, double * S, double * M, int order) {
switch (order) {
  case 2:
    S2M_2(x, y, z, S, M);
    break;
  case 3:
    S2M_3(x, y, z, S, M);
    break;
  case 4:
    S2M_4(x, y, z, S, M);
    break;
  case 5:
    S2M_5(x, y, z, S, M);
    break;
  case 6:
    S2M_6(x, y, z, S, M);
    break;
  case 7:
    S2M_7(x, y, z, S, M);
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
void S2Mc(double x, double y, double z, double * S, double * M, int order) {
switch (order) {
  case 2:
    S2Mc_2(x, y, z, S, M);
    break;
  case 3:
    S2Mc_3(x, y, z, S, M);
    break;
  case 4:
    S2Mc_4(x, y, z, S, M);
    break;
  case 5:
    S2Mc_5(x, y, z, S, M);
    break;
  case 6:
    S2Mc_6(x, y, z, S, M);
    break;
  case 7:
    S2Mc_7(x, y, z, S, M);
    break;
  }
}
void M2Mc(double x, double y, double z, double * M, double * Ms, int order) {
switch (order) {
  case 2:
    M2Mc_2(x, y, z, M, Ms);
    break;
  case 3:
    M2Mc_3(x, y, z, M, Ms);
    break;
  case 4:
    M2Mc_4(x, y, z, M, Ms);
    break;
  case 5:
    M2Mc_5(x, y, z, M, Ms);
    break;
  case 6:
    M2Mc_6(x, y, z, M, Ms);
    break;
  case 7:
    M2Mc_7(x, y, z, M, Ms);
    break;
  }
}
void L2Lc(double x, double y, double z, double * L, double * Ls, int order) {
switch (order) {
  case 2:
    L2Lc_2(x, y, z, L, Ls);
    break;
  case 3:
    L2Lc_3(x, y, z, L, Ls);
    break;
  case 4:
    L2Lc_4(x, y, z, L, Ls);
    break;
  case 5:
    L2Lc_5(x, y, z, L, Ls);
    break;
  case 6:
    L2Lc_6(x, y, z, L, Ls);
    break;
  case 7:
    L2Lc_7(x, y, z, L, Ls);
    break;
  }
}
void L2Pc(double x, double y, double z, double * L, double * F, int order) {
switch (order) {
  case 2:
    L2Pc_2(x, y, z, L, F);
    break;
  case 3:
    L2Pc_3(x, y, z, L, F);
    break;
  case 4:
    L2Pc_4(x, y, z, L, F);
    break;
  case 5:
    L2Pc_5(x, y, z, L, F);
    break;
  case 6:
    L2Pc_6(x, y, z, L, F);
    break;
  case 7:
    L2Pc_7(x, y, z, L, F);
    break;
  }
}
void M2Pc(double x, double y, double z, double * M, double * F, int order) {
switch (order) {
  case 2:
    M2Pc_2(x, y, z, M, F);
    break;
  case 3:
    M2Pc_3(x, y, z, M, F);
    break;
  case 4:
    M2Pc_4(x, y, z, M, F);
    break;
  case 5:
    M2Pc_5(x, y, z, M, F);
    break;
  case 6:
    M2Pc_6(x, y, z, M, F);
    break;
  case 7:
    M2Pc_7(x, y, z, M, F);
    break;
  }
}
void M2Lc(double x, double y, double z, double * M, double * L, int order) {
switch (order) {
  case 2:
    M2Lc_2(x, y, z, M, L);
    break;
  case 3:
    M2Lc_3(x, y, z, M, L);
    break;
  case 4:
    M2Lc_4(x, y, z, M, L);
    break;
  case 5:
    M2Lc_5(x, y, z, M, L);
    break;
  case 6:
    M2Lc_6(x, y, z, M, L);
    break;
  case 7:
    M2Lc_7(x, y, z, M, L);
    break;
  }
}
