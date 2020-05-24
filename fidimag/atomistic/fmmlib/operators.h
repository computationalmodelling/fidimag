#pragma once
#define FMMGEN_MINORDER 2
#define FMMGEN_MAXORDER 8
#define FMMGEN_SOURCEORDER 1
#define FMMGEN_SOURCESIZE 3
#define FMMGEN_OUTPUTSIZE 3
void P2M_2(double x, double y, double z, double q, double * M);
void M2M_2(double x, double y, double z, double * M, double * Ms);
void M2L_2(double x, double y, double z, double * M, double * L);
void L2L_2(double x, double y, double z, double * L, double * Ls);
void L2P_2(double x, double y, double z, double * L, double * F);
void M2P_2(double x, double y, double z, double * M, double * F);
void P2P(double x, double y, double z, double * S, double * F);
void P2M_3(double x, double y, double z, double q, double * M);
void M2M_3(double x, double y, double z, double * M, double * Ms);
void M2L_3(double x, double y, double z, double * M, double * L);
void L2L_3(double x, double y, double z, double * L, double * Ls);
void L2P_3(double x, double y, double z, double * L, double * F);
void M2P_3(double x, double y, double z, double * M, double * F);
void P2M_4(double x, double y, double z, double q, double * M);
void M2M_4(double x, double y, double z, double * M, double * Ms);
void M2L_4(double x, double y, double z, double * M, double * L);
void L2L_4(double x, double y, double z, double * L, double * Ls);
void L2P_4(double x, double y, double z, double * L, double * F);
void M2P_4(double x, double y, double z, double * M, double * F);
void P2M_5(double x, double y, double z, double q, double * M);
void M2M_5(double x, double y, double z, double * M, double * Ms);
void M2L_5(double x, double y, double z, double * M, double * L);
void L2L_5(double x, double y, double z, double * L, double * Ls);
void L2P_5(double x, double y, double z, double * L, double * F);
void M2P_5(double x, double y, double z, double * M, double * F);
void P2M_6(double x, double y, double z, double q, double * M);
void M2M_6(double x, double y, double z, double * M, double * Ms);
void M2L_6(double x, double y, double z, double * M, double * L);
void L2L_6(double x, double y, double z, double * L, double * Ls);
void L2P_6(double x, double y, double z, double * L, double * F);
void M2P_6(double x, double y, double z, double * M, double * F);
void P2M_7(double x, double y, double z, double q, double * M);
void M2M_7(double x, double y, double z, double * M, double * Ms);
void M2L_7(double x, double y, double z, double * M, double * L);
void L2L_7(double x, double y, double z, double * L, double * Ls);
void L2P_7(double x, double y, double z, double * L, double * F);
void M2P_7(double x, double y, double z, double * M, double * F);
void P2M(double x, double y, double z, double q, double * M, int order);
void M2M(double x, double y, double z, double * M, double * Ms, int order);
void M2L(double x, double y, double z, double * M, double * L, int order);
void L2L(double x, double y, double z, double * L, double * Ls, int order);
void L2P(double x, double y, double z, double * L, double * F, int order);
void M2P(double x, double y, double z, double * M, double * F, int order);
