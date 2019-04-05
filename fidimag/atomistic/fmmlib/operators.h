#ifndef OPERATORS_H
#define OPERATORS_H
void P2M_0(double x, double y, double z, double q, double * M);
void M2M_0(double x, double y, double z, double * M, double * Ms);
void M2L_0(double x, double y, double z, double * M, double * L);
void L2L_0(double x, double y, double z, double * L, double * Ls);
void L2P_0(double x, double y, double z, double * L, double * F);
void M2P_0(double x, double y, double z, double * M, double * F);
void P2M_1(double x, double y, double z, double q, double * M);
void M2M_1(double x, double y, double z, double * M, double * Ms);
void M2L_1(double x, double y, double z, double * M, double * L);
void L2L_1(double x, double y, double z, double * L, double * Ls);
void L2P_1(double x, double y, double z, double * L, double * F);
void M2P_1(double x, double y, double z, double * M, double * F);
void P2M_2(double x, double y, double z, double q, double * M);
void M2M_2(double x, double y, double z, double * M, double * Ms);
void M2L_2(double x, double y, double z, double * M, double * L);
void L2L_2(double x, double y, double z, double * L, double * Ls);
void L2P_2(double x, double y, double z, double * L, double * F);
void M2P_2(double x, double y, double z, double * M, double * F);
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
void P2M(double x, double y, double z, double q, double * M, int order);
void M2M(double x, double y, double z, double * M, double * Ms, int order);
void M2L(double x, double y, double z, double * M, double * L, int order);
void L2L(double x, double y, double z, double * L, double * Ls, int order);
void L2P(double x, double y, double z, double * L, double * F, int order);
void M2P(double x, double y, double z, double * M, double * F, int order);
#endif