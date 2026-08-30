#pragma once
#include <stddef.h>
#define FMMGEN_MINORDER 2
#define FMMGEN_MAXORDER 8
#define FMMGEN_SOURCEORDER 1
#define FMMGEN_SOURCESIZE 3
#define FMMGEN_OUTPUTSIZE 3
#define FMMGEN_COMPRESSED 1
/* Array sizes for the *c-suffixed* operators only. The trace-free
   basis drops the C(p+3,3) - (p+1)^2 redundant coefficients, so these
   are (p+1)^2; the unsuffixed operators still take Nterms(order).
   Index by order. */
static const size_t FMMGEN_MULTIPOLESIZE[] = {0, 3, 8, 15, 24, 35, 48, 63, 80};
static const size_t FMMGEN_LOCALSIZE[] = {0, 1, 4, 9, 16, 25, 36, 49, 64};
void S2M_2(double x, double y, double z, double * S, double * M);
void M2M_2(double x, double y, double z, double * M, double * Ms);
void M2L_2(double x, double y, double z, double * M, double * L);
void L2L_2(double x, double y, double z, double * L, double * Ls);
void L2P_2(double x, double y, double z, double * L, double * F);
void M2P_2(double x, double y, double z, double * M, double * F);
void S2Mc_2(double x, double y, double z, double * S, double * M);
void M2Mc_2(double x, double y, double z, double * M, double * Ms);
void L2Lc_2(double x, double y, double z, double * L, double * Ls);
void L2Pc_2(double x, double y, double z, double * L, double * F);
void M2Pc_2(double x, double y, double z, double * M, double * F);
void M2Lc_2(double x, double y, double z, double * M, double * L);
void P2P(double x, double y, double z, double * S, double * F);
void P2P_batch(double tx, double ty, double tz, const double * sx, const double * sy, const double * sz, const double * S, size_t begin, size_t end, double * F);
void S2M_3(double x, double y, double z, double * S, double * M);
void M2M_3(double x, double y, double z, double * M, double * Ms);
void M2L_3(double x, double y, double z, double * M, double * L);
void L2L_3(double x, double y, double z, double * L, double * Ls);
void L2P_3(double x, double y, double z, double * L, double * F);
void M2P_3(double x, double y, double z, double * M, double * F);
void S2Mc_3(double x, double y, double z, double * S, double * M);
void M2Mc_3(double x, double y, double z, double * M, double * Ms);
void L2Lc_3(double x, double y, double z, double * L, double * Ls);
void L2Pc_3(double x, double y, double z, double * L, double * F);
void M2Pc_3(double x, double y, double z, double * M, double * F);
void M2Lc_3(double x, double y, double z, double * M, double * L);
void S2M_4(double x, double y, double z, double * S, double * M);
void M2M_4(double x, double y, double z, double * M, double * Ms);
void M2L_4(double x, double y, double z, double * M, double * L);
void L2L_4(double x, double y, double z, double * L, double * Ls);
void L2P_4(double x, double y, double z, double * L, double * F);
void M2P_4(double x, double y, double z, double * M, double * F);
void S2Mc_4(double x, double y, double z, double * S, double * M);
void M2Mc_4(double x, double y, double z, double * M, double * Ms);
void L2Lc_4(double x, double y, double z, double * L, double * Ls);
void L2Pc_4(double x, double y, double z, double * L, double * F);
void M2Pc_4(double x, double y, double z, double * M, double * F);
void M2Lc_4(double x, double y, double z, double * M, double * L);
void S2M_5(double x, double y, double z, double * S, double * M);
void M2M_5(double x, double y, double z, double * M, double * Ms);
void M2L_5(double x, double y, double z, double * M, double * L);
void L2L_5(double x, double y, double z, double * L, double * Ls);
void L2P_5(double x, double y, double z, double * L, double * F);
void M2P_5(double x, double y, double z, double * M, double * F);
void S2Mc_5(double x, double y, double z, double * S, double * M);
void M2Mc_5(double x, double y, double z, double * M, double * Ms);
void L2Lc_5(double x, double y, double z, double * L, double * Ls);
void L2Pc_5(double x, double y, double z, double * L, double * F);
void M2Pc_5(double x, double y, double z, double * M, double * F);
void M2Lc_5(double x, double y, double z, double * M, double * L);
void S2M_6(double x, double y, double z, double * S, double * M);
void M2M_6(double x, double y, double z, double * M, double * Ms);
void M2L_6(double x, double y, double z, double * M, double * L);
void L2L_6(double x, double y, double z, double * L, double * Ls);
void L2P_6(double x, double y, double z, double * L, double * F);
void M2P_6(double x, double y, double z, double * M, double * F);
void S2Mc_6(double x, double y, double z, double * S, double * M);
void M2Mc_6(double x, double y, double z, double * M, double * Ms);
void L2Lc_6(double x, double y, double z, double * L, double * Ls);
void L2Pc_6(double x, double y, double z, double * L, double * F);
void M2Pc_6(double x, double y, double z, double * M, double * F);
void M2Lc_6(double x, double y, double z, double * M, double * L);
void S2M_7(double x, double y, double z, double * S, double * M);
void M2M_7(double x, double y, double z, double * M, double * Ms);
void M2L_7(double x, double y, double z, double * M, double * L);
void L2L_7(double x, double y, double z, double * L, double * Ls);
void L2P_7(double x, double y, double z, double * L, double * F);
void M2P_7(double x, double y, double z, double * M, double * F);
void S2Mc_7(double x, double y, double z, double * S, double * M);
void M2Mc_7(double x, double y, double z, double * M, double * Ms);
void L2Lc_7(double x, double y, double z, double * L, double * Ls);
void L2Pc_7(double x, double y, double z, double * L, double * F);
void M2Pc_7(double x, double y, double z, double * M, double * F);
void M2Lc_7(double x, double y, double z, double * M, double * L);
void S2M(double x, double y, double z, double * S, double * M, int order);
void M2M(double x, double y, double z, double * M, double * Ms, int order);
void M2L(double x, double y, double z, double * M, double * L, int order);
void L2L(double x, double y, double z, double * L, double * Ls, int order);
void L2P(double x, double y, double z, double * L, double * F, int order);
void M2P(double x, double y, double z, double * M, double * F, int order);
void S2Mc(double x, double y, double z, double * S, double * M, int order);
void M2Mc(double x, double y, double z, double * M, double * Ms, int order);
void L2Lc(double x, double y, double z, double * L, double * Ls, int order);
void L2Pc(double x, double y, double z, double * L, double * F, int order);
void M2Pc(double x, double y, double z, double * M, double * F, int order);
void M2Lc(double x, double y, double z, double * M, double * L, int order);
