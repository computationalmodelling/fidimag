#pragma once
#include<cmath>


void cross(double *output, double * A, double * B);
double dot(double *A, double *B, int n);

double compute_norm(double * a, int n);

void normalise(double * a, int n);
void normalise(double *m, int *pins, int n);

double cross_x(double a0, double a1, double a2, double b0, double b1, double b2);
double cross_y(double a0, double a1, double a2, double b0, double b1, double b2);
double cross_z(double a0, double a1, double a2, double b0, double b1, double b2);
