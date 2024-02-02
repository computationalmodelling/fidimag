/* This file demag_oommf.h is taken from the OOMMF project (oommf/app/oxs/ext/demagcoef.h
downloaded from http://math.nist.gov/oommf/dist/oommf12a5rc_20120928.tar.gz)
with slightly modifications (change OC_REALWIDE to double)
and distributed under this license shown below. */

/* License

OOMMF - Object Oriented MicroMagnetic Framework

The research software provided in this release (“software”) is provided by
NIST as a public service. You may use, copy and distribute copies of the
  software in any medium, provided that you keep intact this entire notice.
  You may improve, modify and create derivative works of the software or any
  portion of the software, and you may copy and distribute such modifications
or works. Modified works should carry a notice stating that you changed the
software and should note the date and nature of any such change. Please
explicitly acknowledge the National Institute of Standards and Technology
as the source of the software.

The software is expressly provided “AS IS.” NIST MAKES NO WARRANTY OF ANY
  KIND, EXPRESS, IMPLIED, IN FACT OR ARISING BY OPERATION OF LAW, INCLUDING,
  WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A
  PARTICULAR PURPOSE, NON-INFRINGEMENT AND DATA ACCURACY. NIST NEITHER
REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE
  UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST
DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE
  SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE
  CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.

You are solely responsible for determining the appropriateness of using
and distributing the software and you assume all risks associated with
  its use, including but not limited to the risks and costs of program
  errors, compliance with applicable laws, damage to or loss of data,
  programs or equipment, and the unavailability or interruption of operation.
This software is not intended to be used in any situation where a failure
could cause risk of injury or damage to property. The software was
developed by NIST employees. NIST employee contributions are not subject
to copyright protection within the United States.

*/

/* FILE: demagcoef.h            -*-Mode: c++-*-
 *
 * Demag coefficients.
 *
 * Constant magnetization demag routines, based on formulae presented in
 * "A Generalization of the Demagnetizing Tensor for Nonuniform
 * Magnetization," by Andrew J. Newell, Wyn Williams, and David
 * J. Dunlop, Journal of Geophysical Research, vol 98, p 9551-9555, June
 * 1993.  This formulae clearly satisfy necessary symmetry and scaling
 * properties, which is not true of the formulae presented in
 * "Magnetostatic Energy Calculations in Two- and Three-Dimensional
 * Arrays of Ferromagnetic Prisms," M. Maicas, E. Lopez, M. CC. Sanchez,
 * C. Aroca and P. Sanchez, IEEE Trans Mag, vol 34, May 1998, p601-607.
 * (Side note: The units in the latter paper are apparently cgs.)  It
 * appears likely that there is an error in the latter paper (attempts
 * to implement the formulae presented there did not produce the proper
 * symmetries), as well as in the older paper, "Magnetostatic
 * Interaction Fields for a Three-Dimensional Array of Ferromagnetic
 * Cubes," Manfred E. Schabes and Amikam Aharoni, IEEE Trans Mag, vol
 * 23, November 1987, p3882-3888.  (Note: The Newell paper deals with
 * uniformly sized rectangular prisms, the Maicas paper allows
 * non-uniformly sized rectangular prisms, and the Schabes paper only
 * considers cubes.)
 *
 *   The kernel here is based on an analytically derived energy, and the
 * effective (discrete) demag field is calculated from the (discrete)
 * energy.
 *
 */

#include <math.h>

double Newell_f(double x, double y, double z);
double Newell_g(double x, double y, double z);

double
CalculateSDA00(double x, double y, double z,
               double dx, double dy, double dz);

inline double
CalculateSDA11(double x, double y, double z,
               double dx, double dy, double dz) { return CalculateSDA00(y, x, z, dy, dx, dz); }

inline double
CalculateSDA22(double x, double y, double z,
               double dx, double dy, double dz) { return CalculateSDA00(z, y, x, dz, dy, dx); }

double
CalculateSDA01(double x, double y, double z,
               double dx, double dy, double dz);

inline double
CalculateSDA02(double x, double y, double z,
               double dx, double dy, double dz) { return CalculateSDA01(x, z, y, dx, dz, dy); }

inline double
CalculateSDA12(double x, double y, double z,
               double dx, double dy, double dz) { return CalculateSDA01(y, z, x, dy, dz, dx); }

double SelfDemagNx(double xsize, double ysize, double zsize);
double SelfDemagNy(double xsize, double ysize, double zsize);
double SelfDemagNz(double xsize, double ysize, double zsize);

double DemagNxxAsymptotic(double x, double y, double z,
                          double dx, double dy, double dz);

double DemagNxyAsymptotic(double x, double y, double z,
                          double dx, double dy, double dz);

double DemagNyyAsymptotic(double x, double y, double z,
                          double dx, double dy, double dz);

double DemagNzzAsymptotic(double x, double y, double z,
                          double dx, double dy, double dz);

double DemagNxzAsymptotic(double x, double y, double z,
                          double dx, double dy, double dz);

double DemagNyzAsymptotic(double x, double y, double z,
                          double dx, double dy, double dz);

////////////////////////////////////////////////////////////////////////////
// Routines to do accurate summation

static int AS_Compare(const void *px, const void *py) {
  // Comparison based on absolute values
  double x = fabs(*((const double *)px));
  double y = fabs(*((const double *)py));
  if (x < y)
    return 1;
  if (x > y)
    return -1;
  return 0;
}

static double
AccurateSum(int n, double *arr) {
  // Order by decreasing magnitude
  qsort(arr, n, sizeof(double), AS_Compare);

  // Add up using doubly compensated summation.  If necessary, mark
  // variables these "volatile" to protect against problems arising
  // from extra precision.  Also, don't expect the compiler to respect
  // order with respect to parentheses at high levels of optimization,
  // i.e., write "u=x; u-=(y-corr)" as opposed to "u=x-(y-corr)".

  double sum, corr, y, u, t, v, z, x, tmp;

  sum = arr[0];
  corr = 0;
  for (int i = 1; i < n; i++) {
    x = arr[i];
    y = corr + x;
    tmp = y - corr;
    u = x - tmp;
    t = y + sum;
    tmp = t - sum;
    v = y - tmp;
    z = u + v;
    sum = t + z;
    tmp = sum - t;
    corr = z - tmp;
  }
  return sum;
}
