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


double Newell_f(double x,double y,double z);
double Newell_g(double x,double y,double z);

double
CalculateSDA00(double x, double y, double z,
		double dx,double dy,double dz);

inline double
CalculateSDA11(double x, double y, double z,
	       double dx,double dy,double dz)
{ return CalculateSDA00(y,x,z,dy,dx,dz); }

inline double
CalculateSDA22(double x, double y, double z,
	       double dx,double dy,double dz)
{ return CalculateSDA00(z,y,x,dz,dy,dx); }


double
CalculateSDA01(double x, double y, double z,
	       double dx,double dy,double dz);

inline double
CalculateSDA02(double x, double y, double z,
	       double dx,double dy,double dz)
{ return CalculateSDA01(x,z,y,dx,dz,dy); }

inline double
CalculateSDA12(double x, double y, double z,
	       double dx,double dy,double dz)
{ return CalculateSDA01(y,z,x,dy,dz,dx); }


double SelfDemagNx(double xsize,double ysize,double zsize);
double SelfDemagNy(double xsize,double ysize,double zsize);
double SelfDemagNz(double xsize,double ysize,double zsize);

double DemagNxxAsymptotic(double x, double y, double z,
                            double dx,double dy,double dz);

double DemagNxyAsymptotic(double x, double y, double z,
                            double dx,double dy,double dz);

double DemagNyyAsymptotic(double x, double y, double z,
                            double dx,double dy,double dz);

double DemagNzzAsymptotic(double x, double y, double z,
                            double dx,double dy,double dz);

double DemagNxzAsymptotic(double x, double y, double z,
                            double dx,double dy,double dz);

double DemagNyzAsymptotic(double x, double y, double z,
                            double dx,double dy,double dz);
