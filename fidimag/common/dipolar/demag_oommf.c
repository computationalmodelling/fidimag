/* This file demag_oommf.c is taken from the OOMMF project (oommf/app/oxs/ext/demagcoef.cc 
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


#include <math.h>
#include <stdlib.h>
#include "dipolar.h"


/* FILE: demagcoef.cc            -*-Mode: c++-*-
 *
 * Demag coefficients.
 *
 * Constant magnetization demag routines, based on formulae presented in
 * "A Generalization of the Demagnetizing Tensor for Nonuniform
 * Magnetization," by Andrew J. Newell, Wyn Williams, and David
 * J. Dunlop, Journal of Geophysical Research - Solid Earth, vol 98,
 * p 9551-9555, June 1993.  This formulae clearly satisfy necessary
 * symmetry and scaling properties, which is not true of the formulae
 * presented in "Magnetostatic Energy Calculations in Two- and
 * Three-Dimensional Arrays of Ferromagnetic Prisms," M. Maicas,
 * E. Lopez, M. CC. Sanchez, C. Aroca and P. Sanchez, IEEE Trans Mag,
 * vol 34, May 1998, p601-607.  (Side note: The units in the latter
 * paper are apparently cgs.)  It appears likely that there is an error
 * in the latter paper (attempts to implement the formulae presented
 * there did not produce the proper symmetries), as well as in the older
 * paper, "Magnetostatic Interaction Fields for a Three-Dimensional
 * Array of Ferromagnetic Cubes," Manfred E. Schabes and Amikam Aharoni,
 * IEEE Trans Mag, vol 23, November 1987, p3882-3888.  (Note: The Newell
 * paper deals with uniformly sized rectangular prisms, the Maicas paper
 * allows non-uniformly sized rectangular prisms, and the Schabes paper
 * only considers cubes.)
 *
 *   The kernel here is based on an analytically derived energy, and the
 * effective (discrete) demag field is calculated from the (discrete)
 * energy.
 *
 * Check values: Below are several values for Nxx and Nxy, calculated
 *  using Maple 7.0 with 100 decimal digits precision, rounded to 50
 *  digits for display.  Normal double precision IEEE floating point
 *  provides approximately 16 decimal digits of accuracy, and so-called
 *  "quadruple" precision provides about 34 decimal digits.
 *
 *
 *  x  y  z  dx dy dz |  Nxx(x,y,z,dx,dy,dz)
 * -------------------+-------------------------------------------------------
 *  0  0  0  50 10  1 |  0.021829576458713811627717362556500594396802771830582
 *  0  0  0   1  1  1 |  0.33333333333333333333333333333333333333333333333333
 *  0  0  0   1  1  2 |  0.40084192360558096752690050789034014000452668298259
 *  0  0  0   2  1  1 |  0.19831615278883806494619898421931971999094663403481
 *  0  0  0   1  2  3 |  0.53879030592371444784959040590642585972177691027128
 *  0  0  0   2  1  3 |  0.27839171603589255540904462483920117770716945564364
 *  0  0  0   3  2  1 |  0.18281797804039299674136496925437296257105363408509
 *  1  0  0   1  1  1 | -0.13501718054449526838713434911401361334238669929852
 *  0  1  0   1  1  1 |  0.067508590272247634193567174557006806671193349649259
 *  1  1  0   1  2  3 | -0.083703755298020462084677435631518487669613050161980
 *  1  1  1   1  1  1 |  0
 *  1  1  1   1  2  3 | -0.056075776617493854142226134670956166201885395511639
 *  1  2  3   1  2  3 |  0.0074263570277919738841364667668809954237071479183522
 * 10  1  1   1  2  3 | -0.00085675752896240944969766580856030369571736381174932
 * 10  4  6   1  2  3 | -0.00025381260722622800624859078080421562302790329870984
 *  3  6  9   1  2  3 |  0.00027042781311956323573639739731014558846864626288393
 *  6  6  6   1  2  3 | -0.000017252712652486259473939673630209925239022411008957
 *
 *  x  y  z  dx dy dz |  Nxy(x,y,z,dx,dy,dz)
 * -------------------+-------------------------------------------------------
 *  0  0  0  50 10  1 |  0
 *  0  0  0   1  1  1 |  0
 *  0  0  0   1  1  2 |  0
 *  0  0  0   2  1  1 |  0
 *  1  0  0   1  1  1 |  0
 *  0  1  0   1  1  1 |  0
 *  1  1  0   1  2  3 | -0.077258075615212400146921495217230818857603260305042
 *  1  1  1   1  1  1 | -0.016062127810508233979724830686189874772059681376565
 *  1  1  1   1  2  3 | -0.060966146490263272608967587158170018091418469887162
 *  1  2  3   1  2  3 | -0.0088226536707711039322880437795490754270432698781039
 * 10  1  1   1  2  3 | -0.00012776400247172360221892601892504762520639604467656
 * 10  4  6   1  2  3 | -0.00020004764005741154294387738750612412053502841766241
 *  3  6  9   1  2  3 | -0.00015720240165711869024193166157368874130207143569916
 *  6  6  6   1  2  3 | -0.00043908646098482886546108269881031774163900540796564
 *
 */

////////////////////////////////////////////////////////////////////////////
// Routines to do accurate summation

static int AS_Compare(const void* px,const void* py)
{
  // Comparison based on absolute values
  double x=fabs(*((const double *)px));
  double y=fabs(*((const double *)py));
  if(x<y) return 1;
  if(x>y) return -1;
  return 0;
}


static double
AccurateSum(int n,double *arr)
{
  // Order by decreasing magnitude
  qsort(arr,n,sizeof(double),AS_Compare);

  // Add up using doubly compensated summation.  If necessary, mark
  // variables these "volatile" to protect against problems arising
  // from extra precision.  Also, don't expect the compiler to respect
  // order with respect to parentheses at high levels of optimization,
  // i.e., write "u=x; u-=(y-corr)" as opposed to "u=x-(y-corr)".

  double sum,corr,y,u,t,v,z,x,tmp;

  sum=arr[0]; corr=0;
  for(int i=1;i<n;i++) {
    x=arr[i];
    y=corr+x;
    tmp=y-corr;
    u=x-tmp;
    t=y+sum;
    tmp=t-sum;
    v=y-tmp;
    z=u+v;
    sum=t+z;
    tmp=sum-t;
    corr=z-tmp;
  }
  return sum;
}


////////////////////////////////////////////////////////////////////////////
// Routines to calculate kernel coefficients
// See Newell et al. for details. The code below follows the
// naming conventions in that paper.

double SelfDemagNx(double x,double y,double z)
{ // Here Hx = -Nxx.Mx (formula (16) in Newell).
  // Note: egcs-2.91.57 on Linux/x86 with -O1 mangles this
  //  function (produces NaN's) unless we manually group terms.

  if(x<=0.0 || y<=0.0 || z<=0.0) return 0.0;
  if(x==y && y==z) return 1./3.;  // Special case: cube

  double xsq=x*x,ysq=y*y,zsq=z*z;
  double diag=sqrt(xsq+ysq+zsq);
  double arr[15];

  double mpxy = (x-y)*(x+y);
  double mpxz = (x-z)*(x+z);

  arr[0] = -4*(2*xsq*x-ysq*y-zsq*z);
  arr[1] =  4*(xsq+mpxy)*sqrt(xsq+ysq);
  arr[2] =  4*(xsq+mpxz)*sqrt(xsq+zsq);
  arr[3] = -4*(ysq+zsq)*sqrt(ysq+zsq);
  arr[4] = -4*diag*(mpxy+mpxz);

  arr[5] = 24*x*y*z*atan(y*z/(x*diag));
  arr[6] = 12*(z+y)*xsq*log(x);

  arr[7] = 12*z*ysq*log((sqrt(ysq+zsq)+z)/y);
  arr[8] = -12*z*xsq*log(sqrt(xsq+zsq)+z);
  arr[9] = 12*z*mpxy*log(diag+z);
  arr[10] = -6*z*mpxy*log(xsq+ysq);

  arr[11] =  12*y*zsq*log((sqrt(ysq+zsq)+y)/z);
  arr[12] = -12*y*xsq*log(sqrt(xsq+ysq)+y);
  arr[13] =  12*y*mpxz*log(diag+y);
  arr[14] =  -6*y*mpxz*log(xsq+zsq);

  double Nxx = AccurateSum(15,arr)/(12*WIDE_PI*x*y*z);
  return Nxx;
}

double SelfDemagNy(double xsize,double ysize,double zsize)
{ return SelfDemagNx(ysize,zsize,xsize); }

double SelfDemagNz(double xsize,double ysize,double zsize)
{ return SelfDemagNx(zsize,xsize,ysize); }


double
Newell_f(double x,double y,double z)
{ // There is mucking around here to handle case where imports
  // are near zero.  In particular, asinh(t) is written as
  // log(t+sqrt(1+t)) because the latter appears easier to
  // handle if t=y/x (for example) as x -> 0.

 // This function is even; the fabs()'s just simplify special case handling.
  x=fabs(x); double xsq=x*x;
  y=fabs(y); double ysq=y*y;
  z=fabs(z); double zsq=z*z;

  double R=xsq+ysq+zsq;
  if(R<=0.0) return 0.0;
  else       R=sqrt(R);

  // f(x,y,z)
  double piece[8];
  int piececount=0;
  if(z>0.) { // For 2D grids, half the calls from F1 have z==0.
    double temp1,temp2,temp3;
    piece[piececount++] = 2*(2*xsq-ysq-zsq)*R;
    if((temp1=x*y*z)>0.)
      piece[piececount++] = -12*temp1*atan2(y*z,x*R);
    if(y>0. && (temp2=xsq+zsq)>0.) {
      double dummy = log(((y+R)*(y+R))/temp2);
      piece[piececount++] = 3*y*zsq*dummy;
      piece[piececount++] = -3*y*xsq*dummy;
    }
    if((temp3=xsq+ysq)>0.) {
      double dummy = log(((z+R)*(z+R))/temp3);
      piece[piececount++] = 3*z*ysq*dummy;
      piece[piececount++] = -3*z*xsq*dummy;
    }
  } else {
    // z==0
    if(x==y) {
      //const double K = 2*sqrt(static_cast<double>(2.0))
      //  -6*log(1+sqrt(static_cast<double>(2.0)));
      double K = -2.4598143973710680537922785014593576970294;
      piece[piececount++] = K*xsq*x;
    } else {
      piece[piececount++] = 2*(2*xsq-ysq)*R;
      if(y>0. && x>0.)
	piece[piececount++] = -6*y*xsq*log((y+R)/x);
    }
  }

  return AccurateSum(piececount,piece)/12.;
}

double
CalculateSDA00(double x,double y,double z,
	       double dx,double dy,double dz)
{ // This is Nxx in Newell's paper
  double result=0.;
  if(x==0. && y==0. && z==0.) {
    // Self demag term.  The base routine can handle x==y==z==0,
    // but this should be more accurate.
    result = SelfDemagNx(dx,dy,dz);
  } else {
    // Simplified (collapsed) formula based on Newell's paper.
    // This saves about half the calls to f().  There is still
    // quite a bit of redundancy from one cell site to the next,
    // but as this is an initialization-only issue speed shouldn't
    // be too critical.
    double arr[27];
    arr[ 0] = -1*Newell_f(x+dx,y+dy,z+dz);
    arr[ 1] = -1*Newell_f(x+dx,y-dy,z+dz);
    arr[ 2] = -1*Newell_f(x+dx,y-dy,z-dz);
    arr[ 3] = -1*Newell_f(x+dx,y+dy,z-dz);
    arr[ 4] = -1*Newell_f(x-dx,y+dy,z-dz);
    arr[ 5] = -1*Newell_f(x-dx,y+dy,z+dz);
    arr[ 6] = -1*Newell_f(x-dx,y-dy,z+dz);
    arr[ 7] = -1*Newell_f(x-dx,y-dy,z-dz);

    arr[ 8] =  2*Newell_f(x,y-dy,z-dz);
    arr[ 9] =  2*Newell_f(x,y-dy,z+dz);
    arr[10] =  2*Newell_f(x,y+dy,z+dz);
    arr[11] =  2*Newell_f(x,y+dy,z-dz);
    arr[12] =  2*Newell_f(x+dx,y+dy,z);
    arr[13] =  2*Newell_f(x+dx,y,z+dz);
    arr[14] =  2*Newell_f(x+dx,y,z-dz);
    arr[15] =  2*Newell_f(x+dx,y-dy,z);
    arr[16] =  2*Newell_f(x-dx,y-dy,z);
    arr[17] =  2*Newell_f(x-dx,y,z+dz);
    arr[18] =  2*Newell_f(x-dx,y,z-dz);
    arr[19] =  2*Newell_f(x-dx,y+dy,z);

    arr[20] = -4*Newell_f(x,y-dy,z);
    arr[21] = -4*Newell_f(x,y+dy,z);
    arr[22] = -4*Newell_f(x,y,z-dz);
    arr[23] = -4*Newell_f(x,y,z+dz);
    arr[24] = -4*Newell_f(x+dx,y,z);
    arr[25] = -4*Newell_f(x-dx,y,z);

    arr[26] =  8*Newell_f(x,y,z);

    result = AccurateSum(27,arr)/(4*WIDE_PI*dx*dy*dz);
  }
  return result;
}


double
Newell_g(double x,double y,double z)
{ // There is mucking around here to handle case where imports
  // are near zero.  In particular, asinh(t) is written as
  // log(t+sqrt(1+t)) because the latter appears easier to
  // handle if t=y/x (for example) as x -> 0.

  double result_sign=1.0;
  if(x<0.0) result_sign *= -1.0;  if(y<0.0) result_sign *= -1.0;
  x=fabs(x); y=fabs(y); z=fabs(z);  // This function is even in z and
  /// odd in x and y.  The fabs()'s simplify special case handling.

  double xsq=x*x,ysq=y*y,zsq=z*z;
  double R=xsq+ysq+zsq;
  if(R<=0.0) return 0.0;
  else       R=sqrt(R);

  // g(x,y,z)
  double piece[7];
  int piececount=0;
  piece[piececount++] = -2*x*y*R;;
  if(z>0.) { // For 2D grids, 1/3 of the calls from CalculateSDA01 have z==0.
    piece[piececount++] = -z*zsq*atan2(x*y,z*R);
    piece[piececount++] = -3*z*ysq*atan2(x*z,y*R);
    piece[piececount++] = -3*z*xsq*atan2(y*z,x*R);

    double temp1,temp2,temp3;
    if((temp1=xsq+ysq)>0.)
      piece[piececount++] = 6*x*y*z*log((z+R)/sqrt(temp1));

    if((temp2=ysq+zsq)>0.)
      piece[piececount++] = y*(3*zsq-ysq)*log((x+R)/sqrt(temp2));

    if((temp3=xsq+zsq)>0.)
      piece[piececount++] = x*(3*zsq-xsq)*log((y+R)/sqrt(temp3));

  } else {
    // z==0.
    if(y>0.) piece[piececount++] = -y*ysq*log((x+R)/y);
    if(x>0.) piece[piececount++] = -x*xsq*log((y+R)/x);
  }

  return result_sign*AccurateSum(piececount,piece)/6.;
}

double
CalculateSDA01(double x,double y,double z,
	       double l,double h,double e)
{ // This is Nxy*(4*PI*tau) in Newell's paper.

  // Simplified (collapsed) formula based on Newell's paper.
  // This saves about half the calls to g().  There is still
  // quite a bit of redundancy from one cell site to the next,
  // but as this is an initialization-only issue speed shouldn't
  // be too critical.
  double arr[27];

  arr[ 0] = -1*Newell_g(x-l,y-h,z-e);
  arr[ 1] = -1*Newell_g(x-l,y-h,z+e);
  arr[ 2] = -1*Newell_g(x+l,y-h,z+e);
  arr[ 3] = -1*Newell_g(x+l,y-h,z-e);
  arr[ 4] = -1*Newell_g(x+l,y+h,z-e);
  arr[ 5] = -1*Newell_g(x+l,y+h,z+e);
  arr[ 6] = -1*Newell_g(x-l,y+h,z+e);
  arr[ 7] = -1*Newell_g(x-l,y+h,z-e);

  arr[ 8] =  2*Newell_g(x,y+h,z-e);
  arr[ 9] =  2*Newell_g(x,y+h,z+e);
  arr[10] =  2*Newell_g(x,y-h,z+e);
  arr[11] =  2*Newell_g(x,y-h,z-e);
  arr[12] =  2*Newell_g(x-l,y-h,z);
  arr[13] =  2*Newell_g(x-l,y+h,z);
  arr[14] =  2*Newell_g(x-l,y,z-e);
  arr[15] =  2*Newell_g(x-l,y,z+e);
  arr[16] =  2*Newell_g(x+l,y,z+e);
  arr[17] =  2*Newell_g(x+l,y,z-e);
  arr[18] =  2*Newell_g(x+l,y-h,z);
  arr[19] =  2*Newell_g(x+l,y+h,z);

  arr[20] = -4*Newell_g(x-l,y,z);
  arr[21] = -4*Newell_g(x+l,y,z);
  arr[22] = -4*Newell_g(x,y,z+e);
  arr[23] = -4*Newell_g(x,y,z-e);
  arr[24] = -4*Newell_g(x,y-h,z);
  arr[25] = -4*Newell_g(x,y+h,z);

  arr[26] =  8*Newell_g(x,y,z);

  return AccurateSum(27,arr)/(4*WIDE_PI*l*h*e);
  /// factor Nxy from Newell's paper.
}

double DemagNxxAsymptotic(double x, double y, double z,
                            double dx,double dy,double dz)
{ // Asymptotic approximation to Nxx term.
  double R2 = x*x + y*y + z*z;

  if(R2<=0.0) {
    // Asymptotic expansion doesn't apply for R==0.  Fall back
    // to self-demag calculation.
    return SelfDemagNx(dx,dy,dz);
  }

  double R  = sqrt(R2);

  double sx2 = x*x/R2;
  double sy2 = y*y/R2;
  double sz2 = z*z/R2;

  double sx4 = sx2*sx2;
  double sy4 = sy2*sy2;
  double sz4 = sz2*sz2;

  double sx6 = sx4*sx2;
  double sy6 = sy4*sy2;
  double sz6 = sz4*sz2;

  double dx2 = dx*dx;
  double dy2 = dy*dy;
  double dz2 = dz*dz;

  double dx4 = dx2*dx2;
  double dy4 = dy2*dy2;
  double dz4 = dz2*dz2;

  double term3 = 2*sx2 - sy2 - sz2;

  double term5 = 0.0;
  if(dx2!=dy2 || dx2!=dz2 || dy2!=dz2) { // Non-cube case
    double a1 =   8*dx2  -  4*dy2  -  4*dz2;
    double a2 = -24*dx2  + 27*dy2  -  3*dz2;
    double a3 = -24*dx2  -  3*dy2  + 27*dz2;
    double a4 =   3*dx2  -  4*dy2  +  1*dz2;
    double a5 =   6*dx2  -  3*dy2  -  3*dz2;
    double a6 =   3*dx2  +  1*dy2  -  4*dz2;
    term5 = a1*sx4 + a2*sx2*sy2 + a3*sx2*sz2
      + a4*sy4 + a5*sy2*sz2 + a6*sz4;
    term5 *= 0.25;
  }

  double term7 = 0.0;
  {
    double b1  =   32*dx4  -  40*dx2*dy2  -  40*dx2*dz2   +  12*dy4   +  10*dy2*dz2   +  12*dz4;
    double b2  = -240*dx4  + 580*dx2*dy2  +  20*dx2*dz2   - 202*dy4   -  75*dy2*dz2   +  22*dz4;
    double b3  = -240*dx4  +  20*dx2*dy2  + 580*dx2*dz2   +  22*dy4   -  75*dy2*dz2   - 202*dz4;
    double b4  =  180*dx4  - 505*dx2*dy2  +  55*dx2*dz2   + 232*dy4   -  75*dy2*dz2   +   8*dz4;
    double b5  =  360*dx4  - 450*dx2*dy2  - 450*dx2*dz2   - 180*dy4   + 900*dy2*dz2   - 180*dz4;
    double b6  =  180*dx4  +  55*dx2*dy2  - 505*dx2*dz2   +   8*dy4   -  75*dy2*dz2   + 232*dz4;
    double b7  =  -10*dx4  +  30*dx2*dy2  -   5*dx2*dz2   -  16*dy4   +  10*dy2*dz2   -   2*dz4;
    double b8  =  -30*dx4  +  55*dx2*dy2  +  20*dx2*dz2   +   8*dy4   -  75*dy2*dz2   +  22*dz4;
    double b9  =  -30*dx4  +  20*dx2*dy2  +  55*dx2*dz2   +  22*dy4   -  75*dy2*dz2   +   8*dz4;
    double b10 =  -10*dx4  -   5*dx2*dy2  +  30*dx2*dz2   -   2*dy4   +  10*dy2*dz2   -  16*dz4;

    term7 = b1*sx6 + b2*sx4*sy2 + b3*sx4*sz2 + b4*sx2*sy4 + b5*sx2*sy2*sz2
      + b6*sx2*sz4 + b7*sy6 + b8*sy4*sz2 + b9*sy2*sz4 + b10*sz6;
    term7 *= (1.0/16.0);
  }

  double Nxx = (-dx*dy*dz/(4*WIDE_PI)) * (((term7/R2 + term5)/R2 + term3)/(R2*R));
  // Error should be of order 1/R^9

  return Nxx;
}

double DemagNxyAsymptotic(double x, double y, double z,
                            double dx,double dy,double dz)
{ // Asympotic approximation to Nxy term.
  double R2 = x*x + y*y + z*z;

  if(R2<=0.0) {
    // Asymptotic expansion doesn't apply for R==0.  Fall back
    // to self-demag calculation.
    return (0.0);
  }

  double R  = sqrt(R2);

  double sx2 = x*x/R2;
  double sy2 = y*y/R2;
  double sz2 = z*z/R2;

  double sx4 = sx2*sx2;
  double sy4 = sy2*sy2;
  double sz4 = sz2*sz2;

  double dx2 = dx*dx;
  double dy2 = dy*dy;
  double dz2 = dz*dz;

  double dx4 = dx2*dx2;
  double dy4 = dy2*dy2;
  double dz4 = dz2*dz2;

  double term3 = 3;

  double term5 = 0.0;
  if(dx2!=dy2 || dx2!=dz2 || dy2!=dz2) { // Non-cube case
    double a1 =   4*dx2  -  3*dy2  -  1*dz2;
    double a2 =  -3*dx2  +  4*dy2  -  1*dz2;
    double a3 =  -3*dx2  -  3*dy2  +  6*dz2;
    term5 = a1*sx2 + a2*sy2 + a3*sz2;
    term5 *= (5.0/4.0);
  }

  double term7 = 0.0;
  {
    double b1  =   16*dx4  -  30*dx2*dy2  -  10*dx2*dz2   +  10*dy4   +   5*dy2*dz2   +   2*dz4;
    double b2  =  -40*dx4  + 105*dx2*dy2  -   5*dx2*dz2   -  40*dy4   -   5*dy2*dz2   +   4*dz4;
    double b3  =  -40*dx4  -  15*dx2*dy2  + 115*dx2*dz2   +  20*dy4   -  35*dy2*dz2   -  32*dz4;
    double b4  =   10*dx4  -  30*dx2*dy2  +   5*dx2*dz2   +  16*dy4   -  10*dy2*dz2   +   2*dz4;
    double b5  =   20*dx4  -  15*dx2*dy2  -  35*dx2*dz2   -  40*dy4   + 115*dy2*dz2   -  32*dz4;
    double b6  =   10*dx4  +  15*dx2*dy2  -  40*dx2*dz2   +  10*dy4   -  40*dy2*dz2   +  32*dz4;
    term7 = b1*sx4 + b2*sx2*sy2 + b3*sx2*sz2
      + b4*sy4 + b5*sy2*sz2 + b6*sz4;
    term7 *= (7.0/16.0);
  }

  double Nxy = (-dx*dy*dz*x*y/(4*WIDE_PI*R2)) * (((term7/R2 + term5)/R2 + term3)/(R2*R));
  // Error should be of order 1/R^9

  return Nxy;
}


double DemagNyyAsymptotic(double x, double y, double z,
                            double dx,double dy,double dz)
{
  return DemagNxxAsymptotic(y,x,z,dy,dx,dz);
}

double DemagNzzAsymptotic(double x, double y, double z,
                            double dx,double dy,double dz)
{
  return DemagNxxAsymptotic(z,y,x,dz,dy,dx);
}

double DemagNxzAsymptotic(double x, double y, double z,
                            double dx,double dy,double dz)
{
  return DemagNxyAsymptotic(x,z,y,dx,dz,dy);
}

double DemagNyzAsymptotic(double x, double y, double z,
                            double dx,double dy,double dz)
{
  return DemagNxyAsymptotic(y,z,x,dy,dz,dx);
}
