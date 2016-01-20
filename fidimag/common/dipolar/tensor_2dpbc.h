

//------------------------------------------------------------------------------
inline double RR(double x,double y,double z)
{return x*x+y*y+z*z;}

inline double R(double x,double y,double z)
{return sqrt(x*x+y*y+z*z);}

//------------------------------------------------------------------------------	       

//------------------------------------------------------------------------------
//about Nxxinf

inline double 
fxx(double x,double y,double z)
{
	double t1,t2,t3,t5;
	t1 = x*x;
	t2 = y*y;
	t3 = z*z;
	t5 = sqrt(t1 + t2 + t3);
	return y * x / (t1 + t3) / t5;
}


inline double Nxxinf(double x,double y,double z,double X0,double Y0) {
	return fxx(x+X0,y-Y0,z)+fxx(x-X0,y+Y0,z)-fxx(x+X0,y+Y0,z)-fxx(x-X0,y-Y0,z);
}

//------------------------------------------------------------------------------
//about Nyyinf

inline double Nyyinf(double x,double y,double z,double X0,double Y0){
	return Nxxinf(y,x,z,Y0,X0);
}


//------------------------------------------------------------------------------
//about Nzzinf

inline double fzz(double x,double y,double z){
	double t1,t2,t3,t5;
	t1 = x*x;
	t2 = y*y;
	t3 = z*z;
	t5 = sqrt(t1 + t2 + t3);
	return  y*x*(t1+t2+2*t3)/((t1 + t3)*(t2+t3)*t5);
}


inline double 
Nzzinf(double x,double y,double z,double X0,double Y0)
{
	return fzz(x+X0,y+Y0,z)+fzz(x-X0,y-Y0,z)-fzz(x+X0,y-Y0,z)-fzz(x-X0,y+Y0,z);
}



//------------------------------------------------------------------------------
//about Nxyinf
inline double fxy(double x,double y,double z){
	return 1.0/R(x,y,z);
}


inline double Nxyinf(double x,double y,double z,double X0,double Y0){
	return (fxy(x-X0,y-Y0,z)+fxy(x+X0,y+Y0,z)-fxy(x-X0,y+Y0,z)-fxy(x+X0,y-Y0,z));
}


//------------------------------------------------------------------------------
//about Nxzinf
inline double fxz(double x,double y,double z){
	double t2,t3,t4,t6;
	t2=x*x;
	t3=y*y;
	t4=z*z;
	t6 = sqrt(t2 + t3 + t4);
	return y * z / (t6 * (t2 + t4));
}


inline double Nxzinf(double x,double y,double z,double X0,double Y0){
	return fxz(x+X0,y-Y0,z)+fxz(x-X0,y+Y0,z)-fxz(x+X0,y+Y0,z)-fxz(x-X0,y-Y0,z);
}

//------------------------------------------------------------------------------
//about Nyzinf

inline double Nyzinf(double x,double y,double z,double X0,double Y0){
	return Nxzinf(y,x,z,Y0,X0);
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//about Nxxdipole
inline double 
DemagNxxDipolar(double x,double y,double z)
{
double  t1 = x*x;
double  t3 = y*y;
double  t4 = z*z;
double  t6 = t1 + t3 + t4;;
double  t7 = t6*t6;
double  t8 = sqrt(t6);
return -(2.0 * t1 - t3 - t4) / (t8 * t7);
}

inline double 
DemagNxyDipolar(double x,double y,double z)
{
double  t6 = RR(x,y,z);
double  t7 = t6*t6;
double  t8 = sqrt(t6);
return -3.0*x*y/(t8*t7);
}
//-----------------------------------------------------------------------------
double DemagTensorNormal(enum Type_Nij comp,double x,double y,double z,double a,double b,double c);

double DemagTensorAsymptotic(enum Type_Nij comp,double x,double y,double z,double a,double b,double c);

double DemagTensorDipolar(enum Type_Nij comp,double x,double y,double z);

double DemagTensorInfinite(enum Type_Nij comp,double x,double y,double z,double X0,double Y0);

