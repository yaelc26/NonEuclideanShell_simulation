#ifndef main_h_
#define main_h_


#include <cassert>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>

/* Mathematical constants */
const double    Pi         = 3.14159265358979323846;
const double    E          = 2.71828182845904509;
const double    EulerGamma = 0.577215664901532843;

/* min and max */
#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

inline double minmod(double a, double b)
{
	if (a*b < 0) return 0;
	if (fabs(a)<fabs(b)) return a;
	return b;
}

#endif
