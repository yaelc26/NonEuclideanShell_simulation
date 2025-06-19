/*
	SpecialFunctions.h

*/

#ifndef _SpecialFunctions_H_
#define _SpecialFunctions_H_	

#include "Main.H"


class SpecialFunctions 
	{
	public:
		
		/* The LogGamma function */
		static double LogGamma(double x)
		{
			if (x<10)
				return LogGamma(x+1) - log(x);
			else
				return (x-0.5)*log(x) - x + 0.5*log(2*Pi) + 1.0/(12.0*x) -
				1.0/(360*pow(x,3)) + 1.0/(1260*pow(x,5)) - 1.0/(1680*pow(x,7)); 
		}
		
		/* The Gamma function */
		static double Gamma(double x)   {return exp(LogGamma(x));}
		
	};


#endif
