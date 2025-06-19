/*
 *  Tensor2D.H
 *  RKLibrary
 *
 *  Created by Raz Kupferman on 8/22/09.
 *  Copyright 2009 The Hebrew University. All rights reserved.
 *
 * This class implements a 2-by-2 real-valued matrix.
 * The reason to design a specific class rather than using a Matrix<Real>
 * is to accelerate calculations and provide functions that are specific 
 * to tensors in 2D fluid mechanics. In particular, we implement functions
 * required by the log-C formations of viscoelasticity.
 */

#ifndef _Tensor2D_H
#define _Tensor2D_H


class Tensor2D
	{
	public:
		
		double m_xx, m_xy, m_yx, m_yy;
		static double s_JacobiTolerance;
		
	public:
		
		/* Default constructors */
		Tensor2D() : m_xx(0), m_xy(0), m_yx(0), m_yy(0) {}
		
		/* Standard constructor with values of components */
		Tensor2D(double a_xx, double a_xy, double a_yx, double a_yy) :
		m_xx(a_xx), m_xy(a_xy), m_yx(a_yx), m_yy(a_yy) 
		{}
		
		/* Copy constructor */
		Tensor2D(const Tensor2D &a_Tensor) : 
		m_xx(a_Tensor.m_xx), m_xy(a_Tensor.m_xy), m_yx(a_Tensor.m_yx), m_yy(a_Tensor.m_yy) 
		{}
		
		/* Static functions of unit and zero */
		static Tensor2D Zero() 
		{
			return Tensor2D(0,0,0,0);
		}
		
		static Tensor2D Unit() 
		{
			return Tensor2D(1,0,0,1);
		}
		
		/* Assignments */
		void operator =(const Tensor2D &a_Tensor)
		{	
			m_xx = a_Tensor.m_xx;
			m_xy = a_Tensor.m_xy;
			m_yx = a_Tensor.m_yx;
			m_yy = a_Tensor.m_yy;
		}			
		
		/* Additon assignment */
		void operator+=(const Tensor2D &a_Tensor)
		{	
			m_xx += a_Tensor.m_xx;
			m_xy += a_Tensor.m_xy;
			m_yx += a_Tensor.m_yx;
			m_yy += a_Tensor.m_yy;
		}			
		
		/* Subtraction assignment */
		void operator-=(const Tensor2D &a_Tensor)
		{	
			m_xx -= a_Tensor.m_xx;
			m_xy -= a_Tensor.m_xy;
			m_yx -= a_Tensor.m_yx;
			m_yy -= a_Tensor.m_yy;
		}			
		
		/* Multiply assignment: note that this is a matric multiply */
		void operator*=(const Tensor2D &a_Tensor)
		{
			double Rxx = m_xx*a_Tensor.m_xx + m_xy*a_Tensor.m_yx;
			double Rxy = m_xx*a_Tensor.m_xy + m_xy*a_Tensor.m_yy;
			double Ryx = m_yx*a_Tensor.m_xx + m_yy*a_Tensor.m_yx;
			double Ryy = m_yx*a_Tensor.m_xy + m_yy*a_Tensor.m_yy;
			m_xx = Rxx;
			m_xy = Rxy;
			m_yx = Ryx;
			m_yy = Ryy;
		}
		
		/* Scale by a scalar */
		void operator*=(double a_val)
		{
			m_xx *= a_val;
			m_xy *= a_val;
			m_yx *= a_val;
			m_yy *= a_val;
		}
		
		/* Divide by a scalar */
		void operator/=(double a_val)
		{
			m_xx /= a_val;
			m_xy /= a_val;
			m_yx /= a_val;
			m_yy /= a_val;
		}
		
		/* Trace */
		double trace() const 
		{
			return m_xx + m_yy;
		}
		
		/* Determinant */
		double det() const 
		{	
			return m_xx*m_yy - m_xy*m_yx;
		}
		
		/* streaming to an ostream */
		friend std::ostream& operator<<(std::ostream &a_os, const Tensor2D &a_Tensor)
		{
			a_os << a_Tensor.m_xx << " " << a_Tensor.m_xy << std::endl;
			a_os << a_Tensor.m_yx << " " << a_Tensor.m_yy << std::endl;
			return a_os;
		}
		
		
		/* Arithmetic functions */
		friend Tensor2D operator-(const Tensor2D &a_Tensor)
		{
			return Tensor2D(-a_Tensor.m_xx, 
							-a_Tensor.m_xy, 
							-a_Tensor.m_yx, 
							-a_Tensor.m_yy);
		}
		
		friend Tensor2D operator+(const Tensor2D &a_Tensor1, const Tensor2D &a_Tensor2)
		{
			Tensor2D ret(a_Tensor1);
			ret += a_Tensor2;
			return ret;
		}
		
		friend Tensor2D operator-(const Tensor2D &a_Tensor1, const Tensor2D &a_Tensor2)
		{
			Tensor2D ret(a_Tensor1);
			ret -= a_Tensor2;
			return ret;
		}
		
		
		friend Tensor2D operator*(const Tensor2D &a_Tensor1, const Tensor2D &a_Tensor2)
		{
			Tensor2D ret(a_Tensor1);
			ret *= a_Tensor2;
			return ret;
		}
		
		friend Tensor2D operator*(double a_val, const Tensor2D &a_Tensor)
		{
			return Tensor2D(a_val*a_Tensor.m_xx, 
							a_val*a_Tensor.m_xy, 
							a_val*a_Tensor.m_yx, 
							a_val*a_Tensor.m_yy);
		}
		
		friend Tensor2D operator*(const Tensor2D &a_Tensor, double a_val)
		{
			return Tensor2D(a_val*a_Tensor.m_xx, 
							a_val*a_Tensor.m_xy, 
							a_val*a_Tensor.m_yx, 
							a_val*a_Tensor.m_yy);
		}
		
		friend Tensor2D operator/(const Tensor2D &a_Tensor, double a_val)
		{
			return Tensor2D(a_Tensor.m_xx/a_val, 
							a_Tensor.m_xy/a_val, 
							a_Tensor.m_yx/a_val, 
							a_Tensor.m_yy/a_val);
		}
		
		
		/* Transpose */
		friend Tensor2D transpose(const Tensor2D &a_Tensor)
		{
			return Tensor2D(a_Tensor.m_xx, 
							a_Tensor.m_yx, 
							a_Tensor.m_xy, 
							a_Tensor.m_yy);
		}
		
		/* The square root of a matrix (assumes symmetry!!!) */
		friend Tensor2D sqrt(const Tensor2D &a_Tensor)
		{
			Tensor2D A = a_Tensor;
			Tensor2D J = Tensor2D::Unit();
			JacobiIterations(A,J);
			A = Tensor2D(sqrt(A.m_xx), 0,0,sqrt(A.m_yy));
			return J*A*transpose(J);
		}
		
		/* The matrix logarithm (assumes symmetry!!!) */ 
		friend Tensor2D log (const Tensor2D &a_Tensor)
		{
			Tensor2D A = a_Tensor;
			Tensor2D J = Tensor2D::Unit();
			JacobiIterations(A,J);
			A = Tensor2D(log(A.m_xx), 0,0,log(A.m_yy));
			return J*A*transpose(J);			
		}
		
		/* The matrix exponential (assumes symmetry!!!) */
		friend Tensor2D exp (const Tensor2D &a_Tensor)
		{
			Tensor2D A = a_Tensor;
			Tensor2D J = Tensor2D::Unit();
			JacobiIterations(A,J);
			A = Tensor2D(exp(A.m_xx), 0,0,exp(A.m_yy));
			return J*A*transpose(J);			
		}
		
		/* Helper function: Jacobi iterarions */
		friend void JacobiIterations(Tensor2D &a_Tensor, Tensor2D &a_J)
		{
			double JacobiTolerance = 1.e-12;
			while (fabs(a_Tensor.m_xy) > JacobiTolerance)
			{
				double tau = 0.5*(a_Tensor.m_xx - a_Tensor.m_yy)/a_Tensor.m_xy;
				double t   = 1.0/(fabs(tau) + sqrt(1+tau*tau));
				t = (tau>0) ? t : -t; 
				double c = 1/sqrt(1+t*t);
				double s = c*t;
				Tensor2D R(c,-s,s,c);
				a_Tensor = transpose(R)*a_Tensor*R;
				a_J = a_J*R;
			}
		}
		
	};


		 
#endif
		 
