/*
 *  TinyMatrix.H
 *  RKLibrary
 *
 *  Created by Raz Kupferman on 9/9/09.
 *  Copyright 2009 The Hebrew University. All rights reserved.
 *
 */

/*
 This class implements a simple square matrix of fixed size, along with the 
 most basic arithmetics. It is a template class with parametrized
 type and size. There is no reference counting.
*/

#ifndef _TINYMATRIX_H_
#define _TINYMATRIX_H_

#include "Main.H"

template <class T, int SIZE> class TinyMatrix;
template <class T, int SIZE> TinyMatrix<T,SIZE> operator+(const TinyMatrix<T,SIZE>&, const TinyMatrix<T,SIZE>&);
template <class T, int SIZE> TinyMatrix<T,SIZE> operator-(const TinyMatrix<T,SIZE>&, const TinyMatrix<T,SIZE>&);
template <class T, int SIZE> TinyVector<T,SIZE> operator*(const TinyMatrix<T,SIZE>&, const TinyVector<T,SIZE>&);
template <class T, int SIZE> TinyMatrix<T,SIZE> operator*(const TinyMatrix<T,SIZE>&, const TinyMatrix<T,SIZE>&);
template <class T, int SIZE> void               luDecompose(TinyMatrix<T,SIZE>& a_lu, 
															TinyVector<int,SIZE>& a_pivot);
template <class T, int SIZE> TinyVector<T,SIZE> luSolve(const TinyMatrix<T,SIZE>& a_lu, 
														const TinyVector<int,SIZE>& a_pivot,
														const TinyVector<T,SIZE>& a_rhs);

template <class T, int SIZE> 
class TinyMatrix
	{
	public:
		
		/* friends functions */
		friend TinyMatrix<T,SIZE> operator+ <>(const TinyMatrix<T,SIZE>&, const TinyMatrix<T,SIZE>&);
		friend TinyMatrix<T,SIZE> operator- <>(const TinyMatrix<T,SIZE>&, const TinyMatrix<T,SIZE>&);
		friend TinyVector<T,SIZE> operator* <>(const TinyMatrix<T,SIZE>&, const TinyVector<T,SIZE>&);
		friend TinyMatrix<T,SIZE> operator* <>(const TinyMatrix<T,SIZE>&, const TinyMatrix<T,SIZE>&);

		/* Default constructor */
		TinyMatrix() 
		{
			for (int i=0; i<SIZE; i++) 
				for (int j=0; j<SIZE; j++) 
					m_data[i][j] = T();
		}
		
		/* Copy constructor */
		TinyMatrix(const TinyMatrix<T,SIZE>& a_rhs)
		{
			for (int i=0; i<SIZE; i++) 
				for (int j=0; j<SIZE; j++) 
					m_data[i][j] = a_rhs.m_data[i][j];
		}
		
		/* Assignment */
		void operator=(const TinyMatrix<T,SIZE>& a_rhs)
		{
			for (int i=0; i<SIZE; i++) 
				for (int j=0; j<SIZE; j++) 
					m_data[i][j] = a_rhs.m_data[i][j];
		}
		
		/* Increment/Decrement/Scale */
		void operator+=(const TinyMatrix<T,SIZE>& a_rhs)
		{
			for (int i=0; i<SIZE; i++) 
				for (int j=0; j<SIZE; j++) 
					m_data[i][j] += a_rhs.m_data[i][j];
		}

		void operator-=(const TinyMatrix<T,SIZE>& a_rhs)
		{
			for (int i=0; i<SIZE; i++) 
				for (int j=0; j<SIZE; j++) 
					m_data[i][j] -= a_rhs.m_data[i][j];
		}

		void scale(T a_value)
		{
			for (int i=0; i<SIZE; i++) 
				for (int j=0; j<SIZE; j++) 
					m_data[i][j] *= a_value;
		}
		
		/* Set value to zero */
		void setToZero()
		{
			for (int i=0; i<SIZE; i++) 
				for (int j=0; j<SIZE; j++) 
					m_data[i][j] = 0;
		}
		
		/* Determinant */
		T det() const
		{
			if (SIZE==2)
			{
				return m_data[0][0]*m_data[1][1] - m_data[0][1]*m_data[1][0];
			}
			else if (SIZE==3)
			{
				return  
				m_data[0][0]*m_data[1][1]*m_data[2][2] - 
				m_data[0][0]*m_data[1][2]*m_data[2][1] + 
				m_data[0][1]*m_data[1][2]*m_data[2][0] - 
				m_data[0][1]*m_data[1][0]*m_data[2][2] + 
				m_data[0][2]*m_data[1][0]*m_data[2][1] - 
				m_data[0][2]*m_data[1][1]*m_data[2][0];
			}
		}
		
		/* Trace */
		T trace() const 
		{
			T ret=0;
			for (int i=0; i<SIZE; i++) ret += m_data[i][i];
			return ret;			
		}
		
		/* Transpose */
		TinyMatrix<T,SIZE> transpose() const
		{
			TinyMatrix<T,SIZE> ret;
			for (int i=0; i<SIZE; i++) 
				for (int j=0; j<SIZE; j++) 
					ret.m_data[i][j] = m_data[j][i];
			return ret;
		}
		
		/* Inverse matrix */
		TinyMatrix<T,SIZE> inverse() const
		{
			TinyMatrix<T,SIZE> ret;
			if (SIZE==2)
			{
				double               d = this->det();
				ret.m_data[0][0] =  m_data[1][1];
				ret.m_data[0][1] = -m_data[1][0];
				ret.m_data[1][0] = -m_data[0][1];
				ret.m_data[1][1] =  m_data[0][0];
				ret.scale(1.0/d);
				return ret;
			}
			else if (SIZE==3)
			{
				double               d = this->det();
				ret.m_data[0][0] = (m_data[1][1]*m_data[2][2] - m_data[1][2]*m_data[2][1]);
				ret.m_data[0][1] = (m_data[2][1]*m_data[0][2] - m_data[2][2]*m_data[0][1]);
				ret.m_data[0][2] = (m_data[0][1]*m_data[1][2] - m_data[0][2]*m_data[1][1]);
				ret.m_data[1][0] = (m_data[1][2]*m_data[2][0] - m_data[1][0]*m_data[2][2]);
				ret.m_data[1][1] = (m_data[2][2]*m_data[0][0] - m_data[2][0]*m_data[0][2]);
				ret.m_data[1][2] = (m_data[0][2]*m_data[1][0] - m_data[0][0]*m_data[1][2]);
				ret.m_data[2][0] = (m_data[1][0]*m_data[2][1] - m_data[1][1]*m_data[2][0]);
				ret.m_data[2][1] = (m_data[2][0]*m_data[0][1] - m_data[2][1]*m_data[0][0]);
				ret.m_data[2][2] = (m_data[0][0]*m_data[1][1] - m_data[0][1]*m_data[1][0]);
				ret.scale(1.0/d);
				return ret;
			}
		}
				
		/* array-like access to elements */
		T  operator()(int a_index1, int a_index2) const {return m_data[a_index1][a_index2];}
		T& operator()(int a_index1, int a_index2)       {return m_data[a_index1][a_index2];}
		
		/* access to pointer */
		T*       getPointer()       {return m_data;}
		const T* getPointer() const {return m_data;}
		
	private:
		
		T m_data[SIZE][SIZE];
	};


/* Functions */
template <class T, int SIZE> 
void luDecompose(TinyMatrix<T,SIZE>& a_lu, TinyVector<int,SIZE>& a_pivot) 
{	
	const double TINY=1.0e-40;
	int          i,imax,j,k;
	double       big,temp;
	double       d = 1.0;
	TinyVector<T,SIZE> vv;
	
	for (i=0;i<SIZE;i++) 
	{
		big=0.0;
		for (j=0;j<SIZE;j++)
			if ((temp=fabs(a_lu(i,j))) > big) big=temp;
		if (big == 0.0) std::cout << "Singular matrix in luDecompose" << std::endl;
		vv(i)=1.0/big;
	}
	for (k=0;k<SIZE;k++) {
		big=0.0;
		for (i=k;i<SIZE;i++) {
			temp=vv(i)*fabs(a_lu(i,k));
			if (temp > big) {
				big  = temp;
				imax = i;
			}
		}
		if (k != imax) 
		{
			for (j=0;j<SIZE;j++) 
			{
				temp         = a_lu(imax,j);
				a_lu(imax,j) = a_lu(k,j);
				a_lu(k,j)    = temp;
			}
			d = -d;
			vv(imax) = vv(k);
		}
		a_pivot(k) = imax;
		if (a_lu(k,k) == 0.0) a_lu(k,k)=TINY;
		for (i=k+1;i<SIZE;i++) 
		{
			a_lu(i,k) = a_lu(i,k) / a_lu(k,k);
			temp = a_lu(i,k);
			for (j=k+1;j<SIZE;j++)  a_lu(i,j) -= temp*a_lu(k,j);
		}
	}
}

template <class T, int SIZE> 
TinyVector<T,SIZE> luSolve(const TinyMatrix<T,SIZE>& a_lu, 
						   const TinyVector<int,SIZE>& a_pivot,
						   const TinyVector<T,SIZE>& a_rhs)
{
	TinyVector<T,SIZE> x;
	
	int i,ii=0,ip,j;
	double sum;
	for (i=0;i<SIZE;i++) x(i) = a_rhs(i);
	for (i=0;i<SIZE;i++) 
	{
		ip    = a_pivot(i);
		sum   = x(ip);
		x(ip) = x(i);
		if (ii != 0)
			for (j=ii-1;j<i;j++) 
				sum -= a_lu(i,j)*x(j);
		else if (sum != 0.0) ii=i+1;
		x(i)=sum;
	}
	for (i=SIZE-1;i>=0;i--) 
	{
		sum = x(i);
		for (j=i+1;j<SIZE;j++) sum -= a_lu(i,j)*x(j);
		x(i)=sum/a_lu(i,i);
	}	
	return x;
}


/* Friend functions */

/* Addition */
template <class T, int SIZE> 
inline TinyMatrix<T,SIZE> operator+(const TinyMatrix<T,SIZE>& a_mat1, const TinyMatrix<T,SIZE>& a_mat2)
{
	TinyMatrix<T,SIZE> ret;
	for (int i=0; i<SIZE; i++) 
		for (int j=0; j<SIZE; j++) 
			ret.m_data[i][j] = a_mat1.m_data[i][j] + a_mat2.m_data[i][j];
	return ret;
}

/* Subtraction */
template <class T, int SIZE> 
inline TinyMatrix<T,SIZE> operator-(const TinyMatrix<T,SIZE>& a_mat1, const TinyMatrix<T,SIZE>& a_mat2)
{
	TinyMatrix<T,SIZE> ret;
	for (int i=0; i<SIZE; i++) 
		for (int j=0; j<SIZE; j++) 
			ret.m_data[i][j] = a_mat1.m_data[i][j] - a_mat2.m_data[i][j];
	return ret;
}

/* Vector multiplication */
template <class T, int SIZE> 
inline TinyVector<T,SIZE> operator*(const TinyMatrix<T,SIZE>& a_mat, const TinyVector<T,SIZE>& a_vec)
{
	TinyVector<T,SIZE> ret;
	for (int i=0; i<SIZE; i++) 
		for (int j=0; j<SIZE; j++) 
			ret.m_data[i] += a_mat.m_data[i][j]*a_vec.m_data[j];
	return ret;
}

/* Matrix multiplication */
template <class T, int SIZE> 
inline TinyMatrix<T,SIZE> operator*(const TinyMatrix<T,SIZE>& a_mat1, const TinyMatrix<T,SIZE>& a_mat2)
{
	TinyMatrix<T,SIZE> ret;
	for (int i=0; i<SIZE; i++) 
		for (int j=0; j<SIZE; j++) 
			for (int k=0; k<SIZE; k++)
					ret.m_data[i][j] += a_mat1.m_data[i][k]*a_mat2.m_data[k][j];
	return ret;
}

/* Print the matrix to stream */
template <class T, int SIZE> 
std::ostream& operator<<(std::ostream &a_os, const TinyMatrix<T, SIZE>& a_mat) 
{
	a_os << "[";
	for (int i=0; i<SIZE; i++) 
	{
		for (int j=0; j<SIZE; j++)
			a_os << a_mat(i,j) << " ";
		a_os << std::endl;
	}
	a_os << "]" << std::endl;
	return a_os;
}
// Scalar multiplication (right multiplication)
template <class T, int SIZE>
inline TinyMatrix<T,SIZE> operator*(const TinyMatrix<T,SIZE>& a_mat, double scalar)
{
    TinyMatrix<T,SIZE> ret;
    for (int i=0; i<SIZE; i++)
        for (int j=0; j<SIZE; j++)
            ret(i,j) = a_mat(i,j) * scalar;
    return ret;
}

// Scalar multiplication (left multiplication)
template <class T, int SIZE>
inline TinyMatrix<T,SIZE> operator*(double scalar, const TinyMatrix<T,SIZE>& a_mat)
{
    return a_mat * scalar;
}

#endif
