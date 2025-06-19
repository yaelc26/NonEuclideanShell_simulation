/*
 *  TinyVector.H
 *  RKLibrary
 *
 *  Created by Raz Kupferman on 9/9/09.
 *  Copyright 2009 The Hebrew University. All rights reserved.
 *
 */

/*
 This class implements a simple vector of fixed size, along with the 
 most basic vector arithmetics. It is a template class with parametrized
 type and size. There is no reference counting.
*/

#ifndef _TINYVECTOR_H_
#define _TINYVECTOR_H_

#include "Main.H"

template <class T, int SIZE> class TinyVector;
template <class T, int SIZE> TinyVector<T,SIZE> operator+(const TinyVector<T,SIZE>&, const TinyVector<T,SIZE>&);
template <class T, int SIZE> TinyVector<T,SIZE> operator-(const TinyVector<T,SIZE>&, const TinyVector<T,SIZE>&);
template <class T, int SIZE> T                  innerProduct(const TinyVector<T,SIZE>&, const TinyVector<T,SIZE>&);
template <class T, int SIZE> T                  norm(const TinyVector<T,SIZE>&);
template <class T, int SIZE> TinyVector<T,SIZE> CrossProduct(const TinyVector<T,SIZE>&, const TinyVector<T,SIZE>&); 

template <class T, int SIZE> class TinyMatrix;
template <class T, int SIZE> TinyVector<T,SIZE> operator*(const TinyMatrix<T,SIZE>&, const TinyVector<T,SIZE>&);

template <class T, int SIZE> 
class TinyVector
	{
	public:
		
		/* friends functions */
		friend TinyVector<T,SIZE> operator+    <>(const TinyVector<T,SIZE>&, const TinyVector<T,SIZE>&);
		friend TinyVector<T,SIZE> operator-    <>(const TinyVector<T,SIZE>&, const TinyVector<T,SIZE>&);
		friend TinyVector<T,SIZE> operator*    <>(const TinyMatrix<T,SIZE>&, const TinyVector<T,SIZE>&);
		friend T                  innerProduct <>(const TinyVector<T,SIZE>&, const TinyVector<T,SIZE>&);
		friend TinyVector<T,SIZE> CrossProduct <>(const TinyVector<T,SIZE>&, const TinyVector<T,SIZE>&);


		/* Default constructor */
		TinyVector() 
		{
			for (int i=0; i<SIZE; i++) m_data[i]=0;
		}
		
		/* Copy constructor */
		TinyVector(const TinyVector<T,SIZE>& a_rhs)
		{
			for (int i=0; i<SIZE; i++) m_data[i] = a_rhs.m_data[i];
		}
		
		/* Assignment */
		void operator=(const TinyVector<T,SIZE>& a_rhs)
		{
			for (int i=0; i<SIZE; i++) m_data[i] = a_rhs.m_data[i];
		}
		
		/* Increment/Decrement/Scale */
		void operator+=(const TinyVector<T,SIZE>& a_rhs)
		{
			for (int i=0; i<SIZE; i++) m_data[i] += a_rhs.m_data[i];
		}

		void operator-=(const TinyVector<T,SIZE>& a_rhs)
		{
			for (int i=0; i<SIZE; i++) m_data[i] -= a_rhs.m_data[i];
		}

		void scale(T a_value)
		{
			for (int i=0; i<SIZE; i++) m_data[i] *= a_value;
		}
		
		/* Set value to zero */
		void setToZero()
		{
			for (int i=0; i<SIZE; i++) m_data[i] = 0;
		}
		
		/* Access to elements (Matlab-like) */
		T  operator()(int a_index) const {return m_data[a_index];}
		T& operator()(int a_index)       {return m_data[a_index];}
		
		/* access to pointer */
		T*       getPointer()       {return m_data;}
		const T* getPointer() const {return m_data;}
		
		/* norm */
		T norm() const {return sqrt(innerProduct(*this,*this));}
				
	private:
		
		T m_data[SIZE];
	};


/* ======================================================================================== */
/* Friend functions                                                                         */
/* ======================================================================================== */

/* Addition */
template <class T, int SIZE> 
inline TinyVector<T,SIZE> operator+(const TinyVector<T,SIZE>& a_vec1, const TinyVector<T,SIZE>& a_vec2)
{
	TinyVector<T,SIZE> ret;
	for (int i=0; i<SIZE; i++) ret.m_data[i] = a_vec1.m_data[i] + a_vec2.m_data[i];
	return ret;
}

/* Subtraction */
template <class T, int SIZE> 
inline TinyVector<T,SIZE> operator-(const TinyVector<T,SIZE>& a_vec1, const TinyVector<T,SIZE>& a_vec2)
{
	TinyVector<T,SIZE> ret;
	for (int i=0; i<SIZE; i++) ret.m_data[i] = a_vec1.m_data[i] - a_vec2.m_data[i];
	return ret;
}

/* Inner product */
template <class T, int SIZE> 
inline T innerProduct(const TinyVector<T,SIZE>& a_vec1, const TinyVector<T,SIZE>& a_vec2)
{
	T ret=0;
	for (int i=0; i<SIZE; i++) ret += a_vec1.m_data[i]*a_vec2.m_data[i];
	return ret;
}

/* Cross product (only if SIZE=3) */
template <class T, int SIZE> 
inline TinyVector<T, SIZE> 
CrossProduct(const TinyVector<T, SIZE>& a_vec1, const TinyVector<T, SIZE>& a_vec2)
{
	assert(SIZE==3);
	TinyVector<double, 3> ret;
	ret.m_data[0] = a_vec1.m_data[1]*a_vec2.m_data[2] - a_vec1.m_data[2]*a_vec2.m_data[1];
	ret.m_data[1] = a_vec1.m_data[2]*a_vec2.m_data[0] - a_vec1.m_data[0]*a_vec2.m_data[2];
	ret.m_data[2] = a_vec1.m_data[0]*a_vec2.m_data[1] - a_vec1.m_data[1]*a_vec2.m_data[0];
	return ret;
}

/* Print the vector to stream */
template <class T, int SIZE> 
std::ostream& operator<<(std::ostream &a_os, const TinyVector<T, SIZE>& a_vec) 
{
	a_os << "[";
	for (int i=0; i<SIZE; i++) a_os << a_vec(i) << " ";
	a_os << "]" << std::endl;
	return a_os;
}

#endif
