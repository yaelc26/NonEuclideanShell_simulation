/*
 *  VectorIterator.H
 *  RKLibrary
 *
 *  Created by Raz Kupferman on 9/6/09.
 *  Copyright 2009 The Hebrew University. All rights reserved.
 *
 */


#ifndef _VECTORITERATOR_H
#define _VECTORITERATOR_H

#include "Main.H"


template <class T> class Vector;
template <class T> class Matrix;
template <class T> class ConstVectorIterator;
template <class T> class VectorIterator;

/* We define both const and non-const vector iterators. This iterator can be used both 
 for Vector<T> and Matrix<T>, but we define MatrixIterator and ConstMatrixIterator<T> to
 have a clean code */


template <class T>
class ConstVectorIterator
	{
	public:
		
		/* Forbid default constructor */
		ConstVectorIterator();
		
		/* Constructor with pointers */
		ConstVectorIterator(T* a_start, T* a_current, T* a_end)
		{
			m_start   = a_start;
			m_current = a_current;
			m_end     = a_end;
		}
		
		/* Copy constructor */
		ConstVectorIterator(const ConstVectorIterator<T> &a_vit)
		{
			m_start   = a_vit.m_start;
			m_current = a_vit.m_current;
			m_end     = a_vit.m_end;
		}
		
		/* Inrementation */
		void operator++(int)  {m_current++;}
		
		/* Check whether we got to the end of the vector */
		bool end() const    {return (m_current==m_end);}
		
		/* Dereferencing */
		T        operator*() const {return *m_current;}
		
	protected:
		
		T*  m_start;
		T*  m_current;
		T*  m_end;
	};


/* ********************************************************************************** */
template <class T>
class VectorIterator : public ConstVectorIterator<T>
{
public:
	
	/* Forbid default constructor */
	VectorIterator();
	
	/* Constructor with pointers */
	VectorIterator(T* a_start, T* a_current, T* a_end) : 
	ConstVectorIterator<T>(a_start, a_current, a_end) 
	{}
	
	/* Copy constructor */
	VectorIterator(const ConstVectorIterator<T> &a_vit) :
	ConstVectorIterator<T>(a_vit.m_start, a_vit.a_current, a_vit.a_end) 
	{}
	
	/* Non-constant dereferencing */
	// T& operator*() const {return *ConstVectorIterator<T>::m_current;}
	T& operator*() const {return *(this->m_current);}
	
};

/* ********************************************************************************** */
template <class T>
class ConstMatrixIterator : public ConstVectorIterator<T>
{
public:
	
	/* Forbid default constructor */
	ConstMatrixIterator();
	
	/* Constructor with pointers */
	ConstMatrixIterator(T* a_start, T* a_current, T* a_end) : 
	ConstVectorIterator<T>(a_start, a_current, a_end) 
	{}
	
	/* Copy constructor */
	ConstMatrixIterator(const ConstMatrixIterator<T> &a_vit) :
	ConstVectorIterator<T>(a_vit.m_start, a_vit.m_current, a_vit.m_end) 
	{}
	
};
/* ********************************************************************************** */
template <class T>
class MatrixIterator : public VectorIterator<T>
{
public:
	
	/* Forbid default constructor */
	MatrixIterator();
	
	/* Constructor with pointers */
	MatrixIterator(T* a_start, T* a_current, T* a_end) : 
	VectorIterator<T>(a_start, a_current, a_end) 
	{}
	
	/* Copy constructor */
	MatrixIterator(const MatrixIterator<T> &a_vit) :
	VectorIterator<T>(a_vit.m_start, a_vit.m_current, a_vit.m_end) 
	{}
	
};



#endif
