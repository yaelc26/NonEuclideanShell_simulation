/*
 Vector.h
 
 A templated reference-counted Vector.
 Allows arithmetics and simple function evaluations.
 */

#ifndef _Vector_h_
#define _Vector_h_

#include "Main.H"
#include "RefCounter.H"
#include "VectorIterator.H"


/* Non-member arithmetics and function evaluations */
template <class T> class Vector;
template <class T> Vector<T> operator- (const Vector<T> &);
template <class T> Vector<T> operator+ (const Vector<T> &, const Vector<T> &);
template <class T> Vector<T> operator- (const Vector<T> &, const Vector<T> &);
template <class T> Vector<T> operator* (const Vector<T> &, const Vector<T> &);
template <class T> Vector<T> operator/ (const Vector<T> &, const Vector<T> &);
template <class T> Vector<T> operator+ (const Vector<T> &, T);
template <class T> Vector<T> operator- (const Vector<T> &, T);
template <class T> Vector<T> operator* (const Vector<T> &, T);
template <class T> Vector<T> operator/ (const Vector<T> &, T);
template <class T> Vector<T> operator+ (T, const Vector<T> &);
template <class T> Vector<T> operator- (T, const Vector<T> &);
template <class T> Vector<T> operator* (T, const Vector<T> &);
template <class T> Vector<T> operator/ (T, const Vector<T> &);
template <class T> Vector<T> applyFunction (const Vector<T> &, T (*)(const T));
template <class T> Vector<T> sin (const Vector<T> &);
template <class T> Vector<T> cos (const Vector<T> &);
template <class T> Vector<T> tan (const Vector<T> &);
template <class T> Vector<T> exp (const Vector<T> &);
template <class T> Vector<T> log (const Vector<T> &);
template <class T> Vector<T> sqrt(const Vector<T> &);
template <class T> Vector<T> pow2(const Vector<T> &);
template <class T> bool EqualSize(const Vector<T> &, const Vector<T> &);


template <class T>
class Vector
	{
	public:
		
		/* Non-member arithmetics and function evaluations */
		friend Vector<T> operator- <>(const Vector<T> &);
		friend Vector<T> operator+ <>(const Vector<T> &, const Vector<T> &);
		friend Vector<T> operator- <>(const Vector<T> &, const Vector<T> &);
		friend Vector<T> operator* <>(const Vector<T> &, const Vector<T> &);
		friend Vector<T> operator/ <>(const Vector<T> &, const Vector<T> &);
		friend Vector<T> operator+ <>(const Vector<T> &, T);
		friend Vector<T> operator- <>(const Vector<T> &, T);
		friend Vector<T> operator* <>(const Vector<T> &, T);
		friend Vector<T> operator/ <>(const Vector<T> &, T);
		friend Vector<T> operator+ <>(T, const Vector<T> &);
		friend Vector<T> operator- <>(T, const Vector<T> &);
		friend Vector<T> operator* <>(T, const Vector<T> &);
		friend Vector<T> operator/ <>(T, const Vector<T> &);
		
		friend Vector<T> applyFunction <>(const Vector<T> &, T (*)(const T));
		friend Vector<T> sin <>(const Vector<T> &);
		friend Vector<T> cos <>(const Vector<T> &);
		friend Vector<T> tan <>(const Vector<T> &);
		friend Vector<T> exp <>(const Vector<T> &);
		friend Vector<T> log <>(const Vector<T> &);
		friend Vector<T> sqrt <>(const Vector<T> &);
		friend Vector<T> pow2 <>(const Vector<T> &);
		
		/* constructors (reference counted) */
		Vector();        
		Vector(int i);
		Vector(const Vector<T> &); 
		
		/* destructor */
		~Vector();
		
		/* access */
		T   operator()(int i)        const   {return m_ptr[i];}
		T&  operator()(int i)                {return m_ptr[i];}
		
		/* return the size */
		int length() const;

		
		/* pointwise assignments */
		void operator= (const Vector<T> &);
		void operator+=(const Vector<T> &);
		void operator-=(const Vector<T> &);
		void operator*=(const Vector<T> &);
		void operator/=(const Vector<T> &);
		
		void operator= (T);   /* data not copied. reference counting! */
		void operator+=(T);
		void operator-=(T);
		void operator*=(T);
		void operator/=(T);
		
		/* clone the Vector */
		Vector<T> clone() const;
		
		/* return the pointer to the data */
		T*       getPointer()           {return m_ptr;}
		const T* getPointer() const     {return m_ptr;}
		
		/* Returns vector iterators (const and non-const) */
		ConstVectorIterator<T> begin() const 
		{
			return ConstVectorIterator<T>(m_ptr, m_ptr, m_ptr + m_length);
		}
		
		VectorIterator<T> begin()  
		{
			return VectorIterator<T>(m_ptr, m_ptr, m_ptr + m_length);
		}
		
		/* various functions */
		T minimum()  const;
		T maximum()  const;
		T sum()      const;
		T mean()     const;
		T norm1()    const;
		T norm2()    const;
		T norminf()  const;
		
	private:
		
		int        m_length;
		T*         m_ptr;
		RefCounter m_cnt;
		
	};



/* Member functions */

/* Empty constructor */
template <class T>
Vector<T>::Vector() : 
m_length(0), 
m_ptr(NULL),
m_cnt()                       
{
}

/* Constructor with size */
template <class T>
Vector<T>::Vector(int a_length) : 
m_length(a_length), 
m_ptr(new T[m_length]),
m_cnt()
{
	(*this) = T();
}

/* Copy constructor */
template <class T>
Vector<T>::Vector(const Vector<T> &rhs) : 
m_length(rhs.m_length), 
m_ptr(rhs.m_ptr),
m_cnt(rhs.m_cnt) 
{
}

/* Destructor */
template <class T>
Vector<T>::~Vector()  
{
	if (m_cnt.is_unique() && m_ptr != NULL) delete [] m_ptr;
}

/* Return the length */
template <class T>
int Vector<T>::length() const 
{
	return m_length;
}

/* Make a brand new copy */
template <class T>
Vector<T> Vector<T>::clone() const
{
	Vector<T> ret(m_length);
	ConstVectorIterator<T> cvit = begin();
	VectorIterator<T>      vit = ret.begin();
	
	for(; !vit.end(); vit++, cvit++)
		*vit = *cvit;
	
	return ret;
}

/* Assignment operator; data not copied; only pointer */
template <class T>
void Vector<T>::operator=(const Vector<T> &rhs)
{
	if (m_ptr != rhs.m_ptr)
	{
		if (m_cnt.is_unique() && m_ptr != NULL) delete [] m_ptr;
		m_ptr    = rhs.m_ptr;
		m_cnt    = rhs.m_cnt;
		m_length = rhs.m_length;
	}
}

/* Pointwise addition */
template <class T>
void Vector<T>::operator+=(const Vector<T> &rhs)
{
	assert(EqualSize(*this,rhs));
	ConstVectorIterator<T> cvit = rhs.begin();
	VectorIterator<T>      vit  = begin();
	
	for(; !vit.end(); vit++, cvit++)
		*vit += *cvit;
}

/* Pointwise subtraction */
template <class T>
void Vector<T>::operator-=(const Vector<T> &rhs)
{
	assert(EqualSize(*this,rhs));
	ConstVectorIterator<T> cvit = rhs.begin();
	VectorIterator<T>      vit  = begin();
	
	for(; !vit.end(); vit++, cvit++)
		*vit -= *cvit;
}


/* Pointwise multiplication */
template <class T>
void Vector<T>::operator*=(const Vector<T> &rhs)
{
	assert(EqualSize(*this,rhs));
	ConstVectorIterator<T> cvit = rhs.begin();
	VectorIterator<T>      vit  = begin();
	
	for(; !vit.end(); vit++, cvit++)
		*vit *= *cvit;
}

/* Pointwise divistion */
template <class T>
void Vector<T>::operator/=(const Vector<T> &rhs)
{
	assert(EqualSize(*this,rhs));
	ConstVectorIterator<T> cvit = rhs.begin();
	VectorIterator<T>      vit  = begin();
	
	for(; !vit.end(); vit++, cvit++)
		*vit /= *cvit;
}

/* Constant assignment */
template <class T>
void Vector<T>::operator=(T a_value)
{
	for(VectorIterator<T> vit  = begin(); !vit.end(); vit++)
		*vit = a_value;
}

/* Constant increment */
template <class T>
void Vector<T>::operator+=(T a_value)
{
	for(VectorIterator<T> vit  = begin(); !vit.end(); vit++)
		*vit += a_value;
}

/* Ccnstant decrement */ 
template <class T>
void Vector<T>::operator-=(T a_value)
{
	for(VectorIterator<T> vit  = begin(); !vit.end(); vit++)
		*vit -= a_value;
}

/* Constant multiplication */
template <class T>
void Vector<T>::operator*=(T a_value)
{
	for(VectorIterator<T> vit  = begin(); !vit.end(); vit++)
		*vit *= a_value;
}

/* Constant division */
template <class T>
void Vector<T>::operator/=(T a_value)
{
	for(VectorIterator<T> vit  = begin(); !vit.end(); vit++)
		*vit /= a_value;
}

/* Smallest element */
template <class T>
T Vector<T>::minimum() const
{
	T ret = m_ptr[0];
	for(ConstVectorIterator<T> cvit = begin(); !cvit.end(); cvit++)
		ret = min(ret,*cvit);
	return ret;
}

/* Largest element */
template <class T>
T Vector<T>::maximum() const
{
	T ret = m_ptr[0];
	for(ConstVectorIterator<T> cvit = begin(); !cvit.end(); cvit++)
		ret = max(ret,*cvit);
	return ret;
}

/* Sum of elements */
template <class T>
T Vector<T>::sum() const
{
	T ret = 0;
	for(ConstVectorIterator<T> cvit = begin(); !cvit.end(); cvit++)
		ret += *cvit;
	return ret;
}

/* Mean value */
template <class T>
T Vector<T>::mean() const
{
	return sum()/double(m_length);
}

/* L1 norm */
template <class T>
T Vector<T>::norm1() const
{
	T ret = 0;
	for(ConstVectorIterator<T> cvit = begin(); !cvit.end(); cvit++)
		ret += fabs(*cvit);
	return ret;
}

/* L2 norm */
template <class T>
T Vector<T>::norm2() const
{
	T ret = 0;
	for(ConstVectorIterator<T> cvit = begin(); !cvit.end(); cvit++)
		ret += (*cvit) * (*cvit);
	return sqrt(ret);
}

/* Linfinity norm */
template <class T>
T Vector<T>::norminf() const
{
	T ret = 0;
	for(ConstVectorIterator<T> cvit = begin(); !cvit.end(); cvit++)
		ret = max(ret, fabs(*cvit));;
	return ret;
}


/* FRIEND FUNCTIONS */


/* Additive inverse */
template <class T>
Vector<T> operator-(const Vector<T> &a_vector)
{
	Vector<T> ret(a_vector.length());
	ret -= a_vector;
	return ret;
}

/* Binary addition */
template <class T>
Vector<T> operator+(const Vector<T> &a_vector, const Vector<T> &a_vector2)
{
	assert(EqualSize(a_vector,a_vector2));
	Vector<T> ret = a_vector.clone();
	ret += a_vector2;
	return ret;
}

/* Binary subtraction */
template <class T>
Vector<T> operator-(const Vector<T> &a_vector, const Vector<T> &a_vector2)
{
	assert(EqualSize(a_vector,a_vector2));
	Vector<T> ret = a_vector.clone();
	ret -= a_vector2;
	return ret;
}


/*Binary multiplication */
template <class T>
Vector<T> operator*(const Vector<T> &a_vector, const Vector<T> &a_vector2)
{
	assert(EqualSize(a_vector,a_vector2));
	Vector<T> ret = a_vector.clone();
	ret *= a_vector2;
	return ret;
}

/* Binary division */
template <class T>
Vector<T> operator/(const Vector<T> &a_vector, const Vector<T> &a_vector2)
{
	assert(EqualSize(a_vector,a_vector2));
	Vector<T> ret = a_vector.clone();
	ret /= a_vector2;
	return ret;
}

/* Addition of vector and scalar */
template <class T>
Vector<T> operator+(const Vector<T> &a_vector, T a_value)
{
	Vector<T> ret = a_vector.clone();
	ret += a_value;
	return ret;
}

/* Subtraction of vector and scalar */
template <class T>
Vector<T> operator-(const Vector<T> &a_vector, T a_value)
{
	Vector<T> ret = a_vector.clone();
	ret -= a_value;
	return ret;
}

/* Multiplication of vector and scalar */
template <class T>
Vector<T> operator*(const Vector<T> &a_vector, T a_value)
{
	Vector<T> ret = a_vector.clone();
	ret *= a_value;
	return ret;
}

/* Division of vector and scalar */
template <class T>
Vector<T> operator/(const Vector<T> &a_vector, T a_value)
{
	Vector<T> ret = a_vector.clone();
	ret /= a_value;
	return ret;
}

/* Addition of scalar and vector */
template <class T>
Vector<T> operator+(T a_value, const Vector<T> &a_vector)
{
	Vector<T> ret = a_vector.clone();
	ret += a_value;
	return ret;
}

/* Subtraction of scalar and vector */
template <class T>
Vector<T> operator-(T a_value, const Vector<T> &a_vector)
{
	Vector<T> ret(a_vector.length());
	ret  = a_value;
	ret -= a_vector;
	return ret;
}

/* Multiplication of scalar and vector */
template <class T>
Vector<T> operator*(T a_value, const Vector<T> &a_vector)
{
	Vector<T> ret = a_vector.clone();
	ret *= a_value;
	return ret;
}

/* Division of scalar and  vector */
template <class T>
Vector<T> operator/(T a_value, const Vector<T> &a_vector)
{
	Vector<T> ret(a_vector.length());
	ret  = a_value;
	ret /= a_vector;
	return ret;
}

/* Unary function application */
template <class T>
Vector<T> applyFunction(const Vector<T> &a_vector, T (*funPtr)(const T))
{
	Vector<T> ret(a_vector.length());
	
	VectorIterator<T>      vit  = ret.begin();
	ConstVectorIterator<T> cvit = a_vector.begin();
	
	for (; !vit.end(); vit++, cvit++)
		*vit = (*funPtr)(*cvit);
	
	return ret;
}

/* sin() funciton */
template <class T>
Vector<T> sin(const Vector<T> &a_vector)
{
	return applyFunction(a_vector, &sin);
}

/* cos() function */
template <class T>
Vector<T> cos(const Vector<T> &a_vector)
{
	return applyFunction(a_vector, &cos);
}

/* tan() function */
template <class T>
Vector<T> tan(const Vector<T> &a_vector)
{
	return applyFunction(a_vector, &tan);
}

/* exp() function */
template <class T>
Vector<T> exp(const Vector<T> &a_vector)
{
	return applyFunction(a_vector, &exp);
}

/* log() function */
template <class T>
Vector<T> log(const Vector<T> &a_vector)
{
	return applyFunction(a_vector, &log);
}

/* sqrt() function */
template <class T>
Vector<T> sqrt(const Vector<T> &a_vector)
{
	return applyFunction(a_vector, &sqrt);
}

/* Square function */
template <class T>
Vector<T> pow2(const Vector<T> &a_vector)
{
	Vector<T> ret(a_vector.length());

	VectorIterator<T>      vit  = ret.begin();
	ConstVectorIterator<T> cvit = a_vector.begin();
	
	for (; !vit.end(); vit++, cvit++)
		*vit = (*cvit)*(*cvit);
	return ret;
}

/* Check that two vectors are of the same size */
template <class T> 
bool EqualSize(const Vector<T> &a_vec1, const Vector<T> &a_vec2)
{
	return (a_vec1.length() == a_vec2.length());
}


#endif
