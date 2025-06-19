/*
 Matrix.h
 
 A templated reference-counted matrix.
 */

#ifndef _Matrix_h_
#define _Matrix_h_

#include "Main.H"
#include "RefCounter.H"
#include "VectorIterator.H"


template <class T> class Matrix;
template <class T> Matrix<T> operator- (const Matrix<T> &);
template <class T> Matrix<T> operator+ (const Matrix<T> &, const Matrix<T> &);
template <class T> Matrix<T> operator- (const Matrix<T> &, const Matrix<T> &);
template <class T> Matrix<T> operator* (const Matrix<T> &, const Matrix<T> &);
template <class T> Matrix<T> operator/ (const Matrix<T> &, const Matrix<T> &);
template <class T> Matrix<T> operator+ (const Matrix<T> &, T);
template <class T> Matrix<T> operator- (const Matrix<T> &, T);
template <class T> Matrix<T> operator* (const Matrix<T> &, T);
template <class T> Matrix<T> operator/ (const Matrix<T> &, T);
template <class T> Matrix<T> operator+ (T, const Matrix<T> &);
template <class T> Matrix<T> operator- (T, const Matrix<T> &);
template <class T> Matrix<T> operator* (T, const Matrix<T> &);
template <class T> Matrix<T> operator/ (T, const Matrix<T> &);
template <class T> Matrix<T> applyFunction(const Matrix<T> &, T (*)(const T));
template <class T> Matrix<T> sin (const Matrix<T> &);
template <class T> Matrix<T> cos (const Matrix<T> &);
template <class T> Matrix<T> tan (const Matrix<T> &);
template <class T> Matrix<T> exp (const Matrix<T> &);
template <class T> Matrix<T> log (const Matrix<T> &);
template <class T> Matrix<T> sqrt (const Matrix<T> &);
template <class T> Matrix<T> pow2 (const Matrix<T> &);
template <class T> bool EqualSize(const Matrix<T> &, const Matrix<T> &);

template <class T>
class Matrix
	{
	public:
		
		typedef MatrixIterator<T>        iterator;
		typedef ConstMatrixIterator<T>   const_iterator;
				
		friend Matrix<T> operator- <>(const Matrix<T> &);
		friend Matrix<T> operator+ <>(const Matrix<T> &, const Matrix<T> &);
		friend Matrix<T> operator- <>(const Matrix<T> &, const Matrix<T> &);
		friend Matrix<T> operator* <>(const Matrix<T> &, const Matrix<T> &);
		friend Matrix<T> operator/ <>(const Matrix<T> &, const Matrix<T> &);
		friend Matrix<T> operator+ <>(const Matrix<T> &, T);
		friend Matrix<T> operator- <>(const Matrix<T> &, T);
		friend Matrix<T> operator* <>(const Matrix<T> &, T);
		friend Matrix<T> operator/ <>(const Matrix<T> &, T);
		friend Matrix<T> operator+ <>(T, const Matrix<T> &);
		friend Matrix<T> operator- <>(T, const Matrix<T> &);
		friend Matrix<T> operator* <>(T, const Matrix<T> &);
		friend Matrix<T> operator/ <>(T, const Matrix<T> &);
		friend Matrix<T> applyFunction <>(const Matrix<T> &, T (*)(const T));
		friend Matrix<T> sin <>(const Matrix<T> &);
		friend Matrix<T> cos <>(const Matrix<T> &);
		friend Matrix<T> tan <>(const Matrix<T> &);
		friend Matrix<T> exp <>(const Matrix<T> &);
		friend Matrix<T> log <>(const Matrix<T> &);
		friend Matrix<T> sqrt <>(const Matrix<T> &);
		friend Matrix<T> pow2 <>(const Matrix<T> &);
		
		/* constructors (reference counting) */
		Matrix();        
		Matrix(int i, int j);
		Matrix(const Matrix<T> &rhs); 
		
		/* destructor */
		~Matrix();
		
		/* access */
		T     operator()(int i, int j) const   {return m_ptr[i + m_rows*j];}
		T     operator()(int i)        const   {return m_ptr[i];}
		T&    operator()(int i, int j)         {return m_ptr[i + m_rows*j];}
		T&    operator()(int i)                {return m_ptr[i];}
		
		/* return the size */
		int rows() const;
		int cols() const;
		
		/* pointwise assignments */
		void operator= (const Matrix<T> &rhs);
		void operator+=(const Matrix<T> &rhs);
		void operator-=(const Matrix<T> &rhs);
		void operator*=(const Matrix<T> &rhs);
		void operator/=(const Matrix<T> &rhs);
		
		void operator= (T val);   /* data not copied. reference counting! */
		void operator+=(T val);
		void operator-=(T val);
		void operator*=(T val);
		void operator/=(T val);
		
		/* clone the matrix */
		Matrix<T> clone() const;
		
		/* return the pointer to the data */
		T*       getPointer()           {return m_ptr;}
		const T* getPointer() const     {return m_ptr;}
		
		/* Returns matrix iterators (const and non-const) */
		ConstMatrixIterator<T> begin() const 
		{
			return ConstMatrixIterator<T>(m_ptr, m_ptr, m_ptr + m_rows*m_cols);
		}
		
		MatrixIterator<T> begin()  
		{
			return MatrixIterator<T>(m_ptr, m_ptr, m_ptr + m_rows*m_cols);
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
		
		int        m_rows, m_cols;
		T*         m_ptr;
		RefCounter m_cnt;
		
	};

/*
 Matrix.cpp
 
 A templated reference-counted matrix.
 */


/* Member function */

/* Default constructor */
template <class T>
Matrix<T>::Matrix() : 
m_rows(0), 
m_cols(0), 
m_ptr(NULL),
m_cnt()                       
{
}

/* Constrctor with size */
template <class T>
Matrix<T>::Matrix(int i, int j) : 
m_rows(i), 
m_cols(j), 
m_ptr(new T[m_rows*m_cols]),
m_cnt()
{
	(*this) = T();
}

/* Copy constructor */
template <class T>
Matrix<T>::Matrix(const Matrix<T> &rhs) : 
m_rows(rhs.m_rows), 
m_cols(rhs.m_cols), 
m_ptr(rhs.m_ptr),
m_cnt(rhs.m_cnt) 
{
}

/* Destructor */
template <class T>
Matrix<T>::~Matrix()  
{
	if (m_cnt.is_unique() && m_ptr != NULL) delete [] m_ptr;
}

/* Returns the nuber of rows */ 
template <class T>
int Matrix<T>::rows() const 
{
	return m_rows;
}

/* Returns the number of columns */
template <class T>
int Matrix<T>::cols() const 
{
	return m_cols;
}  

/* Make a brand new copy */
template <class T>
Matrix<T> Matrix<T>::clone() const
{
	Matrix<T> ret(m_rows,m_cols);
	T* ret_m_ptr = ret.getPointer();
	for (int i=0; i<m_rows*m_cols; i++)
		ret_m_ptr[i] = m_ptr[i];    
	return ret;
}

/* Assignment operator */
template <class T>
void Matrix<T>::operator=(const Matrix<T> &rhs)
{
	if (m_ptr != rhs.m_ptr)
	{
		if (m_cnt.is_unique() && m_ptr != NULL) delete [] m_ptr;
		m_ptr = rhs.m_ptr;
		m_cnt = rhs.m_cnt;
		m_rows = rhs.m_rows;
		m_cols = rhs.m_cols;
	}
}

/* Increment by a matrix */
template <class T>
void Matrix<T>::operator+=(const Matrix<T> &rhs)
{
	assert(EqualSize(*this,rhs));
	for (int i=0; i<m_rows*m_cols; i++)
		m_ptr[i] += rhs.m_ptr[i];    
}

/* Decrement by a matrix */
template <class T>
void Matrix<T>::operator-=(const Matrix<T> &rhs)
{
	assert(EqualSize(*this,rhs));
	for (int i=0; i<m_rows*m_cols; i++)
		m_ptr[i] -= rhs.m_ptr[i];    
}

/* Multiplication by a matrix */
template <class T>
void Matrix<T>::operator*=(const Matrix<T> &rhs)
{
	assert(EqualSize(*this,rhs));
	for (int i=0; i<m_rows*m_cols; i++)
		m_ptr[i] *= rhs.m_ptr[i];    
}

/* Division by matrix */
template <class T>
void Matrix<T>::operator/=(const Matrix<T> &rhs)
{
	assert(EqualSize(*this,rhs));
	for (int i=0; i<m_rows*m_cols; i++)
		m_ptr[i] /= rhs.m_ptr[i];    
}

/* Assignment of a scalar */
template <class T>
void Matrix<T>::operator=(T val)
{
	for (int i=0; i<m_rows*m_cols; i++)
		m_ptr[i] = val;    
}

/* Increment by a scalar */
template <class T>
void Matrix<T>::operator+=(T val)
{
	for (int i=0; i<m_rows*m_cols; i++)
		m_ptr[i] += val;    
}

/* Decrement by a scalar */
template <class T>
void Matrix<T>::operator-=(T val)
{
	for (int i=0; i<m_rows*m_cols; i++)
		m_ptr[i] -= val;    
}

/* Multiplication by a scalar */
template <class T>
void Matrix<T>::operator*=(T val)
{
	for (int i=0; i<m_rows*m_cols; i++)
		m_ptr[i] *= val;    
}

/* Division by a scalar */
template <class T>
void Matrix<T>::operator/=(T val)
{
	for (int i=0; i<m_rows*m_cols; i++)
		m_ptr[i] /= val;    
}

/* Returns the smallest element */
template <class T>
T Matrix<T>::minimum() const
{
	T ret = m_ptr[0];
	for (int i=1; i<m_rows*m_cols; i++)
		ret = min(ret,m_ptr[i]);
	return ret;
}

/* Returns the largest element */
template <class T>
T Matrix<T>::maximum() const
{
	T ret = m_ptr[0];
	for (int i=1; i<m_rows*m_cols; i++)
		ret = max(ret,m_ptr[i]);
	return ret;
}

/* Returns the sum of the elements */
template <class T>
T Matrix<T>::sum() const
{
	T ret = 0;
	for (int i=0; i<m_rows*m_cols; i++)
		ret += m_ptr[i];
	return ret;
}

/* Returns the mean value */
template <class T>
T Matrix<T>::mean() const
{
	return sum()/double(m_rows)/double(m_cols);
}

/* Returns the L1 nurm */
template <class T>
T Matrix<T>::norm1() const
{
	T ret = 0;
	for (int i=0; i<m_rows*m_cols; i++)
		ret += fabs(m_ptr[i]);
	return ret;
}

/* Returns the L2 norm */
template <class T>
T Matrix<T>::norm2() const
{
	T ret = 0;
	for (int i=0; i<m_rows*m_cols; i++)
		ret += m_ptr[i]*m_ptr[i];
	return sqrt(ret);
}

/* Returns the Linfty norm */
template <class T>
T Matrix<T>::norminf() const
{
	T ret = 0;
	for (int i=0; i<m_rows*m_cols; i++)
		ret = max(ret, fabs(m_ptr[i]));;
	return ret;
}

/* Friend functions */


/* Additive inverse */
template <class T>
Matrix<T> operator-(const Matrix<T> &a_matrix)
{
	Matrix<T> ret(a_matrix.rows(),a_matrix.cols());
	ret -= a_matrix;
	return ret;
}

/* Binary addition */
template <class T>
Matrix<T> operator+(const Matrix<T> &a_matrix1, const Matrix<T> &a_matrix2)
{
	assert(EqualSize(a_matrix1,a_matrix2));
	Matrix<T> ret = a_matrix1.clone();
	ret += a_matrix2;
	return ret;
}

/* Binary subtraction */
template <class T>
Matrix<T> operator-(const Matrix<T> &a_matrix1, const Matrix<T> &a_matrix2)
{
	assert(EqualSize(a_matrix1,a_matrix2));
	Matrix<T> ret = a_matrix1.clone();
	ret -= a_matrix2;
	return ret;
}


/*Binary multiplication */
template <class T>
Matrix<T> operator*(const Matrix<T> &a_matrix1, const Matrix<T> &a_matrix2)
{
	assert(EqualSize(a_matrix1,a_matrix2));
	Matrix<T> ret = a_matrix1.clone();
	ret *= a_matrix2;
	return ret;
}

/* Binary division */
template <class T>
Matrix<T> operator/(const Matrix<T> &a_matrix1, const Matrix<T> &a_matrix2)
{
	assert(EqualSize(a_matrix1,a_matrix2));
	Matrix<T> ret = a_matrix1.clone();
	ret /= a_matrix2;
	return ret;
}

/* Addition of vector and scalar */
template <class T>
Matrix<T> operator+(const Matrix<T> &a_matrix, T a_value)
{
	Matrix<T> ret = a_matrix.clone();
	ret += a_value;
	return ret;
}

/* Subtraction of vector and scalar */
template <class T>
Matrix<T> operator-(const Matrix<T> &a_matrix, T a_value)
{
	Matrix<T> ret = a_matrix.clone();
	ret -= a_value;
	return ret;
}

/* Multiplication of vector and scalar */
template <class T>
Matrix<T> operator*(const Matrix<T> &a_matrix, T a_value)
{
	Matrix<T> ret = a_matrix.clone();
	ret *= a_value;
	return ret;
}

/* Division of vector and scalar */
template <class T>
Matrix<T> operator/(const Matrix<T> &a_matrix, T a_value)
{
	Matrix<T> ret = a_matrix.clone();
	ret /= a_value;
	return ret;
}

/* Addition of scalar and vector */
template <class T>
Matrix<T> operator+(T a_value, const Matrix<T> &a_matrix)
{
	Matrix<T> ret = a_matrix.clone();
	ret += a_value;
	return ret;
}

/* Subtraction of scalar and vector */
template <class T>
Matrix<T> operator-(T a_value, const Matrix<T> &a_matrix)
{
	Matrix<T> ret(a_matrix.rows(),a_matrix.cols());
	ret  = a_value;
	ret -= a_matrix;
	return ret;
}

/* Multiplication of scalar and vector */
template <class T>
Matrix<T> operator*(T a_value, const Matrix<T> &a_matrix)
{
	Matrix<T> ret = a_matrix.clone();
	ret *= a_value;
	return ret;
}

/* Division of scalar and  vector */
template <class T>
Matrix<T> operator/(T a_value, const Matrix<T> &a_matrix)
{
	Matrix<T> ret(a_matrix.rows(),a_matrix.cols());
	ret  = a_value;
	ret /= a_matrix;
	return ret;
}

/* Unary function application */
template <class T>
Matrix<T> applyFunction(const Matrix<T> &a_matrix, T (*funPtr)(const T))
{
	Matrix<T> ret(a_matrix.rows(),a_matrix.cols());
	
	MatrixIterator<T>      vit  = ret.begin();
	ConstMatrixIterator<T> cvit = a_matrix.begin();
	
	for (; !vit.end(); vit++, cvit++)
		*vit = (*funPtr)(*cvit);
	
	return ret;
}

/* sin() funciton */
template <class T>
Matrix<T> sin(const Matrix<T> &a_matrix)
{
	return applyFunction(a_matrix, &sin);
}

/* cos() function */
template <class T>
Matrix<T> cos(const Matrix<T> &a_matrix)
{
	return applyFunction(a_matrix, &cos);
}

/* tan() function */
template <class T>
Matrix<T> tan(const Matrix<T> &a_matrix)
{
	return applyFunction(a_matrix, &tan);
}

/* exp() function */
template <class T>
Matrix<T> exp(const Matrix<T> &a_matrix)
{
	return applyFunction(a_matrix, &exp);
}

/* log() function */
template <class T>
Matrix<T> log(const Matrix<T> &a_matrix)
{
	return applyFunction(a_matrix, &log);
}

/* sqrt() function */
template <class T>
Matrix<T> sqrt(const Matrix<T> &a_matrix)
{
	return applyFunction(a_matrix, &sqrt);
}

/* Square function */
template <class T>
Matrix<T> pow2(const Matrix<T> &a_matrix)
{
	Matrix<T> ret(a_matrix.rows(),a_matrix.cols());
	for (int i=0; i<a_matrix.m_rows*a_matrix.m_cols; i++)
		ret.m_ptr[i] = a_matrix.m_ptr[i]*a_matrix.m_ptr[i];
	return ret;
}

/* Check that two matrices are of the same size */
template <class T> 
bool EqualSize(const Matrix<T> &a_vec1, const Matrix<T> &a_vec2)
{
	return (a_vec1.rows() == a_vec2.rows()) && (a_vec1.cols() == a_vec2.cols());
}

#endif
