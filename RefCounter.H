#ifndef REFCOUNTER_H
#define REFCOUNTER_H

/*
  This concrete class implements a "reference counter".
  A RefCounter is a counter that counts how many
  copies of itself exist, and destroys itself when this
  number is zero. Reference counters are commonly used as
  members in larger objects in order to reduce the number
  of time large structures need to be copied.
  
  Examples:
  
  RefCounter cnt1;        // *** New counter set to 1.
  RefCounter cnt2 = cnt1; // *** both are equal to 2.
  RefCounter cnt3;        // *** New counter set to 1;
  cnt3 = cnt1;            // *** All are equal to 3. 
*/

#include <Main.H>

class RefCounter
{
private:

    int* m_ptr;     /* The pointer to the counter */
  
public:

    /* Default constructor */
    RefCounter() : m_ptr(new int(1)) {}
	
	/* Copy constructor */
	RefCounter(const RefCounter &rhs) : m_ptr(rhs.m_ptr)   {increment();}  
	
    /* Desctructor: decrement the counter */
    ~RefCounter()   {decrement();}
  
    /* Assignment: increment RHS and decrement LHS */
	void operator=(const RefCounter &rhs)
	{ 
		decrement();
		m_ptr = rhs.m_ptr;
		increment();
	}
	
    /* is_unique() is TRUE if the counter equals 1 */
	bool is_unique() const    {return (*m_ptr)==1;}
	
private: 

    /* Helper functions */  
    void increment()   {(*m_ptr)++;}
	
	void decrement()
	{
		if ((*m_ptr)==1)
		{  
			delete m_ptr;  
			m_ptr = NULL; 
		}
		else
			--(*m_ptr);
	}
            
};


#endif

