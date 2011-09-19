// intarray2.h -- integer array with two indices, starting at offset.

#ifndef _intarray2_h
#define _intarray2_h

#include <iostream>
//#define BOUNDS		/* Define this if you want bounds checking */
#include <tarray1.h>

class IntArray2 
    {
    IntArray1 array1;
    int size1,size2;
    int off1,off2;
    void make(int,int,int,int);		// Real Resize/Constructor
    
  public:
    IntArray2()			// Default Constructor
	{ size1 = 0; size2 = 0; off1 = 0; off2 = 0;}	
    IntArray2(const IntArray2&);	// Copy constructor
    IntArray2(int s1,int s2)	// Construct with a certain size, starting at 1
	{ make(1,s1,1,s2); }
    IntArray2(int of1,int lim1,int of2,int lim2)	// Construct with offset
	{ make(of1,lim1-of1+1,of2,lim2-of2+1); }
    ~IntArray2() { make(0,0,0,0); }	// Destructor
    int Size1() const { return size1; } // Information functions
    int Lower1() const { return off1; }
    int Upper1() const { return off1+size1-1; }
    int Size2() const { return size2; } // Information functions
    int Lower2() const { return off2; }
    int Upper2() const { return off2+size2-1; }
    int* Store() const { return (int *) array1.Store(); }
    void ReDimension(int s1,int s2)	// Redimension
	{ make(1,s1,1,s2); }
    // Redimension with offset
    void ReDimension(int of1,int lim1,int of2,int lim2)
	{ make(of1,lim1-of1+1,of2,lim2-of2+1); }
    void operator=(const IntArray2&); // Assignment
    void CopyDestroy(IntArray2&);	   // Copy, grabbing storage
    void operator=(int);	   // Assignment to an integer value
    
    int& operator()(int,int);	// Access an element
    int operator()(int,int) const;
    friend std::ostream& operator<<(std::ostream&, const IntArray2&);
};

// Access an element

inline int& IntArray2::operator()(int i,int j)
    {
#ifdef BOUNDS
    void error(const char*);
    if (i < off1 || i >= off1+size1 || j < off2 || j >= off2+size2)
      error("IntArray2: Index out of bounds");
#endif
    return array1[(i-off1)*size2 + j - off2];
    }

inline int IntArray2::operator()(int i,int j) const
    {
#ifdef BOUNDS
    void error(const char*);
    if (i < off1 || i >= off1+size1 || j < off2 || j >= off2+size2)
      error("IntArray2: Index out of bounds");
#endif
    return array1((i-off1)*size2 + j - off2);
    }

#endif
