// array2.h -- integer array with two indices, starting at offset.

#ifndef ARRAY2
#define ELEMENT int		/* type of element */
#define ARRAY2 IntArray2	/* name of class */
#endif

#include <iostream.h>
//#define BOUNDS		/* Define this if you want bounds checking */

class ARRAY2 {
    int size1,size2;
    ELEMENT * store;
    int off1,off2;
    void make(int,int,int,int);		// Real Resize/Constructor
    
  public:
    ARRAY2()			// Default Constructor
	{ store = 0; size1 = 0; size2 = 0; off1 = 0; off2 = 0;}	
    ARRAY2(const ARRAY2&);	// Copy constructor
    ARRAY2(int s1,int s2)	// Construct with a certain size, starting at 1
	{ store = 0 ; make(1,s1,1,s2); }
    ARRAY2(int of1,int lim1,int of2,int lim2)	// Construct with offset
	{ store = 0; make(of1,lim1-of1+1,of2,lim2-of2+1); }
    ~ARRAY2() { make(0,0,0,0); }	// Destructor
    int Size1() const { return size1; } // Information functions
    int Lower1() const { return off1; }
    int Upper1() const { return off1+size1-1; }
    int Size2() const { return size2; } // Information functions
    int Lower2() const { return off2; }
    int Upper2() const { return off2+size2-1; }
    ELEMENT* Store() const { return store; }
    void ReDimension(int s1,int s2)	// Redimension
	{ make(1,s1,1,s2); }
    // Redimension with offset
    void ReDimension(int of1,int lim1,int of2,int lim2)
	{ make(of1,lim1-of1+1,of2,lim2-of2+1); }
    void operator=(const ARRAY2&); // Assignment
    void CopyDestroy(ARRAY2&);	   // Copy, grabbing storage
    void operator=(ELEMENT);	   // Assignment to an integer value
    
    ELEMENT& operator()(int,int);	// Access an element
    ELEMENT operator()(int,int) const;
    friend ostream& operator<<(ostream&, const ARRAY2&);
};

// Access an element

inline ELEMENT& ARRAY2::operator()(int i,int j)
{
#ifdef BOUNDS
    void error(const std::string&);
    if (i < off1 || i >= off1+size1 || j < off2 || j >= off2+size2)
      error("ARRAY2: Index out of bounds");
#endif
    return store[(i-off1)*size2 + j - off2];
}

inline ELEMENT ARRAY2::operator()(int i,int j) const
{
#ifdef BOUNDS
    void error(const std::string&);
    if (i < off1 || i >= off1+size1 || j < off2 || j >= off2+size2)
      error("ARRAY2: Index out of bounds");
#endif
    return store[(i-off1)*size2 + j - off2];
}
