// int array2.cc -- Code for IntArray2 class, an array of ints with two indices

#include "intarray2.h"

void error(const std::string& s);

// Constructor/Resize

void IntArray2::make(int o1,int s1,int o2,int s2)
    {
    if (s1 < 0 || s2 < 0)
	error("IntArray2::make: Bad size");
    size1 = s1;
    size2 = s2;
    off1 = o1;
    off2 = o2;
    array1->ReSize(s1*s2);
    }

// Copy constructor

IntArray2::IntArray2(const IntArray2& V)
    {
    size1 = V.size1;
    size2 = V.size2;
    off1 = V.off1;
    off2 = V.off2;
    array1 = V.array1;
    }

// Assignment

void IntArray2::operator=(const IntArray2& V)
    {
    size1 = V.size1;
    size2 = V.size2;
    off1 = V.off1;
    off2 = V.off2;
    array1 = V.array1;
    }

// Assignment, grabbing storage

void IntArray2::CopyDestroy(IntArray2& V)
    {
    size1 = V.size1;
    size2 = V.size2;
    off1 = V.off1;
    off2 = V.off2;
    array1 = V.array1;
    }

// Assignment to an int value

void IntArray2::operator=(int value)
    {
    array1 = value;
    }

// Output an array

std::ostream& operator<<(std::ostream& s, const IntArray2& V)
    {
    //long f = s.flags();

    s << "\n";
    int k = 0;
    for (int i = 0; i < V.size1; i++)
	{
	for (int j = 0; j < V.size2; j++)
	    s << V.array1(k++) << " ";
	s << "\n";
	}

    s << std::endl;
    //s.flags(f);
    return s;
    }
