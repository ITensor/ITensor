// array2.cc -- Code for ARRAY2 class, an array of ELEMENTs with two indices

//#include "array2.h" /* Use cat in make to make this file */

void error(const string& s);

// Constructor/Resize

void ARRAY2::make(int o1,int s1,int o2,int s2)
{
    if (store != 0 && size1*size2 == 0)
      error("ARRAY2::make Internal inconsistency");
    if (store != 0 && o1 == off1 && o2 == off2 && s1 == size1 && s2 == size2)
      return;
    if (s1 < 0 || s2 < 0) error("ARRAY2::make: Bad size");
    
    if (store != 0) {
	delete[] store;
	store = 0;
    }

    if (s1*s2 > 0) {
	store = new ELEMENT[s1*s2];
    }
    size1 = s1;
    size2 = s2;
    off1 = o1;
    off2 = o2;

//    if (store != 0 && size1*size2 == 0)
//      error("ARRAY2::make Internal inconsistency on leaving make");
}

// Copy constructor

ARRAY2::ARRAY2(const ARRAY2& V)
{
    store = 0;
    make(V.off1,V.size1,V.off2,V.size2);
    for (int i = 0 ; i < size1*size2 ; i++)
      store[i] = V.store[i];
}

// Assignment

void ARRAY2::operator=(const ARRAY2& V)
{
    make(V.off1,V.size1,V.off2,V.size2);
    for (int i = 0 ; i < size1*size2 ; i++)
      store[i] = V.store[i];
}

// Assignment, grabbing storage

void ARRAY2::CopyDestroy(ARRAY2& V)
{
    make(0,0,0,0);
    size1 = V.size1; V.size1 = 0;
    size2 = V.size2; V.size2 = 0;
    off1 = V.off1; V.off1 = 0;
    off2 = V.off2; V.off2 = 0;
    store = V.store; V.store = 0;
}

// Assignment to an ELEMENT value

void ARRAY2::operator=(ELEMENT value)
{
    if (store == 0) error("ARRAY2::operator=(ELEMENT): array not dimensioned");
    for (int i = 0 ; i < size1*size2 ; i++)
      store[i] = value;
}

// Output an array

ostream& operator<<(ostream& s, const ARRAY2& V)
{
    long f = s.flags();

    int k = 0;
    for (int i = 0 ; i < V.size1 ; i++) {
	for (int j = 0 ; j < V.size2 ; j++)
	  s << V.store[k++] << " ";
	s << "\n";
    }

    s << endl;  s.flags(f); return s;
}

