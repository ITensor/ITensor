//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//

// hash.h -- Steve's hash table class

#ifndef _hash_h
#define _hash_h 1

typedef double Real;

class Hash                 	// Hash table
{
	void copy(const Hash &);
	void destroy();
public:
	~Hash() {destroy();}
	Hash(const Hash &other)
		{copy(other);}
	Hash& operator=(const Hash &other)
		{
		if(&other != this)
			{ destroy(); copy(other); }
		return *this;
		}

	int length;
	int* index;
	int* values;

	Hash() {length = 0; values = 0; index = 0;}

	Hash(int l);

        void Resize(int l);

        int Length() { return length; }
    
	void putin(int,int);

	int get(int, int&) const;

	int remove(int);

	int bytesused() const {return 12 + length * 8;}
};

class FirstNext			// Uses Hash to locate things in a list
{
	void copy(const FirstNext &);
	void destroy();
public:
	~FirstNext() {destroy();}
	FirstNext(const FirstNext &other)
		{copy(other);}
	FirstNext& operator=(const FirstNext &other)
		{
		if(&other != this)
			{ destroy(); copy(other); }
		return *this;
		}

	int length;
	Hash first;
	int* next;
	int currentindex;

	FirstNext() {length = 0; next = 0; currentindex = -1;}

	FirstNext(int l);

	void putin(int,int);

	int getfirst(int, int&);

	int getnext(int&);

	int get(int, int&, int&);

	int bytesused() const {return first.bytesused() + 12 + length*4;}
};

#endif
