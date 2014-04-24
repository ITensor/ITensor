//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//

// hash.cc -- Members for Steve's Hash and FirstNext classes

#include "hash.h"
//#include <math>
#include <iostream>


using namespace std;

void error(const string&);

void Hash::copy(const Hash &other)
	{
	length = other.length;
	values = new int[length];
	index = new int[length];
	int i;
	for(i=0; i < length; i++)
		{
		index[i]=other.index[i];
		values[i]=other.values[i];
		}
	}
void Hash::destroy() 
	{
	delete [] values; values = 0;
	delete [] index; index = 0;
	}

Hash::Hash(int l)
	{
	length = l;
	values = new int[length];
	index = new int[length];
	int i;
	for(i=0; i < length; i++)
		index[i]=-1;
	}

// Resize -- added by RMN

void Hash::Resize(int l)
{
    if (l == length) return;
    destroy();
    length = l;
    if (length > 0) {
	values = new int[length];
	index = new int[length];
	int i;
	for(i=0; i < length; i++)
	  index[i]=-1;
    }
}    

/*
inline unsigned hashmap(int key) 
	{
	Real x = (key*2.23457e-3+2.61803398875) * key;
	x = fmod(x,1.0e9);
	return (unsigned) x;
	}
*/
inline unsigned hashmap(int key)
    {
    unsigned int uk = key,r1,r2,r3,r4,r5;
    r1 = 72345219u + uk * 233878732u + 777u * uk * uk + (uk>>8);
    r2 = 123456789u + r1 * 7654321u + (r1>>8);
    r3 = 9234569u + r2 * 12345671u + uk * r1 + (r2>>16);
    r4 = ((uk+r3) >> 8) * 5163151u + r2 * r3 + (r3>>8);
    r5 = ((314159265u + 7u * r1 + r4 * 755332u)>>16) + r3 * uk * r2;
    
    return uk + r1 + r2 + r3 + r4 + r5 + (r5>>16);
    }


void Hash::putin(int key, int val)
	{
	unsigned modkey = hashmap(key) % length;
	int i;
	for(i=0; i < 30; i++)
		{
		if(index[modkey] == -1)		// This space unused.
			{
			index[modkey] = key;
			values[modkey] = val;
			return;
			}
		modkey = (modkey - i - 1 + length) % length;
		}
	cout << endl << "length is " << length << endl;
	cout << "key, original hashmap, original modkey, modkey, val are "
	     << key << " " << hashmap(key) << " " << hashmap(key)%length
	     << " " << modkey << " " << val << endl;
	for(i=0; i < length; i++)
		cout << i << " " << index[i] << " " << hashmap(index[i]) 
		     << " " << hashmap(index[i])%length << " " 
		     << values[i] << endl;
	error("trouble in hash: putin");
	}

int Hash::get(int key, int& val) const
	{
	unsigned modkey = hashmap(key) % length;
	int i;
	for(i=0; i < 30; i++)
		{
		if(index[modkey] == key)
			{ val = values[modkey];  return 1; }
		if(index[modkey] == -1) return 0;
		modkey = (modkey - i - 1 + length) % length;
		}
	cout << "Warning in hash::get " << endl;
	for(i=21; i < 200; i++)
		{
		cout << i << endl;
		if(index[modkey] == key)
			{ val = values[modkey];  return 1; }
		if(index[modkey] == -1) return 0;
		modkey = (modkey - i - 1 + length) % length;
		}
	error("deadly trouble in hash: get");
	return 0;
	}

int Hash::remove(int key)
	{
	unsigned modkey = hashmap(key) % length;
	int i;
	for(i=0; i < 30; i++)
		{
		if(index[modkey] == key)
			{ index[modkey] = -1; return 1; }
		if(index[modkey] == -1) return 0;
		modkey = (modkey - i - 1 + length) % length;
		}
	error("trouble in hash: remove");
	return 0;
	}

void FirstNext::copy(const FirstNext &other)
	{
	length = other.length;
	next = new int[length];
	int i;
	for(i=0; i < length; i++)
		next[i]=other.next[i];
	first = other.first;
	currentindex = other.currentindex;
	}
void FirstNext::destroy() 
	{
	delete [] next;
	}

FirstNext::FirstNext(int l)
	{
	length = l;
	next = new int[length];
	first = Hash(3*length); 
	int i;
	for(i=0; i < length; i++)
		next[i]=-1;
	currentindex = -1;
	}

void FirstNext::putin(int key, int index)
	{
	int j;
	if(first.get(key,j))  //j is now the index of first identical key
		{
	// trace way through previous occurences of key
		while(next[j] != -1) 
			j = next[j];
		next[j] = index;
		}
	else
		{ first.putin(key,index); }
	}

int FirstNext::getfirst(int key, int& index)
	{
	if(first.get(key,index))
		{
		currentindex = index;
		return 1;
		}
	else
		{
		currentindex = -1;
		return 0;
		}
	}

int FirstNext::getnext(int& index)
	{
	if(currentindex == -1) return 0;
	index = currentindex;
	currentindex = next[currentindex];
	if(currentindex != -1)
		{
		index = currentindex;
		return 1;
		}
	else
		return 0;
	}

int FirstNext::get(int key, int& index, int& newkey)
	{
	if(newkey == 1)
		{
		newkey = 0;
		return getfirst(key,index);
		}
	return getnext(index);
	}
