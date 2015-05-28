//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_READWRITE_H_
#define __ITENSOR_READWRITE_H_

#include <fstream>
#include <unistd.h>
#include "itensor/types.h"
#include "itensor/util/error.h"

namespace itensor {

bool inline
fileExists(const std::string& fname)
    {
    std::ifstream file(fname.c_str());
    return file.good();
    }

// Can overload read and write functions with
// signatures below for objects such as std::vector
// 
// For classes having member read/write functions, can 
// leave external read/write overloads undefined
// and the following template overloads will 
// be called
//

//Here we have to use a struct to implement the read(istream,T)
//function because function templates cannot be partially specialized
template<typename T, bool isPod = std::is_pod<T>::value>
struct DoRead
    {
    DoRead(std::istream& s, T& obj)
        {
        obj.read(s);
        }
    };
template<typename T>
struct DoRead<T, true>
    {
    DoRead(std::istream& s, T& val)
        {
        s.read((char*) &val, sizeof(val));
        }
    };
template<typename T>
void
read(std::istream& s, T& val)
    {
    DoRead<T>(s,val);
    }

template<typename T, typename... CtrArgs>
T
read(std::istream& s, CtrArgs&&... args)
    {
    T t(std::forward<CtrArgs>(args)...);
    DoRead<T>(s,t);
    return t;
    }

template<typename T, bool isPod = std::is_pod<T>::value>
struct DoWrite
    {
    DoWrite(std::ostream& s, const T& obj)
        {
        obj.write(s);
        }
    };
template<typename T>
struct DoWrite<T, true>
    {
    DoWrite(std::ostream& s, const T& val)
        {
        s.write((char*) &val, sizeof(val));
        }
    };
template<typename T>
void
write(std::ostream& s, const T& val)
    {
    DoWrite<T>(s,val);
    }

template<typename T>
void
write(std::ostream& s, const std::vector<T>& vec)
    {
    auto size = vec.size();
    s.write((char*)&size,sizeof(size));
    s.write((char*)vec.data(),sizeof(T)*size);
    }

template<typename T>
void
read(std::istream& s, std::vector<T>& vec)
    {
    auto size = vec.size(); //will overwrite
    s.read((char*)&size,sizeof(size));
    vec.resize(size);
    s.read((char*)vec.data(),sizeof(T)*size);
    }

void inline
read(std::istream& s, Cplx& z)
    {
    auto &r = reinterpret_cast<Real(&)[2]>(z)[0];
    auto &i = reinterpret_cast<Real(&)[2]>(z)[1];
    s.read((char*)&r,sizeof(r));
    s.read((char*)&i,sizeof(i));
    }

void inline
write(std::ostream& s, const Cplx& z)
    {
    auto &r = reinterpret_cast<const Real(&)[2]>(z)[0];
    auto &i = reinterpret_cast<const Real(&)[2]>(z)[1];
    s.write((char*)&r,sizeof(r));
    s.write((char*)&i,sizeof(i));
    }

template<class T> 
void
readFromFile(const std::string& fname, T& t) 
    { 
    std::ifstream s(fname.c_str(),std::ios::binary);
    if(!s.good()) 
        throw ITError("Couldn't open file \"" + fname + "\" for reading");
    read(s,t); 
    s.close(); 
    }


template<class T, typename... InitArgs>
T
readFromFile(const std::string& fname, InitArgs&&... iargs)
    { 
    std::ifstream s(fname.c_str(),std::ios::binary); 
    if(!s.good()) 
        throw ITError("Couldn't open file \"" + fname + "\" for reading");
    T t(std::forward<InitArgs>(iargs)...);
    read(s,t); 
    s.close(); 
    return t;
    }


template<class T> 
void
writeToFile(const std::string& fname, const T& t) 
    { 
    std::ofstream s(fname.c_str(),std::ios::binary); 
    if(!s.good()) 
        throw ITError("Couldn't open file \"" + fname + "\" for writing");
    write(s,t); 
    s.close(); 
    }

//Given a prefix (e.g. pfix == "mydir")
//and an optional location (e.g. locn == "/var/tmp/")
//creates a temporary directory and returns its name
//without a trailing slash
//(e.g. /var/tmp/mydir_SfqPyR)
std::string inline
mkTempDir(const std::string& pfix,
          const std::string& locn = "./")
    {
    //Construct dirname
    std::string dirname = locn;
    if(dirname[dirname.length()-1] != '/')
        dirname += '/';
    //Add prefix and template string of X's for mkdtemp
    dirname += pfix + "_XXXXXX";

    //Create C string version of dirname
    auto cstr = std::unique_ptr<char[]>(new char[dirname.size()+1]);
    strcpy(cstr.get(),dirname.c_str());

    //Call mkdtemp
    char* retval = mkdtemp(cstr.get());
    //Check error condition
    if(retval == NULL) throw ITError("mkTempDir failed");

    //Prepare return value
    std::string final_dirname(retval);

    return final_dirname;
    }

} // namespace itensor

#endif
