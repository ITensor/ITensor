//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "index.h"

namespace itensor {

using std::string;
using std::stringstream;


struct TInfo
    {
    IndexType t = NullIndex;
    const char* s = "";
    int n = 0;
    TInfo(IndexType t_, const char* s_) : t(t_), s(s_), n(0) { }
    TInfo(IndexType t_, const char* s_, int n_) : t(t_), s(s_), n(n_) { }
    };

static constexpr int NIType = 13;

using TInfoArray = std::array<TInfo,NIType>;

#define REGISTER_ITYPE(X) TInfo(X,#X)

TInfoArray&
tinfo()
    {
    auto makeTInfoArr = []()
        {
        TInfoArray a =
            {{
            REGISTER_ITYPE(Link),
            REGISTER_ITYPE(Site),
            REGISTER_ITYPE(All),
            REGISTER_ITYPE(NullIndex),
            REGISTER_ITYPE(Atype),
            REGISTER_ITYPE(Btype),
            REGISTER_ITYPE(Ctype),
            REGISTER_ITYPE(Dtype),
            REGISTER_ITYPE(Xtype),
            REGISTER_ITYPE(Ytype),
            REGISTER_ITYPE(Ztype),
            REGISTER_ITYPE(Wtype),
            REGISTER_ITYPE(Vtype) 
            }};
        for(size_t j = 1; j <= a.size(); ++j)
            a[j-1].n = j;
        return a;
        };
    static auto a = makeTInfoArr();
    return a;
    }

int 
IndexTypeToInt(IndexType it)
    {
    for(auto& el : tinfo())
        if(el.t == it)
            {
            return el.n;
            }
    Error("IndexTypeToInt: IndexType not recognized");
    return -1;
    }

IndexType 
IntToIndexType(int i)
    {
    for(auto& el : tinfo())
        if(el.n == i)
            {
            return el.t;
            }
    printfln("No IndexType value defined for i=%d\n",i);
    Error("Undefined IntToIndexType value");
    return Link;
    }

const char*
indexTypeName(IndexType it)
    {
    for(auto& el : tinfo())
        if(el.t == it)
            {
            return el.s;
            }
    Error("Undefined indexTypeName");
    return "";
    }

std::ostream& 
operator<<(std::ostream& s, IndexType it)
    { 
    s << indexTypeName(it);
    return s; 
    }

string 
putprimes(string s, int plev)
    { 
    stringstream str;
    str << s;
    if(plev < 0) Error("Negative prime level");
    if(plev > 3)
        {
        //str << "[" << plev << "']";
        str << "^" << plev;
        }
    else
        {
        for(int i = 1; i <= plev; ++i) str << "\'";
        }
    return str.str();
    }

string 
nameindex(IndexType it, int plev)
    { 
    return putprimes(indexTypeName(it),plev); 
    }

string 
nameint(const string& f, int n)
    { 
    stringstream ss; 
    ss << f << n; 
    return ss.str(); 
    }

//
// class Index
//


string Index::
name() const  { return putprimes(sname_,primelevel_); }

Index& Index::
primeLevel(int plev) 
    { 
    primelevel_ = plev; 
#ifdef DEBUG
    if(primelevel_ < 0)
        Error("Negative primeLevel");
#endif
    return *this;
    }


Index& Index::
mapprime(int plevold, int plevnew, IndexType type)
    {
    if(primelevel_ == plevold)
        {
        if(type == All || type == this->type())
            {
            primelevel_ = plevnew;
#ifdef DEBUG
            if(primelevel_ < 0)
                {
                Error("Negative primeLevel");
                }
#endif
            }
        }
    return *this;
    }


Index& Index::
prime(IndexType type, int inc)
    {
    if(type == this->type() || type == All)
        {
        primelevel_ += inc;
#ifdef DEBUG
        if(primelevel_ < 0)
            {
            Error("Increment led to negative primeLevel");
            }
#endif
        }
    return *this;
    }

void Index::
write(std::ostream& s) const 
    { 
    if(!bool(*this)) Error("Index::write: Index is default initialized");

    s.write((char*) &primelevel_,sizeof(primelevel_));

    const auto t = IndexTypeToInt(type_);
    s.write((char*) &t,sizeof(t));

    s.write((char*) &(id_),sizeof(id_));

    s.write((char*) &(m_),sizeof(m_));

    const int nlength = sname_.length();
    s.write((char*) &nlength,sizeof(nlength));

    s.write(sname_.c_str(),nlength+1);
    }

Index& Index::
read(std::istream& s)
    {
    s.read((char*) &primelevel_,sizeof(primelevel_));
#ifdef DEBUG
    if(primelevel_ < 0)
        {
        Error("Negative primeLevel");
        }
#endif

    int t; 
    s.read((char*) &t,sizeof(t));
    type_ = IntToIndexType(t);

    s.read((char*) &id_, sizeof(id_));

    s.read((char*) &m_,sizeof(m_));

    int nlength; 
    s.read((char*) &nlength,sizeof(nlength));

    auto newname = std::unique_ptr<char[]>(new char[nlength+1]);
    s.read(newname.get(),nlength+1);
    sname_ = string(newname.get()); 

    return *this;
    }

std::ostream& 
operator<<(std::ostream& s, const Index& t)
    {
    if(t.name() != "" && t.name() != " ") s << t.name();
    return s << "(" << nameindex(t.type(),t.primeLevel()) 
             << "," << (t.id() % 10000) << "):" << t.m();
    }

IndexVal::
IndexVal() 
    : 
    val(0) 
    { }

IndexVal::
IndexVal(const Index& index_, long val_) 
    : 
    index(index_),
    val(val_)
    { 
#ifdef DEBUG
    if(!index) Error("IndexVal initialized with default initialized Index");
    //Can also use IndexVal's to indicate prime increments:
    //if(val_ < 1 || val_ > index.m())
    //    {
    //    println("val = ",val_);
    //    println("index = ",index);
    //    Error("val out of range");
    //    }
#endif
    }

bool IndexVal::
operator==(const IndexVal& other) const
    {
    return (index == other.index && val == other.val);
    }


string
showm(const Index& I) { return nameint("m=",I.m()); }

std::ostream& 
operator<<(std::ostream& s, const IndexVal& iv)
    { 
    const Index& ii = iv.index;
    return s << "IndexVal: val = " << iv.val << ", ind = " << ii << "\n"; 
    }

}; //namespace itensor

