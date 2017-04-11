#ifndef __ITENSOR_PARALLEL_H
#define __ITENSOR_PARALLEL_H
#include "mpi.h"
#include <sstream>
#include <vector>
#include <type_traits>

#define DEFAULT_BUFSIZE 500000

namespace itensor {

//
// doRead is a templated function which uses the read method
// of an object to read it from a stream, unless the object
// is plain-old-data (pod) in which case it calls the stream's
// read method directly.
//
// Usage: 
//  XType x;
//  doRead(some_input_stream,x);
//
template<typename T, bool isPod = std::is_pod<T>::value>
struct DoRead
    {
    DoRead(std::istream& data, T& obj)
        {
        obj.read(data);
        }
    };
template<typename T>
struct DoRead<T, true>
    {
    DoRead(std::istream& data, T& val)
        {
        data.read((char*) &val, sizeof(val));
        }
    };
template<typename T>
void
doRead(std::istream& data, T& val)
    {
    DoRead<T>(data,val);
    }

//
// doWrite is a templated functionwhich uses the write method
// of an object to write it from a stream, unless the object
// is plain-old-data (pod) in which case it calls the stream's
// write method directly.
//
// Usage: 
//  XType x = ...;
//  doWrite(some_output_stream,x);
//
template<typename T, bool isPod = std::is_pod<T>::value>
struct DoWrite
    {
    DoWrite(std::ostream& data, const T& obj)
        {
        obj.write(data);
        }
    };
template<typename T>
struct DoWrite<T, true>
    {
    DoWrite(std::ostream& data, const T& val)
        {
        data.write((char*) &val, sizeof(val));
        }
    };
template<typename T>
void
doWrite(std::ostream& data, const T& val)
    {
    DoWrite<T>(data,val);
    }

class Environment
    {
    int rank_,
        nnodes_;
    mutable std::vector<char> buffer;
    public:

    Environment(int argc, char* argv[],
                const Args& args = Global::args());

    ~Environment() { MPI_Finalize(); }

    //
    // Information about nodes
    //

    int 
    rank() const { return rank_; }

    int 
    nnodes() const { return nnodes_; }

    bool 
    firstNode() const { return rank_ == 0; }
    bool 
    lastNode() const { return rank_ == nnodes_-1; }

    int 
    lnode() const { return (rank_ == 0 ? nnodes_-1 : rank_-1); }
    int 
    rnode() const { return (rank_ == nnodes_-1 ? 0 : rank_+1); }

    //
    // Communication and flow control
    //

    template <class T>
    void 
    broadcast(T& obj) const;
    void 
    broadcast(std::stringstream& data) const;

    void 
    barrier() const { MPI_Barrier(MPI_COMM_WORLD); }

    void 
    abort(int code) const { MPI_Abort(MPI_COMM_WORLD,code); }

    double 
    sum(double r) const;

    private:

    //Make this class non-copyable
    void operator=(const Environment&);
    Environment(const Environment&);

    };

class MailBox
    {
    const Environment* env_;
    MPI_Comm com;
    int other_node_;
    char flag_;
    MPI_Request req_;
    MPI_Status rstatus_;
    std::string sdata;
    std::vector<char> rbuffer;
    int tag_;
    public:

    MailBox();

    MailBox(const Environment& env, 
            int other_node,
            const Args& args = Global::args());

    ~MailBox();

    //
    // Accessor methods
    //

    explicit operator bool() const { return env_ != nullptr; }

    int 
    rank() const { checkValid(); return env_->rank(); }

    int 
    nnodes() const { checkValid(); return env_->nnodes(); }

    int
    tag() const { return tag_; }

    //
    // Communication methods
    //

    template <class T>
    void
    receive(T& obj);
    template <class T, typename... Args>
    T
    receive(Args&&... args);
    void 
    receive(std::stringstream& data);

    template <class T> void 
    send(const T& obj);
    void 
    send(const std::stringstream& data);

    template <class T> 
    void 
    broadcast(T& obj) const { checkValid(); env_->broadcast(obj); }

    private:

    void
    checkValid() const
        {
        if(env_==nullptr) throw std::runtime_error("MailBox object is default initialized.");
        }

    void 
    listenForFlag()
        { 
        MPI_Irecv(&flag_,1,MPI_CHAR,other_node_,flagTag(),com,&req_); 
        }

    static int new_tag(const Environment& env, int other_node)
        {
        static std::vector<int> tag(env.nnodes(),0);
        tag.at(other_node) += 3;
        return tag.at(other_node);
        }

    int
    flagTag() const { return tag_+2; }

    int
    sizeTag() const { return tag_+1; }

    }; //class MailBox

inline Environment::
Environment(int argc, char* argv[],
            const Args& args)
    : buffer(args.getInt("Bufsize",DEFAULT_BUFSIZE))
    { 
    MPI_Init(&argc,&argv); 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_); 
    MPI_Comm_size(MPI_COMM_WORLD,&nnodes_); 
    }

template <class T>
void inline Environment::
broadcast(T& obj) const
    {
    if(nnodes_ == 1) return;
    const int root = 0;
    std::stringstream datastream;
    if(rank_ == root) doWrite(datastream,obj);
    broadcast(datastream);
    if(rank_ != root) doRead(datastream,obj);
    }

void inline Environment::
broadcast(std::stringstream& data) const
    { 
    if(nnodes_ == 1) return;
    const int root = 0;
    const int shift = 2;
    if(rank_ == root)
        {
        int size = data.str().length(); 
        int quo = size/buffer.size(), 
            rem = size%buffer.size();
        MPI_Bcast(&quo,1,MPI_INT,root,MPI_COMM_WORLD);
        MPI_Bcast(&rem,1,MPI_INT,root,MPI_COMM_WORLD);
        std::cout << "Doing broadcast (quo,rem) = (" << quo << "," << rem << ")" << std::endl;
        for(int q = 0; q < quo; ++q)
            { 
            MPI_Bcast(const_cast<char*>(data.str().data())+q*buffer.size(),buffer.size(),MPI_CHAR,root,MPI_COMM_WORLD); 
            }
        MPI_Bcast(const_cast<char*>(data.str().data())+quo*buffer.size(),rem+shift,MPI_CHAR,root,MPI_COMM_WORLD);
        }
    else
        {
        int quo=0,rem=0;
        MPI_Bcast(&quo,1,MPI_INT,root,MPI_COMM_WORLD);
        MPI_Bcast(&rem,1,MPI_INT,root,MPI_COMM_WORLD);
        for(int q = 0; q < quo; ++q)
            { 
            MPI_Bcast(&buffer.front(),buffer.size(),MPI_CHAR,root,MPI_COMM_WORLD); 
            data.write(&buffer.front(),buffer.size());
            }
        MPI_Bcast(&buffer.front(),rem+shift,MPI_CHAR,root,MPI_COMM_WORLD);
        data.write(&buffer.front(),rem+shift);
        }
    }

double inline Environment::
sum(double r) const
    {
    if(nnodes_ == 1) return r;
    double res = 0;
    MPI_Reduce(&r,&res,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    return res;
    }


//
// MailBox
//


inline MailBox::
MailBox()
    :
    env_(nullptr)
    {
    }

inline MailBox::
MailBox(const Environment& env, 
        int other_node,
        const Args& args)
    : 
    env_(&env), 
    com(MPI_COMM_WORLD), 
    other_node_(other_node), 
    flag_('f'),
    rbuffer(args.getInt("Bufsize",DEFAULT_BUFSIZE)),
    tag_(new_tag(env,other_node))
    { 
    if(other_node_ >= env_->nnodes())
        { 
        std::cout << "\n\nNode " << env_->rank() << ": other_node = " << other_node_ << " out of range." << std::endl;
        throw std::runtime_error("other_node out of range"); 
        }
    //Initiate flag receive request
    listenForFlag();
    }
inline MailBox::
~MailBox() 
    { 
    if(env_) MPI_Cancel(&req_); 
    }

void inline MailBox::
receive(std::stringstream& data)
    {
    checkValid();
    MPI_Wait(&req_,&rstatus_); 

    int msize = 0;
    MPI_Recv(&msize,1,MPI_INT,other_node_,sizeTag(),com,&rstatus_);

    int quo = msize/rbuffer.size(), 
        rem = msize%rbuffer.size();
    for(int q = 0; q < quo; ++q)
        {
        MPI_Recv(&rbuffer.front(),rbuffer.size(),MPI_CHAR,other_node_,tag(),com,&rstatus_);
        data.write(&rbuffer.front(),rbuffer.size());
        }
    MPI_Recv(&rbuffer.front(),rem,MPI_CHAR,other_node_,tag(),com,&rstatus_);
    data.write(&rbuffer.front(),rem);

    //Reset flag_
    listenForFlag();
    }

template <class T>
void MailBox::
receive(T& obj)
    { 
    std::stringstream data; 
    receive(data); 
    doRead(data,obj);
    }

template <class T, typename... Args>
T MailBox::
receive(Args&&... args)
    { 
    std::stringstream data; 
    receive(data); 
    T obj(std::forward<Args>(args)...);
    doRead(data,obj);
    return obj;
    }


void inline MailBox::
send(const std::stringstream& data)
    {
    checkValid();
    sdata.assign(data.str());
    int msize = sdata.length();
    int quo = msize/rbuffer.size(), 
        rem = msize%rbuffer.size();

    MPI_Send(&flag_,1,MPI_CHAR,other_node_,flagTag(),com);
    MPI_Send(&msize,1,MPI_INT,other_node_,sizeTag(),com);

    auto datap = const_cast<char*>(sdata.data());
    for(int q = 0; q < quo; ++q)
        {
        MPI_Send(datap+q*rbuffer.size(),rbuffer.size(),MPI_CHAR,other_node_,tag(),com);
        }
    MPI_Send(datap+quo*rbuffer.size(),rem,MPI_CHAR,other_node_,tag(),com);
    }


template <class T> 
void inline MailBox::
send(const T& obj)
    {
    std::stringstream data; 
    doWrite(data,obj);
    send(data); 
    }

} //namespace itensor

#endif
