#ifndef __ITENSOR_PARALLEL_H
#define __ITENSOR_PARALLEL_H
#include "mpi.h"
#include <sstream>
#include <vector>
#include <type_traits>
#include "itensor/util/readwrite.h"
#include "itensor/util/args.h"

#define DEFAULT_BUFSIZE 500000

namespace itensor {

class Environment
    {
    int rank_,
        nnodes_;
    mutable std::vector<char> buffer;
    public:

    Environment(int argc, char* argv[],
                Args const& args = Args::global());

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

    template <class T, class... Rest>
    void 
    broadcast(T & obj, Rest &... rest) const;

    void 
    broadcast(std::stringstream& data) const;

    template <typename T>
    void
    scatterVector(std::vector<T> &v);

    void 
    barrier() const { MPI_Barrier(MPI_COMM_WORLD); }

    void 
    abort(int code) const { MPI_Abort(MPI_COMM_WORLD,code); }

    double 
    sum(double r) const;

    template <typename T>
    T
    sum(T &obj) const;

    template <typename T>
    T
    allSum(T &obj) const;

    private:

    //Make this class non-copyable
    void operator=(Environment const&);
    Environment(Environment const&);

    };

class MailBox
    {
    Environment const* env_;
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

    MailBox(Environment const& env, 
            int other_node,
            Args const& args = Args::global());

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
    send(T const& obj);
    void 
    send(std::stringstream const& data);

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

    static int new_tag(Environment const& env, int other_node)
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
            Args const& args)
    : buffer(args.getInt("Bufsize",DEFAULT_BUFSIZE))
    { 
    MPI_Init(&argc,&argv); 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_); 
    MPI_Comm_size(MPI_COMM_WORLD,&nnodes_); 
    }

template <class T>
void Environment::
broadcast(T & obj) const
    {
    if(nnodes_ == 1) return;
    const int root = 0;
    std::stringstream datastream;
    if(rank_ == root) write(datastream,obj);
    broadcast(datastream);
    if(rank_ != root) read(datastream,obj);
    }

template <class T, class... Rest>
void Environment::
broadcast(T & obj, Rest &... rest) const
    {
    broadcast(obj);
    broadcast(rest...);
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

template <typename T>
void Environment::
scatterVector(std::vector<T> &v)
    {
    if(nnodes_ == 1) return;
    const int root = 0;
    
    long mySize;

    if(rank() == root) 
        { 
        auto n = v.size();
        long blockSizes[nnodes_];
        long blockSize = n / nnodes_;

        for(int i = 0; i < nnodes_; i++)
            blockSizes[i] = blockSize;
        if (n % nnodes_ != 0)
            for(int i = 0; i < (n % nnodes_); i++)
            blockSizes[i]++;

        MPI_Scatter(blockSizes, 1, MPI_LONG, &mySize, 1, MPI_LONG, root, MPI_COMM_WORLD);

        auto sum = blockSizes[0];
        for (int i = 1; i < nnodes_; ++i)
            {
            MailBox mailbox(*this, i);
            for (int j = 0; j < blockSizes[i]; ++j)
                {
                mailbox.send(v[sum + j]);
                }
            sum += blockSizes[i];
            }
        v.resize(mySize);
        }
    else
        {
        MPI_Scatter(NULL, 1, MPI_LONG, &mySize, 1, MPI_LONG, root, MPI_COMM_WORLD);   
        v.resize(mySize);
        MailBox mailbox(*this, root);
        for (int i = 0; i < mySize; ++i) { mailbox.receive(v[i]); }
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

template <typename T>
T Environment::
sum(T &obj) const
    {
    if(nnodes_ == 1) return obj;
    const int root = 0;

    T result = obj;
    if(rank() == 0)
        {
        for (int i = 1; i < nnodes_; ++i)
            {
            MailBox mailbox(*this, i);
            T temp;
            mailbox.receive(temp);
            result += temp;
            }
        }
    else
        {
        MailBox mailbox(*this, root);
        mailbox.send(obj);
        }
    return result;
    }

template <typename T>
T inline Environment::
allSum(T &obj) const
    {
    if(nnodes_ == 1) return obj;
    T result = sum(obj);
    broadcast(result);
    return result;
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
MailBox(Environment const& env, 
        int other_node,
        Args const& args)
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
    read(data,obj);
    }

template <class T, typename... Args>
T MailBox::
receive(Args&&... args)
    { 
    std::stringstream data; 
    receive(data); 
    T obj(std::forward<Args>(args)...);
    read(data,obj);
    return obj;
    }


void inline MailBox::
send(std::stringstream const& data)
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
send(T const& obj)
    {
    std::stringstream data; 
    write(data,obj);
    send(data); 
    }

} //namespace itensor

#endif
