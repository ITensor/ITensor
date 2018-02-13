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

class Environment;

void
parallelDebugWait(Environment const& env);

template <class T>
void
broadcast(Environment const& env, T & obj);

template <class T, class... Rest>
void 
broadcast(Environment const& env, T & obj, Rest &... rest);

template <typename T>
void 
scatterVector(Environment const& env, std::vector<T> &v);

double
sum(Environment const& env, double r);

template <typename T>
T
sum(Environment const& env, T &obj);

template <typename T>
T
allSum(Environment const& env, T &obj);

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

    void 
    broadcast(std::stringstream& data) const;

    void 
    barrier() const { MPI_Barrier(MPI_COMM_WORLD); }

    void 
    abort(int code) const { MPI_Abort(MPI_COMM_WORLD,code); }

    private:

    //Make this class non-copyable
    void operator=(Environment const&);
    Environment(Environment const&);

    public:

    //Deprecated methods: kept here for backwards compatibility

    template <class T>
    void 
    broadcast(T& obj) const;

    template <class T, class... Rest>
    void 
    broadcast(T & obj, Rest &... rest) const;

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

void inline
parallelDebugWait(Environment const& env)
    {
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    for(int n = 0; n < env.nnodes(); ++n)
        {
        if(env.rank() == n) printfln("Node %d PID %d on %s",env.rank(),getpid(),hostname);
        env.barrier();
        }
#ifdef DEBUG
    //
    // If compiled in debug mode, require file GO to exist before
    // proceeding. This is to make it possible to attach a debugger
    // such as gdb to each thread before the calculation starts.
    //
    if(env.nnodes() > 1)
        {
        while(!fileExists("GO")) sleep(2);
        printfln("Process %d found file GO, exiting wait loop",env.rank());
        }
    env.barrier();
    if(env.firstNode()) system("rm -f GO");
#endif
    } 

inline Environment::
Environment(int argc, char* argv[],
            Args const& args)
    : buffer(args.getInt("Bufsize",DEFAULT_BUFSIZE))
    { 
    MPI_Init(&argc,&argv); 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_); 
    MPI_Comm_size(MPI_COMM_WORLD,&nnodes_); 
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
        //std::cout << "Doing broadcast (quo,rem) = (" << quo << "," << rem << ")" << std::endl;
        for(int q = 0; q < quo; ++q)
            { 
            //printfln("Node %d q = %d/%d",rank_,q,(quo-1));
            MPI_Bcast(const_cast<char*>(data.str().data())+q*buffer.size(),buffer.size(),MPI_CHAR,root,MPI_COMM_WORLD); 
            }
        MPI_Bcast(const_cast<char*>(data.str().data())+quo*buffer.size(),rem+shift,MPI_CHAR,root,MPI_COMM_WORLD);
        }
    else
        {
        int quo=0,rem=0;
        MPI_Bcast(&quo,1,MPI_INT,root,MPI_COMM_WORLD);
        MPI_Bcast(&rem,1,MPI_INT,root,MPI_COMM_WORLD);
        //printfln("Node %d got quo,rem = %d,%d",rank_,quo,rem);
        for(int q = 0; q < quo; ++q)
            { 
            //printfln("Node %d q = %d/%d",rank_,q,(quo-1));
            MPI_Bcast(&buffer.front(),buffer.size(),MPI_CHAR,root,MPI_COMM_WORLD); 
            data.write(&buffer.front(),buffer.size());
            }
        MPI_Bcast(&buffer.front(),rem+shift,MPI_CHAR,root,MPI_COMM_WORLD);
        data.write(&buffer.front(),rem+shift);
        }
    //printfln("%d reached end of broadcast(stringstream)",rank_);
    }

template <class T>
void
broadcast(Environment const& env, T & obj)
    {
    if(env.nnodes() == 1) return;
    const int root = 0;
    std::stringstream datastream;
    if(env.rank() == root) itensor::write(datastream,obj);
    env.broadcast(datastream);
    if(env.rank() != root) itensor::read(datastream,obj);
    }

template <class T, class... Rest>
void 
broadcast(Environment const& env, T & obj, Rest &... rest)
    {
    broadcast(env,obj);
    broadcast(env,rest...);
    }

template <class T>
void Environment::
broadcast(T& obj) const 
    { 
    itensor::broadcast<T>(*this,obj); 
    }

template <class T, class... Rest>
void Environment::
broadcast(T & obj, Rest &... rest) const
    { 
    itensor::broadcast(*this,obj,rest...); 
    }

template <typename T>
void 
scatterVector(Environment const& env, std::vector<T> &v)
    {
    if(env.nnodes() == 1) return;
    const int root = 0;
    

    if(env.firstNode())
        { 
        auto nnodes = env.nnodes();
        auto n = v.size();
        auto blockSizes = std::vector<long>(nnodes);
        long blockSize = n / nnodes;

        for(int i = 0; i < nnodes; i++) blockSizes[i] = blockSize;

        if (n % nnodes != 0)
            {
            for(int i = 0; i < (n % nnodes); i++) ++blockSizes[i];
            }

        long mySize = 0l;
        MPI_Scatter(blockSizes.data(),1,MPI_LONG,&mySize,1,MPI_LONG,root,MPI_COMM_WORLD);

        auto itp = blockSizes[0];
        for (int i = 1; i < nnodes; ++i)
            {
            MailBox mailbox(env,i);
            mailbox.send(std::vector<T>(v.begin()+itp,v.begin()+itp+blockSizes[i]));
            itp += blockSizes[i];
            }
        v.resize(mySize);
        }
    else
        {
        long mySize = 0l;
        MPI_Scatter(NULL,1,MPI_LONG,&mySize,1,MPI_LONG,root,MPI_COMM_WORLD);   
        v.resize(mySize);
        MailBox mailbox(env,root);
        mailbox.receive(v);
        }
    }

template <typename T>
void 
gatherVector(Environment const& env, std::vector<T> &v)
    { 
    if(env.nnodes() == 1) return;
    const int root = 0;
    
    if(env.firstNode())
        { 
        for(int i = 1; i < env.nnodes(); i++)
            {
            MailBox mailbox(env,i);
            std::vector<T> tmp;
            mailbox.receive(tmp);
            v.insert(v.end(), tmp.begin(), tmp.end());
            }
        }
    else
        {
        MailBox mailbox(env,root);
        mailbox.send(v);
        }
    }

double inline
sum(Environment const& env, double r)
    {
    if(env.nnodes() == 1) return r;
    double res = 0;
    MPI_Reduce(&r,&res,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    return res;
    }

template <typename T>
T
sum(Environment const& env, T &obj)
    {
    if(env.nnodes() == 1) return obj;
    const int root = 0;

    T res = obj;
    if(env.rank() == 0)
        {
        for (int i = 1; i < env.nnodes(); ++i)
            {
            MailBox mailbox(env,i);
            T tmp;
            mailbox.receive(tmp);
            res += tmp;
            }
        }
    else
        {
        MailBox mailbox(env,root);
        mailbox.send(obj);
        }
    return res;
    }

template <typename T>
T
allSum(Environment const& env, T &obj)
    {
    if(env.nnodes() == 1) return obj;
    T result = sum(env,obj);
    broadcast(env,result);
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
    itensor::read(data,obj);
    }

template <class T, typename... Args>
T MailBox::
receive(Args&&... args)
    { 
    std::stringstream data; 
    receive(data); 
    T obj(std::forward<Args>(args)...);
    itensor::read(data,obj);
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
    itensor::write(data,obj);
    send(data); 
    }

} //namespace itensor

#endif
