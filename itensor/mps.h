#ifndef __MPS_H
#define __MPS_H
#include "iqcombiner.h"
#include "model.h"

enum Direction { Fromright, Fromleft, Both, None };

static const Real DefaultLogRef = 1.3E-14;

Vector do_denmat_Real(const ITensor& nA, ITensor& A, ITensor& B, Real cutoff,int minm, int maxm, Direction dir);
Vector do_denmat_Real(const IQTensor& nA, IQTensor& A, IQTensor& B, Real cutoff, int minm,int maxm, Direction dir);
Vector do_denmat_Real(const vector<IQTensor>& nA, const IQTensor& A, const IQTensor& B, IQTensor& U,
	Real cutoff,Real lognormfac,int minm,int maxm, Direction dir, bool donormalize, bool do_relative_cutoff);

template<class Tensor, class TensorSet>
Real doDavidson(Tensor& phi, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH, int niter, int debuglevel, Real errgoal);

extern Real truncerror, svdtruncerr;

template<class Tensor, class IndexT>
IndexT index_in_common(const Tensor& A, const Tensor& B, IndexType t)
{
    for(int j = 1; j <= A.r(); ++j)
    {
        const IndexT& I = A.index(j);
        if(I.type() == t && B.hasindex(I)) { return I; }
    }
    return IndexT();
}
inline Index index_in_common(const ITensor& A, const ITensor& B, IndexType t)
{ return index_in_common<ITensor,Index>(A,B,t); }
inline IQIndex index_in_common(const IQTensor& A, const IQTensor& B, IndexType t)
{ return index_in_common<IQTensor,IQIndex>(A,B,t); }

namespace {
int collapseCols(const Vector& Diag, Matrix& M)
{
    int nr = Diag.Length(), nc = int(Diag.sumels());
    assert(nr != 0);
    if(nc == 0) return nc;
    M = Matrix(nr,nc); M = 0;
    int c = 0;
    for(int r = 1; r <= nr; ++r)
    if(Diag(r) == 1) { M(r,++c) = 1; }
    return nc;
}
}

inline void convertToIQ(const BaseModel& model, const vector<ITensor>& A, vector<IQTensor>& qA, QN totalq = QN(), Real cut = 1E-12)
{
    const int N = A.size()-1;
    qA.resize(A.size());
    const bool is_mpo = A[1].hasindex(model.si(1).primed());
    const int Dim = model.dim();
    const int PDim = (is_mpo ? Dim : 1);

    vector<IQIndex> linkind(N);

    typedef map<QN,Vector>::value_type qD_vt;
    map<QN,Vector> qD; //Diags of compressor matrices by QN
    typedef map<QN,vector<ITensor> >::value_type qt_vt;
    map<QN,vector<ITensor> > qt; //ITensor blocks by QN
    typedef map<QN,ITensor>::value_type qC_vt;
    map<QN,ITensor> qC; //Compressor ITensors by QN
    ITensor block;
    vector<ITensor> nblock;
    vector<inqn> iq;

    QN q;

    qC[totalq] = ITensor(); //Represents Virtual index
    //First value of prev_q below set to totalq

    const int show_s = 0;

    Index bond, prev_bond;
    for(int s = 1; s <= N; ++s)
    {
        qD.clear(); qt.clear();
        if(s > 1) prev_bond = index_in_common(A[s-1],A[s],Link);
        if(s < N) bond = index_in_common(A[s],A[s+1],Link);

        if(s == show_s) { PrintDat(A[s]); }

        foreach(const qC_vt& x, qC) {
        const QN& prev_q = x.first; const ITensor& comp = x.second; 
        for(int n = 1; n <= Dim;  ++n)
        for(int u = 1; u <= PDim; ++u)
        {
            q = (is_mpo ? prev_q+model.si(s).qn(n)-model.si(s).qn(u) : prev_q-model.si(s).qn(n));

            //For the last site, only keep blocks 
            //compatible with specified totalq i.e. q=0 here
            if(s == N && q != QN()) continue;

            //Set Site indices of A[s] and its previous Link Index
            block = A[s];
            if(s != 1) block *= conj(comp);
            block *= model.si(s)(n);
            if(is_mpo) block *= model.siP(s)(u);

            //Initialize D Vector (D records which values of
            //the right Link Index to keep for the current QN q)
            int count = qD.count(q);
            Vector& D = qD[q];
            if(count == 0) { D.ReDimension(bond.m()); D = 0; }

            if(s == show_s)
            {
                cerr << format("For n = %d\n")%n;
                cerr << format("Got a block with norm %.10f\n")%block.norm();
                cerr << format("bond.m() = %d\n")%bond.m();
                PrintDat(block);
                if(s != 1) PrintDat(comp);
            }

            bool keep_block = false;
            if(s == N) keep_block = true;
            else
            {
                if(bond.m() == 1 && block.norm() != 0) { D = 1; keep_block = true; }
                else
                {
                    ITensor summed_block;
                    if(s==1) summed_block = block;
                    else
                    {
                        //Here we sum over the previous link index
                        //which is already ok, analyze the one to the right
                        assert(comp.r()==2);
                        Index new_ind = (comp.index(1)==prev_bond ? comp.index(2) : comp.index(1));
                        summed_block = ITensor(new_ind,1) * block;
                    }
                    //cerr << format("s = %d, bond=")%s << bond << "\n";
                    //summed_block.print("summed_block");

                    Real rel_cut = -1;
                    for(int j = 1; j <= bond.m(); ++j)
                    { rel_cut = max(fabs(summed_block.val1(j)*exp(summed_block.logfac())),rel_cut); }
                    assert(rel_cut >= 0);
                    //Real rel_cut = summed_block.norm()/summed_block.vec_size();
                    rel_cut *= cut;
                    //cerr << "rel_cut == " << rel_cut << "\n";

                    if(rel_cut > 0)
                    for(int j = 1; j <= bond.m(); ++j)
                    if(fabs(summed_block.val1(j)*exp(summed_block.logfac())) > rel_cut) 
                    { D(j) = 1; keep_block = true; }
                }
            } //else (s != N)

            //if(!keep_block && q == totalq) { D(1) = 1; keep_block = true; }

            if(keep_block)
            {
                qD[q] = D;

                if(is_mpo) 
                {
                block.addindex1(conj(model.si(s)(n).index()));
                block.addindex1(model.siP(s)(u).index());
                }
                else { block.addindex1(model.si(s)(n).index()); }

                qt[q].push_back(block);

                if(s==show_s)
                {
                block.print("Kept block",ShowData);
                cerr << "D = " << D << "\n";
                }
            }
        }}

        qC.clear();

        foreach(const qt_vt& x, qt)
        {
            const vector<ITensor>& blks = x.second;
            if(blks.size() != 0)
            {
                q = x.first; 
                if(s == N) 
                { foreach(const ITensor& t, blks) nblock.push_back(t); }
                else
                {
                    Matrix M; int mm = collapseCols(qD[q],M);
                    if(s==show_s)
                    {
                        cerr << format("Adding block, mm = %d\n")%mm;
                        q.print("q");
                        cerr << "qD[q] = " << qD[q] << "\n";
                        cerr << "M = \n" << M << "\n";
                        int count = 0;
                        foreach(const ITensor& t, blks) 
                        t.print((format("t%02d")%(++count)).str(),ShowData);
                    }
                    //string qname = (format("ql%d(%+d:%d:%s)")%s%q.sz()%q.Nf()%(q.Nfp() == 0 ? "+" : "-")).str();
                    string qname = (format("ql%d(%+d:%d)")%s%q.sz()%q.Nf()).str();
                    Index qbond(qname,mm);
                    ITensor compressor(bond,qbond,M);
                    foreach(const ITensor& t, blks) nblock.push_back(t * compressor);
                    iq.push_back(inqn(qbond,q));
                    qC[q] = compressor;
                }
            }
        }

        if(s != N) 
        { 
            if(iq.empty()) 
            {
                cerr << "At site " << s << "\n";
                Error("convertToIQ: no compatible QNs to put into Link.");
            }
            linkind[s] = IQIndex(nameint("qL",s),iq); iq.clear(); 
        }
        if(s == 1)
            qA[s] = (is_mpo ? IQTensor(conj(model.si(s)),model.siP(s),linkind[s]) : IQTensor(model.si(s),linkind[s]));
        else if(s == N)
            qA[s] = (is_mpo ? IQTensor(conj(linkind[s-1]),conj(model.si(s)),model.siP(s)) 
                                    : IQTensor(conj(linkind[s-1]),model.si(s)));
        else
            qA[s] = (is_mpo ? IQTensor(conj(linkind[s-1]),conj(model.si(s)),model.siP(s),linkind[s]) 
                                    : IQTensor(conj(linkind[s-1]),model.si(s),linkind[s]));

        foreach(const ITensor& nb, nblock) { qA[s] += nb; } nblock.clear();

        if(s==show_s)
        {
        qA[s].print((format("qA[%d]")%s).str(),ShowData);
        Error("Stopping");
        }

    } //for loop over s

    IQIndex Center("Center",Index("center",1,Virtual),totalq,In);
    qA[1].addindex1(Center);
}

class InitState
{
    int N;
    vector<IQIndexVal> state;
    typedef IQIndexVal (*SetFuncPtr)(int);
public:
    int NN() const { return N; }

    InitState(int nsite) : N(nsite), state(N+1) { }
    InitState(int nsite,SetFuncPtr setter) : N(nsite), state(N+1) 
    { set_all(setter); }

    void set_all(SetFuncPtr setter)
    { for(int j = 1; j <= N; ++j) GET(state,j-1) = (*setter)(j); }

    IQIndexVal& operator()(int i) { return GET(state,i-1); }
    const IQIndexVal& operator()(int i) const { return GET(state,i-1); }
    operator vector<IQIndexVal>() const { return state; }
};

namespace Internal {

template <class Tensor>
class MPS
{
public:
    typedef Tensor TensorT;
    typedef typename Tensor::IndexT IndexT;
    typedef typename Tensor::IndexValT IndexValT;
    typedef BaseModel ModelT;
protected:
    int N;
    vector<Tensor> A;
    int left_orth_lim,right_orth_lim;
    const ModelT* model_;

    void new_tensors(vector<ITensor>& A_)
    {
        vector<Index> a(N+1);
        for(int i = 1; i <= N; ++i)
        { a[i] = Index(nameint("a",i)); }
        A_[1] = ITensor(si(1),a[1]);
        for(int i = 2; i < N; i++)
        { A_[i] = ITensor(conj(a[i-1]),si(i),a[i]); }
        A_[N] = ITensor(conj(a[N-1]),si(N));
    }

    void random_tensors(vector<ITensor>& A_)
    { new_tensors(A_); for(int i = 1; i <= N; ++i) A_[i].Randomize(); }

    void random_tensors(vector<IQTensor>& A_) { }

    void init_tensors(vector<ITensor>& A_, const InitState& initState)
    { new_tensors(A_); for(int i = 1; i <= N; ++i) A_[i](initState(i))=1; }

    void init_tensors(vector<IQTensor>& A_, const InitState& initState)
    {
        vector<QN> qa(N+1); //qn[i] = qn on i^th bond
        for(int i = 1; i <= N; ++i) { qa[0] -= initState(i).qn()*In; }
        IQIndex Center("Center",Index("center",1,Virtual),qa[0],In);

        //Taking OC to be at the leftmost site,
        //compute the QuantumNumbers of all the Links.
        for(int i = 1; i <= N; ++i)
        {
            //Taking the divergence to be zero,solve for qa[i]
            qa[i] = Out*(-qa[i-1]*In - initState(i).qn());
        }

        vector<IQIndex> a(N+1);
        for(int i = 1; i <= N; ++i)
        { a[i] = IQIndex(nameint("L",i),Index(nameint("l",i)),qa[i]); }
        A_[1] = IQTensor(si(1),a[1]); A_[1].addindex1(Center); A_[1](initState(1))=1;
        for(int i = 2; i < N; ++i)
        { 
        A_[i] = IQTensor(conj(a[i-1]),si(i),a[i]); 
        A_[i](initState(i))=1;
        }
        A_[N] = IQTensor(conj(a[N-1]),si(N)); A_[N](initState(N))=1;
    }

    typedef pair<typename vector<Tensor>::const_iterator,typename vector<Tensor>::const_iterator> const_range_type;
public:
    int minm,maxm;
    Real cutoff;

    //Accessor Methods ------------------------------

    int NN() const { return N;}
    int right_lim() const { return right_orth_lim; }
    int left_lim() const { return left_orth_lim; }
    IQIndex si(int i) const { return model_->si(i); }
    IQIndex siP(int i) const { return model_->siP(i); }
    typedef typename vector<Tensor>::const_iterator AA_it;
    const pair<AA_it,AA_it> AA() const { return make_pair(A.begin()+1,A.end()); }
    const Tensor& AA(int i) const { return GET(A,i); }
    const ModelT& model() const { return *model_; }
    Tensor& AAnc(int i) //nc means 'non const'
    { 
        if(i <= left_orth_lim) left_orth_lim = i-1;
        if(i >= right_orth_lim) right_orth_lim = i+1;
        return GET(A,i); 
    }
    Tensor& setU(int i, Direction dir) //set unitary
    {
        if(dir == Fromleft)
        {
        if(i == left_orth_lim+1) left_orth_lim = i;
        else Error("MPS::setU: left_orth_lim not at i-1");
        }
        else if(dir == Fromright)
        {
        if(i == right_orth_lim-1) right_orth_lim = i;
        else Error("MPS::setU: right_orth_lim not at i-1");
        }
        return GET(A,i);
    }
    bool is_null() const { return (model_==0); }
    bool is_not_null() const { return (model_!=0); }

    Tensor bondTensor(int b) const { Tensor res = A.at(b) * A.at(b+1); return res; }

    //MPS: Constructors --------------------------------------------

    MPS() : N(0), model_(0), minm(1), maxm(MAX_M), cutoff(MAX_CUT) {}

    MPS(const ModelT& mod_,int maxmm = MAX_M, Real cut = MAX_CUT) 
		: N(mod_.NN()), A(mod_.NN()+1),left_orth_lim(0),right_orth_lim(mod_.NN()),
        model_(&mod_), minm(1), maxm(maxmm), cutoff(cut)
	{ random_tensors(A); }

    MPS(const ModelT& mod_,const InitState& initState,int maxmm = MAX_M, Real cut = MAX_CUT) 
		: N(mod_.NN()),A(mod_.NN()+1),left_orth_lim(0),right_orth_lim(2),
        model_(&mod_), minm(1), maxm(maxmm), cutoff(cut)
	{ init_tensors(A,initState); }

    MPS(const ModelT& model, istream& s) { read(model,s); }

    virtual ~MPS() { }

    void read(const ModelT& model, istream& s)
    {
        model_ = &model;
        N = model_->NN();
        A.resize(N);
        for(int j = 1; j <= N; ++j) A[j].read(s);
        s.read((char*) &left_orth_lim,sizeof(left_orth_lim));
        s.read((char*) &right_orth_lim,sizeof(right_orth_lim));
        s.read((char*) &minm,sizeof(minm));
        s.read((char*) &maxm,sizeof(maxm));
        s.read((char*) &cutoff,sizeof(cutoff));
    }

    void write(ostream& s) const
    {
        for(int j = 1; j <= N; ++j) A[j].write(s);
        s.write((char*) &left_orth_lim,sizeof(left_orth_lim));
        s.write((char*) &right_orth_lim,sizeof(right_orth_lim));
        s.write((char*) &minm,sizeof(minm));
        s.write((char*) &maxm,sizeof(maxm));
        s.write((char*) &cutoff,sizeof(cutoff));
    }


    //MPS: operators ------------------------------------------------------

    MPS& operator*=(Real a)
    {
        AAnc(left_orth_lim+1) *= a;
        return *this;
    }
    inline MPS operator*(Real r) const { MPS res(*this); res *= r; return res; }
    friend inline MPS operator*(Real r, MPS res) { res *= r; return res; }

    MPS& operator+=(const MPS& oth);
    inline MPS operator+(MPS res) const { res += *this; return res; }
    //friend inline MPS operator+(MPS A, const MPS& B) { A += B; return A; }
    inline MPS operator-(MPS res) const { res *= -1; res += *this; return res; }

    //MPS: index methods --------------------------------------------------

    void mapprime(int oldp, int newp, PrimeType pt = primeBoth)
	{ for(int i = 1; i <= N; ++i) A[i].mapprime(oldp,newp,pt); }

    void primelinks(int oldp, int newp)
	{ for(int i = 1; i <= N; ++i) A[i].mapprime(oldp,newp,primeLink); }

    void noprimelink()
	{ for(int i = 1; i <= N; ++i) A[i].noprime(primeLink); }

    IndexT LinkInd(int i) const { return index_in_common(A[i],A[i+1],Link); }
    IndexT RightLinkInd(int i) const { assert(i<NN()); return index_in_common(AA(i),AA(i+1),Link); }
    IndexT LeftLinkInd(int i)  const { assert(i>1); return index_in_common(AA(i),AA(i-1),Link); }

    //MPS: orthogonalization methods -------------------------------------

    virtual void doSVD(int i, const Tensor& AA, Direction dir, bool preserve_shape = false)
	{
        /*
        if(!preserve_shape) 
        {
            const int i1 = (dir == Fromleft ? i : i+1);
            const int i2 = (dir == Fromleft ? i-1 : i+2);
            const int boundary = (dir == Fromleft ? 1 : NN()-1);

            vector<IndexT> keep_indices;
            IndexT site = A.at(i1).findtype(Site); site.noprime();
            keep_indices.push_back(site);

            int j = A.at(i1).findindexn(site.primed());
            if(j > -1) //part of an MPO
            { keep_indices.push_back(site.primed()); }

            if(i != boundary) 
            {
            foreach(const IndexT& I, AA.indexn())
            if(A.at(i2).hasindexn(I)) keep_indices.push_back(I);
            foreach(const IndexT& I, AA.index1())
            if(A.at(i2).hasindex1(I)) keep_indices.push_back(I);
            }
        }
        //otherwise the MPS will retain the same index structure
        */

        tensorSVD(AA,GET(A,i),GET(A,i+1),cutoff,minm,maxm,dir,0);
        truncerror = svdtruncerr;

        if(dir == Fromleft)
        {
            if(left_orth_lim == i-1 || i == 1) left_orth_lim = i;
            if(right_orth_lim < i+2) right_orth_lim = i+2;
        }
        else
        {
            if(left_orth_lim > i-1) left_orth_lim = i-1;
            if(right_orth_lim == i+2 || i == N-1) right_orth_lim = i+1;
        }
	}

    //Move the orthogonality center to site i (left_orth_lim = i-1, right_orth_lim = i+1)
    void position(int i, bool preserve_shape = false)
	{
        //if(is_null()) { left_orth_lim = i-1; right_orth_lim = i+1; return; }
        if(is_null()) Error("position: MPS is null");
        while(left_orth_lim < i-1)
        {
            if(left_orth_lim < 0) left_orth_lim = 0;
            Tensor WF = AA(left_orth_lim+1) * AA(left_orth_lim+2);
            doSVD(left_orth_lim+1,WF,Fromleft,preserve_shape);
        }
        while(right_orth_lim > i+1)
        {
            if(right_orth_lim > N+1) right_orth_lim = N+1;
            Tensor WF = AA(right_orth_lim-2) * AA(right_orth_lim-1);
            doSVD(right_orth_lim-2,WF,Fromright,preserve_shape);
        }
	}

    bool is_ortho() const { return (left_orth_lim + 1 == right_orth_lim - 1); }

    int ortho_center() const 
    { 
        if(!is_ortho()) Error("MPS: orthogonality center not well defined.");
        return (left_orth_lim + 1);
    }

    //Checks if A[i] is left (left == true) or right (left == false) orthogonalized
    bool checkOrtho(int i, bool left) const
    {
        IndexT link = (left ? RightLinkInd(i) : LeftLinkInd(i));
        Tensor A = AA(i);
        Tensor Ac = conj(A); Ac.primeind(link,4);

        Tensor rho = A * Ac;

        Tensor Delta(link,link.primed(4),1);

        Tensor Diff = rho - Delta;

        Vector diff(Diff.Length());
        Diff.AssignToVec(diff);

        Real threshold = 1E-13;
        if(Norm(diff) < threshold) return true;

        //Print any helpful debugging info here:
        cerr << "checkOrtho: on line " << __LINE__ << " of mps.h," << endl;
        cerr << "checkOrtho: Tensor at position " << i << " failed to be " << (left ? "left" : "right") << " ortho." << endl;
        cerr << "checkOrtho: Norm(diff) = " << boost::format("%E") % Norm(diff) << endl;
        cerr << "checkOrtho: Error threshold set to " << boost::format("%E") % threshold << endl;
        //-----------------------------

        return false;
    }
    bool checkRightOrtho(int i) const { return checkOrtho(i,false); }
    bool checkLeftOrtho(int i) const { return checkOrtho(i,true); }
    
    bool checkOrtho() const
    {
        for(int i = 1; i <= left_orth_lim; ++i)
        if(!checkLeftOrtho(i))
        {
            cerr << "checkOrtho: A[i] not left orthogonal at site i=" << i << endl;
            return false;
        }

        for(int i = NN(); i >= right_orth_lim; --i)
        if(!checkRightOrtho(i))
        {
            cerr << "checkOrtho: A[i] not right orthogonal at site i=" << i << endl;
            return false;
        }
        return true;
    }

    void getCenter(int j, Direction dir, Tensor& lambda, bool do_signfix = false)
    {
        getCenterMatrix(AAnc(j),(dir == Fromleft ? RightLinkInd(j) : LeftLinkInd(j)),cutoff,minm,maxm,lambda,nameint("c",j));

        if(dir == Fromleft) if(left_orth_lim == j-1 || j == 1) left_orth_lim = j;
        else if(right_orth_lim == j+1 || j == N) right_orth_lim = j;

        if(do_signfix) Error("MPS::getCenter: do_signfix not implemented.");
    }

    template<class TensorSet>
    Real bondDavidson(int b, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH, int niter, int debuglevel, Direction dir, Real errgoal=1E-4)
    {
        if(b-1 > left_orth_lim)
        {
            cerr << format("MPS::bondDavidson: b=%d, Lb=%d\n")%b%left_orth_lim;
            Error("MPS::bondDavidson: b > left_orth_lim");
        }
        if(b+2 < right_orth_lim)
        {
            cerr << format("MPS::bondDavidson: b+1=%d, Rb=%d\n")%(b+1)%right_orth_lim;
            Error("MPS::bondDavidson: b+1 < right_orth_lim");
        }
        Tensor phi = GET(A,b); phi *= GET(A,b+1);
        Real En = doDavidson(phi,mpoh,LH,RH,niter,debuglevel,errgoal);
        doSVD(b,phi,dir);
        return En;
    }

    template<class TensorSet, class OpTensorSet>
    void projectOp(int j, Direction dir, const TensorSet& P, const OpTensorSet& Op, TensorSet& res) const
    {
        if(res.size() != Op.size()) res.resize(Op.size());
        const TensorSet& nP = (P.size() == Op.size() ? P : TensorSet(Op.size()));
        for(unsigned int n = 0; n < Op.size(); ++n) projectOp(j,dir,GET(nP,n),Op[n],GET(res,n));
    }

    template<class OpTensor>
    void projectOp(int j, Direction dir, const Tensor& P, const OpTensor& Op, Tensor& res) const
    {
        if(dir==Fromleft && j > left_orth_lim) 
        { cerr << format("projectOp: from left j > left_orth_lim (j=%d,left_orth_lim=%d)\n")%j%left_orth_lim, Error(""); }
        if(dir==Fromright && j < right_orth_lim) 
        { cerr << format("projectOp: from left j < right_orth_lim (j=%d,right_orth_lim=%d)\n")%j%right_orth_lim, Error(""); }

        res = (P.is_null() ? AA(j) : P * AA(j));
        res *= Op; res *= conj(primed(AA(j)));
    }


    void applygate(const Tensor& gate)
	{
        Tensor AA = A[left_orth_lim+1] * A[left_orth_lim+2] * gate;
        AA.noprime();
        doSVD(left_orth_lim+1,AA,Fromleft);
	}

    Real normalize()
    {
        Real norm2,im;
        psiphi(*this,*this,norm2,im);
        Real norm = sqrt(norm2);
        *this *= 1.0/norm;
        return norm;
    }

    bool is_complex() const
    { return A[left_orth_lim+1].is_complex(); }

    friend inline ostream& operator<<(ostream& s, const MPS& M)
    {
        s << "\n";
        for(int i = 1; i <= M.NN(); ++i) s << M.AA(i) << "\n";
        return s;
    }

    void print(string name = "",Printdat pdat = HideData) const 
    { printdat = (pdat==ShowData); cerr << "\n" << name << " =\n" << *this << "\n"; printdat = false; }

    //-----------------------------------------------------------------
    //IQMPS specific methods

    template <class IQMPSType> 
    void convertToIQ(IQMPSType& iqpsi, QN totalq = QN(), Real cut = 1E-12) const
    {
        assert(model_ != 0);
        const ModelT& sst = *model_;

        iqpsi = IQMPSType(sst,maxm,cutoff);

        if(!A[1].hasindex(si(1))) Error("convertToIQ: incorrect primelevel for conversion");
        bool is_mpo = A[1].hasindex(si(1).primed());
        const int Dim = si(1).m();
        const int PDim = (is_mpo ? Dim : 1);

        vector<IQIndex> linkind(N);

        typedef map<QN,Vector>::value_type qD_vt;
        map<QN,Vector> qD; //Diags of compressor matrices by QN
        typedef map<QN,vector<ITensor> >::value_type qt_vt;
        map<QN,vector<ITensor> > qt; //ITensor blocks by QN
        typedef map<QN,ITensor>::value_type qC_vt;
        map<QN,ITensor> qC; //Compressor ITensors by QN
        ITensor block;
        vector<ITensor> nblock;
        vector<inqn> iq;

        QN q;

        qC[totalq] = ITensor(); //Represents Virtual index
        //First value of prev_q below set to totalq

        const int show_s = 0;

        Index bond, prev_bond;
        for(int s = 1; s <= N; ++s)
        {
            if(debug3) cerr << "s = " << s << "\n";

            qD.clear(); qt.clear();
            if(s > 1) prev_bond = LinkInd(s-1); 
            if(s < N) bond = LinkInd(s);

            if(s == show_s) 
            {
                PrintDat(A[s]);
            }

            foreach(const qC_vt& x, qC) {
            const QN& prev_q = x.first; const ITensor& comp = x.second; 
            for(int n = 1; n <= Dim;  ++n)
            for(int u = 1; u <= PDim; ++u)
            {
                q = (is_mpo ? prev_q+si(s).qn(n)-si(s).qn(u) : prev_q-si(s).qn(n));

                //For the last site, only keep blocks 
                //compatible with specified totalq i.e. q=0 here
                if(s == N && q != QN()) continue;

                //Set Site indices of A[s] and its previous Link Index
                block = A[s];
                if(s != 1) block *= conj(comp);
                block *= si(s)(n);
                if(is_mpo) block *= siP(s)(u);

                //Initialize D Vector (D records which values of
                //the right Link Index to keep for the current QN q)
                int count = qD.count(q);
                Vector& D = qD[q];
                if(count == 0) { D.ReDimension(bond.m()); D = 0; }

                if(s == show_s)
                {
                    cerr << format("For n = %d\n")%n;
                    cerr << format("Got a block with norm %.10f\n")%block.norm();
                    cerr << format("bond.m() = %d\n")%bond.m();
                    PrintDat(block);
                    if(s != 1) PrintDat(comp);
                }

                bool keep_block = false;
                if(s == N) keep_block = true;
                else
                {
                    if(bond.m() == 1 && block.norm() != 0) { D = 1; keep_block = true; }
                    else
                    {
                        ITensor summed_block;
                        if(s==1) summed_block = block;
                        else
                        {
                            //Here we sum over the previous link index
                            //which is already ok, analyze the one to the right
                            assert(comp.r()==2);
                            Index new_ind = (comp.index(1)==prev_bond ? comp.index(2) : comp.index(1));
                            summed_block = ITensor(new_ind,1) * block;
                        }
                        //cerr << format("s = %d, bond=")%s << bond << "\n";
                        //summed_block.print("summed_block");

                        Real rel_cut = -1;
                        for(int j = 1; j <= bond.m(); ++j)
                        { rel_cut = max(fabs(summed_block.val1(j)*exp(summed_block.logfac())),rel_cut); }
                        assert(rel_cut >= 0);
                        //Real rel_cut = summed_block.norm()/summed_block.vec_size();
                        rel_cut *= cut;
                        //cerr << "rel_cut == " << rel_cut << "\n";

                        if(rel_cut > 0)
                        for(int j = 1; j <= bond.m(); ++j)
                        if(fabs(summed_block.val1(j)*exp(summed_block.logfac())) > rel_cut) 
                        { D(j) = 1; keep_block = true; }
                    }
                } //else (s != N)

                //if(!keep_block && q == totalq) { D(1) = 1; keep_block = true; }

                if(keep_block)
                {
                    qD[q] = D;

                    if(is_mpo) 
                    {
                    block.addindex1(conj(si(s)(n).index()));
                    block.addindex1(siP(s)(u).index());
                    }
                    else { block.addindex1(si(s)(n).index()); }

                    qt[q].push_back(block);

                    if(s==show_s)
                    {
                    block.print("Kept block",ShowData);
                    cerr << "D = " << D << "\n";
                    }
                }
            }}

            qC.clear();

            foreach(const qt_vt& x, qt)
            {
                const vector<ITensor>& blks = x.second;
                if(blks.size() != 0)
                {
                    q = x.first; 
                    if(s == N) 
                    { foreach(const ITensor& t, blks) nblock.push_back(t); }
                    else
                    {
                        Matrix M; int mm = collapseCols(qD[q],M);
                        if(s==show_s)
                        {
                            cerr << format("Adding block, mm = %d\n")%mm;
                            q.print("q");
                            cerr << "qD[q] = " << qD[q] << "\n";
                            cerr << "M = \n" << M << "\n";
                            int count = 0;
                            foreach(const ITensor& t, blks) 
                            t.print((format("t%02d")%(++count)).str(),ShowData);
                        }
                        //string qname = (format("ql%d(%+d:%d:%s)")%s%q.sz()%q.Nf()%(q.Nfp() == 0 ? "+" : "-")).str();
                        string qname = (format("ql%d(%+d:%d)")%s%q.sz()%q.Nf()).str();
                        Index qbond(qname,mm);
                        ITensor compressor(bond,qbond,M);
                        foreach(const ITensor& t, blks) nblock.push_back(t * compressor);
                        iq.push_back(inqn(qbond,q));
                        qC[q] = compressor;
                    }
                }
            }

            if(s != N) 
            { 
                if(iq.empty()) 
                {
                    cerr << "At site " << s << "\n";
                    Error("convertToIQ: no compatible QNs to put into Link.");
                }
                linkind[s] = IQIndex(nameint("qL",s),iq); iq.clear(); 
            }
            if(s == 1)
            {
                iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(si(s)),siP(s),linkind[s]) : IQTensor(si(s),linkind[s]));
                IQIndex Center("Center",Index("center",1,Virtual),totalq,In);
                iqpsi.AAnc(1).addindex1(Center);
            }
            else if(s == N)
            {
                iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(linkind[s-1]),conj(si(s)),siP(s)) 
                                        : IQTensor(conj(linkind[s-1]),si(s)));
            }
            else
            {
                iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(linkind[s-1]),conj(si(s)),siP(s),linkind[s]) 
                                        : IQTensor(conj(linkind[s-1]),si(s),linkind[s]));
            }

            foreach(const ITensor& nb, nblock) { iqpsi.AAnc(s) += nb; } nblock.clear();

            if(0) //try to get this working ideally
            if(!is_mpo && s > 1) 
            {
                IQTensor AA = iqpsi.bondTensor(s-1);
                iqpsi.doSVD(s-1,AA,Fromleft);
            }

            if(s==show_s)
            {
            iqpsi.AA(s).print((format("qA[%d]")%s).str(),ShowData);
            Error("Stopping");
            }

        } //for loop over s

        assert(check_QNs(iqpsi));

    } //void convertToIQ(IQMPSType& iqpsi) const

}; //class MPS<Tensor>

namespace {
void plussers(const Index& l1, const Index& l2, Index& sumind, ITensor& first, ITensor& second)
{
    sumind = Index(sumind.rawname(),l1.m()+l2.m(),sumind.type());
    first = ITensor(l1,sumind,1);
    second = ITensor(l2,sumind);
    for(int i = 1; i <= l2.m(); ++i) second(l2(i),sumind(l1.m()+i)) = 1;
}

void plussers(const IQIndex& l1, const IQIndex& l2, IQIndex& sumind, IQTensor& first, IQTensor& second)
{
    map<Index,Index> l1map, l2map;
    vector<inqn> iq;
    foreach(const inqn& x, l1.iq())
	{
        Index ii = x.index;
        Index jj(ii.name(),ii.m(),ii.type());
        l1map[ii] = jj;
        iq.push_back(inqn(jj,x.qn));
	}
    foreach(const inqn& x, l2.iq())
	{
        Index ii = x.index;
        Index jj(ii.name(),ii.m(),ii.type());
        l2map[ii] = jj;
        iq.push_back(inqn(jj,x.qn));
	}
    sumind = IQIndex(sumind,iq);
    first = IQTensor(conj(l1),sumind);
    foreach(const inqn& x, l1.iq())
	{
        Index il1 = x.index;
        Index s1 = l1map[il1];
        ITensor t(il1,s1,1.0);
        first += t;
	}
    second = IQTensor(conj(l2),sumind);
    foreach(const inqn& x, l2.iq())
	{
        Index il2 = x.index;
        Index s2 = l2map[il2];
        ITensor t(il2,s2,1.0);
        second += t;
	}
}
} //anonymous namespace

template <class Tensor>
MPS<Tensor>& MPS<Tensor>::operator+=(const MPS<Tensor>& other)
{
    primelinks(0,4);

    vector<Tensor> first(N), second(N);
    for(int i = 1; i < N; ++i)
    {
        IndexT l1 = this->RightLinkInd(i);
        IndexT l2 = other.RightLinkInd(i);
        IndexT r(l1.rawname());
        plussers(l1,l2,r,first[i],second[i]);
    }

    AAnc(1) = AA(1) * first[1] + other.AA(1) * second[1];
    for(int i = 2; i < N; ++i)
    {
        AAnc(i) = conj(first[i-1]) * AA(i) * first[i] + conj(second[i-1]) * other.AA(i) * second[i];
    }
    AAnc(N) = conj(first[N-1]) * AA(N) + conj(second[N-1]) * other.AA(N);

    noprimelink();

    position(N);
    position(1);

    return *this;
}

} //namespace Internal
typedef Internal::MPS<ITensor> MPS;
typedef Internal::MPS<IQTensor> IQMPS;

inline void diag_denmat(const ITensor& rho, Real cutoff, int minm, int maxm, ITensor& nU, Vector& D, Real logrefnorm = DefaultLogRef)
{
    assert(rho.r() == 2);
    Index active = rho.index(1); active.noprime();

    if(logrefnorm != DefaultLogRef) rho.normlogto(logrefnorm);

    //Do the diagonalization
    Index ri = rho.index(1); ri.noprime();
    Matrix R,U; rho.toMatrix11(ri,ri.primed(),R);
    R *= -1.0; EigenValues(R,D,U); D *= -1.0;

    //Truncate
    Real err = 0.0;
    int mp = D.Length();
    while(mp > maxm || (err+D(mp) < cutoff*D(1) && mp > minm)) err += D(mp--);
    svdtruncerr = (D(1) == 0 ? 0.0 : err/D(1));
    D.ReduceDimension(mp); lastd = D;
    Index newmid(active.rawname(),mp,active.type());
    nU = ITensor(active,newmid,U.Columns(1,mp));
}

inline void diag_denmat(const IQTensor& rho, Real cutoff, int minm, int maxm, IQTensor& nU, Vector& eigs_kept, Real logrefnorm = DefaultLogRef)
{
    IQIndex active = rho.finddir(Out);
    assert(active.primelevel == 0);

    vector<Matrix> mmatrix(rho.iten_size());
    vector<Vector> mvector(rho.iten_size());
    vector<Real> alleig;

    Real maxlogfac = 0;
    if(logrefnorm == DefaultLogRef) 
    {
        foreach(const ITensor& t, rho.itensors())
        { maxlogfac = max(maxlogfac,t.logfac()); }
    }
    else { maxlogfac = logrefnorm; }

    assert(maxlogfac < 100);

    //1. Diagonalize each ITensor within rho
    int itenind = 0;
    foreach(const ITensor& t, rho.itensors())
	{
        assert(t.index(1).noprime_equals(t.index(2)));
        //if(!t.index(1).noprime_equals(t.index(2)))
        //{ Print(rho); Print(t); Error("Non-symmetric ITensor in density matrix"); }

        t.normlogto(maxlogfac);

        Matrix &U = GET(mmatrix,itenind);
        Vector &d = GET(mvector,itenind);

        //Diag ITensors within rho
        int n = t.index(1).m();
        Matrix M(n,n);
        //Real lfac; t.toMatrix11(t.index(2),t.index(1),M,lfac);
        t.toMatrix11(t.index(2),t.index(1),M);

        M *= -1;
        EigenValues(M,d,U);
        d *= -1;

        for(int j = 1; j <= n; ++j) 
        { alleig.push_back(d(j)); }
        ++itenind;
	}

    //2. Truncate eigenvalues

    //Sort all eigenvalues from smallest to largest
    //irrespective of quantum numbers
    sort(alleig.begin(),alleig.end());

    //Truncate
    Real e1 = max(alleig.back(),1.0e-60), docut = 0;
    svdtruncerr = 0;
    int m = 0, mkeep = (int)alleig.size();
    if(mkeep > minm)
    for(; m < (int)alleig.size(); m++, mkeep--)
    if(((svdtruncerr += GET(alleig,m)/e1) > cutoff && mkeep <= maxm) || mkeep <= minm)
    { 
        docut = (m > 0 ?  (GET(alleig,m-1) + GET(alleig,m))/2 : 0);
        svdtruncerr -= GET(alleig,m)/e1;
        break; 
    }
    //cerr << "\nDiscarded " << m << " states in diag_denmat\n";
    m = (int)alleig.size()-m;

    assert(m <= maxm); 
    assert(m < 20000);

    //3. Construct orthogonalized IQTensor nU
    vector<ITensor> terms; terms.reserve(rho.iten_size());
    vector<inqn> iq; iq.reserve(rho.iten_size());
    itenind = 0;
    foreach(const ITensor& t, rho.itensors())
	{
        const Vector& thisD = GET(mvector,itenind);
        int this_m = 1;
        for(; this_m <= thisD.Length(); ++this_m)
        if(thisD(this_m) < docut) { break; }

        if(mkeep == 0 && thisD.Length() >= 1) // zero mps, just keep one arb state
        { this_m = 2; mkeep = 1; docut = 1; }
        --this_m; //since for loop overshoots by 1

        if(this_m == 0) { ++itenind; continue; }

        Index nm("qlink",this_m);
        Index act = t.index(1).deprimed();
        iq.push_back(inqn(nm,active.qn(act)));

        Matrix UU = GET(mmatrix,itenind).Columns(1,this_m);

        ITensor term(act,nm); term.fromMatrix11(act,nm,UU); 
        //term.setlogfac(maxlogfac);
        terms.push_back(term);

        ++itenind;
	}
    IQIndex newmid("qlink",iq,In);
    nU = IQTensor(active,newmid);
    foreach(const ITensor& t, terms) nU += t;

    eigs_kept.ReDimension(m);
    for(int i = 1; i <= m; ++i) eigs_kept(i) = GET(alleig,alleig.size()-i);
    lastd = eigs_kept;
} //void diag_denmat

template<class Tensor>
Vector tensorSVD(const Tensor& AA, Tensor& A, Tensor& B, Real cutoff, int minm, int maxm, Direction dir, Real logrefnorm = DefaultLogRef)
{
    typedef typename Tensor::IndexT IndexT;
    typedef typename Tensor::CombinerT CombinerT;

    if(AA.vec_size() == 0) Error("tensorSVD(Tensor): input tensor had zero size.");

    IndexT mid = index_in_common(A,B,Link);
    if(mid.is_null()) mid = IndexT("mid");

    Tensor& to_orth = (dir==Fromleft ? A : B);
    Tensor& newoc   = (dir==Fromleft ? B : A);

    CombinerT comb;

    int unique_link = 0; //number of Links unique to to_orth
    for(int j = 1; j <= to_orth.r(); ++j) { 
    const IndexT& I = to_orth.index(j);
    if(!(newoc.hasindex(I) || I == Tensor::ReImIndex || I.type() == Virtual))
    {
        if(I.type() == Link) ++unique_link;
        comb.addleft(I);
    }}

    //Check if we're at the edge
    if(unique_link == 0)
    {
        comb.init(mid.rawname());
        assert(comb.check_init());
        comb.product(AA,newoc);
        to_orth = comb; to_orth.conj();
        Vector eigs_kept(comb.right().m()); eigs_kept = 1.0/comb.right().m();
        return eigs_kept; 
    }

    //Apply combiner
    comb.doCondense(true);
    comb.init(mid.rawname());
    Tensor AAc; comb.product(AA,AAc);

    const IndexT& active = comb.right();

    Tensor rho;
    if(AAc.is_complex())
    {
        Tensor re,im;
        AAc.SplitReIm(re,im);
        rho = re; rho.conj(); rho.primeind(active);
        rho *= re;
        im *= conj(primeind(im,active));
        rho += im;
    }
    else { Tensor AAcc = conj(AAc); AAcc.primeind(active); rho = AAc*AAcc; }
    assert(rho.r() == 2);

    Tensor U; Vector eigs_kept;
    diag_denmat(rho,cutoff,minm,maxm,U,eigs_kept,logrefnorm);

    comb.conj();
    comb.product(U,to_orth);
    newoc = conj(U) * AAc;

    return eigs_kept;
}

/*
Vector do_denmat_Real(const vector<IQTensor>& nA, const IQTensor& A, const IQTensor& B, IQTensor& U,
	Real cutoff, int minm,int maxm, Direction dir, bool donormalize, bool do_relative_cutoff)
{
    // Make a density matrix that is summed over the nA
    IQIndex mid = index_in_common(A,B,Link);
    if(mid.is_null()) mid = IQIndex("mid");

    int num_states = nA.size();
    if(num_states == 0) Error("zero size in do_denmat_Real(vector<IQTensor>)");

    const IQTensor& to_orth = (dir==Fromleft ? A : B);
    const IQTensor& newoc   = (dir==Fromleft ? B : A);

    int unique_link = 0;
    //Create combiner
    IQCombiner comb;
    foreach(const IQIndex& i, to_orth.iqinds())
    if(!(newoc.hasindex(i) || i == IQIndReIm || i.type() == Virtual))
    { 
        if(i.type() == Link) ++unique_link;
        comb.addleft(i); 
    }

    //Check if we're at the edge
    //bool edge_case = (to_orth.num_index(Link) <= 1 ? true : false);
    bool edge_case = (unique_link == 0 ? true : false);

    if(edge_case)
    {
        //Handle the right-edge/Fromright and left-edge/Fromleft
        //cases by simply turning the appropriate IQCombiner into
        //an IQTensor and using it as the new edge tensor

        //IQIndex newmid(mid.rawname(),Link,Out);

        comb.init(mid.rawname());

        U = comb; U.conj(); 

        Vector eigs_kept(comb.right().m()); eigs_kept = 1.0/comb.right().m();
        return eigs_kept; 
    }

    //Combine
    //IQIndex c(mid.rawname());
    //comb.init(c);
    comb.doCondense(false);
    comb.init(mid.rawname());
    IQIndex c = comb.right();
    vector<IQTensor> nnA; nnA.reserve(nA.size());
    foreach(const IQTensor& iqt, nA) nnA.push_back(iqt * comb);

    //Condense
    vector<IQTensor> ncA; ncA.reserve(nnA.size());
    IQIndex active("Condensed");
    Condenser cond(c,active);
    foreach(const IQTensor& iqt, nnA) ncA.push_back(iqt * cond);

    IQIndex activep(primeBoth,active,4);
    IQTensor rho;
    if(ncA.front().hasindex(IQIndReIm))
    {
        //Error("not doing complex 873249827");
        foreach(const IQTensor& iqt, ncA)
        {
            IQTensor iqtre,iqtim;
            iqt.SplitReIm(iqtre,iqtim);
            IQTensor iqtreconj = conj(iqtre), iqtimconj = conj(iqtim);
            iqtreconj.ind_inc_prime(active,4);
            //cout << "k = " << k << endl;
            //cout << "ncA[k] is " << ncA[k] << endl;
            //cout << "ncAconj is " << ncAconj << endl;
            IQTensor r = iqtre * iqtreconj;
            iqtimconj.ind_inc_prime(active,4);
            r += iqtim * iqtimconj;
            //cout << "r is " << r << endl;
            if(num_states == 1) rho = r;
            else                rho += r;
            //cout << "now rho is " << rho << endl;
        }
    }
    else
    {
        foreach(const IQTensor& iqt, ncA)
        {
            IQTensor iqtconj = conj(iqt);
            iqtconj.ind_inc_prime(active,4);
            IQTensor r = iqt * iqtconj;
            if(num_states == 1) rho = r;
            else                rho += r;
        }
    }
    if(rho.iten_size() == 0)
	{
        Error("rho iten_size is 0!!!");
	}

    rho *= 1.0/num_states;

    IQTensor nU; Vector eigs_kept;
    diag_denmat(rho,cutoff,minm,maxm,nU,eigs_kept);

    U = nU * cond * conj(comb);

    return eigs_kept;

} //void do_denmat_Real(const vector<IQTensor>& nA, ... )
*/

/* getCenterMatrix:
 * 
 *                    s                   s
 *                    |                   |
 * Decomposes A = -<--A-->- bond into -<--U-<-- -<--Lambda-->- bond
 *
 * A is replaced with the unitary U and Lambda is diagonal.
 * If A is the OC of an MPS, Lambda will contain the Schmidt weights. 
 *
 */
inline void getCenterMatrix(ITensor& A, const Index& bond, Real cutoff,int minm, int maxm, ITensor& Lambda, string newbondname = "")
{
    //Create combiner
    Combiner comb;
    foreach(const Index& i, A.indexn())
    if(!(i == bond || i == IndReIm || i.type() == Virtual))
    { 
        comb.addleft(i); 
    }
    foreach(const Index& i, A.index1())
    if(!(i == bond || i == IndReIm || i.type() == Virtual))
    { 
        comb.addleft(i); 
    }
    comb.init("combined");
    Index active = comb.right();

    //Apply combiner....
    //comb.init(active);
    ITensor Ac = comb * A;

    ITensor rho;
    if(Ac.is_complex())
    {
        ITensor re,im;
        Ac.SplitReIm(re,im);
        ITensor rec = conj(re), imc = conj(im);
        rec.primeind(active);
        rho = re * rec;
        imc.primeind(active);
        rho += im * imc;
    }
    else rho = Ac * primeind(Ac,active);
    assert(rho.r() == 2);

    //Diagonalize & truncate the density matrix
    ITensor Uc; Vector D; diag_denmat(rho,cutoff,minm,maxm,Uc,D);

    Lambda = conj(Uc) * Ac;
    A = Uc * comb; //should be conj(comb) with arrows

    assert(A.checkDim());
    assert(Lambda.checkDim());
}

inline bool check_QNs(const MPS& psi) { return true; }

inline bool check_QNs(const IQMPS& psi)
{
    const int N = psi.NN();
    //Check Link arrows

    //Get the orthogonality center (based on location of Virtual IQIndex)
    int center = -1;
    for(int i = 1; i <= N; ++i) 
    {
        if(psi.AA(i).hastype(Virtual)) { center = i; break; } 
        //else if(i == N) { cerr << "check_QNs: couldn't find a Virtual IQIndex.\n"; return false; }
    }

    if(center > 0)
    {
        //Check arrows from left edge
        for(int i = 1; i < center; ++i)
        {
            if(psi.RightLinkInd(i).dir() != In) 
            {
                cerr << format("check_QNs: At site %d to the left of the OC, Right side Link not pointing In\n")%i;
                return false;
            }
            if(i > 1)
            {
                if(psi.LeftLinkInd(i).dir() != Out) 
                {
                    cerr << format("check_QNs: At site %d to the left of the OC, Left side Link not pointing Out\n")%i;
                    return false;
                }
            }
        }

        //Check arrows from right edge
        for(int i = N; i > center; --i)
        {
            if(i < N)
            if(psi.RightLinkInd(i).dir() != Out) 
            {
                cerr << format("check_QNs: At site %d to the right of the OC, Right side Link not pointing Out\n")%i;
                return false;
            }
            if(psi.LeftLinkInd(i).dir() != In) 
            {
                cerr << format("check_QNs: At site %d to the right of the OC, Left side Link not pointing In\n")%i;
                return false;
            }
        }
    }
    //Done checking arrows

    //Check IQTensors
    for(int i = 1; i <= N; ++i)
    {
        if(!check_QNs(psi.AA(i)))
        {
            cerr << "check_QNs: IQTensor AA(" << i << ") had non-zero divergence.\n";
            return false;
        }
    }
    return true;
}

inline QN total_QN(const IQMPS& psi)
{
    assert(psi.AA(psi.ortho_center()).has_virtual());
    return psi.AA(psi.ortho_center()).virtualQN();
}

template <class MPSType>
void psiphi(const MPSType& psi, const MPSType& phi, Real& re, Real& im)  // <psi | phi>
{
    const int N = psi.NN();
    if(N != phi.NN()) Error("psiphi: mismatched N");

    typename MPSType::TensorT L = phi.AA(1) * conj(primeind(psi.AA(1),psi.LinkInd(1))); 
    for(int i = 2; i < psi.NN(); ++i) L = L * phi.AA(i) * conj(primelink(psi.AA(i)));
    L = L * phi.AA(N);

    Dot(primeind(psi.AA(N),psi.LinkInd(N-1)),L,re,im);
}
template <class MPSType>
Real psiphi(const MPSType& psi, const MPSType& phi) //Re[<psi|phi>]
{
    Real re, im;
    psiphi(psi,phi,re,im);
    if(im != 0) cerr << "Real psiphi: WARNING, dropping non-zero imaginary part of expectation value.\n";
    return re;
}

//Computes an MPS which has the same overlap with psi_basis as psi_to_fit,
//but which differs from psi_basis only on the first site, and has same index
//structure as psi_basis. Result is stored to psi_to_fit on return.
inline void fitWF(const IQMPS& psi_basis, IQMPS& psi_to_fit)
{
    if(!psi_basis.is_ortho()) Error("fitWF: psi_basis must be orthogonolized.");
    if(psi_basis.ortho_center() != 1) Error("fitWF: psi_basis must be orthogonolized to site 1.");
    psi_to_fit.position(1);

    const IQMPS& psib = psi_basis;
    IQMPS& psif = psi_to_fit;
    IQMPS res = psib;

    int N = psib.NN();
    if(psif.NN() != N) Error("fitWF: sizes of wavefunctions must match.");

    IQTensor A = psif.AA(N) * conj(psib.AA(N));
    for(int n = N-1; n > 1; --n)
    {
        A = conj(psib.AA(n)) * A;
        A = psif.AA(n) * A;
    }
    A = psif.AA(1) * A;

    res.AAnc(1) = A;

    assert(check_QNs(res));

    psi_to_fit = res;
}

//Template method for efficiently summing a set of MPS's or MPO's (or any class supporting operator+=)
template <typename MPSType>
void sum(const vector<MPSType>& terms, MPSType& res, Real cut = MAX_CUT, int maxm = MAX_M)
{
    int Nt = terms.size();
    if(Nt == 1) 
    {
        res = terms[0];
        res.cutoff = cut; res.maxm = maxm;
        return;
    }
    else if(Nt == 2)
	{ 
        res = terms[0];
        res.cutoff = cut; res.maxm = maxm;
        //cerr << format("Before +=, cutoff = %.1E, maxm = %d\n")%(res.cutoff)%(res.maxm);
        res += terms[1];
        return;
    }
    else if(Nt > 2)
	{
        //Add all MPS's in pairs
        vector<MPSType> terms2(2), nterms; nterms.reserve(Nt/2);
        for(int n = 0; n < Nt-1; n += 2)
        {
            terms2[0] = terms[n]; terms2[1] = terms[n+1];
            sum(terms2,res,cut,maxm);
            nterms.push_back(res);
        }
        if(Nt%2 == 1) nterms.push_back(terms.back());
        //Recursively call sum again
        sum(nterms,res,cut,maxm);
        return;
	}
    return;
} // void sum(const vector<MPSType>& terms, Real cut, int maxm, MPSType& res)


#ifdef THIS_IS_MAIN
Real truncerror = 0.0;
Real svdtruncerr = 0.0;
#endif

#endif
