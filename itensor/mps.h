#ifndef __MPS_H
#define __MPS_H
#include "combiner.h"
#include "iqcombiner.h"
#include "model.h"

enum Direction { Fromright, Fromleft, Both, None };

static const Real DefaultLogRef = 1.3E-14;

Vector do_denmat_Real(const ITensor& nA, ITensor& A, ITensor& B, Real cutoff,int minm, int maxm, Direction dir);
Vector do_denmat_Real(const IQTensor& nA, IQTensor& A, IQTensor& B, Real cutoff, int minm,int maxm, Direction dir);
Vector do_denmat_Real(const vector<IQTensor>& nA, const IQTensor& A, const IQTensor& B, IQTensor& U,
	Real cutoff,Real lognormfac,int minm,int maxm, Direction dir, bool donormalize, bool do_relative_cutoff);
void getCenterMatrix(ITensor& A, const Index& bond, Real cutoff,int minm, int maxm, ITensor& Lambda, string newbondname = "");
void plussers(const Index& l1, const Index& l2, Index& sumind, ITensor& first, ITensor& second);
void plussers(const IQIndex& l1, const IQIndex& l2, IQIndex& sumind, IQTensor& first, IQTensor& second);

template<class Tensor, class TensorSet>
Real doDavidson(Tensor& phi, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH, int niter, int debuglevel, Real errgoal);

extern Real truncerror, svdtruncerr;

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
        A_[1] = IQTensor(si(1),a[1],Center); A_[1](initState(1))=1;
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
    bool is_null() const { return (model_==0); }
    bool is_not_null() const { return (model_!=0); }

    Tensor bondTensor(int b) const { Tensor res = A.at(b) * A.at(b+1); return res; }

    //Useful for iterating over A's in a foreach loop
    //Do foreach(const I[Q]Tensor& A, psi.all_As()) { ... }
    const_range_type all_As() const { return make_pair(A.begin(),A.end()-1); }

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

    MPS(const ModelT& mod_, istream& s) : N(mod_.NN()), A(mod_.NN()+1), model_(&mod_)
    {
        for(int j = 1; j <= N; ++j) A[j].read(s);
        s.read((char*) &left_orth_lim,sizeof(left_orth_lim));
        s.read((char*) &right_orth_lim,sizeof(right_orth_lim));
        s.read((char*) &minm,sizeof(minm));
        s.read((char*) &maxm,sizeof(maxm));
        s.read((char*) &cutoff,sizeof(cutoff));
    }

    void write(ostream& s) const
    {
        //s.write((char*) &N,sizeof(N));
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
        A[left_orth_lim+1] *= a;
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

    void doSVD(int i, const Tensor& AA, Direction dir, bool preserve_shape = false)
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

        tensorSVD(AA,A[i],A[i+1],cutoff,minm,maxm,dir);
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


    void applygate(Tensor gate)
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
                        { rel_cut = max(fabs(summed_block.val1(j)),rel_cut); }
                        assert(rel_cut >= 0);
                        //Real rel_cut = summed_block.norm()/summed_block.vec_size();
                        rel_cut *= cut;
                        //cerr << "rel_cut == " << rel_cut << "\n";

                        if(rel_cut > 0)
                        for(int j = 1; j <= bond.m(); ++j)
                        if(fabs(summed_block.val1(j)) > rel_cut) 
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
                iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(si(s)),siP(s),linkind[s]) : IQTensor(si(s),linkind[s]));
            else if(s == N)
                iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(linkind[s-1]),conj(si(s)),siP(s)) 
                                        : IQTensor(conj(linkind[s-1]),si(s)));
            else
                iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(linkind[s-1]),conj(si(s)),siP(s),linkind[s]) 
                                        : IQTensor(conj(linkind[s-1]),si(s),linkind[s]));

            foreach(const ITensor& nb, nblock) { iqpsi.AAnc(s) += nb; } nblock.clear();

            //if(s < 5)
            //{
            //    Print(iqpsi.AA(s));
            //}
            //else Error("Stopping.");

            if(s==show_s)
            {
            iqpsi.AA(s).print((format("qA[%d]")%s).str(),ShowData);
            Error("Stopping");
            }

        } //for loop over s

        IQIndex Center("Center",Index("center",1,Virtual),totalq,In);
        iqpsi.AAnc(1).addindex1(Center);

        assert(check_QNs(iqpsi));

    } //void convertToIQ(IQMPSType& iqpsi) const

}; //class MPS<Tensor>

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
        AAnc(i) = first[i-1] * AA(i) * first[i] + second[i-1] * other.AA(i) * second[i];
    }
    AAnc(N) = first[N-1] * AA(N) + second[N-1] * other.AA(N);

    noprimelink();

    position(N);
    position(1);

    return *this;
}


} //namespace Internal
typedef Internal::MPS<ITensor> MPS;
typedef Internal::MPS<IQTensor> IQMPS;

void convertToIQ(const BaseModel& model, const vector<ITensor>& A, vector<IQTensor>& qA, QN totalq = QN(), Real cut = 1E-12)
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
                    { rel_cut = max(fabs(summed_block.val1(j)),rel_cut); }
                    assert(rel_cut >= 0);
                    //Real rel_cut = summed_block.norm()/summed_block.vec_size();
                    rel_cut *= cut;
                    //cerr << "rel_cut == " << rel_cut << "\n";

                    if(rel_cut > 0)
                    for(int j = 1; j <= bond.m(); ++j)
                    if(fabs(summed_block.val1(j)) > rel_cut) 
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

namespace Internal {

template<class Tensor>
class MPO
{
public:
    typedef Tensor TensorT;
    typedef typename Tensor::IndexT IndexT;
    typedef typename Tensor::IndexValT IndexValT;
    typedef typename Tensor::CombinerT CombinerT;
private:
    int N;
    vector<Tensor> A;
    int Lb,Rb;
    const BaseModel* model_;
    Real lref_;
public:
    int minm,maxm;
    Real cutoff;

    operator MPO<IQTensor>()
    { 
        MPO<IQTensor> res(*model_,maxm,cutoff,lref_); 
        res.minm = minm;
        convertToIQ(*model_,A,res.A);
        return res; 
    }

    //Accessor Methods ------------------------------

    int NN() const { return N;}
    // Norm of psi^2 = 1 = norm = sum of denmat evals. 
    // This translates to Tr{Adag A} = norm.  
    // Ref. norm is Tr{1} = d^N, d = 2 S=1/2, d = 4 for Hubbard, etc
    Real lref() const { return lref_; }
    void lref(Real val) 
	{  if(val == 0) { Error("bad lref"); } lref_ = val; }

    const BaseModel& model() const { return *model_; }
    IQIndex si(int i) const { return model_->si(i); }
    IQIndex siP(int i) const { return model_->siP(i); }

    int right_lim() const { return Rb; }
    int left_lim() const { return Lb; }

    const Tensor& AA(int i) const { return GET(A,i); }
    Tensor& AAnc(int i) //nc means 'non const'
    { 
        if(i <= Lb) Lb = i-1;
        if(i >= Rb) Rb = i+1;
        return GET(A,i); 
    }
    bool is_null() const { return (model_==0); }
    bool is_not_null() const { return (model_!=0); }

    Tensor bondTensor(int b) const { Tensor res = A.at(b) * A.at(b+1); return res; }

    MPO() : N(0), lref_(DefaultLogRef) { }
    MPO(const BaseModel& model, int maxm_ = MAX_M, Real cutoff_ = MAX_CUT, Real _lref = DefaultLogRef) 
    : N(model.NN()), A(N+1), Lb(0), Rb(N), model_(&model), lref_(_lref), minm(1), maxm(maxm_), cutoff(cutoff_)
	{ 
        if(_lref == 0) Error("MPO<Tensor>: Setting lref_ to zero");
        if(_lref == DefaultLogRef) lref_ = model.NN() * log(2.0); 
	}

    MPO(BaseModel& model, istream& s) : N(model.NN()), A(N+1), model_(&model)
    {
        for(int j = 1; j <= N; ++j) A[j].read(s);
        s.read((char*) &Lb,sizeof(Lb));
        s.read((char*) &Rb,sizeof(Rb));
        s.read((char*) &lref_,sizeof(lref_));
        s.read((char*) &minm,sizeof(minm));
        s.read((char*) &maxm,sizeof(maxm));
        s.read((char*) &cutoff,sizeof(cutoff));
    }

    void write(ostream& s) const
    {
        for(int j = 1; j <= N; ++j) A[j].write(s);
        s.write((char*) &Lb,sizeof(Lb));
        s.write((char*) &Rb,sizeof(Rb));
        s.write((char*) &lref_,sizeof(lref_));
        s.write((char*) &minm,sizeof(minm));
        s.write((char*) &maxm,sizeof(maxm));
        s.write((char*) &cutoff,sizeof(cutoff));
    }

    //MPO: index methods --------------------------------------------------

    void mapprime(int oldp, int newp, PrimeType pt = primeBoth)
	{ for(int i = 1; i <= N; ++i) A[i].mapprime(oldp,newp,pt); }

    void primelinks(int oldp, int newp)
	{ for(int i = 1; i <= N; ++i) A[i].mapprime(oldp,newp,primeLink); }

    void noprimelink()
	{ for(int i = 1; i <= N; ++i) A[i].noprime(primeLink); }

    IndexT LinkInd(int i) const { return index_in_common(A[i],A[i+1],Link); }
    IndexT RightLinkInd(int i) const { assert(i<NN()); return index_in_common(AA(i),AA(i+1),Link); }
    IndexT LeftLinkInd(int i)  const { assert(i>1); return index_in_common(AA(i),AA(i-1),Link); }

    void primeall()	// sites i,i' -> i',i'';  link:  l -> l'
	{
        for(int i = 1; i <= this->NN(); i++)
        {
            AAnc(i).mapprime(0,1,primeLink);
            AAnc(i).mapprime(1,2,primeSite);
            AAnc(i).mapprime(0,1,primeSite);
        }
	}

    void doSVD(int i, const Tensor& AA, Direction dir, bool preserve_shape = false)
	{
        tensorSVD(AA,A[i],A[i+1],cutoff,minm,maxm,dir,lref_);
        truncerror = svdtruncerr;

        if(dir == Fromleft)
        {
            if(Lb == i-1 || i == 1) Lb = i;
            if(Rb < i+2) Rb = i+2;
        }
        else
        {
            if(Lb > i-1) Lb = i-1;
            if(Rb == i+2 || i == N-1) Rb = i+1;
        }
	}

    //Move the orthogonality center to site i (left_orth_lim = i-1, right_orth_lim = i+1)
    void position(int i, bool preserve_shape = false)
	{
        if(is_null()) Error("position: MPS is null");
        while(Lb < i-1)
        {
            if(Lb < 0) Lb = 0;
            Tensor WF = AA(Lb+1) * AA(Lb+2);
            doSVD(Lb+1,WF,Fromleft,preserve_shape);
        }
        while(Rb > i+1)
        {
            if(Rb > N+1) Rb = N+1;
            Tensor WF = AA(Rb-2) * AA(Rb-1);
            doSVD(Rb-2,WF,Fromright,preserve_shape);
        }
	}

    bool is_ortho() const { return (Lb + 1 == Rb - 1); }

    int ortho_center() const 
    { 
        if(!is_ortho()) Error("MPS: orthogonality center not well defined.");
        return (Lb + 1);
    }

    bool is_complex() const
    { return A[Lb+1].is_complex(); }

    friend inline ostream& operator<<(ostream& s, const MPO& M)
    {
        s << "\n";
        for(int i = 1; i <= M.NN(); ++i) s << M.AA(i) << "\n";
        return s;
    }

    void print(string name = "",Printdat pdat = HideData) const 
    { printdat = (pdat==ShowData); cerr << "\n" << name << " =\n" << *this << "\n"; printdat = false; }

private:
    friend class MPO<ITensor>;
    friend class MPO<IQTensor>;
}; //class MPO<Tensor>
} //namespace Internal
typedef Internal::MPO<ITensor> MPO;
typedef Internal::MPO<IQTensor> IQMPO;


namespace Internal {

template<class Tensor>
class MPOSet
{
    int N, size_;
    vector<vector<const Tensor*> > A;
public:
    typedef vector<Tensor> TensorT;

    MPOSet() : N(-1), size_(0) { }

    MPOSet(const MPS<Tensor>& Op1) 
    : N(-1), size_(0) 
    { include(Op1); }

    MPOSet(const MPS<Tensor>& Op1, 
           const MPS<Tensor>& Op2) 
    : N(-1), size_(0) 
    { include(Op1); include(Op2); }

    MPOSet(const MPS<Tensor>& Op1, 
           const MPS<Tensor>& Op2,
           const MPS<Tensor>& Op3) 
    : N(-1), size_(0) 
    { include(Op1); include(Op2); include(Op3); }

    MPOSet(const MPS<Tensor>& Op1, 
           const MPS<Tensor>& Op2,
           const MPS<Tensor>& Op3, 
           const MPS<Tensor>& Op4) 
    : N(-1), size_(0) 
    { include(Op1); include(Op2); include(Op3); include(Op4); }

    void include(const MPS<Tensor>& Op)
    {
        if(N < 0) { N = Op.NN(); A.resize(N+1); }
        for(int n = 1; n <= N; ++n) GET(A,n).push_back(&(Op.AA(n))); 
        ++size_;
    }

    int NN() const { return N; }
    int size() const { return size_; }
    const vector<const Tensor*>& AA(int j) const { return GET(A,j); }
    const vector<Tensor> bondTensor(int b) const
    { vector<Tensor> res = A[b] * A[b+1]; return res; }

}; //class Internal::MPOSet

} //namespace Internal
typedef Internal::MPOSet<ITensor> MPOSet;
typedef Internal::MPOSet<IQTensor> IQMPOSet;


namespace Internal {

template <class Tensor>
class HamBuilder
{
public:
    typedef Tensor TensorT;
    typedef typename Tensor::IndexT IndexT;
    typedef BaseModel SiteSetT;
private:
    const SiteSetT& iss;
    const int N;
    vector<IndexT> currentlinks;
public:

    int NN() const { return N; }
    IndexT si(int i) const { return iss.si(i); }

    HamBuilder(const SiteSetT& iss_) : iss(iss_), N(iss.NN()), currentlinks(N) { }

    void newlinks(vector<Index>& currentlinks)
	{
        currentlinks.resize(N);
        static int ver = 0; ++ver;
        for(int i = 1; i < N; i++)
        {
            stringstream ss;
            ss << "h" << ver << "-" << i;
            currentlinks[i] = IndexT(ss.str(),1);
        }
	}

    void newlinks(vector<IQIndex>& currentlinks)
	{
        currentlinks.resize(N);
        static int ver = 0; ++ver;
        for(int i = 1; i < N; i++)
        {
            stringstream ss;
            ss << "h" << ver << "-" << i;
            currentlinks[i] = IndexT(ss.str());
        }
	}

    Tensor unit(int i) const { return Tensor(si(i),si(i).primed(),1); }

    void getidentity(Real factor, MPO<ITensor>& res)
	{
        newlinks(currentlinks);
        res = MPO<ITensor>(iss,res.maxm,res.cutoff);
        res.AAnc(1) = unit(1); res.AAnc(1).addindex1(GET(currentlinks,1));
        res.AAnc(1) *= factor;
        res.AAnc(N) = unit(N); res.AAnc(N).addindex1(GET(currentlinks,N-1));
        for(int i = 2; i < N; ++i)
        {
            res.AAnc(i) = unit(i); 
            res.AAnc(i).addindex1(GET(currentlinks,i-1));
            res.AAnc(i).addindex1(GET(currentlinks,i));
        }
	}

    void getMPO(Real factor, int i, Tensor op, MPO<ITensor>& res)
	{
        getidentity(1,res);
        res.AAnc(i) = op;
        if(i > 1) res.AAnc(i).addindex1(GET(currentlinks,i-1));
        if(i < N) res.AAnc(i).addindex1(GET(currentlinks,i));
        res *= factor;
	}

    void getMPO(Real factor, int i1, Tensor op1, int i2, Tensor op2, MPO<ITensor>& res)
	{
        if(i1 == i2) Error("HamBuilder::getMPO: i1 cannot equal i2.");
        getMPO(1,i2,op2,res);
        res.AAnc(i1) = op1;
        if(i1 > 1) res.AAnc(i1).addindex1(GET(currentlinks,i1-1));
        if(i1 < N) res.AAnc(i1).addindex1(GET(currentlinks,i1));
        res *= factor;
	}

    template <typename Iterable1, typename Iterable2>
    void getMPO(Real factor, Iterable1 sites, Iterable2 ops, MPO<ITensor>& res)
	{
        for(int i = 0; i < (int) sites.size(); ++i)
        for(int j = 0; j < (int) sites.size(); ++j)
        {
            if(i == j) continue;
            if(sites[i] == sites[j]) Error("HamBuilder::getMPO: all sites should be unique.");
        }

        if(sites.size() != ops.size()) Error("HamBuilder::getMPO: need same number of sites as ops.");

        getidentity(1,res);
        for(int i = 0; i < (int) sites.size(); ++i)
        {
            const int s = sites[i];
            res.AAnc(s) = GET(ops,i);
            assert(GET(ops,i).hasindex(si(sites[i])));
            if(s > 1) res.AAnc(s).addindex1(GET(currentlinks,s-1));
            if(s < N) res.AAnc(s).addindex1(GET(currentlinks,s));
        }
        res *= factor;
	}

};

} //namespace Internal
typedef Internal::HamBuilder<ITensor> HamBuilder;
//typedef Internal::HamBuilder<IQTensor> IQHamBuilder;

class MPOBuilder
{
protected:
    const SiteSet& sst;
    const int Ns;
public:

    MPOBuilder(const SiteSet& sst_) : sst(sst_), Ns(sst_.NN()) { }

    ITensor makeLedge(const Index& L) const
    {
        ITensor res(L); res(L(L.m())) = 1;
        return res;
    }

    ITensor makeRedge(const Index& R) const
    {
        ITensor res(R); res(R(1)) = 1;
        return res;
    }
};

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
    if(active.primelevel != 0)
    {
        Print(rho.index(1));
        Print(rho.index(2));
        Print(active);
    }
    assert(active.primelevel == 0);

    vector<Matrix> mmatrix(rho.iten_size());
    vector<Vector> mvector(rho.iten_size());
    vector<Real> alleig;

    Real maxlogfac = -1E20;
    if(logrefnorm == DefaultLogRef) 
    {
        foreach(const ITensor& t, rho.itensors())
        { maxlogfac = max(maxlogfac,t.logfac()); }
    }
    else { maxlogfac = 2.0*logrefnorm; }

    //1. Diagonalize each ITensor within rho
    int itenind = 0;
    foreach(const ITensor& t, rho.itensors())
	{
        //assert(t.index(1).noprime_equals(t.index(2)));
        if(!t.index(1).noprime_equals(t.index(2)))
        { Print(rho); Print(t); Error("Non-symmetric ITensor in density matrix"); }

        t.normlogto(maxlogfac);

        Matrix &U = mmatrix[itenind];
        Vector &d = mvector[itenind];

        //Diag ITensors within rho
        int n = t.index(1).m();
        Matrix M(n,n);
        t.toMatrix11(t.index(2),t.index(1),M);

        M *= -1.0;
        EigenValues(M,d,U);
        d *= -1.0;

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
    if(((svdtruncerr += alleig[m]/e1) > cutoff && mkeep <= maxm) || mkeep <= minm)
    { 
        docut = (m > 0 ?  (alleig[m-1] + alleig[m]) * 0.5 : 0);
        svdtruncerr -= alleig[m]/e1;
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
        const Vector& thisD = mvector[itenind];
        int this_m = 1;
        for(; this_m <= thisD.Length(); ++this_m)
        if(thisD(this_m) < docut) { break; }
        --this_m; //since for loop overshoots by 1

        if(this_m == 0) { ++itenind; continue; }

        if(mkeep == 0 && thisD.Length() >= 1) // zero mps, just keep one arb state
        { this_m = 2; mkeep = 1; docut = 1; }

        Index nm("qlink",this_m);
        Index act = t.index(1).deprimed();
        iq.push_back(inqn(nm,active.qn(act)));

        Matrix UU = mmatrix[itenind].Columns(1,this_m);

        assert(act.primelevel == 0);
        assert(active.hasindex(act));
        assert(act.m() == UU.Nrows());

        terms.push_back(ITensor(act,nm,UU));
        ++itenind;
	}
    //Print(iq.size());
    //Print(iq);
    //cerr << "Got here\n";
    IQIndex newmid("qlink",iq,In);
    nU = IQTensor(active,newmid);
    foreach(const ITensor& t, terms) nU += t;

    eigs_kept.ReDimension(m);
    for(int i = 1; i <= m; ++i) eigs_kept(i) = alleig[alleig.size()-i];
    lastd = eigs_kept;
} //void diag_denmat

template<class Tensor>
Vector tensorSVD(const Tensor& AA, Tensor& A, Tensor& B, Real cutoff, int minm, int maxm, Direction dir, Real logrefnorm = DefaultLogRef)
{
    typedef typename Tensor::IndexT IndexT;
    typedef typename Tensor::CombinerT CombinerT;

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

    to_orth = U * conj(comb);
    newoc   = conj(U) * AAc;

    return eigs_kept;
}

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
                cerr << "check_QNs: Right side Link not pointing In to left of orthog. center at site " << i << endl;
                cerr << "AA(" << i << ") = " << psi.AA(i) << endl;
                return false;
            }
            if(i > 1)
            {
                if(psi.LeftLinkInd(i).dir() != Out) 
                {
                    cerr << "check_QNs: Left side Link not pointing Out to left of orthog. center at site " << i << endl;
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
                cerr << "check_QNs: Right side Link not pointing Out on right side of orthog. center at site " << i << endl;
                return false;
            }
            if(psi.LeftLinkInd(i).dir() != In) 
            {
                cerr << "check_QNs: Left side Link not pointing In on right side of orthog. center at site " << i << endl;
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

template <class MPSType, class MPOType>
void psiHphi(const MPSType& psi, const MPOType& H, const MPSType& phi, Real& re, Real& im) //<psi|H|phi>
{
    typedef typename MPSType::TensorT Tensor;
    const int N = H.NN();

    if(psi.NN() != phi.NN() || psi.NN() != N) Error("psiHphi: mismatched N");

    Tensor L = phi.AA(1) * H.AA(1) * conj(primed(psi.AA(1)));
    for(int i = 2; i < N; ++i) { L = L * phi.AA(i) * H.AA(i) * conj(primed(psi.AA(i))); }
    L = L * phi.AA(N) * H.AA(N);

    Dot(primed(psi.AA(N)),L,re,im);
}
template <class MPSType, class MPOType>
Real psiHphi(const MPSType& psi, const MPOType& H, const MPSType& phi) //Re[<psi|H|phi>]
{
    Real re, im;
    psiHphi(psi,H,phi,re,im);
    if(im != 0) cerr << format("\nReal psiHphi: WARNING, dropping non-zero (im = %.5f) imaginary part of expectation value.\n")%im;
    return re;
}

void psiHphi(const MPS& psi, const MPO& H, const ITensor& LB, const ITensor& RB, const MPS& phi, Real& re, Real& im);
Real psiHphi(const MPS& psi, const MPO& H, const ITensor& LB, const ITensor& RB, const MPS& phi); //Re[<psi|H|phi>]

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
int specialdebug1 = 0, specialdebug2 = 0;
#endif

#endif
