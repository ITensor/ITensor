//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "autompo.h"
#include <algorithm>

using std::find;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::array;
using std::pair;
using std::make_pair;

namespace itensor {

bool
isReal(const Complex& z) { return z.imag() == 0; }

bool
isApproxReal(const Complex& z, Real epsilon = 1E-12) { return fabs(z.imag()) < epsilon; }

SiteTerm::
SiteTerm() : i(-1) { }

SiteTerm::
SiteTerm(const std::string& op_,
         int i_)
    :
    op(op_),
    i(i_)
    { }

bool SiteTerm::
operator==(const SiteTerm& other) const
    {
    return (op == other.op && i == other.i);
    }

bool SiteTerm::
isFermionic() const
    {
    if(!op.empty() && op.front() == 'C') return true;
    return false;
    }

std::string OpString(const SiteTermProd &prod)
    {
    if(prod.empty())
        return "";

    std::string opstr;
    for(auto t = prod.begin(); t != prod.end()-1; t++)
        opstr += t->op + "*";
    opstr += prod.back().op;
    return opstr;
    }
    
bool IsFermionic(const SiteTermProd &prod)
    {
    std::string opstr = OpString(prod);
    auto num = std::count(opstr.begin(), opstr.end(), 'C');
    return num % 2;
    }

string
fermionicTerm(const std::string& op)
    {
    static array<pair<string,string>,6>
           rewrites =
           {{
           make_pair("Cdagup","Adagup"),
           make_pair("Cup","Aup"),
           make_pair("Cdagdn","Adagdn*Fup"),
           make_pair("Cdn","Adn*Fup"),
           make_pair("C","A"),
           make_pair("Cdag","Adag")
           }};
    for(auto& p : rewrites)
        {
        if(p.first == op) return p.second;
        }
    return op;
    }
    
void RewriteFermionic(SiteTermProd &prod, bool isleftFermionic)
    {
    if(prod.empty())
        Error("Empty product in RewriteFermionic is not expected.");    
    
    int i = prod.front().i;
    for(const SiteTerm &t : prod)
        if(t.i != i)
            Error("Multi-site product in RewriteFermionic is not expected.");    

    // Rewrite a fermionic single site product using the Jordan-Wigner string            
    bool isSiteFermionic = IsFermionic(prod);
    if(isSiteFermionic)
        {
        for(SiteTerm &t : prod)
            if(t.isFermionic())
                t.op = fermionicTerm(t.op);
        }
    
    // Add a FermiPhase operator at the end if the product of operators
    // to the left (including this site) is fermionic
    if((isleftFermionic && !isSiteFermionic) || (!isleftFermionic && isSiteFermionic))
        prod.emplace_back("F", i);         
    }

QN QuantumNumber(const SiteSet &sites, const SiteTermProd &prod)
    {
    QN qn;
    for(const SiteTerm &t: prod)
        {
        IQTensor op = sites.op(t.op,t.i);
        qn += -div(op);                    
        }
    return qn;
    }
    
bool SiteTermProd::operator<(const SiteTermProd &other) const
    { 
    if(size() != other.size())
        return size() < other.size();
        
    for(unsigned i=0; i<size(); i++)
        if(at(i).i != other.at(i).i)
            return at(i).i < other.at(i).i;
        else if(at(i).op != other.at(i).op)
            return at(i).op < other.at(i).op;
            
    return false; // the two are equal
    }
    
SiteTermProd mult(const SiteTermProd &first, const SiteTermProd &second)
    {
    SiteTermProd prod = first;
    prod.insert( prod.end(), second.begin(), second.end() );
    return prod;
    }    

Term& Term::
operator*=(Real x)
    {
    coef *= x;
    return *this;
    }

Term& Term::
operator*=(Complex x)
    {
    coef *= x;
    return *this;
    }
    
Term Term::
operator*(Real x) const
    {
    Term t(x*this->coef, this->ops);
    return t;
    }

Term Term::
operator*(Complex x) const
    {
    Term t(x*this->coef, this->ops);
    return t;
    }

void TermSum::operator+=(const Term &term)
    {
    // Check if the Term already appears in the sum and if yes, then only add the coefficients
    SiteTermProd ops = term.ops;
    auto isEqualProd = [&ops](const Term &t) { return t.ops == ops; };
    auto it = find_if(sum.begin(), sum.end(), isEqualProd);
    if(it != sum.end())
        it->coef += term.coef;
    else
        sum.push_back(term);
    }    

void HTerm::
add(const std::string& op,
    int i,
    Real x)
    {
    //The following ensures operators remain
    //in site order within the vector "ops"
    auto it = ops.begin();
    while(it != ops.end() && it->i <= i) ++it;
    
    SiteTerm t(op,i);
    
    // If the operator is fermionic and being inserted in between existing operators 
    // need to check if an extra minus is required
    if(it != ops.end() && t.isFermionic())
    {
        SiteTermProd rightOps(it, ops.end());
        if(IsFermionic(rightOps))
            coef *= -1;
    }
    
    coef *= x;
    ops.insert(it,t);
    }

bool HTerm::
startsOn(int i) const 
    { 
    if(ops.empty()) Error("No operators in HTerm");
    return first().i == i; 
    }

bool HTerm::
endsOn(int i) const 
    { 
    if(ops.empty()) Error("No operators in HTerm");
    return last().i == i; 
    }

bool HTerm::
contains(int i) const 
    { 
    if(ops.empty()) Error("No operators in HTerm");
    return i >= first().i && i <= last().i; 
    }

bool HTerm::
proportionalTo(const HTerm& other) const
    {  
    return ops == other.ops;
    }
    
void AutoMPO::
add(const HTerm& t)
    { 
    if(abs(t.coef) == 0)
        return; 
        
    terms_[t.ops] += t.coef;
    }


AutoMPO::Accumulator::
Accumulator(AutoMPO* pa_, 
            Real x_)
    :
    pa(pa_),
    state(New),
    coef(x_)
    {}

AutoMPO::Accumulator::
Accumulator(AutoMPO* pa_, 
            Complex x_)
    :
    pa(pa_),
    state(New),
    coef(x_)
    {}

AutoMPO::Accumulator::
Accumulator(AutoMPO* pa_)
    : 
    Accumulator(pa_,1)
    {}


AutoMPO::Accumulator::
Accumulator(AutoMPO* pa_, 
            const char* op_)
    :
    pa(pa_),
    state(Op),
    coef(1),
    op(op_)
    {}

AutoMPO::Accumulator::
Accumulator(AutoMPO* pa_, 
            const std::string& op_)
    :
    pa(pa_),
    state(Op),
    coef(1),
    op(op_)
    {}


AutoMPO::Accumulator::
~Accumulator()
    {
    if(state==Op) Error("Invalid input to AutoMPO (missing site number?)");
    term *= coef;
    pa->add(term);
    }
    

AutoMPO::Accumulator& AutoMPO::Accumulator::
operator,(Real x)
    {
    coef *= x;
    return *this;
    }

AutoMPO::Accumulator& AutoMPO::Accumulator::
operator,(Complex x)
    {
    coef *= x;
    return *this;
    }

AutoMPO::Accumulator& AutoMPO::Accumulator::
operator,(int i)
    {
    if(state==Op)
        {
        term.add(op,i);
        state = New;
        op = "";
        }
    else
        {
        coef *= Real(i);
        }
    return *this;
    }

AutoMPO::Accumulator& AutoMPO::Accumulator::
operator,(const char* op_)
    {
    if(state == New)
        {
        op = op_;
        state = Op;
        }
    else
        {
        Error("Invalid input to AutoMPO (two strings in a row?)");
        }
    return *this;
    }

AutoMPO::Accumulator& AutoMPO::Accumulator::
operator,(const std::string& op_)
    {
    if(state == New)
        {
        op = op_;
        state = Op;
        }
    else
        {
        Error("Invalid input to AutoMPO (two strings in a row?)");
        }
    return *this;
    }
    
bool IQMPOMatElement::operator==(const IQMPOMatElement &other) const
    {
        return rowqn == other.rowqn && colqn == other.colqn && 
                row == other.row && col == other.col && 
                val == other.val;
    }

ComplexMatrix::ComplexMatrix(const std::vector<CoefMatElement> &M)
    {
    int nr = 0, nc = 0;
    bool isComplex = false;
    
    for(const CoefMatElement &elem : M)
        {
        nr = max(nr, elem.ind.row);
        nc = max(nc, elem.ind.col);
        if(!isReal(elem.val))
            isComplex = true;
        }
    
    Re.Enlarge(nr, nc);
    if(isComplex)
        Im.Enlarge(nr, nc);
        
    for(const CoefMatElement &elem : M)
        {
        Re(elem.ind.row,elem.ind.col) = elem.val.real();
        if(!isReal(elem.val))
            Im(elem.ind.row,elem.ind.col) = elem.val.imag();
        }
    }
    
Complex ComplexMatrix::operator() (int i, int j) const
    {
        Complex re, im;
        if(Im.Storage())
            im = Complex(0, Im(i,j));
        if(Re.Storage())
            re = Complex(Re(i,j),0);
        return re+im;
    }

void AutoMPO::DecomposeTerm(int n, const SiteTermProd &ops, 
                            SiteTermProd &left, SiteTermProd &onsite, SiteTermProd &right) const
    {
    auto isOnSiteOrOnTheRight = [&n](const SiteTerm &t) {return t.i >= n;};
    auto startOfOnSite = find_if(ops.begin(), ops.end(), isOnSiteOrOnTheRight);
    
    auto isOnTheRight = [&n](const SiteTerm &t) {return t.i > n;};
    auto startOfRightPart = find_if(startOfOnSite, ops.end(), isOnTheRight);

    left = SiteTermProd(ops.begin(), startOfOnSite);
    onsite = SiteTermProd(startOfOnSite, startOfRightPart);
    right = SiteTermProd(startOfRightPart, ops.end());
    }  

// Returns a 1-based index of the SiteTermProd ops in the vector
// If ops is not in the vector adds it is added
int AutoMPO::PosInVec(const SiteTermProd &ops, std::vector<SiteTermProd> &vec) const
    {   
    auto it = find(vec.begin(), vec.end(), ops);
    if (it != vec.end())
        return it - vec.begin() + 1;
    vec.push_back(ops);
    return vec.size();
    }

// Translate terms map into a vector for backwards compatibility
std::vector<HTerm> AutoMPO::terms() const
    {
        std::vector<HTerm> hterms;
        hterms.reserve(terms_.size());
    
        for(const std::pair<SiteTermProd, Complex> &t : terms_)
            hterms.emplace_back(HTerm(t.second, t.first));
        
        return hterms;
    }

// Construct left & right partials and the ceofficients matrix on each link as well as the temporary MPO
void AutoMPO::PartitionHTerms(std::vector<PartitionByQN> &part, std::vector<MPOSparseMatrix> &tempMPO) const
    {
    std::vector<HTerm> hterms = terms();
    for(const HTerm &ht : hterms)
        for(int n = ht.first().i; n <= ht.last().i; n++)
            {
            SiteTermProd left, onsite, right;
            DecomposeTerm(n, ht.ops, left, onsite, right);
            
            QN lqn = QuantumNumber(sites_, left);
            QN sqn = QuantumNumber(sites_, onsite);
            
            int j,k,l;

            // part.at(i) is the partition at the link between sites i+1 and i+2
            // i.e. part.at(0) is the partition at the link between sites 1 and 2
            // and part.at(N-2) is the partition at the link between sites N-1 and N
            // for site n the link on the left is part.at(n-2) and the link on the right is part.at(n-1)
            if(left.empty())
                {
                j=0;
                if(right.empty()) // on site term
                    k = 0;
                else // term starting on site n
                    k = PosInVec(right, part.at(n-1)[sqn].right);
                }
            else
                {
                if(right.empty()) // term ending on site n
                    {
                    k = 0;
                    j = PosInVec(onsite, part.at(n-2)[lqn].right);
                    }
                else
                    {
                    j = PosInVec(mult(onsite,right), part.at(n-2)[lqn].right);
                    k = PosInVec(right, part.at(n-1)[lqn+sqn].right);
                    }
                l = PosInVec(left, part.at(n-2)[lqn].left);
                part.at(n-2)[lqn].Coeff.emplace_back(CoefMatElement(MatIndex(l, j), ht.coef));
                }
                
            // Place the coefficient of the HTerm when the term starts
            Complex c = j==0 ? ht.coef : 1;
            
            bool leftF = IsFermionic(left);
            if(onsite.empty())
                {
                if(leftF)
                    onsite.emplace_back(SiteTerm("F",n));
                else 
                    onsite.emplace_back(SiteTerm("Id",n));
                }
            else
                RewriteFermionic(onsite, leftF);
            
            IQMPOMatElement elem(lqn, lqn+sqn, j, k, Term(c, onsite));
            auto el = find(tempMPO.at(n-1).begin(), tempMPO.at(n-1).end(), elem);
            if(el == tempMPO.at(n-1).end())
                tempMPO.at(n-1).push_back(elem);
            }    
            
#ifdef SHOW_AUTOMPO

    println("Left and Right Partials:");
    for(unsigned n=0; n<part.size(); n++)
        {
        for(const auto &pqn : part.at(n))
            {
            const Partition &p = pqn.second;

            println("Left, QN = ", pqn.first);
            for(const SiteTermProd &prod : p.left)
                println(prod);
            println("Right, QN = ", pqn.first);
            for(const SiteTermProd &prod : p.right)
                println(prod);
            println("Coef, QN = ", pqn.first);
            for(const CoefMatElement &elem : p.Coeff)
                println(elem.ind.row,',',elem.ind.col,'\t',elem.val);
            }
        println("=========================================");
        }
        
    println();
    println("TempMPO Elements:");
    for(unsigned n=0; n<tempMPO.size(); n++)
        {
        for(const IQMPOMatElement &elem: tempMPO.at(n))
            println(elem.rowqn,',',elem.row,'\t',elem.colqn,',',elem.col,'\t',elem.val.coef,'\t',elem.val.ops);
        println("=========================================");
        }

#endif   
    }

// SVD the coefficients matrix on each link and construct the compressed MPO matrix
void AutoMPO::CompressMPO(const std::vector<PartitionByQN> &part, const std::vector<MPOSparseMatrix> &tempMPO,
                        std::vector<MPOMatrix> &finalMPO, std::vector<IQIndex> &links, 
                        bool isExpH = false, Complex tau = 0) const
    {
    const int N = sites_.N();
    
    // For the MPO matrix on site n we need the SVD on both the previous link and the following link
    
    std::map<QN, ComplexMatrix> V_n, V_npp;
    std::map<QN, int> d_n, d_npp; 	// num of non-zero singular values for each QN block
    int d_n_tot = 0, d_npp_tot = 0;
    int max_d = 0;
    
    const QN ZeroQN;
    
    int d0 = isExpH ? 1 : 2;
    
    d_n[ZeroQN] = 0;
    d_n_tot = d0;
    
    // qn block offset in the compressed MPO
    std::map<QN, int> qnstart_n, qnstart_npp;
    qnstart_n[ZeroQN] = d0; 
        
    vector<IndexQN> inqn;
    inqn.emplace_back(Index("hl0_0",d_n_tot),ZeroQN);
    links.at(0) = IQIndex("Hl0",inqn);

    for(int n=1; n<=N; n++)
        {
        if(n==N || part.at(n-1).empty())
            d_npp[ZeroQN] = 0;        
        else    
            for(const std::pair<QN, Partition> &v : part.at(n-1) )
                {
                QN qn = v.first;
                Partition p = v.second;
                
                // Convert the coefficients of the partition to a dense Matrix                
                ComplexMatrix C(p.Coeff);
                
                Vector D;
                if(C.isComplex())
                    {                    
                    ComplexMatrix U;
                    SVD(C.Re, C.Im, U.Re, U.Im, D, V_npp[qn].Re, V_npp[qn].Im);
                    }
                else
                    {
                    Matrix U;                
                    SVD(C.Re, U, D, V_npp[qn].Re);
                    }

                auto isApproxZero = [&](const Real &val){ return fabs(val) < epsilon_; };
                auto firstApproxZero = std::find_if(D.begin(), D.end(), isApproxZero);
                d_npp[qn]=firstApproxZero - D.begin();
                }

        qnstart_npp[ZeroQN] = d0;
        
        inqn.clear();
        int count = 0;

        d_npp_tot = d_npp[ZeroQN]+d0;
        
        // Make sure zero QN is first in the list of indices
        inqn.emplace_back(Index(format("hl%d_%d",n,count++),d_npp_tot),ZeroQN);        
        
        for(const std::pair<QN, int> &d : d_npp)
            {
            if(d.first == ZeroQN)
                continue;   // was already taken care of
            
            qnstart_npp[d.first] = d_npp_tot;
            d_npp_tot += d.second;
            
            inqn.emplace_back(Index(format("hl%d_%d",n,count++),d.second),d.first);
            }
            
        links.at(n) = IQIndex(nameint("Hl",n),inqn);

#ifdef SHOW_AUTOMPO        
        println("Size of MPO matrix for site ", n, " is (", d_n_tot, ", ", d_npp_tot, ")");
#endif       
       
        // Construct the compressed MPO
        
        finalMPO.at(n-1).resize(d_n_tot);
        for(auto &v : finalMPO.at(n-1))
            v.resize(d_npp_tot);
        
        SiteTermProd opId;
        opId.emplace_back("Id", n);
        
        finalMPO.at(n-1).at(0).at(0) += Term(1, opId);
        if(!isExpH)
            finalMPO.at(n-1).at(1).at(1) += Term(1, opId);
            
        Complex Zero(0,0);
        for(const IQMPOMatElement &elem: tempMPO.at(n-1))
            {
            int k = elem.row;
            int l = elem.col;
            
            // if constructing ExpH multiply by tau when term is starting (k=0)
            Term t = (isExpH && k==0) ? elem.val*(-tau) : elem.val;
            int rowOffset = isExpH ? 0 : 1;
    
            if(l==0 && k==0)	// on-site terms
                finalMPO.at(n-1).at(rowOffset).at(0) += t;
            else if(k==0)  	// terms starting on site n
                {
                for(int j=1; j<=d_npp[elem.colqn]; j++)
                    if(V_npp[elem.colqn](j,l) != Zero) // 1-based access of matrix elements
                        finalMPO.at(n-1).at(rowOffset).at(qnstart_npp[elem.colqn]+j-1) += Term(t.coef*conj(V_npp[elem.colqn](j,l)), t.ops);
                        
                }
            else if(l==0) 	// terms ending on site n
                {
                for(int i=1; i<=d_n[elem.rowqn]; i++)
                    if(V_n[elem.rowqn](i,k) != Zero) // 1-based access of matrix elements
                        finalMPO.at(n-1).at(qnstart_n[elem.rowqn]+i-1).at(0) += Term(t.coef*V_n[elem.rowqn](i,k), t.ops);
                }
            else 
                {
                for(int i=1; i<=d_n[elem.rowqn]; i++)
                    for(int j=1; j<=d_npp[elem.colqn]; j++) 
                        if( (V_n[elem.rowqn](i,k) != Zero)  && (V_npp[elem.colqn](j,l) != Zero) ) // 1-based access of matrix elements
                            finalMPO.at(n-1).at(qnstart_n[elem.rowqn]+i-1).at(qnstart_npp[elem.colqn]+j-1) += Term(t.coef*V_n[elem.rowqn](i,k)*conj(V_npp[elem.colqn](j,l)), t.ops);
                }
            }

        // Store SVD computed at this step for next link
        V_n = V_npp;
        d_n = d_npp;
        d_n_tot = d_npp_tot;
        qnstart_n = qnstart_npp;
        
        V_npp.clear();
        d_npp.clear();
        qnstart_npp.clear();
        
        max_d = max(max_d, d_n_tot);
        }
        
        println("Maximal dimension of the MPO is ", max_d);
        
#ifdef SHOW_AUTOMPO
        for(int n=1; n<=N; n++)
            {
            for(unsigned r = 0; r < finalMPO.at(n-1).size(); ++r, println())
                for(unsigned c = 0; c < finalMPO.at(n-1).at(r).size(); ++c)
                    print(finalMPO.at(n-1).at(r).at(c), "\t");
            println("=========================================");
            }
#endif            
    }

IQMPO AutoMPO::ConstructMPOTensors(const std::vector<MPOMatrix> &finalMPO, 
                                    const std::vector<IQIndex> &links, bool isExpH = false) const
    {
    
    IQMPO H(sites_);
    
    const int N = sites_.N();
    int min_n = isExpH ? 1 : 2;
    
    for(int n=1; n<=N; n++)
        {
        int nr = finalMPO.at(n-1).size();
        int nc = n == N ? min_n : finalMPO.at(n).size();
        
        auto &row = links.at(n-1),
             &col = links.at(n);

        H.Anc(n) = IQTensor(dag(sites_(n)),prime(sites_(n)),dag(row),col);

        for(int r = 1; r <= nr; ++r)
            for(int c = 1; c <= nc; ++c)
                for(const Term &term : finalMPO.at(n-1).at(r-1).at(c-1).sum)
                    {
                    if(abs(term.coef) < 1E-12)
                        continue;
                        
                    IQTensor op = sites_.op(term.ops.front().op, n);
                    for(auto it = term.ops.begin()+1; it != term.ops.end(); it++)
                        op = multSiteOps(op, sites_.op(it->op, n));

                    if(isReal(term.coef))        
                        H.Anc(n) += term.coef.real() * op * row(r) * col(c);       
                    else
                        H.Anc(n) += term.coef * op * row(r) * col(c);
                    }
        }
    
    H.Anc(1) *= IQTensor(links.at(0)(min_n));
    H.Anc(N) *= IQTensor(dag(links.at(N))(1));   
    
    return H;
    }

IQMPO AutoMPO::ConstructMPOUsingSVD() const
    {
    const int N = sites_.N();
    
    std::vector<PartitionByQN> part(N-1);   // There are N-1 links between N sites
    std::vector<MPOSparseMatrix> tempMPO(N);

    PartitionHTerms(part, tempMPO);        
    
    std::vector<MPOMatrix> finalMPO(N);
    std::vector<IQIndex> links(N+1);

    CompressMPO(part, tempMPO, finalMPO, links);

    return ConstructMPOTensors(finalMPO, links);
    }

IQMPO AutoMPO::toExpHUsingSVD_ZW1(Complex tau) const
    {
    const int N = sites_.N();
    
    std::vector<PartitionByQN> part(N-1); // There are N-1 links between N sites
    std::vector<MPOSparseMatrix> tempMPO(N);

    PartitionHTerms(part, tempMPO);        
    
    std::vector<MPOMatrix> finalMPO(N);
    std::vector<IQIndex> links(N+1);
    
    CompressMPO(part, tempMPO, finalMPO, links, /*isExpH*/ true, tau);

    return ConstructMPOTensors(finalMPO, links, /*isExpH*/ true);
    }

/*
MPO convention:
===============
For each link of the MPO, define a set of bases 
that describe the terms of the Hamiltonian
corresponding to the left "half" of the MPO.
The terms include "IL", which means the product
of identities to the left, and "HL", the sum of
all terms entirely contained on the left.

Usually these two special terms occupy positions 1 and two,
respectively.

The rest of the bases are each site term on the left that
is connected to something on the right.

So for neighbor and next neighbor, operator pair A B, 
coefs t1 and t2, on site n, the MPO matrix is:
n-1             n
      1111   HL  11A1  111A  <== bases
1111   1     0     0    A         
HL     0     1     0    0   
11A1   0    t2 B   0    0   
111A   0    t1 B   1    0   

For neighbor and next neighbor, operator pair A B and B A, t1 and t2
site n:
n-1             n
      1111  HL    11A1 11B1  111A  111B 
1111   1     0     0     0     A     B  
HL     0     1     0     0     0     0  
11A1   0    t2 B   0     0     0     0  
11B1   0    t2 A   0     0     0     0  
111A   0    t1 B   1     0     0     0  
111B   0    t1 A   0     1     0     0  

F == fermiPhase, i.e. F = (-1)^(# of fermions of either type of spin)
Then we make c and cdagger both have F's going off to the left.

Fermion operator rewriting convention:

//
//Spinless fermions
//

Cdag_i C_j  = (F_1 F_2 F_3 ... F_{i-1})^2 (Adag_i F_i) F_{i+1} ... A_j
            = Adag_i F_{i+1} ... A_j

C_i Cdag_j = (A_i F_i) F_{i+1} ... Adag_j

//
//Fermions with spin
//

Cdagup_i Cup_j  = (F_1 F_2 F_3 ... )^2 (Adagup_i F_i) F_{i+1} ... Aup_j
                = (Adagup_i F_i) F_{i+1} ... Aup_j //cancel squared F operators

Cup_i Cdagup_j = (Aup_i F_i) F_{i+1} ... Adagup_j

Cdagdn_i Cdn_j  = (Adagdn_i F_i) F_{i+1} ... Fup_j Adn_j 
                = - Adagdn_i F_{i+1} ... Fup_j Adn_j     //use Adagdn_i * F_i = -Adagdn_i
                = Adagdn_i F_{i+1} ... Fup_j Fdn_j Adn_j //use Adn_j = -Fdn_j*Adn_j
                = Adagdn_i F_{i+1} ... (F_j Adn_j)       //combine Fup_j*Fdn_j = F_j (definition)

Cdn_i Cdagdn_j = (Adn_i F_i) F_{i+1} ... Fup_j Adagdn_j
               = - Adn_i F_{i+1} ... Fup_j Adagdn_j      //use Adn_i*F_i = -Adn_i
               = Adn_i F_{i+1} ... Fup_j Fdn_j Adagdn_j  //use Adagdn_j = -Fdn_j*Adagdn_j
               = Adn_i F_{i+1} ... (F_j Adagdn_j)        //combined Fup_j*Fdn_j = F_j (definition)


*/

string
startTerm(const std::string& op)
    {
    static array<pair<string,string>,6>
           rewrites =
           {{
           make_pair("Cdagup","Adagup*F"),
           make_pair("Cup","Aup*F"),
           make_pair("Cdagdn","Adagdn"),
           make_pair("Cdn","Adn"),
           make_pair("C","A*F"),
           make_pair("Cdag","Adag")
           }};
    for(auto& p : rewrites)
        {
        if(p.first == op) return p.second;
        }
    return op;
    }

string
endTerm(const std::string& op)
    {
    static array<pair<string,string>,6>
           rewrites =
           {{
           make_pair("Cup","Aup"),
           make_pair("Cdagup","Adagup"),
           make_pair("Cdn","F*Adn"),
           make_pair("Cdagdn","F*Adagdn"),
           make_pair("C","A"),
           make_pair("Cdag","Adag")
           }};
    for(auto& p : rewrites)
        {
        if(p.first == op) return p.second;
        }
    return op;
    }

struct SiteQN
    {
    SiteTerm st;
    QN q;
    SiteQN() { }
    SiteQN(const SiteTerm& st_,
           const QN& q_)
        :
        st(st_),
        q(q_)
        {
        }
    };

void
plusAppend(std::string& s, const std::string& a)
    {
    if(s.size() == 0 || s == "0") s = a;
    else 
        {
        s += "+";
        s += a;
        }
    }
    
template<>
IQMPO
toMPO<IQTensor>(const AutoMPO& am,
                const Args& args)
    {
    const SiteSet& sites = am.sites();
    IQMPO H(sites);
    const int N = sites.N();
    const auto checkqn = args.getBool("CheckQNs",true);
    
    std::vector<HTerm> terms = am.terms();

    for(auto& t : terms)
        if(t.Nops() > 2) 
            Error("Only at most 2-operator terms allowed for AutoMPO conversion to MPO/IQMPO when SVD option is set to false.");
        else if(t.Nops() > 1)
            {
                if(t.first().i == t.last().i)
                    Error("AutoMPO: operators must be placed on distinct sites when SVD option is set to false.");
            }

    //Special SiteTerm objects indicating either
    //a string of identities coming from the first
    //site of the system or the completed Hamitonian
    //for the left-hand side of the system
    SiteTerm IL("IL",0),
             HL("HL",0);

    vector<vector<SiteQN>> basis(N+1);
    for(int n = 0; n < N; ++n)
        basis.at(n).emplace_back(IL,QN());
    for(int n = 1; n <= N; ++n)
        basis.at(n).emplace_back(HL,QN());

    const QN Zero;

    //Fill up the basis array at each site with 
    //the unique operator types occurring on the site
    //and starting a string of operators (i.e. first op of an HTerm)
    for(const auto& ht : terms)
    for(int n = ht.first().i; n < ht.last().i; ++n)
        {
        auto& bn = basis.at(n);
        auto test_has_first = [&ht](const SiteQN& sq){ return sq.st == ht.first(); };
        bool has_first = (find_if(bn.cbegin(),bn.cend(),test_has_first) != bn.end());
        if(!has_first) 
            {
            auto Op = sites.op(ht.first().op,ht.first().i);
            if(checkqn)
                bn.emplace_back(ht.first(),-div(Op));
            else
                bn.emplace_back(ht.first(),Zero);
            }
        }

    if(checkqn)
        {
        auto qn_comp = [&Zero](const SiteQN& sq1,const SiteQN& sq2)
                       {
                       //First two if statements are to artificially make
                       //the default-constructed Zero QN come first in the sort
                       if(sq1.q == Zero && sq2.q != Zero) return true;
                       else if(sq2.q == Zero && sq1.q != Zero) return false;
                       return sq1.q < sq2.q;
                       };
        //Sort bond "basis" elements by quantum number sector:
        for(auto& bn : basis) std::sort(bn.begin(),bn.end(),qn_comp);
        }

    int max_d = 0;
    
    vector<IQIndex> links(N+1);
    vector<IndexQN> inqn;
    for(int n = 0; n <= N; n++)
        {
        auto& bn = basis.at(n);
        inqn.clear();
        QN currq = bn.front().q;
        int currm = 0;
        int count = 0;
        for(auto& sq : bn)
            {
            if(sq.q == currq)
                {
                ++currm;
                }
            else
                {
                inqn.emplace_back(Index(format("hl%d_%d",n,count++),currm),currq);
                currq = sq.q;
                currm = 1;
                }
            }
        inqn.emplace_back(Index(format("hl%d_%d",n,count++),currm),currq);

        links.at(n) = IQIndex(nameint("Hl",n),inqn);
        
        max_d = max(max_d, links.at(n).m());

        //if(n <= 2 or n == N)
        //    {
        //    println("basis for site ",n);
        //    for(size_t l = 0; l < bn.size(); ++l) printfln("%d %s %s",l,bn.at(l).st,bn.at(l).q);
        //    println();
        //    printfln("IQIndex for site %d:\n%s",n,links.at(n));
        //    }
        }
        
        println("Maximal dimension is ", max_d);

#ifdef SHOW_AUTOMPO
    static string ws[100][100];
#endif

    //Create arrays indexed by lattice sites.
    //For lattice site "j", ht_by_n[j] contains
    //all HTerms (operator strings) which begin on,
    //end on, or cross site "j"
    vector<vector<HTerm>> ht_by_n(N+1);
    for(const HTerm& ht : terms) 
    for(const auto& st : ht.ops)
        {
        ht_by_n.at(st.i).push_back(ht);
        }

    for(int n = 1; n <= N; n++)
        {
        auto& bn1 = basis.at(n-1);
        auto& bn  = basis.at(n);

        auto& W = H.Anc(n);
        auto &row = links.at(n-1),
             &col = links.at(n);

        W = IQTensor(dag(sites(n)),prime(sites(n)),dag(row),col);

        for(int r = 0; r < row.m(); ++r)
        for(int c = 0; c < col.m(); ++c)
            {
            auto& rst = bn1.at(r).st;
            auto& cst = bn.at(c).st;

#ifdef SHOW_AUTOMPO
            ws[r][c] = "0";
#endif
            auto rc = IQTensor(dag(row)(r+1)) * IQTensor(col(c+1));

            //Start a new operator string
            if(cst.i == n && rst == IL)
                {
                //Call startTerm to handle fermionic cases with Jordan-Wigner strings
                auto op = startTerm(cst.op);
                W += sites.op(op,n) * rc;
#ifdef SHOW_AUTOMPO
                ws[r][c] = op;
#endif
                }

            //Add identity "string" connecting operator
            //strings of more than two sites in length
            if(cst == rst)
                {
                if(cst.isFermionic())
                    W += sites.op("F",n) * rc;
                else
                    W += sites.op("Id",n) * rc;
#ifdef SHOW_AUTOMPO
                if(cst.isFermionic()) ws[r][c] = "F";
                else                 ws[r][c] = "1";
#endif
                }

            //End operator strings
            if(cst == HL)
                {
                for(const auto& ht : ht_by_n.at(n))
                if(rst == ht.first() && ht.last().i == n)
                    {
                    auto op = endTerm(ht.last().op);
                    W += ht.coef * sites.op(op,n) * rc;
#ifdef SHOW_AUTOMPO
                    ws[r][c] = op;
#endif
                    }
                }

            //Include on-site operators
            if(rst == IL && cst == HL)
                {
                for(const auto& ht : ht_by_n.at(n))
                if(ht.first().i == ht.last().i)
                    {
#ifdef SHOW_AUTOMPO
                    if(isApproxReal(ht.coef))
                        ws[r][c] = format("%.2f %s",ht.coef.real(),ht.first().op);
                    else
                        ws[r][c] = format("%.2f %s",ht.coef,ht.first().op);
#endif
                    W += ht.coef * sites.op(ht.first().op,n) * rc;
                    }
                }

            }

#ifdef SHOW_AUTOMPO
        if(n <= 10 or n == N)
            {
            for(int r = 0; r < row.m(); ++r, println())
            for(int c = 0; c < col.m(); ++c)
                {
                print(ws[r][c],"\t");
                if(ws[r][c].length() < 8 && c == 1) 
                print("\t");
                }
            println("=========================================");
            }
#endif
        }

    H.Anc(1) *= IQTensor(links.at(0)(1));
    H.Anc(N) *= IQTensor(dag(links.at(N))(1));

    //checkQNs(H);

    return H;
    }

template<>
MPO
toMPO<ITensor>(const AutoMPO& a,
               const Args& args)
    {
    static Args checkqn("CheckQNs",false);
    IQMPO res = toMPO<IQTensor>(a,args+checkqn);
    return res.toMPO();
    }
    
AutoMPO::operator MPO() const
    { 
    if(svd_) 
        {
        IQMPO H = ConstructMPOUsingSVD(); 
        return H.toMPO();
        }
    else 
        return toMPO<ITensor>(*this); 
    }

AutoMPO::operator IQMPO() const
    { 
    if(svd_) 
        return ConstructMPOUsingSVD(); 
    else 
        return toMPO<IQTensor>(*this); 
    }


IQMPO
toExpH_ZW1(const AutoMPO& am,
           Complex tau,
           const Args& args)
    {
        
    if(am.usingSVD())
        return am.toExpHUsingSVD_ZW1(tau);
        
    const SiteSet& sites = am.sites();
    IQMPO H(sites);
    const int N = sites.N();

    std::vector<HTerm> terms = am.terms();

    for(auto& t : terms)
    if(t.Nops() > 2) 
        {
        Error("Only at most 2-operator terms allowed for AutoMPO conversion to MPO/IQMPO");
        }

    bool is_complex = fabs(tau.imag()) > fabs(1E-12*tau.real());

    //Special SiteTerm objects indicating either
    //a string of identities coming from the first
    //site of the system or the completed Hamitonian
    //for the left-hand side of the system
    SiteTerm IL("IL",0);

    vector<vector<SiteQN>> basis(N+1);
    for(int n = 0; n <= N; ++n)
        basis.at(n).emplace_back(IL,QN());

    //Fill up the basis array at each site with 
    //the unique operator types occurring on the site
    //and starting a string of operators (i.e. first op of an HTerm)
    for(const auto& ht : terms)
    for(int n = ht.first().i; n < ht.last().i; ++n)
        {
        auto& bn = basis.at(n);
        auto test = [&ht](const SiteQN& sq){ return sq.st == ht.first(); };
        bool has_first = (std::find_if(bn.cbegin(),bn.cend(),test) != bn.end());
        if(!has_first) 
            {
            auto Op = sites.op(ht.first().op,ht.first().i);
            bn.emplace_back(ht.first(),-div(Op));
            }
        }

    const QN Zero;
    auto qn_comp = [&Zero](const SiteQN& sq1,const SiteQN& sq2)
                   {
                   //First two if statements are to artificially make
                   //the default-constructed Zero QN come first in the sort
                   if(sq1.q == Zero && sq2.q != Zero) return true;
                   else if(sq2.q == Zero && sq1.q != Zero) return false;
                   return sq1.q < sq2.q;
                   };
    //Sort bond "basis" elements by quantum number sector:
    for(auto& bn : basis) std::sort(bn.begin(),bn.end(),qn_comp);

    vector<IQIndex> links(N+1);
    vector<IndexQN> inqn;
    for(int n = 0; n <= N; n++)
        {
        auto& bn = basis.at(n);
        inqn.clear();
        QN currq = bn.front().q;
        int currm = 0;
        int count = 0;
        for(auto& sq : bn)
            {
            if(sq.q == currq)
                {
                ++currm;
                }
            else
                {
                inqn.emplace_back(Index(format("hl%d_%d",n,count++),currm),currq);
                currq = sq.q;
                currm = 1;
                }
            }
        inqn.emplace_back(Index(format("hl%d_%d",n,count++),currm),currq);

        links.at(n) = IQIndex(nameint("Hl",n),inqn);

        //if(n <= 2 or n == N)
        //    {
        //    println("basis for site ",n);
        //    for(size_t l = 0; l < bn.size(); ++l) printfln("%d %s %s",l,bn.at(l).st,bn.at(l).q);
        //    println();
        //    printfln("IQIndex for site %d:\n%s",n,links.at(n));
        //    }
        }

#ifdef SHOW_AUTOMPO
    static string ws[100][100];
#endif

    //Create arrays indexed by lattice sites.
    //For lattice site "j", ht_by_n[j] contains
    //all HTerms (operator strings) which begin on,
    //end on, or cross site "j"
    vector<vector<HTerm>> ht_by_n(N+1);
    for(const HTerm& ht : terms) 
    for(const auto& st : ht.ops)
        {
        ht_by_n.at(st.i).push_back(ht);
        }

    for(int n = 1; n <= N; n++)
        {
        auto& bn1 = basis.at(n-1);
        auto& bn  = basis.at(n);

        auto& W = H.Anc(n);
        auto &row = links.at(n-1),
             &col = links.at(n);

        W = IQTensor(dag(sites(n)),prime(sites(n)),dag(row),col);

        for(int r = 0; r < row.m(); ++r)
        for(int c = 0; c < col.m(); ++c)
            {
            auto& rst = bn1.at(r).st;
            auto& cst = bn.at(c).st;

#ifdef SHOW_AUTOMPO
            ws[r][c] = "0";
#endif
            auto rc = IQTensor(dag(row)(r+1)) * IQTensor(col(c+1));

            //Start a new operator string
            if(cst.i == n && rst == IL)
                {
#ifdef SHOW_AUTOMPO
                ws[r][c] = "(-t)*" + cst.op;
#endif
                auto opname = startTerm(cst.op);
                auto op = sites.op(opname,n) * rc;
                if(is_complex) op *= (-tau);
                else           op *= (-tau.real());
                W += op;
                }

            //Add identity "string" connecting operator
            //strings of more than two sites in length
            if(cst == rst)
                {
#ifdef SHOW_AUTOMPO
                if(cst.isFermionic()) plusAppend(ws[r][c],"F");
                else                 plusAppend(ws[r][c],"1");
#endif
                if(cst.isFermionic())
                    W += sites.op("F",n) * rc;
                else
                    W += sites.op("Id",n) * rc;
                }

            //End operator strings
            if(cst == IL)
                {
                for(const auto& ht : ht_by_n.at(n))
                if(rst == ht.first() && ht.last().i == n)
                    {
#ifdef SHOW_AUTOMPO
                if(isApproxReal(ht.coef))
                    ws[r][c] = format("%.2f %s",ht.coef.real(),cst.op);
                else
                    ws[r][c] = format("%.2f %s",ht.coef,cst.op);
#endif
                    W += ht.coef * sites.op(endTerm(ht.last().op),n) * rc;
                    }
                }

            //Include on-site operators
            if(rst == IL && cst == IL)
                {
                for(const auto& ht : ht_by_n.at(n))
                if(ht.first().i == ht.last().i)
                    {
#ifdef SHOW_AUTOMPO
                    if(isApproxReal(ht.coef))
                        plusAppend(ws[r][c],format("(-t*%.2f)*%s",ht.coef.real(),ht.first().op));
                    else
                        plusAppend(ws[r][c],format("(-t*%.2f)*%s",ht.coef,ht.first().op));
#endif
                    auto op = ht.coef * sites.op(ht.first().op,n) * rc;
                    if(is_complex) op *= (-tau);
                    else           op *= (-tau.real());
                    W += op;
                    }
                }

            }

#ifdef SHOW_AUTOMPO
        if(n <= 10 or n == N)
            {
            for(int r = 0; r < row.m(); ++r, println())
            for(int c = 0; c < col.m(); ++c)
                {
                print(ws[r][c],"\t");
                if(ws[r][c].length() < 8 && c == 1) 
                print("\t");
                }
            println("=========================================");
            }
#endif
        }

    H.Anc(1) *= IQTensor(links.at(0)(1));
    H.Anc(N) *= IQTensor(dag(links.at(N))(1));

    //checkQNs(H);

    return H;
    }

template<>
IQMPO
toExpH<IQTensor>(const AutoMPO& a,
         Complex tau,
         const Args& args)
    {
    auto approx = args.getString("Approx","ZW1");
    IQMPO res;
    if(approx == "ZW1")
        {
        res = toExpH_ZW1(a,tau,args);
        }
    else
        {
        Error(format("Unknown approximation Approx=\"%s\"",approx));
        }
    return res;
    }

template<>
MPO
toExpH<ITensor>(const AutoMPO& a,
                Complex tau,
                const Args& args)
    {
    IQMPO res = toExpH<IQTensor>(a,tau,args);
    return res.toMPO();
    }

std::ostream& 
operator<<(std::ostream& s, const SiteTerm& t)
    {
    s << format("%s(%d)",t.op,t.i);
    return s;
    }

std::ostream& 
operator<<(std::ostream& s, const SiteTermProd& ops)
    {
    const char* pfix = "";
    for(const auto& st : ops) 
        {
        s << format("%s%s(%d)",pfix,st.op,st.i);
        pfix = " ";
        }
    return s;
    }
    
std::ostream& 
operator<<(std::ostream& s, const TermSum& tSum)
    {
    const char* pfix = "";     
    for(const Term &t : tSum.sum)
        {        
        if(abs(t.coef-1.0) > 1E-12) 
            if(isReal(t.coef))
                {
                if (t.coef.real() > 0)
                    s << pfix;
                s << format("%s%f ",pfix,t.coef.real());
                }
            else
                format("%s%f ",pfix,t.coef);
        else
            s << pfix;
            
        s << t.ops;
        pfix = "+";
        }
    return s;
    }

std::ostream& 
operator<<(std::ostream& s, const HTerm& t)
    {
    if(abs(t.coef-1.0) > 1E-12) 
        s << (isReal(t.coef) ? format("%f ",t.coef.real()) : format("%f ",t.coef));
    s << t.ops;
    return s;
    }

std::ostream& 
operator<<(std::ostream& s, const AutoMPO& a)
    {
    s << "AutoMPO:\n";
    for(const auto& t : a.terms()) s << t << "\n";
    return s;
    }


}
