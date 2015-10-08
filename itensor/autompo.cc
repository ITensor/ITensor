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

    
void RewriteFermionic(SiteTermProd &prod, bool isleftFermionic)
    {
    if(prod.empty())
        Error("Empty product in RewriteFermionic is not expected.");    
    
    int i = prod.front().i;
    for(const SiteTerm &t : prod)
        if(t.i != i)
            Error("Multi-site product in RewriteFermionic is not expected.");    
            
    bool isSiteFermionic = IsFermionic(prod);
    if(isSiteFermionic)
        {
        for(SiteTerm &t : prod)
            if(t.isFermionic())
                t.op = fermionicTerm(t.op);
        }
    
    if((isleftFermionic && !isSiteFermionic) || (!isleftFermionic && isSiteFermionic))
        prod.emplace_back("F", i);         
    }

SiteTermProd mult(const SiteTermProd &first, const SiteTermProd &second)
    {
    SiteTermProd prod = first;
    prod.insert( prod.end(), second.begin(), second.end() );
    return prod;
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

HTerm& HTerm::
operator*=(Real x)
    {
    coef *= x;
    return *this;
    }

HTerm& HTerm::
operator*=(Complex x)
    {
    coef *= x;
    return *this;
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
    
    // Check if a proportional term already exists
    auto isProportional = [&t](const HTerm &ht) {return ht.proportionalTo(t); };
    auto it = find_if(terms_.begin(), terms_.end(), isProportional);
    if(it == terms_.end())
        terms_.push_back(t); 
    else
        it->coef += t.coef;
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
    
Complex ComplexMatrix::operator() (int i, int j) const
    {
        Complex re, im;
        if(Im.Storage())
            im = Complex(0, Im(i,j));
        if(Re.Storage())
            re = Complex(Re(i,j),0);
        return re+im;
    }
    
// i,j are 1-based    
void ComplexMatrix::insert(int i, int j, Complex val) 
    {
   // TODO: Determine the size in advance instead of resizing every time
    if(i > Re.Nrows() || j > Re.Ncols())
        {
        int newNRows = max(i, Re.Nrows());
        int newNCols = max(j, Re.Ncols());
        Re.Enlarge(newNRows, newNCols);
        if(Im.Storage() || !isReal(val))
            Im.Enlarge(newNRows, newNCols);
        }
    
    Re(i,j) = val.real();
    if(!isReal(val))
        Im(i,j) = val.imag();    
    }   

// Note that n is 1-based site index
void AutoMPO::AddToTempMPO(int n, const MatElement &elem)
    {
        auto el = find(tempMPO_.at(n-1).begin(), tempMPO_.at(n-1).end(), elem);
        if(el == tempMPO_.at(n-1).end())
            tempMPO_.at(n-1).push_back(elem);
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
int AutoMPO::AddToVec(const SiteTermProd &ops, std::vector<SiteTermProd> &vec)
    {   
    auto it = find(vec.begin(), vec.end(), ops);
    if (it != vec.end())
        return it - vec.begin() + 1;
    vec.push_back(ops);
    return vec.size();
    }
    
void AutoMPO::AddToMPO(int n, MatIndex ind, const Index &row, const Index &col, const TermSum &tSum)
    {
        
    for(const Term &term : tSum.sum)
        {
        if(fabs(term.coef) < 1E-12)
            continue;
            
        clock_t dt1 = clock();
        IQTensor op = sites_.op(term.ops.front().op, n);
        for(auto it = term.ops.begin()+1; it != term.ops.end(); it++)
            op = multSiteOps(op, sites_.op(it->op, n));
        dt1 = clock()-dt1;
        dt1_ += dt1;
        
        clock_t dt2 = clock();            
        if(isReal(term.coef))        
            H_.Anc(n) += term.coef.real() * op * row(ind.row) * col(ind.col);       
        else
            H_.Anc(n) += term.coef * op * row(ind.row) * col(ind.col);
        dt2 = clock()-dt2;
        dt2_ += dt2;
        }
    }

void AutoMPO::ConstructMPOUsingSVD()
    {
    if(H_.valid())
        return;
        
    printfln("Constructing MPO using SVD");
        
    const int N = sites_.N();
    
    H_ = MPO(sites_);
    
    tempMPO_.resize(N);
    finalMPO_.resize(N);
    leftPart_.resize(N);
    rightPart_.resize(N);
    Coeff_.resize(N);
    
    clock_t t = clock();
        
    // Construct left & right partials and the ceofficients matrix on each link 
    // as well as the temporary MPO
    for(const HTerm &ht : terms_)
        for(int n = ht.first().i; n <= ht.last().i; n++)
            {
            SiteTermProd left, onsite, right;
            DecomposeTerm(n, ht.ops, left, onsite, right);
#ifdef SHOW_AUTOMPO            
            println(ht, ", n=",n,": ", left, ", ", onsite, ", ", right);
#endif
            int j,k,l;
            
            if(left.empty())
                {
                j=0;
                if(right.empty()) // on site term
                    k = 0;
                else // term starting on site n
                    k = AddToVec(right, rightPart_.at(n));
                }
            else
                {
                if(right.empty()) // term ending on site n
                    {
                    k = 0;
                    j = AddToVec(onsite, rightPart_.at(n-1));
                    }
                else
                    {
                    j = AddToVec(mult(onsite,right), rightPart_.at(n-1));
                    k = AddToVec(right, rightPart_.at(n));
                    }
                l = AddToVec(left, leftPart_.at(n-1));
                Coeff_.at(n-1).insert(l,j, ht.coef);
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
            MatElement elem(MatIndex(j, k), Term(c, onsite));
            AddToTempMPO(n, elem);
            }
    
    t = clock() - t;
    println("It took ", ((float)t)/CLOCKS_PER_SEC, " seconds to construct the temporary MPO");
    
#ifdef SHOW_AUTOMPO
    println("TempMPO Elements:");
    for(int n=1; n<=N; n++)
        {
        for(const MatElement &elem: tempMPO_.at(n-1))
            println(elem.ind.row,',',elem.ind.col,'\t',elem.term.coef,'\t',elem.term.ops);
        println("=========================================");
        }

    println("Left and Right Partials:");
    for(int n=0; n<N; n++)
    {
        println("Left:");
        for(const SiteTermProd &prod : leftPart_.at(n))
            println(prod);
        println("Right:");
        for(const SiteTermProd &prod : rightPart_.at(n))
            println(prod);

         println("=========================================");
    }
    
    println("Left-Right Coefficients:");        
    for(int n=0; n<N; n++)
    {
         println(Coeff_.at(n).Re);
         if(Coeff_.at(n).isComplex())
             println(Coeff_.at(n).Im);
         println("=========================================");
    }
#endif            

    // SVD Coeff matrix on each link and construct the final MPO matrix for each site
    // Note that for the MPO matrix on site n we need the SVD on both the previous link and the following link
    ComplexMatrix V_n, V_npp;
    int d_n = 0, d_npp = 0; 	// d_n = num of non-zero singular values of Coeff_[n]
    int max_d = 0;

    t = clock();
    
    for(int n=1; n<=N; n++)
        {
            
        // for 1st site there are no left partials and Coeff_[n-1] is empty but Coeff_[n] is not
        // for Nth site Coeff_[n] is empty but Coeff_[n-1] is not
        if(n==N || Coeff_.at(n).isEmpty())
            d_npp = 0;
        else
            {
            clock_t svdt = clock();
            
            Vector D;
            if(Coeff_.at(n).isComplex())
                {                    
                ComplexMatrix U;
                SVD(Coeff_.at(n).Re, Coeff_.at(n).Im, U.Re, U.Im, D, V_npp.Re, V_npp.Im);
                }
            else
                {
                Matrix U;                
                SVD(Coeff_.at(n).Re, U, D, V_npp.Re);
                }

            svdt = clock() - svdt;
            println("SVD for site ", n, " took ", ((float)svdt)/CLOCKS_PER_SEC, " seconds");

            Real epsilon = 1E-12;
            auto isApproxZero = [&epsilon](const Real &val){ return fabs(val) < epsilon; };
            auto firstApproxZero = std::find_if(D.begin(), D.end(), isApproxZero);
            d_npp=firstApproxZero - D.begin();
            }
            
        println("Num of elements in temporary MPO for site ", n, " is ", tempMPO_.at(n-1).size());
        
        clock_t sitet = clock();
        
        finalMPO_.at(n-1).resize(d_n+2);
        for(auto &v : finalMPO_.at(n-1))
            v.resize(d_npp+2);
            
        println("Size of MPO matrix for site ", n, " is (", d_n+2, ", ", d_npp+2, ")");
        
        SiteTermProd opId;
        opId.emplace_back("Id", n);
        
        finalMPO_.at(n-1).at(0).at(0) += Term(1, opId);
        finalMPO_.at(n-1).at(1).at(1) += Term(1, opId);
            
        Complex Zero(0,0);
        for(const MatElement &elem: tempMPO_.at(n-1))
            {
            int k = elem.ind.row;
            int l = elem.ind.col;
    
            if(l==0 && k==0)	// on-site terms
                finalMPO_.at(n-1).at(1).at(0) += elem.term;
            else if(k==0)  	// terms starting on site n
                {
                for(int j=1; j<=d_npp; j++)
                    if(V_npp(j,l) != Zero) // 1-based access of matrix elements
                        finalMPO_.at(n-1).at(1).at(1+j) += Term(elem.term.coef*V_npp(j,l), elem.term.ops);
                        
                }
            else if(l==0) 	// terms ending on site n
                {
                for(int i=1; i<=d_n; i++)
                    if(V_n(i,k) != Zero) // 1-based access of matrix elements
                        finalMPO_.at(n-1).at(1+i).at(0) += Term(elem.term.coef*V_n(i,k), elem.term.ops);
                }
            else 
                {
                for(int i=1; i<=d_n; i++)
                    for(int j=1; j<=d_npp; j++) 
                        if( (V_n(i,k) != Zero)  && (V_npp(j,l) != Zero) ) // 1-based access of matrix elements
                            finalMPO_.at(n-1).at(1+i).at(1+j) += Term(elem.term.coef*V_n(i,k)*V_npp(j,l), elem.term.ops);
                }
            }
            sitet = clock() - sitet;
            println("It took ", ((float)sitet)/CLOCKS_PER_SEC, " seconds to construct the final MPO for site ", n);

        // Store SVD computed at this step for next link
        V_n = V_npp;
        d_n = d_npp;        
        max_d = max(max_d, d_n);
        }
        
        t = clock() - t;
        println("It took ", ((float)t)/CLOCKS_PER_SEC, " seconds to construct the final MPO");
        
        println("Maximal dimension of MPO is ", max_d+2);
        
#ifdef SHOW_AUTOMPO
        for(int n=1; n<=N; n++)
            {
            for(unsigned r = 0; r < finalMPO_.at(n-1).size(); ++r, println())
                for(unsigned c = 0; c < finalMPO_.at(n-1).at(r).size(); ++c)
                    print(finalMPO_.at(n-1).at(r).at(c), "\t");
            println("=========================================");
            }
#endif    

        t = clock();
        
        vector<Index> links(N+1);
        links.at(0) = Index(nameint("Hl",0),finalMPO_.at(0).size());

        for(int n=1; n<=N; n++)
            {
            int nr = finalMPO_.at(n-1).size();
            int nc = n == N ? 2 : finalMPO_.at(n).size();
            
            links.at(n) = Index(nameint("Hl",n),nc);
            auto &row = links.at(n-1),
                 &col = links.at(n);

            H_.Anc(n) = ITensor(sites_.si(n),sites_.siP(n),row,col);
                
            clock_t sitet = clock();
            
            for(int r = 1; r <= nr; ++r)
                for(int c = 1; c <= nc; ++c)
                    AddToMPO(n, MatIndex(r,c), row, col, finalMPO_.at(n-1).at(r-1).at(c-1));    

            sitet = clock() - sitet;
            println("It took ", ((float)sitet)/CLOCKS_PER_SEC, " seconds to construct the MPO matrix for site ", n);

        }
        
    t = clock() - t;
    println("It took ", ((float)t)/CLOCKS_PER_SEC, " seconds to construct the MPO matrices");
    
    println("Operators multiplication took ", ((float)dt1_)/CLOCKS_PER_SEC, " seconds");
    println("Operators addition took ", ((float)dt2_)/CLOCKS_PER_SEC, " seconds");

    
    H_.Anc(1) *= ITensor(links.at(0)(2));
    H_.Anc(N) *= ITensor(links.at(N)(1));
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

    for(auto& t : am.terms())
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
    for(const auto& ht : am.terms())
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
    for(const HTerm& ht : am.terms()) 
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


IQMPO
toExpH_ZW1(const AutoMPO& am,
           Complex tau,
           const Args& args)
    {
    const SiteSet& sites = am.sites();
    IQMPO H(sites);
    const int N = sites.N();

    for(auto& t : am.terms())
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
    for(const auto& ht : am.terms())
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
    for(const HTerm& ht : am.terms()) 
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
