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
SiteTerm() : i(-1), coef(0) { }

SiteTerm::
SiteTerm(const std::string& op_,
         int i_,
         Real coef_)
    :
    op(op_),
    i(i_),
    coef(coef_)
    { }

bool SiteTerm::
operator==(const SiteTerm& other) const
    {
    return (op == other.op && i == other.i && abs(coef-other.coef) < 1E-12);
    }

bool SiteTerm::
proportialTo(const SiteTerm& other) const
    {
    return (op == other.op && i == other.i);
    }

bool
isFermionic(const SiteTerm& st)
    {
    if(!st.op.empty() && st.op.front() == 'C') return true;
    return false;
    }

std::string SiteTermProd::
opStr() const
    {
        std::string opstr;
        for(auto t = ops.begin(); t != ops.end()-1; t++)
            opstr += t->op + "*";
        opstr += ops.back().op;
        return opstr;
    }

SiteTermProd SiteTermProd::
operator*(const SiteTermProd &other) const
    {
    SiteTermProd prod(*this);
    for(SiteTerm op : other.ops)
        prod.ops.push_back(op);
    return prod;
    }
    
bool SiteTermProd::operator==(const SiteTermProd &other) const
    {
    if(ops.size() != other.ops.size())
        return false;

    for(size_t n = 0; n < ops.size(); ++n)
        if(ops[n] != other.ops.at(n)) 
            {
            return false;
            }
    return true;
    }
    
void SiteTermSum::
operator+=(const std::pair<Complex, SiteTermProd> &term)
    {
    opSum.push_back(term);
    }    


HTerm::
HTerm() { }

HTerm::
HTerm(const std::string& op1_,
      int i1_,
      Real x_)
    { 
    add(op1_,i1_,x_);
    }

HTerm::
HTerm(const std::string& op1_,
      int i1_,
      const std::string& op2_,
      int i2_,
      Real x_)
    { 
    add(op1_,i1_,x_);
    add(op2_,i2_);
    }

void HTerm::
add(const std::string& op,
    int i,
    Real x)
    {
    //The following ensures operators remain
    //in site order within the vector "ops"
    auto it = opProd.ops.begin();
    while(it != opProd.ops.end() && it->i < i) ++it;
    if(it!= opProd.ops.end() && it->i == i) Error("AutoMPO: operators must be placed on distinct sites");
    opProd.ops.emplace(it,op,i,x);
    }

bool HTerm::
startsOn(int i) const 
    { 
    if(opProd.ops.empty()) Error("No operators in HTerm");
    return first().i == i; 
    }

bool HTerm::
endsOn(int i) const 
    { 
    if(opProd.ops.empty()) Error("No operators in HTerm");
    return last().i == i; 
    }

bool HTerm::
contains(int i) const 
    { 
    if(opProd.ops.empty()) Error("No operators in HTerm");
    return i >= first().i && i <= last().i; 
    }

Complex HTerm::
coef() const
    {
    if(Nops() == 0) return 0;
    Complex c = 1;
    for(const auto& op : opProd.ops) c *= op.coef;
    return c;
    }

HTerm& HTerm::
operator*=(Real x)
    {
    if(Nops() == 0) Error("No operators in HTerm");
    opProd.ops.front().coef *= x;
    return *this;
    }

HTerm& HTerm::
operator*=(Complex x)
    {
    if(Nops() == 0) Error("No operators in HTerm");
    opProd.ops.front().coef *= x;
    return *this;
    }

bool HTerm::
operator==(const HTerm& other) const
    {  
    if(opProd != other.opProd) 
        return false;

    return true;
    }

bool HTerm::
operator!=(const HTerm& other) const
    {
    return !operator==(other);
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


//TODO:
// o Add support for > 2 site operators
// o Add support for long-range (exponentially-decaying type) operator strings
// o Add support for fermionic operator strings

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

//#define SHOW_AUTOMPO

// Note that n is 1-based site index
void AutoMPO::AddToTempMPO(int n, const MatElement &elem)
    {
        auto el = find(tempMPO_.at(n-1).begin(), tempMPO_.at(n-1).end(), elem);
        if(el == tempMPO_.at(n-1).end())
            tempMPO_.at(n-1).push_back(elem);
    }

void AutoMPO::DecomposeTerm(int n, const SiteTermProd &term, 
                            SiteTermProd &left, SiteTermProd &onsite, SiteTermProd &right) const
    {
    auto isOnSiteOrOnTheRight = [&n](const SiteTerm &t) {return t.i >= n;};
    auto startOfOnSite = find_if(term.ops.begin(), term.ops.end(), isOnSiteOrOnTheRight);
    
    auto isOnTheRight = [&n](const SiteTerm &t) {return t.i > n;};
    auto startOfRightPart = find_if(startOfOnSite, term.ops.end(), isOnTheRight);

    left.ops = std::vector<SiteTerm>(term.ops.begin(), startOfOnSite);
    onsite.ops = std::vector<SiteTerm>(startOfOnSite, startOfRightPart);
    right.ops = std::vector<SiteTerm>(startOfRightPart, term.ops.end());
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
    
void AutoMPO::AddToMPO(int n, Complex coeff, MatIndex ind, 
                        const Index &row, const Index &col, const SiteTermProd &prod)
    {
    std::string opstr = prod.opStr();

    if(isReal(coeff))        
        {
        H_.Anc(n) += coeff.real() * sites_.op(opstr, n) * row(ind.first) * col(ind.second);
        
#ifdef SHOW_AUTOMPO
        std::string str = format("%.2f %s",coeff.real(), opstr);
        if(!mpoStr_[ind.first-1][ind.second-1].empty())
            mpoStr_[ind.first-1][ind.second-1] += "+" + str;
        else 
            mpoStr_[ind.first-1][ind.second-1] = str;
#endif
        }        
    else
        // TODO: SHOW_AUTOMPO
        H_.Anc(n) += coeff * sites_.op(opstr, n) * row(ind.first) * col(ind.second);

    }


void AutoMPO::ConstructMPOUsingSVD()
    {
    if(H_.valid())
        return;
        
    printfln("Constructing MPO using SVD");
        
    const int N = sites_.N();
    
    H_ = MPO(sites_);
    
    tempMPO_.resize(N);
    leftPart_.resize(N);
    rightPart_.resize(N);
    Coeff_.resize(N);
        
    // Construct left & right partials and the ceofficients matrix on each link 
    // as well as the temporary MPO
    for(const HTerm &ht : terms_)
        for(int n = ht.first().i; n <= ht.last().i; n++)
            {
            SiteTermProd left, onsite, right;
            DecomposeTerm(n, ht.opProd, left, onsite, right);
            int j,k,l;
            
            if(left.ops.empty())
                {
                j=0;
                if(right.ops.empty()) // on site term
                    k = 0;
                else // term starting on site n
                    k = AddToVec(right, rightPart_.at(n));
                }
            else
                {
                if(right.ops.empty()) // term ending on site n
                    {
                    k = 0;
                    j = AddToVec(onsite, rightPart_.at(n-1));
                    }
                else
                    {
                    j = AddToVec(onsite*right, rightPart_.at(n-1));
                    k = AddToVec(right, rightPart_.at(n));
                    }
                l = AddToVec(left, leftPart_.at(n-1));
                // TODO: Determine the size in advance instead of resizing every time
                if(l > Coeff_.at(n-1).Nrows() || j > Coeff_.at(n-1).Ncols())
                    Coeff_.at(n-1).Enlarge(max(l,Coeff_.at(n-1).Nrows()), max(j, Coeff_.at(n-1).Ncols()));
                // TODO: Complex coef
                Coeff_.at(n-1).el(l-1, j-1) = ht.coef().real();
                }
                
            // Place the coefficient of the HTerm when the term starts
            Complex c = j==0 ? ht.coef() : 1;
            
            if(onsite.ops.empty())
                // TODO: Handle Fermions
                onsite.ops.emplace_back(SiteTerm("Id",n));
            MatElement elem({j, k}, c, onsite);
            AddToTempMPO(n, elem);
            }
            
#ifdef SHOW_AUTOMPO
    for(int n=1; n<=N; n++)
        {
        for(const MatElement &elem: tempMPO_.at(n-1))
            {
            MatIndex ind = std::get<0>(elem);
            println(ind.first,',',ind.second,'\t',std::get<1>(elem),'\t',std::get<2>(elem));
            }
        println("=========================================");
        }


#endif            

    // SVD Coeff matrix on each link and construct the final MPO matrix for each site
    // Note that for the MPO matrix on site n we need the SVD on both the previous link and the following link
    Matrix V_n, V_npp;
    int d_n = 0, d_npp = 0; 	// d_n = num of non-zero singular values of Coeff_[n]

    vector<Index> links(N+1);
    links.at(0) = Index(nameint("Hl",0),d_n+2);

    for(int n=1; n<=N; n++)
        {
            
        // for 1st site there are no left partials and Coeff_[n-1] is empty but Coeff_[n] is not
        // for Nth site Coeff_[n] is empty but Coeff_[n-1] is not
        if(n==N || !Coeff_.at(n).Store())
            d_npp = 0;
        else
            {
            Matrix U;
            Vector D;
            
            SVD(Coeff_.at(n), U, D, V_npp);
            auto firstZero = find(D.begin(), D.end(), 0);
            d_npp=firstZero - D.begin();
            D.ReduceDimension(d_npp);
            }
            
        links.at(n) = Index(nameint("Hl",n),d_npp+2);
        auto &row = links.at(n-1),
             &col = links.at(n);

        H_.Anc(n) = ITensor(sites_.si(n),sites_.siP(n),row,col);

        H_.Anc(n) += sites_.op("Id",n) * row(1) * col(1);
        H_.Anc(n) += sites_.op("Id",n) * row(2) * col(2);

#ifdef SHOW_AUTOMPO
        mpoStr_[0][0] = "1";
        mpoStr_[1][1] = "1";
#endif        

        for(const MatElement &elem: tempMPO_.at(n-1))
            {
            MatIndex ind = std::get<0>(elem);
            int k = ind.first;
            int l = ind.second;
            Complex c = std::get<1>(elem);            
            if(l==0 && k==0)	// on-site terms
                AddToMPO(n, c, {2,1}, row, col, std::get<2>(elem));
            else if(k==0)  	// terms starting on site n
                {
                for(int j=1; j<=d_npp; j++)
                    if(V_npp(j,l) != 0) // 1-based access of matrix elements
                        AddToMPO(n, c*V_npp(j,l), {2,2+j}, row, col, std::get<2>(elem));
                }
            else if(l==0) 	// terms ending on site n
                {
                for(int i=1; i<=d_n; i++)
                    if(V_n(i,k) != 0) // 1-based access of matrix elements
                        AddToMPO(n, c*V_n(i,k), {2+i,1}, row, col,std::get<2>(elem));
                }
            else 
                {
                for(int i=1; i<=d_n; i++)
                    for(int j=1; j<=d_npp; j++) 
                        if( (V_n(i,k) != 0) && (V_npp(j,l) != 0) ) // 1-based access of matrix elements
                            AddToMPO(n, c*V_n(i,k)*V_npp(j,l), {2+i,2+j}, row, col, std::get<2>(elem));
                }
            }

#ifdef SHOW_AUTOMPO
        if(n <= 10 or n == N)
            {
            for(int r = 0; r < d_n+2; ++r, println())
            for(int c = 0; c < d_npp+2; ++c)
                {
                print(mpoStr_[r][c],"\t");
                if(mpoStr_[r][c].length() < 8 && c == 1) 
                print("\t");
                // reset
                mpoStr_[r][c] = "";
                }
            println("=========================================");
            }
#endif    
        // Store SVD computed at this step for next link
        V_n = V_npp;
        d_n = d_npp;        
        }
    
    H_.Anc(1) *= ITensor(links.at(0)(2));
    H_.Anc(N) *= ITensor(links.at(N)(1));
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
        {
        Error("Only at most 2-operator terms allowed for AutoMPO conversion to MPO/IQMPO");
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
    for(const auto& st : ht.opProd.ops)
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
                W += cst.coef * sites.op(op,n) * rc;
#ifdef SHOW_AUTOMPO
                if(isApproxReal(cst.coef))
                    ws[r][c] = format("%.2f %s",cst.coef.real(),op);
                else
                    ws[r][c] = format("%.2f %s",cst.coef,op);
#endif
                }

            //Add identity "string" connecting operator
            //strings of more than two sites in length
            if(cst == rst)
                {
                /*
                int found = 0;
                for(const auto& ht : ht_by_n.at(n))
                    {
                    if(ht.first() == rst &&
                       ht.first().i != n && 
                       ht.last().i  != n)
                        {
                        for(const auto& st : ht.ops)
                        if(st.i == n)
                            {
                            found += 1;
#ifdef SHOW_AUTOMPO
                            ws[r][c] = format("%.2f %s",st.coef,st.op);
#endif
                            W += st.coef * sites.op(st.op,n) * rc;
                            }
                        }
                    }
                //if(found == 0)
                    */

                if(isFermionic(cst))
                    W += sites.op("F",n) * rc;
                else
                    W += sites.op("Id",n) * rc;
#ifdef SHOW_AUTOMPO
                if(isFermionic(cst)) ws[r][c] = "F";
                else                 ws[r][c] = "1";
#endif
                //if(found > 1)
                //    {
                //    println("Warning: found > 1 at site ",n);
                //    PAUSE
                //    }
                }

            //End operator strings
            if(cst == HL)
                {
                for(const auto& ht : ht_by_n.at(n))
                if(rst == ht.first() && ht.last().i == n)
                    {
                    auto op = endTerm(ht.last().op);
                    W += ht.last().coef * sites.op(op,n) * rc;
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
                    if(isApproxReal(ht.first().coef))
                        ws[r][c] = format("%.2f %s",ht.first().coef.real(),ht.first().op);
                    else
                        ws[r][c] = format("%.2f %s",ht.first().coef,ht.first().op);
#endif
                    W += ht.first().coef * sites.op(ht.first().op,n) * rc;
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
    for(const auto& st : ht.opProd.ops)
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
                if(isApproxReal(cst.coef))
                    ws[r][c] = format("(-t*%.2f)*%s",cst.coef.real(),cst.op);
                else
                    ws[r][c] = format("(-t*%.2f)*%s",cst.coef,cst.op);
#endif
                auto opname = startTerm(cst.op);
                auto op = cst.coef * sites.op(opname,n) * rc;
                if(is_complex) op *= (-tau);
                else           op *= (-tau.real());
                W += op;
                }

            //Add identity "string" connecting operator
            //strings of more than two sites in length
            if(cst == rst)
                {
#ifdef SHOW_AUTOMPO
                if(isFermionic(cst)) plusAppend(ws[r][c],"F");
                else                 plusAppend(ws[r][c],"1");
#endif
                if(isFermionic(cst))
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
                    ws[r][c] = ht.last().op;
#endif
                    W += ht.last().coef * sites.op(endTerm(ht.last().op),n) * rc;
                    }
                }

            //Include on-site operators
            if(rst == IL && cst == IL)
                {
                for(const auto& ht : ht_by_n.at(n))
                if(ht.first().i == ht.last().i)
                    {
#ifdef SHOW_AUTOMPO
                    if(isApproxReal(ht.first().coef))
                        plusAppend(ws[r][c],format("(-t*%.2f)*%s",ht.first().coef.real(),ht.first().op));
                    else
                        plusAppend(ws[r][c],format("(-t*%.2f)*%s",ht.first().coef,ht.first().op));
#endif
                    auto op = ht.first().coef * sites.op(ht.first().op,n) * rc;
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
    if(isReal(t.coef))
        s << format("%f * %s(%d)",t.coef.real(),t.op,t.i);
    else
        s << format("%f * %s(%d)",t.coef,t.op,t.i);
    return s;
    }


std::ostream& 
operator<<(std::ostream& s, const SiteTermProd& prod)
    {
    const char* pfix = "";
    for(const auto& st : prod.ops) 
        {
        s << format("%s%s(%d)",pfix,st.op,st.i);
        pfix = " ";
        }
    return s;
    }
    
std::ostream& 
operator<<(std::ostream& s, const SiteTermSum& sum)
    {
    const char* pfix = "";     
    for(const std::pair<Complex, SiteTermProd> &t : sum.opSum)
        {        
        if(abs(t.first-1.0) > 1E-12) 
            if(isReal(t.first))
                {
                if (t.first.real() > 0)
                    s << pfix;
                s << format("%s%f ",pfix,t.first.real());
                }
            else
                format("%s%f ",pfix,t.first);
        else
            s << pfix;
            
        s << t.second;
        pfix = "+";
        }
    return s;
    }

std::ostream& 
operator<<(std::ostream& s, const HTerm& t)
    {
    if(abs(t.coef()-1.0) > 1E-12) 
        s << (isReal(t.coef()) ? format("%f ",t.coef().real()) : format("%f ",t.coef()));
    s << t.opProd;
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
