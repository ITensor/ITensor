//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include <algorithm>
#include "itensor/util/print_macro.h"
#include "itensor/mps/autompo.h"

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
isReal(const Cplx& z) { return z.imag() == 0; }

bool
isApproxReal(const Cplx& z, Real epsilon = 1E-12) { return std::fabs(z.imag()) < epsilon; }

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
    ops.emplace_back(op,i,x);
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

Complex HTerm::
coef() const
    {
    if(Nops() == 0) return 0;
    Complex c = 1;
    for(const auto& op : ops) c *= op.coef;
    return c;
    }

HTerm& HTerm::
operator*=(Real x)
    {
    if(Nops() == 0) Error("No operators in HTerm");
    ops.front().coef *= x;
    return *this;
    }

HTerm& HTerm::
operator*=(Complex x)
    {
    if(Nops() == 0) Error("No operators in HTerm");
    ops.front().coef *= x;
    return *this;
    }

bool HTerm::
operator==(const HTerm& other) const
    {
    if(Nops() != other.Nops()) return false;

    for(size_t n = 0; n <= ops.size(); ++n)
    if(ops[n] != other.ops.at(n)) 
        {
        return false;
        }

    return true;
    }

bool HTerm::
operator!=(const HTerm& other) const
    {
    return !operator==(other);
    }

void
sort(HTerm & ht)
    {
    if(ht.ops.size() <= 1) return;

    auto op = [&ht](size_t n)->SiteTerm& { return ht.ops.at(n); };

    //print("Before sorting, ops are:");
    //for(auto n : range(ht.ops.size()))
    //    {
    //    printf(" %s*%s_%d",op(n).coef,op(n).op,op(n).i);
    //    }
    //println();

    //Do bubble sort: O(n^2) but allows making 
    //pair-wise comparison for fermion signs
    bool did_swap = true;
    while(did_swap)
        {
        did_swap = false;
        for(auto n : range(ht.ops.size()-1))
            {
            if(op(n).i == op(n+1).i) 
                {
                Error("AutoMPO: cannot put two operators on same site in a single term");
                }
            if(op(n).i > op(n+1).i) 
                {
                std::swap(op(n),op(n+1));
                did_swap = true;
                if(isFermionic(op(n)) && isFermionic(op(n+1)))
                    {
                    op(n+1).coef *= -1;
                    }
                }
            }
        }
    //print("After sorting, ops are:");
    //for(auto n : range(ht.ops.size()))
    //    {
    //    printf(" %s*%s_%d",op(n).coef,op(n).op,op(n).i);
    //    }
    //println();
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

template<typename Tensor>
MPOt<Tensor>
toMPOImpl(AutoMPO const& am,
          Args const& args)
    {
    using IndexT = typename Tensor::index_type;
    auto checkqn = args.getBool("CheckQN",true);

    SiteSet const& sites = am.sites();
    auto H = MPOt<Tensor>(sites);
    auto N = sites.N();

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
    //(unique including their coefficient)
    //and starting a string of operators (i.e. first op of an HTerm)
    for(auto& ht : am.terms())
    for(auto n : range(ht.first().i,ht.last().i))
        {
        auto& bn = basis.at(n);
        auto test_has_first = [&ht](SiteQN const& sq){ return sq.st == ht.first(); };
        bool has_first = (stdx::find_if(bn,test_has_first) != bn.end());
        if(!has_first) 
            {
            auto Op = sites.op(ht.first().op,ht.first().i);
            if(checkqn)
                {
                bn.emplace_back(ht.first(),-div(Op));
                }
            else
                {
                bn.emplace_back(ht.first(),Zero);
                }
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

    auto links = vector<IndexT>(N+1);
    vector<IndexQN> inqn;
    for(int n = 0; n <= N; ++n)
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

        links.at(n) = IQIndex(nameint("Hl",n),std::move(inqn));

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
    for(auto& ht : am.terms()) 
    for(auto& st : ht.ops)
        {
        ht_by_n.at(st.i).push_back(ht);
        }

    for(auto n : range1(N))
        {
        auto& bn1 = basis.at(n-1);
        auto& bn  = basis.at(n);

        auto& W = H.Anc(n);
        auto &row = links.at(n-1),
             &col = links.at(n);

        W = Tensor(dag(sites(n)),prime(sites(n)),dag(row),col);

        for(auto r : range(row.m()))
        for(auto c : range(col.m()))
            {
            auto& rst = bn1.at(r).st;
            auto& cst = bn.at(c).st;


#ifdef SHOW_AUTOMPO
            ws[r][c] = "0";
#endif
            auto rc = setElt(dag(row)(r+1)) * setElt(col(c+1));

            //Start a new operator string
            if(cst.i == n && rst == IL)
                {
                //Call startTerm to handle fermionic cases with Jordan-Wigner strings
                auto op = startTerm(cst.op);
                //if(Global::debug1())
                //    {
                //    println("\nAttempting to add the following");
                //    PrintData(sites.op(op,n));
                //    printfln("cst.coef = %f",cst.coef);
                //    PrintData(cst.coef * sites.op(op,n));
                //    auto tmp = cst.coef * sites.op(op,n) * rc;
                //    PrintData(tmp);
                //    PrintData(W);
                //    EXIT
                //    }
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
                //Check if operator is an ending operator
                for(const auto& ht : ht_by_n.at(n))
                if(rst == ht.first() && ht.last().i == n)
                    {
                    auto op = endTerm(ht.last().op);
                    W += ht.last().coef * sites.op(op,n) * rc;
#ifdef SHOW_AUTOMPO
                    ws[r][c] = op;
                    auto coef = ht.last().coef;
                    if(isApproxReal(coef))
                        ws[r][c] = format("%.2f %s",coef.real(),op);
                    else
                        ws[r][c] = format("%.2f %s",coef,op);
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

    H.Anc(1) *= setElt(links.at(0)(1));
    H.Anc(N) *= setElt(dag(links.at(N))(1));

    //checkQNs(H);

    return H;
    }

template<>
MPO 
toMPO(AutoMPO const& am, Args const& args) 
    { 
    return toMPOImpl<ITensor>(am,{args,"CheckQN",false});
    }
template<>
IQMPO 
toMPO(AutoMPO const& am, Args const& args) 
    { 
    return toMPOImpl<IQTensor>(am,args);
    }

//template<>
//MPO
//toMPO<ITensor>(const AutoMPO& a,
//               const Args& args)
//    {
//    auto checkqn = Args("CheckQNs",false);
//    auto res = toMPO<IQTensor>(a,args+checkqn);
//    return res.toMPO();
//    }


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

    bool is_complex = std::fabs(tau.imag()) > std::fabs(1E-12*tau.real());

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

        links.at(n) = IQIndex(nameint("Hl",n),std::move(inqn));

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
            auto rc = setElt(dag(row)(r+1)) * setElt(col(c+1));

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

    H.Anc(1) *= setElt(links.at(0)(1));
    H.Anc(N) *= setElt(dag(links.at(N))(1));

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
operator<<(std::ostream& s, const HTerm& t)
    {
    const char* pfix = "";
    if(abs(t.coef()-1.0) > 1E-12) 
        s << (isReal(t.coef()) ? format("%f ",t.coef().real()) : format("%f ",t.coef()));
    for(const auto& st : t.ops) 
        {
        s << format("%s%s(%d)",pfix,st.op,st.i);
        pfix = " ";
        }
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
