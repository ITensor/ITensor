//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "autompo.h"

using std::find;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace itensor {

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

If F == fermiPhase, i.e. F = (-1)^(# of fermions)
Then we make c and cdagger both have F's going off to the left,
towards site.   So in c(1) cdag(2),  c(1) has an F in front of
it, which evaluates to -1 since c requres
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

IQMPO
toIQMPO(const AutoMPO& am,
        const OptSet& opts)
    {
    const SiteSet& sites = am.sites();
    IQMPO H(sites);
    const int N = sites.N();

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

        if(n <= 2 or n == N)
            {
            println("basis for site ",n);
            for(size_t l = 0; l < bn.size(); ++l) printfln("%d %s %s",l,bn.at(l).st,bn.at(l).q);
            println();
            printfln("IQIndex for site %d:\n%s",n,links.at(n));
            }
        }

    static string ws[100][100];

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

            ws[r][c] = "0";
            auto rc = IQTensor(dag(row)(r+1)) * IQTensor(col(c+1));

            //Start a new operator string
            if(cst.i == n && rst == IL)
                {
                ws[r][c] = format("%.2f %s",cst.coef,cst.op);
                W += cst.coef * sites.op(cst.op,n) * rc;
                }

            //Add identity "string" connecting operator
            //strings of more than two sites in length
            if(cst == rst)
                {
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
                            ws[r][c] = format("%.2f %s",st.coef,st.op);
                            W += st.coef * sites.op(st.op,n) * rc;
                            }
                        }
                    }

                if(found == 0)
                    {
                    ws[r][c] = "1";
                    W += sites.op("Id",n) * rc;
                    }

                if(found > 1)
                    {
                    println("Warning: found > 1 at site ",n);
                    PAUSE
                    }
                }

            //End operator strings
            if(cst == HL)
                {
                for(const auto& ht : ht_by_n.at(n))
                if(rst == ht.first() && ht.last().i == n)
                    {
                    ws[r][c] = ht.last().op;
                    W += ht.last().coef * sites.op(ht.last().op,n) * rc;
                    }
                }

            //Include on-site operators
            if(rst == IL && cst == HL)
                {
                for(const auto& ht : ht_by_n.at(n))
                if(ht.first().i == ht.last().i)
                    {
                    ws[r][c] = format("%.2f %s",ht.first().coef,ht.first().op);
                    W += ht.first().coef * sites.op(ht.first().op,n) * rc;
                    }
                }

            }

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
        }

    H.Anc(1) *= IQTensor(links.at(0)(1));
    H.Anc(N) *= IQTensor(dag(links.at(N))(1));

    //checkQNs(H);

    return H;
    }

MPO
toMPO(const AutoMPO& a,
      const OptSet& opts)
    {
    IQMPO res = toIQMPO(a,opts);
    return res.toMPO();
    }


IQMPO
toIQExpH_ZW1(const AutoMPO& am,
             Complex tau,
             const OptSet& opts)
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

        if(n <= 2 or n == N)
            {
            println("basis for site ",n);
            for(size_t l = 0; l < bn.size(); ++l) printfln("%d %s %s",l,bn.at(l).st,bn.at(l).q);
            println();
            printfln("IQIndex for site %d:\n%s",n,links.at(n));
            }
        }

    static string ws[100][100];

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

            ws[r][c] = "0";
            auto rc = IQTensor(dag(row)(r+1)) * IQTensor(col(c+1));

            //Start a new operator string
            if(cst.i == n && rst == IL)
                {
                ws[r][c] = format("(-t*%.2f)*%s",cst.coef,cst.op);
                auto op = cst.coef * sites.op(cst.op,n) * rc;
                if(is_complex) op *= (-tau);
                else           op *= (-tau.real());
                W += op;
                }

            //Add identity "string" connecting operator
            //strings of more than two sites in length
            if(cst == rst)
                {
                plusAppend(ws[r][c],"1");
                W += sites.op("Id",n) * rc;
                }

            //End operator strings
            if(cst == IL)
                {
                for(const auto& ht : ht_by_n.at(n))
                if(rst == ht.first() && ht.last().i == n)
                    {
                    ws[r][c] = ht.last().op;
                    W += ht.last().coef * sites.op(ht.last().op,n) * rc;
                    }
                }

            //Include on-site operators
            if(rst == IL && cst == IL)
                {
                for(const auto& ht : ht_by_n.at(n))
                if(ht.first().i == ht.last().i)
                    {
                    plusAppend(ws[r][c],format("(-t*%.2f)*%s",ht.first().coef,ht.first().op));
                    auto op = ht.first().coef * sites.op(ht.first().op,n) * rc;
                    if(is_complex) op *= (-tau);
                    else           op *= (-tau.real());
                    W += op;
                    }
                }

            }

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
        }

    H.Anc(1) *= IQTensor(links.at(0)(1));
    H.Anc(N) *= IQTensor(dag(links.at(N))(1));

    //checkQNs(H);

    return H;
    }

IQMPO
toIQExpH(const AutoMPO& a,
         Complex tau,
         const OptSet& opts)
    {
    auto approx = opts.getString("Approx","ZW1");
    IQMPO res;
    if(approx == "ZW1")
        {
        res = toIQExpH_ZW1(a,tau,opts);
        }
    else
        {
        Error(format("Unknown approximation Approx=\"%s\"",approx));
        }
    return res;
    }

MPO
toExpH(const AutoMPO& a,
       Complex tau,
       const OptSet& opts)
    {
    IQMPO res = toIQExpH(a,tau,opts);
    return res.toMPO();
    }

std::ostream& 
operator<<(std::ostream& s, const SiteTerm& t)
    {
    s << format("%f * %s(%d)",t.coef,t.op,t.i);
    return s;
    }


std::ostream& 
operator<<(std::ostream& s, const HTerm& t)
    {
    const char* pfix = "";
    if(fabs(t.coef()-1) > 1E-12) s << format("%f ",t.coef());
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


};
