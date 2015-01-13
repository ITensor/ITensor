#define THIS_IS_MAIN
#include "core.h"
#include "svdalgs.h"
#include "matrix.h"
#include "indent.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>

#define rNum()    ( (double)rand() / (1.0+(double)RAND_MAX) )

using namespace std;
using namespace itensor;

const Real pi = atan2(0,-1);
const Real err = 10e-10;
const int maxm = 4;
const Real printcutoff = 1E-3;

class myMPS 
    {
    public:

    void
    init(const vector<Index>);

    int
    N();

    int
    C();
    
    ITensor
    D();

    Index
    S(int);

    Index
    L(int);

    ITensor
    V(int);

    ITensor
    getPsi();
    
    Real 
    getNorm();
    
    void
    shiftCenter(int);

    private:

    int n_; // Number of sites
    int c_; // Center
    vector<ITensor> V_; // Tensors
    vector<Index> s_;
    vector<Index> l_;
    ITensor D_;

    };

void
myMPS::init(const vector<Index> s)
    {
    n_ = s.size()-1;
    s_ = s;
    c_ = 1;

    vector<Index> l(n_+1);
    for(int i = 1; i <= n_; i++)
        l[i] = Index(nameint("l",i),1,Link);
    l_ = l;
    
    D_ = ITensor(l_[1],l_[2],1.0);

    vector<ITensor> V(n_+1);
    V[1] = ITensor(s_[1],l_[1]);
    V[1](s_[1](1),l_[1](1)) = 1.0;
    for(int i = 2; i < n_; i++)
        {
        V[i] = ITensor(s_[i],l_[i],l_[i+1]);
        V[i](s_[i](1),l_[i](1),l_[i+1](1)) = 1.0;
        }
    V[n_] = ITensor(s_[n_],l_[n_]);
    V[n_](s[n_](1),l_[n_](1)) = 1.0;
    V_ = V;
    
    }

int
myMPS::N()
    {
    return n_;
    }

int
myMPS::C()
    {
    return c_;
    }

ITensor
myMPS::D()
    {
    return D_;
    }

Index
myMPS::S(int i)
    {
    return s_[i];
    }

Index
myMPS::L(int i)
    {
    return l_[i];
    }

ITensor
myMPS::V(int i)
    {
    return V_[i];
    }

ITensor
myMPS::getPsi()
    {
    ITensor Psi = V_[1];
    for(int i = 2; i <= c_; i++)
        Psi = Psi*V_[i];
    Psi = Psi*D_;
    for(int i = c_+1; i <= n_; i++)
        Psi = Psi*V_[i];
    return Psi;
    }

Real
myMPS::getNorm()
    {
    ITensor Norm = V_[1]*prime(conj(V_[1]),Link);
    for(int i = 2; i <= c_; i++)
        Norm = Norm*V_[i]*prime(conj(V_[i]),Link);
    Norm = Norm*D_;
    for(int i = c_+1; i <= n_; i++)
        Norm = Norm*V_[i]*prime(conj(V_[i]),Link);
    return Norm.toReal();
    }

void
myMPS::shiftCenter(int shift)
    {
    ITensor psi;
    OptSet opts;
    opts.add("Cutoff",err);
    opts.add("Maxm",maxm);

    if(shift > 0)
        {
        for(int k = 1; k <= shift; k++)
            {
            c_++;
            V_[c_] = D_*V_[c_];
            psi = V_[c_]*V_[c_+1]; 
            svd(psi,V_[c_],D_,V_[c_+1],opts);
            l_[c_] = commonIndex(V_[c_],D_);
            l_[c_+1] = commonIndex(V_[c_+1],D_);
            }
        }
    else if(shift < 0)
        {
        for(int k = 1; k <= -shift; k++)
            {
            c_--;
            V_[c_+1] = V_[c_+1]*D_;
            psi = V_[c_]*V_[c_+1];
            svd(psi,V_[c_],D_,V_[c_+1],opts);
            l_[c_] = commonIndex(V_[c_],D_);
            l_[c_+1] = commonIndex(V_[c_+1],D_);
            }
        }
    }

typedef pair<Real,int> mypair;
bool comparator ( const mypair& l, const mypair& r)
    { return l.first > r.first; }

ITensor
rot(Index s1, Index s2, Real theta, string type)
    {
    Index sP1 = prime(s1);
    Index sP2 = prime(s2);
    ITensor rres(s1,s2,sP1,sP2), ires(s1,s2,sP1,sP2);

    if(type == "14")
        {
        rres(s1(1),s2(1),sP1(1),sP2(1)) = cos(theta/2);
        rres(s1(1),s2(2),sP1(1),sP2(2)) = cos(theta/2);
        rres(s1(2),s2(1),sP1(2),sP2(1)) = cos(theta/2);
        rres(s1(2),s2(2),sP1(2),sP2(2)) = cos(theta/2);
        ires(s1(2),s2(2),sP1(1),sP2(1)) = sin(theta/2);
        ires(s1(1),s2(2),sP1(2),sP2(1)) = -sin(theta/2);
        ires(s1(2),s2(1),sP1(1),sP2(2)) = -sin(theta/2);
        ires(s1(1),s2(1),sP1(2),sP2(2)) = sin(theta/2);
        }
    if(type == "13")
        {
        rres(s1(1),s2(1),sP1(1),sP2(1)) = cos(theta/2);
        rres(s1(1),s2(2),sP1(1),sP2(2)) = cos(theta/2);
        rres(s1(2),s2(1),sP1(2),sP2(1)) = cos(theta/2);
        rres(s1(2),s2(2),sP1(2),sP2(2)) = cos(theta/2);
        rres(s1(2),s2(2),sP1(1),sP2(1)) = sin(theta/2);
        rres(s1(1),s2(2),sP1(2),sP2(1)) = -sin(theta/2);
        rres(s1(2),s2(1),sP1(1),sP2(2)) = sin(theta/2);
        rres(s1(1),s2(1),sP1(2),sP2(2)) = -sin(theta/2);
        }
    if(type == "12")
        {
        rres(s1(1),s2(1),sP1(1),sP2(1)) = cos(theta/2);
        rres(s1(1),s2(2),sP1(1),sP2(2)) = cos(theta/2);
        rres(s1(2),s2(1),sP1(2),sP2(1)) = cos(theta/2);
        rres(s1(2),s2(2),sP1(2),sP2(2)) = cos(theta/2);
        ires(s1(1),s2(1),sP1(1),sP2(1)) = sin(theta/2);
        ires(s1(1),s2(2),sP1(1),sP2(2)) = sin(theta/2);
        ires(s1(2),s2(1),sP1(2),sP2(1)) = -sin(theta/2);
        ires(s1(2),s2(2),sP1(2),sP2(2)) = -sin(theta/2);
        }
    if(type == "24")
        {
        rres(s1(1),s2(1),sP1(1),sP2(1)) = cos(theta/2);
        rres(s1(1),s2(2),sP1(1),sP2(2)) = cos(theta/2);
        rres(s1(2),s2(1),sP1(2),sP2(1)) = cos(theta/2);
        rres(s1(2),s2(2),sP1(2),sP2(2)) = cos(theta/2);
        rres(s1(2),s2(2),sP1(1),sP2(1)) = -sin(theta/2);
        rres(s1(1),s2(2),sP1(2),sP2(1)) = -sin(theta/2);
        rres(s1(2),s2(1),sP1(1),sP2(2)) = sin(theta/2);
        rres(s1(1),s2(1),sP1(2),sP2(2)) = sin(theta/2);
        }
    if(type == "23")
        {
        rres(s1(1),s2(1),sP1(1),sP2(1)) = cos(theta/2);
        rres(s1(1),s2(2),sP1(1),sP2(2)) = cos(theta/2);
        rres(s1(2),s2(1),sP1(2),sP2(1)) = cos(theta/2);
        rres(s1(2),s2(2),sP1(2),sP2(2)) = cos(theta/2);
        ires(s1(2),s2(2),sP1(1),sP2(1)) = sin(theta/2);
        ires(s1(1),s2(2),sP1(2),sP2(1)) = sin(theta/2);
        ires(s1(2),s2(1),sP1(1),sP2(2)) = sin(theta/2);
        ires(s1(1),s2(1),sP1(2),sP2(2)) = sin(theta/2);
        }

    return rres + Complex_i*ires;
    }

//
// Block diagonalize a real anti-symmetric matrix A
// Order the eigenvalues from largest to smallest
//
void
BlockDiagonalize(const MatrixRef& A, Matrix& Atilde, Matrix& W)
    {
    int n = A.Nrows()/2;
    SchurDecomp(A, Atilde, W);

    // Sort Schur vectors
    for(int j = 1; j <= n; j++)
        {
        if(Atilde(2*j-1,2*j) < 0)
            {
            Vector Vtemp = W.Row(2*j-1);
            W.Row(2*j-1) = W.Row(2*j);
            W.Row(2*j) = Vtemp;
            }
        }

    Atilde = W*A*W.t();

    vector<mypair> eigs(n);
    vector<int> order(n+1);
    for(int j = 1; j <= n; j++)
        {
        eigs[j-1].first = Atilde(2*j-1,2*j);
        eigs[j-1].second = j;
        }

    sort(eigs.begin(), eigs.end(), comparator);

    for(int j = 1; j <= n; j++)
        {
        order[j] = eigs[j-1].second;
        }

    Matrix Wtemp = W;

    for(int j = 1; j <= n; j++)
        {
        W.Row(2*j-1) = Wtemp.Row(2*order[j]-1);
        W.Row(2*j) = Wtemp.Row(2*order[j]);
        }

    // Make Schur vectors columns of W
    W = W.t();

    Atilde = W.t()*A*W;
    }

Real
GetGate(Matrix& Gamma, Matrix& WB, int col, int j1, int j2)
    {
    int b = WB.Ncols()/2;
    int n = Gamma.Ncols()/2;
    Matrix rot(2,2);
    Real theta;
    Vector WBj,
           Gammaj;

    // Get \theta_{j1,j2}
    // cout << "Get \\theta_{" << j1 << "," << j2 << "}" << endl;
    theta = atan2(-WB(j2,col),WB(j1,col));
    //cout << "\\theta_{" << j1 << "," << j2 << "} = " << theta << endl;
    rot(1,1) = rot(2,2) = cos(theta);
    rot(1,2) = sin(theta);
    rot(2,1) = -sin(theta);
    // Swap j1+1, j2
    WBj = WB.Row(j1+1);
    WB.Row(j1+1) = WB.Row(j2);
    WB.Row(j2) = WBj;
    Gammaj = Gamma.Column(j1+1);
    Gamma.Column(j1+1) = Gamma.Column(j2);
    Gamma.Column(j2) = Gammaj;
    Gammaj = Gamma.Row(j1+1);
    Gamma.Row(j1+1) = Gamma.Row(j2);
    Gamma.Row(j2) = Gammaj;
    // Rotate j1, j1+1
    WB.SubMatrix(j1,j1+1,1,2*b) = rot.t()*WB.SubMatrix(j1,j1+1,1,2*b);
    Gamma.SubMatrix(j1,j1+1,1,2*n) = rot.t()*Gamma.SubMatrix(j1,j1+1,1,2*n);
    Gamma.SubMatrix(1,2*n,j1,j1+1) = Gamma.SubMatrix(1,2*n,j1,j1+1)*rot;
    // Swap back j1+1, j2
    WBj = WB.Row(j1+1);
    WB.Row(j1+1) = WB.Row(j2);
    WB.Row(j2) = WBj;
    Gammaj = Gamma.Column(j1+1);
    Gamma.Column(j1+1) = Gamma.Column(j2);
    Gamma.Column(j2) = Gammaj;
    Gammaj = Gamma.Row(j1+1);
    Gamma.Row(j1+1) = Gamma.Row(j2);
    Gamma.Row(j2) = Gammaj;

    //cout << "W_B = " << endl;
    //cout << WB;
    
    return theta;
    }

void
ApplyGate(Matrix& A, Real theta, int j1, int j2)
    {
    int n = A.Ncols()/2;
    Matrix rot(2,2);
    Vector Aj;

    rot(1,1) = rot(2,2) = cos(theta);
    rot(1,2) = sin(theta);
    rot(2,1) = -sin(theta);
    // Swap j1+1, j2
    Aj = A.Column(j1+1);
    A.Column(j1+1) = A.Column(j2);
    A.Column(j2) = Aj;
    Aj = A.Row(j1+1);
    A.Row(j1+1) = A.Row(j2);
    A.Row(j2) = Aj;
    // Rotate j1, j1+1
    A.SubMatrix(j1,j1+1,1,2*n) = rot.t()*A.SubMatrix(j1,j1+1,1,2*n);
    A.SubMatrix(1,2*n,j1,j1+1) = A.SubMatrix(1,2*n,j1,j1+1)*rot;
    // Swap back j1+1, j2
    Aj = A.Column(j1+1);
    A.Column(j1+1) = A.Column(j2);
    A.Column(j2) = Aj;
    Aj = A.Row(j1+1);
    A.Row(j1+1) = A.Row(j2);
    A.Row(j2) = Aj;

    }
ITensor
op(Index s, string name)
    {
    Index sP = prime(s);

    ITensor res(s,sP);

    if(name == "Id")
        {
        res(s(1),sP(1)) = 1.0;
        res(s(2),sP(2)) = 1.0;
        }
    else if(name == "a")
        res(s(1),sP(2)) = 1.0;
    else if(name == "aT")
        res(s(2),sP(1)) = 1.0;
    else if(name == "n")
        res(s(2),sP(2)) = 1.0;

    return res;
    }

void
GMPS(Matrix phi, int N_up, int N_par, vector<ITensor>& VB, int&c, ITensor& D, vector<Index>& s, vector<Index>& l, int b)
    {
    int N_sites = phi.Nrows();
    int N_dn = N_par - N_up;
    int n = N_sites;

    // Form Gamma matrix
    Matrix phi_up;
    Matrix phi_dn;
    Matrix Lambda_up(N_sites,N_sites);
    Matrix Lambda_dn(N_sites,N_sites);
    Lambda_up = 0.0, Lambda_dn = 0.0;

    if(N_up > 0)
        {
        phi_up = phi.SubMatrix(1,N_sites,1,N_up);
        Lambda_up = phi_up*phi_up.t();
        }
    if(N_dn > 0)
        {
        phi_dn = phi.SubMatrix(1,N_sites,N_up+1,N_par);
        Lambda_dn = phi_dn*phi_dn.t();
        }

    Matrix Gamma(2*n,2*n);
    Gamma = 0.0;
    Matrix Lambdare(n,n);
    Lambdare = Lambda_up;
    Matrix Lambdaim(n,n);
    Lambdaim = 0.0;

    Matrix BlockId(2*n,2*n);
    BlockId = 0;

    for(int j = 1; j <= n; j++)
        BlockId(2*j-1,2*j) = 1.0;
    BlockId -= BlockId.t();

    Gamma.SubMatrix(1,n,1,n) = Lambdaim;
    Gamma.SubMatrix(1,n,n+1,2*n) = -Lambdare;
    Gamma.SubMatrix(n+1,2*n,1,n) = Lambdare;
    Gamma.SubMatrix(n+1,2*n,n+1,2*n) = Lambdaim;

    Matrix Id(2*n,2*n);
    Id = 1.0;
    Matrix P(2*n,2*n);
    P = 0.0;
    for(int j = 1; j <= n; j++)
        P.Row(j) = Id.Row(2*j-1);
    for(int j = 1; j <= n; j++)
        P.Row(n+j) = Id.Row(2*j);

    Gamma = P.t() * Gamma * P;

    Gamma = BlockId+2.0*Gamma;

    cout << Gamma;

    // Output file
    ofstream fangles, fpositions;
    fangles.open("angles.txt");
    fpositions.open("positions.txt");

    int ngates = (b-1)*(2*n-b)/2;
    Vector angles(5*ngates);
    angles = 0.0;
    int jangle = 0;

    vector<int> begpos(5*ngates+1);
    vector<int> endpos(5*ngates+1);

    Vector nu(n);
    nu = 0.0;
    int site = 1;

    Matrix W;
    Matrix Gammatilde;
    BlockDiagonalize(Gamma,Gammatilde,W);

    Matrix Gammaorig = Gamma;
    Matrix Gammatildeorig = Gammatilde;

    Matrix Gammaappr = Gamma;

    int ntemp = n;
    int btemp = b;
    while(ntemp >= btemp)
        {
        Matrix GammaB = Gamma.SubMatrix(1,2*btemp,1,2*btemp);
        Matrix GammaBtilde(2*btemp,2*btemp);
        Matrix WB(2*btemp,2*btemp);
    
        BlockDiagonalize(GammaB,GammaBtilde,WB);

        for(int j = btemp-1; j >= 1; j--)
            {
            int start = 2*j-1;
            Real theta14 = GetGate(Gamma,WB,1,start,start+3);
            ApplyGate(Gammaappr,theta14,2*site-2+start,2*site-2+start+3);
            fpositions << 2*site-2+start << " " << 2*site-2+start+3 << endl;
            fangles << theta14 << endl;
            jangle++;
            angles(jangle) = theta14;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+3;

            Real theta13 = GetGate(Gamma,WB,1,start,start+2);
            ApplyGate(Gammaappr,theta13,2*site-2+start,2*site-2+start+2);
            fpositions << 2*site-2+start << " " << 2*site-2+start+2 << endl;
            fangles << theta13 << endl;
            jangle++;
            angles(jangle) = theta13;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+2;

            Real theta12 = GetGate(Gamma,WB,1,start,start+1);
            ApplyGate(Gammaappr,theta12,2*site-2+start,2*site-2+start+1);
            fpositions << 2*site-2+start << " " << 2*site-2+start+1 << endl;
            fangles << theta12 << endl;
            jangle++;
            angles(jangle) = theta12;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+1;

            Real theta24 = GetGate(Gamma,WB,2,start+1,start+3);
            ApplyGate(Gammaappr,theta24,2*site-2+start+1,2*site-2+start+3);
            fpositions << 2*site-2+start+1 << " " << 2*site-2+start+3 << endl;
            fangles << theta24 << endl;
            jangle++;
            angles(jangle) = theta24;
            begpos[jangle] = 2*site-2+start+1;
            endpos[jangle] = 2*site-2+start+3;

            Real theta23 = GetGate(Gamma,WB,2,start+1,start+2);
            ApplyGate(Gammaappr,theta23,2*site-2+start+1,2*site-2+start+2);
            fpositions << 2*site-2+start+1 << " " << 2*site-2+start+2 << endl;
            fangles << theta23 << endl;
            jangle++;
            angles(jangle) = theta23;
            begpos[jangle] = 2*site-2+start+1;
            endpos[jangle] = 2*site-2+start+2;
            }

        nu(site) = Gamma(1,2);
        site++;

        Gamma = Gamma.SubMatrix(3,2*ntemp,3,2*ntemp);
        ntemp -= 1;

        }
    
    btemp -= 1;

    // Get gate for sites n-b to n-2
    while(ntemp > 1)
        {
        Matrix GammaB = Gamma.SubMatrix(1,2*btemp,1,2*btemp);
        Matrix GammaBtilde(2*btemp,2*btemp);
        Matrix WB(2*btemp,2*btemp);
    
        BlockDiagonalize(GammaB,GammaBtilde,WB);

        for(int j = btemp-1; j >= 1; j--)
            {
            int start = 2*j-1;
            Real theta14 = GetGate(Gamma,WB,1,start,start+3);
            ApplyGate(Gammaappr,theta14,2*site-2+start,2*site-2+start+3);
            fpositions << 2*site-2+start << " " << 2*site-2+start+3 << endl;
            fangles << theta14 << endl;
            jangle++;
            angles(jangle) = theta14;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+3;

            Real theta13 = GetGate(Gamma,WB,1,start,start+2);
            ApplyGate(Gammaappr,theta13,2*site-2+start,2*site-2+start+2);
            fpositions << 2*site-2+start << " " << 2*site-2+start+2 << endl;
            fangles << theta13 << endl;
            jangle++;
            angles(jangle) = theta13;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+2;

            Real theta12 = GetGate(Gamma,WB,1,start,start+1);
            ApplyGate(Gammaappr,theta12,2*site-2+start,2*site-2+start+1);
            fpositions << 2*site-2+start << " " << 2*site-2+start+1 << endl;
            fangles << theta12 << endl;
            jangle++;
            angles(jangle) = theta12;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+1;

            Real theta24 = GetGate(Gamma,WB,2,start+1,start+3);
            ApplyGate(Gammaappr,theta24,2*site-2+start+1,2*site-2+start+3);
            fpositions << 2*site-2+start+1 << " " << 2*site-2+start+3 << endl;
            fangles << theta24 << endl;
            jangle++;
            angles(jangle) = theta24;
            begpos[jangle] = 2*site-2+start+1;
            endpos[jangle] = 2*site-2+start+3;

            Real theta23 = GetGate(Gamma,WB,2,start+1,start+2);
            ApplyGate(Gammaappr,theta23,2*site-2+start+1,2*site-2+start+2);
            fpositions << 2*site-2+start+1 << " " << 2*site-2+start+2 << endl;
            fangles << theta23 << endl;
            jangle++;
            angles(jangle) = theta23;
            begpos[jangle] = 2*site-2+start+1;
            endpos[jangle] = 2*site-2+start+2;
            }

        nu(site) = Gamma(1,2);

        site++;

        Gamma = Gamma.SubMatrix(3,2*ntemp,3,2*ntemp);
        ntemp -= 1;
        btemp -= 1;

        }
   
    nu(site) = Gamma(1,2);

    int parity = 1; // Default even parity

    if(Gamma(1,2) > 0)
        {
        cout << "Ground state in even subspace" << endl;
        cout << "Project last site onto \\ket{0}" << endl;
        cout << endl;
        parity = 1;
        }
    else
        {
        cout << "Ground state in odd subspace" << endl;
        cout << "Project last site onto \\ket{1}" << endl;
        cout << endl;
        parity = 2;
        }

    if(parity == 2)
        {
        Vector vtemp = Gammatilde.Column(2*n-1);
        Gammatilde.Column(2*n-1) = Gammatilde.Column(2*n);
        Gammatilde.Column(2*n) = vtemp;
        vtemp = Gammatilde.Row(2*n-1);
        Gammatilde.Row(2*n-1) = Gammatilde.Row(2*n);
        Gammatilde.Row(2*n) = vtemp;
        nu(n) = -nu(n);
        }

    Matrix Gammaappr2 = Gammatilde;
    for(int j = 1; j <= 5*ngates; j++)
        ApplyGate(Gammaappr2,-angles(5*ngates+1-j),begpos[5*ngates+1-j],endpos[5*ngates+1-j]);

    // Calculate entanglement
    int L = n/2;
    Vector nuL(L), nuLappr(L);
    Matrix GammaL = Gammaorig.SubMatrix(1,2*L,1,2*L);
    Matrix GammaLappr = Gammaappr2.SubMatrix(1,2*L,1,2*L);
    Matrix GammaLP = GammaL;
    Matrix GammaLPappr = GammaLappr;
    Matrix WL(2*L,2*L), WLappr(2*L,2*L);
    BlockDiagonalize(GammaL, GammaLP, WL);
    BlockDiagonalize(GammaLappr, GammaLPappr, WLappr);
    Real SL = 0.0;
    Real SLappr = 0.0;
    for(int j = 1; j <= L; j++)
        {
        nuL(j) = GammaLP(2*j-1,2*j);
        nuLappr(j) = GammaLPappr(2*j-1,2*j);
        Real x = 0.5*(1.0 + nuL(j));
        Real xappr = 0.5*(1.0 + nuLappr(j));
        if(x < 1.0-1E-10)
            SL -= x*log2(x) + (1-x)*log2(1-x);
        if(xappr < 1.0-1E-10)
            SLappr -= xappr*log2(xappr) + (1-xappr)*log2(1-xappr);
        }

    // Output
    // Exact
    cout << "n = " << n << endl;
    cout << "S_L = " << SL << endl;
    cout << endl;

    // Approximate
    cout << "b = " << b << endl;
    int Chimax = 1 << (b-1);
    cout << "Chi_max = 2^(b-1) = " << Chimax << endl;

    Vector epserr(n);
    Real epsmax = 0.0;
    epserr = 0.0;
    for(int j = 1; j <= n; j++)
        {
        epserr(j) = 1.0 - nu(j);
        if(epserr(j) > epsmax)
            epsmax = epserr(j);
        }
    cout << "eps_max = " << epsmax << endl;
    cout << "log(eps_max) = " << log10(epsmax) << endl;

    cout << "S_L appr = " << SLappr << endl;
    cout << endl;

    Vector nocc(n);
    for(int j = 1; j <= n; j++)
        nocc(j) = 0.5*(1-Gammaorig(2*j-1,2*j));

    Vector noccappr(n);
    for(int j = 1; j <= n; j++)
        noccappr(j) = 0.5*(1-Gammaappr2(2*j-1,2*j));

    fangles.close();
    fpositions.close();

    //
    // Create MPS from gates
    //

    vector<ITensor> V(ngates+1);
    
    // Start with product state
    for(int j = 1; j <= n; j++)
        l[j] = Index(nameint("l",j),1,Link);

    int jgates = 1;

    VB[1] = ITensor(s[1],l[1]);
    VB[1](s[1](1),l[1](1)) = 1.0;

    for(int j = 2; j < n; j++)
        {
        VB[j] = ITensor(s[j],l[j-1],l[j]);
        VB[j](s[j](1),l[j-1](1),l[j](1)) = 1.0;
        }
    
    D = ITensor(l[n-1],l[n],1.0);
    VB[n] = ITensor(s[n],l[n]);
    VB[n](s[n](parity),l[n](1)) = 1.0;

    // Get gates
    for(int jB = b; jB <= n; jB++)
        {
        for(int j = 1; j <= b-1; j++)
            {
            Real theta14 = angles(5*(jgates-1)+1);
            Real theta13 = angles(5*(jgates-1)+2);
            Real theta12 = angles(5*(jgates-1)+3);
            Real theta24 = angles(5*(jgates-1)+4);
            Real theta23 = angles(5*(jgates-1)+5);

            ITensor Vtemp;
            int site = jB-j;

            V[jgates] = rot(s[site],s[site+1],theta14,"14");
            Vtemp = rot(s[site],s[site+1],theta13,"13");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta12,"12");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta24,"24");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta23,"23");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);

            jgates++;
            }
        }

    // Get last gates

    for(int btemp = b-1; btemp > 1; btemp--)
        {
        for(int j = 1; j <= btemp-1; j++)
            {
            Real theta14 = angles(5*(jgates-1)+1);
            Real theta13 = angles(5*(jgates-1)+2);
            Real theta12 = angles(5*(jgates-1)+3);
            Real theta24 = angles(5*(jgates-1)+4);
            Real theta23 = angles(5*(jgates-1)+5);

            ITensor Vtemp;
            int site = n-j;

            V[jgates] = rot(s[site],s[site+1],theta14,"14");
            Vtemp = rot(s[site],s[site+1],theta13,"13");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta12,"12");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta24,"24");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta23,"23");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);

            jgates++;
            }
        }

    // Gate order
    vector<int> order(ngates+1);
    
    int jg = ngates;
    for(int row = b-1; row >=1; row--)
        {
        int ord = b-row;
        int btemp = b;
        for(int j = 1; j <= n-row; j++)
            {
            order[jg] = ord;
            if(j > n-b) { btemp--; }
            ord = ord + btemp - 1;
            jg--;
            }
        }

    int jgate = 1;
    
    // Apply gates to product state

    // First row

    OptSet opts;
    opts.add("Cutoff",err);
    opts.add("Maxm",maxm);

    // Keep track of time
    clock_t t;
    t = clock();

    for(int start = 1; start < b; start++)
        {
        ITensor Vtemp;
        for(int j = n-1; j > start; j--)
            {
            // Apply gate
            Vtemp = prime(VB[j],Site)*D*prime(VB[j+1],Site)*V[order[jgate]];
            jgate++;
            svd(Vtemp,VB[j],D,VB[j+1],opts);
            l[j] = commonIndex(VB[j],D);
            l[j+1] = commonIndex(D,VB[j+1]);
            // Move center left
            Vtemp = VB[j-1]*VB[j]*D;
            svd(Vtemp,VB[j-1],D,VB[j],opts);
            l[j-1] = commonIndex(VB[j-1],D);
            l[j] = commonIndex(D,VB[j]);
            }
        Vtemp = prime(VB[start],Site)*D*prime(VB[start+1],Site)*V[order[jgate]];
        jgate++;
        svd(Vtemp,VB[start],D,VB[start+1],opts);
        l[start] = commonIndex(VB[start],D);
        l[start+1] = commonIndex(D,VB[start+1]);

        // Move center back to end
        for(int j = start+1; j < n; j++)
            {
            // Move center right
            Vtemp = D*VB[j]*VB[j+1];
            ITensor VBjtemp = D*VB[j];
            svd(Vtemp,VBjtemp,D,VB[j+1],opts);
            VB[j] = VBjtemp;
            l[j] = commonIndex(VB[j],D);
            l[j+1] = commonIndex(D,VB[j+1]);
            }

        }

    t = clock() - t;
    Real time = ((Real)t)/CLOCKS_PER_SEC;
    cout << "time (s) = " << time << endl;
    cout << "log(time) = " << log10(time) << endl;
    }

void
GMPS(Matrix Gamma, vector<ITensor>& VB, int&c, ITensor& D, vector<Index>& s, vector<Index>& l, int b)
    {
    int n = Gamma.Nrows()/2;
    // Output file
    ofstream fangles, fpositions;
    fangles.open("angles.txt");
    fpositions.open("positions.txt");

    int ngates = (b-1)*(2*n-b)/2;
    Vector angles(5*ngates);
    angles = 0.0;
    int jangle = 0;

    vector<int> begpos(5*ngates+1);
    vector<int> endpos(5*ngates+1);

    Vector nu(n);
    nu = 0.0;
    int site = 1;

    Matrix W;
    Matrix Gammatilde;
    BlockDiagonalize(Gamma,Gammatilde,W);

    Matrix Gammaorig = Gamma;
    Matrix Gammatildeorig = Gammatilde;

    Matrix Gammaappr = Gamma;

    int ntemp = n;
    int btemp = b;
    while(ntemp >= btemp)
        {
        Matrix GammaB = Gamma.SubMatrix(1,2*btemp,1,2*btemp);
        Matrix GammaBtilde(2*btemp,2*btemp);
        Matrix WB(2*btemp,2*btemp);
    
        BlockDiagonalize(GammaB,GammaBtilde,WB);

        for(int j = btemp-1; j >= 1; j--)
            {
            int start = 2*j-1;
            Real theta14 = GetGate(Gamma,WB,1,start,start+3);
            ApplyGate(Gammaappr,theta14,2*site-2+start,2*site-2+start+3);
            fpositions << 2*site-2+start << " " << 2*site-2+start+3 << endl;
            fangles << theta14 << endl;
            jangle++;
            angles(jangle) = theta14;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+3;

            Real theta13 = GetGate(Gamma,WB,1,start,start+2);
            ApplyGate(Gammaappr,theta13,2*site-2+start,2*site-2+start+2);
            fpositions << 2*site-2+start << " " << 2*site-2+start+2 << endl;
            fangles << theta13 << endl;
            jangle++;
            angles(jangle) = theta13;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+2;

            Real theta12 = GetGate(Gamma,WB,1,start,start+1);
            ApplyGate(Gammaappr,theta12,2*site-2+start,2*site-2+start+1);
            fpositions << 2*site-2+start << " " << 2*site-2+start+1 << endl;
            fangles << theta12 << endl;
            jangle++;
            angles(jangle) = theta12;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+1;

            Real theta24 = GetGate(Gamma,WB,2,start+1,start+3);
            ApplyGate(Gammaappr,theta24,2*site-2+start+1,2*site-2+start+3);
            fpositions << 2*site-2+start+1 << " " << 2*site-2+start+3 << endl;
            fangles << theta24 << endl;
            jangle++;
            angles(jangle) = theta24;
            begpos[jangle] = 2*site-2+start+1;
            endpos[jangle] = 2*site-2+start+3;

            Real theta23 = GetGate(Gamma,WB,2,start+1,start+2);
            ApplyGate(Gammaappr,theta23,2*site-2+start+1,2*site-2+start+2);
            fpositions << 2*site-2+start+1 << " " << 2*site-2+start+2 << endl;
            fangles << theta23 << endl;
            jangle++;
            angles(jangle) = theta23;
            begpos[jangle] = 2*site-2+start+1;
            endpos[jangle] = 2*site-2+start+2;
            }

        nu(site) = Gamma(1,2);
        site++;

        Gamma = Gamma.SubMatrix(3,2*ntemp,3,2*ntemp);
        ntemp -= 1;

        }
    
    btemp -= 1;

    // Get gate for sites n-b to n-2
    while(ntemp > 1)
        {
        Matrix GammaB = Gamma.SubMatrix(1,2*btemp,1,2*btemp);
        Matrix GammaBtilde(2*btemp,2*btemp);
        Matrix WB(2*btemp,2*btemp);
    
        BlockDiagonalize(GammaB,GammaBtilde,WB);

        for(int j = btemp-1; j >= 1; j--)
            {
            int start = 2*j-1;
            Real theta14 = GetGate(Gamma,WB,1,start,start+3);
            ApplyGate(Gammaappr,theta14,2*site-2+start,2*site-2+start+3);
            fpositions << 2*site-2+start << " " << 2*site-2+start+3 << endl;
            fangles << theta14 << endl;
            jangle++;
            angles(jangle) = theta14;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+3;

            Real theta13 = GetGate(Gamma,WB,1,start,start+2);
            ApplyGate(Gammaappr,theta13,2*site-2+start,2*site-2+start+2);
            fpositions << 2*site-2+start << " " << 2*site-2+start+2 << endl;
            fangles << theta13 << endl;
            jangle++;
            angles(jangle) = theta13;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+2;

            Real theta12 = GetGate(Gamma,WB,1,start,start+1);
            ApplyGate(Gammaappr,theta12,2*site-2+start,2*site-2+start+1);
            fpositions << 2*site-2+start << " " << 2*site-2+start+1 << endl;
            fangles << theta12 << endl;
            jangle++;
            angles(jangle) = theta12;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+1;

            Real theta24 = GetGate(Gamma,WB,2,start+1,start+3);
            ApplyGate(Gammaappr,theta24,2*site-2+start+1,2*site-2+start+3);
            fpositions << 2*site-2+start+1 << " " << 2*site-2+start+3 << endl;
            fangles << theta24 << endl;
            jangle++;
            angles(jangle) = theta24;
            begpos[jangle] = 2*site-2+start+1;
            endpos[jangle] = 2*site-2+start+3;

            Real theta23 = GetGate(Gamma,WB,2,start+1,start+2);
            ApplyGate(Gammaappr,theta23,2*site-2+start+1,2*site-2+start+2);
            fpositions << 2*site-2+start+1 << " " << 2*site-2+start+2 << endl;
            fangles << theta23 << endl;
            jangle++;
            angles(jangle) = theta23;
            begpos[jangle] = 2*site-2+start+1;
            endpos[jangle] = 2*site-2+start+2;
            }

        nu(site) = Gamma(1,2);

        site++;

        Gamma = Gamma.SubMatrix(3,2*ntemp,3,2*ntemp);
        ntemp -= 1;
        btemp -= 1;

        }
   
    nu(site) = Gamma(1,2);

    int parity = 1; // Default even parity

    if(Gamma(1,2) > 0)
        {
        cout << "Ground state in even subspace" << endl;
        cout << "Project last site onto \\ket{0}" << endl;
        cout << endl;
        parity = 1;
        }
    else
        {
        cout << "Ground state in odd subspace" << endl;
        cout << "Project last site onto \\ket{1}" << endl;
        cout << endl;
        parity = 2;
        }

    if(parity == 2)
        {
        Vector vtemp = Gammatilde.Column(2*n-1);
        Gammatilde.Column(2*n-1) = Gammatilde.Column(2*n);
        Gammatilde.Column(2*n) = vtemp;
        vtemp = Gammatilde.Row(2*n-1);
        Gammatilde.Row(2*n-1) = Gammatilde.Row(2*n);
        Gammatilde.Row(2*n) = vtemp;
        nu(n) = -nu(n);
        }

    Matrix Gammaappr2 = Gammatilde;
    for(int j = 1; j <= 5*ngates; j++)
        ApplyGate(Gammaappr2,-angles(5*ngates+1-j),begpos[5*ngates+1-j],endpos[5*ngates+1-j]);

    // Calculate entanglement
    int L = n/2;
    Vector nuL(L), nuLappr(L);
    Matrix GammaL = Gammaorig.SubMatrix(1,2*L,1,2*L);
    Matrix GammaLappr = Gammaappr2.SubMatrix(1,2*L,1,2*L);
    Matrix GammaLP = GammaL;
    Matrix GammaLPappr = GammaLappr;
    Matrix WL(2*L,2*L), WLappr(2*L,2*L);
    BlockDiagonalize(GammaL, GammaLP, WL);
    BlockDiagonalize(GammaLappr, GammaLPappr, WLappr);
    Real SL = 0.0;
    Real SLappr = 0.0;
    for(int j = 1; j <= L; j++)
        {
        nuL(j) = GammaLP(2*j-1,2*j);
        nuLappr(j) = GammaLPappr(2*j-1,2*j);
        Real x = 0.5*(1.0 + nuL(j));
        Real xappr = 0.5*(1.0 + nuLappr(j));
        if(x < 1.0-1E-10)
            SL -= x*log2(x) + (1-x)*log2(1-x);
        if(xappr < 1.0-1E-10)
            SLappr -= xappr*log2(xappr) + (1-xappr)*log2(1-xappr);
        }

    // Output
    // Exact
    cout << "n = " << n << endl;
    cout << "S_L = " << SL << endl;
    cout << endl;

    // Approximate
    cout << "b = " << b << endl;
    int Chimax = 1 << (b-1);
    cout << "Chi_max = 2^(b-1) = " << Chimax << endl;

    Vector epserr(n);
    Real epsmax = 0.0;
    epserr = 0.0;
    for(int j = 1; j <= n; j++)
        {
        epserr(j) = 1.0 - nu(j);
        if(epserr(j) > epsmax)
            epsmax = epserr(j);
        }
    cout << "eps_max = " << epsmax << endl;
    cout << "log(eps_max) = " << log10(epsmax) << endl;

    cout << "S_L appr = " << SLappr << endl;
    cout << endl;

    Vector nocc(n);
    for(int j = 1; j <= n; j++)
        nocc(j) = 0.5*(1-Gammaorig(2*j-1,2*j));

    Vector noccappr(n);
    for(int j = 1; j <= n; j++)
        noccappr(j) = 0.5*(1-Gammaappr2(2*j-1,2*j));

    fangles.close();
    fpositions.close();

    //
    // Create MPS from gates
    //

    vector<ITensor> V(ngates+1);
    
    // Start with product state
    for(int j = 1; j <= n; j++)
        l[j] = Index(nameint("l",j),1,Link);

    int jgates = 1;

    VB[1] = ITensor(s[1],l[1]);
    VB[1](s[1](1),l[1](1)) = 1.0;

    for(int j = 2; j < n; j++)
        {
        VB[j] = ITensor(s[j],l[j-1],l[j]);
        VB[j](s[j](1),l[j-1](1),l[j](1)) = 1.0;
        }
    
    D = ITensor(l[n-1],l[n],1.0);
    VB[n] = ITensor(s[n],l[n]);
    VB[n](s[n](parity),l[n](1)) = 1.0;

    // Get gates
    for(int jB = b; jB <= n; jB++)
        {
        for(int j = 1; j <= b-1; j++)
            {
            Real theta14 = angles(5*(jgates-1)+1);
            Real theta13 = angles(5*(jgates-1)+2);
            Real theta12 = angles(5*(jgates-1)+3);
            Real theta24 = angles(5*(jgates-1)+4);
            Real theta23 = angles(5*(jgates-1)+5);

            ITensor Vtemp;
            int site = jB-j;

            V[jgates] = rot(s[site],s[site+1],theta14,"14");
            Vtemp = rot(s[site],s[site+1],theta13,"13");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta12,"12");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta24,"24");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta23,"23");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);

            jgates++;
            }
        }

    // Get last gate

    for(int btemp = b-1; btemp > 1; btemp--)
        {
        for(int j = 1; j <= btemp-1; j++)
            {
            Real theta14 = angles(5*(jgates-1)+1);
            Real theta13 = angles(5*(jgates-1)+2);
            Real theta12 = angles(5*(jgates-1)+3);
            Real theta24 = angles(5*(jgates-1)+4);
            Real theta23 = angles(5*(jgates-1)+5);

            ITensor Vtemp;
            int site = n-j;

            V[jgates] = rot(s[site],s[site+1],theta14,"14");
            Vtemp = rot(s[site],s[site+1],theta13,"13");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta12,"12");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta24,"24");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta23,"23");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);

            jgates++;
            }
        }

    // Gate order
    vector<int> order(ngates+1);
    
    int jg = ngates;
    for(int row = b-1; row >=1; row--)
        {
        int ord = b-row;
        int btemp = b;
        for(int j = 1; j <= n-row; j++)
            {
            order[jg] = ord;
            if(j > n-b) { btemp--; }
            ord = ord + btemp - 1;
            jg--;
            }
        }

    int jgate = 1;
    
    // Apply gates to product state

    // First row

    OptSet opts;
    opts.add("Cutoff",err);
    opts.add("Maxm",maxm);

    // Keep track of time
    clock_t t;
    t = clock();

    for(int start = 1; start < b; start++)
        {
        ITensor Vtemp;
        for(int j = n-1; j > start; j--)
            {
            // Apply gate
            Vtemp = prime(VB[j],Site)*D*prime(VB[j+1],Site)*V[order[jgate]];
            jgate++;
            svd(Vtemp,VB[j],D,VB[j+1],opts);
            l[j] = commonIndex(VB[j],D);
            l[j+1] = commonIndex(D,VB[j+1]);
            // Move center left
            Vtemp = VB[j-1]*VB[j]*D;
            svd(Vtemp,VB[j-1],D,VB[j],opts);
            l[j-1] = commonIndex(VB[j-1],D);
            l[j] = commonIndex(D,VB[j]);
            }
        Vtemp = prime(VB[start],Site)*D*prime(VB[start+1],Site)*V[order[jgate]];
        jgate++;
        svd(Vtemp,VB[start],D,VB[start+1],opts);
        l[start] = commonIndex(VB[start],D);
        l[start+1] = commonIndex(D,VB[start+1]);

        // Move center back to end
        for(int j = start+1; j < n; j++)
            {
            // Move center right
            Vtemp = D*VB[j]*VB[j+1];
            ITensor VBjtemp = D*VB[j];
            svd(Vtemp,VBjtemp,D,VB[j+1],opts);
            VB[j] = VBjtemp;
            l[j] = commonIndex(VB[j],D);
            l[j+1] = commonIndex(D,VB[j+1]);
            }

        }

    t = clock() - t;
    Real time = ((Real)t)/CLOCKS_PER_SEC;
    cout << "time (s) = " << time << endl;
    cout << "log(time) = " << log10(time) << endl;

    }

Real
GMPS(Matrix Hre, Matrix Him, Matrix Gre, Matrix Gim, vector<ITensor>& VB, int&c, ITensor& D, vector<Index>& s, vector<Index>& l, int b)
    {
    // Number of sites
    // Number of Majorana modes = 2n
    int n = Hre.Nrows();

    // Form matrix [A_{j,k}] (Hamiltonian in Majorana basis)
    // H = \frac{i}{2} \sum_{j,k=1}^{2n} c_j A_{j,k} c_k 
    Matrix A(2*n,2*n),
           Atilde(2*n,2*n),
           W(2*n,2*n);
    A = 0; Atilde = 0; W = 0;

    A.SubMatrix(1,n,1,n) = Him + Gim;
    A.SubMatrix(1,n,n+1,2*n) = Hre + Gre;
    A.SubMatrix(n+1,2*n,1,n) = -Hre + Gre;
    A.SubMatrix(n+1,2*n,n+1,2*n) = Him - Gim;

    //cout << "A = " << endl;
    //cout << A;

    Matrix Id(2*n,2*n);
    Id = 1.0;

    //cout << "Id = " << endl;
    //cout << Id;
    
    Matrix P(2*n,2*n);
    P = 0.0;
   
    for(int j = 1; j <= n; j++)
        P.Row(j) = Id.Row(2*j-1);
    for(int j = 1; j <= n; j++)
        P.Row(n+j) = Id.Row(2*j);

    cout << "P = " << endl;
    cout << P;
    
    A = P.t() * A * P;
    //cout << "A = " << endl;
    //cout << A;

    // Get \tilde{A} and W such that
    // A = W \tilde{A} W^T
    BlockDiagonalize(A,Atilde,W);
    
    /*
    cout << "A = " << endl;
    cout << A;
    cout << "\\tilde{A} = " << endl;
    cout << Atilde;
    cout << "W = " << endl;
    cout << W;
    cout << "\\det(W) = " << endl;
    cout << Determinant(W) << endl;
    cout << endl;
    */

    Vector eps(n);
    for(int j = 1; j <= n; j++)
        eps(j) = Atilde(2*j-1,2*j);

    // Get \Gamma matrix
    // \Gamma = W \tilde{\Gamma} W^T
    Matrix Gammatilde(2*n,2*n),
           Gamma(2*n,2*n);
    Gammatilde = 0; Gamma = 0;

    for(int j = 1; j <= n; j++)
        Gammatilde(2*j-1,2*j) = 1.0;
    Gammatilde -= Gammatilde.t();

    Gamma = W*Gammatilde*W.t();

    //cout << "\\tilde{\\Gamma} = " << endl;
    //cout << Gammatilde;
    cout << "\\Gamma = W \\tilde{\\Gamma} W^T = " << endl;
    cout << Gamma;
    
    //cout << A;
    //cout << Hre;

    // TEBD algorithm
    // Submatrix in block B of \Gamma
  
    // Output file
    ofstream fangles, fpositions;
    fangles.open("angles.txt");
    fpositions.open("positions.txt");

    int ngates = (b-1)*(2*n-b)/2;
    //cout << "ngates = " << ngates << endl;
    //cout << endl;
    Vector angles(5*ngates);
    angles = 0.0;
    int jangle = 0;

    vector<int> begpos(5*ngates+1);
    vector<int> endpos(5*ngates+1);

    Vector nu(n);
    nu = 0.0;
    int site = 1;

    Matrix Aorig = A;
    Matrix Gammaorig = Gamma;
    Matrix Atildeorig = Atilde;
    Matrix Gammatildeorig = Gammatilde;

    Matrix Gammaappr = Gamma;
    Matrix Aappr = A;

    int ntemp = n;
    int btemp = b;
    while(ntemp >= btemp)
        {
        Matrix GammaB = Gamma.SubMatrix(1,2*btemp,1,2*btemp);
        Matrix GammaBtilde(2*btemp,2*btemp);
        Matrix WB(2*btemp,2*btemp);
    
        BlockDiagonalize(GammaB,GammaBtilde,WB);
        //cout << "\\Gamma_B = " << endl;
        //cout << GammaB;
        //cout << "\\tilde{\\Gamma}_B = " << endl;
        //cout << GammaBtilde;
        //cout << "W_B = " << endl;
        //cout << WB;
        //cout << "\\det(W_B) = " << endl;
        //cout << Determinant(WB) << endl;
        //cout << endl;

        for(int j = btemp-1; j >= 1; j--)
            {
            int start = 2*j-1;
            Real theta14 = GetGate(Gamma,WB,1,start,start+3);
            ApplyGate(Gammaappr,theta14,2*site-2+start,2*site-2+start+3);
            fpositions << 2*site-2+start << " " << 2*site-2+start+3 << endl;
            ApplyGate(Aappr,theta14,2*site-2+start,2*site-2+start+3);
            fangles << theta14 << endl;
            jangle++;
            angles(jangle) = theta14;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+3;

            Real theta13 = GetGate(Gamma,WB,1,start,start+2);
            ApplyGate(Gammaappr,theta13,2*site-2+start,2*site-2+start+2);
            fpositions << 2*site-2+start << " " << 2*site-2+start+2 << endl;
            ApplyGate(Aappr,theta13,2*site-2+start,2*site-2+start+2);
            fangles << theta13 << endl;
            jangle++;
            angles(jangle) = theta13;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+2;

            Real theta12 = GetGate(Gamma,WB,1,start,start+1);
            ApplyGate(Gammaappr,theta12,2*site-2+start,2*site-2+start+1);
            fpositions << 2*site-2+start << " " << 2*site-2+start+1 << endl;
            ApplyGate(Aappr,theta12,2*site-2+start,2*site-2+start+1);
            fangles << theta12 << endl;
            jangle++;
            angles(jangle) = theta12;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+1;

            Real theta24 = GetGate(Gamma,WB,2,start+1,start+3);
            ApplyGate(Gammaappr,theta24,2*site-2+start+1,2*site-2+start+3);
            fpositions << 2*site-2+start+1 << " " << 2*site-2+start+3 << endl;
            ApplyGate(Aappr,theta24,2*site-2+start+1,2*site-2+start+3);
            fangles << theta24 << endl;
            jangle++;
            angles(jangle) = theta24;
            begpos[jangle] = 2*site-2+start+1;
            endpos[jangle] = 2*site-2+start+3;

            Real theta23 = GetGate(Gamma,WB,2,start+1,start+2);
            ApplyGate(Gammaappr,theta23,2*site-2+start+1,2*site-2+start+2);
            fpositions << 2*site-2+start+1 << " " << 2*site-2+start+2 << endl;
            ApplyGate(Aappr,theta23,2*site-2+start+1,2*site-2+start+2);
            fangles << theta23 << endl;
            jangle++;
            angles(jangle) = theta23;
            begpos[jangle] = 2*site-2+start+1;
            endpos[jangle] = 2*site-2+start+2;
            }

        //cout << Gamma;
        //GammaP.SubMatrix(2*site-1,2*n,2*site-1,2*n) = Gamma;
        //AP.SubMatrix(2*site-1,2*n,2*site-1,2*n) = APtemp;
        //cout << "Gammaappr = " << endl;
        //cout << Gammaappr;
        //cout << "Aappr = " << endl;
        //cout << Aappr;

        nu(site) = Gamma(1,2);
        site++;

        //cout << "\\Gamma_B = " << endl;
        //cout << Gamma.SubMatrix(1,2*btemp,1,2*btemp) << endl;

        Gamma = Gamma.SubMatrix(3,2*ntemp,3,2*ntemp);
        //APtemp = APtemp.SubMatrix(3,2*ntemp,3,2*ntemp);
        ntemp -= 1;

        //cout << "\\Gamma = " << endl;
        //cout << Gamma << endl;
        }
    
    //cout << "GammaP = " << endl;
    //cout << GammaP;

    btemp -= 1;

    // Get gate for sites n-b to n-2
    while(ntemp > 1)
        {
        Matrix GammaB = Gamma.SubMatrix(1,2*btemp,1,2*btemp);
        Matrix GammaBtilde(2*btemp,2*btemp);
        Matrix WB(2*btemp,2*btemp);
    
        BlockDiagonalize(GammaB,GammaBtilde,WB);
        //cout << "\\Gamma_B = " << endl;
        //cout << GammaB;
        //cout << "\\tilde{\\Gamma}_B = " << endl;
        //cout << GammaBtilde;
        //cout << "W_B = " << endl;
        //cout << WB;
        //cout << "\\det(W_B) = " << endl;
        //cout << Determinant(WB) << endl;
        //cout << endl;

        for(int j = btemp-1; j >= 1; j--)
            {
            int start = 2*j-1;
            Real theta14 = GetGate(Gamma,WB,1,start,start+3);
            ApplyGate(Gammaappr,theta14,2*site-2+start,2*site-2+start+3);
            fpositions << 2*site-2+start << " " << 2*site-2+start+3 << endl;
            ApplyGate(Aappr,theta14,2*site-2+start,2*site-2+start+3);
            fangles << theta14 << endl;
            jangle++;
            angles(jangle) = theta14;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+3;

            Real theta13 = GetGate(Gamma,WB,1,start,start+2);
            ApplyGate(Gammaappr,theta13,2*site-2+start,2*site-2+start+2);
            fpositions << 2*site-2+start << " " << 2*site-2+start+2 << endl;
            ApplyGate(Aappr,theta13,2*site-2+start,2*site-2+start+2);
            fangles << theta13 << endl;
            jangle++;
            angles(jangle) = theta13;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+2;

            Real theta12 = GetGate(Gamma,WB,1,start,start+1);
            ApplyGate(Gammaappr,theta12,2*site-2+start,2*site-2+start+1);
            fpositions << 2*site-2+start << " " << 2*site-2+start+1 << endl;
            ApplyGate(Aappr,theta12,2*site-2+start,2*site-2+start+1);
            fangles << theta12 << endl;
            jangle++;
            angles(jangle) = theta12;
            begpos[jangle] = 2*site-2+start;
            endpos[jangle] = 2*site-2+start+1;

            Real theta24 = GetGate(Gamma,WB,2,start+1,start+3);
            ApplyGate(Gammaappr,theta24,2*site-2+start+1,2*site-2+start+3);
            fpositions << 2*site-2+start+1 << " " << 2*site-2+start+3 << endl;
            ApplyGate(Aappr,theta24,2*site-2+start+1,2*site-2+start+3);
            fangles << theta24 << endl;
            jangle++;
            angles(jangle) = theta24;
            begpos[jangle] = 2*site-2+start+1;
            endpos[jangle] = 2*site-2+start+3;

            Real theta23 = GetGate(Gamma,WB,2,start+1,start+2);
            ApplyGate(Gammaappr,theta23,2*site-2+start+1,2*site-2+start+2);
            fpositions << 2*site-2+start+1 << " " << 2*site-2+start+2 << endl;
            ApplyGate(Aappr,theta23,2*site-2+start+1,2*site-2+start+2);
            fangles << theta23 << endl;
            jangle++;
            angles(jangle) = theta23;
            begpos[jangle] = 2*site-2+start+1;
            endpos[jangle] = 2*site-2+start+2;
            }

        //GammaP.SubMatrix(2*site-1,2*n,2*site-1,2*n) = Gamma;
        //AP.SubMatrix(2*site-1,2*n,2*site-1,2*n) = APtemp;
        //cout << "Gammaappr = " << endl;
        //cout << Gammaappr;
        //cout << "Aappr = " << endl;
        //cout << Aappr;

        nu(site) = Gamma(1,2);

        site++;
        //cout << "\\Gamma_B = " << endl;
        //cout << Gamma.SubMatrix(1,2*btemp,1,2*btemp) << endl;

        Gamma = Gamma.SubMatrix(3,2*ntemp,3,2*ntemp);
        //APtemp = APtemp.SubMatrix(3,2*ntemp,3,2*ntemp);
        ntemp -= 1;
        btemp -= 1;

        //cout << "\\Gamma = " << endl;
        //cout << Gamma << endl;
        }
   
    nu(site) = Gamma(1,2);

    int parity = 1; // Default even parity

    if(Gamma(1,2) > 0)
        {
        cout << "Ground state in even subspace" << endl;
        cout << "Project last site onto \\ket{0}" << endl;
        cout << endl;
        parity = 1;
        }
    else
        {
        cout << "Ground state in odd subspace" << endl;
        cout << "Project last site onto \\ket{1}" << endl;
        cout << endl;
        parity = 2;
        }

    //cout << "Atilde = " << endl;
    //cout << Atilde;

    //cout << "AP = " << endl;
    //cout << AP;

    if(parity == 2)
        {
        Vector vtemp = Gammatilde.Column(2*n-1);
        Gammatilde.Column(2*n-1) = Gammatilde.Column(2*n);
        Gammatilde.Column(2*n) = vtemp;
        vtemp = Gammatilde.Row(2*n-1);
        Gammatilde.Row(2*n-1) = Gammatilde.Row(2*n);
        Gammatilde.Row(2*n) = vtemp;
        nu(n) = -nu(n);
        }

    /*
    cout << "Gammaappr = " << endl;
    cout << Gammaappr;
    for(int j = 1; j <= 2*n; j++)
        {
        Vector vtemp = Gammaappr.Column(j);
        Real norm = 0.0;
        for(int k = 1; k <= 2*n; k++)
            {
            norm += vtemp(k);
            }
        cout << norm << endl;
        }
    */

    //cout << "Gammatilde = " << endl;
    //cout << Gammatilde;
    
    /*
    Matrix GammaPP(2*n,2*n);
    GammaPP = 0.0;
    for(int j = 1; j <= n; j++)
        GammaPP(2*j-1,2*j) = nu(j);
    GammaPP -= GammaPP.t();
    
    cout << "GammaPP = " << endl;
    cout << GammaPP;
    */

    Matrix Gammaappr2 = Gammatilde;
    for(int j = 1; j <= 5*ngates; j++)
        ApplyGate(Gammaappr2,-angles(5*ngates+1-j),begpos[5*ngates+1-j],endpos[5*ngates+1-j]);

    // Calculate entanglement
    int L = n/2;
    Vector nuL(L), nuLappr(L);
    Matrix GammaL = Gammaorig.SubMatrix(1,2*L,1,2*L);
    Matrix GammaLappr = Gammaappr2.SubMatrix(1,2*L,1,2*L);
    Matrix GammaLP = GammaL;
    Matrix GammaLPappr = GammaLappr;
    Matrix WL(2*L,2*L), WLappr(2*L,2*L);
    BlockDiagonalize(GammaL, GammaLP, WL);
    BlockDiagonalize(GammaLappr, GammaLPappr, WLappr);
    Real SL = 0.0;
    Real SLappr = 0.0;
    for(int j = 1; j <= L; j++)
        {
        nuL(j) = GammaLP(2*j-1,2*j);
        nuLappr(j) = GammaLPappr(2*j-1,2*j);
        Real x = 0.5*(1.0 + nuL(j));
        Real xappr = 0.5*(1.0 + nuLappr(j));
        if(x < 1.0-1E-10)
            SL -= x*log2(x) + (1-x)*log2(1-x);
        if(xappr < 1.0-1E-10)
            SLappr -= xappr*log2(xappr) + (1-xappr)*log2(1-xappr);
        }

    // Calculate energy
    Real Egs = 0.25*Trace(Atildeorig*Gammatildeorig) + 0.5*Trace(Hre);
    
    // Output
    // Exact
    cout << "n = " << n << endl;
    //cout << "E = 1/4 Tr(Atilde \\Gammatilde) + 1/2 Tr(Hre) = " << Egs << endl;
    cout << "E = " << Egs << endl;
    cout << "S_L = " << SL << endl;
    cout << endl;

    // Approximate
    cout << "b = " << b << endl;
    int Chimax = 1 << (b-1);
    cout << "Chi_max = 2^(b-1) = " << Chimax << endl;

    //cout << "Number of sites = " << site << endl;
    //cout << "\\nu_k = " << endl;
    //cout << nu;
    Vector epserr(n);
    Real epsmax = 0.0;
    epserr = 0.0;
    for(int j = 1; j <= n; j++)
        {
        epserr(j) = 1.0 - nu(j);
        if(epserr(j) > epsmax)
            epsmax = epserr(j);
        }
    //cout << "\\epserr_k = " << endl;
    //cout << epserr;
    cout << "eps_max = " << epsmax << endl;
    cout << "log(eps_max) = " << log10(epsmax) << endl;

    //cout << "Energy:" << endl;
    Real Eappr = 0.25*Trace(Aappr*Gammatilde) + 0.5*Trace(Hre);
    
    //cout << "Approximate energy = 0.25*Trace(Aappr*Gammatilde) + 0.5*Trace(Hre) = " << Eappr << endl;
    cout << "Eappr = " << Eappr << endl;
    Real delE = fabs(Egs-Eappr)/fabs(Egs);
    cout << "|delE / E| = " << delE << endl;
    cout << "log|delE / E| = " << log10(delE) << endl;
    cout << "S_L appr = " << SLappr << endl;
    cout << endl;

    //cout << "Occupation:" << endl;
    //cout << "Exact occupation:" << endl;
    Vector nocc(n);
    for(int j = 1; j <= n; j++)
        nocc(j) = 0.5*(1-Gammaorig(2*j-1,2*j));
    //cout << "<n_j> exact = " << endl;
    //cout << nocc;

    //cout << "Appr occupation:" << endl;
    Vector noccappr(n);
    for(int j = 1; j <= n; j++)
        noccappr(j) = 0.5*(1-Gammaappr2(2*j-1,2*j));
    //cout << "<n_j> appr = " << endl;
    //cout << noccappr;

    //cout << "Entanglement:" << endl;
    //cout << "Exact entanglement:" << endl;
    //cout << "nuL = " << endl;
    //cout << nuL;
    //cout << "nuLappr = " << endl;
    //cout << nuLappr;

    /*
    Matrix Gammatildeappr = Gammaorig;
    for(int j = 1; j <= 5*ngates; j++)
        ApplyGate(Gammatildeappr,angles(j),begpos[j],endpos[j]);
    for(int j = 1; j <= n; j++)
        cout << Gammatildeappr(2*j-1,2*j) << endl;
    */

    fangles.close();
    fpositions.close();

    //
    // Create MPS from gates
    //

    vector<ITensor> V(ngates+1);
    
    // Start with product state
    for(int j = 1; j <= n; j++)
        l[j] = Index(nameint("l",j),1,Link);

    int jgates = 1;

    VB[1] = ITensor(s[1],l[1]);
    VB[1](s[1](1),l[1](1)) = 1.0;

    for(int j = 2; j < n; j++)
        {
        VB[j] = ITensor(s[j],l[j-1],l[j]);
        VB[j](s[j](1),l[j-1](1),l[j](1)) = 1.0;
        }
    
    D = ITensor(l[n-1],l[n],1.0);
    VB[n] = ITensor(s[n],l[n]);
    VB[n](s[n](parity),l[n](1)) = 1.0;

    // Product state wavefunction
    /*
    ITensor psi = VB[1];
    for(int j = 2; j < n; j++)
        {
        psi = psi*VB[j];    
        }
    psi = psi*D;
    psi = psi*VB[n];
    PrintDat(psi);
    */

    // Get gates
    for(int jB = b; jB <= n; jB++)
        {
        for(int j = 1; j <= b-1; j++)
            {
            Real theta14 = angles(5*(jgates-1)+1);
            Real theta13 = angles(5*(jgates-1)+2);
            Real theta12 = angles(5*(jgates-1)+3);
            Real theta24 = angles(5*(jgates-1)+4);
            Real theta23 = angles(5*(jgates-1)+5);

            ITensor Vtemp;
            int site = jB-j;

            V[jgates] = rot(s[site],s[site+1],theta14,"14");
            Vtemp = rot(s[site],s[site+1],theta13,"13");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta12,"12");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta24,"24");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta23,"23");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            //PrintDat(V[jgates]);

            jgates++;
            }
        }

    // Get last gate

    for(int btemp = b-1; btemp > 1; btemp--)
        {
        for(int j = 1; j <= btemp-1; j++)
            {
            Real theta14 = angles(5*(jgates-1)+1);
            Real theta13 = angles(5*(jgates-1)+2);
            Real theta12 = angles(5*(jgates-1)+3);
            Real theta24 = angles(5*(jgates-1)+4);
            Real theta23 = angles(5*(jgates-1)+5);

            ITensor Vtemp;
            int site = n-j;

            V[jgates] = rot(s[site],s[site+1],theta14,"14");
            Vtemp = rot(s[site],s[site+1],theta13,"13");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta12,"12");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta24,"24");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            Vtemp = rot(s[site],s[site+1],theta23,"23");
            V[jgates] = V[jgates]*prime(Vtemp);
            V[jgates].mapprime(2,1);
            //PrintDat(V[jgates]);

            jgates++;
            }
        }
    //cout << "Number of gates = " << jgates-1 << endl;

    // Gate order
    vector<int> order(ngates+1);
    
    int jg = ngates;
    for(int row = b-1; row >=1; row--)
        {
        int ord = b-row;
        int btemp = b;
        for(int j = 1; j <= n-row; j++)
            {
            order[jg] = ord;
            if(j > n-b) { btemp--; }
            ord = ord + btemp - 1;
            jg--;
            }
        }

    //for(int j = 1; j <= ngates; j++)
    //    cout << "ord[" << j << "] = " << order[j] << "\n" << endl;
    int jgate = 1;
    
    // Apply gates to product state

    // First row
    //int start = 1;

    OptSet opts;
    opts.add("Cutoff",err);
    opts.add("Maxm",maxm);

    // Keep track of time
    clock_t t;
    t = clock();

    for(int start = 1; start < b; start++)
        {
        //cout << "Row " << start << endl;
        //cout << endl;
        ITensor Vtemp;
        for(int j = n-1; j > start; j--)
            {
            // Apply gate
            //cout << "j = " << j << endl;
            Vtemp = prime(VB[j],Site)*D*prime(VB[j+1],Site)*V[order[jgate]];
            jgate++;
            svd(Vtemp,VB[j],D,VB[j+1],opts);
            l[j] = commonIndex(VB[j],D);
            l[j+1] = commonIndex(D,VB[j+1]);
            // Move center left
            Vtemp = VB[j-1]*VB[j]*D;
            svd(Vtemp,VB[j-1],D,VB[j],opts);
            l[j-1] = commonIndex(VB[j-1],D);
            l[j] = commonIndex(D,VB[j]);
            }
        Vtemp = prime(VB[start],Site)*D*prime(VB[start+1],Site)*V[order[jgate]];
        jgate++;
        svd(Vtemp,VB[start],D,VB[start+1],opts);
        l[start] = commonIndex(VB[start],D);
        l[start+1] = commonIndex(D,VB[start+1]);

        // Measure occupation
        //ITensor dj(s[start]);
        //dj(s[start](2)) = 1.0;
        //ITensor nj = dj*VB[start]*D;
        //PrintDat((nj*dag(nj)).toComplex());
        //cout << "n_" << start << " = " << 2.0/(n+1)*sin(pi*start/(n+1))*sin(pi*start/(n+1)) << endl;
        //cout << endl;

        // Move center back to end
        for(int j = start+1; j < n; j++)
            {
            // Move center right
            Vtemp = D*VB[j]*VB[j+1];
            ITensor VBjtemp = D*VB[j];
            svd(Vtemp,VBjtemp,D,VB[j+1],opts);
            VB[j] = VBjtemp;
            l[j] = commonIndex(VB[j],D);
            l[j+1] = commonIndex(D,VB[j+1]);
            // Measure occupation
            //dj = ITensor(s[j]);
            //dj(s[j](2)) = 1.0;
            //nj = dj*VB[j]*D;
            //PrintDat((nj*dag(nj)).toComplex());
            //cout << "n_" << j << " = " << 2.0/(n+1)*sin(pi*j/(n+1))*sin(pi*j/(n+1)) << endl;
            //cout << endl;
            }

        // Measure occupation
        //dj = ITensor(s[n]);
        //dj(s[n](2)) = 1.0;
        //nj = dj*VB[n]*D;
        //PrintDat((nj*dag(nj)).toComplex());
        //cout << "n_" << n << " = " << 2.0/(n+1)*sin(pi*n/(n+1))*sin(pi*n/(n+1)) << endl;
        //cout << endl;
        }

    t = clock() - t;
    Real time = ((Real)t)/CLOCKS_PER_SEC;
    cout << "time (s) = " << time << endl;
    cout << "log(time) = " << log10(time) << endl;

    // Get wavefunction
    /*
    psi = VB[1];
    for(int j = 2; j < n; j++)
        {
        psi = psi*VB[j];    
        }
    psi = psi*D;
    psi = psi*VB[n];
    PrintDat(psi);
    */

    return Egs;
    }

void
shiftCenter(vector<ITensor>& A, int& c, ITensor& D, vector<Index> s, vector<Index>& l, int shift)
    {
    ITensor psi;
    OptSet opts;
    opts.add("Cutoff",err);
    opts.add("Maxm",maxm);

    if(shift > 0)
        {
        for(int k = 1; k <= shift; k++)
            {
            c++;
            A[c] = D*A[c];
            psi = A[c]*A[c+1]; 
            svd(psi,A[c],D,A[c+1],opts);
            l[c] = commonIndex(A[c],D);
            l[c+1] = commonIndex(A[c+1],D);
            }
        }
    else if(shift < 0)
        {
        for(int k = 1; k <= -shift; k++)
            {
            c--;
            A[c+1] = A[c+1]*D;
            psi = A[c]*A[c+1];
            svd(psi,A[c],D,A[c+1],opts);
            l[c] = commonIndex(A[c],D);
            l[c+1] = commonIndex(A[c+1],D);
            }
        }
    }

Vector
getOcc(vector<ITensor>& A, int& c, ITensor& D, vector<Index> s, vector<Index>& l)
    {
    shiftCenter(A,c,D,s,l,1-c);

    int n = s.size() - 1;
    Vector res(n);
    res = 0.0;

    for(int j = 1; j <= n; j++)
        {
        ITensor dj(s[j]);
        dj(s[j](2)) = 1.0;
        ITensor nj = dj*A[j]*D;
        Complex occ = (nj*dag(nj)).toComplex();
        res(j) = occ.real();
        if(j < n-1) { shiftCenter(A,c,D,s,l,1); }
        //cout << "res(" << j << ") = " << res(j) << endl;
        //cout << endl;
        }

    //cout << "n_" << j << " = " << 2.0/(n+1)*sin(pi*j/(n+1))*sin(pi*j/(n+1)) << endl;
    //cout << endl;
    return res;
    }
            
Matrix
makePotential(int nx, int ny, int natomx, int natomy, Real mu, Real V0, string potential)
    {
    Matrix V(ny,nx);
    // Atom spacing
    Real ax = (Real)(nx-1)/natomx;
    Real ay = (Real)(ny-1)/natomy;
    for(int jx = 1; jx <= nx; jx++)
        {
        for(int jy = 1; jy <= ny; jy++)
            {
            V(jy,jx) = mu;
            if(potential == "periodic")
                {
                Real thetax = 0.0;
                Real thetay = 0.0;
                if(ax > 0.0)
                    thetax = (jx-1)*pi/ax;
                if(ay > 0.0)
                    thetay = (jy-1)*pi/ay;
                Real cx = cos(thetax);
                Real cy = cos(thetay);
                Real cx2 = cx*cx;
                Real cy2 = cy*cy;
                V(jy,jx) += V0*(cx2+cy2);
                }
            else if(potential == "harmonic")
                {
                for(int kx = 1; kx <= natomx; kx++)
                    {
                    for(int ky = 1; ky <= natomy; ky++)
                        {
                        Real delx = (Real)jx - ax*((Real)kx-0.5);
                        Real dely = (Real)jy - ay*((Real)ky-0.5);
                        V(jy,jx) += V0*(delx*delx+dely*dely);
                        }
                    }
                }
            else if(potential == "pseudoatom")
                {
                for(int kx = 1; kx <= natomx; kx++)
                    {
                    for(int ky = 1; ky <= natomy; ky++)
                        {
                        Real delx02 = 0.01*sqrt(ax*ax+ay*ay);
                        Real delx = (Real)jx - ax*((Real)kx-0.5);
                        Real dely = (Real)jy - ay*((Real)ky-0.5);
                        Real delx2 = delx*delx + dely*dely;
                        V(jy,jx) -= V0/sqrt(1.0 + delx2/delx02);
                        }
                    }
                }
            }
        }
    return V;
    }

void
makeMPO(vector<ITensor>& H, vector<Index> s, vector<Index>& hl, Matrix V, Real phi)
    {
    int nx = V.Ncols();
    int ny = V.Nrows();
    int k = 2*(ny+1);
    int n = nx*ny;
    
    Real t = 1.0;

    for(int j = 0; j <= n; j++)
        hl[j] = Index(nameint("hl",j),k,Link);
    
    int j = 1;
    for(int jx = 1; jx <= nx; ++jx)
        {
        for(int jy = 1; jy <= ny; ++jy)
            {
            int j = (jx - 1)*ny + jy;
            ITensor& W = H[j];
            Index& row = hl[j-1];
            Index& col = hl[j];

            W = ITensor(s[j],prime(s[j]),row,col);
            
            W += op(s[j],"Id") * row(1) * col(1);
          
            W += V(jy,jx) * op(s[j],"n") * row(k) * col(1);

            W += op(s[j],"aT") * row(2) * col(1);
            W += op(s[j],"a") * row(k/2+1) * col(1);
           
            for(int b = 1; b <= ny-1; b++)
                {
                W -= (op(s[j],"Id")-2.0*op(s[j],"n")) * row(b+2) * col(b+1);
                W -= (op(s[j],"Id")-2.0*op(s[j],"n")) * row(k/2+b+1) * col(k/2+b); 
                }

            // nearest neighbor
            if(jy < ny)
                {
                W -= 0.5*t*cos(2.0*(jx-1)*pi*phi) * op(s[j],"a") * row(k) * col(2);
                W -= Complex_i*0.5*t*sin(2.0*(jx-1)*pi*phi) * op(s[j],"a") * row(k) * col(2);
                W -= 0.5*t*cos(2.0*(jx-1)*pi*phi) * op(s[j],"aT") * row(k) * col(k/2+1);
                W += Complex_i*0.5*t*sin(2.0*(jx-1)*pi*phi) * op(s[j],"aT") * row(k) * col(k/2+1);
                }

            // next nearest neighbor
            if(jx < nx)
                {
                W -= 0.5*t * op(s[j],"a") * row(k) * col(k/2);
                W -= 0.5*t * op(s[j],"aT") * row(k) * col(k-1);
                }

            W += op(s[j],"Id") * row(k) * col(k);
            }
        }

    ITensor vL(hl[0]);
    ITensor vR(hl[n]);

    vL(hl[0](k)) = 1.0;
    vR(hl[n](1)) = 1.0;

    H[1] *= vL;
    H[n] *= vR;
    }

Complex
getEnergy(vector<ITensor> H, vector<ITensor>& A, int& c, ITensor& D, vector<Index> s, vector<Index>& l)
    {
    int n = s.size() - 1;
    ITensor En;
    shiftCenter(A,c,D,s,l,1-c);
    En = H[1]*A[1]*D*dag(prime(A[1]))*prime(D);
    for(int j = 2; j <= n; j++)
        En = En*H[j]*A[j]*dag(prime(A[j]));
    return En.toComplex();
    }

Complex
innerProd(vector<ITensor> A1, int& c1, ITensor& D1, vector<ITensor>& A2, int& c2, ITensor& D2, vector<Index> s, vector<Index>& l1, vector<Index>& l2)
    {
    int n = s.size() - 1;
    ITensor prod;
    shiftCenter(A1,c1,D1,s,l1,1-c1);
    shiftCenter(A2,c2,D2,s,l2,1-c2);
    prod = A1[1]*D1*dag(prime(A2[1],Link))*prime(D2);
    for(int j = 2; j <= n; j++)
        prod = prod*A1[j]*dag(prime(A2[j],Link));
    return prod.toComplex();
    }

Real
DMRG(vector<ITensor>& H, vector<ITensor>& A, int&c, ITensor& D, vector<Index>& s, vector<Index>& l, int nsweeps, string verbose)
    {
    cout << "Begin DMRG" << endl;
    cout << endl;

    shiftCenter(A,c,D,s,l,1-c);
    
    int N = s.size()-1;

    Real energy = NAN;
    Real S = NAN;
    Real Vtot = NAN;
    Vector V;
    int M;

    LocalOp<ITensor> Heff;

    ITensor Atemp;
    Index linktemp;

    std::vector<ITensor> Hblock (N+2);

    // Initialize blocks

    Hblock[N] = dag(prime(A[N])) * H[N] * A[N];

    for(int j = N-1; j > 2; --j)
        {
        Hblock[j] = Hblock[j+1] * dag(prime(A[j])) * H[j] * A[j]; 
        }

    int dir = 0;
    int n = 0;
    Real totaltime = 0.0;

    for(int sw = 1; sw <= nsweeps; ++sw)
        {
        // Keep track of time
        clock_t t;
        t = clock();
        for(int k = 1; k < N; ++k)
            {

            if(dir==0) { n = k; }
            else { n = N - k; }

            Heff.update(H[n],H[n+1],Hblock[n-1],Hblock[n+2]);
            ITensor phi = A[n]*D*A[n+1];
            energy = davidson(Heff,phi);
            
            OptSet opts;
            opts.add("Cutoff",err);
            opts.add("Maxm",maxm);
            
            svd(phi,A[n],D,A[n+1],opts);

            // Get entanglement information
            if(n == N/2)
                {
                V = D.diag();
                M = V.Length();
                Vtot = 0.0;
                S = 0.0;
                for(int k = 1; k <= M; ++k)
                    {
                    Real lamk = V(k)*V(k);
                    Vtot += lamk;
                    S -= lamk*log2(lamk);
                    }
                }

            l[n] = commonIndex(A[n],D);
            l[n+1] = commonIndex(D,A[n+1]);

            if(dir==0)
                {
                if(n==1) 
                    { 
                    Hblock[1] = H[1] * A[1] * dag(prime(A[1])); 
                    }
                else 
                    { 
                    Hblock[n] = Hblock[n-1] * H[n] * A[n] * dag(prime(A[n])); 
                    }
                if(n < N-1)
                    {
                    shiftCenter(A,c,D,s,l,1);
                    }
                }
            else
                {
                if(n==N-1) 
                    {
                    Hblock[N] = dag(prime(A[N])) * H[N] * A[N];
                    }
                else 
                    {
                    Hblock[n+1] = Hblock[n+2] * dag(prime(A[n+1])) * H[n+1] * A[n+1];
                    }
                if(n > 1)
                    {
                    shiftCenter(A,c,D,s,l,-1);
                    }
                }
            }
        t = clock() - t;
        Real time = ((Real)t)/CLOCKS_PER_SEC;
        totaltime += time;
        if(verbose == "verbose")
            {
            cout << "Sweep " << sw << endl;
            cout << "Energy = " << energy << endl;
            cout << "States kept = " << M << endl;
            cout << "Smallest singular value = " << V[M-1] << endl;
            cout << "Vtot = " << Vtot << endl;
            cout << "S = " << S << endl;
            cout << "time (s) = " << time << endl;
            cout << "log(time) = " << log10(time) << endl;
            cout << endl;
            }
        if(dir==0) { dir = 1; }
        else { dir = 0; }
        }
    cout << "DMRG complete" << endl;
    cout << "Final information:" << endl;
    cout << "Sweep " << nsweeps << endl;
    cout << "Energy = " << energy << endl;
    cout << "States kept = " << M << endl;
    cout << "Smallest singular value = " << V[M-1] << endl;
    cout << "Vtot = " << Vtot << endl;
    cout << "S = " << S << endl;
    cout << "Total time = " << totaltime << endl;
    cout << endl;

    return energy;
    }

void
makeHam(Matrix& Hre, Matrix& Him, Matrix V, Real t, Real phi)
    {
    int nx = V.Ncols();
    int ny = V.Nrows();
    int n = nx*ny;

    // Create hopping Hamiltonian
    Hre = 0.0;
    Him = 0.0;
    for(int jx = 1; jx <= nx; jx++)
	    {
        for(int jy = 1; jy <= ny; jy++)
            {
            int j = (jx-1)*ny + jy;
            Hre(j,j) += V(jy,jx);
            if(jy < ny && j < n)
                {
                Hre(j,j+1) = Hre(j+1,j) = -0.5*t*cos(2.0*(jx-1)*pi*phi);
                Him(j,j+1) = -0.5*t*sin(2.0*(jx-1)*pi*phi);
                Him(j+1,j) = -Him(j,j+1);
                }
            if(jx < nx)
                Hre(j,j+ny) = Hre(j+ny,j) = -0.5*t;
            }
        }
    }

void
MPOGS(const Matrix& H, vector<ITensor>& V, int&c, ITensor& D, vector<Index>& s, vector<Index>& l)
    {
    cout << "MPO method" << endl;
    int n = H.Nrows();
    cout << "H = " << endl;
    cout << H;
    Vector eps; 
    Matrix U;
    EigenValues(H,eps,U);
    cout << "eps = " << endl;
    cout << eps;
    cout << "U = " << endl;
    cout << U;
    
    int nF = 0;
    for(int j = 1; j <= n; j++)
        if(eps(j) < 0.0) { nF++; }
    cout << "nF = " << nF;
    cout << endl;
    
    // Start with empty product state
    for(int j = 1; j <= n; j++)
        l[j] = Index(nameint("l",j),1,Link);

    V[1] = ITensor(s[1],l[1]);
    V[1](s[1](1),l[1](1)) = 1.0;
    D = ITensor(l[1],l[2],1.0);

    for(int j = 2; j < n; j++)
        {
        V[j] = ITensor(s[j],l[j],l[j+1]);
        V[j](s[j](1),l[j](1),l[j+1](1)) = 1.0;
        }
    
    V[n] = ITensor(s[n],l[n]);
    V[n](s[n](1),l[n](1)) = 1.0;
    
    PrintDat(D);
    PrintDat(V[n]);

    for(int j = 1; j <= nF; j++)
        {
        Vector Uj = U.Column(j);
        cout << Uj;

        vector<ITensor> A(n+1);
        vector<Index> lA(n+1);
        }

    }

void
halfK(Matrix& phi, Real& w, Real& O, Matrix& invO_matrix_up, Matrix& invO_matrix_dn, Matrix Proj_k_half, Matrix Phi_T, int N_up, int N_par)
    {
    phi = Proj_k_half*phi;
    
    int N_site = phi.Nrows();
   
    // Edit (overlap)

    Matrix Phi_T_up;
    Matrix phi_up;
    Matrix Phi_T_dn;
    Matrix phi_dn;
    Real detinvO_matrix_up = 1.0;
    Real detinvO_matrix_dn = 1.0;

    if(N_up > 0)
        {
        Phi_T_up = Phi_T.SubMatrix(1,N_site,1,N_up);
        phi_up = phi.SubMatrix(1,N_site,1,N_up);
        invO_matrix_up = Inverse(Phi_T_up.t()*phi_up);
        cout << invO_matrix_up;
        detinvO_matrix_up = Determinant(invO_matrix_up);
        }
    
    if(N_par-N_up > 0)
        {
        Phi_T_dn = Phi_T.SubMatrix(1,N_site,N_up+1,N_par);
        phi_dn = phi.SubMatrix(1,N_site,N_up+1,N_par);
        invO_matrix_dn = Inverse(Phi_T_dn.t()*phi_dn);
        detinvO_matrix_dn = Determinant(invO_matrix_dn);
        }
    
    Real O_new = 1.0/(detinvO_matrix_up*detinvO_matrix_dn);
    Real O_ratio=O_new/O;

    if(O_ratio > 0)
        {
        O=O_new;
        w=w*O_ratio;
        }
    else
        {
        w = 0;
        }

    }

void
V(Vector& phi, Vector phi_T, int N_up, int N_par, Real& O, Real& w, Matrix& invO_matrix_up, Matrix& invO_matrix_dn, Matrix aux_fld)
    {
    Vector Gii(2);
    Matrix RR(2,2);
    Gii = 0.0;
    RR = 0.0;

    // Edit (overlap)

    Vector temp1_up(N_up);
    Vector temp2_up(N_up);
    Vector temp1_dn(N_par-N_up);
    Vector temp2_dn(N_par-N_up);
    temp1_up = 0.0;
    temp2_up = 0.0;
    temp1_dn = 0.0;
    temp2_dn = 0.0;

    if(N_up > 0)
        {
        temp1_up = phi.SubVector(1,N_up)*invO_matrix_up;
        temp2_up = invO_matrix_up*phi_T.SubVector(1,N_up);
        Gii(1) = temp1_up*phi_T.SubVector(1,N_up);
        }

    if(N_par-N_up > 0)
        {
        temp1_dn = phi.SubVector(N_up+1,N_par)*invO_matrix_dn;
        temp2_dn = invO_matrix_dn*phi_T.SubVector(N_up+1,N_par);
        Gii(2) = temp1_dn*phi_T.SubVector(N_up+1,N_par);
        }

    RR(1,1) = (aux_fld(1,1)-1.0)*Gii(1) + 1.0;
    RR(1,2) = (aux_fld(1,2)-1.0)*Gii(1) + 1.0;
    RR(2,1) = (aux_fld(2,1)-1.0)*Gii(2) + 1.0;
    RR(2,2) = (aux_fld(2,2)-1.0)*Gii(2) + 1.0;

    Vector O_ratio_temp(2);
    O_ratio_temp(1) = RR(1,1)*RR(2,1);
    O_ratio_temp(2) = RR(1,2)*RR(2,2);

    Vector O_ratio_temp_real(2);
    O_ratio_temp_real = 0.0;
    if(O_ratio_temp(1) > 0.0)
        O_ratio_temp_real(1) = O_ratio_temp(1);
    if(O_ratio_temp(2) > 0.0)
        O_ratio_temp_real(2) = O_ratio_temp(2);
    Real sum_O_ratio_temp_real = O_ratio_temp_real(1)+O_ratio_temp_real(2);

    if(sum_O_ratio_temp_real <= 0)
        w=0;

    if(w > 0)
        {
        w=w*0.5*sum_O_ratio_temp_real;
        
        int x_spin;
        if(O_ratio_temp_real(1)/sum_O_ratio_temp_real >= rNum())
            x_spin=1;
        else
            x_spin=2;
        if(N_up > 0)
            phi.SubVector(1,N_up)=phi.SubVector(1,N_up)*aux_fld(1,x_spin);
        if(N_par-N_up > 0)
            phi.SubVector(N_up+1,N_par)=phi.SubVector(N_up+1,N_par)*aux_fld(2,x_spin);
    
        O=O*O_ratio_temp(x_spin);
        invO_matrix_up=invO_matrix_up+(1-aux_fld(1,x_spin))/RR(1,x_spin)*temp2_up*temp1_up;
        invO_matrix_dn=invO_matrix_dn+(1-aux_fld(2,x_spin))/RR(2,x_spin)*temp2_dn*temp1_dn;
        }
    }

Real
measure(Matrix H_k, Matrix phi, Matrix Phi_T, Matrix invO_matrix_up, Matrix invO_matrix_dn, int N_up, int N_par, Real U)
    {
    Real e = 0.0;
    int N_sites = phi.Nrows();

    // Edit (energy)

    Vector diag_G_up(N_sites);
    diag_G_up = 0.0;
    Vector diag_G_dn(N_sites);
    diag_G_dn = 0.0;
    Matrix G_up(N_sites,N_sites);
    Matrix G_dn(N_sites,N_sites);
    G_up = 0.0;
    G_dn = 0.0;

    if(N_up > 0)
        {
        Matrix temp_up = phi.SubMatrix(1,N_sites,1,N_up)*invO_matrix_up;
        G_up = temp_up*(Phi_T.SubMatrix(1,N_sites,1,N_up)).t();
        diag_G_up = G_up.Diagonal();
        }

    if(N_par - N_up > 0)
        {
        Matrix temp_dn = phi.SubMatrix(1,N_sites,N_up+1,N_par)*invO_matrix_dn;
        G_dn = temp_dn*(Phi_T.SubMatrix(1,N_sites,N_up+1,N_par)).t();
        diag_G_dn = G_dn.Diagonal();
        }
   
    Real n_int = diag_G_up*diag_G_dn;
    Real potentialEnergy = n_int*U;
    
    Real kineticEnergy = 0.0;
    for(int i = 1; i <= N_sites; i++)
        {
        for(int j = 1; j <= N_sites; j++)
            {
            kineticEnergy += H_k(i,j)*(G_up(i,j)+G_dn(i,j));
            }
        }
    
    e = potentialEnergy + kineticEnergy;

    return e;
    }

void
stepwlk(vector<Matrix>& phi, int N_wlk, int N_sites, Vector& w, Vector& O, Real& E, Real& W, Matrix H_k, Matrix Proj_k_half, int flag_mea, Matrix Phi_T, int N_up, int N_par, Real U, Real fac_norm, Matrix aux_fld)
    {
    Vector e(N_wlk);
    e = 0.0;
    Matrix invO_matrix_up;
    Matrix invO_matrix_dn;

    for(int i_wlk = 1; i_wlk <= N_wlk; i_wlk++)
        {
        Matrix Phi = phi[i_wlk];
        if(w(i_wlk) > 0)
            {
            w(i_wlk) = w(i_wlk)*exp(fac_norm);
            halfK(Phi, w(i_wlk), O(i_wlk), invO_matrix_up, invO_matrix_dn, Proj_k_half, Phi_T, N_up, N_par);
            if(w(i_wlk) > 0)
                {
                for(int j_site = 1; j_site <= N_sites; j_site++)
                    {
                    if(w(i_wlk) > 0)
                        {
                        Vector Phi_site = Phi.Row(j_site);
                        V(Phi_site, Phi_T.Row(j_site), N_up, N_par, O(i_wlk), w(i_wlk), invO_matrix_up, invO_matrix_dn, aux_fld);
                        Phi.Row(j_site) = Phi_site;
                        }
                    }
                }
            if(w(i_wlk) > 0)
                {
                halfK(Phi, w(i_wlk), O(i_wlk), invO_matrix_up, invO_matrix_dn, Proj_k_half, Phi_T, N_up, N_par);
                if(w(i_wlk) > 0)
                    {
                    if(flag_mea == 1)
                        {
                        e(i_wlk) = measure(H_k, Phi, Phi_T,  invO_matrix_up, invO_matrix_dn, N_up, N_par, U); 
                        }
                    }
                }
            }
        phi[i_wlk] = Phi;
        }

    // Compute the ensemble's total energy and weight if measurement took place
    if(flag_mea == 1)
        {
        for(int i_wlk = 1; i_wlk <= N_wlk; i_wlk++)
            {
            if(w(i_wlk) > 0)
                {
                E = E + e(i_wlk)*w(i_wlk);
                W = W + w(i_wlk);
                }
            }
        }
    }

void
stblz(vector<Matrix>& Phi, int N_wlk, Vector& O, int N_up, int N_par)
    {
    int N_sites = Phi[1].Nrows();
    int N_dn = N_par - N_up;

    for(int i_wlk = 1; i_wlk <= N_wlk; i_wlk++)
        {
        Matrix phi = Phi[i_wlk];
        Real det_R_up = 1.0;
        Real det_R_dn = 1.0;

        if(N_up > 0)
            {
            Matrix Q_up;
            Matrix R_up;
            Matrix phi_up(N_sites,N_sites);
            phi_up = 0.0;
            phi_up.SubMatrix(1,N_sites,1,N_up) = phi.SubMatrix(1,N_sites,1,N_up);
            QRDecomp(phi_up,Q_up,R_up);
            phi.SubMatrix(1,N_sites,1,N_up) = Q_up.SubMatrix(1,N_sites,1,N_up);
            det_R_up = Determinant(R_up.SubMatrix(1,N_up,1,N_up));
            }
        if(N_dn > 0)
            {
            Matrix Q_dn;
            Matrix R_dn;
            Matrix phi_dn(N_sites,N_sites);
            phi_dn = 0.0;
            phi_dn.SubMatrix(1,N_sites,1,N_dn) = phi.SubMatrix(1,N_sites,N_up+1,N_par);
            QRDecomp(phi_dn,Q_dn,R_dn);
            phi.SubMatrix(1,N_sites,N_up+1,N_par) = Q_dn.SubMatrix(1,N_sites,1,N_dn);
            det_R_dn = Determinant(R_dn.SubMatrix(1,N_dn,1,N_dn));
            }
        Phi[i_wlk] = phi;
        O(i_wlk) = O(i_wlk)/det_R_up/det_R_dn;
        }
    }

void
pop_cntrl(vector<Matrix>& Phi, Vector& w, Vector& O, int N_wlk, int N_sites, int N_par)
    {
    vector<Matrix> new_Phi(N_wlk+1);
    Matrix zeros(N_sites, N_par);
    zeros = 0.0;
    for(int i = 1; i <= N_wlk; i++)
        new_Phi[i] = zeros;

    Vector new_O(N_wlk);
    new_O = 0.0;
    
    Real sum_w = 0.0;
    for(int i = 1; i <= N_wlk; i++)
        sum_w += w(i);
    Real d = 1.0*N_wlk/sum_w;
    
    sum_w = -rNum();
    int n_wlk = 0;

    for(int i_wlk = 1; i_wlk <= N_wlk; i_wlk++)
        {
        sum_w += w(i_wlk)*d;
        int n = ceil(sum_w);
        for(int j = n_wlk+1; j <= n; j++)
            {
            new_Phi[j] = Phi[i_wlk];
            new_O(j) = O(i_wlk);
            }
        n_wlk = n;
        }

    Phi = new_Phi;
    O = new_O;

    w = Vector(N_wlk,1.0);
    }

int main(int argc, char* argv[])
    {
   
    cout << setprecision(16);


    int Lx = 4;
    int Ly = 1;

    int N_up = 2;
    int N_dn = 0;

    Real U = 0.0;
    Real tx = 1.0;
    Real ty = 1.0;

    Real deltau = 0.01;
    int N_wlk = 1;
    int N_blksteps = 1;
    int N_eqblk = 1;
    int N_blk = 1;
    int itv_modsvd = 1;
    int itv_pc = 1;
    int itv_Em = 1;
    
    int N_sites = Lx*Ly;
    int N_par = N_up + N_dn;

    Matrix H_k(N_sites,N_sites);
    H_k = 0.0;

    int r = 0;
    for(int iy = 1; iy <= Ly; iy++)
        {
        for(int jx = 1; jx <= Lx; jx++)
            {
            r++;
            if(Lx != 1)
                {
                if(jx == 1)
                    {
                    H_k(r,r+1) = H_k(r,r+1) - tx;
                    }
                else if(jx == Lx)
                    {
                    H_k(r,r-1) = H_k(r,r-1) - tx;
                    }
                else
                    {
                    H_k(r,r-1) = -tx;
                    H_k(r,r+1) = -tx;
                    }
                }

            if(Ly != 1)
                {
                if(iy == 1)
                    {
                    H_k(r,r+Lx) = H_k(r,r-Lx) - ty;
                    }
                else if(iy == Ly)
                    {
                    H_k(r,r-Lx) = H_k(r,r-Lx) - ty;
                    }
                else
                    {
                    H_k(r,r-Lx) = -ty;
                    H_k(r,r-Lx) = -ty;
                    }
                }
            }
        }
   
    Matrix Proj_k_half = Exp(-0.5*deltau*H_k);

    int n = N_sites;

    Vector E_nonint;
    Matrix psi_nonint;

    EigenValues(H_k, E_nonint, psi_nonint);

    Matrix Phi_T(N_sites,N_par);
    if(N_up > 0)
        {
        Phi_T.SubMatrix(1,N_sites,1,N_up) = psi_nonint.SubMatrix(1,N_sites,1,N_up);
        }
    if(N_dn > 0)
        {
        Phi_T.SubMatrix(1,N_sites,N_up+1,N_par) = psi_nonint.SubMatrix(1,N_sites,1,N_dn);;
        }

    vector<ITensor> V_T(n+1);
    ITensor D_T;
    int c_T = n-1;
    vector<Index> s(n+1);
    vector<Index> l_T(n+1);

    int d = 2;
    for(int j = 1; j <= n; j++)
        s[j] = Index(nameint("s",j),d,Site);

    cout << "Creating myMPS\n" << endl;
    myMPS MPS_T;
    MPS_T.init(s);
    cout << "Sites = " << MPS_T.N() << endl;
    cout << "Center = " << MPS_T.C() << endl;
    cout << "Site 1 = " << MPS_T.S(1) << endl;
    cout << "Link 1 = " << MPS_T.L(1) << endl;
    PrintDat(MPS_T.D());
    PrintDat(MPS_T.V(1));
    MPS_T.shiftCenter(1);
    cout << "Center = " << MPS_T.C() << endl;
    PrintDat(MPS_T.V(2));
    PrintDat(MPS_T.D());

    //PrintDat(MPS_T.V(1)*MPS_T.D()*MPS_T.V(2)*MPS_T.V(3)*MPS_T.V(4));
    PrintDat(MPS_T.getPsi());
    cout << "Norm = " << MPS_T.getNorm() << endl;

    int b = 4;
    if(b > n)
        b = n;
    cout << "b = " << b << endl;
    int expb = 1 << (b-1);
    cout << "Chi_max = 2^(b-1) = " << expb << endl;
    cout << endl;

    GMPS(Phi_T,N_up,N_par,V_T,c_T,D_T,s,l_T,b);
    ITensor Psi = V_T[1]*D_T*V_T[2]*V_T[3]*V_T[4];
    //ITensor Psi = V_T[1]*D_T*V_T[2];
    cout << Psi;

    vector<Index> hl(n+1);
    vector<ITensor> H(n+1);

    Matrix V(Lx,Ly);
    V = 0.0;
    Real phi = 0.0;

    makeMPO(H,s,hl,V,phi);
    
    ITensor Ham = H[1]*H[2]*H[3]*H[4];
    //ITensor Ham = H[1]*H[2];
    cout << 2.0*Psi*Ham*prime(conj(Psi));

    Real E_K = 0.0;

    for(int i = 1; i <= N_up; i++)
        E_K += E_nonint(i);
    for(int i = 1; i <= N_dn; i++)
        E_K += E_nonint(i);

    Vector n_r_up(N_sites);
    Vector n_r_dn(N_sites);
    n_r_up = 0.0;
    n_r_dn = 0.0;

    // Edit (energy)

    if(N_up > 0)
        {
        Matrix Lambda_up = Phi_T.SubMatrix(1,N_sites,1,N_up)*(Phi_T.SubMatrix(1,N_sites,1,N_up)).t();
        n_r_dn = Lambda_up.Diagonal();
        }
    if(N_dn > 0)
        {
        Matrix Lambda_dn = Phi_T.SubMatrix(1,N_sites,N_up+1,N_par)*(Phi_T.SubMatrix(1,N_sites,N_up+1,N_par)).t();
        n_r_dn = Lambda_dn.Diagonal();
        }

    Real E_V = U*n_r_up*n_r_dn;
    Real E_T = E_K + E_V;

    // Assemble the initial population of walkers
    vector<Matrix> Phi(N_wlk+1);
    for(int i = 1; i <= N_wlk; i++)
        {
        Phi[i] = Phi_T;
        }

    Vector w(N_wlk);
    Vector O(N_wlk);

    w = 1.0, O = 1.0;

    Vector E_blk(N_blk);
    Vector W_blk(N_blk);

    E_blk = 0.0, W_blk = 0.0;

    // initialize auxiliary field constants
    Real fac_norm = (E_T-0.5*U*N_par)*deltau;
    Real gamma = acosh(exp(0.5*deltau*U));
    Matrix aux_fld(2,2);
    aux_fld = 0.0;

    for(int i = 1; i <= 2; i++)
        {
        for(int j = 1; j <= 2; j++)
            {
            if((i+j)%2 == 0)
                aux_fld(i,j)=exp(gamma);
            else
                aux_fld(i,j)=exp(-gamma);
            }
        }
   
    srand(time(NULL));

    int flag_mea = 0;
    Real E = 0.0;
    Real W = 0.0;

    // Equilibration phase
    for(int i_blk = 1; i_blk <= N_eqblk; i_blk++)
        {
        for(int j_step = 1; j_step <= N_blksteps; j_step++)
            {
            stepwlk(Phi, N_wlk, N_sites, w, O, E, W, H_k, Proj_k_half, flag_mea, Phi_T, N_up, N_par, U, fac_norm, aux_fld);
            if(j_step % itv_modsvd == 0)
                {
                stblz(Phi, N_wlk, O, N_up, N_par);
                }
            if(j_step % itv_pc == 0)
                {
                pop_cntrl(Phi, w, O, N_wlk, N_sites, N_par);
                }
            }
        }
    
    // Measurement phase
    for(int i_blk = 1; i_blk <= N_blk; i_blk++)
        {
        for(int j_step = 1; j_step <= N_blksteps; j_step++)
            {
            if(j_step % itv_Em == 0)
                flag_mea = 1;
            else
                flag_mea = 0;
            stepwlk(Phi, N_wlk, N_sites, w, O, E_blk(i_blk), W_blk(i_blk), H_k, Proj_k_half, flag_mea, Phi_T, N_up, N_par, U, fac_norm, aux_fld);
            if (j_step % itv_modsvd == 0)
                {
                stblz(Phi, N_wlk, O, N_up, N_par);
                }
            if (j_step % itv_pc == 0)
                pop_cntrl(Phi, w, O, N_wlk, N_sites, N_par);
            if (j_step % itv_Em ==0)
                fac_norm = (E_blk(i_blk)/W_blk(i_blk)-0.5*U*N_par)*deltau;
            }
        E_blk(i_blk)=E_blk(i_blk)/W_blk(i_blk);
        }
    
    Real E_ave = 0.0;
    for(int i = 1; i <= N_blk; i++)
        E_ave += E_blk(i);
    E_ave = E_ave/N_blk;
    Real E_err = 0.0;
    if(N_blk > 1)
        {
        for(int i = 1; i <= N_blk; i++)
            E_err += (E_blk(i) - E_ave)*(E_blk(i) - E_ave);
        E_err = sqrt(E_err) / sqrt(N_blk - 1);
        E_err = E_err / sqrt(N_blk);
        }

    cout << "E_ave = " << E_ave << endl;
    cout << "E_err = " << E_err << endl;
    cout << endl;
    
    /*
    // Number of sites
    int mx = 2;
    int my = 0;
    int m = mx+my;
    int nx = 1 << mx;
    int ny = 1 << my;
    int n = 1 << m;
    cout << "mx = " << mx << endl;
    cout << "my = " << my << endl;
    cout << "m = mx + my = " << m << endl;
    cout << "nx = 2^mx = " << nx << endl;
    cout << "ny = 2^my = " << ny << endl;
    cout << "n = 2^m = " << n << endl;
    cout << endl;

    Real mu = 0.0,
         t = 1.0,
         Delta = 0.0,
         phi = 0.0,
         V0 = 0.0;
    int natomx = 1,
        natomy = 1;

    //string potential = "pseudoatom";
    string potential = "periodic";

    Matrix V = makePotential(nx,ny,natomx,natomy,mu,V0,potential);
    //cout << V << endl;

    Matrix Hre(n,n),
           Him(n,n),
           Gre(n,n),
           Gim(n,n);
    Hre = 0, Him = 0, Gre = 0, Gim = 0;
    
    makeHam(Hre,Him,V,t,phi);
    
    vector<ITensor> VB(n+1);
    ITensor D;
    int c = n-1;
    vector<Index> s(n+1);
    vector<Index> l(n+1);

    int d = 2;
    for(int j = 1; j <= n; j++)
        s[j] = Index(nameint("s",j),d,Site);

    int b = 3;
    cout << "b = " << b << endl;
    int expb = 1 << (b-1);
    cout << "Chi_max = 2^(b-1) = " << expb << endl;
    cout << endl;

    Real Egs = GMPS(Hre,Him,Gre,Gim,VB,c,D,s,l,b);

    int L = n/2;
    shiftCenter(VB, c, D, s, l, L-c);
    Vector DL = D.diag();
    int chiL = DL.Length();
    Real normL = 0.0;
    Real SL = 0.0;
    for(int j = 1; j <= chiL; j++)
        {
        Real lambdaj = DL(j)*DL(j);
        normL += lambdaj;
        SL -= lambdaj*log2(lambdaj);
        }

    vector<Index> hl(n+1);
    vector<ITensor> H(n+1);

    makeMPO(H,s,hl,V,phi);

    Real EGMPS = getEnergy(H,VB,c,D,s,l).real();
    Real delE = fabs(EGMPS-Egs)/fabs(Egs);

    cout << "Chi = " << chiL << endl;
    cout << "EGMPS = " << EGMPS << endl;
    cout << "|Del E / E| = " << delE << endl;
    cout << "log|delE / E| = " << log10(delE) << endl;
    cout << "1-norm = " << 1.0-normL << endl;
    cout << "S_L = " << SL << endl;
    cout << endl;

    vector<ITensor> VBP(n+1);
    for(int j = 1; j <= n; j++)
        {
        VBP[j] = VB[j];
        //VBP[j].randomize();
        }
    ITensor DP = D;
    int cP = c;
    vector<Index> lP(n+1);
    for(int j = 1; j <= n; j++)
        lP[j] = l[j];
    //shiftCenter(VBP,cP,DP,s,lP,n-1-cP);
    //shiftCenter(VBP,cP,DP,s,lP,1-cP);
    
    int nsweeps = 2;
    Real EDMRG = DMRG(H,VBP,cP,DP,s,lP,nsweeps,"verbose");
    Real delEDMRG = fabs(EDMRG-Egs)/fabs(Egs); 

    cout << "EDMRG = " << EDMRG << endl;
    cout << "|Del E / E| = " << delEDMRG << endl;
    cout << "log|delE / E| = " << log10(delEDMRG) << endl;
    */

    }


