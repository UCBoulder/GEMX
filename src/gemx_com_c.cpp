//ACTUAL GEMX_COM FILE IN C
#include "gemx_com_c.hpp"
//Array Declarations
CArray2D<int> ileft;
CArray2D<int> jleft;
CArray2D<int> iright;
CArray2D<int> jright;

CArray2D<double> ke;
CArray2D<double> nos;

CArray2D<double> efle; //1st index 1-based
CArray2D<double> pfle; //1st index 1-based
CArray2D<double> pfl; //1st index 1-based, but +1 in ftn
CArray2D<double> efl;   //1st index 1-based

CArray3D<double> phi_k;
CArray3D<double> dphidr;
CArray3D<double> dphi_kdr;
CArray3D<double> d2phidr2;
CArray3D<double> d2phi_kdr2;
CArray3D<double> d2phidz2;
CArray3D<double> d2phi_kdz2;
CArray3D<double> dphidz;
CArray3D<double> dphi_kdz;
CArray3D<double> OPPphi;
CArray3D<double> OPPphik;
CArray3D<double> l_hand;
CArray3D<double> r_hand;

CArray2D<double> den2d1; 
CArray2D<double> den2d2; 
CArray2D<double> dden2d;

CArray2D<double> bmag;
CArray2D<double> gbtor;
CArray2D<double> gbx;
CArray2D<double> gbz;

CArray2D<double> gn0i;
CArray2D<double> gn0e;
CArray2D<double> gt0i;
CArray2D<double> gt0e;
CArray2D<double> xforw;
CArray2D<double> zforw;
CArray2D<double> zbackw;
CArray2D<double> xbackw;

CArray2D<double> gcpnex;
CArray2D<double> gcpnez;
CArray2D<double> gcptex;
CArray2D<double> gcptez;

CArray2D<double> gnuobx;
CArray2D<double> gnuoby;
CArray2D<double> gupae0;

CArray3D<double> rho;
CArray3D<double> phi;//!,den_pre(0:imx,0:jmx,0:kmx),dden(0:imx,0:jmx,0:kmx))

CArray3D<double> ex;
CArray3D<double> ez; //note: +1 since arrays are 0-imx with imx included in ftn
CArray3D<double> ezeta;

CArray3D<double> delbx; 
CArray3D<double> delbz; 
CArray3D<double> delby;

CArray3D<double> apar;
CArray3D<double> dene;

CArray3D<double> upar;
CArray3D<double> jpar;
CArray3D<double> upars;
CArray3D<double> phis;
CArray3D<double> denes;
CArray3D<double> apars;

CArray3D<double> dnedx;
CArray3D<double> dnedy;
CArray3D<double> dupadx;
CArray3D<double> dupady;

CArray3D<double> rk_hand;

MPI_Comm TUBE_COMM, GRID_COMM, PETSC_COMM;
int imx, jmx, kmx, mmx;
int numprocs;
int last,myid, cnt , ierr;
int PADE, CST, weightscheme,modes,filtering_iterations, cold_start,checkpoint;

int nmx,nsmx,nsubd=8,ntube=4,petsc_color,petsc_rank,iBoltzmann,globle_integer=0,eBoltzmann,eAdiabatic,iterations,dbg;
int rand_table[10007];
    char outname[71];
    double endtm,begtm,pstm;
    double starttm, lasttm, tottm;
      double start_total_tm, end_total_tm, start_integ_tm, end_integ_tm, start_ppush_tm, end_ppush_tm, start_cpush_tm, end_cpush_tm;
      double total_tm = 0.0, integ_tm = 0.0, ppush_tm = 0.0, cpush_tm = 0.0;
//      imx,jmx,kmx = max no. of grid pts in x,y,z
// !    mmx         = max no. of particles
// !    nmx         = max. no. of time steps
// !    nsmx        = max. no. of species (including tracer particles
int* mm = nullptr; int *tmm = nullptr; int *lr = nullptr;
double *mims = nullptr; double *q = nullptr;
int timestep, iez;
int iseed;
//double *time = nullptr; //causing issues, don't want to deal with it

double dx,dz,dzeta,dt,totvol,n0,tcurr;
double etaohm;
double lx,lz;
int nm,nsm,ncurr,iflr,ifield_solver,ntracer,i3D,icollision;
double cut,amp,tor,amie,emass,qel,rneu;
int iput,iget,ision,isham,peritr,iadi;
int idg;
const int master = 0;

//    Various diagnostic arrays and scalars
//    plotting constants

    int nplot, xnplt;

double vcut;
int nonlin,nonline,iflut,ifluid,ipara;
std::complex<double> IU;

double *xg = nullptr;
double *zg = nullptr;
double *jac = nullptr;

//          particle array declarations
double *mu = nullptr; //note: be careful here - these are 1-based index in c++, need account for in code
double *x2 = nullptr;
double *zeta2 = nullptr;
double *z2 = nullptr;
double *u2 = nullptr;
double *x3 = nullptr;
double *zeta3 = nullptr;
double *z3 = nullptr;
double *u3 = nullptr;
double *w2 = nullptr;
double *w3 = nullptr;
double *gw = nullptr;

double *fe = nullptr;
double *te = nullptr;
double *rmsphi = nullptr;
double *rmsez = nullptr;
double *rmsapa = nullptr;
double *avewi = nullptr;
double *vol = nullptr;

CArray4D<double> den;



//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void new_gemx_com(){
    mm = new int[nsmx+1]; std::fill(mm, mm+nsmx+1, 0);
    tmm = new int[nsmx+1]; std::fill(tmm, tmm+nsmx, 0);
    lr = new int[nsmx+1]; std::fill(lr, lr+nsmx+1, 0);
    mims = new double[nsmx+1]; std::fill(mims, mims+nsmx+1, 0.0);
    q = new double[nsmx]; std::fill(q, q+nsmx, 0.0);

    //time = new double[nmx+1];

    xg = new double[imx+1]; std::fill(xg, xg+imx+1, 0.0);
    zg = new double[jmx+1]; std::fill(zg, zg+jmx+1, 0.0);
    jac = new double[imx+1]; std::fill(jac, jac+imx+1, 0.0);

    //          particle array declarations
    mu = new double[mmx]; std::fill(mu, mu+mmx, 0.0);//note: be careful here - these are 1-based index in c++, need account for in code
    x2 = new double[mmx]; std::fill(x2, x2+mmx, 0.0);
    zeta2 = new double[mmx]; std::fill(zeta2, zeta2+mmx, 0.0);
    z2 = new double[mmx]; std::fill(z2, z2+mmx, 0.0);
    u2 = new double[mmx]; std::fill(u2, u2+mmx, 0.0);
    x3 = new double[mmx]; std::fill(x3, x3+mmx, 0.0);
    zeta3 = new double[mmx]; std::fill(zeta3, zeta3+mmx, 0.0);
    z3 = new double[mmx]; std::fill(z3, z3+mmx, 0.0);
    u3 = new double[mmx]; std::fill(u3, u3+mmx, 0.0);
    w2 = new double[mmx]; std::fill(w2, w2+mmx, 0.0);
    w3 = new double[mmx]; std::fill(w3, w3+mmx, 0.0);
    gw = new double[mmx]; std::fill(gw,gw+mmx, 0.0);

//      variables for tracing a grid (i,j) along field lien to the neighboring planes
    ileft.resize(imx+1, jmx+1); 
    jleft.resize(imx+1, jmx+1);
    iright.resize(imx+1, jmx+1);
    jright.resize(imx+1, jmx+1);

//  energy diagnostic arrays
    ke.resize(nsmx, nmx+1); //careful - 1st index 1 based in ftn (next one too)
    fe = new double[nmx+1]; std::fill(fe, fe+nmx+1, 0.0);
    te = new double[nmx+1]; std::fill(te, te+nmx+1, 0.0);
    rmsphi = new double[nmx+1]; std::fill(rmsphi, rmsphi+nmx+1, 0.0);
    rmsez = new double[nmx+1]; std::fill(rmsez, rmsez+nmx+1, 0.0);
    rmsapa = new double[nmx+1]; std::fill(rmsapa, rmsapa+nmx+1, 0.0);
    avewi = new double[nmx+1]; std::fill(avewi, avewi+nmx+1, 0.0);
    nos.resize(nsmx, nmx+1); 

    //  flux diagnostics
    vol = new double[nsubd]; std::fill(vol, vol+nsubd, 0.0);//careful, 1-indexed
    efle.resize(nsubd, nmx+1); //1st index 1-based
    pfle.resize(nsubd, nmx+1); //1st index 1-based
    pfl.resize(nsmx+1, nmx+1); //1st index 1-based, but +1 in ftn
    efl.resize(nsmx, nmx+1);   //1st index 1-based

    
    
        
    //Boltzmann Electron subroutine arrays for Newton solve
    phi_k.resize(imx+1, jmx+1, kmx+1);
    dphidr.resize(imx+1, jmx+1, kmx+1);
    dphi_kdr.resize(imx+1, jmx+1, kmx+1);
    d2phidr2.resize(imx+1, jmx+1, kmx+1);
    d2phi_kdr2.resize(imx+1, jmx+1, kmx+1);
    d2phidz2.resize(imx+1, jmx+1, kmx+1);
    d2phi_kdz2.resize(imx+1, jmx+1, kmx+1);
    dphidz.resize(imx+1, jmx+1, kmx+1);
    dphi_kdz.resize(imx+1, kmx+1, kmx+1);
    OPPphi.resize(imx+1, jmx+1, kmx+1);
    OPPphik.resize(imx+1, jmx+1, kmx+1);
    l_hand.resize(imx+1, jmx+1, kmx+1);
    r_hand.resize(imx+1, jmx+1, kmx+1);

    /*The following arrays didn't have a label*/

    //2D Arrays
    den2d1.resize(imx+1, jmx+1); 
    den2d2.resize(imx+1, jmx+1); 
    dden2d.resize(imx+1, jmx+1);

    bmag.resize(imx+1, jmx+1);
    gbtor.resize(imx+1, jmx+1);
    gbx.resize(imx+1, jmx+1);
    gbz.resize(imx+1, jmx+1);

    gn0i.resize(imx+1, jmx+1);
    gn0e.resize(imx+1, jmx+1);
    gt0i.resize(imx +1, jmx+1);
    gt0e.resize(imx+1, jmx+1);
    xforw.resize(imx+1, jmx+1);
    zforw.resize(imx+1, jmx+1);
    zbackw.resize(imx+1, jmx+1);
    xbackw.resize(imx+1, jmx+1);

    gcpnex.resize(imx+1, jmx+1);
    gcpnez.resize(imx+1, jmx+1);
    gcptex.resize(imx+1, jmx+1);
    gcptez.resize(imx+1, jmx+1);

    gnuobx.resize(imx+1, jmx+1);
    gnuoby.resize(imx+1, jmx+1);
    gupae0.resize(imx+1, jmx+1);

        // Calder Edit
        // real(8),dimension(:,:),allocatable :: phiavg
        // Calder Edit End

    //3D Arrays
    rk_hand.resize(imx+1, jmx+1, kmx+1);
    rho.resize(imx+1, jmx+1, kmx+1);
    phi.resize(imx+1, jmx+1, kmx+1);//!,den_pre(0:imx,0:jmx,0:kmx),dden(0:imx,0:jmx,0:kmx))
    //allocate(phiavg(0:imx,0:jmx))!,0:kmx)) !Calder Edit

    ex.resize(imx+1, jmx+1, kmx+1);
    ez.resize(imx+1, jmx+1, kmx+1); //note: +1 since arrays are 0-imx with imx included in ftn
    ezeta.resize(imx+1, jmx+1, kmx+1);

    delbx.resize(imx+1, jmx+1, kmx+1); 
    delbz.resize(imx+1, jmx+1, kmx+1); 
    delby.resize(imx+1, jmx+1, kmx+1);

    apar.resize(imx+1, jmx+1, kmx+1);
    dene.resize(imx+1, jmx+1, kmx+1);

    upar.resize(imx+1, jmx+1, kmx+1);
    jpar.resize(imx+1, jmx+1, kmx+1);
    upars.resize(imx+1, jmx+1, kmx+1);
    phis.resize(imx+1, jmx+1, kmx+1);
    denes.resize(imx+1, jmx+1, kmx+1);
    apars.resize(imx+1, jmx+1, kmx+1);

    dnedx.resize(imx+1, jmx+1, kmx+1);
    dnedy.resize(imx+1, jmx+1, kmx+1);
    dupadx.resize(imx+1, jmx+1, kmx+1);
    dupady.resize(imx+1, jmx+1, kmx+1);

    den.resize(2, imx+1, jmx+1, kmx+1); //1st index 1 based - careful
}

//cleans 1d arrays in COM
void cleanupCom(){
    delete[] mm;
    delete[] mims;
    //delete[] time;
    delete[] xg;
    delete[] zg;
    delete[] jac;
    delete[] fe;
    delete[] te;
    delete[] rmsphi;
    delete[] rmsez;
    delete[] rmsapa;
    delete[] avewi;
    delete[] vol;
    delete[] mu;
    delete[] x2;
    delete[] zeta2;
    delete[] z2;
    delete[] u2;
    delete[] x3;
    delete[] zeta3;
    delete[] z3;
    delete[] u3;
    delete[] w2;
    delete[] w3;
    delete[] gw;
}