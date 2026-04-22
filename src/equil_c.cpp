#include "equil_c.hpp"
//As a note to everyone, all arrays are made to be inclusive of final index like fortran. This means a 2D array of size nx+1, nz+1 has a final index of nx,nz. This is how fortran does it

double mimp = 2, chgi = 1;
double betaVal,rmaj0,a,q0,r0,q0p,q0abs,shat0;
double phi_diag,phi_diag_freq,weight_diag;
double dR,dth,mu0,e,proton;
int nr=200,nr2=100,ntheta=200,isgnf=1,isgnq=-1,isupae0=0,tor_n;
double psi_max=0.31, psi_min=-1.0 ,R_min=1.0, Z_min=-1.5, Z_internal=-1.2, psi_div=0.305,psi_a=0.311647;

// GEM-X
//integer :: cont=259002 !Calder Edit

int nzeta = 64; 
int nx=299, nz=299; //indexing problems
// int nx=999, nz=999;
double zctr;
double  dxeq, xdim, xctr, zdim, dzeq;
double pi,pi2;

CArray2D<double> b0, b0x, b0z, b0zeta,dbdx,dbdz, c2_over_vA2, q_grid;
CArray2D<double> t0i,t0e,xn0i,xn0e,captix,captex,capnix,capnex,captiz,captez,capniz,capnez;
CArray2D<double> t0D,xn0D;

double *psi = nullptr;
double *psip = nullptr;
double *sf = nullptr;
double *vpari = nullptr;
double *vparip = nullptr;
double *zeff = nullptr;
double *nue0 = nullptr;
double *phinc = nullptr;
double *phincp = nullptr;
double *er = nullptr;
double *upari = nullptr;
double *Rgrid = nullptr;
double *Zgrid = nullptr;

//careful, next declaration first index 1 based in ftn
CArray2D<double> t0s,xn0s,capts,capns,vpars,vparsp,psi_p,mask,mask2,mask3,mask4;
double bu,tu,nu,xu,frequ,vu,eru;

//     for including bstar effects (this was the only one allocated)
CArray2D<double> bdcrvb;
CArray2D<double> rho_i;
CArray2D<double> dpsi_dr;
CArray2D<double> dpsi_dz;
/*double *psip2 = new double[?] //this array isn't allocated, keeping just in case
also these arrays are not allocated (2D):  curvbz,srbr,srbz,thbr,thbz,prsrbr,prsrbz,pthsrbr,pthsrbz*/

// for equilibrium current 
CArray2D<double> upae0,nuob,dnuobdr,dnuobdt;
CArray3D<double> curlb; //last index 1 based in ftn - cPP:(0,1,2) vs ftn:(1,2,3)

// for phi average Calder Edit
// int num_lines = 80817; //NON-CBC VERSION
int num_lines = 44816;    //CBC VERSION
// int num_lines = 78427;
int line;
CArray2D<double> phiavg;
double *psitab;
double *weight00;
double *weight01;
double *weight10;
double *weight11;
double *jacobian;
double *deno;
int *gindex;
int *iarray;
int *jarray;
int *priv;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void new_equil_c(){
    using namespace std;
    int i;
    double omegau;
    
    //allocating arrays
    psitab = new double[num_lines+1]; std::fill(psitab, psitab+num_lines+1, 0.0);
    weight00 = new double[num_lines+1]; std::fill(weight00, weight00+num_lines+1, 0.0);
    weight01 = new double[num_lines+1]; std::fill(weight01, weight01+num_lines+1, 0.0);
    weight10 = new double[num_lines+1]; std::fill(weight10, weight10+num_lines+1, 0.0);
    weight11 = new double[num_lines+1]; std::fill(weight11, weight11+num_lines+1, 0.0);
    jacobian = new double[num_lines+1]; std::fill(jacobian, jacobian+num_lines+1, 0.0);
    deno = new double[num_lines+1];std::fill(deno, deno+num_lines+1, 0.0);
    gindex = new int[num_lines+1]; std::fill(gindex, gindex+num_lines+1, 0);
    iarray = new int[num_lines+1]; std::fill(iarray, iarray+num_lines+1, 0);
    jarray = new int[num_lines+1]; std::fill(jarray, jarray+num_lines+1, 0);
    priv = new int[num_lines+1]; std::fill(priv, priv+num_lines+1, 0);

    psi = new double[nr+1]; std::fill(psi, psi+nr+1, 0.0);
    psip = new double[nr+1]; std::fill(psip, psip+nr+1, 0.0);
    sf = new double[nr+1]; std::fill(sf, sf+nr+1, 0.0);
    vpari = new double[nr+1]; std::fill(vpari, vpari+nr+1, 0.0);
    vparip = new double[nr+1]; std::fill(vparip, vparip+nr+1, 0.0);
    zeff = new double[nr+1]; std::fill(zeff, zeff+nr+1, 0.0);
    nue0 = new double[nr+1]; std::fill(nue0, nue0+nr+1, 0.0);
    phinc = new double[nr+1]; std::fill(phinc, phinc+nr+1, 0.0);
    phincp = new double[nr+1]; std::fill(phincp, phincp+nr+1, 0.0);
    er = new double[nr+1]; std::fill(er, er+nr+1, 0.0);
    upari = new double[nr+1]; std::fill(upari, upari+nr+1, 0.0);
    Rgrid = new double[nx+1]; std::fill(Rgrid, Rgrid+nx+1, 0.0);
    Zgrid = new double[nz+1]; std::fill(Zgrid, Zgrid+nz+1, 0.0);

    b0.resize(nx+1,nz+1), b0x.resize(nx+1,nz+1), b0z.resize(nx+1,nz+1),
                 b0zeta.resize(nx+1,nz+1),dbdx.resize(nx+1,nz+1),dbdz.resize(nx+1,nz+1),
                 c2_over_vA2.resize(nx+1,nz+1),q_grid.resize(nx+1,nz+1);
    t0i.resize(nx+1,nz+1),t0e.resize(nx+1,nz+1),xn0i.resize(nx+1,nz+1),xn0e.resize(nx+1,nz+1),
                 captix.resize(nx+1,nz+1),captex.resize(nx+1,nz+1),capnix.resize(nx+1,nz+1),capnex.resize(nx+1,nz+1),
                 captiz.resize(nx+1,nz+1),captez.resize(nx+1,nz+1),capniz.resize(nx+1,nz+1),capnez.resize(nx+1,nz+1);

    t0D.resize(nx+1,nz+1),xn0D.resize(nx+1,nz+1);

    t0s.resize(5,nr+1),xn0s.resize(5,nr+1),capts.resize(5,nr+1),capns.resize(5,nr+1),
                 vpars.resize(5,nr+1),vparsp.resize(5,nr+1),psi_p.resize(nx+1,nz+1),mask.resize(nx+1,nz+1),
                 mask2.resize(nx+1,nz+1),mask3.resize(nx+1,nz+1),mask4.resize(nx+1,nz+1);

    bdcrvb.resize(nx+1,nz+1);
    
    // for equilibrium current 
    upae0.resize(nr+1,ntheta+1),nuob.resize(nr+1,ntheta+1),dnuobdr.resize(nr+1,ntheta+1),
                    dnuobdt.resize(nr+1,ntheta+1);
    curlb.resize(nx+1, nz+1, 3); //last index 1 based in ftn - cPP:(0,1,2) vs ftn:(1,2,3)

    phiavg.resize(nx+1, nz+1);
    rho_i.resize(nx+1, nz+1);
    dpsi_dr.resize(nx+1,nz+1);
    dpsi_dz.resize(nx+1,nz+1);

    //global equilibrium data 
    //open R.dat and input into Rgrid
    read1D("R.dat", Rgrid, 0);
    //open Z.dat and input into Zgrid
    read1D("Z.dat", Zgrid, 0);
    //open psi_p.dat and input into psi_p
    read2D("psi_p.dat", psi_p, nx, nz);

    // !Calder Edits Start
    // !character(len=100) :: line_buffer
    // ! Open the dataset file
    //read jacodata.dat into gindex, psitab, iarray, jarray, weight00, wheight01, wheight10, wehight11, jacobian, deno, priv
    std::string delimeter = "	";
    std::string line;
    std::string num;
    std::ifstream file;
    file.open("jacodata_cbc.dat");
    int index = 0;
    int currI = 0;

    while(std::getline(file, line)) {
        std::stringstream str(line);
        while(getline(str, num, delimeter[0])) {
            int res = index % 11;
            switch (res)
            {
            case 0: 
            gindex[currI] = stod(num); 
            break;
            case 1: 
            psitab[currI] = stod(num); 
            break;
            case 2: 
            iarray[currI] = stod(num); 
            break;
            case 3: 
            jarray[currI] = stod(num); 
            break;
            case 4: 
            weight00[currI] = stod(num); 
            break;
            case 5: 
            weight10[currI] = stod(num); 
            break;
            case 6: 
            weight01[currI] = stod(num); 
            break;
            case 7: 
            weight11[currI] = stod(num);
             break;
            case 8: 
            jacobian[currI] = stod(num); 
            break;
            case 9: 
            deno[currI] = stod(num); 
            break;
            case 10: 
            priv[currI] = stod(num); 
            currI+=1; 
            break;
            default: break;
            }
            index+=1;
        }
    }
    //Calder Edits Finish 



//         open(unit=11, file = 'psi_test',status='unknown',action='write')
//      do i=0,nx
//          write(11,*) psi_p
//       end do
//       close(11) 



    xdim = Rgrid[nx] - Rgrid[0];
    zdim = Zgrid[nz] - Zgrid[0];
    dR = (Rgrid[nx]-Rgrid[0])/nx;
//    dx=(Rgrid(nx)-Rgrid(0))/nx
//    dZ=(Zgrid(nz)-Zgrid(0))/nz

//  open BR.dat and iput into b0x
    read2D("BR.dat", b0x, nx, nz);

//    open(unit=11, file = 'debug.dat',status='unknown',action='write')
//       write(11,*) b0x(5,0), 'dR='dR
//       close(11)

//  open Bz.dat and input into b0z
    read2D("Bz.dat", b0z, nx, nz);
//  open Bt.dat and input into b0zeta
    read2D("Bt.dat", b0zeta, nx, nz);
//       B0x=0
//       B0z=0  

    for(int i = 0; i <= nx; ++i){
        for(int j = 0; j <= nz; ++j){
            b0(i,j) = sqrt((b0x(i,j)*b0x(i,j))+(b0z(i,j)*b0z(i,j))+(b0zeta(i,j)*b0zeta(i,j)));
        }
    }
    
    e = 1.6e-19;
    mu0 = 1.25663706212e-6;
    proton = 1.67e-27;
    // bu = 1.98562;

    bu = 1; //?

    tu = 1000*e;
    omegau = e*bu/proton;
    frequ = omegau;

    // vu = sqrt(tu/proton);
    vu = 1; //?
    // xu = proton*vu/(e*bu);
    xu = 1; //?
    nu = 2.5e19;
    betaVal = 4*3.14159*1e-7*nu*tu/(bu*bu); //is this the whole thing squared or bu squared? Double check if ever used
    
//     assign T, n profiles 
    // for(int i = 0; i <= nx; ++i){
    //     for(int j = 0; j <= nz; ++j){
    //         t0i(i,j) = 1.*tu;
    //         t0e(i,j) = 1.*tu;
    //         xn0i(i,j) = 1.*nu;
    //         xn0e(i,j) = 1.*nu;
    //     }
    // }
    //   Calder Edit: realistic profilies for ITG runs
    //   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //    open(unit=10, file = 'Profiles/ti0_profile.dat',status='old',action='read')
    //    read(10,*) t0i
    //    close(10)

    //    open(unit=10, file = 'Profiles/te0_profile.dat',status='old',action='read')
    //    read(10,*) t0e
    //    close(10)

    //    open(unit=10, file = 'Profiles/ni0_profile.dat',status='old',action='read')
    //    read(10,*) xn0i
    //    close(10)

    //    open(unit=10, file = 'Profiles/ne0_profile.dat',status='old',action='read')
    //    read(10,*) xn0e
    //    close(10) 
    
//   Manufactured Temperature Profiles: 0.5*(-np.tanh(5*np.array(profiles['psinorm'])-2.5)+1)*np.max(np.array(profiles['ti']))
    // GORLER PROFILES
    // read2D("ni0_gorler.dat", xn0i, nx, nz);
    // read2D("ne0_gorler.dat", xn0e, nx, nz);
    // read2D("ti0_gorler.dat", t0i, nx, nz); //problem likley here, so likely also with others
    // read2D("te0_gorler.dat", t0e, nx, nz);
    read2D("ni0_wpqh.dat", xn0i, nx, nz);
    read2D("ne0_wpqh.dat", xn0e, nx, nz);
    read2D("ti0_wpqh.dat", t0i, nx, nz);
    read2D("te0_wpqh.dat", t0e, nx, nz);
    // read2D("ni0_wpqh_ideal.dat", xn0i, nx, nz);
    // read2D("ne0_wpqh_ideal.dat", xn0e, nx, nz);
    // read2D("ti0_wpqh_ideal.dat", t0i, nx, nz);
    // read2D("te0_wpqh_ideal.dat", t0e, nx, nz);
    // printf("Testing");
    //PUT INTO SI UNITS
    for(int i = 0; i <= nx; ++i) {
        for(int j = 0; j <=nz; ++j) {
            t0i(i,j) = t0i(i,j)*tu;
            t0e(i,j) = t0e(i,j)*tu;
            xn0i(i,j) = xn0i(i,j)*1e19;
            xn0e(i,j) = xn0e(i,j)*1e19;

            t0D(i,j) = 0.01*t0i(i,j);
            xn0D(i,j) = xn0i(i,j)*0.01;
        }
    }

    //Put into SI units
    // for(int i = 0; i <= nx; ++i){
    //     for(int j = 0; j <= nz; ++j){
    //         t0i(i,j) = t0i(i,j)*tu;
    //         t0e(i,j) = t0e(i,j)*tu;
    //         // xn0i = xn0i*1.0e21
    //         // xn0e = xn0e*1.0e21
    //         // xn0i = xn0e 
    //     }
    // }   
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    // //open ne0.dat input into xn0e
    // read2D("ne0.dat", xn0e, nx, nz);
    // //open ti0.dat input into t0i
    // read2D("ti0.dat", t0i, nx, nz);

    // UNIFORM PROFILES
    // for(int i = 0; i <= nx; ++i){
    //     for(int j = 0; j <= nz; ++j){
    //         xn0e(i,j) = xn0e(i,j)*nu;
    //         t0i(i,j) = t0i(i,j)*tu;
    //     }
    // }
    // xn0i = xn0e; 

    for(int i = 0; i <= nx; ++i){
        for(int j = 0; j <= nz; ++j){
            q_grid(i,j) = 2.52 * pow(sqrt(pow(Rgrid[i] - Rgrid[imx/2],2) + pow(Zgrid[j]-Zgrid[jmx/2],2))/0.6012,2) -0.16*(sqrt(pow(Rgrid[i] - Rgrid[imx/2],2) + pow(Zgrid[j]-Zgrid[jmx/2],2)))/0.6012 + 0.86;
            if(xn0e(i,j)<1e-9){
                c2_over_vA2(i,j)=mu0*2*proton*0.01*xn0e(nx/2,nz/2)/(b0(i,j)*b0(i,j))*(vu*vu);
            } else {
                c2_over_vA2(i,j)=mu0*2*proton*xn0e(i,j)/(b0(i,j)*b0(i,j))*(vu*vu);//2*Rgrid(i)**2/(Rgrid(0)+Rgrid(nx))**2
            }
            rho_i(i,j) = sqrt(t0i(i,j)/(2*proton))*(2*proton)/(e*b0(i,j));
            // if (c2_over_vA2(i,j)<1.0e-19) c2_over_vA2(i,j)=0.01*mu0*2*proton*xn0e(i/2,j/2)/(b0(i/2,j/2)**2)*vu**2
        }
    }

    mask.Clear();
    //c2_over_vA2=1
    cout << psi_min << endl;
    for(int i = 2; i <= nx-2; ++i){ 
        for(int j = 2; j <= nz-2; ++j){ 
//          c2_over_vA2(i,j)=mu0*2*proton*xn0e(i,j)/(b0(i,j)**2)*vu**2!2*Rgrid(i)**2/(Rgrid(0)+Rgrid(nx))**2
//          write(*,*) c2_over_vA2(i,j)
            if(psi_p(i,j)<psi_max && psi_p(i,j)>psi_min && (Zgrid[j]>Z_internal || (Zgrid[j] && psi_p(i,j)>psi_div)) && Rgrid[i]>R_min && Zgrid[j]<1.2){
               // if(psi_p(i,j)<0.3 && Zgrid(j)>-1.08){}
            //    if(psi_p(i,j)>0.21 && Zgrid[j]<-1.2)
                mask(i,j) = 1;
                if (psi_p(i,j)<0.21 && Zgrid[j]<-1.2){
                    mask(i,j) = 0;
                }
            } else {
                mask(i,j) = 0;
            }
        }
    }

    for(int i = 2; i <= nx-2; ++i){
        for(int j = 2; j <= nz-2; ++j){
            if((mask(i+1,j) * mask(i-1,j) * mask(i,j+1) * mask(i,j-1)) == 0){ //idea: bitwise & can make this faster technically
                mask2(i,j) = 0;
            } else {
                mask2(i,j) = 2;
            }   
        }
    }

    for(int i = 2; i <= nx-2; ++i){
        for(int j = 2; j <= nz-2; ++j){
            if((mask2(i+1,j) * mask2(i-1,j) * mask2(i,j+1) * mask2(i,j-1))==0){ //idea: bitwise & can make this faster technically
                mask3(i,j) = 0;
            } else {
                mask3(i,j) = 3;
            }   
        }
    }

    for(int i = 2; i <= nx-2; ++i){
        for(int j = 2; j <= nz-2; ++j){
            if((mask3(i+1,j) * mask3(i-1,j) * mask3(i,j+1) * mask3(i,j-1))==0){ //idea: bitwise & can make this faster technically
                mask4(i,j) = 0;
            } else {
                mask4(i,j) = 4;
            }   
        }
    }

    if(myid==0){ 
        ofstream file;
        file.open("test_1_over_vA2");
         for(int j = 0; j <= jmx; ++j)  {
            for(int i = 0; i <= imx; ++i) {
               file << c2_over_vA2(i,j) << "    ";
            }
            file << "\n";
         }
         file.close();
        
        file.open("mask");
         for(int j = 0; j <= jmx; ++j)  {
            for(int i = 0; i <= imx; ++i) {
               file << mask(i,j) << "    ";
            }
            file << "\n";
         }
         file.close();

        file.open("mask2");
         for(int j = 0; j <= jmx; ++j)  {
            for(int i = 0; i <= imx; ++i) {
               file << mask2(i,j) << "    ";
            }
            file << "\n";
         }
         file.close();

         file.open("mask3");
         for(int j = 0; j <= jmx; ++j)  {
            for(int i = 0; i <= imx; ++i) {
               file << mask3(i,j) << "    ";
            }
            file << "\n";
         }
         file.close();

         file.open("mask4");
         for(int j = 0; j <= jmx; ++j)  {
            for(int i = 0; i <= imx; ++i) {
               file << mask4(i,j) << "    ";
            }
            file << "\n";
         }
         file.close();
    }

    //Normalization
    for(int i = 0; i <= nx; ++i){
        for(int j = 0; j <= nz; ++j){
            b0(i,j) = b0(i,j)/bu;
            b0x(i,j) = b0x(i,j)/bu;
            b0z(i,j) = b0z(i,j)/bu;
            b0zeta(i,j) = b0zeta(i,j)/bu;
        }
    }

//  cout << "betaU=" << beta << "\n";

    pi = atan(1.0)*4.0;
    pi2 = 2*pi;
    rmaj0 = 1000.;
    a = 360.;

//  xctr = a*1.5
//  zctr = 0.

    // xctr = 0.5*(Rgrid[nx]+Rgrid[0])/xu;
    // zctr = 0.5*(Zgrid[nz]+Zgrid[0])/xu;
    xctr = Rgrid[158]; // WPQH mode specific
    zctr = Zgrid[148]; // WPQH mode specific
    xdim = xdim/xu;
    zdim = zdim/xu;
//  xdim = a*2;
//  zdim = a*3
    dxeq = xdim/nx;
    dzeq = zdim/nz;

    dxeq = abs(Rgrid[1]-Rgrid[0]);
    dzeq = abs(Zgrid[1]-Zgrid[0]);
    
    for(int i = 1; i <= nx-1; ++i){
        for(int j = 1; j <= nz-1; ++j){
            curlb(i,j,0)=-((-b0zeta(i,j+1)/b0(i,j+1)+b0zeta(i,j-1)/b0(i,j-1))/(2*dzeq));
            curlb(i,j,1)=-(1/Rgrid[i]*(Rgrid[i+1]*b0zeta(i+1,j)/b0(i+1,j)-Rgrid[i-1]*b0zeta(i-1,j)/b0(i-1,j))/(2*dxeq));
            curlb(i,j,2)=-((b0x(i,j+1)/b0(i,j+1)-b0x(i,j-1)/b0(i,j-1))/(2*dzeq)-(b0z(i+1,j)/b0(i+1,j)-b0z(i-1,j)/b0(i-1,j))/(2*dxeq));
//            write(*,*) curlb(i,j,2)+(b0zeta(i,j)/b0(i,j)/Rgrid(i)*xu+(b0zeta(i+1,j)/b0(i+1,j)-b0zeta(i-1,j)/b0(i-1,j))/(2*dxeq))
        }
    }

    for(int i = 1; i <= nx-1; ++i){
        for(int j = 1; j <= nz-1; ++j){
            bdcrvb(i,j)=-1/b0(i,j)*(b0x(i,j)*(-b0zeta(i,j+1)/b0(i,j+1)+b0zeta(i,j-1)/b0(i,j-1))/(2*dzeq)
                                              +b0z(i,j)/Rgrid[i]*(Rgrid[i+1]*b0zeta(i+1,j)/b0(i+1,j)-Rgrid[i-1]*b0zeta(i-1,j)/b0(i-1,j))/(2*dxeq)
                                              +b0zeta(i,j)*((b0x(i,j+1)/b0(i,j+1)-b0x(i,j-1)/b0(i,j-1))/(2*dzeq)-(b0z(i+1,j)/b0(i+1,j)-b0z(i-1,j)/b0(i-1,j))/(2*dxeq)));
        }
    }

//   Gradients of T and n profiles for global simulation. Flux-tube parameters done at the end
    // captix.Clear();
    // captex.Clear();
    // capnix.Clear();
    // capnex.Clear();
    // captiz.Clear(); //shouldn't need this since array filled with zeros upon initialization
    // captez.Clear();
    // capniz.Clear();
    // capnez.Clear();

    // dspi_dr.Clear();
    // dpsi_dz.Clear();
    
    //-CALDER EDIT 05/13/2025 ("<" specifically used - not a mistake)----
    for(int i = 1; i < nx; ++i) {
        for(int j = 1; j < nz; ++j) {
            captix(i,j) = (-1/t0i(i,j))*(t0i(i+1,j)-t0i(i-1,j))/(2*dxeq);
            captex(i,j) = (-1/t0e(i,j))*(t0e(i+1,j)-t0e(i-1,j))/(2*dxeq);
            captiz(i,j) = (-1/t0i(i,j))*(t0i(i,j+1)-t0i(i,j-1))/(2*dzeq);
            captez(i,j) = (-1/t0e(i,j))*(t0e(i,j+1)-t0e(i,j-1))/(2*dzeq);
            capnix(i,j) = (-1/xn0i(i,j))*(xn0i(i+1,j)-xn0i(i-1,j))/(2*dxeq);
            capnex(i,j) = (-1/xn0e(i,j))*(xn0e(i+1,j)-xn0e(i-1,j))/(2*dxeq);
            capniz(i,j) = (-1/xn0i(i,j))*(xn0i(i,j+1)-xn0i(i,j-1))/(2*dzeq);
            capnez(i,j) = (-1/xn0e(i,j))*(xn0e(i,j+1)-xn0e(i,j-1))/(2*dzeq);

            dpsi_dr(i,j) = (psi_p(i+1,j)-psi_p(i-1,j))/(2*dxeq);
            dpsi_dz(i,j) = (psi_p(i,j+1)-psi_p(i,j-1))/(2*dzeq);
        }
    }
    //-------------------------------------------------------------------

    //dbdx.Clear();
    //dbdz.Clear();
    for(int i = 1; i < nx; ++i){
        for(int j = 1; j < nz; ++j){
            dbdx(i,j) = (b0(i+1,j)-b0(i-1,j))/(2*dxeq);
            dbdz(i,j) = (b0(i,j+1)-b0(i,j-1))/(2*dzeq); 
        }
    }

    if(myid == 0){
        //open file to write xu, omegau, vu
        //do later, really not as important
        //format();
    }
}

//delete's 1D array pointers in EQUIL
void cleanUpEquil(){
    delete[] psi;
    delete[] psip;
    delete[] sf;
    delete[] vpari;
    delete[] vparip;
    delete[] zeff;
    delete[] nue0;
    delete[] phinc;
    delete[] phincp;
    delete[] er;
    delete[] upari;
    delete[] Rgrid;
    delete[] Zgrid;
    delete[] psitab;
    delete[] weight00;
    delete[] weight01;
    delete[] weight10;
    delete[] weight11;
    delete[] jacobian;
    delete[] deno;
    delete[] gindex;
    delete[] iarray;
    delete[] jarray;
    delete[] priv;
}