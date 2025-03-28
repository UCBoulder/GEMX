#include "gemx_com_externs.h"
#include "equil_externs.h"

//2D Pointers t0 use in c
Array2D<double> b0_c;
Array2D<double> t0i_c; 
Array2D<double> xn0i_c; 
Array2D<double> ileft_c;
Array2D<double> xbackw_c;
Array2D<double> zbackw_c;
Array2D<double> jleft_c;
Array2D<double> iright_c;
Array2D<double> xforw_c;
Array2D<double> jright_c;
Array2D<double> zforw_c;
Array2D<double> b0zeta_c;
Array2D<double> phiavg_c;
Array2D<double> t0s_c;
Array2D<double> xn0s_c;
Array2D<double> capts_c;
Array2D<double> capns_c;
Array2D<double> vpars_c;
Array2D<double> vparsp_c;
Array2D<double> psi_p_c;
Array2D<double> mask_c;
Array2D<double> mask2_c;
Array2D<double> mask3_c;
Array2D<double> mask4_c;
Array2D<double> c2_over_vA2_c;
Array2D<double> xn0e_c;
Array2D<double> t0e_c;
Array2D<double> den2d2_c;
Array2D<double> dden2d_c;
Array2D<double> den2d1_c;
Array2D<double> dbdx_c;
Array2D<double> dbdz_c;
Array2D<double> b0x_c;
Array2D<double> b0z_c;
Array2D<double> captix_c;
Array2D<double> captiz_c;
Array2D<double> capnix_c;
Array2D<double> capniz_c;

//3D pointer too be used in C
Array3D<double> phi_c;
Array3D<double> ex_c;
Array3D<double> ez_c;
Array3D<double> ezeta_c;
Array3D<double> phi_k_c;
Array3D<double> dphidr_c;
Array3D<double> dphi_kdr_c;
Array3D<double> dphidz_c;
Array3D<double> dphi_kdz_c;
Array3D<double> d2phidr2_c;
Array3D<double> d2phi_kdr2_c;
Array3D<double> d2phidz2_c;
Array3D<double> d2phi_kdz2_c;
Array3D<double> OPPphi_c;
Array3D<double> OPPphik_c;
Array3D<double> l_hand_c;
Array3D<double> r_hand_c;
Array3D<double> upar_c;
Array3D<double> apars_c;
Array3D<double> apar_c;
Array3D<double> jpar_c;
Array3D<double> rho_c;
Array3D<double> dene_c;
Array3D<double> curlb_c;
Array3D<double> delbx_c;
Array3D<double> delby_c;
Array3D<double> delbz_c;

//4D Arrays
Array4D<double> den_c;

void new_gemx_com_c_(){
    //call to functions that will calculate location of data in fortran arrays
    Allocate2dPointerArrays_gemx_com();
    Allocate3dPointerArrays_gemx_com();
    Allocate4dPointerArrays_gemx_com();
}

//Allocation of 2D and 3D arrays below
void Allocate2dPointerArrays_gemx_com(){
    b0_c.CreateArray2D(b0_ptr,nx+1,nz+1);
    t0i_c.CreateArray2D(t0i_ptr,nx+1,nz+1);
    xn0i_c.CreateArray2D(xn0i_ptr, nx+1, nz+1);
    ileft_c.CreateArray2D(ileft_ptr,imx+1,jmx+1);
    xbackw_c.CreateArray2D(xbackw_ptr,imx+1,jmx+1);
    zbackw_c.CreateArray2D(zbackw_ptr,imx+1,jmx+1);
    jleft_c.CreateArray2D(jleft_ptr,imx+1,jmx+1);
    iright_c.CreateArray2D(iright_ptr,imx+1,jmx+1);
    xforw_c.CreateArray2D(xforw_ptr,imx+1,jmx+1);
    jright_c.CreateArray2D(jright_ptr,imx+1,jmx+1);
    zforw_c.CreateArray2D(zforw_ptr,imx+1,jmx+1);
    b0zeta_c.CreateArray2D(b0zeta_ptr, nx+1, nz+1);
    t0s_c.CreateArray2D(t0s_ptr, 5, nr+1);
    xn0s_c.CreateArray2D(xn0s_ptr,5, nr+1);
    capts_c.CreateArray2D(capts_ptr,5,nr+1);
    capns_c.CreateArray2D(capns_ptr,5,nr+1);
    vpars_c.CreateArray2D(vpars_ptr,5,nr+1);
    vparsp_c.CreateArray2D(vparsp_ptr,5,nr+1);
    psi_p_c.CreateArray2D(psi_p_ptr, nx+1, nz+1);
    mask_c.CreateArray2D(mask_ptr, nx+1, nz+1);
    mask2_c.CreateArray2D(mask2_ptr, nx+1, nz+1);
    mask3_c.CreateArray2D(mask3_ptr, nx+1, nz+1);
    mask4_c.CreateArray2D(mask4_ptr, nx+1, nz+1);
    //phiavg_c.CreateArray2D(phiavg_ptr, nx, nz); //All instances currently being passed to function
    c2_over_vA2_c.CreateArray2D(c2_over_va2_ptr, nx+1, nz+1);
    xn0e_c.CreateArray2D(xn0e_ptr, nx+1, nz+1);
    t0e_c.CreateArray2D(t0e_ptr, nx+1, nz+1);
    den2d2_c.CreateArray2D(den2d2_ptr, imx+1, jmx+1);
    dden2d_c.CreateArray2D(dden2d_ptr, imx+1, jmx+1);
    den2d1_c.CreateArray2D(den2d1_ptr, imx+1, jmx+1);
    dbdx_c.CreateArray2D(dbdx_ptr, nx+1, nz+1);
    dbdz_c.CreateArray2D(dbdz_ptr, nx+1, nz+1);
    b0x_c.CreateArray2D(b0x_ptr, nx+1, nz+1);
    b0z_c.CreateArray2D(b0z_ptr, nx+1, nz+1);
    captix_c.CreateArray2D(captix_ptr, nx+1, nz+1);
    captiz_c.CreateArray2D(captiz_ptr, nx+1, nz+1);
    capnix_c.CreateArray2D(capnix_ptr, nx+1, nz+1);
    capniz_c.CreateArray2D(capniz_ptr, nx+1, nz+1);
}

void Allocate3dPointerArrays_gemx_com(){
    //phi_c.CreateArray3D(phi_ptr, imx, jmx, kmx); //All instances currently being passed to function
    ex_c.CreateArray3D(ex_ptr, imx+1, jmx+1, kmx+1);
    ez_c.CreateArray3D(ez_ptr, imx+1, jmx+1, kmx+1);
    ezeta_c.CreateArray3D(ezeta_ptr, imx+1, jmx+1, kmx+1);
    phi_k_c.CreateArray3D(phi_k_ptr, imx+1, jmx+1, kmx+1);
    dphidr_c.CreateArray3D(dphidr_ptr, imx+1, jmx+1, kmx+1);
    dphi_kdr_c.CreateArray3D(dphi_kdr_ptr, imx+1, jmx+1 ,kmx+1);
    dphidz_c.CreateArray3D(dphidz_ptr, imx+1, jmx+1, kmx+1);
    dphi_kdz_c.CreateArray3D(dphi_kdz_ptr, imx+1, jmx+1, kmx+1);
    d2phidr2_c.CreateArray3D(d2phidr2_ptr, imx+1, jmx+1, kmx+1);
    d2phi_kdr2_c.CreateArray3D(d2phi_kdr2_ptr, imx+1, jmx+1, kmx+1);
    d2phidz2_c.CreateArray3D(d2phidz2_ptr, imx+1, jmx+1, kmx+1);
    d2phi_kdz2_c.CreateArray3D(d2phi_kdz2_ptr, imx+1, jmx+1, kmx+1);
    OPPphi_c.CreateArray3D(oppphi_ptr, imx+1, jmx+1, kmx+1);
    OPPphik_c.CreateArray3D(oppphik_ptr, imx+1, jmx+1, kmx+1);
    l_hand_c.CreateArray3D(l_hand_ptr, imx+1, jmx+1, kmx+1);
    r_hand_c.CreateArray3D(r_hand_ptr, imx+1, jmx+1, kmx+1);
    upar_c.CreateArray3D(upar_ptr, imx+1, jmx+1, kmx+1);
    apars_c.CreateArray3D(apars_ptr, imx+1, jmx+1, kmx+1);
    apar_c.CreateArray3D(apar_ptr, imx+1, jmx+1, kmx+1);
    jpar_c.CreateArray3D(jpar_ptr, imx+1, jmx+1, kmx+1);
    rho_c.CreateArray3D(rho_ptr, imx+1, jmx+1, kmx+1);
    dene_c.CreateArray3D(dene_ptr, imx+1, jmx+1, jmx+1);
    curlb_c.CreateArray3D(curlb_ptr, nx+1, nz+1, 3); 
    delbx_c.CreateArray3D(delbx_ptr, imx+1, jmx+1, kmx+1);
    delby_c.CreateArray3D(delby_ptr, imx+1, jmx+1, kmx+1);
    delbz_c.CreateArray3D(delbz_ptr, imx+1, jmx+1, kmx+1);
}

void Allocate4dPointerArrays_gemx_com(){
    den_c.CreateArray4D(den_ptr, 2, imx+1, jmx+1, kmx+1);
}