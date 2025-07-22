#include "outd_c.hpp"
using namespace std;

void outd_c_(const int &n) { 
    ofstream outTracer;
    ofstream testNe;
    ofstream testPhi;
    if(myid == 0) {
        printf("timestep = %d\n", timestep);
    }

    outTracer.open("./out/tracer.out", ios::app); 
    for(int m = 0; m < ntracer; ++m) {
        outTracer << "            " << timestep << "            " << m << "    " << setprecision(16) << (x3[m])*xu+Rgrid[0] << "        " << (z3[m])*xu+Zgrid[0] << "\n";
    }
    outTracer.close();

    testNe.open("./out/testne");
    if(testNe.is_open()) {
        for(int i = 0; i <= imx; ++i) {
            for(int j = 0; j <= jmx; ++j) {  
                testNe << fixed << setprecision(16) << dene(i,j,0) << "    ";
            }
            testNe << "\n";
        }
        testNe.close();
    } else {
        cout << "Warning testne failed to open/ wasn't created" << "\n";
    }

    testPhi.open("./out/testphi");
    if(testPhi) {
        for(int i = 0; i <= imx; ++i) {
            for(int j = 0; j <= jmx; ++j) {
                testPhi << fixed << setprecision(16) << phi(i,j,0) <<  "    "; 
            }
            testPhi << "\n";
        }
        testPhi.close();
    } else {
        cout << "Warning testphi failed to open/ wasn't created" << "\n";
    }


//  if((n % nplot) == 0){                       //no point in check since everything commented at this time
// !SEP	   call phixy(phi(:,:,:),'phixy',31,n)  //if we use, we would update the arrays in C to take Array3D objects
// !SEP	   call phixz(phi(:,:,:),'phixz',32,n)

// !SEP	   call phixy(apar(:,:,:),'apaxy',36,n)
// !SEP       call phixz(apar(:,:,:),'apaxz',37,n)
//    }

}

// void phixy(Array3D<double> grd, string fl, int unt, int n) { //TODO - later

// }

// void phixz(Array3D<double> grd, string fl, int unt, int n) { 

// }