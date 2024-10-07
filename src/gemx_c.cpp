#include "gemx.hpp"
#include "fcnt.hpp"
#include <iostream>
#include <cmath>
#include <chrono>
using namespace std;

void parperp_c_(double* vpar,double* vperp2,int& m,int& cnt,int& MyId){
    double r1 = 0.0;
    double r2 = 0.0;
    double t = 0.0;

    const double c0 = 2.515517;
    const double c1 = 0.802853;
    const double c2 = 0.010328;

    const double d1 = 1.432788;
    const double d2 = 0.189269;
    const double d3 = 0.001308;
    double temp = 0.0;
    int iflag = 1;

    //const auto start = std::chrono::high_resolution_clock::now(); //timing stuff
    r1 = revers_c(m+MyId*cnt, 7); //sets values for r1 and r2 (random numbers)
    r2 = revers_c(m+MyId*cnt, 11); 

        //.....quiet start---see denavit pf '71(?) & abramowitz hand book
        //.....fibonacci start---see denavit comm. pla. phy. & con. fus. '81
        // warning: we have g1=1 in the x-direction. This surpresses all odd
        //          modes in the x-direction!!!

    
    
    if (r1 <= 0.5) goto jump;

    r1 = 1.0 - r1;
    iflag = -1;
    
    jump:
    if(r1 >= 1.0e-6){
        t = sqrt(log(1.0/(r1*r1)));
    }
    else{
        t = 5.0;
        cout << "parperp2 warning m= " << m << endl;
    }

    *vpar = t-(c0+c1*t+c2*(t*t))/(1.+d1*t+d2*(t*t)+d3*(t*t*t));
    *vpar = temp*iflag;

    *vperp2 = -2.0*log(r2); 
    
    //timing stuff
    // const auto end = std::chrono::high_resolution_clock::now();
    // const std::chrono::duration<double> diff = end - start;
    // printf("%.10f\n",diff);

    return;
}