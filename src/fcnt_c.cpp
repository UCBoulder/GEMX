#include "fcnt.hpp"

#include <chrono>
#include <iostream>

//Function used in gemx.f90
double revers_c(const int& num, const int& n){
   //const auto start = std::chrono::high_resolution_clock::now();
   double rev = 0.0;
   double power = 1.0;
   int inum = num;
   int iquot = 0;
   int irem = 0;

   while(inum > 0){
      iquot = int(inum/(n));
      irem = inum - n*iquot;
      power = power/(n);
      rev = rev + irem*power;
      inum = iquot;
   }
   return rev; 

   /* Timing
   const auto end = std::chrono::high_resolution_clock::now();
   const std::chrono::duration<double> diff = end - start;
   printf("%.10f\n",diff);
   */
}

