#include "fcnt.hpp"
#include <iostream>

double revers_c_(const int& num, const int& n){
   double rev = 0.0;
   double power = 1.0;
   int inum = num;
   int iquot = 0;
   int irem = 0;

   while(inum > 0){
      iquot = int(inum/n);
      irem = inum - n*iquot;
      power = power/n;
      rev = rev + irem*power;
      inum = iquot;
   }
   return rev; 
}