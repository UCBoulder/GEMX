#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "gemx_com_c.hpp"
#include "equil_c.hpp"
#pragma once

void outd_c_(const int &n);
// void phixy(Array3D<double> grd, string fl, int unt, int n);
// void phixz(Array3D<double> grd, string fl, int unt, int n); //TOTO - weird error here - compiler doesn;t like int. Remove if magic happens