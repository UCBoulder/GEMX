#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "MultiArraysC.hpp"

#pragma once

void read1D(std::string fname, double arr[], int dflag);
void read2D(std::string fname, CArray2D<double> &arr, int x, int y);