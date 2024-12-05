#pragma once

// void revers_c_(const int& num, const int& n, double* const r); //double* const r (use as argument if calling from fortran)

extern "C"
{
  double revers_c_(const int& num, const int& n); // double* revers (use as argument if calling from fortran)
}
