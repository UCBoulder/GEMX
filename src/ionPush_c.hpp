#include "gemx_com_c.hpp"
#include "equil_c.hpp"
#include <cmath>
#include <iostream>
#include <chrono>

#include <mpi.h>
#include <nvToolsExt.h>
#pragma once

void ppush_c_(const int &n);
void cpush_c_(const int &n);
inline void updateDeviceData();
inline void updateHostData();
