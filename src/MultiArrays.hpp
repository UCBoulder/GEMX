#pragma once

#include <cstddef>

// 2D array class with column-major order
template<typename T>
class Array2D {
public:
   Array2D(){}

   void CreateArray2D(const T* const data, const std::size_t x, const std::size_t y){
      data_ = const_cast<T*>(data);
      x_ = x+1; //for fortran arrays starting at 0, edge cases need to be included for max value. if array starts at 1,1 need to subtract one from index
      y_ = y+1;
   }

   // Column-major access
   T& operator()(const std::size_t i, const std::size_t j) {
      return data_[j * x_ + i];  // Column-major order
   }

private:
   T* data_ = nullptr;
   std::size_t x_ = 0;
   std::size_t y_ = 0;
};

// 3D array class with column-major order
template<typename T>
class Array3D {
public:
   Array3D(){}

   void CreateArray3D(const T* const data, const std::size_t x, const std::size_t y, const std::size_t z){
      data_ = const_cast<T*>(data);
      x_ = x+1;
      y_ = y+1;
      z_ = z+1;
   }

   // Column-major access
   inline T& operator()(const std::size_t i, const std::size_t j, const std::size_t k) {
      return data_[(k * y_ + j) * x_ + i];  // Column-major order
   }

private:
   T* data_ = nullptr;
   std::size_t x_ = 0;
   std::size_t y_ = 0;
   std::size_t z_ = 0;
};