#pragma once

#include <cstddef>
#include <stdio.h>
#include <algorithm>

// 2D array class with column-major order
template<typename T>
class Array2D {
public:
   Array2D(){}

   void CreateArray2D(const T* const data, const std::size_t x, const std::size_t y){
      data_ = const_cast<T*>(data);
      x_ = x; //for fortran arrays starting at 0, edge cases need to be included for max value. if array starts at 1,1 need to subtract one from index
      y_ = y;
   }

   void Clear(){
      if (data_ != nullptr){std::fill(data_,data_+x_*y_,0);}
   }

   // Column-major access
   T& operator()(const std::size_t i, const std::size_t j) {
      return data_[j * x_ + i];  // Column-major order
   }

   inline Array2D& operator=(const Array2D &arr) {
      for(auto i = 0; i < size_; ++i){
         data_[i] = arr.data_[i];
      }
      return *this;
   }

private:
   T* data_ = nullptr;
   std::size_t x_ = 0;
   std::size_t y_ = 0;
   std::size_t size_ = 0;
};

// 3D array class with column-major order
template<typename T>
class Array3D {
public:
   Array3D(){}

   void CreateArray3D(const T* const data, const std::size_t x, const std::size_t y, const std::size_t z){
      data_ = const_cast<T*>(data);
      x_ = x; 
      y_ = y;
      z_ = z;
      size_ = x * y * z;
   }

   void Clear(){
      if (data_ != nullptr){std::fill(data_,data_+x_*y_*z_,0);}
   }

   inline T& operator()(const std::size_t i, const std::size_t j, const std::size_t k) {
      return data_[(k * y_ + j) * x_ + i];  // Column-major order
   }

   inline Array3D& operator=(const Array3D &arr) {
      for(auto i = 0; i < size_; ++i){
         data_[i] = arr.data_[i];
      }
      return *this;
   }

private:
   T* data_ = nullptr;
   std::size_t x_ = 0;
   std::size_t y_ = 0;
   std::size_t z_ = 0;
   std::size_t size_ = 0;
};

// 4D array class with column-major order
template<typename T>
class Array4D {
   public:
   Array4D(){}

   void CreateArray4D(const T* const data, const std::size_t x, const std::size_t y, const std::size_t z, const std::size_t q){
      data_ = data_ = const_cast<T*>(data);
      x_ = x;
      y_ = y;
      z_ = z;
      q_ = q;
      size_ = x * y * z * q;
   }

   void Clear(){
      if (data_ != nullptr){std::fill(data_,data_+x_*y_*z_*q_,0);}
   }

   // Column-major access
   inline T& operator()(const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t l) {
      //return data_[i + x_ * (j + y_ * (k + z_ * l))];  // Column-major order
      return data_[((l * z_ + k) * y_ + j) * x_ + i];
   }

   inline Array4D& operator=(const Array4D &arr) {
      for(auto i = 0; i < size_; ++i){
         data_[i] = arr.data_[i];
      }
      return *this;
   }
      
   private:
   T* data_ = nullptr;
   std::size_t x_ = 0;
   std::size_t y_ = 0;
   std::size_t z_ = 0;
   std::size_t q_ = 0; 
   std::size_t size_ = 0;
};