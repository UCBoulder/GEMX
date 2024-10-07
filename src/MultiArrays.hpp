#pragma once

#include <cstddef>

// 2D array class with column-major order
template<typename T>
class Array2D {
public:
    Array2D(const T* const data, const std::size_t x, const std::size_t y)
        : data_(const_cast<T*>(data)), x_(x), y_(y) {}

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
    Array3D(const T* const data, const std::size_t x, const std::size_t y, const std::size_t z)
        : data_(const_cast<T*>(data)), x_(x), y_(y), z_(z) {}

    // Column-major access
    T& operator()(const std::size_t i, const std::size_t j, const std::size_t k) {
        return data_[(k * y_ + j) * x_ + i];  // Column-major order
    }

private:
    T* data_ = nullptr;
    std::size_t x_ = 0;
    std::size_t y_ = 0;
    std::size_t z_ = 0;
};