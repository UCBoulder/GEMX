#pragma once

#include <algorithm>
#include <cstddef>

template<typename T>
class CArray4D {
    public:
        CArray4D() {
            data_ = nullptr;
            x_ = 0;
            y_ = 0;
            z_ = 0;
            q_ = 0;
            size_ = 0;
        }

        ~CArray4D() {
            data_ = NULL;
            delete[] data_;
        }

        //x + y*D1 + z*D1*D2 + t*D1*D2*D3
    inline T& operator()(const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t l) {
        return data_[i * (x_*y_*z_) + j * (y_*z_) + (k * z_) + l];
    }

    void resize(size_t xsize, size_t ysize, size_t zsize, size_t qsize){
        size_ = xsize*ysize*zsize*qsize;
        data_ = new T[size_];
        std::fill(data_, data_ + (size_), 0);
    }

    void Clear() {
        std::fill(data_, data_+size_, T{});//T{} instead of 0 because of complex types where 0 is real
    }

    inline CArray4D& operator=(const CArray4D &arr){
        for(auto i = 0; i < size_; ++i){
            data_[i] = arr.data_[i];
        }
        return *this;
    }

    T* start() {
        return data_;
    }

    private:
    T* data_ = nullptr;
    std::size_t size_ = 0;
    std::size_t x_ = 0;
    std::size_t y_ = 0;
    std::size_t z_ = 0;
    std::size_t q_ = 0;
};

template<typename T>
class CArray3D {
public:
    CArray3D() {
        data_ = nullptr;
        x_ = 0;
        y_ = 0;
        z_ = 0;
        size_ = 0;
    }

    ~CArray3D() {
        data_ = NULL;
        delete[] data_;
    }

    void resize(size_t xsize, size_t ysize, size_t zsize) { //int flag (idea for future, allow dynamic re-allocation without deleting previous data)
        x_ = xsize;
        y_ = ysize;
        z_ = zsize;
        if(data_) {
            delete[] data_;
        }
        size_ = x_*y_*z_;
        data_ = new T[xsize * ysize * zsize];
        std::fill(data_, data_ + (xsize * ysize * zsize), 0);
    }

    void CreateArray3D(const T* const data, const std::size_t x, const std::size_t y, const std::size_t z) { //legacy?
       data_ = const_cast<T*>(data);
       x_ = x;
       y_ = y;
       z_ = z;
       size_ = x_ * y_ * z_;
    }

    void Clear() {
        std::fill(data_, data_+size_, T{});//T{} instead of 0 because of complex types where 0 is real
    }

    int getX(){
        return x_;
    }

    int getY(){
        return y_;
    }

    int getZ(){
        return z_;
    }

    size_t getSize() {
        return size_;
    }

    T* start() {
        return data_;
    }

    inline T& operator()(const std::size_t i, const std::size_t j, const std::size_t k) {
        return data_[(i * y_ + j) * z_ + k];
    }

    inline CArray3D& operator=(const CArray3D &arr){
        for(auto i = 0; i < size_; ++i){
            data_[i] = arr.data_[i];
        }
        return *this;
    }


private:
    T* data_ = nullptr;
    std::size_t size_ = 0;
    std::size_t x_ = 0;
    std::size_t y_ = 0;
    std::size_t z_ = 0;
};

template<typename T>
class CArray2D{
    public:
    CArray2D() {
        x_ = 0;
        y_ = 0;
        this->data_ = nullptr;
        size_ = 0;
    }

    ~CArray2D(){
        x_ = 0;
        y_ = 0;
        delete[] this->data_;
    }

    void resize(size_t xsize, size_t ysize) {
        x_ = xsize;
        y_ = ysize;
        if(data_) {
            delete[] data_;
        }
        size_ = x_*y_;
        data_ = new T[size_];
        std::fill(data_, data_ + (size_), 0);
    }

    void Clear() {
        std::fill(data_, data_+size_, T{});
    }

    T* start() {
        return data_;
    }

    inline T& operator()(const std::size_t i, const std::size_t j) {
        return data_[i * y_ + j];
    }

    inline CArray2D& operator=(const CArray2D &arr){
        for(auto i = 0; i < size_; ++i){
            data_[i] = arr.data_[i];
        }
        return *this;
    }

    private:
    T* data_ = nullptr;
    size_t x_ = 0;
    size_t y_ = 0;
    size_t size_ = 0;
};