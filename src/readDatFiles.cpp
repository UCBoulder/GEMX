#include "readDatFiles.hpp"

// Read values from a 1D array .dat file and store into a 1D array in C++
void read1D(std::string fname, double arr[]){
    std::ifstream file;
    std::string line;
    int i = 0;
    file.open(fname);
    while(getline(file, line)) {
        arr[i] = stod(line); 
        ++i;    
    }
    file.close();
}
// Read values from 2D array .dat file and store into 2D array in C++ (bounds inclusive)
void read2D(std::string fname, CArray2D<double> arr, int x, int y){
    int size = (x+1)*(y+1);
    double tempArr[size]; 
    std::ifstream file;
    std::string line;
    std::string num;
    std::string delimeter = "	";
    int i = 0;
    file.open(fname);
    //Break Down lines in file to store values in 1D array
    while(getline(file, line)){
        std::stringstream str(line);
        while(getline(str, num, delimeter[0])){
            tempArr[i] = stod(num);
            i+=1;
        }
    }
    file.close();
    //use vlaues stored in 1D array to translate to 2D array
    int r = 0; //indexing term for 1D array
    for(int i = 0; i <= x; ++i){
        for(int j = 0; j <= y; ++j){
            arr(i,j) = tempArr[r];
            r++;
        }
    }
}