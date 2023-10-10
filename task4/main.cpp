#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>
#include <map>
#include <functional>
#include <chrono>

#include <mpi.h>

class Timer {
public:

    Timer()
    :   _begin(std::chrono::high_resolution_clock::now())
    {}

    double elapsed() const
    {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::microseconds>(end - _begin).count() / 1000000.0;
    }

private:

    std::chrono::high_resolution_clock::time_point _begin;
};

int main(int argc, char** argv) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <A filename> <B filename> <C filename> <mode> [iterations]\n";
        return 1;
    }

    int P_row = std::atoi(argv[1]),
        P_col = std::atoi(argv[2]);

    const char *A_filename = argv[3],
               *B_filename = argv[4],
               *C_filename = argv[5];

    int iterations = 1;
    if (argc > 6)
        iterations = std::atoi(argv[6]);
    
    MPI_Init(NULL, NULL);
    


    return 0;
}