#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>
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
        std::cerr << "Usage: " << argv[0] << " <P_row> <P_col> <A filename> <B filename> <C filename>\n";
        return 1;
    }

    MPI_Init(NULL, NULL);

    int P_row = std::atoi(argv[1]),
        P_col = std::atoi(argv[2]);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    assert(size == P_row * P_col);

    const char *A_filename = argv[3],
               *B_filename = argv[4],
               *C_filename = argv[5];

    // Matrix file
    std::ifstream A_file(A_filename);
    assert(A_file);

    // Matrix dimentions header
    struct {
        int32_t m, n;
    } matrix_header;
    A_file.read((char*)&matrix_header, sizeof(matrix_header));

    // Cut world communicator if there are more processes than matrix elements in the same axis

    P_row = std::min(P_row, matrix_header.m);
    P_col = std::min(P_col, matrix_header.n);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm world;
    MPI_Comm_split(MPI_COMM_WORLD, rank < P_row * P_col, 0, &world);

    // Get rid of out-of-bounds processes
    if (rank >= P_row * P_col) {
        MPI_Finalize();
        return 0;
    }

    MPI_Comm_rank(world, &rank);

    // Get per-axis ranks and row communicators for later row-wise reducing

    int y_rank = rank / P_col, x_rank;
    MPI_Comm row_comm;
    MPI_Comm_split(world, y_rank, 0, &row_comm);
    MPI_Comm_rank(row_comm, &x_rank);

    // Matrix preprocessing //

    // Block dimentions
    int height = matrix_header.m / P_row + (y_rank < matrix_header.m % P_row),
        width  = matrix_header.n / P_col + (x_rank < matrix_header.n % P_col);

    // Original matrix offsets
    int y_begin = matrix_header.m / P_row * y_rank + std::min(y_rank, matrix_header.m % P_row),
        x_begin = matrix_header.n / P_col * x_rank + std::min(x_rank, matrix_header.n % P_col);

    // Block reading
    std::vector<int32_t> A(height * width);
    for (int i = 0; i < height; ++i) {
        A_file.seekg(sizeof(matrix_header) + (matrix_header.n * (y_begin + i) + x_begin) * sizeof(int32_t));
        A_file.read((char*)&A[i * width], width * sizeof(int32_t));
    }

    A_file.close();

    // Vector preprocessing //

    // Vector file
    std::ifstream B_file(B_filename);
    assert(B_file);

    // Vector dimention header
    struct {
        int32_t n;
    } vector_header;
    B_file.read((char*)&vector_header, sizeof(vector_header));
    assert(vector_header.n == matrix_header.n);

    // Vector part reading
    std::vector<int32_t> B(width);
    B_file.seekg(sizeof(vector_header.n) + x_begin * sizeof(int32_t));
    B_file.read((char*)B.data(), B.size() * sizeof(int32_t));

    B_file.close();

    // Result preparation //

    std::vector<int32_t> C_partial(height);

    // Performing matvec //

    Timer timer;

    for (int32_t y = 0; y < height; ++y) {
        int32_t cache = 0;
        for (int32_t x = 0; x < width; ++x) {
            cache += A[y * width + x] * B[x];
        }
        C_partial[y] = cache;
    }

    // Collecting result //

    std::vector<int32_t> C;
    if (x_rank == 0)
        C.resize(height);
        
    MPI_Reduce(C_partial.data(), C.data(), C_partial.size(), MPI_INT32_T, MPI_SUM, 0, row_comm);

    // Time estimation //

    float local_time = timer.elapsed(), max_time;
    MPI_Reduce(&local_time, &max_time, 1, MPI_FLOAT, MPI_MAX, 0, world);
    if (rank == 0)
        std::cout << "Elapsed time: " << max_time << " s\n";

    // Saving results to file //

    std::ofstream C_file;
    if (rank == 0) {
        C_file.open(C_filename, std::ios_base::out | std::ios_base::trunc);
        C_file.write((char*)&matrix_header.m, sizeof(matrix_header.m));
    } else if (x_rank == 0)
        C_file.open(C_filename, std::ios_base::out);

    MPI_Barrier(world);

    if (x_rank == 0) {
        C_file.seekp(sizeof(matrix_header.m) + y_begin * sizeof(int32_t));
        C_file.write((char*)C.data(), C.size() * sizeof(int32_t));
    }

    C_file.close();

    MPI_Finalize();

    return 0;
}
