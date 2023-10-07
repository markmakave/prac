#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>
#include <map>
#include <functional>
#include <chrono>
#include <papi.h>

#define PAPI_CHECK(x) { int ret = x; if (ret != PAPI_OK) std::cerr << "PAPI Error on line " << __LINE__ << ": " << PAPI_strerror(ret) << '\n'; }

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

enum STRATEGY { IJK, IKJ, KIJ, JIK, JKI, KJI };

template <enum STRATEGY>
void matmul(
    const int32_t* A,
    const int32_t* B,
          int32_t* C,
    const int n
);

template <>
void matmul<IJK>(
    const int32_t* A,
    const int32_t* B,
          int32_t* C,
    const int n
) {
    for (int32_t i = 0; i < n; ++i) {
        for (int32_t j = 0; j < n; ++j) {
            int32_t cache = 0;
            for (int32_t k = 0; k < n; ++k) {
                cache += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = cache;
        }
    }
}

template <>
void matmul<IKJ>(
    const int32_t* A,
    const int32_t* B,
          int32_t* C,
    const int n
) {
    std::fill(C, C + n * n, 0);
    for (int32_t i = 0; i < n; ++i) {
        for (int32_t k = 0; k < n; ++k) {
            int32_t cache = A[i * n + k];
            for (int32_t j = 0; j < n; ++j) {
                C[i * n + j] += cache * B[k * n + j];
            }
        }
    }
}

template <>
void matmul<KIJ>(
    const int32_t* A,
    const int32_t* B,
          int32_t* C,
    const int n
) {
    std::fill(C, C + n * n, 0);
    for (int32_t k = 0; k < n; ++k) {
        for (int32_t i = 0; i < n; ++i) {
            int32_t cache = A[i * n + k];
            for (int32_t j = 0; j < n; ++j) {
                C[i * n + j] += cache * B[k * n + j];
            }
        }
    }
}

template <>
void matmul<JIK>(
    const int32_t* A,
    const int32_t* B,
          int32_t* C,
    const int n
) {
    for (int32_t j = 0; j < n; ++j) {
        for (int32_t i = 0; i < n; ++i) {
            int32_t cache = 0;
            for (int32_t k = 0; k < n; ++k) {
                cache += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = cache;
        }
    }
}

template <>
void matmul<JKI>(
    const int32_t* A,
    const int32_t* B,
          int32_t* C,
    const int n
) {
    std::fill(C, C + n * n, 0);
    for (int32_t j = 0; j < n; ++j) {
        for (int32_t k = 0; k < n; ++k) {
            int32_t cache = B[k * n + j];
            for (int32_t i = 0; i < n; ++i) {
                C[i * n + j] += A[i * n + k] * cache;
            }
        }
    }
}

template <>
void matmul<KJI>(
    const int32_t* A,
    const int32_t* B,
          int32_t* C,
    const int n
) {
    std::fill(C, C + n * n, 0);
    for (int32_t k = 0; k < n; ++k) {
        for (int32_t j = 0; j < n; ++j) {
            int32_t cache = B[k * n + j];
            for (int32_t i = 0; i < n; ++i) {
                C[i * n + j] += A[i * n + k] * cache;
            }
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <A filename> <B filename> <C filename> <mode>\n";
        return 1;
    }

    const char *A_filename = argv[1],
               *B_filename = argv[2],
               *C_filename = argv[3],
               *mode = argv[4];
    
    auto read_matrix = [](const char* filename) -> std::vector<int32_t> {
        std::ifstream f(filename);
        assert(f);

        int32_t n;
        f.read((char*)&n, sizeof(n));

        std::vector<int32_t> m(n*n);
        f.read((char*)m.data(), n*n*sizeof(int32_t));

        return m;
    };

    auto write_matrix = [](const std::vector<int32_t>& M, const char* filename) {
        std::ofstream f(filename);
        assert(f);
        
        int32_t n = std::sqrt(M.size());
        f.write((const char*)&n, sizeof(n));
        f.write((const char*)M.data(), M.size() * sizeof(int32_t));
    };

    std::vector<int32_t> A = read_matrix(A_filename),
                         B = read_matrix(B_filename);
    assert(A.size() == B.size());

    std::vector<int32_t> C(A.size());
    std::fill(C.begin(), C.end(), 0);

    int n = std::sqrt(A.size());
    
    std::map<std::string, void(*)(const int32_t*, const int32_t*, int32_t*, const int)> launch{
        { "ijk", matmul<IJK> },
        { "ikj", matmul<IKJ> },
        { "kij", matmul<KIJ> },
        { "jik", matmul<JIK> },
        { "jki", matmul<JKI> },
        { "kji", matmul<KJI> }
    };

    int ret = PAPI_library_init(PAPI_VER_CURRENT);
    assert(ret == PAPI_VER_CURRENT);

    int event_set = PAPI_NULL;
    PAPI_CHECK(PAPI_create_eventset(&event_set));
    PAPI_CHECK(PAPI_add_event(event_set, PAPI_L1_DCM));
    PAPI_CHECK(PAPI_add_event(event_set, PAPI_L2_));

    auto& entry = launch[mode];

    long long cm[2];
    PAPI_CHECK(PAPI_start(event_set));
    entry(A.data(), B.data(), C.data(), n);
    PAPI_CHECK(PAPI_stop(event_set, cm));

    std::cout << "Mode: " << mode << "\nL1 misses: " << cm[0] << "\nL2 misses: " << cm[1] << '\n';

    write_matrix(C, C_filename);

    return 0;
}