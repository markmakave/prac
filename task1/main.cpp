#include <iostream>
#include <vector>
#include <thread>
#include <chrono>
#include <mutex>
#include <cmath>

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

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <num_segments> <num_threads>" << std::endl;
        return 1;
    }

    int num_segments = std::atoi(argv[1]), num_threads = std::atoi(argv[2]);
    int iterations = num_segments / num_threads;
    double dx = 1.0 / num_segments;

    double pi;
    std::mutex pi_mutex;

    auto job = [&](int thread_id) {
        double start = thread_id * iterations * dx;

        double local_pi = 0;
        for (int i = 0; i < iterations; ++i)
            local_pi += 4.0 / (1.0 + std::pow(start + dx * i, 2));

        std::lock_guard<std::mutex> lock(pi_mutex);
        pi += local_pi;
    };

    Timer timer;

    enum { ITERATIONS = 1 };
    for (int iteration = 0; iteration < ITERATIONS; ++iteration) {
        pi = 0;

        std::vector<std::thread> thread_pool;
        for (int i = 0; i < num_threads; ++i)
            thread_pool.emplace_back(job, i);

        for (auto& thread : thread_pool)
            thread.join();

        pi *= dx;
    }

    std::cout << pi << "\nElapsed time: " << timer.elapsed() / ITERATIONS << " s\n";

    return 0;
}
