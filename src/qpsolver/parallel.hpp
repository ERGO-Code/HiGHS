#ifndef __PARALLEL_HPP__
#define __PARALLEL_HPP__

#include <algorithm>
#include <thread>
#include <functional>
#include <vector>
#include <mutex>

enum class PARALLELISM_SETTING {
    NONE,
    OMP,
    BUILTIN
};

/// @param[in] nb_elements : size of your for loop
/// @param[in] functor(start, end) :
/// your function processing a sub chunk of the for loop.
/// "start" is the first index to process (included) until the index "end"
/// (excluded)
/// @code
///     for(int i = start; i < end; ++i)
///         computation(i);
/// @endcode
/// @param use_threads : enable / disable threads.
///
///
static
void parallel_for(unsigned nb_elements,
                  std::function<void (int start, int end)> functor,
                  PARALLELISM_SETTING mode = PARALLELISM_SETTING::NONE)
{
    // -------
    unsigned nb_threads_hint = std::thread::hardware_concurrency();
    unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);

    unsigned batch_size = nb_elements / nb_threads;
    unsigned batch_remainder = nb_elements % nb_threads;

    std::vector< std::thread > my_threads(nb_threads);
    int start;
    switch (mode) {
        case PARALLELISM_SETTING::NONE:
            functor( 0, nb_elements );
            break;
        case PARALLELISM_SETTING::BUILTIN:
            for(unsigned i = 0; i < nb_threads; ++i)
            {
                start = i * batch_size;
                my_threads[i] = std::thread(functor, start, start+batch_size);
            }
            start = nb_threads * batch_size;
            functor( start, start+batch_remainder);

            std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));
            break;
        case PARALLELISM_SETTING::OMP:
        #pragma omp parallel for
            for(unsigned i = 0; i < nb_elements; ++i){
                functor( i, i+1 );
            }
            break;
    }
}

static void parallel_for_obo(int nb_elements,
                  std::function<void (int idx)> functor) {
    // -------
    unsigned nb_threads_hint = std::thread::hardware_concurrency();
    unsigned nb_threads = 1;//nb_threads_hint == 0 ? 8 : (nb_threads_hint);

    std::vector< std::thread > my_threads(nb_threads);
    
    int assigned = 0;
    std::vector<int> current (nb_threads, 0);
    std::mutex lock;
    for (int i=0; i<nb_threads; i++) {
        my_threads[i] = std::thread([&](int tid)
        {
            while (true) {
                lock.lock();
                if (assigned < nb_elements) {
                    current[tid] = assigned++;
                    lock.unlock();
                } else {
                    lock.unlock();
                    break;
                }
                
                functor(current[tid]);
            }
        }, i);
    }
    std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));
}

static void parallel_for_frac(int nb_elements,
                  std::function<void (int start, int end)> functor) {
    // -------
    unsigned nb_threads_hint = std::thread::hardware_concurrency();
    unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);

    std::vector< std::thread > my_threads(nb_threads);
    
    int assigned = 0;
    std::vector<int> current_start(nb_threads, 0);
    std::vector<int> current_end(nb_threads, 0);
    std::mutex lock;
    for (int i=0; i<nb_threads; i++) {
        my_threads[i] = std::thread([&](int tid) {
            while (true) {
                lock.lock();
                if (assigned < nb_elements) {
                    int left = nb_elements - assigned;
                    int assign = left / (3* nb_threads) + 1;
                    current_start[tid] = assigned;
                    assigned = std::min(nb_elements, assigned + assign);
                    current_end[tid] = assigned;
                    lock.unlock();
                } else {
                    lock.unlock();
                    break;
                }
                
                functor(current_start[tid], current_end[tid]);
            }
        }, i);
    }
    std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));
}

#endif
