
#include <omp.h>

#include <iostream>

#include "catch.hpp"
#include "matrix_multiplication.hpp"
#include "parallel/HighsParallel.h"

using namespace highs;

// const bool dev_run = false;

int64_t fib_sequential(const int64_t n) {
  if (n <= 1) return 1;
  return fib_sequential(n - 1) + fib_sequential(n - 2);
}

int64_t fib(const int64_t n) {
  if (n <= 20) return fib_sequential(n);

  int64_t n1;
  parallel::spawn([&]() { n1 = fib(n - 1); });
  int64_t n2 = fib(n - 2);
  parallel::sync();

  // printf("fib(%ld) = %ld + %ld = %ld\n", n, n1, n2, n1 + n2);
  return n1 + n2;
}

int64_t fib_omp(const int64_t n) {
  if (n <= 30) return fib_sequential(n);

  int64_t n1;
#pragma omp task shared(n1)
  n1 = fib_omp(n - 1);

  int64_t n2 = fib_omp(n - 2);

#pragma omp taskwait

  // printf("fib(%ld) = %ld + %ld = %ld\n", n, n1, n2, n1 + n2);
  return n1 + n2;
}

// matrix_multiplication_omp
// reference: https://computing.llnl.gov/tutorials/openMP/samples/C/omp_mm.c
void matrix_multiplication_omp(unsigned nthreads) {
  omp_set_num_threads(nthreads);

  int i, j, k;

#pragma omp parallel for private(i, j)
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; j++) {
      a[i][j] = i + j;
    }
  }

#pragma omp parallel for private(i, j)
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; j++) {
      b[i][j] = i * j;
    }
  }

#pragma omp parallel for private(i, j)
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; j++) {
      c[i][j] = 0;
    }
  }

#pragma omp parallel for private(i, j, k)
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; j++) {
      for (k = 0; k < N; k++) {
        c[i][j] += a[i][k] * b[k][j];
      }
    }
  }

  // int edge;

  //#pragma omp parallel shared(a, b, c, nthreads) private(i, j, k)
  //{
  //  #pragma omp single private(i, j)
  //  for(i = 0; i<N; i++) {
  //    #pragma omp task private(j) firstprivate(i) depend(out: edge)
  //    for (j=0; j<N; j++)
  //      a[i][j]= i+j;
  //  }

  //  #pragma omp single private(i, j)
  //  for(i = 0; i<N; i++) {
  //    #pragma omp task private(j) firstprivate(i) depend(out: edge)
  //    for (j=0; j<N; j++)
  //      b[i][j]= i*j;
  //  }

  //  #pragma omp single private(i, j)
  //  for(i = 0; i<N; i++) {
  //    #pragma omp task private(j) firstprivate(i) depend(out: edge)
  //    for (j=0; j<N; j++)
  //      c[i][j]= 0;
  //  }

  //  #pragma omp single private(i, j)
  //  for(i = 0; i<N; i++) {
  //    #pragma omp task private(j, k) firstprivate(i) depend(in: edge)
  //    for(j=0; j<N; j++) {
  //      for (k=0; k<N; k++) {
  //        c[i][j] += a[i][k] * b[k][j];
  //      }
  //    }
  //  }
  //}

  // std::cout << reduce_sum() << std::endl;
}

void matrix_multiplication_highs(unsigned nthreads) {
  //#pragma omp parallel for private(i, j)
  parallel::for_each(0, N, [&](HighsInt start, HighsInt end) {
    for (int i = start; i < end; ++i) {
      for (int j = 0; j < N; j++) {
        a[i][j] = i + j;
      }
    }
  });

  parallel::for_each(0, N, [&](HighsInt start, HighsInt end) {
    for (int i = start; i < end; ++i) {
      for (int j = 0; j < N; j++) {
        b[i][j] = i * j;
      }
    }
  });

  parallel::for_each(0, N, [&](HighsInt start, HighsInt end) {
    for (int i = start; i < end; ++i) {
      for (int j = 0; j < N; j++) {
        c[i][j] = 0;
      }
    }
  });

  parallel::for_each(0, N, [&](HighsInt start, HighsInt end) {
    for (int i = start; i < end; ++i) {
      for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {
          c[i][j] += a[i][k] * b[k][j];
        }
      }
    }
  });
}

std::chrono::microseconds measure_time_omp(unsigned num_threads) {
  auto beg = std::chrono::high_resolution_clock::now();
  matrix_multiplication_omp(num_threads);
  auto end = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(end - beg);
}

std::chrono::microseconds measure_time_highs(unsigned num_threads) {
  auto beg = std::chrono::high_resolution_clock::now();
  matrix_multiplication_highs(num_threads);
  auto end = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(end - beg);
}

int N = 0;
double **a = nullptr, **b = nullptr, **c = nullptr;

void matrix_multiplication(const std::string& model, const unsigned num_threads,
                           const unsigned num_rounds) {
  std::cout << std::setw(12) << "size" << std::setw(12) << "runtime"
            << std::endl;

  for (int i = 128; i <= 1024; i += 32) {
    N = i;

    allocate_matrix();

    double runtime{0.0};

    for (unsigned j = 0; j < num_rounds; ++j) {
      if (model == "highs") {
        runtime += measure_time_highs(num_threads).count();
      } else if (model == "omp") {
        runtime += measure_time_omp(num_threads).count();
      } else
        assert(false);
    }

    std::cout << std::setw(12) << N << std::setw(12)
              << runtime / num_rounds / 1e3 << std::endl;

    deallocate_matrix();
  }
}

TEST_CASE("MatrixMultHighs", "[parallel]") {
  parallel::initialize_scheduler(16);
  std::cout << "\nhighs workstealing for loop:" << std::endl;
  matrix_multiplication("highs", parallel::num_threads(), 10);
}

TEST_CASE("MatrixMultOmp", "[parallel]") {
  std::cout << "\nomp for loop:" << std::endl;
  matrix_multiplication("omp", 16, 10);
}

TEST_CASE("FibonacciTasksHighs", "[parallel]") {
  parallel::initialize_scheduler(16);
  int64_t result = fib(50);

  REQUIRE(result == 20365011074);
  // fib 46
  // REQUIRE(result == 2971215073);
}

TEST_CASE("FibonacciTasksOmp", "[parallel]") {
  int64_t result;
#pragma omp parallel num_threads(16)
  {
#pragma omp single
    { result = fib_omp(50); }
  }

  REQUIRE(result == 20365011074);
  // fib 46
  // REQUIRE(result == 2971215073);
}
