// The Computer Language Benchmarks Game
// https://salsa.debian.org/benchmarksgame-team/benchmarksgame/
//
// Original C contributed by Sebastien Loisel
// Conversion to C++ by Jon Harrop
// OpenMP parallelize by The Anh Tran
// Add SSE by The Anh Tran
// Additional SSE optimization by Krzysztof Jakubowski
// Converted to AVX by Tomas Wain

// g++ -pipe -O3  -fomit-frame-pointer -march=native -fopenmp -mavx2 \
// ./spectralnorm.cpp-7.c++ -o spectralnorm.gpp-7.c++.o && \
// g++ spectralnorm.gpp-7.c++.o -o spectralnorm.gpp-7.gpp_run -fopenmp

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <sys/sysinfo.h>
#include <sched.h>
#include <omp.h>
#include <immintrin.h>

// compute values of A 2 at a time instead of 1
static inline __m128d EvalA(__m128d i, __m128d j) {
    __m128d ip1   = i + _mm_set1_pd(1.0);
    __m128d ipj   = i + j;
    __m128d ipjp1 = ip1 + j;
    __m128d a     = ipj * ipjp1 * _mm_set1_pd(0.5) + ip1;
    return  a;
}

void EvalATimesU(double *__restrict__ src, double *__restrict__ dst,
                 int begin, int end, int length) {
    for(int i = begin; i < end; i += 2) {
      __m128d src0 = _mm_set1_pd(src[0]);
      __m128d ti   = _mm_set1_pd(i) + _mm_setr_pd(0.0, 1.0);
      __m128d eA   = EvalA(ti, _mm_set1_pd(0.0));
      __m128d sum  = src0 / eA;
      for(int j = 1; j < length; j++) {
        __m128d srcj = _mm_set1_pd(src[j]); 
        __m128d idx = eA + ti + _mm_set1_pd(j);
        eA = idx;
        sum += srcj / idx;
      }
      _mm_storeu_pd(dst + i, sum);
    }
}

void EvalAtTimesU(double *__restrict__ src, double *__restrict__ dst,
                  int begin, int end, int length) {
    for(int i = begin; i < end; i += 2) {
      __m128d src0 = _mm_set1_pd(src[0]);
      __m128d ti   = _mm_set1_pd(i) + _mm_setr_pd(1.0, 2.0);    
      __m128d eAt  = EvalA(_mm_set1_pd(0.0), (ti - _mm_set1_pd(1.0)));
      __m128d sum  = src0 / eAt;      	
      for(int j = 1; j < length; j++) {
        __m128d srcj = _mm_set1_pd(src[j]); 
        __m128d idx = eAt + ti + _mm_set1_pd(j);
        eAt = idx;
        sum += srcj / idx;
      }
      _mm_storeu_pd(dst + i, sum);
    }
}

void EvalAtATimesU(double *src, double *dst, double *tmp,
                   int begin, int end, int N) {
    EvalATimesU (src, tmp, begin, end, N);
    #pragma omp barrier
    EvalAtTimesU(tmp, dst, begin, end, N);
    #pragma omp barrier
}

int GetThreadCount() {
  return get_nprocs();
}

double SpectralNorm(int N) {
  __attribute__((aligned(32))) double u[N+3], v[N+3], tmp[N+3];

    double vBv = 0.0;
    double vv = 0.0;

#pragma omp parallel default(shared) num_threads(GetThreadCount())
    {
        // this block will be executed by NUM_THREADS
        // variable declared in this block is private for each thread
        int threadid = omp_get_thread_num();
        int threadcount = omp_get_num_threads();
        int chunk = N / threadcount;

        // calculate each thread's working range [r1 .. r2) => static schedule
        int begin = threadid * chunk;
        int end = (threadid < (threadcount -1)) ? (begin + chunk) : N;

        for(int i = begin; i < end; i++)
            u[i] = 1.0;
        #pragma omp barrier

        for (int ite = 0; ite < 10; ++ite) {
            EvalAtATimesU(u, v, tmp, begin, end, N);
            EvalAtATimesU(v, u, tmp, begin, end, N);
        }

        double sumvb = 0.0, sumvv = 0.0;
        for (int i = begin; i < end; i++) {
            sumvv += v[i] * v[i];
            sumvb += u[i] * v[i];
        }

        #pragma omp critical
        {
            vBv += sumvb;
            vv += sumvv;
        }
    }
    
    return sqrt(vBv / vv);
}

int main(int argc, char *argv[]) {
    int N = ((argc >= 2) ? atoi(argv[1]) : 5500);
    printf("%.9f\n", SpectralNorm(N));
    return 0;
}
