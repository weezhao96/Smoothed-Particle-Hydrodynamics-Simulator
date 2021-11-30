/* linalg.h
    Wrapper for Fortran BLAS routine.
*/

#pragma once

#define F77NAME(x) x##_

extern "C" {

// BLAS Level 1
    // x <- alpha * x
    void F77NAME(dscal)(const int& n,
                        const double& alpha,const double* x, const int& incx);

    // y <- alpha * x + y
    void F77NAME(daxpy)(const int& n,
                        const double& alpha, const double* x, const int& incx,
                        const double* y, const int& incy);

    // z <- x' * y
    double F77NAME(ddot)(const int& n,
                         const double* x, const int& incx,
                         const double* y, const int& incy);

    // z <- sum(x)
    double F77NAME(dasum)(const int& n, const double* x, const int& incx);

    // y <- x
    void F77NAME(dcopy)(const int& n,
                        const double* x, const int& incx,
                        const double* y, const int& incy);

// BLAS Level 3
    // C <- alpha * A * B + beta * C
    void F77NAME(dgemm)(const char& transa, const char& transb,
                        const int& m, const int& n, const int& k,
                        const double& alpha, const double* A, const int& lda,
                        const double* B, const int& ldb,
                        const double& beta, const double* C, const int& ldc);
}