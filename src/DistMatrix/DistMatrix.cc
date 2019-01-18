// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "DistMatrix.h"
#include "BlacsContext.h"
#include "MGmol_MPI.h"

#include <cassert>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <string.h>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef SCALAPACK
#include "blacs.h"
#endif

#ifdef __KCC
#define TEMP_DECL /**/
#else
#define TEMP_DECL template <>
#endif

#ifdef ADD_
#define idamax idamax_
#define isamax isamax_
#define pdamax pdamax_
#define psamax psamax_
#define dswap dswap_
#define indxl2g indxl2g_
#endif

extern "C"
{
#ifdef SCALAPACK
    // PBLAS
    void pdswap(
        Pint, double*, Pint, Pint, Pint, Pint, double*, Pint, Pint, Pint, Pint);
    void psswap(
        Pint, float*, Pint, Pint, Pint, Pint, float*, Pint, Pint, Pint, Pint);
    void pdsymm(Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint, Pint, Pint,
        Pdouble, Pint, Pint, Pint, Pdouble, double* const, Pint, Pint, Pint);
    void pssymm(Pchar, Pchar, Pint, Pint, Pfloat, Pfloat, Pint, Pint, Pint,
        Pfloat, Pint, Pint, Pint, Pfloat, float* const, Pint, Pint, Pint);
    void pdgemm(Pchar, Pchar, Pint, Pint, Pint, Pdouble, Pdouble, Pint, Pint,
        Pint, Pdouble, Pint, Pint, Pint, Pdouble, double* const, Pint, Pint,
        Pint);
    void pdgemv(Pchar, Pint, Pint, Pdouble, Pdouble, Pint, Pint, Pint, Pdouble,
        Pint, Pint, Pint, Pint, Pdouble, double* const, Pint, Pint, Pint, Pint);
    void pdsymv(Pchar, Pint, Pdouble, Pdouble, Pint, Pint, Pint, Pdouble, Pint,
        Pint, Pint, Pint, Pdouble, double* const, Pint, Pint, Pint, Pint);
    void psgemm(Pchar, Pchar, Pint, Pint, Pint, Pfloat, Pfloat, Pint, Pint,
        Pint, Pfloat, Pint, Pint, Pint, Pfloat, float* const, Pint, Pint, Pint);
    void psgemv(Pchar, Pint, Pint, Pfloat, Pfloat, Pint, Pint, Pint, Pfloat,
        Pint, Pint, Pint, Pint, Pfloat, float* const, Pint, Pint, Pint, Pint);
    void pssymv(Pchar, Pint, Pfloat, Pfloat, Pint, Pint, Pint, Pfloat, Pint,
        Pint, Pint, Pint, Pfloat, float* const, Pint, Pint, Pint, Pint);
    void pdsyrk(Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint, Pint, Pint,
        Pdouble, double* const, Pint, Pint, Pint);
    void pssyrk(Pchar, Pchar, Pint, Pint, Pfloat, Pfloat, Pint, Pint, Pint,
        Pfloat, float* const, Pint, Pint, Pint);
    void pdtran(Pint, Pint, Pdouble, Pdouble, Pint, Pint, Pint, Pdouble,
        double* const, Pint, Pint, Pint);
    void pstran(Pint, Pint, Pfloat, Pfloat, Pint, Pint, Pint, Pfloat,
        float* const, Pint, Pint, Pint);
    void pdtrmm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint,
        Pint, Pint, double*, Pint, Pint, Pint);
    void pstrmm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pfloat, Pfloat, Pint,
        Pint, Pint, float*, Pint, Pint, Pint);
    void pdtrsm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint,
        Pint, Pint, double*, Pint, Pint, Pint);
    void pstrsm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pfloat, Pfloat, Pint,
        Pint, Pint, float*, Pint, Pint, Pint);
    void pdamax(Pint, double*, int*, double*, Pint, Pint, Pint, Pint);
    void psamax(Pint, float*, int*, float*, Pint, Pint, Pint, Pint);

    // SCALAPACK
    double pdlatra(Pint, Pdouble, Pint, Pint, Pint);
    double pslatra(Pint, Pfloat, Pint, Pint, Pint);
    double pdlaset(
        char*, int*, int*, double*, double*, double*, int*, int*, int*);
    float pslaset(char*, int*, int*, float*, float*, float*, int*, int*, int*);
    double pdlacp3(Pint, Pint, double*, Pint, double*, Pint, Pint, Pint, Pint);
    float pslacp3(Pint, Pint, float*, Pint, float*, Pint, Pint, Pint, Pint);
    double pdlange(Pchar, int*, int*, double*, int*, int*, int*, double*);
    float pslange(Pchar, int*, int*, float*, int*, int*, int*, float*);
    void pdtrtrs(Pchar, Pchar, Pchar, Pint, Pint, Pdouble, Pint, Pint, Pint,
        double*, Pint, Pint, Pint, int*);
    void pstrtrs(Pchar, Pchar, Pchar, Pint, Pint, Pfloat, Pint, Pint, Pint,
        float*, Pint, Pint, Pint, int*);
    void pdgemr2d(Pint, Pint, Pdouble, Pint, Pint, Pint, double* const, Pint,
        Pint, Pint, Pint);
    void psgemr2d(Pint, Pint, Pfloat, Pint, Pint, Pint, float* const, Pint,
        Pint, Pint, Pint);
    void pigemr2d(
        Pint, Pint, Pint, Pint, Pint, Pint, int*, Pint, Pint, Pint, Pint);
    void pdpotrf(Pchar, int*, double*, int*, int*, int*, int*);
    void pspotrf(Pchar, int*, float*, int*, int*, int*, int*);
    void pdpotrs(Pchar, int*, int*, double*, int*, int*, int*, double*, int*,
        int*, int*, int*);
    void pspotrs(Pchar, int*, int*, float*, int*, int*, int*, float*, int*,
        int*, int*, int*);
    void pdgetrf(int*, int*, double*, int*, int*, int*, int*, int*);
    void psgetrf(int*, int*, float*, int*, int*, int*, int*, int*);
    void pdgetrs(Pchar, int*, int*, double*, int*, int*, int*, int*, double*,
        int*, int*, int*, int*);
    void psgetrs(Pchar, int*, int*, float*, int*, int*, int*, int*, float*,
        int*, int*, int*, int*);
    void pdpotri(Pchar, int*, double*, int*, int*, int*, int*);
    void pspotri(Pchar, int*, float*, int*, int*, int*, int*);
    void pdtrtri(Pchar, Pchar, int*, double*, int*, int*, int*, int*);
    void pdpocon(Pchar, int*, double*, int*, int*, int*, double*, double*,
        double*, int*, int*, int*, int*);
    void pspocon(Pchar, int*, float*, int*, int*, int*, float*, float*, float*,
        int*, int*, int*, int*);
    void pdsygst(Pint, Pchar, Pint, double*, Pint, Pint, Pint, Pdouble, Pint,
        Pint, Pint, double*, int*);
    void pssygst(Pint, Pchar, Pint, float*, Pint, Pint, Pint, Pfloat, Pint,
        Pint, Pint, float*, int*);
    void pdsyev(Pchar, Pchar, int*, double*, int*, int*, int*, double*, double*,
        int*, int*, int*, double*, int*, int*);
    void pssyev(Pchar, Pchar, int*, float*, int*, int*, int*, float*, float*,
        int*, int*, int*, float*, int*, int*);
    void pdgesvd(Pchar, Pchar, int*, int*, double*, int*, int*, int*, double*,
        double*, int*, int*, int*, double*, int*, int*, int*, double*, int*,
        int*);
    void psgesvd(Pchar, Pchar, int*, int*, float*, int*, int*, int*, float*,
        float*, int*, int*, int*, float*, int*, int*, int*, float*, int*, int*);
    // SCALAPACK TOOLS
    int numroc(Pint, Pint, Pint, Pint, Pint);
    int indxl2g(Pint, Pint, Pint, Pint, Pint);
#endif
    // BLAS1
    void daxpy(Pint, Pdouble, Pdouble, Pint, double*, Pint);
    void dcopy(Pint, Pdouble, Pint, double* const, Pint);
    //  int idamax(Pint,double*,Pint);
    void dswap(Pint, double*, Pint, double*, Pint);

    // BLAS2
    void dgemv(Pchar, Pint, Pint, Pdouble, Pdouble, Pint, Pdouble, Pint,
        Pdouble, double*, Pint);

    void dsymv(Pchar, Pint, Pdouble, Pdouble, Pint, Pdouble, Pint, Pdouble,
        double*, Pint);

    // BLAS3
    void dsymm(Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint, Pdouble, Pint,
        Pdouble, double*, Pint);
    void dgemm(Pchar, Pchar, Pint, Pint, Pint, Pdouble, Pdouble, Pint, Pdouble,
        Pint, Pdouble, double*, Pint);
    void dsyrk(Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint, Pdouble,
        double* const, Pint);
    void dtrmm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint,
        double* const, Pint);
    void dtrsm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint,
        double* const, Pint);
    // LAPACK
    double dlange(Pchar, int*, int*, double*, int*, double*);
    void dtrtrs(
        Pchar, Pchar, Pchar, Pint, Pint, Pdouble, Pint, double*, Pint, int*);
    void dpotrf(Pchar, int*, double*, int*, int*);
    void dpotrs(Pchar, int*, int*, double*, int*, double*, int*, int*);
    void dgetrf(int*, int*, double*, int*, int*, int*);
    void dgetrs(Pchar, int*, int*, double*, int*, int*, double*, int*, int*);
    void dpotri(Pchar, int*, double*, int*, int*);
    void dtrtri(Pchar, Pchar, int*, double*, int*, int*);
    void dpocon(
        Pchar, int*, double*, int*, double*, double*, double*, int*, int*);
    void dsygst(int*, Pchar, int*, double*, int*, Pdouble, Pint, int*);
    void dsyev(Pchar, Pchar, int*, double*, int*, double*, double*, int*, int*);
    void dgesvd(Pchar, Pchar, int*, int*, double*, int*, double*, double*, int*,
        double*, int*, double*, int*, int*);
}

namespace dist_matrix
{

#ifndef SCALAPACK
int numroc(int* a, int* b, int* c, int* d, int* e) { return *a; }
int indxl2g(Pint indxloc, Pint nb, Pint iproc, Pint isrcproc, Pint nprocs)
{
    return *indxloc;
}
#endif

template <class T>
DistMatrix<T>::DistMatrix(const string& name)
    : object_name_(name),
      comm_global_(default_bc_->comm_global()),
      bc_(*default_bc_),
      m_(0),
      n_(0),
      mb_(0),
      nb_(0)
{
    assert(default_bc_ != 0);
}

template <class T>
DistMatrix<T>::DistMatrix(const string& name, const BlacsContext& bc)
    : object_name_(name),
      comm_global_(bc.comm_global()),
      bc_(bc),
      m_(0),
      n_(0),
      mb_(0),
      nb_(0)
{
}

// Construct a DistMatrix of dimensions m,n
template <class T>
DistMatrix<T>::DistMatrix(
    const string& name, const BlacsContext& bc, const int m, const int n)
    : object_name_(name), comm_global_(bc.comm_global()), bc_(bc)
{
    resize(m, n, distmatrix_def_block_size_, distmatrix_def_block_size_);
}

template <class T>
DistMatrix<T>::DistMatrix(const string& name, const int m, const int n)
    : object_name_(name),
      comm_global_(default_bc_->comm_global()),
      bc_(*default_bc_)
{
    resize(m, n, distmatrix_def_block_size_, distmatrix_def_block_size_);
}

// Construct a DistMatrix of dimensions m,n
template <class T>
DistMatrix<T>::DistMatrix(const string& name, const BlacsContext& bc,
    const int m, const int n, const int mb, const int nb)
    : object_name_(name), comm_global_(bc.comm_global()), bc_(bc)
{
    resize(m, n, mb, nb);
}

// Construct a DistMatrix of dimensions m,n
template <class T>
DistMatrix<T>::DistMatrix(
    const string& name, const int m, const int n, const int mb, const int nb)
    : object_name_(name),
      comm_global_(default_bc_->comm_global()),
      bc_(*default_bc_)
{
    resize(m, n, mb, nb);
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
DistMatrix<T>::DistMatrix(const DistMatrix<T>& a)
    : bc_(a.bc_), comm_global_(a.bc_.comm_global())
{
    object_name_ = a.object_name_ + "_duplicate";
    resize(a.m(), a.n(), a.mb(), a.nb());
    *this = a;
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::resize(const int m, const int n, const int mb, const int nb)
{
    assert(m > 0);
    assert(n > 0);
    assert(mb > 0);
    assert(nb > 0);
    m_ = m;
    n_ = n;
#ifdef SCALAPACK
    mb_ = min(mb, m);
    nb_ = min(nb, n);
#else
    mb_ = m;
    nb_ = n;
#endif
    ictxt_       = bc_.ictxt();
    nprow_       = bc_.nprow();
    npcol_       = bc_.npcol();
    myrow_       = bc_.myrow();
    mycol_       = bc_.mycol();
    active_      = myrow_ >= 0;
    int isrcproc = 0;
    mloc_        = numroc(&m_, &mb_, &myrow_, &isrcproc, &nprow_);
    nloc_        = numroc(&n_, &nb_, &mycol_, &isrcproc, &npcol_);

    // set leading dimension of val array to mloc_;
    lld_ = mloc_;
    if (lld_ == 0) lld_ = 1;

    // total and local number of blocks
    mblocks_ = 0;
    nblocks_ = 0;
    if (active_)
    {
        mblocks_ = (mloc_ + mb_ - 1) / mb_;
        nblocks_ = (nloc_ + nb_ - 1) / nb_;
    }

    m_incomplete_ = mloc_ % mb_ != 0;
    n_incomplete_ = nloc_ % nb_ != 0;

    desc_[0] = 1;
    desc_[1] = ictxt_;
    desc_[2] = m_;
    desc_[3] = n_;
    desc_[4] = mb_;
    desc_[5] = nb_;
    desc_[6] = 0;
    desc_[7] = 0;
    desc_[8] = lld_;

    // allocate local space
    if (active_)
    {
        size_ = mloc_ * nloc_;
        val_  = new T[size_];
    }
    else
    {
        size_ = 0;
        val_  = NULL;
    }
    clear();
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::clear(void)
{
    if (active_) memset(&val_[0], 0, size_ * sizeof(T));
}

////////////////////////////////////////////////////////////////////////////////
// Compute the global index of a distributed matrix entry pointed to by the
// local index of the process
template <class T>
int DistMatrix<T>::indxl2grow(const int indxloc) const
{
    const int isrcproc = 0;
    const int findx    = indxloc + 1;
    return indxl2g(&findx, &mb_, &myrow_, &isrcproc, &nprow_) - 1;
}

template <class T>
int DistMatrix<T>::indxl2gcol(const int indxloc) const
{
    const int isrcproc = 0;
    const int findx    = indxloc + 1;
    return indxl2g(&findx, &nb_, &mycol_, &isrcproc, &npcol_) - 1;
}

////////////////////////////////////////////////////////////////////////////////
// set diagonal of local blocks of DistMatrix to global diag_values
template <class T>
void DistMatrix<T>::setDiagonal(const vector<T>& diag_values)
{
    assert((int)diag_values.size() == m_);

    if (active_)
        for (int i = 0; i < mloc_; i++)
        {
            const int ig = indxl2grow(i);
            for (int j = 0; j < nloc_; j++)
            {
                const int jg = indxl2gcol(j);
                if (ig == jg) val_[i + j * mloc_] = diag_values[ig];
            }
        }
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
double DistMatrix<T>::dot(const DistMatrix<T>& x) const
{
    assert(ictxt_ == x.ictxt());
    double sum  = 0.;
    double tsum = 0.;
    if (active_)
    {
        assert(m_ == x.m());
        assert(n_ == x.n());
        assert(mloc_ == x.mloc());
        assert(nloc_ == x.nloc());
        //    int ione=1;
        assert(x.size_ == size_);
        tsum = MPdot(size_, &val_[0], &x.val_[0]);
        //    tsum=ddot(&size_, &val_[0], &ione, &x.val_[0], &ione);
    }
#ifdef SCALAPACK
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&tsum, &sum, 1, MPI_SUM);
//  MPI_Allreduce(&tsum, &sum, 1, MPI_DOUBLE, MPI_SUM, comm_global_ );
#else
    sum = tsum;
#endif
    return sum;
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
double DistMatrix<T>::nrm2() const
{
    double dott = this->dot(*this);
    return sqrt(dott);
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::scal(const double alpha)
{
    if (alpha == 1.0) return;
    if (alpha == 0.0)
    {
        clear();
        return;
    }
    //  int     ione = 1;
    if (active_) MPscal(size_, alpha, &val_[0]);
    //  if ( active_ ) dscal(&size_, &alpha, &val_[0], &ione);
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::axpy(const double alpha, const DistMatrix<T>& x)
{
    assert(ictxt_ == x.ictxt());
    //  int ione=1;
    assert(m_ == x.m());
    assert(n_ == x.n());
    assert(mloc_ == x.mloc_);
    assert(nloc_ == x.nloc_);
    assert(size_ == x.size_);
    if (active_) MPaxpy(size_, alpha, &x.val_[0], &val_[0]);
    //  if( active_ ) daxpy(&size_, &alpha, &x.val_[0], &ione, &val_[0], &ione);
}

////////////////////////////////////////////////////////////////////////////////
// identity: initialize matrix to identity
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::identity(void)
{
#ifdef SCALAPACK
    char uplo    = 'n';
    double alpha = 0.0;
    double beta  = 1.0;
    int ione     = 1;
    if (active_)
        pdlaset(&uplo, &m_, &n_, &alpha, &beta, &val_[0], &ione, &ione, desc_);
#else
    clear();
    for (int i = 0; i < min(m_, n_); i++)
        val_[i + mloc_ * i] = 1.0;
#endif
}

TEMP_DECL
void DistMatrix<float>::identity(void)
{
#ifdef SCALAPACK
    char uplo   = 'n';
    float alpha = 0.0;
    float beta  = 1.0;
    int ione    = 1;
    if (active_)
        pslaset(&uplo, &m_, &n_, &alpha, &beta, &val_[0], &ione, &ione, desc_);
#else
    clear();
    for (int i = 0; i < min(m_, n_); i++)
        val_[i + mloc_ * i] = 1.0;
#endif
}

////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
DistMatrix<double>& DistMatrix<double>::operator=(const DistMatrix<double>& src)
{
    if (this == &src) return *this;

    // cout << " DistMatrix<double>::operator= for "<<object_name_<<" and
    // "<<a.object_name_ << endl;
    if (src.m() != m_ || src.n() != n_)
    {
        cerr << " DistMatrix<double>::operator=: size inconsistency" << endl;
        exit(1);
    }

#ifdef SCALAPACK
    if (src.mb_ == mb_ && src.nb_ == nb_ && src.bc_ == bc_)
#endif
    {
        // simple copy
        if (active_)
        {
            for (int i = 0; i < 9; i++)
            {
                // cout << " desc_[" << i << "] = " << desc_[i] << " / "
                //     << " src.desc_[" << i << "] = " << src.desc_[i] << endl;
                assert(desc_[i] == src.desc_[i]);
            }
            // cout << " DistMatrix<double>::operator=: simple copy" << endl;
            memcpy(&val_[0], &src.val_[0], size_ * sizeof(double));
        }
    }
#ifdef SCALAPACK
    else
    {
        // redistribute using function pdgemr2d
        // (see file REDIST/SRC/pdgemr.c in ScaLAPACK)
        int ione = 1;
        // get a context englobing at least all processors included
        // in either *this context or src context
        BlacsContext bc(comm_global_);
        int gictxt = bc.ictxt();
        pdgemr2d_tm_.start();
        pdgemr2d(&m_, &n_, src.val_, &ione, &ione, src.desc(), &val_[0], &ione,
            &ione, desc_, &gictxt);
        pdgemr2d_tm_.stop();
    }
#endif
    return *this;
}

TEMP_DECL
DistMatrix<float>& DistMatrix<float>::operator=(const DistMatrix<float>& src)
{
    if (this == &src) return *this;

    // cout << " DistMatrix<float>::operator= for "<<object_name_<<" and
    // "<<a.object_name_ << endl;
    if (src.m() != m_ || src.n() != n_)
    {
        cerr << " DistMatrix<float>::operator=: size inconsistency" << endl;
        exit(1);
    }

#ifdef SCALAPACK
    if (src.mb_ == mb_ && src.nb_ == nb_ && src.bc_ == bc_)
#endif
    {
        // simple copy
        if (active_)
        {
            for (int i = 0; i < 9; i++)
            {
                // cout << " desc_[" << i << "] = " << desc_[i] << " / "
                //     << " src.desc_[" << i << "] = " << src.desc_[i] << endl;
                assert(desc_[i] == src.desc_[i]);
            }
            // cout << " DistMatrix<double>::operator=: simple copy" << endl;
            memcpy(&val_[0], &src.val_[0], size_ * sizeof(float));
        }
    }
#ifdef SCALAPACK
    else
    {
        // redistribute using function pdgemr2d
        // (see file REDIST/SRC/pdgemr.c in ScaLAPACK)
        int ione = 1;
        // get a context englobing at least all processors included
        // in either *this context or src context
        BlacsContext bc(comm_global_);
        int gictxt = bc.ictxt();
        pdgemr2d_tm_.start();
        psgemr2d(&m_, &n_, src.val_, &ione, &ione, src.desc(), &val_[0], &ione,
            &ione, desc_, &gictxt);
        pdgemr2d_tm_.stop();
    }
#endif
    return *this;
}
////////////////////////////////////////////////////////////////////////////////
// Copy a distributed matrix src into another in a block at position (ib,jb)
TEMP_DECL
DistMatrix<double>& DistMatrix<double>::assign(
    const DistMatrix<double>& src, const int ib, const int jb)
{
    if (this == &src) return *this;

    int m = src.m();
    int n = src.n();

    if (ib + m > m_ || jb + n > n_)
    {
        cerr << " DistMatrix<double>::assign: size inconsistency" << endl;
        exit(1);
    }

#ifdef SCALAPACK
    if (src.bc_ == bc_ && src.m_ == m_ && src.n_ == n_ && src.mb_ == mb_
        && src.nb_ == nb_)
#endif
    {
        // simple copy
        if (active_)
        {
            for (int j = 0; j < n; j++)
            {
                memcpy(&val_[(j + jb) * m_ + ib], &src.val_[j * m],
                    m * sizeof(double));
            }
        }
    }
#ifdef SCALAPACK
    else
    {
        // redistribute
        // get global context index
        int pib  = ib + 1;
        int pjb  = jb + 1;
        int ione = 1;
        BlacsContext bc(comm_global_);
        int gictxt = bc.ictxt();
        assert(m + ib <= m_);
        assert(n + jb <= n_);
        pdgemr2d(&m, &n, src.val_, &ione, &ione, src.desc(), val_, &pib, &pjb,
            desc_, &gictxt);
    }
#endif
    return *this;
}

TEMP_DECL
DistMatrix<float>& DistMatrix<float>::assign(
    const DistMatrix<float>& src, const int ib, const int jb)
{
    if (this == &src) return *this;

    int m = src.m();
    int n = src.n();

    if (ib + m > m_ || jb + n > n_)
    {
        cerr << " DistMatrix<float>::assign: size inconsistency" << endl;
        exit(1);
    }

#ifdef SCALAPACK
    if (src.bc_ == bc_ && src.m_ == m_ && src.n_ == n_ && src.mb_ == mb_
        && src.nb_ == nb_)
#endif
    {
        // simple copy
        if (active_)
        {
            for (int j = 0; j < n; j++)
            {
                memcpy(&val_[(j + jb) * m_ + ib], &src.val_[j * m],
                    m * sizeof(float));
            }
        }
    }
#ifdef SCALAPACK
    else
    {
        // redistribute
        // get global context index
        int pib  = ib + 1;
        int pjb  = jb + 1;
        int ione = 1;
        BlacsContext bc(comm_global_);
        int gictxt = bc.ictxt();
        assert(m + ib <= m_);
        assert(n + jb <= n_);
        psgemr2d(&m, &n, src.val_, &ione, &ione, src.desc(), val_, &pib, &pjb,
            desc_, &gictxt);
    }
#endif
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
// set value of diagonal or off-diagonal elements to a constant
// uplo=='u': set strictly upper part to x
// uplo=='l': set strictly lower part to x
// uplo=='d': set diagonal to x
////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::setVal(const char uplo, const T xx)
{
    if (active_)
    {
        if (uplo == 'l' || uplo == 'L')
        {
            // initialize strictly lower part
            for (int li = 0; li < mblocks_; li++)
            {
                for (int lj = 0; lj < nblocks_; lj++)
                {
                    for (int ii = 0; ii < mbs(li); ii++)
                    {
                        for (int jj = 0; jj < nbs(lj); jj++)
                        {
                            if (i(li, ii) > j(lj, jj))
                                val_[(ii + li * mb_) + (jj + lj * nb_) * mloc_]
                                    = xx;
                        }
                    }
                }
            }
        }
        else if (uplo == 'u' || uplo == 'U')
        {
            // initialize strictly upper part
            for (int li = 0; li < mblocks_; li++)
            {
                for (int lj = 0; lj < nblocks_; lj++)
                {
                    for (int ii = 0; ii < mbs(li); ii++)
                    {
                        for (int jj = 0; jj < nbs(lj); jj++)
                        {
                            if (i(li, ii) < j(lj, jj))
                                val_[(ii + li * mb_) + (jj + lj * nb_) * mloc_]
                                    = xx;
                        }
                    }
                }
            }
        }
        else if (uplo == 'd' || uplo == 'D')
        {
            // initialize diagonal elements
            if (active_)
            {
                // loop through all local blocks (ll,mm)
                for (int ll = 0; ll < mblocks(); ll++)
                {
                    for (int mm = 0; mm < nblocks(); mm++)
                    {
                        // check if block (ll,mm) has diagonal elements
                        int imin = i(ll, 0);
                        int imax = imin + mbs(ll) - 1;
                        int jmin = j(mm, 0);
                        int jmax = jmin + nbs(mm) - 1;
                        // cout << " process (" << myrow_ << "," << mycol_ <<
                        // ")"
                        // << " block (" << ll << "," << mm << ")"
                        // << " imin/imax=" << imin << "/" << imax
                        // << " jmin/jmax=" << jmin << "/" << jmax << endl;

                        if ((imin <= jmax) && (imax >= jmin))
                        {
                            // block (ll,mm) holds diagonal elements
                            int idiagmin = max(imin, jmin);
                            int idiagmax = min(imax, jmax);

                            // cout << " process (" << myrow_ << "," << mycol_
                            // << ")"
                            // << " holds diagonal elements " << idiagmin << "
                            // to " << idiagmax << " in block (" << ll << "," <<
                            // mm << ")" << endl;

                            for (int ii = idiagmin; ii <= idiagmax; ii++)
                            {
                                // access element (ii,ii)
                                int jj                  = ii;
                                int iii                 = ll * mb_ + x(ii);
                                int jjj                 = mm * nb_ + y(jj);
                                val_[iii + mloc_ * jjj] = xx;
                            }
                        }
                    }
                }
            }
        }
        else
        {
            cerr << " DistMatrix<double>::set: invalid argument" << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// matrix-vector multiplication
// this = alpha*op(A)*b+beta*this
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::gemv(const char transa, const double alpha,
    const DistMatrix<double>& a, const DistMatrix<double>& b, const double beta)
{
    assert(ictxt_ == a.ictxt());
    assert(ictxt_ == b.ictxt());

    if (active_)
    {
        int m, n;
        if (transa == 'N' || transa == 'n')
        {
            m = a.m();
            n = a.n();
            assert(a.m() == m_);
        }
        else
        {
            m = a.n();
            n = a.m();
            assert(a.n() == m_);
        }

#ifdef SCALAPACK
        int ione = 1;
        pdgemv(&transa, &m, &n, &alpha, &a.val_[0], &ione, &ione, a.desc(),
            &b.val_[0], &ione, &ione, b.desc(), &ione, &beta, &val_[0], &ione,
            &ione, desc_, &ione);
#else
        int ione = 1;
        dgemv(&transa, &m, &n, &alpha, &a.val_[0], &a.lld_, &b.val_[0], &ione,
            &beta, &val_[0], &ione);
#endif
    }
}

TEMP_DECL
void DistMatrix<float>::gemv(const char transa, const float alpha,
    const DistMatrix<float>& a, const DistMatrix<float>& b, const float beta)
{
    assert(ictxt_ == a.ictxt());
    assert(ictxt_ == b.ictxt());

    if (active_)
    {
        int m, n;
        if (transa == 'N' || transa == 'n')
        {
            m = a.m();
            n = a.n();
            assert(a.m() == m_);
        }
        else
        {
            m = a.n();
            n = a.m();
            assert(a.n() == m_);
        }

#ifdef SCALAPACK
        int ione = 1;
        psgemv(&transa, &m, &n, &alpha, &a.val_[0], &ione, &ione, a.desc(),
            &b.val_[0], &ione, &ione, b.desc(), &ione, &beta, &val_[0], &ione,
            &ione, desc_, &ione);
#else
        int ione = 1;
        sgemv(&transa, &m, &n, &alpha, &a.val_[0], &a.lld_, &b.val_[0], &ione,
            &beta, &val_[0], &ione);
#endif
    }
}

TEMP_DECL
void DistMatrix<double>::symv(const char uplo, const double alpha,
    const DistMatrix<double>& a, const DistMatrix<double>& b, const double beta)
{
    assert(ictxt_ == a.ictxt());
    assert(ictxt_ == b.ictxt());

    if (active_)
    {
        assert(a.n() == m_);

#ifdef SCALAPACK
        int ione = 1;
        pdsymv(&uplo, &m_, &alpha, &a.val_[0], &ione, &ione, a.desc(),
            &b.val_[0], &ione, &ione, b.desc(), &ione, &beta, &val_[0], &ione,
            &ione, desc_, &ione);
#else
        int ione = 1;
        dsymv(&uplo, &m_, &alpha, &a.val_[0], &a.lld_, &b.val_[0], &ione, &beta,
            &val_[0], &ione);
#endif
    }
}

TEMP_DECL
void DistMatrix<float>::symv(const char uplo, const float alpha,
    const DistMatrix<float>& a, const DistMatrix<float>& b, const float beta)
{
    assert(ictxt_ == a.ictxt());
    assert(ictxt_ == b.ictxt());

    if (active_)
    {
        assert(a.n() == m_);

#ifdef SCALAPACK
        int ione = 1;
        pssymv(&uplo, &m_, &alpha, &a.val_[0], &ione, &ione, a.desc(),
            &b.val_[0], &ione, &ione, b.desc(), &ione, &beta, &val_[0], &ione,
            &ione, desc_, &ione);
#else
        int ione = 1;
        ssymv(&uplo, &m_, &alpha, &a.val_[0], &a.lld_, &b.val_[0], &ione, &beta,
            &val_[0], &ione);
#endif
    }
}

////////////////////////////////////////////////////////////////////////////////
// matrix multiplication
// this = alpha*op(A)*op(B)+beta*this
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::gemm(const char transa, const char transb,
    const double alpha, const DistMatrix<double>& a,
    const DistMatrix<double>& b, const double beta)
{
    assert(ictxt_ == a.ictxt());
    assert(ictxt_ == b.ictxt());

    if (active_)
    {
        int m, n, k;
        if (transa == 'N' || transa == 'n')
        {
            m = a.m();
            k = a.n();
            assert(a.m() == m_);
        }
        else
        {
            m = a.n();
            k = a.m();
            assert(a.n() == m_);
        }
        if (transb == 'N' || transb == 'n')
        {
            n = b.n();
            assert(k == b.m());
        }
        else
        {
            n = b.m();
            assert(k == b.n());
        }

#ifdef SCALAPACK
        int ione = 1;
        pdgemm(&transa, &transb, &m, &n, &k, &alpha, &a.val_[0], &ione, &ione,
            a.desc(), &b.val_[0], &ione, &ione, b.desc(), &beta, &val_[0],
            &ione, &ione, desc_);
#else
        dgemm(&transa, &transb, &m, &n, &k, &alpha, &a.val_[0], &a.lld_,
            &b.val_[0], &b.lld_, &beta, &val_[0], &lld_);
#endif
    }
}

TEMP_DECL
void DistMatrix<float>::gemm(const char transa, const char transb,
    const float alpha, const DistMatrix<float>& a, const DistMatrix<float>& b,
    const float beta)
{
    assert(ictxt_ == a.ictxt());
    assert(ictxt_ == b.ictxt());

    if (active_)
    {
        int m, n, k;
        if (transa == 'N' || transa == 'n')
        {
            m = a.m();
            k = a.n();
            assert(a.m() == m_);
        }
        else
        {
            m = a.n();
            k = a.m();
            assert(a.n() == m_);
        }
        if (transb == 'N' || transb == 'n')
        {
            n = b.n();
            assert(k == b.m());
        }
        else
        {
            n = b.m();
            assert(k == b.n());
        }

#ifdef SCALAPACK
        int ione = 1;
        psgemm(&transa, &transb, &m, &n, &k, &alpha, &a.val_[0], &ione, &ione,
            a.desc(), &b.val_[0], &ione, &ione, b.desc(), &beta, &val_[0],
            &ione, &ione, desc_);
#else
        sgemm(&transa, &transb, &m, &n, &k, &alpha, &a.val_[0], &a.lld_,
            &b.val_[0], &b.lld_, &beta, &val_[0], &lld_);
#endif
    }
}
////////////////////////////////////////////////////////////////////////////////
// symmetrix_matrix * matrix multiplication
// this = beta * this + alpha * a * b
// this = beta * this + alpha * b * a
// with "a" symetric
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::symm(const char side, const char uplo,
    const double alpha, const DistMatrix<double>& a,
    const DistMatrix<double>& b, const double beta)
{
    assert(ictxt_ == a.ictxt());
    assert(ictxt_ == b.ictxt());

    if (active_)
    {
        assert(a.n() == a.m());
        if (side == 'L' || side == 'l')
        {
            assert(m_ == a.m());
            assert(n_ == b.n());
            assert(a.n() == b.m());
        }
        else
        {
            assert(n_ == a.m());
            assert(m_ == b.n());
            assert(a.m() == b.n());
        }

#ifdef SCALAPACK
        int ione = 1;
        pdsymm(&side, &uplo, &m_, &n_, &alpha, &a.val_[0], &ione, &ione,
            a.desc(), &b.val_[0], &ione, &ione, b.desc(), &beta, &val_[0],
            &ione, &ione, desc_);
#else
        dsymm(&side, &uplo, &m_, &n_, &alpha, &a.val_[0], &a.lld_, &b.val_[0],
            &b.lld_, &beta, &val_[0], &lld_);
#endif
    }
}

TEMP_DECL
void DistMatrix<float>::symm(const char side, const char uplo,
    const float alpha, const DistMatrix<float>& a, const DistMatrix<float>& b,
    const float beta)
{
    assert(ictxt_ == a.ictxt());
    assert(ictxt_ == b.ictxt());

    if (active_)
    {
        assert(a.n() == a.m());
        if (side == 'L' || side == 'l')
        {
            assert(m_ == a.m());
            assert(n_ == b.n());
            assert(a.n() == b.m());
        }
        else
        {
            assert(n_ == a.m());
            assert(m_ == b.n());
            assert(a.m() == b.n());
        }

#ifdef SCALAPACK
        int ione = 1;
        pssymm(&side, &uplo, &m_, &n_, &alpha, &a.val_[0], &ione, &ione,
            a.desc(), &b.val_[0], &ione, &ione, b.desc(), &beta, &val_[0],
            &ione, &ione, desc_);
#else
        ssymm(&side, &uplo, &m_, &n_, &alpha, &a.val_[0], &a.lld_, &b.val_[0],
            &b.lld_, &beta, &val_[0], &lld_);
#endif
    }
}

////////////////////////////////////////////////////////////////////////////////
// Compute a matrix-matrix product for a real triangular
// matrix or its transpose.
// *this = alpha op(A) * *this    if side=='l'
// *this = alpha * *this * op(A)  if side=='r'
// where op(A) = A or trans(A)
// alpha is a scalar, *this is an m by n matrix, and A is a unit or non-unit,
// upper- or lower-triangular matrix.
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::trmm(const char side, const char uplo,
    const char trans, const char diag, const double alpha,
    const DistMatrix<double>& a)
{
    if (active_)
    {
        assert(a.m_ == a.n_);
        if (side == 'L' || side == 'l')
        {
            assert(a.n_ == m_);
        }
        else
        {
            assert(a.n_ == n_);
        }
#ifdef SCALAPACK
        int ione = 1;
        pdtrmm(&side, &uplo, &trans, &diag, &m_, &n_, &alpha, &a.val_[0], &ione,
            &ione, a.desc_, &val_[0], &ione, &ione, desc_);
#else
        dtrmm(&side, &uplo, &trans, &diag, &m_, &n_, &alpha, &a.val_[0], &a.m_,
            &val_[0], &m_);
#endif
    }
}

TEMP_DECL
void DistMatrix<float>::trmm(const char side, const char uplo, const char trans,
    const char diag, const float alpha, const DistMatrix<float>& a)
{
    if (active_)
    {
        assert(a.m_ == a.n_);
        if (side == 'L' || side == 'l')
        {
            assert(a.n_ == m_);
        }
        else
        {
            assert(a.n_ == n_);
        }
#ifdef SCALAPACK
        int ione = 1;
        pstrmm(&side, &uplo, &trans, &diag, &m_, &n_, &alpha, &a.val_[0], &ione,
            &ione, a.desc_, &val_[0], &ione, &ione, desc_);
#else
        strmm(&side, &uplo, &trans, &diag, &m_, &n_, &alpha, &a.val_[0], &a.m_,
            &val_[0], &m_);
#endif
    }
}
////////////////////////////////////////////////////////////////////////////////
// solve one of the matrix equations
// op( A )*X = alpha*(*this), or X*op( A ) = alpha*(*this)
// where alpha is a scalar, X and (*this) are m by n matrices, A is a unit,
// or non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
// op( A ) = A   or   op( A ) = A'.
//
// The matrix X is overwritten on (*this).
//
// op( A )*X = alpha*(*this)     if side=='l'
// X*op( A ) = alpha*(*this)     if side=='r'
// where op(A) = A or trans(A)
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::trsm(const char side, const char uplo,
    const char trans, const char diag, const double alpha,
    const DistMatrix<double>& a)
{
    if (active_)
    {
        assert(a.m_ == a.n_);
        if (side == 'L' || side == 'l')
        {
            assert(a.n_ == m_);
        }
        else
        {
            assert(a.n_ == n_);
        }
#ifdef SCALAPACK
        int ione = 1;
        pdtrsm(&side, &uplo, &trans, &diag, &m_, &n_, &alpha, &a.val_[0], &ione,
            &ione, a.desc_, &val_[0], &ione, &ione, desc_);
#else
        dtrsm(&side, &uplo, &trans, &diag, &m_, &n_, &alpha, &a.val_[0], &a.m_,
            &val_[0], &m_);
#endif
    }
}

TEMP_DECL
void DistMatrix<float>::trsm(const char side, const char uplo, const char trans,
    const char diag, const float alpha, const DistMatrix<float>& a)
{
    if (active_)
    {
        assert(a.m_ == a.n_);
        if (side == 'L' || side == 'l')
        {
            assert(a.n_ == m_);
        }
        else
        {
            assert(a.n_ == n_);
        }
#ifdef SCALAPACK
        int ione = 1;
        pstrsm(&side, &uplo, &trans, &diag, &m_, &n_, &alpha, &a.val_[0], &ione,
            &ione, a.desc_, &val_[0], &ione, &ione, desc_);
#else
        strsm(&side, &uplo, &trans, &diag, &m_, &n_, &alpha, &a.val_[0], &a.m_,
            &val_[0], &m_);
#endif
    }
}
////////////////////////////////////////////////////////////////////////////////
// Solves a triangular system of the form A * X = B or
// A**T * X = B, where A is a triangular matrix of  order  N,
// and  B  is an N-by-NRHS matrix.
// Output in B.
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::trtrs(const char uplo, const char trans,
    const char diag, DistMatrix<double>& b) const
{
    int info;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        int ione = 1;
        pdtrtrs(&uplo, &trans, &diag, &m_, &b.n_, &val_[0], &ione, &ione, desc_,
            &b.val_[0], &ione, &ione, b.desc_, &info);
#else
        dtrtrs(&uplo, &trans, &diag, &m_, &b.n_, &val_[0], &m_, &b.val_[0],
            &b.m_, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::trtrs, info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
    }
}

TEMP_DECL
void DistMatrix<float>::trtrs(const char uplo, const char trans,
    const char diag, DistMatrix<float>& b) const
{
    int info;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        int ione = 1;
        pstrtrs(&uplo, &trans, &diag, &m_, &b.n_, &val_[0], &ione, &ione, desc_,
            &b.val_[0], &ione, &ione, b.desc_, &info);
#else
        strtrs(&uplo, &trans, &diag, &m_, &b.n_, &val_[0], &m_, &b.val_[0],
            &b.m_, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::trtrs, info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
// Computes the Cholesky factorization of a square real
// symmetric positive definite distributed matrix
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
int DistMatrix<double>::potrf(char uplo)
{
    potrf_tm_.start();

    int info = 0;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        int ione = 1;
        pdpotrf(&uplo, &m_, &val_[0], &ione, &ione, desc_, &info);
#else
        dpotrf(&uplo, &m_, &val_[0], &m_, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::potrf, object " << object_name_
                 << ", info=" << info << endl;
        }
    }

    potrf_tm_.stop();

    return info;
}

TEMP_DECL
int DistMatrix<float>::potrf(char uplo)
{
    potrf_tm_.start();

    int info = 0;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        int ione = 1;
        pspotrf(&uplo, &m_, &val_[0], &ione, &ione, desc_, &info);
#else
        spotrf(&uplo, &m_, &val_[0], &m_, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::potrf, object " << object_name_
                 << ", info=" << info << endl;
        }
    }

    potrf_tm_.stop();

    return info;
}
////////////////////////////////////////////////////////////////////////////////
// Computes an LU factorization of a general m by n real
// distributed matrix using partial pivoting with row interchanges
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::getrf(vector<int>& ipiv)
{
    int info;
    if (active_)
    {
#ifdef SCALAPACK
        int ione = 1;
        ipiv.resize(mb_ + mloc_);
        pdgetrf(&m_, &n_, &val_[0], &ione, &ione, desc_, &ipiv[0], &info);
#else
        ipiv.resize(n_);
        dgetrf(&m_, &n_, &val_[0], &m_, &ipiv[0], &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::getrf, info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
    }
}

TEMP_DECL
void DistMatrix<float>::getrf(vector<int>& ipiv)
{
    int info;
    if (active_)
    {
#ifdef SCALAPACK
        int ione = 1;
        ipiv.resize(mb_ + mloc_);
        psgetrf(&m_, &n_, &val_[0], &ione, &ione, desc_, &ipiv[0], &info);
#else
        ipiv.resize(n_);
        sgetrf(&m_, &n_, &val_[0], &m_, &ipiv[0], &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::getrf, info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
// Solve a system of linear equations A*X = B with a symmetric posi-
//  tive definite matrix A using the Cholesky factorization
// A = U**T*U or A = L*L**T computed by potrf
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::potrs(char uplo, DistMatrix<double>& b)
{
    int info;
    if (active_)
    {
        assert(m_ == n_);
        assert(b.m_ == m_);
        assert(b.n_ <= n_);
#ifdef SCALAPACK
        int ione = 1;
        pdpotrs(&uplo, &m_, &b.n_, &val_[0], &ione, &ione, desc_, &b.val_[0],
            &ione, &ione, b.desc_, &info);
#else
        dpotrs(&uplo, &m_, &b.n_, &val_[0], &m_, &b.val_[0], &b.m_, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::potrs, object " << object_name_
                 << ", info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
    }
}

TEMP_DECL
void DistMatrix<float>::potrs(char uplo, DistMatrix<float>& b)
{
    int info;
    if (active_)
    {
        assert(m_ == n_);
        assert(b.m_ == m_);
        assert(b.n_ <= n_);
#ifdef SCALAPACK
        int ione = 1;
        pspotrs(&uplo, &m_, &b.n_, &val_[0], &ione, &ione, desc_, &b.val_[0],
            &ione, &ione, b.desc_, &info);
#else
        spotrs(&uplo, &m_, &b.n_, &val_[0], &m_, &b.val_[0], &b.m_, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::potrs, object " << object_name_
                 << ", info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
// Solve a system of linear equations A * X = B or A' * X = B with a general
// N-by-N matrix A using the LU factorization computed by getrf()
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::getrs(
    char trans, DistMatrix<double>& b, vector<int>& ipiv)
{
    int info;
    if (active_)
    {
        assert(m_ == n_);
        assert(b.m_ == m_);
        assert(b.n_ <= n_);
#ifdef SCALAPACK
        int ione = 1;
        assert((int)ipiv.size() == (mb_ + mloc_));
        pdgetrs(&trans, &m_, &b.n_, &val_[0], &ione, &ione, desc_, &ipiv[0],
            &b.val_[0], &ione, &ione, b.desc_, &info);
#else
        assert((int)ipiv.size() == m_);
        dgetrs(&trans, &m_, &b.n_, &val_[0], &m_, &ipiv[0], &b.val_[0], &b.m_,
            &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::getrs, object " << object_name_
                 << ", info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
    }
}

TEMP_DECL
void DistMatrix<float>::getrs(
    char trans, DistMatrix<float>& b, vector<int>& ipiv)
{
    int info;
    if (active_)
    {
        assert(m_ == n_);
        assert(b.m_ == m_);
        assert(b.n_ <= n_);
#ifdef SCALAPACK
        int ione = 1;
        assert((int)ipiv.size() == (mb_ + mloc_));
        psgetrs(&trans, &m_, &b.n_, &val_[0], &ione, &ione, desc_, &ipiv[0],
            &b.val_[0], &ione, &ione, b.desc_, &info);
#else
        assert((int)ipiv.size() == m_);
        sgetrs(&trans, &m_, &b.n_, &val_[0], &m_, &ipiv[0], &b.val_[0], &b.m_,
            &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::getrs, object " << object_name_
                 << ", info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
// Compute  the inverse of a real symmetric positive definite matrix
// using the Cholesky factorization A = U**T*U or A = L*L**T computed
// by DistMatrix<double>::potrf
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
int DistMatrix<double>::potri(char uplo)
{
    potri_tm_.start();

    int info = 0;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        int ione = 1;
        pdpotri(&uplo, &m_, &val_[0], &ione, &ione, desc_, &info);
#else
        dpotri(&uplo, &m_, &val_[0], &m_, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::potri, object " << object_name_
                 << ", info=" << info << endl;
        }
    }
    potri_tm_.stop();

    return info;
}

TEMP_DECL
int DistMatrix<float>::potri(char uplo)
{
    potri_tm_.start();

    int info = 0;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        int ione = 1;
        pspotri(&uplo, &m_, &val_[0], &ione, &ione, desc_, &info);
#else
        spotri(&uplo, &m_, &val_[0], &m_, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::potri, object " << object_name_
                 << ", info=" << info << endl;
        }
    }
    potri_tm_.stop();

    return info;
}

// DTRTRI - computes the inverse of a real upper or lower triangular matrix A
TEMP_DECL
int DistMatrix<double>::trtri(char uplo, char diag)
{
    trtri_tm_.start();

    int info = 0;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        int ione = 1;
        pdtrtri(&uplo, &diag, &m_, &val_[0], &ione, &ione, desc_, &info);
#else
        dtrtri(&uplo, &diag, &m_, &val_[0], &m_, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::trtri, object " << object_name_
                 << ", info=" << info << endl;
        }
    }
    trtri_tm_.stop();

    return info;
}

////////////////////////////////////////////////////////////////////////////////
// returns the value of the 1-norm, or the Frobenius norm,
// or the infinity norm, or the element of largest absolute value of a
// distributed matrix
//
// max(abs(A(i,j))) if ty== 'm' or 'M'
// 1-norm           if ty == '1'
// infinity-norm    if ty=='I' or 'i'
// Frobenius norm   if ty=='F' or 'f'
TEMP_DECL
double DistMatrix<double>::norm(char ty)
{
    int lwork       = 0;
    double norm_val = 1.;
    if (active_)
    {
#ifdef SCALAPACK
        int ione = 1;
        if (ty == 'I' || ty == 'i') lwork = mloc_;
        if (ty == '1') lwork = nloc_;

        vector<double> work(lwork);

        norm_val
            = pdlange(&ty, &m_, &n_, &val_[0], &ione, &ione, desc_, &work[0]);
#else
        if (ty == 'I' || ty == 'i') lwork = m_;
        vector<double> work(lwork);
        norm_val = dlange(&ty, &m_, &n_, &val_[0], &m_, &work[0]);
#endif
    }
#ifdef USE_MPI
    MPI_Bcast(&norm_val, 1, MPI_DOUBLE, 0, comm_global_);
#endif
    return norm_val;
}

TEMP_DECL
double DistMatrix<float>::norm(char ty)
{
    int lwork      = 0;
    float norm_val = 1.;
    if (active_)
    {
#ifdef SCALAPACK
        int ione = 1;
        if (ty == 'I' || ty == 'i') lwork = mloc_;
        if (ty == '1') lwork = nloc_;

        vector<float> work(lwork);

        norm_val = (double)pslange(
            &ty, &m_, &n_, &val_[0], &ione, &ione, desc_, &work[0]);
#else
        if (ty == 'I' || ty == 'i') lwork = m_;
        vector<float> work(lwork);
        norm_val     = (double)slange(&ty, &m_, &n_, &val_[0], &m_, &work[0]);
#endif
    }
#ifdef USE_MPI
    MPI_Bcast(&norm_val, 1, MPI_FLOAT, 0, comm_global_);
#endif
    return norm_val;
}

////////////////////////////////////////////////////////////////////////////////
// estimate the reciprocal of the condition number (in the 1-norm) of a
// real symmetric positive definite matrix given by its Cholesky factorization
// U or L such that A  = U**T*U or A = L*L**T
// (as computed by DistMatrix<double>::potrf)
// anorm: 1-norm of matrix A
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
double DistMatrix<double>::pocon(char uplo, double anorm)
{
    int info;
    double rcond = 1.;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        int ione     = 1;
        int lwork    = 2 * mloc_ + 3 * nloc_ + nb_ * min(npcol_, nprow_);
        lwork        = max(lwork, 1);
        int liwork   = max(mloc_, 1);
        double* work = new double[lwork];
        int* iwork   = new int[liwork];
        pdpocon(&uplo, &m_, &val_[0], &ione, &ione, desc_, &anorm, &rcond, work,
            &lwork, iwork, &liwork, &info);
        if (info != 0)
        {
            cerr << "mloc_=" << mloc_ << ", nloc_=" << nloc_ << endl;
            cerr << "npcol_=" << npcol_ << ", nprow_=" << nprow_ << endl;
            cerr << "lwork=" << lwork << ", but should be at least " << work[0]
                 << endl;
            cerr << "liwork=" << liwork << ", but should be at least "
                 << iwork[0] << endl;
        }
        delete[] iwork;
        delete[] work;
#else
        double* work = new double[3 * m_];
        int* iwork   = new int[m_];
        dpocon(&uplo, &m_, &val_[0], &m_, &anorm, &rcond, work, iwork, &info);
        delete[] iwork;
        delete[] work;
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::pocon, info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
    }
    return rcond;
}

TEMP_DECL
double DistMatrix<float>::pocon(char uplo, float anorm)
{
    int info;
    float rcond = 1.;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        int ione    = 1;
        int lwork   = 2 * mloc_ + 3 * nloc_ + nb_ * min(npcol_, nprow_);
        lwork       = max(lwork, 1);
        int liwork  = max(mloc_, 1);
        float* work = new float[lwork];
        int* iwork  = new int[liwork];
        pspocon(&uplo, &m_, &val_[0], &ione, &ione, desc_, &anorm, &rcond, work,
            &lwork, iwork, &liwork, &info);
        if (info != 0)
        {
            cerr << "mloc_=" << mloc_ << ", nloc_=" << nloc_ << endl;
            cerr << "npcol_=" << npcol_ << ", nprow_=" << nprow_ << endl;
            cerr << "lwork=" << lwork << ", but should be at least " << work[0]
                 << endl;
            cerr << "liwork=" << liwork << ", but should be at least "
                 << iwork[0] << endl;
        }
        delete[] iwork;
        delete[] work;
#else
        float* work = new float[3 * m_];
        int* iwork  = new int[m_];
        spocon(&uplo, &m_, &val_[0], &m_, &anorm, &rcond, work, iwork, &info);
        delete[] iwork;
        delete[] work;
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::pocon, info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
    }
    return (double)rcond;
}
////////////////////////////////////////////////////////////////////////////////
// symmetric rank k update
// this = beta * this + alpha * A * A^T  (trans=='n')
// this = beta * this + alpha * A^T * A  (trans=='t')
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::syrk(
    char uplo, char trans, double alpha, DistMatrix<double>& a, double beta)
{
    assert(ictxt_ == a.ictxt());
    assert(n_ == m_); // *this must be a square matrix

    if (active_)
    {
        int n, k;
        if (trans == 'N' || trans == 'n')
        {
            n = m_;
            k = a.n();
        }
        else
        {
            n = m_;
            k = a.m();
        }

#ifdef SCALAPACK
        int ione = 1;
        pdsyrk(&uplo, &trans, &n, &k, &alpha, &a.val_[0], &ione, &ione,
            a.desc(), &beta, &val_[0], &ione, &ione, desc_);
#else
        dsyrk(&uplo, &trans, &n, &k, &alpha, &a.val_[0], &a.m_, &beta, &val_[0],
            &m_);
#endif
    }
}

TEMP_DECL
void DistMatrix<float>::syrk(
    char uplo, char trans, float alpha, DistMatrix<float>& a, float beta)
{
    assert(ictxt_ == a.ictxt());
    assert(n_ == m_); // *this must be a square matrix

    if (active_)
    {
        int n, k;
        if (trans == 'N' || trans == 'n')
        {
            n = m_;
            k = a.n();
        }
        else
        {
            n = m_;
            k = a.m();
        }

#ifdef SCALAPACK
        int ione = 1;
        pssyrk(&uplo, &trans, &n, &k, &alpha, &a.val_[0], &ione, &ione,
            a.desc(), &beta, &val_[0], &ione, &ione, desc_);
#else
        ssyrk(&uplo, &trans, &n, &k, &alpha, &a.val_[0], &a.m_, &beta, &val_[0],
            &m_);
#endif
    }
}
////////////////////////////////////////////////////////////////////////////////
// getsub: *this = sub(A)
// copy submatrix A(ia:ia+m, ja:ja+n) into *this;
// *this and A may live in different contexts,
// Assumes that the context of *this is a subcontext of the context of A
TEMP_DECL
void DistMatrix<double>::getsub(
    const DistMatrix<double>& a, int m, int n, int ia, int ja)
{
#ifdef SCALAPACK
    int iap = ia + 1;
    int jap = ja + 1;
    assert(n <= n_);
    assert(n <= a.n());
    assert(m <= m_);
    assert(m <= a.m());
    int ione = 1;
    BlacsContext bc(comm_global_);
    int gictxt = bc.ictxt();
    // int gictxt = ictxt_;
    pdgemr2d(&m, &n, a.val_, &iap, &jap, a.desc_, val_, &ione, &ione, desc_,
        &gictxt);
#else
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
            val_[i + j * m_] = a.val_[(i + ia) + (j + ja) * a.m()];
#endif
}

TEMP_DECL
void DistMatrix<float>::getsub(
    const DistMatrix<float>& a, int m, int n, int ia, int ja)
{
#ifdef SCALAPACK
    int iap = ia + 1;
    int jap = ja + 1;
    assert(n <= n_);
    assert(n <= a.n());
    assert(m <= m_);
    assert(m <= a.m());
    int ione = 1;
    BlacsContext bc(comm_global_);
    int gictxt = bc.ictxt();
    // int gictxt = ictxt_;
    psgemr2d(&m, &n, a.val_, &iap, &jap, a.desc_, val_, &ione, &ione, desc_,
        &gictxt);
#else
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
            val_[i + j * m_] = a.val_[(i + ia) + (j + ja) * a.m()];
#endif
}
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<int>::getsub(
    const DistMatrix<int>& a, int m, int n, int ia, int ja)
{
#ifdef SCALAPACK
    int iap = ia + 1;
    int jap = ja + 1;
    assert(n <= n_);
    assert(n <= a.n());
    assert(m <= m_);
    assert(m <= a.m());
    int ione = 1;
    BlacsContext bc(comm_global_);
    int gictxt = bc.ictxt();
    pigemr2d(&m, &n, a.val_, &iap, &jap, a.desc(), val_, &ione, &ione, desc_,
        &gictxt);
#else
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
            val_[i + j * m_] = a.val_[(i + ia) + (j + ja) * a.m()];
#endif
}

////////////////////////////////////////////////////////////////////////////////
// matrix transpose
// this = alpha * transpose(A) + beta * this
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::transpose(
    const double alpha, const DistMatrix<double>& a, const double beta)
{
    assert(ictxt_ == a.ictxt());

    if (active_)
    {
        assert(a.m() == n_);
        assert(a.n() == m_);

#ifdef SCALAPACK
        int ione = 1;
        pdtran(&m_, &n_, &alpha, &a.val_[0], &ione, &ione, a.desc(), &beta,
            &val_[0], &ione, &ione, desc_);
#else
        scal(beta);
        for (int i = 0; i < m_; i++)
            for (int j = 0; j < i; j++)
            {
                val_[i * m_ + j] += alpha * a.val_[j * m_ + i];
                val_[j * m_ + i] += alpha * a.val_[i * m_ + j];
            }
        for (int i = 0; i < m_; i++)
            val_[i * m_ + i] += alpha * a.val_[i * m_ + i];
#endif
    }
}

TEMP_DECL
void DistMatrix<float>::transpose(
    const float alpha, const DistMatrix<float>& a, const float beta)
{
    assert(ictxt_ == a.ictxt());

    if (active_)
    {
        assert(a.m() == n_);
        assert(a.n() == m_);

#ifdef SCALAPACK
        int ione = 1;
        pstran(&m_, &n_, &alpha, &a.val_[0], &ione, &ione, a.desc(), &beta,
            &val_[0], &ione, &ione, desc_);
#else
        scal(beta);
        for (int i = 0; i < m_; i++)
            for (int j = 0; j < i; j++)
            {
                val_[i * m_ + j] += alpha * a.val_[j * m_ + i];
                val_[j * m_ + i] += alpha * a.val_[i * m_ + j];
            }
        for (int i = 0; i < m_; i++)
            val_[i * m_ + i] += alpha * a.val_[i * m_ + i];
#endif
    }
}
////////////////////////////////////////////////////////////////////////////////
// matrix transpose
// *this = transpose(a)
////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::transpose(const DistMatrix<T>& a)
{
    assert(this != &a);
    transpose(1.0, a, 0.0);
}

////////////////////////////////////////////////////////////////////////////////
// set to zero upper or lower part of a lower or upper triangular matrix
////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::trset(const char uplo)
{
    if (active_)
    {

        assert(m_ == n_);

        if (uplo == 'l' || uplo == 'L')
        {
            for (int li = 0; li < mblocks_; li++)
            {
                for (int lj = 0; lj < nblocks_; lj++)
                {
                    for (int ii = 0; ii < mbs(li); ii++)
                    {
                        for (int jj = 0; jj < nbs(lj); jj++)
                        {
                            if (i(li, ii) < j(lj, jj))
                                val_[(ii + li * mb_) + (jj + lj * nb_) * mloc_]
                                    = 0.0;
                        }
                    }
                }
            }
        }
        else
        {
            for (int li = 0; li < mblocks_; li++)
            {
                for (int lj = 0; lj < nblocks_; lj++)
                {
                    for (int ii = 0; ii < mbs(li); ii++)
                    {
                        for (int jj = 0; jj < nbs(lj); jj++)
                        {
                            if (i(li, ii) > j(lj, jj))
                                val_[(ii + li * mb_) + (jj + lj * nb_) * mloc_]
                                    = 0.0;
                        }
                    }
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::init(const T* const a, const int lda)
{
    if (active_)
    {
        for (int li = 0; li < mblocks_; li++)
        {
            for (int lj = 0; lj < nblocks_; lj++)
            {
                for (int ii = 0; ii < mbs(li); ii++)
                {
                    for (int jj = 0; jj < nbs(lj); jj++)
                    {
                        val_[(ii + li * mb_) + (jj + lj * nb_) * mloc_]
                            = a[i(li, ii) + j(lj, jj) * lda];
                    }
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::initTest()
{
    if (active_)
    {
        for (int li = 0; li < mblocks_; li++)
        {
            for (int lj = 0; lj < nblocks_; lj++)
            {
                for (int ii = 0; ii < mbs(li); ii++)
                {
                    for (int jj = 0; jj < nbs(lj); jj++)
                    {
                        val_[(ii + li * mb_) + (jj + lj * nb_) * mloc_]
                            = 1000 * i(li, ii) + j(lj, jj);
                    }
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
double DistMatrix<T>::sumProdElements(const DistMatrix<T>& a) const
{
    assert(mblocks_ == a.mblocks_);
    assert(nblocks_ == a.nblocks_);
    assert(mb_ == a.mb_);
    assert(nb_ == a.nb_);

    double tsum = 0.;

    if (active_)
    {
        for (int li = 0; li < mblocks_; li++)
        {
            for (int lj = 0; lj < nblocks_; lj++)
            {
                for (int ii = 0; ii < mbs(li); ii++)
                {
                    for (int jj = 0; jj < nbs(lj); jj++)
                    {
                        tsum += (double)val_[(ii + li * mb_)
                                             + (jj + lj * nb_) * mloc_]
                                * (double)a.val_[(ii + li * mb_)
                                                 + (jj + lj * nb_) * mloc_];
                    }
                }
            }
        }
    }

    double sum = 0.;
#ifdef SCALAPACK
    MPI_Allreduce(&tsum, &sum, 1, MPI_DOUBLE, MPI_SUM, comm_global_);
#else
    sum = tsum;
#endif
    return sum;
}

////////////////////////////////////////////////////////////////////////////////
// Constructor for a distributed matrix
// from a duplicated matrix
////////////////////////////////////////////////////////////////////////////////
template <class T>
DistMatrix<T>::DistMatrix(const string& name, const BlacsContext& bc,
    const int lda, const T* const a, const int m, const int n, const int mb,
    const int nb)
    : object_name_(name), comm_global_(bc.comm_global()), bc_(bc)
{
    resize(m, n, mb, nb);
    init(a, lda);
}

template <class T>
DistMatrix<T>::DistMatrix(const string& name, const BlacsContext& bc,
    const int lda, const T* const a, const int m, const int n)
    : object_name_(name), comm_global_(bc.comm_global()), bc_(bc)
{
    resize(m, n, distmatrix_def_block_size_, distmatrix_def_block_size_);
    init(a, lda);
}

////////////////////////////////////////////////////////////////////////////////
//
// Generate a duplicated matrix from a distributed matrix
//
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::matgather(double* a, const int lda) const
{
    matgather_tm_.start();

    assert(m_ == n_);

#ifdef SCALAPACK
    const int ii  = -1;
    const int jj  = -1;
    const int rev = 0;
    const int i   = 1;
    if (active_) pdlacp3(&m_, &i, &val_[0], desc_, a, &lda, &ii, &jj, &rev);

    // int nrpocs=bc_.size();
    int nrpocs     = bc_.nprocs();
    int n_active   = bc_.nprow() * bc_.npcol();
    int n_inactive = nrpocs - n_active;
    int ns         = n_inactive / n_active;
    if (ns * n_active < n_inactive) ns++;
    int mpirc = MPI_SUCCESS;
    int size  = lda * n_;
    int mype  = bc_.mype();
    MPI_Status status;
    // send data from active PEs to inactive PEs
    for (int i = 0; i < ns; i++)
    {
        if (!active_ && mype >= (i + 1) * n_active && mype < (i + 2) * n_active)
        {
            int rpe = mype % n_active; // remote PE
            mpirc
                = MPI_Recv(a, size, MPI_DOUBLE, rpe, 0, comm_global_, &status);
        }
        if (active_ && (mype + (i + 1) * n_active) < nrpocs)
            mpirc = MPI_Send(a, size, MPI_DOUBLE, mype + (i + 1) * n_active, 0,
                comm_global_);
        if (mpirc != MPI_SUCCESS)
            cerr << "DistMatrix<T>::matgather: MPI Send/Recv failed for PE "
                 << mype << endl;
    }
#else
    memcpy(a, &val_[0], lda * n_ * sizeof(double));

#endif
    matgather_tm_.stop();
}

TEMP_DECL
void DistMatrix<float>::matgather(float* a, const int lda) const
{
    matgather_tm_.start();

    assert(m_ == n_);

#ifdef SCALAPACK
    const int ii  = -1;
    const int jj  = -1;
    const int rev = 0;
    const int i   = 1;
    if (active_) pslacp3(&m_, &i, &val_[0], desc_, a, &lda, &ii, &jj, &rev);

    // int nrpocs=bc_.size();
    int nrpocs     = bc_.nprocs();
    int n_active   = bc_.nprow() * bc_.npcol();
    int n_inactive = nrpocs - n_active;
    int ns         = n_inactive / n_active;
    if (ns * n_active < n_inactive) ns++;
    int mpirc = MPI_SUCCESS;
    int size  = lda * n_;
    int mype  = bc_.mype();
    MPI_Status status;
    // send data from active PEs to inactive PEs
    for (int i = 0; i < ns; i++)
    {
        if (!active_ && mype >= (i + 1) * n_active && mype < (i + 2) * n_active)
        {
            int rpe = mype % n_active; // remote PE
            mpirc = MPI_Recv(a, size, MPI_FLOAT, rpe, 0, comm_global_, &status);
        }
        if (active_ && (mype + (i + 1) * n_active) < nrpocs)
            mpirc = MPI_Send(
                a, size, MPI_FLOAT, mype + (i + 1) * n_active, 0, comm_global_);
        if (mpirc != MPI_SUCCESS)
            cerr << "DistMatrix<T>::matgather: MPI Send/Recv failed for PE "
                 << mype << endl;
    }
#else
    memcpy(a, &val_[0], lda * n_ * sizeof(float));

#endif
    matgather_tm_.stop();
}
////////////////////////////////////////////////////////////////////////////////
// Constructor for a diagonal distributed matrix
// from the diagonal elements
////////////////////////////////////////////////////////////////////////////////
template <class T>
DistMatrix<T>::DistMatrix(const string& name, const BlacsContext& bc,
    const T* const dmat, const int m, const int n, const int mb, const int nb)
    : object_name_(name), comm_global_(bc.comm_global()), bc_(bc)
{
    assert(m == n);
    assert(mb == nb);

    resize(m, n, mb, nb);
    setDiagonalValues(dmat);
}

template <class T>
DistMatrix<T>::DistMatrix(const string& name, const BlacsContext& bc,
    const T* const dmat, const int m, const int n)
    : object_name_(name), comm_global_(bc.comm_global()), bc_(bc)
{
    resize(m, n, distmatrix_def_block_size_, distmatrix_def_block_size_);
    setDiagonalValues(dmat);
}

template <class T>
DistMatrix<T>::DistMatrix(
    const string& name, const T* const dmat, const int m, const int n)
    : object_name_(name),
      comm_global_(default_bc_->comm_global()),
      bc_(*default_bc_)
{
    resize(m, n, distmatrix_def_block_size_, distmatrix_def_block_size_);
    setDiagonalValues(dmat);
}

template <class T>
void DistMatrix<T>::setDiagonalValues(const T* const dmat)
{
    if (active_)
    {
        if (myrow_ == mycol_)
        {
            for (int li = 0; li < nblocks_; li++)
            {
                for (int ii = 0; ii < nbs(li); ii++)
                {
                    val_[(ii + li * nb_) * (mloc_ + 1)] = dmat[i(li, ii)];
                }
            }
        }
    }
}

template <class T>
void DistMatrix<T>::getDiagonalValues(T* const dmat) const
{
    memset(dmat, 0, m_ * sizeof(T));

    if (active_)
    {
        if (myrow_ == mycol_)
        {
            for (int li = 0; li < nblocks_; li++)
            {
                for (int ii = 0; ii < nbs(li); ii++)
                {
                    dmat[i(li, ii)] = val_[(ii + li * nb_) * (mloc_ + 1)];
                }
            }
        }
    }
#ifdef USE_MPI
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    T* recvbuf      = new T[m_];
    mmpi.allreduce(dmat, recvbuf, m_, MPI_SUM);
    //  MPI_Allreduce ( dmat, recvbuf, m_,
    //                   MPI_DOUBLE, MPI_SUM, comm_global_ );
    memcpy(dmat, recvbuf, m_ * sizeof(T));
    delete[] recvbuf;
#endif
}

////////////////////////////////////////////////////////////////////////////////
// Compute the trace of the DistMatrix
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
double DistMatrix<double>::trace(void) const
{
    assert(m_ == n_);
    double trace = 0.0;

    if (active_)
    {
#ifdef SCALAPACK
        int ione = 1;
        trace    = pdlatra(&n_, &val_[0], &ione, &ione, desc_);
#else
        int inc = m_ + 1;
        for (int i = 0; i < n_; i++)
            trace += val_[i * inc];
#endif
    }

#ifdef USE_MPI
    MPI_Bcast(&trace, 1, MPI_DOUBLE, 0, comm_global_);
#endif
    return trace;
}

TEMP_DECL
double DistMatrix<float>::trace(void) const
{
    assert(m_ == n_);
    double trace = 0.0;

    if (active_)
    {
#ifdef SCALAPACK
        int ione = 1;
        trace    = (double)pslatra(&n_, &val_[0], &ione, &ione, desc_);
#else
        int inc = m_ + 1;
        for (int i = 0; i < n_; i++)
            trace += val_[i * inc];
#endif
    }

#ifdef USE_MPI
    MPI_Bcast(&trace, 1, MPI_DOUBLE, 0, comm_global_);
#endif
    return trace;
}
////////////////////////////////////////////////////////////////////////////////
// Reduces a real  symmetric-definite  generalized  eigenproblem  to  standard
// form.  If itype = 1, the problem is A*x = lambda*B*x,
// and A (=*this) is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
// If itype = 2 or 3, the problem is A*B*x = lambda*x or
// B*A*x = lambda*x, and *this is overwritten by U*A*U**T or L**T*A*L.
// B must have been previously factorized as U**T*U or L*L**T by
// DistMatrix<double>::dpotrf.
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::sygst(
    int itype, char uplo, const DistMatrix<double>& b)
{
    int info;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        int ione = 1;
        double scale;
        pdsygst(&itype, &uplo, &m_, &val_[0], &ione, &ione, desc_, &b.val_[0],
            &ione, &ione, b.desc_, &scale, &info);
#else
        dsygst(&itype, &uplo, &m_, &val_[0], &m_, &b.val_[0], &b.m_, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::sygst, info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
    }
}

TEMP_DECL
void DistMatrix<float>::sygst(int itype, char uplo, const DistMatrix<float>& b)
{
    int info;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        int ione = 1;
        float scale;
        pssygst(&itype, &uplo, &m_, &val_[0], &ione, &ione, desc_, &b.val_[0],
            &ione, &ione, b.desc_, &scale, &info);
#else
        ssygst(&itype, &uplo, &m_, &val_[0], &m_, &b.val_[0], &b.m_, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::sygst, info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::syev(
    char jobz, char uplo, vector<double>& w, DistMatrix<double>& z)
{
    int info;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        assert(mb_ > 1);
        int ione = 1, izero = 0;
        int nn    = max(mb_, m_);
        int np0   = numroc(&nn, &mb_, &izero, &izero, &nprow_);
        int lwork = 5 * m_ + 3 * np0 + mb_ * nn + m_ * np0;
        lwork *= 2;
        double* work = new double[lwork];
        pdsyev(&jobz, &uplo, &m_, &val_[0], &ione, &ione, desc_, &w[0],
            &z.val_[0], &ione, &ione, z.desc_, work, &lwork, &info);
#else
        int lwork    = 3 * m_;
        double* work = new double[lwork];
        z            = *this;
        dsyev(&jobz, &uplo, &m_, &z.val_[0], &m_, &w[0], work, &lwork, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::syev, info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
        delete[] work;
    }

#ifdef SCALAPACK
    MPI_Bcast(&w[0], m_, MPI_DOUBLE, 0, comm_global_);
#endif
}

TEMP_DECL
void DistMatrix<float>::syev(
    char jobz, char uplo, vector<float>& w, DistMatrix<float>& z)
{
    int info;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        assert(mb_ > 1);
        int ione = 1, izero = 0;
        int nn    = max(mb_, m_);
        int np0   = numroc(&nn, &mb_, &izero, &izero, &nprow_);
        int lwork = 5 * m_ + 3 * np0 + mb_ * nn + m_ * np0;
        lwork *= 2;
        float* work = new float[lwork];
        pssyev(&jobz, &uplo, &m_, &val_[0], &ione, &ione, desc_, &w[0],
            &z.val_[0], &ione, &ione, z.desc_, work, &lwork, &info);
#else
        int lwork   = 3 * m_;
        float* work = new float[lwork];
        z           = *this;
        ssyev(&jobz, &uplo, &m_, &z.val_[0], &m_, &w[0], work, &lwork, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::syev, info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
        delete[] work;
    }

#ifdef SCALAPACK
    MPI_Bcast(&w[0], m_, MPI_FLOAT, 0, comm_global_);
#endif
}
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::gesvd(char jobu, char jobvt, vector<double>& s,
    DistMatrix<double>& u, DistMatrix<double>& vt)
{
    int info;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        assert(mb_ > 1);
        int ione  = 1;
        int lwork = -1;
        double qwork[1];
        pdgesvd(&jobu, &jobvt, &m_, &n_, &val_[0], &ione, &ione, desc_, &s[0],
            &u.val_[0], &ione, &ione, u.desc_, &vt.val_[0], &ione, &ione,
            vt.desc_, &qwork[0], &lwork, &info);
        lwork        = qwork[0];
        double* work = new double[lwork];
        pdgesvd(&jobu, &jobvt, &m_, &n_, &val_[0], &ione, &ione, desc_, &s[0],
            &u.val_[0], &ione, &ione, u.desc_, &vt.val_[0], &ione, &ione,
            vt.desc_, work, &lwork, &info);
#else
        int lwork = -1;
        double qwork[1];
        dgesvd(&jobu, &jobvt, &m_, &n_, &val_[0], &m_, &s[0], &u.val_[0], &m_,
            &vt.val_[0], &n_, qwork, &lwork, &info);
        lwork        = qwork[0];
        double* work = new double[lwork];
        dgesvd(&jobu, &jobvt, &m_, &n_, &val_[0], &m_, &s[0], &u.val_[0], &m_,
            &vt.val_[0], &n_, qwork, &lwork, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::gesvd, info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
        delete[] work;
    }

#ifdef SCALAPACK
    int m = min(m_, n_);
    MPI_Bcast(&s[0], m, MPI_DOUBLE, 0, comm_global_);
#endif
}

TEMP_DECL
void DistMatrix<float>::gesvd(char jobu, char jobvt, vector<float>& s,
    DistMatrix<float>& u, DistMatrix<float>& vt)
{
    int info;
    if (active_)
    {
        assert(m_ == n_);

#ifdef SCALAPACK
        assert(mb_ > 1);
        int ione  = 1;
        int lwork = -1;
        float qwork[1];
        psgesvd(&jobu, &jobvt, &m_, &n_, &val_[0], &ione, &ione, desc_, &s[0],
            &u.val_[0], &ione, &ione, u.desc_, &vt.val_[0], &ione, &ione,
            vt.desc_, &qwork[0], &lwork, &info);
        lwork       = qwork[0];
        float* work = new float[lwork];
        psgesvd(&jobu, &jobvt, &m_, &n_, &val_[0], &ione, &ione, desc_, &s[0],
            &u.val_[0], &ione, &ione, u.desc_, &vt.val_[0], &ione, &ione,
            vt.desc_, work, &lwork, &info);
#else
        int lwork = -1;
        float qwork[1];
        sgesvd(&jobu, &jobvt, &m_, &n_, &val_[0], &m_, &s[0], &u.val_[0], &m_,
            &vt.val_[0], &n_, qwork, &lwork, &info);
        lwork       = qwork[0];
        float* work = new float[lwork];
        sgesvd(&jobu, &jobvt, &m_, &n_, &val_[0], &m_, &s[0], &u.val_[0], &m_,
            &vt.val_[0], &n_, qwork, &lwork, &info);
#endif
        if (info != 0)
        {
            cerr << " DistMatrix::gesvd, info=" << info << endl;
#ifdef USE_MPI
            MPI_Abort(comm_global_, 2);
#else
            exit(2);
#endif
        }
        delete[] work;
    }

#ifdef SCALAPACK
    int m = min(m_, n_);
    MPI_Bcast(&s[0], m, MPI_FLOAT, 0, comm_global_);
#endif
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::sygv(char jobz, char uplo, const char fac, DistMatrix<T>& b,
    vector<T>& w, DistMatrix<T>& z)
{
    // compute Cholesky factorization of b if necessary
    if (fac != 'c') b.potrf(uplo);

    // Transform the generalized eigenvalue problem to a standard form
    sygst(1, uplo, b);

    // solve a standard symmetric eigenvalue problem
    syev('v', uplo, w, z);

    // Get the eigenvectors Z of the generalized eigenvalue problem
    // Solve Z=L**(-T)*U
    b.trtrs(uplo, 't', 'n', z);
}

////////////////////////////////////////////////////////////////////////////////
// get max in absolute value of column j
TEMP_DECL
int DistMatrix<double>::iamax(const int j, double& val)
{
    int indx = -1;
    int incx = 1;
#ifdef SCALAPACK
    int ix          = 1;
    int jx          = y(j) + 1; // C to fortran
    int proc_col    = pc(j);
    int desc_col[9] = { desc_[0], desc_[1], desc_[2], 1, desc_[4], 1, desc_[6],
        proc_col, desc_[8] };

    if (mycol_ == proc_col)
        pdamax(&m_, &val, &indx, &val_[y(j) * desc_[8]], &ix, &jx, desc_col,
            &incx);
    indx--; // fortran to C

    int proot = proc_col;
    MPI_Bcast(&indx, 1, MPI_INT, proot, comm_global_);
    MPI_Bcast(&val, 1, MPI_DOUBLE, proot, comm_global_);
#else
    indx     = idamax(&m_, &val_[0] + m_ * j, &incx) - 1;
    val      = val_[m_ * j + indx];
#endif
    return indx;
}

TEMP_DECL
int DistMatrix<float>::iamax(const int j, float& val)
{
    int indx = -1;
    int incx = 1;
#ifdef SCALAPACK
    int ix          = 1;
    int jx          = y(j) + 1; // C to fortran
    int proc_col    = pc(j);
    int desc_col[9] = { desc_[0], desc_[1], desc_[2], 1, desc_[4], 1, desc_[6],
        proc_col, desc_[8] };

    if (mycol_ == proc_col)
        psamax(&m_, &val, &indx, &val_[y(j) * desc_[8]], &ix, &jx, desc_col,
            &incx);
    indx--; // fortran to C

    int proot = proc_col;
    MPI_Bcast(&indx, 1, MPI_INT, proot, comm_global_);
    MPI_Bcast(&val, 1, MPI_FLOAT, proot, comm_global_);
#else
    indx     = isamax(&m_, &val_[0] + m_ * j, &incx) - 1;
    val      = val_[m_ * j + indx];
#endif
    return indx;
}
////////////////////////////////////////////////////////////////////////////////
TEMP_DECL
void DistMatrix<double>::swapColumns(const int j1, const int j2)
{
    if (j1 == j2) return;
#ifdef SCALAPACK
    int jx   = j1 + 1; // C to Fortran
    int jy   = j2 + 1; // C to Fortran
    int ione = 1;
    int inc  = 1; // indicates swap columns
    if (active_)
        pdswap(&m_, &val_[0], &ione, &jx, desc_, &inc, &val_[0], &ione, &jy,
            desc_, &inc);
#else
    int ione = 1;
    dswap(&m_, &val_[0] + m_ * j1, &ione, &val_[0] + m_ * j2, &ione);
#endif
}

TEMP_DECL
void DistMatrix<float>::swapColumns(const int j1, const int j2)
{
    if (j1 == j2) return;
#ifdef SCALAPACK
    int jx   = j1 + 1; // C to Fortran
    int jy   = j2 + 1; // C to Fortran
    int ione = 1;
    int inc  = 1; // indicates swap columns
    if (active_)
        psswap(&m_, &val_[0], &ione, &jx, desc_, &inc, &val_[0], &ione, &jy,
            desc_, &inc);
#else
    int ione = 1;
    sswap(&m_, &val_[0] + m_ * j1, &ione, &val_[0] + m_ * j2, &ione);
#endif
}
////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::print(ostream& os) const
{
    if (m_ == 0 || n_ == 0) return;
    BlacsContext bcl(comm_global_, 1, 1);
#if 0
  // First version: uses a complete copy of *this on process (0,0)
  DistMatrix<T> t(bcl,m_,n_,mb_,nb_);
  t = *this;
  if ( t.active() )
  {
    // this is done only on pe 0
    os << " print: (m,n)=" << "(" << m_ << "," << n_ << ")" << endl;
    for ( int ii = 0; ii < m_; ii++ )
    {
      for ( int jj = 0; jj < n_; jj++ )
      {
        os << "(" << ii << "," << jj << ")=" << t.val_[ii+t.mloc()*jj] << endl;
      }
    }
  }
#else
    // Copy blocks of <blocksize> columns and print them on process (0,0)
    const int blockmemsize = 32768; // maximum memory size of a block in bytes
    // compute maximum block size: must be at least 1
    int maxbs = max(1, (int)((blockmemsize / sizeof(T)) / m_));
    DistMatrix<T> t("DistMatrix<T>::print", bcl, m_, maxbs);
    int nblocks = n_ / maxbs + ((n_ % maxbs == 0) ? 0 : 1);
    int ia      = 0;
    int ja      = 0;
    for (int jb = 0; jb < nblocks; jb++)
    {
        int blocksize = ((jb + 1) * maxbs > n_) ? n_ % maxbs : maxbs;
        t.getsub(*this, t.m(), blocksize, ia, ja);
        ja += blocksize;
        if (t.active())
        {
            os << scientific;
            // this is done only on pe 0
            for (int jj = 0; jj < blocksize; jj++)
            {
                for (int ii = 0; ii < m_; ii++)
                {
                    os << "(" << ii << "," << jj + jb * maxbs
                       << ")=" << t.val_[ii + t.mloc() * jj] << endl;
                }
            }
        }
    }
#endif
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::printMM(ostream& os) const
{
    if (m_ == 0 || n_ == 0) return;
    BlacsContext bcl(comm_global_, 1, 1);
    // Copy blocks of <blocksize> columns and print them on process (0,0)
    const int blockmemsize = 32768; // maximum memory size of a block in bytes
    // compute maximum block size: must be at least 1
    const int maxbs = max(1, (int)((blockmemsize / sizeof(T)) / m_));
    DistMatrix<T> t("DistMatrix<T>::print", bcl, m_, maxbs);
    const int nblocks = n_ / maxbs + ((n_ % maxbs == 0) ? 0 : 1);
    int ia            = 0;
    int ja            = 0;
    if (t.active())
    {
        os << "%%MatrixMarket matrix array real general" << endl;
        os << m_ << " " << n_ << endl;
    }
    for (int jb = 0; jb < nblocks; jb++)
    {
        const int blocksize = ((jb + 1) * maxbs > n_) ? n_ % maxbs : maxbs;
        t.getsub(*this, t.m(), blocksize, ia, ja);
        ja += blocksize;
        if (t.active()) // this is done only on pe 0
        {
            for (int jj = 0; jj < blocksize; jj++)
            {
                // print one column
                for (int ii = 0; ii < m_; ii++)
                {
                    os << t.val_[ii + t.mloc() * jj] << endl;
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::axpyColumn(
    const int icol, const double alpha, const DistMatrix<T>& x)
{
    const int ni = mloc_ * icol;
    //    int ione=1;
    MPaxpy(mloc_, alpha, &x.val_[ni], &val_[ni]);
    //    daxpy(&mloc_, &alpha, &x.val_[ni], &ione,
    //                          &val_[ni],   &ione);
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::scalColumn(const int icol, const double alpha)
{
    const int ni = mloc_ * icol;
    //    int ione=1;
    MPscal(mloc_, alpha, &val_[ni]);
    //    dscal(&mloc_, &alpha, &val_[ni],   &ione);
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
void DistMatrix<T>::print(
    ostream& os, const int ia, const int ja, const int ma, const int na) const
{
    const int m = min(ma, max(m_ - ia, 0));
    const int n = min(na, max(n_ - ja, 0));
    BlacsContext bcl(comm_global_, 1, 1);
    DistMatrix<T> t("DistMatrix<T>::print", bcl, m, n);

#ifdef USE_MPI
    MPI_Barrier(comm_global_);
#endif
    t.getsub(*this, m, n, ia, ja);
    if (t.active())
    {
        os << setprecision(5);
        for (int ii = 0; ii < m; ii++)
        {
            // this is done only on pe 0
            for (int jj = 0; jj < n; jj++)
            {
                os << t.val_[ii + m * jj] << "\t";
            }
            os << endl;
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
// get value for global index (i,j)
// assuming we are on the right processor to get it!
template <class T>
const T DistMatrix<T>::getGlobalVal(
    const int i, const int j, const bool onpe) const
{
    assert(i < m_);
    assert(j < n_);
    assert(i >= 0);
    assert(j >= 0);

    bool flag = true;
    T val     = 0.;
    if (!onpe)
    {
        const int taski = (i / mb_) % nprow_;
        const int taskj = (j / nb_) % npcol_;
        flag            = (taski == myrow_ && taskj == mycol_);
    }

    if (flag)
    {
        const int ib = i / (nprow_ * mb_);
        const int jb = j / (npcol_ * nb_);
        const int x  = i % mb_;
        const int y  = j % nb_;
        assert(ib * mb_ + x + (jb * nb_ + y) * mloc_ < size_);
        val = val_[ib * mb_ + x + (jb * nb_ + y) * mloc_];
    }

#ifdef USE_MPI
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    if (!onpe)
        mmpi.bcast(
            &val, 1, 0); // MPI_Bcast(&val, 1, MPI_DOUBLE, 0, comm_global_);
#endif
    return val;
}

////////////////////////////////////////////////////////////////////////////////
// set value for global index (i,j)
template <class T>
void DistMatrix<T>::setGlobalVal(
    const int i, const int j, const T val, const bool onpe)
{
    assert(i < m_);
    assert(j < n_);
    assert(i >= 0);
    assert(j >= 0);

    bool flag = true;
    if (!onpe)
    {
        const int taski = (i / mb_) % nprow_;
        const int taskj = (j / nb_) % npcol_;
        flag            = (taski == myrow_ && taskj == mycol_);
    }
    if (flag)
    {

        const int ib = i / (nprow_ * mb_);
        const int jb = j / (npcol_ * nb_);
        const int x  = i % mb_;
        const int y  = j % nb_;
        assert(ib * mb_ + x + (jb * nb_ + y) * mloc_ < size_);

        val_[ib * mb_ + x + (jb * nb_ + y) * mloc_] = val;
    }
}
////////////////////////////////////////////////////////////////////////////////
template <class T>
MPI_Comm DistMatrix<T>::comm_global() const
{
    return comm_global_;
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
ostream& operator<<(ostream& os, const DistMatrix<T>& a)
{
    a.print(os);
    return os;
}

template ostream& operator<<(ostream& os, const DistMatrix<double>&);
template ostream& operator<<(ostream& os, const DistMatrix<float>&);

#ifdef __KCC
// explicit instantiation declaration
template void DistMatrix<double>::clear(void);
template DistMatrix<double>::DistMatrix(
    const BlacsContext&, const int, const int, const int, const int);
template DistMatrix<double>::DistMatrix(const DistMatrix&);
template DistMatrix<double>::DistMatrix(const BlacsContext&,
    const double* const, const int, const int, const int, const int);

template void DistMatrix<double>::resize(
    const int m, const int n, const int mb, const int nb);
template void DistMatrix<double>::init(const double* const a, const int lda);
template void DistMatrix<double>::trset(const char);
template void DistMatrix<double>::print(
    ostream&, const int, const int, const int, const int) const;
template void DistMatrix<double>::print(ostream&) const;
template class vector<int>;
#else
template class DistMatrix<double>;
template class DistMatrix<float>;
#endif

} // namespace
