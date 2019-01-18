// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef _PCG_SOLVER_DIEL_H_
#define _PCG_SOLVER_DIEL_H_

#include "Control.h"
#include "PB.h"
#include "PBh2.h"
#include "PBh4.h"
#include "PBh4M.h"
#include "PBh4MP.h"
#include "PBh6.h"
#include "PBh8.h"

#include <vector>

template <class T, typename T2>
class PCGSolver_Diel
{

private:
    std::vector<pb::Grid*> grid_;
    short lap_type_;
    short bc_[3];
    // operators
    T oper_;
    std::vector<T*> pc_oper_;
    std::vector<pb::GridFunc<T2>*> gf_work_;
    std::vector<pb::GridFunc<T2>*> gf_rcoarse_;
    std::vector<pb::GridFunc<T2>*> gf_newv_;
    // solver params
    int maxiters_;
    double tol_;
    double final_residual_;
    double residual_reduction_;
    // precon params
    short nu1_;
    short nu2_;
    short max_nlevels_;
    short nlevels_;
    void setupPrecon();

public:
    PCGSolver_Diel(T& oper, const short px, const short py, const short pz)
        : oper_(oper)
    {
        maxiters_           = 10; // default
        nu1_                = 2; // default
        nu2_                = 2; // default
        tol_                = 1.e-16;
        max_nlevels_        = 10;
        final_residual_     = -1.;
        residual_reduction_ = -1.;

        // boundary conditions
        bc_[0] = px;
        bc_[1] = py;
        bc_[2] = pz;
        //        fully_periodic_=( (bc_[0]==1) && (bc_[1]==1) && (bc_[2]==1) );
        Control& ct = *(Control::instance());
        lap_type_   = ct.lap_type;
    };

    void setup(const short nu1, const short nu2, const short max_sweeps,
        const double tol, const short max_nlevels)
    {
        maxiters_    = max_sweeps;
        nu1_         = nu1;
        nu2_         = nu2;
        tol_         = tol;
        max_nlevels_ = max_nlevels;
    }

    void clear();

    void preconSolve(pb::GridFunc<T2>& gf_v, const pb::GridFunc<T2>& gf_f,
        const short level = 0);

    bool solve(pb::GridFunc<T2>& gf_phi, pb::GridFunc<T2>& gf_rhs);

    bool solve(pb::GridFunc<T2>& gf_phi, pb::GridFunc<T2>& gf_rhs,
        pb::GridFunc<T2>& gf_rhod, pb::GridFunc<T2>& gf_vks);

    double getFinalResidual() const { return final_residual_; }
    double getResidualReduction() const { return residual_reduction_; }

    // Destructor
    ~PCGSolver_Diel() { clear(); }
};

#endif
