// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id:$
#ifndef OrbitalsExtrapolationFACTORY_H
#define OrbitalsExtrapolationFACTORY_H

#include "OrbitalsExtrapolationOrder2.h"
#include "OrbitalsExtrapolationOrder3.h"
#include "SpreadPenalty.h"

class OrbitalsExtrapolationFactory
{
public:

    static OrbitalsExtrapolation* create(const int type, SpreadPenalty *spread_penalty)
    {
        OrbitalsExtrapolation* orbitals_extrapol;
        switch( type )
        {
            case 1:
                orbitals_extrapol = new OrbitalsExtrapolationOrder2();
                break;
            case 2:
                orbitals_extrapol = new OrbitalsExtrapolationOrder3();
                break;
            default:
                (*MPIdata::serr)<<"OrbitalsExtrapolation* create() --- option invalid:"
                               <<type<<endl;
                exit(2);
        }
        return orbitals_extrapol;
    }
};

#endif