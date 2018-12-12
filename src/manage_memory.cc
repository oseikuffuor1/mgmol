#include "Control.h"
#include "BlockVector.h"
#include "global.h"

// Increase memory slots in BlockVector as needed based on runtime
// options
void increaseMemorySlotsForOrbitals()
{
    Control& ct     = *(Control::instance());

    switch (ct.OuterSolver())
    {
        case OuterSolverType::ABPG:
        {
            // r_k-1, phi_k-1
            BlockVector<ORBDTYPE>::incMaxAllocInstances(2);
            break;
        }
        case OuterSolverType::PolakRibiere:
        {
            // r_k-1, z_k, z_k-1, p_k
            BlockVector<ORBDTYPE>::incMaxAllocInstances(4);
            break;
        }

    }

    switch (ct.WFExtrapolation())
    {
        case WFExtrapolationType::Reversible:
        {
            BlockVector<ORBDTYPE>::incMaxAllocInstances(2);
            break;
        }
        case WFExtrapolationType::Order2:
        {
            BlockVector<ORBDTYPE>::incMaxAllocInstances(1);
            break;
        }
        case WFExtrapolationType::Order3:
        {
            BlockVector<ORBDTYPE>::incMaxAllocInstances(2);
            break;
        }

    }

    for (short i = 1; i < ct.wf_m; i++)
        BlockVector<ORBDTYPE>::incMaxAllocInstances(2);
    if (ct.use_kernel_functions) BlockVector<ORBDTYPE>::incMaxAllocInstances(1);
}

