// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#include "Rho.h"
#include "LocGridOrbitals.h"
#include "Control.h"
#include "Mesh.h"
#include "MPIdata.h"
#include "SubMatrices.h"
#include "ProjectedMatrices.h"
#include "SquareLocalMatrices.h"
#include "mputils.h"
#include "numerical_kernels.h"
using namespace std;

Timer Rho::update_tm_("Rho::update");
Timer Rho::compute_tm_("Rho::compute");

//class ProjectedMatricesInterface;

Rho::Rho()
    : orbitals_type_(-1),
      iterative_index_(-10),
      verbosity_level_(0) // default value
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    myspin_=mmpi.myspin();
    nspin_ =mmpi.nspin();

    //default values for block sizes
    block_functions_=4;
    block_space_    =128;

    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();
    np_=mygrid.size();
    
    rho_.resize(nspin_);
    for(short is=0;is<nspin_;is++){
        rho_[is].resize(np_);
        memset(&rho_[is][0],0,np_*sizeof(RHODTYPE));
    }

    assert( nspin_==1 || nspin_==2 );
    assert( myspin_==0 || myspin_==1 );
};

Rho::~Rho()
{}

void Rho::setup(const int orbitals_type,
                const vector<vector<int> >& orbitals_indexes)
{
    if( verbosity_level_>2 && onpe0 )
        (*MPIdata::sout)<<" Rho::setup()"<<endl;

    orbitals_type_=orbitals_type;

    orbitals_indexes_=orbitals_indexes;
}

    
void Rho::extrapolate()
{
    double minus=-1;
//    int ione=1;
    double two=2.;
    if( rho_minus1_.empty() )
    {
        rho_minus1_.resize(nspin_);
        rho_minus1_[myspin_].resize(np_);
        memcpy(&rho_minus1_[myspin_][0],&rho_[myspin_][0],np_*sizeof(RHODTYPE));
        return;
    }
    RHODTYPE* tmp=new RHODTYPE[np_];
    memcpy(tmp,&rho_[myspin_][0],np_*sizeof(RHODTYPE));

    MPscal(np_,two,&rho_[myspin_][0]);
    MPaxpy(np_,minus,&rho_minus1_[myspin_][0],&rho_[myspin_][0]);

    memcpy(&rho_minus1_[myspin_][0],tmp,np_*sizeof(RHODTYPE));

    delete[] tmp;
}

void Rho::axpyRhoc(const double alpha, RHODTYPE *rhoc)
{
//    int ione=1;
    
    double factor = ( nspin_ > 1 ) ? 0.5*alpha : alpha;
    MPaxpy(np_,factor,&rhoc[0],&rho_[myspin_][0]);
}

void Rho::update(LocGridOrbitals& current_orbitals)
{
    const ProjectedMatricesInterface& proj_matrices( *(current_orbitals.getProjMatrices()) );

    assert( current_orbitals.getIterativeIndex()>=0 );

    update_tm_.start();

    if( verbosity_level_>2 && onpe0 )
        (*MPIdata::sout)<<"Rho::update()"<<endl; 

    const int new_iterative_index
        = ((1+current_orbitals.getIterativeIndex())%100)
         +(proj_matrices.getDMMatrixIndex()%100)*100;

    if( iterative_index_==new_iterative_index ){
        if( onpe0 && verbosity_level_>2)
            (*MPIdata::sout)<<"Rho already up to date, iterative_index_="
                <<iterative_index_
                <<endl;
        return;
    }
    iterative_index_=new_iterative_index;
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        (*MPIdata::sout)<<"Rho::update(), iterative_index_="<<iterative_index_<<endl;
#endif
    
    computeRho(current_orbitals);

    gatherSpin();

    rescaleTotalCharge();
    
    update_tm_.stop();
}

// note: rho can be negative because of added background charge
double Rho::computeTotalCharge()
{
    const int nspin=(int)rho_.size();

    double tcharge = 0.;
    for(short ispin=0;ispin<nspin;ispin++)
    {
        const RHODTYPE* const prho=&rho_[ispin][0];
        for(int idx = 0;idx < np_;idx++)
        {
            assert(prho[idx]<1.e6);
            tcharge += (double)prho[idx];
        }
    } 

    // reduce over spin communicator since charge from other spin is already included
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&tcharge,1,MPI_SUM);

    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();
    tcharge *= mygrid.vel();
#ifdef DEBUG
    if( mmpi.instancePE0() )
        (*MPIdata::sout)<<fixed<<setprecision(8)<<" myspin="<<myspin_<<", Total charge = "<<tcharge<<endl;
#endif

    return tcharge;
}

void Rho::rescaleTotalCharge()
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    // Check total charge
    const double tcharge = computeTotalCharge();
    
    const int nspin=(int)rho_.size();

    Control& ct = *(Control::instance());
    double nel=ct.getNel();
    if( tcharge>0. ){
        double t1=nel/tcharge;

        if( mmpi.PE0() && fabs(t1-1.)>0.001)
        {
            (*MPIdata::sout)<<" Rho::rescaleTotalCharge(), charge = "<<tcharge<<endl;
            (*MPIdata::sout)<<" Rescaling factor: "<<t1<<endl;
            (*MPIdata::sout)<<" Num. electrons: "<<nel<<endl;
        }

        for(int ispin=0;ispin<nspin;ispin++)
            MPscal(np_, t1, &rho_[ispin][0]);
    }
#ifdef DEBUG
    Mesh* mymesh = Mesh::instance();
    // Check total charge again
    double ttcharge = 0.;
    for(int ispin=0;ispin<nspin;ispin++)
    for(int idx = 0;idx < np_;idx++) {
        ttcharge += (double)rho_[ispin][idx];
    } 

    mmpi.allreduce(&ttcharge,1,MPI_SUM);
    const pb::Grid& mygrid  = mymesh->grid();
    ttcharge *= mygrid.vel();

    if( mmpi.instancePE0() )
        (*MPIdata::sout)<<fixed<<setprecision(8)<<" Total charge after rescaling: "<<ttcharge<<endl;
#endif
}

int Rho::setupSubdomainData(const int iloc,
                            const vector<const LocGridOrbitals*>& vorbitals,
                            const ProjectedMatricesInterface* const projmatrices,
                            vector<MATDTYPE>& melements,
                            vector<vector<const ORBDTYPE*> >& vmpsi)
{
    //printWithTimeStamp("Rho::setupSubdomainData()...",cout);
    
    const short norb=(short)vorbitals.size();
    vmpsi.resize(norb);
    
    const vector<int>& loc_indexes( orbitals_indexes_[iloc] );
    const int n_colors=(int)loc_indexes.size();
    
    if( n_colors==0 )return 0;

    vector<int> mycolors;
    mycolors.reserve(n_colors);
    for(int icolor = 0;icolor < n_colors;icolor++)
        if(loc_indexes[icolor]!=-1)
            mycolors.push_back(icolor);
    const int nmycolors=(int)mycolors.size();
    
    for(short j=0;j<norb;j++)
        vmpsi[j].resize(nmycolors);
            
    short j=0;
    for(vector<const LocGridOrbitals*>::const_iterator it =vorbitals.begin();
                                                       it!=vorbitals.end();
                                                     ++it)
    {
        for(int color=0;color<nmycolors;color++)
        {
            vmpsi[j][color]=(*it)->getPsi(mycolors[color],iloc);
            assert( vmpsi[j][color]!=NULL );
        }
        j++;
    }

#ifndef USE_DIS_MAT
    const int numst    =vorbitals[0]->numst();
#else
    SquareLocalMatrices<MATDTYPE>& localX( projmatrices->getLocalX() );
    const MATDTYPE* const localX_iloc=localX.getSubMatrix(iloc);
#endif
    melements.clear();
    melements.resize(nmycolors*nmycolors);
    
    for(int i=0;i<nmycolors;i++)
    {
        const int icolor=mycolors[i];
        if( norb==1 )
        {
#ifdef USE_DIS_MAT
            melements[i*nmycolors+i]=localX_iloc[icolor+n_colors*icolor];
#else
            const int    iist    = loc_indexes[icolor]*numst;
            const int    index_X = loc_indexes[icolor]*(numst+1);
            melements[i*nmycolors+i] = dm_.getVal(index_X);
#endif
        }
        const int jmax = ( norb==1 ) ? i : nmycolors;
        for(int j=0;j<jmax;j++)
        {
            int jcolor=mycolors[j];
#ifdef USE_DIS_MAT
            melements[j*nmycolors+i]=localX_iloc[icolor+n_colors*jcolor];
#else
            melements[j*nmycolors+i]=dm_.getVal(iist + loc_indexes[jcolor]);
#endif
        }                
    }
    
    return nmycolors;
}

void Rho::accumulateCharge(
    const double alpha,
    const short ix_max,
    const ORBDTYPE* const psii,
    const ORBDTYPE* const psij,
    RHODTYPE* const plrho)
{
    for(int ix = 0;ix < ix_max;ix++)
        plrho[ix] += (RHODTYPE)(alpha
                   * (double)psii[ix]
                   * (double)psij[ix]);
}

void Rho::computeRhoSubdomain(const int iloc_init, 
                             const int iloc_end,
                             const LocGridOrbitals& orbitals)
{
    assert( orbitals_type_==0 || orbitals_type_==1 );
    
    compute_tm_.start();
    
    Mesh* mymesh = Mesh::instance();

    const int loc_numpt=mymesh->locNumpt();
    
    RHODTYPE* const prho=&rho_[myspin_][0];
    
    vector<const LocGridOrbitals*> vorbitals;
    vorbitals.push_back(&orbitals);
    
    const double nondiagfactor = 2.;

    for(int iloc=iloc_init;iloc<iloc_end;iloc++)
    {
        const int istart = iloc*loc_numpt;
        
        RHODTYPE* const lrho  = &prho[istart];
        
        vector<MATDTYPE> melements;
        vector< vector<const ORBDTYPE*> > vmpsi;
        
        const int nmycolors=setupSubdomainData(iloc, vorbitals, orbitals.projMatrices(),
                                               melements, vmpsi);
        assert( vmpsi.size()==1 );
        vector<const ORBDTYPE*> mpsi=vmpsi[0];

        const int nblocks_color=nmycolors/block_functions_;            
        const int max_icolor= ( nmycolors%block_functions_ == 0 )
                              ? nmycolors
                              : nblocks_color*block_functions_;
        const int missed_rows=nmycolors-max_icolor;
        //(*MPIdata::sout)<<"max_icolor="<<max_icolor<<endl;
        //(*MPIdata::sout)<<"Rho::computeRhoSubdomain: missed_rows="<<missed_rows<<endl;
        //(*MPIdata::sout)<<"nblocks_color="<<nblocks_color<<endl;

#ifdef _OPENMP    
#pragma omp parallel for
#endif
        for(int idx = 0;idx < loc_numpt;idx+=block_space_)
        {
            const short ix_max=min(block_space_,loc_numpt-idx);
            //RHODTYPE* const plrho=lrho+idx;
            
            // non-diagonal blocks
            for(int icolor=0;icolor<max_icolor;icolor+=block_functions_)
            for(int jcolor=0;jcolor<icolor;    jcolor+=block_functions_){
                
                nonOrthoRhoKernel(icolor,block_functions_,
                                 jcolor,block_functions_,
                                 idx,ix_max,
                                 &melements[0],nmycolors,
                                 mpsi,nondiagfactor,
                                 lrho);
            }
            // finish missing rows ()
            if( missed_rows )
            for(int jcolor=0;jcolor<max_icolor;jcolor+=block_functions_){
                nonOrthoRhoKernel(max_icolor,missed_rows,
                                 jcolor,block_functions_,
                                 idx,ix_max,
                                 &melements[0],nmycolors,
                                 mpsi,nondiagfactor,
                                 lrho);
            }
            
            // diagonal blocks (jcolor=icolor)
            for(int icolor=0;icolor<max_icolor;icolor+=block_functions_)
            {
                nonOrthoRhoKernelDiagonalBlock(icolor,block_functions_,
                                 idx,ix_max,
                                 &melements[0],nmycolors,
                                 mpsi,
                                 lrho);
            }
            // finish missing rows () for diagonal blocks
            if( missed_rows )
            {
                nonOrthoRhoKernelDiagonalBlock(max_icolor,missed_rows,
                                 idx,ix_max,
                                 &melements[0],nmycolors,
                                 mpsi,
                                 lrho);
            }
        }
                       
    }

    compute_tm_.stop();
}

void Rho::computeRhoSubdomainOffDiagBlock(const int iloc_init, 
                                         const int iloc_end,
                                         const vector<const LocGridOrbitals*>& vorbitals,
                                         const ProjectedMatricesInterface* const projmatrices)
{
    assert( orbitals_type_==0 || orbitals_type_==1 );
    assert( vorbitals.size()==2 );
    
    //printWithTimeStamp("Rho::computeRhoSubdomainOffDiagBlock()...",cout);
    
    Mesh* mymesh = Mesh::instance();

    const int loc_numpt=mymesh->locNumpt();
    
    RHODTYPE* const prho=&rho_[myspin_][0];
    
    for(int iloc=iloc_init;iloc<iloc_end;iloc++)
    {
        const int istart = iloc*loc_numpt;
        
        RHODTYPE* const lrho  = &prho[istart];
        
        vector<MATDTYPE> melements;
        vector< vector<const ORBDTYPE*> > vmpsi;
        
        const int nmycolors=setupSubdomainData(iloc, vorbitals, projmatrices,
                                               melements, vmpsi);
        const int nblocks_color=nmycolors/block_functions_;            
        const int max_icolor= ( nmycolors%block_functions_ == 0 )
                              ? nmycolors
                              : nblocks_color*block_functions_;
        const int missed_rows=nmycolors-max_icolor;
        //(*MPIdata::sout)<<"max_icolor="<<max_icolor<<endl;
        //(*MPIdata::sout)<<"Rho::computeRhoSubdomainOffDiagBlock: missed_rows="<<missed_rows<<endl;
        //(*MPIdata::sout)<<"nblocks_color="<<nblocks_color<<endl;
        
        for(int idx = 0;idx < loc_numpt;idx+=block_space_)
        {
            const short ix_max=min(block_space_,loc_numpt-idx);
            RHODTYPE* const plrho=lrho+idx;
            
            for(int istart=0;istart<max_icolor;istart+=block_functions_)
            for(int jstart=0;jstart<max_icolor;jstart+=block_functions_){
                for(short jcolor=jstart;jcolor<jstart+block_functions_;jcolor++){
                    const int jld=jcolor*nmycolors;
                    const ORBDTYPE* const psij=&vmpsi[1][jcolor][idx];
                    for(short icolor=istart;icolor<istart+block_functions_;icolor++){
                        const ORBDTYPE* const psii=&vmpsi[0][icolor][idx];
                        const double alpha=(double)melements[jld+icolor];
                        accumulateCharge(alpha,ix_max,psii,psij,plrho);
                    }
                }
            }
            // finish missing rows and cols
            if( missed_rows ){
                int icolor=max_icolor;
                int imax =nmycolors%block_functions_;
                for(int jcolor=0;jcolor<max_icolor;jcolor+=block_functions_){
                    for(short j=0;j<block_functions_;j++){
                        const int jstart=(jcolor+j)*nmycolors+icolor;
                        const ORBDTYPE* const psij=&vmpsi[1][jcolor+j][idx];
                        for(short i=0;i<imax;i++){
                            const ORBDTYPE* const psii=&vmpsi[0][icolor+i][idx];
                            const double alpha=(double)melements[jstart+i];
                            accumulateCharge(alpha,ix_max,psii,psij,plrho);
                        }
                    }
                }
            }
            if( missed_rows ){
                int jcolor=max_icolor;
                int jmax=nmycolors%block_functions_;
                for(int icolor=0;icolor<max_icolor;icolor+=block_functions_){
                    for(short j=0;j<jmax;j++){
                        const int jstart=(jcolor+j)*nmycolors+icolor;
                        const ORBDTYPE* const psij=&vmpsi[1][jcolor+j][idx];
                        for(short i=0;i<block_functions_;i++){
                            const ORBDTYPE* const psii=&vmpsi[0][icolor+i][idx];
                            const double alpha=(double)melements[jstart+i];
                            accumulateCharge(alpha,ix_max,psii,psij,plrho);
                        }
                    }
                }
            }
            if( missed_rows ){
                int icolor=max_icolor;
                int imax =nmycolors%block_functions_;
                int jcolor=max_icolor;
                int jmax=nmycolors%block_functions_;
                for(short j=0;j<jmax;j++){
                    const int jstart=(jcolor+j)*nmycolors+icolor;
                    const ORBDTYPE* const psij=&vmpsi[1][jcolor+j][idx];
                    for(short i=0;i<imax;i++){
                        const ORBDTYPE* const psii=&vmpsi[0][icolor+i][idx];
                        const double alpha=(double)melements[jstart+i];
                        accumulateCharge(alpha,ix_max,psii,psij,plrho);
                    }
                }
             }
            
        }
                       
    }
}

void Rho::computeRhoSubdomain(const int iloc_init, 
                             const int iloc_end,
                             const LocGridOrbitals& orbitals,
                             const vector<PROJMATDTYPE>& occ)
{
    assert( orbitals_type_==0 || orbitals_type_==1 || orbitals_type_==2 );
    if( verbosity_level_>2 && onpe0 )
        (*MPIdata::sout)<<"Rho::computeRhoSubdomain, diagonal case..."<<endl;
    
    Mesh* mymesh = Mesh::instance();

    const int loc_numpt=mymesh->locNumpt();
    const int n_colors =orbitals_indexes_[0].size();
    
    RHODTYPE* const prho=&rho_[myspin_][0];
    
    for(int iloc=iloc_init;iloc<iloc_end;iloc++)
    {
        const int istart = iloc*loc_numpt;
        vector<int>& loc_indexes=orbitals_indexes_[iloc];

        RHODTYPE* const lrho  = prho+istart;

        // Loop over states and accumulate charge
        for(int icolor = 0;icolor < n_colors;icolor++)
        {
            const ORBDTYPE* const psi=orbitals.getPsi(icolor,iloc);
            const int st=loc_indexes[icolor];
            assert( st<(int)occ.size() );
            if( st>=0 )
            {
                const double t1=2.*(double)occ[st];
                for(int idx = 0;idx < loc_numpt;idx++)
                {
                    const double alpha=(double)psi[idx];        
                    lrho[idx] += (RHODTYPE)(t1 * alpha * alpha);
                }
            }
         
        }
    }
}

void Rho::computeRho(LocGridOrbitals& orbitals)
{
    ProjectedMatricesInterface& proj_matrices( *(orbitals.getProjMatrices()) );

    computeRho(orbitals, proj_matrices);
}

void Rho::computeRho(LocGridOrbitals& orbitals, ProjectedMatricesInterface& proj_matrices)
{
    assert( rho_.size()>0 );
    assert( rho_[myspin_].size()>0 );

    Mesh* mymesh = Mesh::instance();
    const int subdivx=mymesh->subdivx();
    const int loc_numpt=mymesh->locNumpt();
    Control& ct = *(Control::instance());

    memset(&rho_[myspin_][0],0, subdivx*loc_numpt*sizeof(RHODTYPE));

    if( orbitals_type_==0 || (orbitals_type_==2 && ct.fullyOccupied()) )
    {
        vector<PROJMATDTYPE> occ( orbitals.numst() );
        proj_matrices.getOccupations(occ);
        computeRhoSubdomain(0, subdivx, orbitals, occ);
    }
    else
    {
#ifdef USE_DIS_MAT
        proj_matrices.updateSubMatX();
#else
        dm_=proj_matrices.dm();
#endif
        computeRhoSubdomain(0, subdivx, orbitals);
    }
}

void Rho::computeRho(LocGridOrbitals& orbitals1,
                     LocGridOrbitals& orbitals2,
                     const dist_matrix::DistMatrix<DISTMATDTYPE>& dm11,
                     const dist_matrix::DistMatrix<DISTMATDTYPE>& dm12,
                     const dist_matrix::DistMatrix<DISTMATDTYPE>& dm21,
                     const dist_matrix::DistMatrix<DISTMATDTYPE>& dm22)
{
    assert( orbitals_type_==1 );

    Mesh* mymesh = Mesh::instance();
    const int subdivx=mymesh->subdivx();
    const int loc_numpt=mymesh->locNumpt();

    memset(&rho_[myspin_][0],0, subdivx*loc_numpt*sizeof(RHODTYPE));

    // 11 diagonal block
#ifdef USE_DIS_MAT
    ProjectedMatrices* projmatrices1=dynamic_cast<ProjectedMatrices*>( orbitals1.getProjMatrices() );
    projmatrices1->updateSubMatX( dm11 );
#else
    dm_=dm11;
#endif
    computeRhoSubdomain(0, subdivx, orbitals1);

    // 22 diagonal block
#ifdef USE_DIS_MAT
    ProjectedMatrices* projmatrices2=dynamic_cast<ProjectedMatrices*>( orbitals2.getProjMatrices() );
    projmatrices2->updateSubMatX( dm22 );
#else
    dm_=dm22;
#endif
    computeRhoSubdomain(0, subdivx, orbitals2);
    
    // non diagonal blocks
    vector<const LocGridOrbitals*> vorbitals;
    vorbitals.push_back(&orbitals1);
    vorbitals.push_back(&orbitals2);

    dist_matrix::DistMatrix<DISTMATDTYPE> dm(dm12);
    dm.scal(2.); // use symmetry to reduce work
#ifdef USE_DIS_MAT
    projmatrices1->updateSubMatX( dm );
#else
    dm_=dm;
#endif
    computeRhoSubdomainOffDiagBlock(0, subdivx, vorbitals, projmatrices1);    
}

void Rho::computeRho(LocGridOrbitals& orbitals,
                     const dist_matrix::DistMatrix<DISTMATDTYPE>& dm)
{
    assert( orbitals_type_==1 );

    iterative_index_++;

    Mesh* mymesh = Mesh::instance();
    const int subdivx=mymesh->subdivx();
    const int loc_numpt=mymesh->locNumpt();

    memset(&rho_[myspin_][0],0, subdivx*loc_numpt*sizeof(RHODTYPE));

    ProjectedMatrices* projmatrices=dynamic_cast<ProjectedMatrices*>( orbitals.getProjMatrices() );
    projmatrices->updateSubMatX( dm );

    computeRhoSubdomain(0, subdivx, orbitals);
}

void Rho::init(const RHODTYPE* const rhoc)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int     ione=1;
    
    // Initialize the charge density  
    if( verbosity_level_>2 && mmpi.instancePE0() )
        (*MPIdata::sout)<<"Rho: Initialize electronic density value with rhoc"<<endl;
    for(int i=0;i<(int)rho_.size();i++){
        Tcopy(&np_, rhoc, &ione, &rho_[i][0], &ione);
        if(rho_.size()==2)MPscal(np_,0.5,&rho_[i][0]);
    }
    iterative_index_=0;
    
    rescaleTotalCharge();
}

void Rho::initUniform()
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    // Initialize the charge density  
    if( mmpi.instancePE0() )
        (*MPIdata::sout)<<" Initialize electronic density with uniform value"<<endl;
    for(int i=0;i<(int)rho_.size();i++){
        for(int j=0;j<np_;j++)
            rho_[i][j]=1.;
        if(rho_.size()==2)MPscal(np_,0.5,&rho_[i][0]);
    }
    iterative_index_=0;
    
    rescaleTotalCharge();
}

// read rho and potentials form a hdf5 file
int Rho::readRestart(HDFrestart& file)
{
    Control& ct = *(Control::instance());
    if( onpe0 && ct.verbose>0 )
        (*MPIdata::sout)<<"Try to read density"<<endl;
      
    // Read the Density
    file.read_1func_hdf5(&rho_[myspin_][0],"Density");
    iterative_index_=0;

    return 0;
}

template <typename T>
double Rho::dotWithRho(const T* const func)const
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
//    int ione=1;

    double val=MPdot(np_, &rho_[mmpi.myspin()][0], func);
//    double val=ddot(&np_, &rho_[mmpi.myspin()][0], &ione, func, &ione);

    double esum=0.;
    mmpi.allreduce(&val, &esum, 1, MPI_SUM);
    val=esum;

    mmpi.allreduceSpin(&val, &esum, 1, MPI_SUM);
    val=esum;

    return val;
}

void Rho::gatherSpin()
{
    if( nspin_<2 )return;
    
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    //if(mmpi.instancePE0())cout<<"Rho::gatherSpin(), myspin_="<<myspin_<<endl;
    mmpi.exchangeDataSpin(&rho_[ myspin_     ][0],
                          &rho_[(myspin_+1)%2][0],np_);
}

void Rho::printTimers(std::ostream& os)
{
   update_tm_.print(os);
   compute_tm_.print(os);
}

template double Rho::dotWithRho<double>(const double* const func)const;
#ifdef USE_MP
template double Rho::dotWithRho<float>(const float* const func)const;
#endif