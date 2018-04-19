// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "hdf_tools.h"
#include "IonData.h"

#include <iostream>
#include <cassert>
#include <string.h>

using namespace std;

namespace mgmol_tools
{
void write1d(hid_t file_id, string datasetname, vector<int>& data, int length)
{
    assert( file_id>=0 );
    
    if( length==0 )return;
   
    hsize_t dim = length;
 
    // Create the data space for the datasets
    hid_t    dataspace_id = H5Screate_simple(1, &dim, NULL);
    if( dataspace_id<0 )
    {
        cerr<<"write1d(), H5Screate_simple failed!!!"<<endl;
        return;
    }

    // Create dataset
    hid_t dataset_id = H5Dcreate2(file_id, datasetname.c_str(),
                                  H5T_NATIVE_INT, 
                                  dataspace_id,
                                  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    if( dataset_id<0 )
    {
        cerr<<"write1d(), H5Dcreate2 failed!!!"<<endl;
        return;
    }
    H5Sclose(dataspace_id);
    
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, &data[0]);
    if( status<0 )
    {
        cerr<<"write1d(), H5Dwrite failed!!!"<<endl;
        return;
    }
    
    status=H5Dclose(dataset_id);
    if( status<0 )
    {
        cerr<<"write1d(), H5Dclose failed!!!"<<endl;
        return;
    }
}

void write2d(hid_t file_id, string datasetname, vector<int>& data, int* dims)
{
    assert( file_id>=0 );
    
    //cout<<"Write "<<dims[0]<<" atomic numbers..."<<endl;
    
    if( dims[0]==0 )return;
    
    hsize_t dimsm[2]={(hsize_t)dims[0],(hsize_t)dims[1]};

    // Create the data space for the datasets
    hid_t    dataspace_id = H5Screate_simple(2, dimsm, NULL);
    if( dataspace_id<0 )
    {
        cerr<<"write2d(), H5Screate_simple failed!!!"<<endl;
        return;
    }

    // Create dataset
    hid_t dataset_id = H5Dcreate2(file_id, datasetname.c_str(),
                                  H5T_NATIVE_INT, 
                                  dataspace_id,
                                  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    if( dataset_id<0 )
    {
        cerr<<"write2d(), H5Dcreate2 failed!!!"<<endl;
        return;
    }
    H5Sclose(dataspace_id);
    
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, &data[0]);
    if( status<0 )
    {
        cerr<<"write2d(), H5Dwrite failed!!!"<<endl;
        return;
    }
    
    status=H5Dclose(dataset_id);
    if( status<0 )
    {
        cerr<<"write2d(), H5Dclose failed!!!"<<endl;
        return;
    }
}

void write2d(hid_t file_id, string datasetname, vector<unsigned short>& data, int* dims)
{
    assert( file_id>=0 );
    
    //cout<<"Write "<<dims[0]<<" atomic numbers..."<<endl;
    
    if( dims[0]==0 )return;
    
    hsize_t dimsm[2]={(hsize_t)dims[0],(hsize_t)dims[1]};

    // Create the data space for the datasets
    hid_t    dataspace_id = H5Screate_simple(2, dimsm, NULL);
    if( dataspace_id<0 )
    {
        cerr<<"write2d(), H5Screate_simple failed!!!"<<endl;
        return;
    }

    // Create dataset
    hid_t dataset_id = H5Dcreate2(file_id, datasetname.c_str(),
                                  H5T_NATIVE_USHORT, 
                                  dataspace_id,
                                  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    if( dataset_id<0 )
    {
        cerr<<"write2d(), H5Dcreate2 failed!!!"<<endl;
        return;
    }
    H5Sclose(dataspace_id);
    
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_USHORT, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, &data[0]);
    if( status<0 )
    {
        cerr<<"write2d(), H5Dwrite failed!!!"<<endl;
        return;
    }
    
    status=H5Dclose(dataset_id);
    if( status<0 )
    {
        cerr<<"write2d(), H5Dclose failed!!!"<<endl;
        return;
    }
}

void write2d(hid_t file_id, string datasetname, vector<double>& data, int* dims)
{
    assert( file_id>=0 );
    
    //cout<<"Write "<<dims[0]<<" atomic numbers..."<<endl;
    
    if( dims[0]==0 )return;
    
    hsize_t dimsm[2]={(hsize_t)dims[0],(hsize_t)dims[1]};

    // Create the data space for the datasets
    hid_t    dataspace_id = H5Screate_simple(2, dimsm, NULL);
    if( dataspace_id<0 )
    {
        cerr<<"write2d(), H5Screate_simple failed!!!"<<endl;
        return;
    }

    // Create dataset
    hid_t dataset_id = H5Dcreate2(file_id, datasetname.c_str(),
                                  H5T_NATIVE_DOUBLE, 
                                  dataspace_id,
                                  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    if( dataset_id<0 )
    {
        cerr<<"write2d(), H5Dcreate2 failed!!!"<<endl;
        return;
    }
    H5Sclose(dataspace_id);
    
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, 
                             H5S_ALL, H5P_DEFAULT, &data[0]);
    if( status<0 )
    {
        cerr<<"write2d(), H5Dwrite failed!!!"<<endl;
        return;
    }
    
    status=H5Dclose(dataset_id);
    if( status<0 )
    {
        cerr<<"write2d(), H5Dclose failed!!!"<<endl;
        return;
    }
}

void write2d(hid_t file_id, string datasetname, vector<string>& data, int* dims)
{
    assert( file_id>=0 );
    
    // create type for strings of length IonData_MaxStrLength
    hid_t strtype =  H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, IonData_MaxStrLength);
    
    //cout<<"Write "<<dims[0]<<" atomic numbers..."<<endl;
    
    if( dims[0]==0 )return;
    
    hsize_t dimsm[2]={(hsize_t)dims[0],(hsize_t)dims[1]};

    // Create the data space for the datasets
    hid_t    dataspace_id = H5Screate_simple(2, dimsm, NULL);
    if( dataspace_id<0 )
    {
        cerr<<"write2d(), H5Screate_simple failed!!!"<<endl;
        return;
    }

    // Create dataset
    hid_t dataset_id = H5Dcreate2(file_id, datasetname.c_str(),
                                  strtype, 
                                  dataspace_id,
                                  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    if( dataset_id<0 )
    {
        cerr<<"write2d(), H5Dcreate2 failed!!!"<<endl;
        return;
    }
    H5Sclose(dataspace_id);
    
    // First copy the contents of the vector into a temporary container
    vector<FixedLengthString> tc;
    for (vector<string>::const_iterator i = data.begin (), end = data.end ();
                                        i != end;
                                      ++i)
    {
        FixedLengthString t;
        strncpy (t.mystring, i->c_str (), IonData_MaxStrLength);
        tc.push_back (t);
    }
    
    string attname("String_Length");
    hsize_t dimsA[1]={1};
    hid_t dataspaceA_id = H5Screate_simple(1, dimsA, NULL);
    hid_t attribute_id = H5Acreate2(dataset_id, attname.c_str(), H5T_NATIVE_INT, 
                                    dataspaceA_id, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = H5Awrite(attribute_id, H5T_NATIVE_INT, &IonData_MaxStrLength);
    if( status<0 )
    {
        cerr<<"write2d(), Attribute: "<<attname<<" --- H5Awrite failed!!!"<<endl;
    }
    
    status = H5Dwrite(dataset_id, strtype, H5S_ALL, 
                             H5S_ALL, H5P_DEFAULT, &tc[0]);
    if( status<0 )
    {
        cerr<<"write2d(), H5Dwrite failed!!!"<<endl;
        return;
    }
    
    H5Tclose (strtype);
    
    status = H5Sclose(dataspaceA_id);
    status = H5Aclose(attribute_id);
    status=H5Dclose(dataset_id);
    if( status<0 )
    {
        cerr<<"write2d(), H5Dclose failed!!!"<<endl;
        return;
    }
}

void parallelWrite2d(hid_t file_id, string datasetname, vector<int>& data, int* dims, MPI_Comm comm)
{
    assert( file_id>=0 );
    assert( !data.empty() );
    
    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
   
    hsize_t dimsm[2]={(hsize_t)dims[0],(hsize_t)dims[1]};
    hsize_t dimsf[2]={(hsize_t)(dimsm[0]*mpi_size),dimsm[1]};

    hid_t filespace = H5Screate_simple(2, dimsf, NULL);
    if( filespace<0 )
    {
        cerr<<"parallelWrite2d(), H5Screate_simple failed for filespace!!!"<<endl;
        return;
    }
    hid_t memspace  = H5Screate_simple(2, dimsm, NULL);
    if( memspace<0 )
    {
        cerr<<"parallelWrite2d(), H5Screate_simple failed for memspace!!!"<<endl;
        return;
    }


    // Create dataset
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 2, dimsm);
    hid_t dset_id = H5Dcreate2(file_id, datasetname.c_str(),
                               H5T_NATIVE_INT, filespace,
                               H5P_DEFAULT, plist_id, H5P_DEFAULT);
    if( dset_id<0 )
    {
        cerr<<"parallelWrite2d() for dataset "<<datasetname<<", H5Dcreate2() failed!!!"<<endl;
        return;
    }
    H5Pclose(plist_id);
    
    hsize_t offset[2]={mpi_rank*dimsm[0],0};
    hsize_t stride[2]={1,1};
    hsize_t count[2]={1,1};
    hsize_t block[2]={dimsm[0],dimsm[1]};
    
    /* Select hyperslab in the file. */
    herr_t status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
    if( status<0 )
    {
        cerr<<"parallelWrite2d(), H5Sselect_hyperslab() failed!!!"<<endl;
        return;
    }
    
    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &data[0]);
    if( status<0 )
    {
        cerr<<"parallelWrite2d(), H5Dwrite failed!!!"<<endl;
        return;
    }
    
    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
}

void parallelWrite2d(hid_t file_id, string datasetname, vector<unsigned short>& data, int* dims, MPI_Comm comm)
{
    assert( file_id>=0 );
    assert( !data.empty() );
    
    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
   
    hsize_t dimsm[2]={(hsize_t)dims[0],(hsize_t)dims[1]};
    hsize_t dimsf[2]={(hsize_t)(dimsm[0]*mpi_size),dimsm[1]};

    hid_t filespace = H5Screate_simple(2, dimsf, NULL);
    if( filespace<0 )
    {
        cerr<<"parallelWrite2d(), H5Screate_simple failed for filespace!!!"<<endl;
        return;
    }
    hid_t memspace  = H5Screate_simple(2, dimsm, NULL);
    if( memspace<0 )
    {
        cerr<<"parallelWrite2d(), H5Screate_simple failed for memspace!!!"<<endl;
        return;
    }


    // Create dataset
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 2, dimsm);
    hid_t dset_id = H5Dcreate2(file_id, datasetname.c_str(),
                               H5T_NATIVE_USHORT, filespace,
                               H5P_DEFAULT, plist_id, H5P_DEFAULT);
    if( dset_id<0 )
    {
        cerr<<"parallelWrite2d() for dataset "<<datasetname<<", H5Dcreate2() failed!!!"<<endl;
        return;
    }
    H5Pclose(plist_id);
    
    hsize_t offset[2]={mpi_rank*dimsm[0],0};
    hsize_t stride[2]={1,1};
    hsize_t count[2]={1,1};
    hsize_t block[2]={dimsm[0],dimsm[1]};
    
    /* Select hyperslab in the file. */
    herr_t status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
    if( status<0 )
    {
        cerr<<"parallelWrite2d(), H5Sselect_hyperslab() failed!!!"<<endl;
        return;
    }
    
    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
    status = H5Dwrite(dset_id, H5T_NATIVE_USHORT, memspace, filespace, plist_id, &data[0]);
    if( status<0 )
    {
        cerr<<"parallelWrite2d(), H5Dwrite failed!!!"<<endl;
        return;
    }
    
    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
}

void parallelWrite2d(hid_t file_id, string datasetname, vector<double>& data, int* dims, MPI_Comm comm)
{
    assert( file_id>=0 );
    assert( !data.empty() );
    
    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
   
    hsize_t dimsm[2]={(hsize_t)dims[0],(hsize_t)dims[1]};
    hsize_t dimsf[2]={(hsize_t)(dimsm[0]*mpi_size),dimsm[1]};

    hid_t filespace = H5Screate_simple(2, dimsf, NULL);
    if( filespace<0 )
    {
        cerr<<"parallelWrite2d(), H5Screate_simple failed for filespace!!!"<<endl;
        return;
    }
    hid_t memspace  = H5Screate_simple(2, dimsm, NULL);
    if( memspace<0 )
    {
        cerr<<"parallelWrite2d(), H5Screate_simple failed for memspace!!!"<<endl;
        return;
    }


    // Create dataset
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 2, dimsm);
    hid_t dset_id = H5Dcreate2(file_id, datasetname.c_str(),
                               H5T_NATIVE_DOUBLE, filespace,
                               H5P_DEFAULT, plist_id, H5P_DEFAULT);
    if( dset_id<0 )
    {
        cerr<<"parallelWrite2d() for dataset "<<datasetname<<", H5Dcreate2() failed!!!"<<endl;
        return;
    }
    H5Pclose(plist_id);
    
    hsize_t offset[2]={mpi_rank*dimsm[0],0};
    hsize_t stride[2]={1,1};
    hsize_t count[2]={1,1};
    hsize_t block[2]={dimsm[0],dimsm[1]};
    
    /* Select hyperslab in the file. */
    herr_t status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
    if( status<0 )
    {
        cerr<<"parallelWrite2d(), H5Sselect_hyperslab() failed!!!"<<endl;
        return;
    }
    
    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &data[0]);
    if( status<0 )
    {
        cerr<<"parallelWrite2d(), H5Dwrite failed!!!"<<endl;
        return;
    }
    
    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
}

void parallelWrite2d(hid_t file_id, string datasetname, vector<string>& data, int* dims, MPI_Comm comm)
{
    assert( file_id>=0 );
    assert( !data.empty() );
    
    // create type for strings of length IonData_MaxStrLength
    hid_t strtype =  H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, IonData_MaxStrLength);
    
    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
   
    hsize_t dimsm[2]={(hsize_t)dims[0],(hsize_t)dims[1]};
    hsize_t dimsf[2]={(hsize_t)(dimsm[0]*mpi_size),dimsm[1]};

    hid_t filespace = H5Screate_simple(2, dimsf, NULL);
    if( filespace<0 )
    {
        cerr<<"parallelWrite2d(), H5Screate_simple failed for filespace!!!"<<endl;
        return;
    }
    hid_t memspace  = H5Screate_simple(2, dimsm, NULL);
    if( memspace<0 )
    {
        cerr<<"parallelWrite2d(), H5Screate_simple failed for memspace!!!"<<endl;
        return;
    }


    // Create dataset
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 2, dimsm);
    hid_t dset_id = H5Dcreate2(file_id, datasetname.c_str(),
                               strtype, filespace,
                               H5P_DEFAULT, plist_id, H5P_DEFAULT);
    if( dset_id<0 )
    {
        cerr<<"parallelWrite2d() for dataset "<<datasetname<<", H5Dcreate2() failed!!!"<<endl;
        return;
    }
    H5Pclose(plist_id);
    
    hsize_t offset[2]={mpi_rank*dimsm[0],0};
    hsize_t stride[2]={1,1};
    hsize_t count[2]={1,1};
    hsize_t block[2]={dimsm[0],dimsm[1]};
    
    /* Select hyperslab in the file. */
    herr_t status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
    if( status<0 )
    {
        cerr<<"parallelWrite2d(), H5Sselect_hyperslab() failed!!!"<<endl;
        return;
    }
    
    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
    // First copy the contents of the vector into a temporary container
    vector<FixedLengthString> tc;
    for (vector<string>::const_iterator i = data.begin (), end = data.end ();
                                        i != end;
                                      ++i)
    {
        FixedLengthString t;
        strncpy (t.mystring, i->c_str (), IonData_MaxStrLength);
        tc.push_back (t);
    }
    status = H5Dwrite(dset_id, strtype, memspace, filespace, plist_id, &tc[0]);
    if( status<0 )
    {
        cerr<<"parallelWrite2d(), H5Dwrite failed!!!"<<endl;
        return;
    }
    
    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
}

void addAttribute2Dataset(hid_t dset_id, const char* attname, 
                          const vector<double>&  attr_data)
{
    assert( dset_id>-1 );
    
    // Create the data space for the attribute.
    hsize_t dim = (hsize_t)attr_data.size();

    //  Open a dataset attribute.
    hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
    if( dataspace_id<0 )
    {
        cerr<<"H5Screate failed!!!"<<endl;
        return;
    }
    hid_t attribute_id = H5Acreate2(dset_id, attname, H5T_NATIVE_DOUBLE, 
                               dataspace_id,
                               H5P_DEFAULT,H5P_DEFAULT);
    if( attribute_id<0 )
    {
        cerr<<"H5Acreate failed!!!"<<endl;
        return;
    }

    herr_t status = H5Sclose(dataspace_id);
    if( status<0 ) cerr<<"H5Sclose failed!!!"<<endl;
    
    //(*MPIdata::sout)<<"Write attribute "<<attname<<endl;
    status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &attr_data[0]);
    if( status<0 ) cerr<<"H5Awrite failed!!!"<<endl;

    status = H5Aclose(attribute_id);
    if( status<0 ) cerr<<"H5Aclose failed!!!"<<endl;
}

void addAttribute2Dataset(hid_t dset_id, const char* attname, 
                          const vector<int>&  attr_data)
{
    assert( dset_id>-1 );
    
    // Create the data space for the attribute.
    hsize_t dim = (hsize_t)attr_data.size();

    //  Open a dataset attribute.
    hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
    if( dataspace_id<0 )
    {
        cerr<<"H5Screate failed!!!"<<endl;
        return;
    }
    hid_t attribute_id = H5Acreate2(dset_id, attname, H5T_NATIVE_INT, 
                               dataspace_id,
                               H5P_DEFAULT,H5P_DEFAULT);
    if( attribute_id<0 )
    {
        cerr<<"H5Acreate failed!!!"<<endl;
        return;
    }

    herr_t status = H5Sclose(dataspace_id);
    if( status<0 )
    {
        cerr<<"H5Sclose failed!!!"<<endl;
        return;
    }
    
    //(*MPIdata::sout)<<"Write attribute "<<attname<<endl;
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, &attr_data[0]);
    if( status<0 )
    {
        cerr<<"H5Awrite failed!!!"<<endl;
        return;
    }

    status = H5Aclose(attribute_id);
    if( status<0 )
    {
        cerr<<"H5Aclose failed!!!"<<endl;
    }
}

}