// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <vector>
#include <string>
#include <mpi.h>
#include "hdf5.h"

namespace mgmol_tools
{
void write1d(hid_t file_id, std::string datasetname, std::vector<int>& data, int length);
void write2d(hid_t file_id, std::string datasetname, std::vector<int>& data, int* dims);
void write2d(hid_t file_id, std::string datasetname, std::vector<unsigned short>& data, int* dims);
void write2d(hid_t file_id, std::string datasetname, std::vector<double>& data, int* dims);
void write2d(hid_t file_id, std::string datasetname, std::vector<std::string>& data, int* dims);

void parallelWrite2d(hid_t file_id, std::string datasetname, std::vector<int>& data, int* dims, MPI_Comm comm);
void parallelWrite2d(hid_t file_id, std::string datasetname, std::vector<unsigned short>& data, int* dims, MPI_Comm comm);
void parallelWrite2d(hid_t file_id, std::string datasetname, std::vector<double>& data, int* dims, MPI_Comm comm);
void parallelWrite2d(hid_t file_id, std::string datasetname, std::vector<std::string>& data, int* dims, MPI_Comm comm);

void addAttribute2Dataset(hid_t dset_id, const char* attname, 
                          const std::vector<double>&  attr_data);
void addAttribute2Dataset(hid_t dset_id, const char* attname, 
                          const std::vector<int>&  attr_data);
}