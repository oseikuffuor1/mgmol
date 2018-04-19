// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

////////////////////////////////////////////////////////////////////////////////
//
//  Timer.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Timer.h,v 1.6 2010/01/29 01:10:10 jeanluc Exp $

#ifndef TIMER_H
#define TIMER_H

#include <time.h>
#include <sys/time.h>
#include <cstring>
#include <iomanip>
#include <iostream>
#include "mpi.h"

class Timer
{
private:

    clock_t clk_;
    std::string  name_;
    double   t_;
    double   total_cpu_;
    double   total_real_;
    bool    running_;
    int     ncalls_;

    MPI_Comm comm_;

    double cpu()const
    {
      if ( running_ ) 
      {
        return total_cpu_ + ((double)(clock()-clk_))/CLOCKS_PER_SEC;
      }
      else
      {
        return total_cpu_;
      } 
    };

    double gtod(void)const;

    double real()const
    {
      if ( running_ ) 
      {
        return total_real_ + gtod()-t_;
      }
      else
      {
        return total_real_;
      } 
    };

public:

    Timer(const std::string& name, MPI_Comm comm=MPI_COMM_WORLD) :
        name_(name),
        total_cpu_(0.0),
        total_real_(0.0),
        running_(false),
        ncalls_(0),
        comm_(comm)
    {};

    void reset()
    {
        total_cpu_  = 0.0;
        total_real_ = 0.0;
        running_ = false;
        ncalls_  = 0;
    };

    void start()
    {
      clk_ = clock();
      t_   = gtod();
      running_ = true;
      ncalls_++;
    };

    bool running()const { return running_; };

    void stop()
    {
      if ( running_ ) 
      {
        total_cpu_  += ((double)(clock()-clk_))/CLOCKS_PER_SEC;
        total_real_ += gtod()-t_;
        running_ = false;
      }
    };

    void print(std::ostream& os)const;
};

#endif