/* Copyright 2019-2021
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "AnyFFT.H"

#include <AMReX.H>
// #include <AMReX_IntVect.H>
#include <AMReX_REAL.H>

// #include <fftw3.h>
#include "fftx_plan.fftx.precompile.hpp"
#include "device_macros.h"

namespace AnyFFT
{
  FFTplan CreatePlan(// const amrex::IntVect& real_size,
                     const fftx::point_t<3>& real_size,
                     amrex::Real * const real_array,
                     Complex * const complex_array,
                     const direction dir,
                     const int dim)
    {
        FFTplan fft_plan;
        fft_plan.m_plan = fftx::fftx_plan();

        // Initialize fft_plan.m_plan with the vendor fft plan.
        // Swap dimensions: AMReX FAB are Fortran-order but FFTX is C-order
        if (dir == direction::R2C){
            if (dim == 3) {
              // fftx::point_t<3> sz = fftx::point_t<3>({{
              // real_size[2], real_size[1], real_size[0]}});
              fftx::point_t<3> sz = fftx::point_t<3>({{
                    real_size[2], real_size[1], real_size[0]}});
              fft_plan.m_plan.m_tfm_3d_r2c = new fftx::mdprdft<3>(sz);
              // VendorCreatePlanR2C3D(
              // real_size[2], real_size[1], real_size[0],
              // real_array, complex_array, FFTW_ESTIMATE);
            } else if (dim == 2) { // FIXME
              // amrex::Abort("only dim=3 has been implemented");
              std::cout << "only dim=3 has been implemented" << std::endl;
              exit(-1);
            } else {
              // amrex::Abort("only dim=2 and dim=3 have been implemented");
              std::cout << "only dim=2 and dim=3 have been implemented" << std::endl;
            }
        } else if (dir == direction::C2R){
            if (dim == 3) {
              // fftx::point_t<3> sz = fftx::point_t<3>({{
              // real_size[2], real_size[1], real_size[0]}});
              fftx::point_t<3> sz = fftx::point_t<3>({{
                    real_size[2], real_size[1], real_size[0]}});
              fft_plan.m_plan.m_tfm_3d_c2r = new fftx::imdprdft<3>(sz);
              // VendorCreatePlanC2R3D(
              //   real_size[2], real_size[1], real_size[0],
              //   complex_array, real_array, FFTW_ESTIMATE);
            } else if (dim == 2) { // FIXME
              // amrex::Abort("only dim=3 has been implemented");
              std::cout << "only dim=3 has been implemented" << std::endl;
            } else {
              // amrex::Abort("only dim=2 and dim=3 have been implemented. Should be easy to add dim=1.");
              std::cout << "only dim=2 and dim=3 have been implemented. Should be easy to add dim=1." << std::endl;
            }
        }

        // Store meta-data in fft_plan
        fft_plan.m_real_array = real_array;
        fft_plan.m_complex_array = complex_array;
        fft_plan.m_dir = dir;
        fft_plan.m_dim = dim;

        return fft_plan;
    }

    void DestroyPlan(FFTplan& fft_plan)
    {
      // Call destructor on fft_plan.m_plan.
    }

    void Execute(FFTplan& fft_plan)
    {
      // First get sizes of real and complex arrays.
      fftx::point_t<3> rSize, cSize;
      if (fft_plan.m_dir == direction::R2C)
        {
          fftx::mdprdft<3>* tfm_3d_r2c = fft_plan.m_plan.m_tfm_3d_r2c;
          rSize = tfm_3d_r2c->inputSize();
          cSize = tfm_3d_r2c->outputSize();
        }
      else if (fft_plan.m_dir == direction::C2R)
        {
          fftx::imdprdft<3>* tfm_3d_c2r = fft_plan.m_plan.m_tfm_3d_c2r;
          cSize = tfm_3d_c2r->inputSize();
          rSize = tfm_3d_c2r->outputSize();
        }
      else
        {
          std::cout << "direction must be R2C or C2R" << std::endl;
          exit(-1);
        }
      fftx::box_t<3> rDomain(fftx::point_t<3>({{1, 1, 1}}),
                             fftx::point_t<3>({{rSize[0], rSize[1], rSize[2]}}));
      fftx::box_t<3> cDomain(fftx::point_t<3>({{1, 1, 1}}),
                             fftx::point_t<3>({{cSize[0], cSize[1], cSize[2]}}));

      auto real_bytes = rDomain.size() * sizeof(double); // (amrex::Real);
      auto complex_bytes = cDomain.size() * sizeof(Complex);

      // Device buffer to contain inmput array followed by output array.
      char* bufferDevicePtr;
      DEVICE_MALLOC(&bufferDevicePtr, real_bytes + complex_bytes);
      double* realDevicePtr;
      Complex* complexDevicePtr;
      if (fft_plan.m_dir == direction::R2C)
        {
          realDevicePtr = (double*) bufferDevicePtr;
          bufferDevicePtr += real_bytes;
          complexDevicePtr = (Complex*) bufferDevicePtr;
        }
      else if (fft_plan.m_dir == direction::C2R)
        {
          complexDevicePtr = (Complex*) bufferDevicePtr;
          bufferDevicePtr += complex_bytes;
          realDevicePtr = (double*) bufferDevicePtr;
        }

      fftx::array_t<3, double> // fftx::array_t<3, amrex::Real>
        realDeviceArray(fftx::global_ptr<double> // fftx::global_ptr<amrex::Real>
                        (realDevicePtr, 0, 1),
                        rDomain);

      fftx::array_t<3, Complex>
        complexDeviceArray(fftx::global_ptr<Complex>
                           (complexDevicePtr, 0, 1),
                           cDomain);

      if (fft_plan.m_dir == direction::R2C)
        {
          DEVICE_MEM_COPY(realDevicePtr,
                          fft_plan.m_real_array,
                          real_bytes,
                          MEM_COPY_HOST_TO_DEVICE);

          fftx::mdprdft<3>* tfm_3d_r2c = fft_plan.m_plan.m_tfm_3d_r2c;

          tfm_3d_r2c->transform(realDeviceArray, complexDeviceArray);

          DEVICE_MEM_COPY(fft_plan.m_complex_array,
                          complexDevicePtr,
                          complex_bytes,
                          MEM_COPY_DEVICE_TO_HOST);
        }
      else if (fft_plan.m_dir == direction::C2R)
        {
          DEVICE_MEM_COPY(complexDevicePtr,
                          fft_plan.m_complex_array,
                          complex_bytes,
                          MEM_COPY_HOST_TO_DEVICE);

          fftx::imdprdft<3>* tfm_3d_c2r = fft_plan.m_plan.m_tfm_3d_c2r;

          tfm_3d_c2r->transform(complexDeviceArray, realDeviceArray);

          DEVICE_MEM_COPY(fft_plan.m_real_array,
                          realDevicePtr,
                          real_bytes,
                          MEM_COPY_DEVICE_TO_HOST);
        }
      else
        {
          std::cout << "direction must be R2C or C2R" << std::endl;
          exit(-1);
        }
      DEVICE_FREE(bufferDevicePtr);
    }
}
