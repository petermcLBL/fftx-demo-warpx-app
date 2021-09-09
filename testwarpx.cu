#include <stdio.h>

#include "device_macros.h"

// WrapFFTX.cpp includes AnyFFT.H
#include "WrapFFTX.cpp"

#include "fftx3utilities.h"
#include "sizes.h"

template<typename T>
T minSubarray(const T* arr, int lo, int hi)
{
  T val = arr[lo];
  for (int i = lo+1; i <= hi; i++)
    {
      if (arr[i] < val)
        {
          val = arr[i];
        }
    }
  return val;
}


template<typename T>
T maxSubarray(const T* arr, int lo, int hi)
{
  T val = arr[lo];
  for (int i = lo+1; i <= hi; i++)
    {
      if (arr[i] > val)
        {
          val = arr[i];
        }
    }
  return val;
}


template<typename T>
T avgSubarray(const T* arr, int lo, int hi)
{
  T tot = 0.;
  int len = 0;
  for (int i = lo; i <= hi; i++)
    {
      tot += arr[i];
      len++;
    }
  T avg = tot / (len * 1.);
  return avg;
}

void setRand(double& a_val)
{
  a_val = 1. - ((double) rand()) / (double) (RAND_MAX/2);
}

void setRand(std::complex<double>& a_val)
{
  double x, y;
  setRand(x);
  setRand(y);
  a_val = std::complex<double>(x, y);
}

double diffAbs(double a_x,
               double a_y)
{
  double diffNorm = a_x - a_y;
  if (diffNorm < 0.) diffNorm = -diffNorm;
  return diffNorm;
}

double diffAbs(std::complex<double>& a_x,
               std::complex<double>& a_y)
{
  double diffNorm = std::abs(a_x - a_y);
  return diffNorm;
}

DEVICE_FFT_RESULT deviceExecD2Z(DEVICE_FFT_HANDLE a_plan,
                                double* a_in,
                                std::complex<double>* a_out)
{
  return DEVICE_FFT_EXECD2Z(a_plan,
                            (DEVICE_FFT_DOUBLEREAL*) a_in,
                            (DEVICE_FFT_DOUBLECOMPLEX*) a_out);
}


template<typename T_IN, typename T_OUT>
struct deviceTransform
{
  deviceTransform(DEVICE_FFT_TYPE a_tp,
                  int a_dir = 0)
  {
    m_tp = a_tp;
    m_dir = a_dir;
  }
                  
  DEVICE_FFT_TYPE m_tp;

  int m_dir;

  DEVICE_FFT_RESULT exec(DEVICE_FFT_HANDLE a_plan,
                         T_IN* a_in,
                         T_OUT* a_out)
  {
    if (m_tp == DEVICE_FFT_Z2Z)
      {
        return DEVICE_FFT_EXECZ2Z(a_plan,
                                  (DEVICE_FFT_DOUBLECOMPLEX*) a_in,
                                  (DEVICE_FFT_DOUBLECOMPLEX*) a_out,
                                  m_dir);
      }
    else if (m_tp == DEVICE_FFT_D2Z)
      {
        return DEVICE_FFT_EXECD2Z(a_plan,
                                  (DEVICE_FFT_DOUBLEREAL*) a_in,
                                  (DEVICE_FFT_DOUBLECOMPLEX*) a_out);
      }
    else if (m_tp == DEVICE_FFT_Z2D)
      {
        return DEVICE_FFT_EXECZ2D(a_plan,
                                  (DEVICE_FFT_DOUBLECOMPLEX*) a_in,
                                  (DEVICE_FFT_DOUBLEREAL*) a_out);
      }
    else
      {
        return (DEVICE_FFT_RESULT) -1;
      }
  }
};
  

deviceTransform<double, std::complex<double> >
mdprdftDevice(DEVICE_FFT_D2Z);

deviceTransform<std::complex<double>, double>
imdprdftDevice(DEVICE_FFT_Z2D);


template<typename T_IN, typename T_OUT>
void compareSize(const fftx::point_t<3>& a_sz,
                 AnyFFT::direction a_dir,
                 deviceTransform<T_IN, T_OUT>& a_tfmDevice,
                 int a_verbosity)
{
  /*
    Allocate space for arrays, and set input array.
  */
  // const fftx::point_t<3> realSize = a_sz;
  // const fftx::point_t<3> complexSize = fftx::point_t<3>({{a_sz[0]/2 + 1, a_sz[1], a_sz[2]}});

  // This doesn't work. :/
  // const fftx::point_t<3> unit = fftx::point_t<3>::Unit();
  //   fftx::box_t<3> inputDomain(unit, inputSize);
  //   fftx::box_t<3> outputDomain(unit, outputSize);

  // N.B.: in transform, real size is {tfm_sz[0], tfm_sz[1], tfm_sz[2]}
  // and complex size is {tfm_sz[0]/2+1, tfm_sz[1], tfm_sz[2]},
  // where tfm_sz is the size of the real domain within the FFTX transform,
  // which has the reverse of the dimensions in CreatePlan.
  // Because of the reversals, here we define the domains with reversals.
  // The only time we use realDomain or complexDomain here other than to get
  // the number of points in them is when we look up the indices of where
  // there are differences.
  fftx::box_t<3> realDomain(fftx::point_t<3>({{1, 1, 1}}),
                            fftx::point_t<3>({{a_sz[2], a_sz[1], a_sz[0]}}));
  // The only reason 
  fftx::box_t<3> complexDomain(fftx::point_t<3>({{1, 1, 1}}),
                               fftx::point_t<3>({{a_sz[2]/2+1, a_sz[1], a_sz[0]}}));
  
  fftx::array_t<3, double> realArrayHost(realDomain);
  fftx::array_t<3, std::complex<double>> complexArrayHost(complexDomain);
  if (a_dir == AnyFFT::direction::R2C)
    {
      forall([](double(&v), const fftx::point_t<3>& p)
             {
               setRand(v);
             }, realArrayHost);
    }
  else if (a_dir == AnyFFT::direction::C2R)
    {
      forall([](std::complex<double>(&v), const fftx::point_t<3>& p)
             {
               setRand(v);
             }, complexArrayHost);
    }
  else
    {
      std::cout << "direction must be either R2C or C2R" << std::endl;
      exit(-1);
    }
    
  /*
    Define the FFTX transform.
  */
  std::cout << "Creating plan "
            << ((a_dir == AnyFFT::direction::R2C) ? "R2C" : "C2R")
            << " on size " << a_sz << std::endl;
  AnyFFT::FFTplan this_fftx_plan =
    AnyFFT::CreatePlan(// amrex::IntVect(AMREX_D_DECL(a_sz[0], a_sz[1], a_sz[2])),
                       a_sz,
                       realArrayHost.m_data.local(),
                       complexArrayHost.m_data.local(),
                       a_dir,
                       3);

  /*
    Execute the FFTX transform.
  */
  if (this_fftx_plan.m_plan.defined())
    {
      std::cout << "Executing plan" << std::endl;
      AnyFFT::Execute(this_fftx_plan);

      /*
        Define the device transform.
      */
      DEVICE_FFT_HANDLE this_device_plan;
      {
        std::cout << "Defining plan for device transform" << std::endl;
        auto result =
          DEVICE_FFT_PLAN3D(&this_device_plan, a_sz[2], a_sz[1], a_sz[0],
                            a_tfmDevice.m_tp);
        if (result != DEVICE_FFT_SUCCESS)
          {
            std::cout << "deviceFFT plan define failed\n" << std::endl;
            exit(-1);
          }
      }

      /*
        Execute the device transform, and compare result with that from FFTX.
      */
      fftx::box_t<3> inputDomain, outputDomain;
      T_IN* inputHostPtr;
      T_OUT* outputFFTXHostPtr;
      std::string tfmName;
      if (a_dir == AnyFFT::direction::R2C)
        {
          inputDomain = realDomain;
          outputDomain = complexDomain;
          inputHostPtr = (T_IN*) realArrayHost.m_data.local();
          outputFFTXHostPtr = (T_OUT*) complexArrayHost.m_data.local();
          tfmName = this_fftx_plan.m_plan.m_tfm_3d_r2c->name();
        }
      else if (a_dir == AnyFFT::direction::C2R)
        {
          inputDomain = complexDomain;
          outputDomain = realDomain;
          inputHostPtr = (T_IN*) complexArrayHost.m_data.local();
          outputFFTXHostPtr = (T_OUT*) realArrayHost.m_data.local();
          tfmName = this_fftx_plan.m_plan.m_tfm_3d_c2r->name();
        }
      else
        {
          std::cout << "direction must be either R2C or C2R" << std::endl;
          exit(-1);
        }
      
      {
        auto input_bytes = inputDomain.size() * sizeof(T_IN);
        auto output_bytes = outputDomain.size() * sizeof(T_OUT);
        
        char* bufferDevicePtr;
        DEVICE_MALLOC(&bufferDevicePtr, input_bytes + output_bytes);
        T_IN* inputDevicePtr = (T_IN*) bufferDevicePtr;
        bufferDevicePtr += input_bytes;
        T_OUT* outputDevicePtr = (T_OUT*) bufferDevicePtr;
        
        // Have already set realArrayHost.
        DEVICE_MEM_COPY(inputDevicePtr, // dest
                        inputHostPtr, // source
                        input_bytes, // bytes
                        MEM_COPY_HOST_TO_DEVICE); // type
        
        auto result = a_tfmDevice.exec(this_device_plan,
                                       inputDevicePtr,
                                       outputDevicePtr);
        if (result != DEVICE_FFT_SUCCESS)
          {
            std::cout << "deviceFFTExec launch failed\n" << std::endl;
            exit(-1);
          }
        
        auto nptsOutput = outputDomain.size();
        T_OUT* outputHostPtr = new T_OUT[nptsOutput];
        DEVICE_MEM_COPY(outputHostPtr, // dest
                        outputDevicePtr, // source
                        output_bytes, // bytes
                        MEM_COPY_DEVICE_TO_HOST); // type
        DEVICE_FREE(bufferDevicePtr);
        
        // Now find differences between outputHostPtr
        // and complexArrayHost.m_data.local().
        const double tol = 1.e-7;
        bool match = true;
        double maxDiff = 0.;
        for (size_t ind = 0; ind < nptsOutput; ind++)
          {
            T_OUT outputFFTXPoint = outputFFTXHostPtr[ind];
            T_OUT outputDeviceFFTPoint = outputHostPtr[ind];
            double diffAbsPoint = diffAbs(outputFFTXPoint, outputDeviceFFTPoint);
            updateMaxAbs(maxDiff, diffAbsPoint);
            bool matchPoint = (diffAbsPoint < tol);
            if (!matchPoint)
              {
                match = false;
                if (a_verbosity >= 3)
                  {
                    fftx::point_t<3> pt = pointFromPositionBox(ind, outputDomain);
                    // Take the flipped indices because remember we
                    // flipped the dimensions of the domains.
                    std::cout << "error at " << pt.flipped()
                              << ": FFTX " << outputFFTXPoint
                              << ", deviceFFT " << outputDeviceFFTPoint
                              << std::endl;
                  }
              }
          }
        delete[] outputHostPtr;
        if (match)
          {
            printf("YES, results match for %s. Max diff %11.5e\n",
                   tfmName.c_str(), maxDiff);
          }
        else
          {
            printf("NO, results do not match for %s. Max diff %11.5e\n",
                   tfmName.c_str(), maxDiff);
          }
      }
    }
  
  /*
    Destroy the FFTX transform.
  */
  std::cout << "Destroying plan" << std::endl;
  AnyFFT::DestroyPlan(this_fftx_plan);
}


int main(int argc, char* argv[])
{
  printf("Usage:  %s [verbosity=0]\n", argv[0]);
  // printf("verbosity 0 for avg times, 1 for min/max, 2 for all iterations, 3 for errors\n");
  int verbosity = 0;
  if (argc > 1)
    {
      verbosity = atoi(argv[1]);
    }
  printf("Running with verbosity %d\n", verbosity);

  // last entry is { 0, 0, 0 }
  int numentries = sizeof ( AllSizes3 ) / sizeof ( fftx::point_t<3> ) - 1;

  for ( int ind = 0; ind < numentries; ind++ )
    {
      fftx::point_t<3> sz = AllSizes3[ind];

      compareSize(sz, AnyFFT::direction::R2C, mdprdftDevice, verbosity);

      compareSize(sz, AnyFFT::direction::C2R, imdprdftDevice, verbosity);
    }
  
  printf("%s: All done, exiting\n", argv[0]);
  return 0;
}
