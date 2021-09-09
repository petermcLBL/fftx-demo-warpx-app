#ifndef fftx_plan_PRECOMPILE_H
#define fftx_plan_PRECOMPILE_H

#include "fftx_mdprdft_public.h"
#include "fftx_imdprdft_public.h"

#include "fftx3.hpp"
#include "mdprdft.fftx.precompile.hpp"
#include "imdprdft.fftx.precompile.hpp"

/*
  Kludge that stands in for a vendor FFT plan on the model of
  fftw_plan, cufftHandle, rocfft_plan.

  The only transforms used by WarpX are R2C and C2R,
  but they may be 2D or 3D, single- or double-precision.
  Here I have only 3D double-precision.
 */
namespace fftx {
  
  // template <int DIM, typename T_IN, typename T_OUT>
  class fftx_plan
  {
  public:
    fftx_plan()
    {
      m_tfm_3d_r2c = NULL;
      m_tfm_3d_c2r = NULL;
    }

    bool defined()
    {
      bool got_r2c = (m_tfm_3d_r2c != NULL);
      bool got_c2r = (m_tfm_3d_c2r != NULL);
      if (got_r2c && !got_c2r)
        {
          return m_tfm_3d_r2c->defined();
        }
      else if (!got_r2c && got_c2r)
        {
          return m_tfm_3d_c2r->defined();
        }
      else
        { // neither r2c nor c2r is set, or both are set
          return false;
        }
    }

    
    ~fftx_plan()
    {
      if (m_tfm_3d_r2c != NULL)
        {
          delete m_tfm_3d_r2c;
        }
      if (m_tfm_3d_c2r != NULL)
        {
          delete m_tfm_3d_c2r;
        }
    }

    // List of all possible template instantiations.
    mdprdft<3>* m_tfm_3d_r2c;

    imdprdft<3>* m_tfm_3d_c2r;
  };
};

#endif  
