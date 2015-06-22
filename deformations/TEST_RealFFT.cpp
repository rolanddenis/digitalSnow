#include <cmath>
#include <algorithm>

#include <DGtal/base/Common.h>
#include <DGtal/kernel/domains/HyperRectDomain.h>
#include <DGtal/kernel/SpaceND.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include "RealFFT.h" 
#include "VTKWriter.h"

int main()
{
  constexpr typename DGtal::Dimension N = 2;
  using real = double;
  using integer = int;
  using Space = DGtal::SpaceND<N, integer>;
  using Point = typename Space::Point;
  using Domain = DGtal::HyperRectDomain<Space>;
  using Image = DGtal::ImageContainerBySTLVector<Domain, real>;
  using FFT = DGtal::RealFFT<Domain, real>;

  real dT = 5; // Diffusion coefficient

  // Image initialization
  Domain domain{ Point::diagonal(-64), Point::diagonal(63) };
  Image image{ domain };
  for ( auto const& pt : domain )
    image.setValue( pt, pt.norm1() <= 30 ? 1 : 0 );

  FFT fft(image.domain());

  // Copy data
  auto spatial_image = fft.getSpatialImage();
  std::copy( image.cbegin(), image.cend(), spatial_image.begin() );

  // Forward transformation
  fft.forwardFFT(FFTW_ESTIMATE);

  // Convolution
  auto const extent = fft.getSpatialExtent();
  auto freq_image = fft.getFreqImage();
  for ( auto it = freq_image.begin(); it != freq_image.end(); ++it )
    {
      auto const& point = it.getPoint();
      
      double norm2 = 0;
      for ( size_t j = 0; j < Image::dimension; ++j)
        {
          double coord = static_cast<double>(point[j]) / extent[j];
          if ( coord >= 0.5 ) coord -= 1.;
          norm2 += coord*coord;
        }
      const double c = std::exp( -4*M_PI*M_PI*dT*norm2 );

      // New value
      auto const v = *it;
      *it = { c*std::real(v), c*std::imag(v) };
    }

  // Back in spatial space
  fft.backwardFFT(FFTW_ESTIMATE);

  // Store the result
  const size_t n = fft.getSpatialDomain().size();
  std::transform(
      spatial_image.cbegin(),
      spatial_image.cend(),
      image.begin(),
      [n] (real x) { return x/n; }
  );

  // Export
  if ( N == 2 )
    {
      DGtal::VTKWriter<Domain>( "fft_test", image.domain() ) << "data" << image;
    }

  return 0;
}

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */
