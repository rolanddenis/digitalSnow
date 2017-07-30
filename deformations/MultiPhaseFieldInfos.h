#pragma once

#include <vector>
#include <ostream>

template < typename TMultiImage >
struct MultiPhaseFieldInfos
{
  using Value = typename TMultiImage::Value;

  MultiImageInfos<TMultiImage> multiImageInfos; ///< Informations of the multi-image container.
  std::vector<Value> phasePerimeters;           ///< Estimated perimeter of each phases.
  Value morganCost;                             ///< Morgan's cost of the partition.
};

template < typename TMultiImage >
inline
std::ostream&
operator<< ( std::ostream & out, MultiPhaseFieldInfos<TMultiImage> const& object )
{
  out << "Multi-image infos:" << std::endl
      << object.multiImageInfos << std::endl;


  out << std::setprecision(5);

  // Phase perimeters.
  for ( std::size_t i = 0; i < object.phasePerimeters.size(); ++i )
    out << "P" << i << "=" << object.phasePerimeters[i] << " ";
  out << std::endl;

  // Morgan'cost
  out << "Morgan's cost = " << object.morganCost << std::endl;

  return out;
}

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */
