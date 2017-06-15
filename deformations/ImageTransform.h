#pragma once

namespace DGtal
{
    template < typename TDomain, typename TValue > class ImageContainerBySTLVector;
} // namespace DGtal

/// Generic image transformation implementation
template < typename TInputImage,
           typename TOutputImage,
           typename TFunctor >
void ImageTransform( TInputImage const& anInputImage, TOutputImage & anOutputImage, TFunctor aFunctor )
{
    ASSERT( anInputImage.domain().lowerBound() == anOutputImage.domain().lowerBound()
            && anInputImage.domain().upperBound() == anOutputImage.domain().upperBound() );

    // TODO: what if points have different types ??
    for ( auto const& pt : anInputImage.domain() )
        anOutputImage.setValue( pt, aFunctor( pt, anInputImage(pt) ) );
}

/// Image transformation implementation for ImageContainerBySTLVector.
template < typename TInputDomain, typename TInputValue,
           typename TOutputDomain, typename TOutputValue,
           typename TFunctor >
void ImageTransform(
    DGtal::ImageContainerBySTLVector< TInputDomain,  TInputValue > const& anInputImage,
    DGtal::ImageContainerBySTLVector< TOutputDomain, TOutputValue > & anOutputImage,
    TFunctor aFunctor )
{
    ASSERT( anInputImage.domain().lowerBound() == anOutputImage.domain().lowerBound()
            && anInputImage.domain().upperBound() == anOutputImage.domain().upperBound() );

    auto input_it            = anInputImage.begin();
    auto pt_it               = anInputImage.domain().begin();

    // TODO: should we detect lambda parameters or do we trust in compiler optimizations ?
    for ( auto & v : anOutputImage )
        v = aFunctor( *pt_it++, *input_it++ );
}
