#pragma once

namespace DGtal
{
    namespace execution
    {
        struct Sequential {};
        
        template < typename TDomainSplitter >
        struct ParallelOpenMP
        {
            TDomainSplitter domainSplitter;
        };

    } // namespace execution

    // Sequential dispatch
    template < typename TAlgorithm,
               typename... TArgs >
    auto executionPolicyDispatch(
            execution::Sequential,
            TArgs... && theArgs )
    {
        return TAlgorithm::apply( std::forward<TArgs>(theArgs)... );
    }


    template < typename TAlgorithm,
               typename TDomainSplitter,
               typename... TArgs >
    auto executionPolicyDispatch(
            execution::ParallelOpenMP< TDomainSplitter > const& anOpenMP,
            TArgs... && theArgs

} // namespace DGtal
