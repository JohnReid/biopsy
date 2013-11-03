
#ifndef BIO_BIFA_TEST_CASE_H_
#define BIO_BIFA_TEST_CASE_H_

#include <bio/defs.h>
#include <bio/bifa_algorithm.h>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/operators.hpp>

#include <string>
#include <vector>


BIO_NS_START


struct ROCPoint
{
    double specificity;
    double sensitivity;

    ROCPoint( double specificity, double sensitivity );
};
std::ostream &
operator<<( std::ostream & os, const ROCPoint & roc_point );



struct BinaryTestResults
    : boost::addable< BinaryTestResults >
{
    unsigned true_positives;
    unsigned false_positives;
    unsigned true_negatives;
    unsigned false_negatives;

    BinaryTestResults(
        unsigned true_positives = 0,
        unsigned false_positives = 0,
        unsigned true_negatives = 0,
        unsigned false_negatives = 0);

    BinaryTestResults & operator()(bool tested_positive, bool should_be_positive);

    ROCPoint get_roc_point() const;

    BinaryTestResults & operator+=( const BinaryTestResults & rhs );
};



/**
Contains the input for a test case and a list of the factors expected to bind there.
*/ 
struct BiFaTestCase
{
    typedef BiFaTestCase this_t;
    typedef boost::shared_ptr< this_t > ptr_t;
    typedef std::vector< ptr_t > vec_t;
    typedef std::map< double, ROCPoint > threshold_result_map_t; /**< Maps thresholds to test results. */

    std::string name;
    BiFaInput input;
    std::vector< BindingModel * > binders;

    template< typename BinderIt >
    ROCPoint get_roc_point( double threshold, const BiFaOutput & test_output, BinderIt binders_begin, BinderIt binders_end ) const;

    template< typename BinderIt >
    BinaryTestResults get_results( double threshold, const BiFaOutput & test_output, BinderIt binders_begin, BinderIt binders_end ) const;

    /** Examines the result map and calculates next threshold to try. */
    static double calculate_next_threshold(
        threshold_result_map_t::const_iterator begin,
        threshold_result_map_t::const_iterator end,
        double min_threshold = 0.0,
        double max_threshold= 1.0);

};


#if 0
template< typename BinderIt >
BinaryTestResults
BiFaTestCase::get_results( double threshold, const BiFaOutput & test_output, BinderIt binders_begin, BinderIt binders_end ) const
{
    BinaryTestResults results;

    for ( ; binders_end != binders_begin; ++binders_begin)
    {
        //an iterator to the binder score 
        DnaBinder::score_map_t::type::const_iterator binder_score = test_output.binding_probs.find( *binders_begin );

        //is the score greater than the threshold?
        const bool above_threshold = test_output.binding_probs.end() != binder_score && binder_score->score  > threshold;

        //should we find this in the results
        const bool in_test_case = binders.end() != std::find( binders.begin(), binders.end(), *binders_begin );

        //increment statistics
        results( above_threshold, in_test_case );
    }

    return results;
}
#endif


template< typename BinderIt >
ROCPoint
BiFaTestCase::get_roc_point( double threshold, const BiFaOutput & test_output, BinderIt binders_begin, BinderIt binders_end ) const
{
    return make_results( get_statistics( threshold, test_output, binders_begin, binders_end ) );
}




BIO_NS_END

#endif //BIO_BIFA_TEST_CASE_H_

