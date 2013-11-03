/**
@file

Copyright John Reid 2006-2010

*/

#include <boost/python.hpp>
#include "biopsy/python.h"
#include "biopsy/analyse.h"
#include "biopsy/convert_hit_to_bio.h"

#include <bio/svg_match.h>



using namespace boost;
using namespace boost::python;
using namespace std;


namespace biopsy {


void
build_analysis_svg(
    binding_hit::vec_ptr hits,
    binding_hit::vec_ptr max_chain,
    const sequence & seq,
    const BIO_NS::BuildSvgArgs & args
)
{
    //build the svg with all the arguments...
    namespace fs = boost::filesystem;
    biopsy::build_svg(
        fs::path( args.file ),
        args.title,
        seq,
        args.min_threshold,
        *hits,
        args.max_num_factors,
        args.show_labels,
        args.open_file,
        max_chain.get(),
        args.notes,
        args.max_threshold
    );
}


void export_analyse()
{
    using boost::python::arg;


    /**
    Analyse a sequence.
    */
    def(
        "append_sequence",
        append_sequence,
        ( arg( "sequence_vec" ), "sequence" ),
        "Appends a sequence to the sequence vector" );

    def(
        "append_sequences",
        append_sequences,
        ( arg( "sequence_vec" ), "sequence_vec" ),
        "Appends a sequence vector to the sequence vector" );


    def(
        "analyse",
        analyse,
        ( arg( "sequence" ), arg( "threshold" ) = BIOPSY_ANALYSE_THRESHOLD_DEFAULT ),
        "Analyses a sequence" );

    def(
        "analyse_phylo",
        analyse_phylo,
        (
            arg( "centre_sequence" ),
            arg( "phylo_sequences" ),
            arg( "threshold" ) = BIOPSY_ANALYSE_THRESHOLD_DEFAULT ),
        "Analyses a sequence including a phylogenetic comparison" );

    def(
        "analyse_max_chain",
        analyse_max_chain,
        (
            arg( "hit_array" ),
            arg( "max_box_limit" ) = 50000 ),
        "Calculates the maximal chain across the hit vectors." );

    def(
        "score_pssm_on_sequence",
        score_pssm_on_sequence,
        (
            arg( "pssm_name" ),
            arg( "sequence" ),
            arg( "threshold" ) = BIOPSY_ANALYSE_THRESHOLD_DEFAULT,
            arg( "result" ) ),
        "Scores a pssm on a sequence. Returns estimate that pssm binds there." );

    def(
        "score_pssms_on_sequence",
        score_pssms_on_sequence,
        (
            arg( "pssm_names" ),
            arg( "sequence" ),
            arg( "threshold" ) = BIOPSY_ANALYSE_THRESHOLD_DEFAULT ),
        "Scores a pssm on a sequence. Returns hit results." );

    def(
        "biobase_score_pssms_on_sequence",
        biobase_score_pssms_on_sequence,
        (
            arg( "pssm_names" ),
            arg( "sequence" ),
            arg( "threshold" ) = BIOPSY_BIOBASE_SCORE_THRESHOLD_DEFAULT ),
        "Biobase scores for a pssm on a sequence. Returns hit results." );

    def(
        "get_pathway_for_pssm",
        get_pathway_for_pssm,
        (
            arg( "pssm_name" )),
            "Determines the pathway associated with a PSSM." );

    def(
        "score_pssms_on_phylo_sequences",
        score_pssms_on_phylo_sequences,
        (
            arg( "pssm_names" ),
            arg( "sequences" ),
            arg( "threshold" ) = BIOPSY_ANALYSE_THRESHOLD_DEFAULT,
            arg( "phylo_threshold" ) = BIOPSY_ANALYSE_THRESHOLD_DEFAULT,
            arg( "calculate_maximal_chain" ) = true
        ),
        "Scores a pssm on sequences. Centre sequence is the first one. Returns: (hits, max_chain, unadjusted_hits)." );

    to_python_converter< phylo_sequences_result, tupleconverter< phylo_sequences_result > >();

    USING_BIO_NS;
    class_<
        BuildSvgArgs
    >(
        "BuildSvgArgs",
        "Parameters for building SVG file from old style BiFa analysis (pssm_match)",
        init<
            std::string,
            std::string,
            BIO_NS::float_t,
            BIO_NS::float_t,
            unsigned,
            bool,
            bool,
            std::string
        >(
            (
                arg( "file" ) = "matches.svg",
                arg( "title" ) = "BiFa analysis",
                arg( "max_threshold" ) = 0.00f,
                arg( "min_threshold" ) = 0.05f,
                arg( "max_num_factors" ) = 10,
                arg( "show_labels" ) = false,
#ifdef WIN32
                arg( "open_file" ) = true,
#else //WIN32
                arg( "open_file" ) = false,
#endif //WIN32
                arg( "notes" ) = ""
            ),
            "Initialise parameters for running a pssm_match"
        )
    )
        .def_readwrite(
            "file",
            &BuildSvgArgs::file,
            "Filename that SVG will be written to" )
        .def_readwrite(
            "title",
            &BuildSvgArgs::title,
            "Title of SVG" )
        .def_readwrite(
            "max_threshold",
            &BuildSvgArgs::max_threshold,
            "Maximum threshold to display on SVG" )
        .def_readwrite(
            "min_threshold",
            &BuildSvgArgs::min_threshold,
            "Minimum threshold to display on SVG" )
        .def_readwrite(
            "max_num_factors",
            &BuildSvgArgs::max_num_factors,
            "Maximum number of factors to display at top of SVG" )
        .def_readwrite(
            "show_labels",
            &BuildSvgArgs::show_labels,
            "Show labels next to hits" )
        .def_readwrite(
            "open_file",
            &BuildSvgArgs::open_file,
            "Open the SVG in browser" )
        .def_readwrite(
            "notes",
            &BuildSvgArgs::notes,
            "Notes to be added to SVG" )
    ;

    def(
        "build_analysis_svg",
        build_analysis_svg,
        (
            arg( "hits" ),
            arg( "max_chain" ),
            arg( "sequence" ),
            arg( "args" ) = BuildSvgArgs()
        ),
        "Build svg representation of the hits and maximal chain"
    );
}



} //namespace biopsy
