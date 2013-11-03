/**
@file

Copyright John Reid 2006

*/

#include "biopsy/python.h"
#include "biopsy/remo.h"



using namespace boost;
using namespace boost::python;
using namespace boost::python::indexing;
using namespace std;



namespace biopsy {



void export_remo()
{
	using namespace remo;

	class_<
		remome,
		noncopyable,
		remome::ptr
	> remome_class(
		"Remome",
		"A space of remos indexed by which sequences were aligned to find them",
		no_init )
		;


	//put everything else in the Remome class scope
	scope outer = remome_class;


	enum_< region >( "region" )
		.value( "upstream", region_upstream )
		.value( "downstream", region_downstream )
		.value( "gene", region_gene )
		.value( "undefined_region", region_undefined )
		;


	class_<
		location
	>(
		"Location",
		"A location w.r.t. a particular gene and region",
		boost::python::init< int, int >()
	)
		.def( "__repr__", &location::str )
		.def_readwrite( "start", &location::_start )
		.def_readwrite( "end", &location::_end )
		;


	class_<
		ensembl_id
	>(
		"EnsemblId",
		"An ensembl id",
		boost::python::init< std::string, unsigned >()
	)
		.def( "__eq__", &ensembl_id::operator== )
		.def( "__repr__", &ensembl_id::str )
		.def_readwrite( "prefix", &ensembl_id::_prefix )
		.def_readwrite( "num", &ensembl_id::_num )
		.def( "parse", &ensembl_id::parse )
		.staticmethod( "parse" )
		.def( "looks_like", &ensembl_id::looks_like )
		.staticmethod( "looks_like" )
		;
	class_<
		ensembl_id::list,
		ensembl_id::list_ptr
	>( "EnsemblIdList" )
		.def( container_suite< ensembl_id::list >() )
		;


	class_<
		exon
	>(
		"Exon",
		"An exon",
		boost::python::init< location, ensembl_id >()
	)
		.def_readwrite( "location", &exon::_location )
		.def_readwrite( "id", &exon::_id )
		;
	class_< exon::list >( "ExonList" )
		.def( container_suite< exon::list >() )
		;


	class_<
		ensembl_database_id
	>(
		"EnsemblDatabaseId",
		"Identifies an ensembl database",
		init< const species &, int, int, const std::string & >()
	)
		.def( "__eq__", &ensembl_database_id::operator== )
		.def( "__repr__", &ensembl_database_id::str )
		.def_readwrite( "species", &ensembl_database_id::_species )
		.def_readwrite( "software_version", &ensembl_database_id::_software_version )
		.def_readwrite( "ncbi_build", &ensembl_database_id::_ncbi_build )
		.def_readwrite( "build_version", &ensembl_database_id::_build_version )
		.def( "parse", &ensembl_database_id::parse )
		.staticmethod( "parse" )
		.def( "looks_like", &ensembl_database_id::looks_like )
		.staticmethod( "looks_like" )
		;


	class_<
		alignment_sequence_id,
		alignment_sequence_id::ptr,
		noncopyable
	>(
		"AlignmentSequenceId",
		"Identifies a sequence by gene id, transcript id, ensembl database and alignment version number",
		no_init
	)
		.def( "__repr__", &alignment_sequence_id::str )
		.def_readwrite( "gene_id", &alignment_sequence_id::_gene_id )
		.def_readwrite( "transcript_id", &alignment_sequence_id::_transcript_id )
		.def_readwrite( "db_id", &alignment_sequence_id::_db_id )
		.def_readwrite( "version", &alignment_sequence_id::_version )
		.def( "parse", &alignment_sequence_id::parse )
		.staticmethod( "parse" )
		.def( "looks_like", &alignment_sequence_id::looks_like )
		.staticmethod( "looks_like" )
		;
	class_<
		alignment_sequence_id::list,
		alignment_sequence_id::list_ptr,
		noncopyable
	>( "AlignmentSequenceIdList" )
		.def( container_suite< alignment_sequence_id::list >() )
		;


	class_<
		alignment_sequence_info,
		alignment_sequence_info::ptr,
		noncopyable
	>(
		"AlignmentSequenceInfo",
		"Holds information about a sequence, such as position, the exons, ...",
		no_init )
		.def_readwrite( "length", &alignment_sequence_info::_length )
		.def_readwrite( "has_position", &alignment_sequence_info::_has_position )
		.def_readwrite( "position", &alignment_sequence_info::_position )
		.def_readwrite( "exons", &alignment_sequence_info::_exons )
		.def_readwrite( "region", &alignment_sequence_info::_region )
		;


	class_<
		remo_sequence,
		remo_sequence::ptr,
		noncopyable
	>(
		"RemoSequence",
		"A sequence in a remo for one species",
		no_init )
		.def_readwrite( "masked_sequence", &remo_sequence::_masked_sequence )
		.def_readwrite( "unmasked_sequence", &remo_sequence::_unmasked_sequence )
		.def_readwrite( "location", &remo_sequence::_location )
		.def_readwrite( "target_location", &remo_sequence::_target_location )
		.def_readwrite( "conservation", &remo_sequence::_conservation )
		.def_readwrite( "repeat_ratio", &remo_sequence::_repeat_ratio )
		.def_readwrite( "belief", &remo_sequence::_belief )
		;
	class_<
		remo_sequence::list,
		remo_sequence::list_ptr,
		noncopyable
	>( "RemoSequenceList" )
		.def( container_suite< remo_sequence::list >() )
		;


	class_<
		module,
		module::ptr,
		noncopyable
	>(
		"Remo",
		"A regulatory module",
		no_init )
		.def( "get_sequence_ids", &module::get_sequence_ids )
		.def( "get_sequences", &module::get_sequences )
		.def( "get_sequence_for", &module::get_sequence_for )
		;
	class_<
		module::list,
		module::list_ptr,
		noncopyable
	>( "RemoList" )
		.def( container_suite< module::list >() )
		;


	class_<
		aligned_sequence_set,
		aligned_sequence_set::ptr,
		noncopyable
	>(
		"AlignedSequenceSet",
		"Defines the sequences that have been aligned to search for remos",
		no_init )
		.def_readwrite( "centre_sequence", &aligned_sequence_set::_centre_sequence )
		.def( "get_sequence_ids", &aligned_sequence_set::get_sequence_ids )
		.def( "get_sequence_info", &aligned_sequence_set::get_sequence_info )
		.def(
			"get_sequences_for_remo",
			get_sequences_for_remo,
			"Get the list of sequences for a remo. The centre sequence is first." )
		.def(
			"get_remo_id",
			get_remo_id,
			"Get the id for a remo" )
		;
	class_<
		aligned_sequence_set::list,
		noncopyable,
		aligned_sequence_set::list_ptr
	>( "AlignedSequencesList" )
		.def( container_suite< aligned_sequence_set::list >() )
		;


	class_<
		ensembl_id_alignment_map,
		ensembl_id_alignment_map::ptr,
		noncopyable
	>(
		"EnsemblIdAlignmentMap",
		"Maps Ensembl ids to lists of aligned sequences containing them",
		no_init )
		.def( "get_genes", &ensembl_id_alignment_map::get_genes )
		.def( "get_alignments_for", &ensembl_id_alignment_map::get_alignments_for )
		;


	remome_class
		.def( "get_aligned_sequences", &remome::get_aligned_sequences )
		.def( "get_remos_for", &remome::get_remos_for )
		.def( "remove_remos_for", &remome::remove_remos_for )
		.def( "serialise", &remome::serialise )
		.def( "deserialise", &remome::deserialise )
		.staticmethod( "deserialise" )
		.def( "load", &load_remome_from_file )
		.staticmethod( "load" )
		.def( "parse", &parse_remome_from_file )
		.staticmethod( "parse" )
		.def(
			"make_gene_alignment_map",
			&make_gene_alignment_map,
			"Takes a list of aligned sequences and creates a map from gene ids to them" )
		.staticmethod( "make_gene_alignment_map" )
		.def(
			"get_remo_from_id",
			get_remo_from_id,
			"Get the remo from the id" )
		;

	to_python_converter< remo_locator, tupleconverter< remo_locator > >();
}



} //namespace biopsy
