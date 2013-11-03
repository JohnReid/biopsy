/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/site_test_case.h"
#include "bio/biobase_db.h"
#include "bio/biobase_data_traits.h"
#include "bio/remo.h"
#include "bio/environment.h"
#include "bio/serialisable.h"

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable : 4312)
# pragma warning(disable : 4311)
#endif //_MSC_VER

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/filesystem/path.hpp>
using boost::multi_index_container;
using namespace boost::multi_index;
namespace fs = boost::filesystem;

#ifdef _MSC_VER
# pragma warning(pop)
#endif //_MSC_VER

#include <map>
using namespace std;

BIO_NS_START


namespace detail {

/* tags for accessing both sides of a bidirectional map */
struct from{};
struct to{};

/* The class template bidirectional_map wraps the specification
 * of a bidirectional map based on multi_index_container.
 */

template<typename FromType,typename ToType>
struct bidirectional_multimap
{
  typedef std::pair<FromType,ToType> value_type;

  /* A bidirectional multimap can be simulated as a multi_index_container
   * of pairs of (FromType,ToType) with two indices, one
   * for each member of the pair.
   */

  typedef multi_index_container<
    value_type,
    indexed_by<
      ordered_non_unique<
        tag<from>,member<value_type,FromType,&value_type::first> >,
      ordered_non_unique<
        tag<to>,  member<value_type,ToType,&value_type::second> >
    >
  > type;
};


typedef bidirectional_multimap< std::string, Site * >::type multimap_t;
typedef std::pair< multimap_t::const_iterator, multimap_t::const_iterator > multimap_const_range_t;

//visit a site
struct SiteVisitor
{
	multimap_t & map;
	unsigned num_visited;

	SiteVisitor(multimap_t & map)
		: map(map)
		, num_visited(0)
	{
	}

	void report() const
	{
		std::cout << "Visited " << num_visited << " sites\n";
		std::cout << "Found " << map.size() << " ensembl links\n";
	}

	void operator()(Site::map_t::value_type site)
	{
		++num_visited;

		//get the gene link from the description
		int gene_acc = site.second->get_transfac_gene_accession();
		if ( -1 == gene_acc )
		{
			return;
		}
		const TableLink gene_link( GENE_DATA, gene_acc );

		//can we find the gene?
		Gene * gene = BiobaseDb::singleton().get_entry< GENE_DATA >( gene_link );
		if ( 0 == gene )
		{
			return;
		}

		//look for ensembl links 
		for (DatabaseRefVec::const_iterator d = gene->database_refs.begin();
			gene->database_refs.end() != d;
			++d)
		{
			if ( ENSEMBL_DB == d->db )
			{
				map.insert( multimap_t::value_type( BIO_MAKE_STRING( *d ), site.second.get() ) );
			}
		}
	}
};


} //detail

void
SiteTestCases::init_singleton()
{
	deserialise_or_init< true >(
		*this,
		fs::path(
			BioEnvironment::singleton().get_site_test_cases_file()
		),
		boost::bind< void >(
			&SiteTestCases::init_from_biobase,
			_1
		) 
	);
}

void
SiteTestCases::init_from_biobase()
{
	using namespace detail;

	const bool verbose = true;
	const bool masked = true;
	multimap_t map;
	std::set< Site * > sites_found;
	std::set< std::string > genes_found;

	ReMoExtraction::ptr_t extraction;
	fs::path
		remo_extraction_archive(
			BioEnvironment::singleton().get_default_remo_archive_file()
		);
	cout << "Deserialising remo extraction from \"" << remo_extraction_archive._BOOST_FS_NATIVE() << "\"\n";
	extraction = ReMoExtraction::deserialise(remo_extraction_archive);

	std::for_each(
		BiobaseDb::singleton().get_sites().begin(),
		BiobaseDb::singleton().get_sites().end(),
		SiteVisitor( map ) );

	const unsigned num_groups = extraction->sequence_groups.size();
	cout << "Finding sites in " << num_groups << " sequence groups\n";
	cout << (masked ? "Using" : "Not using") << " masked remos\n";

	//for each sequence group
	for (ReMoSequenceGroup::list_t::const_iterator sg = extraction->sequence_groups.begin();
		extraction->sequence_groups.end() != sg;
		++sg)
	{
		try
		{
			//for each sequence in the group
			for (ReMoSequence::list_t::const_iterator s = sg->get()->sequences.begin();
				sg->get()->sequences.end() != s;
				++s)
			{
				//parse the sequence id
				ReMoSequenceId seq_id;
				if ( parse_sequence_id( s->get()->id, seq_id ) )
				{
					//look for the gene in our map
					multimap_const_range_t range = map.equal_range( seq_id.gene_id );
					if ( range.first != range.second )
					{
						//for each time we find the gene in the map
						while ( range.first != range.second )
						{
							//what sequence are we looking for?
							std::string site_sequence =	range.first->second->sequence;
							BIO_TO_UPPER( site_sequence );

							//and its reverse complement
							std::string site_sequence_rev_comp;
							reverse_complement( site_sequence, std::back_inserter( site_sequence_rev_comp ) );
							BIO_TO_UPPER( site_sequence_rev_comp );

							//for each bundle associated with this group
							for (ReMoBundle::map_t::const_iterator bundle = sg->get()->remo_bundles.begin();
								sg->get()->remo_bundles.end() != bundle;
								++bundle)
							{

								//the remo list for the sequence we are interested in
								ReMo::map_t::const_iterator rl = bundle->second->remos.find( s->get()->id );
								if ( bundle->second->remos.end() == rl )
								{
									//std::cout << "Could not find remo list for sequence: " << s->get()->id << "\n";
									continue;
								}
								const ReMo::list_t & remo_list = rl->second;

								//for each part of the remo
								for (ReMo::list_t::const_iterator remo_part = remo_list.begin();
									remo_list.end() != remo_part;
									++remo_part)
								{
									//the remo sequence
									const std::string & remo_sequence = remo_part->get()->get_sequence( masked );

									//for the site sequence and its reverse complement
									for (unsigned i = 0;
										2 != i;
										++i)
									{
										const bool rev_comp = i > 0;
										const std::string & seq = rev_comp ? site_sequence_rev_comp : site_sequence;

										//look for the sequence
										std::string::size_type index = remo_sequence.find( seq );
										if ( std::string::npos != index )
										{
											//we found it
											sites_found.insert( range.first->second );
											genes_found.insert( seq_id.gene_id );

											//generate a test case
											SiteTestCase::ptr_t test_case( new SiteTestCase );
											test_case->name =
												BIO_MAKE_STRING(
													(rev_comp ? "RC; " : "SS; ")
													<< range.first->second->get_name() << "; "
													<< seq << "; "
													<< bundle->first);
											test_case->site = range.first->second->get_link();

											//add each remo list in the bundle
											for (ReMo::map_t::const_iterator rl2 = bundle->second->remos.begin();
												bundle->second->remos.end() != rl2;
												++rl2)
											{
												seq_t sequence;

												//add its sequence to the test case
												ReMo::copy_sequence(
													rl2->second,
													std::back_inserter( sequence ),
													masked );

												//is it the centre sequence?
												if ( rl2 == rl )
												{
													test_case->centre_sequence = sequence;
												}
												else
												{
													test_case->sequences.push_back( sequence );
												}
											}
											push_back( test_case );

											//print the details if required
											if ( verbose )
											{
												std::cout << test_case->name << "\n";
											}
										}
									}
								}
							}

							++range.first;
						}
					}
				}
			}
		}
		catch (const std::exception & ex)
		{
			cerr << "Error: " << ex.what() << endl;
		}
		catch (const string & msg)
		{
			cerr << "Error: " << msg << endl;
		}
		catch (const char * msg)
		{
			cerr << "Error: " << msg << endl;
		}
		catch (...)
		{
			cerr << "Undefined error" << endl;
		}
	}

	cout
		<< "Found " << sites_found.size() << " different sites in "
		<< genes_found.size() << " different genes\n";
}


BIO_NS_END
