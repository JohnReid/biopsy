/**
@file

Copyright John Reid 2007, 2013
*/

#include "bio-pch.h"



#include <bio/hmm_dna.h>
#include <bio/species_file_sets.h>
#include <bio/options.h>
#include <bio/application.h>
#include <bio/environment.h>
#include <bio/serialisable.h>
USING_BIO_NS;

#include <boost/test/execution_monitor.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace boost;


#include <fstream>
using namespace std;


struct HmmTrainerApp : Application
{
	unsigned num_species_hmm_training_seqs;
	unsigned species_hmm_training_seq_length;
	unsigned max_order;
	unsigned max_num_states;
	bool increase_parameters;
	bool want_to_exit;
	//bool want_to_serialise; //now we serialise every iteration
	//size_t report_freq;

	HmmTrainerApp()
		: want_to_exit(false)
		//, want_to_serialise(false)
	{
		get_options().add_options()
			("max_order,o", po::value(&max_order)->default_value(5), "highest order")
			("max_num_states,s", po::value(&max_num_states)->default_value(3), "highest # states")
			("num_seqs,n", po::value(&num_species_hmm_training_seqs)->default_value(2000), "# sequences")
			("seq_length,l", po::value(&species_hmm_training_seq_length)->default_value(100), "sequence length")
			("increase_parameters", po::bool_switch(&increase_parameters)->default_value(false), "increase parameters every iteration")
			//("report_freq,r", po::value(&report_freq)->default_value(3), "how many iterations before printing likelihood")
			;
	}



	void init()
	{
		register_ctrl_handler();

		cout
			<< endl
			<< "Hit Ctrl-BREAK to save current state of HMMs" << endl
			<< "Hit Ctrl-C to save current state of HMMs and exit" << endl
			<< endl;
	}


	bool ctrl_handler(CtrlSignal signal)
	{
		//want_to_serialise = true;
		want_to_exit = (CTRL_BREAK_SIGNAL != signal);

		return true;
	}



	int task()
	{
		cout << "Training HMMs" << endl;

		//default hmm map
		DnaHmmOrderNumStateMap & hmm_map = DnaHmmOrderNumStateMap::singleton();

		for (unsigned num_states = 1; num_states <= max_num_states; ++num_states)
		{
			for (unsigned order = 0; order <= max_order; ++order)
			{
				if (! hmm_map.contains_model(num_states, order))
				{
					cout << "Inserting new model of order " << order << " and with " << num_states << " states\n";
					hmm_map.insert_model(num_states, order, create_dna_model(num_states, order));
				}
			}
		}

		//forever
		while (true)
		{
			//get the random sequences
			cout
				<< "Building " << num_species_hmm_training_seqs << " sequences "
				<< "each of length "  << species_hmm_training_seq_length << endl;

			SequenceCollection::ptr_t sequences = 
				get_random_sequence_collection(
					num_species_hmm_training_seqs,
					species_hmm_training_seq_length);

			//train the hmms
			cout << "Training\n";
			hmm_map.train_all(*sequences);

			//calculate the likelihood of the sequences under each hmm in the map
			for (DnaHmmOrderNumStateMap::model_map_t::const_iterator i = hmm_map.models.begin();
				i != hmm_map.models.end();
				++i)
			{
				cout << "(" << i->first.num_states << "," << i->first.order << "): "
					<< i->second->get_likelihood(*sequences)
					<< endl;
			}

			//serialise the map
			serialise< false >(
				hmm_map,
				fs::path(
					BioEnvironment::singleton().get_species_hmm_file().c_str()
				)
			);

			if (want_to_exit)
			{
				break;
			}

			//increase parameters
			if (increase_parameters)
			{
				num_species_hmm_training_seqs = num_species_hmm_training_seqs * 11 / 10;
				species_hmm_training_seq_length
					= size_t(species_hmm_training_seq_length + std::log((float_t) species_hmm_training_seq_length));
			}
		}

		return 0;
	}
};

int
main(int argc, char * argv [])
{
	return HmmTrainerApp().main(argc, argv);
}
