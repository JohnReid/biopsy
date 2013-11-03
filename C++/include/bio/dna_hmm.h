#ifndef BIO_DNA_HMM_H_
#define BIO_DNA_HMM_H_

#include <bio/defs.h>
#include <bio/sequence.h>
#include <bio/hidden_markov_model.h>
#include <bio/hmm_gen_sequence.h>
#include <bio/chromosomes_file_set.h>

#include <boost/filesystem/path.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/access.hpp>

#include <string>
#include <map>



BIO_NS_START

typedef MarkovState<NucleoCode> DnaMarkovState;
typedef HiddenMarkovModel<NucleoCode> DnaHiddenMarkovModel;
typedef boost::shared_ptr<DnaHiddenMarkovModel> DnaHiddenMarkovModelPtr;
typedef HmmSequenceGenerator<NucleoCode> DnaSequenceGenerator;


/** Maps names to HMMs. */
struct HmmMap
{
private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & map;
    }

public:
    typedef std::map<std::string, DnaHiddenMarkovModelPtr> map_t;
    map_t map;
};

typedef std::map<std::string, ChromosomesFileSet::ptr_t> species_file_sets_t;
extern species_file_sets_t species_file_sets;

extern HmmMap species_hmms; //Hidden markov models trained on the genomes of named species
extern size_t species_hmm_training_seq_length; //The size of the sequences we use to train the species hmms
extern size_t num_species_hmm_training_seqs; //The # of sequences we use to train the species hmms

/** Make sure we have built the file sets of chromosomes for each species. */
void
ensure_species_file_sets_built();


/** Deserialise the species_hmms from disk or build from scratch. */
void
deserialise_or_build_species_hmms();


/** Build the species_hmms from scratch. */
void
build_species_hmms();


/** Deserialise the species_hmms from disk. */
void
deserialise_species_hmms();


/** Serialise the species_hmms to disk. */
void
serialise_species_hmms();


/** Continuously train the species hmms. */
prob_t
train_species_hmms(size_t training_seq_size);


/** Train the species hmms on multiple sequences of the length=species_hmm_training_size. */
prob_t
train_species_hmms_multiple(size_t num_sequences, size_t seq_size);


/** Trains the given hmm on one good sequence of the given length from the file set. */
prob_t
train_hmm_on_file_set(
    DnaHiddenMarkovModel & hmm,
    const ChromosomesFileSet & file_set,
    size_t training_seq_size);

/** Trains the given hmm multiply on given sequences. */
prob_t
train_hmm_multiply_on_files(
    DnaHiddenMarkovModel & hmm,
    FileSeqVec & seqs);

/** Trains the given hmm multiply on several sequences of the given length from the file set. */
prob_t
train_hmm_multiply_on_file_set(
    DnaHiddenMarkovModel & hmm,
    const ChromosomesFileSet & file_set,
    size_t num_seqs,
    size_t seq_length);

void
gen_sequence_from_random_species_dna_hmm(seq_t & seq, size_t seq_length);

BIO_NS_END

#endif //#ifndef BIO_DNA_HMM_H_
