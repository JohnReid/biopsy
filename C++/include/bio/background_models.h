
#ifndef BIO_BACKGROUND_MODELS_H_
#define BIO_BACKGROUND_MODELS_H_

#include <bio/defs.h>
#include <bio/sequence.h>
#include <bio/sequence_collection.h>

#include <boost/shared_ptr.hpp>

BIO_NS_START


/** A model of DNA sequences. */
struct DnaModel
{
    typedef boost::shared_ptr<DnaModel> ptr_t;

    virtual ~DnaModel();

    /** Train on one sequence. */
    virtual void train(const seq_t & sequence) = 0;

    /** Return the per base likelihood under this model of the given sequence. */
    virtual prob_t get_likelihood(const seq_t & sequence) = 0;

    /** Train on multiple sequences concurrently. */
    virtual void train(const SequenceCollection & sequences) = 0;

    /** Return the per base likelihood under this model of the given sequences. */
    virtual prob_t get_likelihood(const SequenceCollection & sequences) = 0;

    /** Generate a random sequence from this model and append it to the sequence. */
    virtual void append_random_sequence(seq_t & seq, unsigned seq_length) const = 0;
};

void
gen_sequence_from_random_species_dna_hmm(seq_t & seq, size_t seq_length);


/**
Hands out references to DNA models.
*/
class DnaModelBroker
{
    //
    // Methods
    //
public:
    /**
    Gets a HMM sequence model with given number of states and order.
    */
    DnaModel::ptr_t get_hmm_model(unsigned num_states, unsigned order);

    /**
    Gets the sequence model used in nMica.
    */
    DnaModel::ptr_t get_nmica_model();

};


BIO_NS_END

#endif // BIO_BACKGROUND_MODELS_H_
