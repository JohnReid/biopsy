
#ifndef BIO_BIFA_H_
#define BIO_BIFA_H_

#include "bio/defs.h"
#include "bio/sequence.h"
#include "bio/binding_model.h"
#include "bio/pssm_match_args.h"


BIO_NS_START


/**
Input to a BiFa algorithm. A centre sequence and a list of sequences conserved in other species.
*/ 
struct BiFaInput
{
	seq_t centre_sequence;
	SeqList conserved_sequences;

	BiFaInput( const seq_t & centre_sequence = "" );

private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & centre_sequence;
        ar & conserved_sequences;
    }		 
};



/**
Output from a particular BiFa algorithm.
*/
struct BiFaOutput
{
	typedef BiFaOutput this_t;
	typedef boost::shared_ptr< this_t > ptr_t;
	typedef std::vector< ptr_t > vec_t;

	/** The hits. */
	BindingModel::hit_set_t hits;
};





/**
A BiFa algorithm.
*/
struct BiFaAlgorithm
{
	typedef BiFaAlgorithm this_t;
	typedef boost::shared_ptr< this_t > ptr_t;
	typedef std::vector< ptr_t > vec_t;

	virtual ~BiFaAlgorithm();

	/**
	Apply the algorithm to the input.
	*/
	virtual BiFaOutput::ptr_t operator()( const BiFaInput & input ) = 0;

	/**
	Get the model for the key (e.g. TableLink).
	*/
	virtual BindingModel * get_model_for( const boost::any & key) = 0;

	/**
	Get all the models applied to the input
	*/
	virtual void fill_model_universe( BindingModel::set_t & universe ) = 0;

	/**
	The name of the model.
	*/
	virtual std::string get_name() const = 0;

	static ptr_t get_default_algorithm();
};



struct MatchBiFaAlgorithm
	: BiFaAlgorithm
{
	/**
	Apply the algorithm to the input.
	*/
	virtual BiFaOutput::ptr_t operator()( const BiFaInput & input );

	/**
	Get the model for the key (e.g. TableLink).
	*/
	virtual BindingModel * get_model_for( const boost::any & key);

	/**
	Get all the models applied to the input
	*/
	virtual void fill_model_universe( BindingModel::set_t & universe );

	/**
	The name of the model.
	*/
	virtual std::string get_name() const;
};




/**
Transfac BiFa algorithm.
*/
struct TransfacBiFaAlgorithm
	: BiFaAlgorithm
{
	bool verbose;
	bool adjust_phylo;
	PssmMatchArgs args;

	TransfacBiFaAlgorithm(
		const PssmMatchArgs & args,
		bool adjust_phylo,
		bool verbose = false);

	virtual ~TransfacBiFaAlgorithm();

	/**
	Apply the algorithm to the input.
	*/
	virtual BiFaOutput::ptr_t operator()( const BiFaInput & input );

	/**
	Get the model for the key (e.g. TableLink).
	*/
	virtual BindingModel * get_model_for( const boost::any & key);

	/**
	Get all the models applied to the input
	*/
	virtual void fill_model_universe( BindingModel::set_t & universe );

	/**
	The name of the model.
	*/
	virtual std::string get_name() const;
};




BIO_NS_END

#endif //BIO_BIFA_H_

