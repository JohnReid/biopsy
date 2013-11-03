#ifndef BIO_SEQUENCE_COLLECTION_H_
#define BIO_SEQUENCE_COLLECTION_H_

#include "bio/defs.h"
#include "bio/sequence.h"


BIO_NS_START

class SequenceCollection
{
public:
	virtual unsigned num_sequences() const = 0;
	virtual const seq_t & get_sequence(unsigned index) const = 0;
	virtual ~SequenceCollection() { };

	typedef boost::shared_ptr<SequenceCollection> ptr_t;
};

class SequenceCollectionList : public SequenceCollection
{
public:
	const SeqList & list;

	SequenceCollectionList(const SeqList & list) : list(list) { }

	virtual unsigned num_sequences() const { return list.size(); }
	virtual const seq_t & get_sequence(unsigned index) const
	{
		SeqList::const_iterator s = list.begin();
		for (unsigned i = 0; index != i; ++i)
		{
			++s;
		}
		return *s;
	}
};

class SequenceCollectionVector : public SequenceCollection
{
	typedef std::vector<seq_t> SeqVector;

	SeqVector sequences;

public:
	unsigned num_sequences() const { return sequences.size(); }

	const seq_t & get_sequence(unsigned index) const { return sequences[index]; }

	void add_sequence(const seq_t & sequence) { sequences.push_back(sequence); }
};




BIO_NS_END


#endif //BIO_SEQUENCE_COLLECTION_H_
