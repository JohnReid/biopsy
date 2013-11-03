

#ifndef BIO_COUNTER_H_
#define BIO_COUNTER_H_

#include "bio/defs.h"

#include <boost/io/ios_state.hpp>

#include <iostream>


BIO_NS_START


/** Counts occurences of keys. Would rather have protected inheritance - boost serialization won't work with that though. */
template <class Key>
struct Counter : public std::map< Key, unsigned >
{
	typedef std::map< Key, unsigned > base_t;
	typedef Counter< Key > this_t;

	using base_t::clear;
	using base_t::begin;
	using base_t::end;
	using base_t::const_iterator;
	using base_t::find;
	using base_t::size;

	//typedef typename base_t::iterator iterator;
	//typedef typename base_t::const_iterator const_iterator;

	/** Sets count to 0 unless it already exists. */
	void initialise(const Key & key)
	{
		insert(value_type(key, 0));
	}

	/** Increment the count by one. */
	void increment(const Key & key)
	{
		typename base_t::iterator i = find(key);
		if (end() == i)
		{
			i = base_t::insert(typename base_t::value_type(key, 0)).first;
		}
		i->second++;
	}

	void operator()(const Key & key)
	{
		increment(key);
	}

	/** what is the total count? */
	unsigned get_total() const
	{
		unsigned result = 0;
		for (typename base_t::const_iterator i = begin();
			end() != i;
			++i)
		{
			result += i->second;
		}
		return result;
	}

	/** What is the maximum count? */
	unsigned get_max() const
	{
		typename base_t::const_iterator max = find_max();
		return end() == max ? 0 : max->second;
	}

	/** Which is the maximum count? */
	typename base_t::const_iterator find_max() const
	{
		typename base_t::const_iterator result = begin();
		for (typename base_t::const_iterator i = begin();
			end() != i;
			++i)
		{
			if (i->second > result->second)
			{
				result = i;
			}
		}
		return result;
	}

	/** How many for this Key. */
	unsigned get_count(Key key) const
	{
		typename base_t::const_iterator i = find(key);
		return end() == i ? 0 : i->second;
	}

	/** Add another set of counts to this one. */
	this_t &
	operator+=(const this_t & rhs)
	{
		for (typename base_t::const_iterator r = rhs.begin();
			rhs.end() != r;
			++r)
		{
			typename base_t::iterator t = find(r->first);
			if (end() != t)
			{
				t->second += r->second;
			}
			else
			{
				base_t::insert(*r);
			}
		}

		return *this;
	}

	void
	print(
		bool sorted_by_count = false,
		std::ostream & os = std::cout,
		const std::string & key_name = "Key",
		bool show_percentages = false,
		unsigned num_output = 0,
		unsigned histogram_width = 0,
		bool totals = false) const
	{
		boost::io::ios_base_all_saver ias( os );

		const unsigned total = get_total();
		const unsigned max = get_max();

		if (0 != histogram_width)
		{
			os << std::string(histogram_width + 1, ' ');
		}
		os << "       #, ";
		if (show_percentages)
		{
			os << "  %, ";
		}
		os << key_name << "\n";

		if (sorted_by_count)
		{
			count_set_t count_set;
			insert_into(count_set);
			unsigned i = 0;
			for (typename count_set_t::const_reverse_iterator c = count_set.rbegin();
				count_set.rend() != c;
				++c)
			{
				if (0 != histogram_width)
				{
					const unsigned bar_size = c->count * histogram_width / max;
					os
						<< std::string(bar_size, '*')
						<< std::string(histogram_width - bar_size + 1, ' ');
				}
				os << std::setw(8) << c->count << ", ";
				if (show_percentages)
				{
					os << std::setw(3) << 100 * c->count / total << ", ";
				}
				os << c->key << "\n";

				++i;
				if (num_output == i)
				{
					break;
				}
			}
		}
		else
		{
			unsigned i = 0;
			for (typename base_t::const_iterator c = begin();
				end() != c;
				++c, ++i)
			{
				if (0 != histogram_width)
				{
					const unsigned bar_size = c->second * histogram_width / max;
					os
						<< std::string(bar_size, '*')
						<< std::string(histogram_width - bar_size + 1, ' ');
				}
				os << std::setw(8) << c->second << ", ";
				if (show_percentages)
				{
					os << std::setw(3) << 100 * c->second / total << ", ";
				}
				os << c->first << "\n";

				++i;
				if (num_output == i)
				{
					break;
				}
			}
		}

		if (totals)
		{
			if (0 != histogram_width)
			{
				os << std::string(histogram_width + 1, ' ');
			}
			os << "   TOTAL, ";
			if (show_percentages)
			{
				os << "     ";
			}
			os << "# Entries" << "\n";

			if (0 != histogram_width)
			{
				os << std::string(histogram_width + 1, ' ');
			}
			os << std::setw(8) << total << ", ";
			if (show_percentages)
			{
				os << "     ";
			}
			os << size() << "\n";
		}
	}

	/** Holds one key and a count. */
	struct Count
	{
		Key key;
		unsigned count;

		Count(const Key & key, unsigned count)
			: key(key)
			, count(count)
		{
		}

		bool operator<(const Count & rhs) const
		{
			return
				count == rhs.count
					? key < rhs.key
					: count < rhs.count;
		}
	};
	typedef std::set< Count > count_set_t;

	void insert_into(count_set_t & count_set) const
	{
		for (typename base_t::const_iterator c = begin();
			end() != c;
			++c)
		{
			count_set.insert(Count(c->first, c->second));
		}
	}


protected:
	friend class boost::serialization::access;
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object< base_t >(*this);
    }
};


BIO_NS_END


#endif //BIO_COUNTER_H_

