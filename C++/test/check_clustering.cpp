


#include <bio/bayesian_hierarchical_clustering.h>
#include <bio/random.h>
USING_BIO_NS

#include <boost/assign.hpp>
#include <boost/assert.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/vector_property_map.hpp>
#include <boost/array.hpp>
using namespace boost;
using namespace boost::assign;
using boost::unit_test::test_suite;

#include <iostream>
using namespace std;

#ifndef _WIN32
# include <stdlib.h>
#endif //_WIN32

typedef boost::array<bool, 4> bool_data_t;
typedef std::vector<bool_data_t> bool_data_vec_t;
typedef bool_data_vec_t::const_iterator bool_data_it;
typedef std::set<bool_data_it> bool_data_set_t;
bool_data_vec_t bool_data;

typedef boost::array<double, 2> real_data_t;
typedef std::vector<real_data_t> real_data_vec_t;
typedef real_data_vec_t::const_iterator real_data_it;
typedef std::set<real_data_it> real_data_set_t;
real_data_vec_t real_data;

void
create_data()
{
	using namespace boost::assign;

	if (bool_data.empty())
	{
		bool_data +=
			list_of( true)( true)( true)( true ),
			list_of( true)(false)(false)(false ),
			list_of(false)( true)(false)(false ),
			list_of(false)(false)( true)(false ),
			list_of(false)(false)(false)( true ),
			list_of(false)( true)( true)( true ),
			list_of( true)(false)( true)( true ),
			list_of( true)( true)(false)( true ),
			list_of( true)( true)( true)(false ),
			list_of(false)(false)(false)(false );
	}

	if (real_data.empty())
	{
		real_data +=
			//list_of( 1.00)( 1.00),
			//list_of( 1.00)( 1.00),
			list_of( 1.00)( 1.00),
			//list_of(-1.00)( 1.00),
			//list_of(-1.00)( 1.00),
			list_of(-1.00)( 1.00),
			//list_of( 1.00)(-1.00),
			//list_of( 1.00)(-1.00),
			list_of( 1.00)(-1.00),
			//list_of(-1.00)(-1.00),
			//list_of(-1.00)(-1.00),
			list_of(-1.00)(-1.00),
			//list_of( 0.00)( 0.00),
			//list_of( 0.00)( 0.00),
			list_of( 0.00)( 0.00);
	}
}

ostream &
operator<<(ostream & os, const bool_data_t & bool_data)
{
	for (bool_data_t::const_iterator d = bool_data.begin();
		bool_data.end() != d;
		++d)
	{
		os << (*d ? "1" : "0");
	}
	return os;
}

ostream &
operator<<(ostream & os, const real_data_t & real_data)
{
	boost::io::ios_base_all_saver ias(os);
	//os.precision(3);
	os.unsetf(ios::adjustfield);
	os.unsetf(ios::floatfield);

	for (real_data_t::const_iterator d = real_data.begin();
		real_data.end() != d;
		++d)
	{
		os << *d << ":";
	}
	return os;
}

struct RandomScorer
{
	typedef double score_t;
	static score_t min_score() { return numeric_limits<score_t>::min(); }

	template <class Cluster>
	score_t
	operator()(
		const Cluster & cluster) const
	{
		return get_uniform_01();
	}

	template <class Cluster>
	score_t
	operator()(
		const Cluster & cluster_1,
		const Cluster & cluster_2) const
	{
		return get_uniform_01();
	}
};

template <class DataIt, class Scorer>
struct DfsVisitor : boost::default_dfs_visitor
{
	typedef typename ClusterTraits<DataIt, Scorer>::data_set_t data_set_t;
	data_set_t & data_set;

	DfsVisitor(data_set_t & data_set) : data_set(data_set) { }

	template < typename Vertex, typename Graph >
	void discover_vertex(Vertex u, const Graph & g)
	{
		if (g[u].first)
		{
			data_set.insert(g[u].second);
		}
	}
};


template <class Clusterer>
void
show_clusters(Clusterer & clusterer)
{
	typedef Clusterer clusterer_t;
	typedef typename clusterer_t::data_it data_it;
	typedef typename clusterer_t::scorer_t scorer_t;

	typename clusterer_t::tree_vertex_it v, end;
	for (tie(v, end) = vertices(clusterer.tree); end != v; ++v)
	{
		typedef DfsVisitor<data_it, scorer_t> dfs_visitor_t;
		typedef typename dfs_visitor_t::data_set_t data_set_t;
		data_set_t data_set;

		depth_first_visit(
			clusterer.tree,
			*v,
			dfs_visitor_t(data_set),
			vector_property_map<typename clusterer_t::tree_vertex>());

#ifdef VERBOSE_CHECKING
		for (typename data_set_t::const_iterator i = data_set.begin();
			data_set.end() != i;
			++i)
		{
			cout << **i << ",";
		}
		cout << endl;
#endif
	}
}

template <class ModelIt, class DataIt>
void
compare_models(
	ModelIt models_begin,
	ModelIt models_end,
	DataIt data_begin,
	DataIt data_end)
{
	boost::io::ios_base_all_saver ias(cout);
	cout.precision(9);
	cout.setf(ios::left, ios::adjustfield);
	cout.setf(ios::fixed, ios::floatfield);

	for (DataIt i = data_begin;
		data_end != i;
		++i)
	{
		typedef std::list<DataIt> data_it_list_t;
		data_it_list_t single_data_set = list_of(i);

		cout << *i << "\t";
		for (ModelIt m = models_begin;
			models_end != m;
			++m)
		{
			cout << "\t" << setw(12) << (*m)(single_data_set.begin(), single_data_set.end());
		}
		cout << endl;

		for (DataIt j = data_begin;
			data_end != j;
			++j)
		{
			data_it_list_t double_data_set = list_of(i)(j);

			cout << *i << "," << *j;
			for (ModelIt m = models_begin;
				models_end != m;
				++m)
			{
				cout << "\t" << setw(10) << (*m)(double_data_set.begin(), double_data_set.end());
			}
			cout << endl;

		}
	}
}

void
check_wishart_normal_prior()
{
	using namespace boost::assign;

	cout << "******* check_wishart_normal_prior()" << endl;

	create_data();

	typedef NormalWishartPrior model_t;
	typedef std::vector<model_t> model_vec_t;
	model_vec_t models(3);

	models[0] = model_t(2);
	models[0].s(0,0) = 1;
	models[0].s(1,0) = 0;
	models[0].s(0,1) = 0;
	models[0].s(1,1) = 1;
	models[0].m(0) = 0;
	models[0].m(1) = 0;
	models[0].r = 0.001;
	models[0].v = 2.0;

	models[1] = models[0];
	models[1].r = .1;

	models[2] = models[0];
	models[2].r = 1;

#ifdef VERBOSE_CHECKING
	compare_models(models.begin(), models.end(), real_data.begin(), real_data.end());
#endif
}

void
check_bernoulli_beta_prior()
{
	using namespace boost::assign;

	cout << "******* check_bernoulli_beta_prior()" << endl;

	create_data();

	typedef BernoulliBetaPrior model_t;
	typedef std::vector<model_t> model_vec_t;
	model_vec_t models(3);

	models[0] = model_t(4);
	models[0].alpha = list_of(1)(1)(1)(1);
	models[0].beta = list_of(1)(1)(1)(1);

	models[1] = models[0];
	models[1].alpha = list_of(1)(1)(1)(10);

	models[2] = models[0];
	models[2].beta = list_of(1)(1)(1)(10);

#ifdef VERBOSE_CHECKING
	compare_models(models.begin(), models.end(), bool_data.begin(), bool_data.end());
#endif
}

void
check_multinomial_dirichlet_prior()
{
	using namespace boost::assign;

	cout << "******* check_multinomial_dirichlet_prior()" << endl;

	create_data();

}

void
check_bayesian_hierarchical_clustering()
{
	using namespace boost::assign;

	cout << "******* check_bayesian_hierarchical_clustering()" << endl;

	create_data();

	if (false) //don't run for now
	{
#ifdef VERBOSE_CHECKING
		cout << "Random clustering" << endl;
#endif

		typedef RandomScorer scorer_t;
		scorer_t scorer;
		typedef TreeBuildingClusterer<bool_data_it, scorer_t> clusterer_t;
		clusterer_t clusterer;
		cluster<bool_data_it, scorer_t, clusterer_t>(
			bool_data.begin(),
			bool_data.end(),
			scorer,
			clusterer);

		show_clusters(clusterer);
	}

	if (false) //don't run for now
	{
#ifdef VERBOSE_CHECKING
		cout << "Bayesian clustering of bool data" << endl;
#endif

		typedef BernoulliBetaPrior model_t;
		typedef BayesianClusterScorer<BernoulliBetaPrior> scorer_t;
		typedef TreeBuildingClusterer<bool_data_it, scorer_t> clusterer_t;
		model_t model(4);
		model.alpha += 0.1,1.0,0.7,1.0;
		model.beta += 1.0,1.5,1.0,2.0;
		clusterer_t clusterer;
		scorer_t scorer(model);
		cluster<bool_data_it, scorer_t, clusterer_t>(
			bool_data.begin(),
			bool_data.end(),
			scorer,
			clusterer);

		show_clusters(clusterer);
	}

	{
#ifdef VERBOSE_CHECKING
		cout << "Bayesian clustering of real data" << endl;
#endif

		typedef NormalWishartPrior model_t;
		typedef BayesianClusterScorer<model_t> scorer_t;
		typedef TreeBuildingClusterer<real_data_it, scorer_t> clusterer_t;
		model_t model(2);
		model.s(0,0) = 1;
		model.s(1,0) = 0;
		model.s(0,1) = 0;
		model.s(1,1) = 1;
		model.m(0) = 0.0;
		model.m(1) = 0.0;
		model.r = 1;
		model.v = 3.0;
		clusterer_t clusterer;
		scorer_t scorer(model);
		cluster<real_data_it, scorer_t, clusterer_t>(
			real_data.begin(),
			real_data.end(),
			scorer,
			clusterer);

		show_clusters(clusterer);
	}
}


void register_clustering_tests(test_suite * test)
{
	test->add(BOOST_TEST_CASE(&check_bayesian_hierarchical_clustering), 0);
	test->add(BOOST_TEST_CASE(&check_multinomial_dirichlet_prior), 0);
	test->add(BOOST_TEST_CASE(&check_wishart_normal_prior), 0);
	test->add(BOOST_TEST_CASE(&check_bernoulli_beta_prior), 0);
}


