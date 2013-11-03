
#ifndef BIO_BAYESIAN_HIERARCHICAL_CLUSTERING_H_
#define BIO_BAYESIAN_HIERARCHICAL_CLUSTERING_H_

#include "bio/defs.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <set>
#include <list>
#include <limits>
#include <iterator>
#include <algorithm>
#include <iostream>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_math.h>

BIO_NS_START

/** @file
See http://www.gatsby.ucl.ac.uk/~heller/bhc.pdf for the motivation for these algorithms.
*/

/** Aggregates a data set with its associated score. */
template <class DataIt, class Score>
struct Cluster
{
    typedef std::set<DataIt> data_set_t;

    data_set_t data_set;
    Score score;
};


/** A Normal-Wishart prior.

Each datum should be an iterator that points to a sequence of reals.

I am not convinced this is implemented correctly...
*/
struct NormalWishartPrior
{
    typedef double real_t;
    typedef boost::numeric::ublas::vector<real_t> vector_t;
    typedef boost::numeric::ublas::matrix<real_t> matrix_t;
    matrix_t s;
    vector_t m;
    real_t r;
    real_t v;
    unsigned get_num_dims() const { return s.size1(); }

    NormalWishartPrior(unsigned num_dimensions = 0)
        : s(num_dimensions,  num_dimensions)
        , m(num_dimensions)
    {
    }

    template <class DataIt>
    real_t
    operator()(
        DataIt data_begin,
        DataIt data_end) const
    {
        using namespace boost::numeric::ublas;

        //if we have no data
        if (data_end == data_begin)
        {
            return real_t(1);
        }

        //calculate the number of data points
        unsigned n = 0; //std::count_if(data_begin, data_end, true);
        for (DataIt di = data_begin;
            data_end != di;
            ++di)
        {
            ++n;
        }

        //put the data set in a matrix
        matrix_t x(get_num_dims(), n);
        unsigned i = 0;
        for (DataIt di = data_begin;
            data_end != di;
            ++di, ++i)
        {
            for (unsigned j = 0;
                get_num_dims() != j;
                ++j)
            {
                x(j,i) = (**di)[j];
            }
        }

        //calculate the sum of the data in each dimension
        vector_t x_sum(get_num_dims());
        for (unsigned j = 0;
            get_num_dims() != j;
            ++j)
        {
            x_sum(j) = 0;
            for (unsigned i = 0;
                n != i;
                ++i)
            {
                x_sum(j) += x(j, i);
            }
        }

        //calculate s_dash - the first term is s
        matrix_t s_dash(s);
        {
            //add X.XT to s_dash
            s_dash += prod(x, trans(x));

            //add m.mT term to s_dash
            s_dash += (r * n / (n + r)) * outer_prod(m, m);

            //add x_sum.x_sumT term to s_dash
            s_dash += (1.0 / (n + r)) * outer_prod(x_sum, x_sum);

            //subtract m.x_sum term from s_dash
            s_dash -= (r / (n + r)) * (outer_prod(m, x_sum) + outer_prod(x_sum, m));
        }

        //std::cout << s << std::endl;
        //std::cout << s_dash << std::endl;
        //std::cout << prod(x, trans(x)) << std::endl;

        //calculate the determinants
        const real_t s_det = get_determinant(s);
        const real_t s_dash_det = get_determinant(s_dash);

        //v_dash
        const real_t v_dash = v + n;

        //the number of dimensions
        const real_t k = get_num_dims();

        //the number of data points
        const real_t N = n;

        //the natural logarithm of the result
        real_t ln_p_D_given_H1 = 0.0;

        ln_p_D_given_H1 -= N * k / 2 * gsl_sf_log(2 * M_PI);
        ln_p_D_given_H1 += k / 2 * gsl_sf_log(r / (N + r));
        ln_p_D_given_H1 += v / 2 * gsl_sf_log(fabs(s_det));
        ln_p_D_given_H1 -= v_dash / 2 * gsl_sf_log(fabs(s_dash_det));
        ln_p_D_given_H1 += v_dash * k / 2 * gsl_sf_log(2);
        ln_p_D_given_H1 -= v * k / 2 * gsl_sf_log(2);

        //add/subtract the pi terms
        for (unsigned d = 0; k != d; ++d)
        {
            ln_p_D_given_H1 +=
                gsl_sf_lngamma((v_dash + 1 - d) / 2)
                - gsl_sf_lngamma((v + 1 - d) / 2);

        }

        return gsl_sf_exp(ln_p_D_given_H1);
    }

    template <class M>
    static real_t
    get_determinant(const M & m)
    {
        BOOST_ASSERT(m.size1() == m.size2());

        matrix_t lu(m);
        //permutation_matrix<unsigned> p(m.size());

        lu_factorize(lu);

        real_t result = 0;
        for (unsigned j = 0;
            lu.size1() != j;
            ++j)
        {
            result += lu(j,j);
        }
        return result;
    }
};




/** A Bernoulli - Beta prior.

Each datum should be an iterator that points to a sequence of bools (or convertible types).
*/
struct BernoulliBetaPrior
{
    typedef double real_t;
    typedef std::vector<real_t> param_vec_t;
    param_vec_t alpha;
    param_vec_t beta;

    BernoulliBetaPrior(unsigned num_dimensions = 0)
        : alpha(num_dimensions)
        , beta(num_dimensions)
    {
    }

    template <class DataIt>
    real_t
    operator()(
        DataIt data_begin,
        DataIt data_end)
    {
        if (data_end == data_begin)
        {
            return real_t(1);
        }

        assert(alpha.size() == beta.size());
        assert((*data_begin)->size() == alpha.size());

        //the log of p(D|H1)
        double ln_p_D_H1 = 0.0;

        unsigned n = 0;
        for (DataIt i = data_begin; data_end != i; ++i)
        {
            ++n;
        }

        //for each dimension
        for (size_t d = 0; alpha.size() != d; ++d)
        {
            //work out m_d
            unsigned m_d = 0;
            for (DataIt i = data_begin; data_end != i; ++i)
            {
                assert((*i)->size() == alpha.size());

                if ((**i)[d])
                {
                    ++m_d;
                }
            }

            ln_p_D_H1 +=
                gsl_sf_lngamma(alpha[d] + beta[d])
                + gsl_sf_lngamma(alpha[d] + m_d)
                + gsl_sf_lngamma(beta[d] + n - m_d)
                - gsl_sf_lngamma(alpha[d])
                - gsl_sf_lngamma(beta[d])
                - gsl_sf_lngamma(alpha[d] + beta[d] + n);
        }

        return gsl_sf_exp(ln_p_D_H1);
    }
};




/** A Multinomial Dirichlet prior.

Each datum should be an iterator that points to a sequence of bools (or convertible types).
*/
struct MultinomialDirichletPrior
{
    typedef double real_t;
    typedef std::vector<real_t> param_vec_t;
    param_vec_t alpha;

    MultinomialDirichletPrior(unsigned num_dimensions = 0)
        : alpha(num_dimensions)
    {
    }

    template <class DataIt>
    real_t
    operator()(
        DataIt data_begin,
        DataIt data_end)
    {
        if (data_end == data_begin)
        {
            return real_t(1);
        }

        assert((*data_begin)->size() == alpha.size());

        //the log of p(D|H1) ie log of the result
        double ln_p_D_H1 = 0.0;

        //TODO - not implemented

        return gsl_sf_exp(ln_p_D_H1);
    }
};






/** The score of a cluster based on a prior model and the score of the two clusters
it is composed of. */
struct BayesianClusterScore
{
    typedef double real_t;

    real_t p_Dk_given_H1k;
    real_t p_Dk_given_Tk;

    BayesianClusterScore(real_t p_Dk_given_H1k = 1, real_t p_Dk_given_Tk = 1)
        : p_Dk_given_H1k(p_Dk_given_H1k)
        , p_Dk_given_Tk(p_Dk_given_Tk)
    {
    }

    bool operator<(const BayesianClusterScore & rhs) const
    {
        if (real_t(0) == p_Dk_given_Tk)
        {
            return false;
        }
        if (real_t(0) == rhs.p_Dk_given_Tk)
        {
            return false;
        }
        return operator()() < rhs();
    }

    real_t operator()() const
    {
        return p_Dk_given_H1k / p_Dk_given_Tk;
    }
};

BIO_NS_END


std::ostream &
operator<<(std::ostream & os, const BIO_NS::BayesianClusterScore & score)
{
    boost::io::ios_all_saver ias(os);
    os.precision(3);

    os << score.p_Dk_given_H1k << "/" << score.p_Dk_given_Tk << "=" << score();
    return os;
}


BIO_NS_START

/** Scores clusters according to http://www.gatsby.ucl.ac.uk/~heller/bhcnew.pdf */
template <class Model>
struct BayesianClusterScorer
{
    typedef double real_t;
    typedef Model model_t;

    model_t & model;
    BayesianClusterScorer(model_t & model)
        : model(model)
    {
    }

    typedef BayesianClusterScore score_t;
    static score_t min_score() { return score_t(-1, 1); }

    template <class Cluster>
    score_t
    operator()(
        const Cluster & cluster) const
    {
        const real_t marginal_likelihood = model(cluster.data_set.begin(), cluster.data_set.end());
        return score_t(marginal_likelihood, marginal_likelihood);
    }

    /** Returns the score of the data set resulting from merging the two data sets. */
    template <class Cluster>
    score_t
    operator()(
        const Cluster & cluster_1,
        const Cluster & cluster_2) const
    {
        //merge the two data sets
        typename Cluster::data_set_t merged_data_set;
        set_union(
            cluster_1.data_set.begin(),
            cluster_1.data_set.end(),
            cluster_2.data_set.begin(),
            cluster_2.data_set.end(),
            std::inserter(merged_data_set, merged_data_set.begin()));

        score_t result;
        const real_t pi_k = 0.5; //?????
        result.p_Dk_given_H1k =
            pi_k * model(merged_data_set.begin(), merged_data_set.end());
        result.p_Dk_given_Tk =
            result.p_Dk_given_H1k
            + (1 - pi_k) * cluster_1.score.p_Dk_given_Tk * cluster_2.score.p_Dk_given_Tk;

        return result;
    }
};

/** Defines some types for given data types and scorers. */
template <class DataIt, class Scorer>
struct ClusterTraits
{
    typedef DataIt data_it;
    typedef Scorer scorer_t;

    typedef typename data_it::value_type data_t;
    typedef std::set<data_it> data_set_t;

    typedef Cluster<data_it, typename scorer_t::score_t> cluster_t;
    typedef typename cluster_t::data_set_t::iterator cluster_data_it;

    typedef std::list<cluster_t> cluster_list_t;
    typedef typename cluster_list_t::iterator cluster_it;
};

/** Clusters the data points according to pairs' scores. */
template <class DataIt, class Scorer, class Output>
void
cluster(
    DataIt data_begin,
    DataIt data_end,
    Scorer & scorer,
    Output & output)
{
    using namespace boost;
    using namespace std;

    typedef ClusterTraits<DataIt, Scorer> traits_t;

    //initialise the list of clusters with one entry for every data point
    typename traits_t::cluster_list_t cluster_list;
    for (DataIt d = data_begin;
        data_end != d;
        ++d)
    {
        typename traits_t::cluster_t cluster;
        cluster.data_set.insert(d);
        cluster.score = scorer(cluster);
        output(cluster_list.insert(cluster_list.begin(), cluster));
    }

    //while we have at least 2 clusters to merge
    while (1 < cluster_list.size())
    {
        //look for the best pair of clusters to merge in the hierarchy
        typename std::pair<typename traits_t::cluster_it, typename traits_t::cluster_it> best_clusters =
            make_pair(cluster_list.end(), cluster_list.end());
        typename Scorer::score_t best_score = Scorer::min_score();

        //for each pair of clusters, c1 and c2
        for (typename traits_t::cluster_it c1 = cluster_list.begin();
            cluster_list.end() != c1;
            ++c1)
        {
            typename traits_t::cluster_it c2 = c1;
            ++c2;
            for ( ;
                cluster_list.end() != c2;
                ++c2)
            {
                //score the pair of clusters
                const typename Scorer::score_t score = scorer(*c1, *c2);

                //is it the best score so far?
                if (best_score < score)
                {
                    //update best
                    best_score = score;
                    best_clusters.first = c1;
                    best_clusters.second = c2;
                }
            }
        }
        assert(best_clusters.first != best_clusters.second);
        assert(cluster_list.end() != best_clusters.first);
        assert(cluster_list.end() != best_clusters.second);

        //create and add the merged cluster
        typename traits_t::cluster_t merged_cluster;
        merged_cluster.score = best_score;
        set_union(
            best_clusters.first->data_set.begin(),
            best_clusters.first->data_set.end(),
            best_clusters.second->data_set.begin(),
            best_clusters.second->data_set.end(),
            inserter(merged_cluster.data_set, merged_cluster.data_set.begin()));

        //let the output know
        output(
            best_clusters.first,
            best_clusters.second,
            cluster_list.insert(cluster_list.begin(), merged_cluster));

        //remove the old data sets and add the merged one
        cluster_list.erase(best_clusters.first);
        cluster_list.erase(best_clusters.second);
    }
}

struct LoggingClusterer
{
    std::ostream & os;

    LoggingClusterer(std::ostream & os) : os(os)
    {
    }

    template <class ClusterIt>
    void
    operator()(
        const ClusterIt & cluster)
    {
        using namespace std;

        os << "Initial cluster: ";
        for (typename ClusterIt::value_type::data_set_t::const_iterator d = cluster->data_set.begin();
            cluster->data_set.end() != d;
            ++d)
        {
            os << **d << ",";
        }
        os << " score:" << cluster->score << endl;
    }

    template <class ClusterIt>
    void
    operator()(
        ClusterIt cluster_1,
        ClusterIt cluster_2,
        ClusterIt merged_cluster)
    {
        using namespace std;

        //print what we are doing
        cout << "Clustering: ";
        for (typename ClusterIt::value_type::data_set_t::const_iterator d = cluster_1->data_set.begin();
            cluster_1->data_set.end() != d;
            ++d)
        {
            os << **d << ",";
        }
        cout << " with ";
        for (typename ClusterIt::value_type::data_set_t::const_iterator d = cluster_2->data_set.begin();
            cluster_2->data_set.end() != d;
            ++d)
        {
            os << **d << ",";
        }
        os << " score:" << merged_cluster->score << endl;
    }
};

template <class DataIt, class Scorer>
struct TreeBuildingClusterer
{
    typedef DataIt data_it;
    typedef Scorer scorer_t;
    typedef ClusterTraits<data_it, scorer_t> traits_t;

    typedef std::pair<bool, typename traits_t::data_it> vertex_property;
    typedef boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::directedS,
        vertex_property,
        typename Scorer::score_t> tree_t;
    typedef typename boost::graph_traits<tree_t>::vertex_iterator tree_vertex_it;
    typedef typename boost::graph_traits<tree_t>::vertex_descriptor tree_vertex;
    tree_t tree;

    /** A map from clusters to tree nodes. */
    typedef std::map<typename traits_t::cluster_t *, tree_vertex> cluster_vertex_map_t;
    cluster_vertex_map_t cluster_vertex_map;

    LoggingClusterer logger;
    bool logging_on;

    TreeBuildingClusterer(bool logging_on = false) : logger(std::cout), logging_on(logging_on)
    {
    }

    void
    reset()
    {
        cluster_vertex_map.clear();
        tree.clear();
    }

    void
    operator()(
        typename traits_t::cluster_it cluster)
    {
        using namespace std;
        using namespace boost;

        if (logging_on)
        {
            logger(cluster);
        }

        cluster_vertex_map[&*cluster] = add_vertex(make_pair(true, *cluster->data_set.begin()), tree);
    }

    void
    operator()(
        typename traits_t::cluster_it cluster_1,
        typename traits_t::cluster_it cluster_2,
        typename traits_t::cluster_it merged_cluster)
    {
        using namespace std;
        using namespace boost;

        if (logging_on)
        {
            logger(cluster_1, cluster_2, merged_cluster);
        }

        //find the vertex descriptors for cluster_1 and cluster_2
        assert(cluster_vertex_map.end() != cluster_vertex_map.find(&*cluster_1));
        assert(cluster_vertex_map.end() != cluster_vertex_map.find(&*cluster_2));

        cluster_vertex_map[&*merged_cluster] = add_vertex(make_pair(false, typename traits_t::data_it()), tree);

        add_edge(cluster_vertex_map[&*merged_cluster], cluster_vertex_map[&*cluster_1], tree);
        add_edge(cluster_vertex_map[&*merged_cluster], cluster_vertex_map[&*cluster_2], tree);
    }
};

BIO_NS_END


#endif //BIO_BAYESIAN_HIERARCHICAL_CLUSTERING_H_
