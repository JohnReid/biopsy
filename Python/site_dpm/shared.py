#
# Copyright John Reid 2008, 2009, 2010
#


"""
Various shared code (options, logging, etc...)
"""

import logging
logging.basicConfig(level=logging.INFO)

from cookbook.workflow import *
import default_options as options
initialise_workflow(options)

from time import time
from itertools import chain, imap, ifilter, islice
import cookbook.cache_decorator, biopsy, numpy, os, itertools, cPickle, sys, shutil, biopsy
from cookbook.simple_logging import create_log_exceptions_decorator as log_exceptions






def global_cached_method(name):
    "Decorator to store output of methods in site dpm data directory (global across all runs of site DPM)."
    return cookbook.cache_decorator.pickled_cached_method(os.path.join(get_site_dpm_data_dir(), '%s.pickle' % name))



def chain_biomart_queries(query, values, max_query_size=150):
    "Takes a list of values to be passed to a query function and splits it into several smaller queries."
    num_queries = len(values) / max_query_size + 1
    return chain(*[query(islice(values, i, None, num_queries)) for i in xrange(num_queries)])



def write_gene_set_with_names(f, genes, ensembl_names):
    """Writes gene set to file object with their names."""
    for g in genes:
        print >> f, '%s,%s' % (g, ensembl_names.get(g, '<unknown>'))


def yield_ensembl_names(gene_ids):
    import biopsy.identifiers.biomart as biomart, csv
    query = biomart.new_query()
    dataset = biomart.add_dataset(query, 'mmusculus_gene_ensembl')
    biomart.add_filter(dataset, 'ensembl_gene_id', ",".join(imap(str, gene_ids)))
    biomart.add_attribute(dataset, 'ensembl_gene_id')
    biomart.add_attribute(dataset, 'external_gene_id')
    for row in biomart.yield_csv_query_results(query):
        yield row[0], row[1]


def get_ensembl_names(gene_ids):
    "@return: A map from the given genes to their external ids."
    logging.info("Querying Ensembl biomart for %d genes' names", len(gene_ids))
    return dict(chain_biomart_queries(yield_ensembl_names, gene_ids))



#
# Set up some directories to write output to
#
_programs_dir = os.path.join(options.output_dir, 'programs')
"the directory to put program specific info into."

def get_programs_dir():
    "@return: The directory to put program specific info into."
    #ensure_dir_exists(_programs_dir)
    return _programs_dir


_summaries_dir = os.path.join(options.output_dir, 'summaries')
"the directory to put DPM summaries into"

def get_summaries_dir():
    "@return: The directory to put DPM summaries into."
    #ensure_dir_exists(_summaries_dir)
    return _summaries_dir


_site_dpm_data_dir = os.path.join(biopsy.get_data_dir(), 'site-dpm')
"The directory where results that are reused across runs are cached."

def get_site_dpm_data_dir():
    "@return: The directory where results that are reused across runs are cached."
    #ensure_dir_exists(_site_dpm_data_dir)
    return _site_dpm_data_dir



try:
    import pylab
except:
    import warnings
    warnings.warn('Could not set matplotlib figure size')
    print sys.exc_info()


@global_cached_method('all-ensembl-names')
def get_all_ensembl_names():
    "@return: All ensembl names."
    from biopsy.identifiers.biomart import quick_query
    logging.info("Querying Ensembl biomart for all mouse genes' names")
    return dict(
        quick_query(
            dataset='mmusculus_gene_ensembl',
            attributes=('ensembl_gene_id', 'external_gene_id')
        )
    )
