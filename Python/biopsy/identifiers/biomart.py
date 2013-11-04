#
# Copyright John Reid 2007,2008,2009
#

"""
Code to access biomart.

www.ebi.ac.uk/biomart/

www.biomart.org/

    query = new_query()
    dataset = add_dataset(query, 'mmusculus_gene_ensembl')
    add_attribute(dataset, 'ensembl_gene_id')
    add_attribute(dataset, 'external_gene_id')
    add_filter(dataset, 'ensembl_gene_id', 'ENSMUSG00000056216')
    for row in csv.reader(execute_query(query), delimiter=','):
        print row

or

    for row in quick_query(
        dataset='mmusculus_gene_ensembl',
        attributes=['ensembl_gene_id', 'external_gene_id'],
        filters=[('ensembl_gene_id', 'ENSMUSG00000056216,ENSMUSG00000034957')],
    ):
        print row

"""

import urllib, csv, sys
import xml.etree.ElementTree as ET

_dataset_config_version = "0.7"

def print_query(query, file):
    'Print the query to the file.'
    ET.ElementTree(query).write(file)

def new_query(
        virtualSchemaName='default',
        formatter='CSV',
        header="0",
        uniqueRows="0",
        count="",
        datasetConfigVersion=_dataset_config_version
):
    'Create a new query.'
    query = ET.Element('Query')
    query.set('virtualSchemaName', virtualSchemaName)
    query.set('formatter', formatter)
    query.set('header', header)
    query.set('uniqueRows', uniqueRows)
    query.set('count', count)
    query.set('datasetConfigVersion', datasetConfigVersion)
    return query

def query_to_string(query):
    return """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
%s""" % ET.tostring(query)

def add_dataset(query, name, interface='default'):
    'Add a dataset to a query'
    dataset = ET.SubElement(query, "Dataset")
    dataset.set('name', name)
    dataset.set('interface', interface)
    return dataset

def add_attribute(dataset, name):
    'Add an attribute to a dataset'
    attribute = ET.SubElement(dataset, "Attribute")
    attribute.set('name', name)
    return attribute

def add_filter(dataset, name, value):
    'Add a filter to a dataset'
    filter = ET.SubElement(dataset, "Filter")
    filter.set('name', name)
    filter.set('value', value)
    return filter

def query_url(query, webservice = 'http://www.biomart.org/biomart/martservice?'):
    return webservice + urllib.urlencode([('query',query_to_string(query))])

max_query_len = 4000

def execute_query(query, webservice = 'http://www.biomart.org/biomart/martservice?'):
    "Execute the query at the given web service"
    url = query_url(query, webservice)
    if len(url) > max_query_len:
        raise RuntimeError('Query URL very long (%d), will probably be rejected.' % len(url))
    return urllib.urlopen(url)

def yield_csv_query_results(query, webservice = 'http://www.biomart.org/biomart/martservice?'):
    "Yield the rows returned from csv query."
    import csv
    for row in csv.reader(execute_query(query), delimiter=','):
        yield row

def split_big_list(iterable, chunk_size=200):
    l = list(iterable)
    for i in xrange(0, len(l), chunk_size):
        yield l[i:i+chunk_size]

def quick_query(dataset, attributes, filters=()):
    """
    Convenience function to perform a query returning the named attributes and using the (name, value) pairs
    in filters.
    """
    query = new_query()
    dataset = add_dataset(query, dataset)
    for attribute in attributes:
        add_attribute(dataset, attribute)
    for name, value in filters:
        add_filter(dataset, name, value)
    for row in csv.reader(execute_query(query), delimiter=','):
        yield row


if '__main__' == __name__:
    _typical_query = """<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName="default" formatter='CSV' header="0" uniqueRows="0" count="" datasetConfigVersion="0.5" >

        <Dataset name = "compara_hsap_mmus_orthologs" interface = "default" >
            <Attribute name = "homol_description" />
            <Attribute name = "hsap_gene_stable_id" />
            <Attribute name = "mmus_gene_stable_id" />
        </Dataset>
    </Query>
    """
    #check we can parse it
    tree = ET.fromstring(_typical_query)

    query = new_query()
    dataset = add_dataset(query, 'mmusculus_gene_ensembl')
    add_attribute(dataset, 'ensembl_gene_id')
    add_attribute(dataset, 'external_gene_id')
    add_filter(dataset, 'ensembl_gene_id', 'ENSMUSG00000056216')
    for row in csv.reader(execute_query(query), delimiter=','):
        print row

    for row in quick_query(
        dataset='mmusculus_gene_ensembl',
        attributes=['ensembl_gene_id', 'external_gene_id'],
        filters=[('ensembl_gene_id', 'ENSMUSG00000056216,ENSMUSG00000034957')],
    ):
        print row
