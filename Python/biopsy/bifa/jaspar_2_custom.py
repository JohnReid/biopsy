#!/usr/bin/env python

#
# Copyright John Reid 2012
#


"""
Converts JASPAR matrix file into BiFa custom PSSM files.
"""

from __future__ import with_statement
import sys, warnings

def parse_jaspar(f):
    return [
        map(int, l.split())
        for l in f
    ]


def jaspar_url(jaspar_id):
    if jaspar_id.startswith('CN'):
        return 'http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=%s&rm=present&collection=CNE' % jaspar_id
    if jaspar_id.startswith('MA'):
        return 'http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=%s&rm=present&collection=CORE' % jaspar_id
    if jaspar_id.startswith('MF'):
        return 'http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=%s&rm=present&collection=FAM' % jaspar_id
    if jaspar_id.startswith('PB'):
        return 'http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=%s&rm=present&collection=PBM' % jaspar_id
    if jaspar_id.startswith('PL'):
        return 'http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=%s&rm=present&collection=PBM_HLH' % jaspar_id
    if jaspar_id.startswith('PH'):
        return 'http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=%s&rm=present&collection=PBM_HOMEO' % jaspar_id
    if jaspar_id.startswith('PF'):
        return 'http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=%s&rm=present&collection=PHYLO_FACTS' % jaspar_id
    if jaspar_id.startswith('POL'):
        return 'http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=%s&rm=present&collection=POLII' % jaspar_id
    if jaspar_id.startswith('SA'):
        return 'http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=%s&rm=present&collection=SPLICE' % jaspar_id
    else:
        return 'http://jaspar.genereg.net/'



prefix = sys.argv[1] # get prefix from command line
collections = dict()

#
# For each PSSM described in stdin
#
with open('jaspar_custom_map.txt', 'w') as map_f:
    for i, line in enumerate(sys.stdin):
        line = line.strip()
        if not line.startswith('>'):
            warnings.warn('Expecting line to start with ">": %s' % line)
        idx = line.find(' ')
        jaspar_id = line[1:idx]
        tf = line[idx+1:]
        custom_id = '%s-%06d' % (prefix, i)
        print >> map_f, '%s : %s : %s' % (custom_id, jaspar_id, tf)
        
        #
        # Work out the collection
        #
        collection = jaspar_id[:2]
        if collection not in collections:
            collections[collection] = open('JASPAR-%s.pssm_set' % collection, 'w')
            print >> collections[collection], '# A collection of JASPAR pssms: %s : http://jaspar.genereg.net/' % collection
        print >> collections[collection], custom_id
        
        #
        # Get the matrix
        #
        matrix = []
        for base in ('A', 'C', 'G', 'T'):
            line = sys.stdin.next()
            if base != line[0]:
                raise ValueError('Expecting first character of matrix row to be "%s": %s' % (base, line))
            start = line.find('[')
            end = line.find(']')
            matrix.append(map(float, line[start+1:end].split()))
        length = len(matrix[0])
        
        #
        # Write the custom PSSM
        #
        custom_filename = '%s.pssm' % custom_id
        with open(custom_filename, 'w') as out:
            print >> out, "NA  %s: %s" % (jaspar_id, tf)
            print >> out, "WI  %d" % length
            print >> out, "PO  %s" % ' '.join('%02d' % (j+1) for j in xrange(length))
            for code, row in (('CA', matrix[0]), ('CC', matrix[1]), ('CG', matrix[2]), ('CT', matrix[3])): 
                print >> out, '%s  %s' % (code, ' '.join('%02f' % entry for entry in row))
            #IU  G  G  A  C  A  C  G  T  G  G  C  N 
            print >> out, "UR  %s" % jaspar_url(jaspar_id)


for f in collections.values():
    f.close()
