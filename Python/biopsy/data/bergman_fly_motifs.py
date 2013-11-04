#
# Copyright John Reid 2010
#

"""
Code to deal with the Bergman curated set of fly motifs.
"""

import os, biopsy.data as D, numpy as N
import xml.etree.ElementTree as ET


def xms_filename():
    "@return: The filename of the XMS file where the motifs are stored."
    return os.path.join(D.data_dir(), "Bergman-Fly-Motifs", "SelexConsensus1.1.xms")


def parse_xms(f):
    "Parse an XMS file. @return: Yield the motifs."
    tree = ET.parse(f)
    root = tree.getroot()
    for motif in root.findall('motif'):
        name = motif.find('name').text
        weightmatrix = motif.find('weightmatrix')
        columns = int(weightmatrix.get('columns'))
        alphabet = weightmatrix.get('alphabet')
        if 'DNA' == alphabet:
            alphabet_size = 4
        matrix = N.zeros((columns, alphabet_size))
        for column in weightmatrix.findall('column'):
            pos = int(column.get('pos'))
            for weight in column.findall('weight'):
                symbol = weight.get('symbol')
                if 'thymine' == symbol:
                    b = 3
                elif 'guanine' == symbol:
                    b = 2
                elif 'cytosine' == symbol:
                    b = 1
                elif 'adenine' == symbol:
                    b = 0
                else:
                    raise RuntimeError('Unrecognized symbol: ' + symbol)
                value = float(weight.text)
                matrix[pos,b] = value
        properties = dict()
        for prop in motif.findall('prop'):
            key = prop.find('key')
            value = prop.find('value')
            assert None != key
            assert None != value
            properties[key.text] = value.text
        threshold = float(motif.find('threshold').text)
        yield name, alphabet, matrix, properties, threshold


def write_as_custom_pssm(f, id_, name, matrix, comments=None, url=None, field_width=3, scale=1):
    """
    Write the motif as a custom PSSM to the file, f.

#
# Drosophila Hunchback from JASPAR
#
ID  DN-000001
NA  D$Hunchback
WI  10
PO  01 02 03 04 05 06 07 08 09 10
CA  01 06 09 04 13 16 16 14 15 09
CC  05 08 03 03 01 00 00 00 01 02
CG  08 02 04 01 00 00 00 02 00 02
CT  02 00 00 08 02 00 00 00 00 03
IU   G  C  A  T  A  A  A  A  A  A
UR  None
    """
    if None != comments:
        print >> f, '#'
        for comment in comments:
            print >> f, '# %s' % comment
        print >> f, '#'

    print >> f, 'ID  %s' % id_
    print >> f, 'NA  %s' % name
    print >> f, 'WI  %s' % len(matrix)
    print >> f, 'PO  %s' % ' '.join('%*d' % (field_width, i+1) for i in xrange(len(matrix)))
    for b, tag in enumerate(('CA', 'CC', 'CG', 'CT')):
        print >> f, '%s  %s' % (tag, ' '.join('%*d' % (field_width, int(v)) for v in matrix[:,b]*scale))
    print >> f, 'UR  %s' % (None != url and url or 'None')


def normalise_matrix(matrix):
    "@return: A normalised version of the argument."
    return (matrix.T / matrix.sum(axis=1)).T


def smooth_matrix_with_pseudo_count(matrix, pseudo_count):
    "@return: A smoothed version the matrix using the given pseudo counts."
    smoothed = matrix + pseudo_count
    return normalise_matrix(smoothed)


def write_matrix_to_file(f, id_, name, alphabet, matrix, properties, threshold, scale=1):
    "Write the matrix to the file in the custom PSSM format."
    comments = [
        'PSSM parsed from set of fly TFs curated by Bergman.'
    ]
    comments.extend('%20s : %s' % (k, v) for k, v in properties.iteritems())
    write_as_custom_pssm(f, id_, name, matrix, comments=comments, scale=scale)




if '__main__' == __name__:

    import sys

    output_dir = '/home/john/Data/custom-pssms'
    pssm_set_tag = 'BG'
    scale = 30

    pssm_set_f = open(os.path.join(output_dir, 'bergman-fly.pssm_set'), 'w')
    print >> pssm_set_f, '#'
    print >> pssm_set_f, '# Set of fly TFs curated by Bergman.'
    print >> pssm_set_f, '# PSSMs were scaled as if there were %d observations.' % scale
    print >> pssm_set_f, '#'
    for i, (name, alphabet, matrix, properties, threshold) in enumerate(parse_xms(open(xms_filename()))):
        id_ = '%s-%06d' % (pssm_set_tag, i+1)
        print id_, name
        print >> pssm_set_f, id_
        f = open(os.path.join(output_dir, '%s.pssm' % id_), 'w')
        properties['Equivalent # observations'] = str(scale)
        write_matrix_to_file(f, id_, name, alphabet, matrix, properties, threshold, scale=scale)
        f.close()
    pssm_set_f.close()
