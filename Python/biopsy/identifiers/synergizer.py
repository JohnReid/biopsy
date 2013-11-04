#
# Copyright John Reid 2009
#

import simplejson as json
import httplib
from itertools import chain

"""
Code to use the Syngergizer web service:

http://llama.med.harvard.edu/synergizer/doc/
"""


server = 'llama.med.harvard.edu'
service_url = '/cgi/synergizer/serv'
headers = {
  'Content-type' : 'application/json',
}




def make_request(method, params):
    """
    Make a request to the Synergizer web service.

    @return: result
    """
    conn = httplib.HTTPConnection(server)
    payload = json.dumps({'method':method, 'params':params, 'id':0})
    conn.request("POST", service_url, body=payload, headers=headers)
    response = conn.getresponse()
    response_data = response.read()
    decoded_response = json.loads(response_data)
    conn.close()
    if decoded_response['error']:
        raise RuntimeError('Syngerizer error: %s' % str(decoded_response['error']['message']))
    return decoded_response['result']


def version():
    """
    @return: The version of the Synergizer service.
    """
    return make_request('version', [])


def available_authorities():
    """
    @return: The available authorities in the Synergizer service.
    """
    return make_request('available_authorities', [])



def available_species(authority):
    """
    @return: The available species for the authority in the Synergizer service.
    """
    return make_request('available_species', [authority])



def available_domains(authority, species):
    """
    @return: The available domains for the authority and species in the Synergizer service.
    """
    return make_request('available_domains', [authority, species])



def available_ranges(authority, species, domain):
    """
    @return: The available ranges for the authority, species and domain in the Synergizer service.
    """
    return make_request('available_ranges', [authority, species, domain])



def translate(authority, species, domain, range, ids):
    """
    @return: The translated ids.
    """
    return make_request(
        'translate',
        [
            {
                'authority' : authority,
                'species'   : species,
                'domain'    : domain,
                'range'     : range,
                'ids'       : ids
            }
        ]
    )


def more_than_one_translation(x):
    return len(x) > 2

def how_many_have_translations(translated):
    count = 0
    for x in those_with_translation(translated):
        count += 1
    return count

def get_translations(translated):
    """@return: A set of the identifiers that were mapped to."""
    return set(chain(*(x[1:] for x in translated if x[1])))

def those_with_translation(translated):
    """Yield those translations that succeeded."""
    for x in translated:
        if x[1]:
            yield x

def yield_mapping(translated):
    """Yield (from, to) tuples where from is an identifier in the domain and to is a sequence of identifiers in the range."""
    for x in those_with_translation(translated):
        yield x[0], x[1:]



if '__main__' == __name__:
    print 'Version: %s' % version()
    print 'Available authorities: %s' % ', '.join(available_authorities())
    print 'Available species: %s' % ', '.join(available_species('ensembl'))
    print 'Available domains: %s' % ', '.join(available_domains('ensembl', 'Homo sapiens'))
    print 'Available ranges: %s' % ', '.join(available_ranges('ensembl', 'Homo sapiens', 'ensembl_gene_id'))
    translated = translate('ensembl', 'Homo sapiens', 'unigene', 'ensembl_gene_id', [ 'Hs.239', 'Hs.715518' ])
    print 'Translated: %s' % ', '.join('->'.join(map(str, x)) for x in translated)
