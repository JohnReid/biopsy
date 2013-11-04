
_base_url = 'http://www.biomart.org/biomart/martservice'

def _attribute_xml( attribute ):
    "Returns xml suitable for inclusion into query"
    return '<Attribute name = "%s" />' % attribute

def _filter_xml( name, value ):
    "Returns xml suitable for inclusion into query"
    return '<Filter name = "%s" value = "%s"/>' % ( name, value )

def _meta_query( fields ):
    import urllib
    return urllib.urlopen(
            _base_url,
            urllib.urlencode( fields )
    )

def _convert_tab_separated( handle ):
    for l in handle:
        l = l.strip()
        if '' != l:
            yield l.strip().split( '\t' )

def registry():
    for f in _convert_tab_separated(
            _meta_query(
                    {
                            'type' : 'registry',
                    }
            )
    ): yield f
#print '\n'.join([str(x) for x in registry()])
#raise


def datasets( mart = 'ensembl' ):
    for f in _convert_tab_separated(
            _meta_query(
                    {
                            'type' : 'datasets',
                            'mart' : mart
                    }
            )
    ): yield f
#print '\n'.join([str(x) for x in datasets( 'sequence' ) if -1 != x[1].find('mmus')])
#raise

def configuration( dataset = 'mmusculus_gene_ensembl' ):
    for f in _convert_tab_separated(
            _meta_query(
                    {
                            'type' : 'configuration',
                            'dataset' : dataset
                    }
            )
    ): yield f
#'\n'.join([str(c) for c in configuration()])
#raise

def attributes( dataset = 'mmusculus_gene_ensembl' ):
    for f in _convert_tab_separated(
            _meta_query(
                    {
                            'type' : 'attributes',
                            'dataset' : dataset
                    }
            )
    ): yield f

def filters( dataset = 'mmusculus_gene_ensembl' ):
    for f in _convert_tab_separated(
            _meta_query(
                    {
                            'type' : 'filters',
                            'dataset' : dataset
                    }
            )
    ): yield f

#for f in filters():
#       if f[0].find('go') != -1:
#               print f

def _query( xml ):
    """Execute query and return result"""
    import urllib
    data = urllib.urlencode( [ ( 'query', xml ) ] )
    return urllib.urlopen(
            _base_url,
            data
    )

_query_xml_template = """
<!DOCTYPE Query>
<Query
        virtualSchemaName = "default"
        Header = "1"
        count = "%s"
        softwareVersion = "0.5" >
        <Dataset name = "%s" interface = "default" >
                %s
                %s
        </Dataset>
</Query>"""

class Query( object ):
    """A biomart query"""

    def __init__(
            self,
            filters = [ ],
            attributes = [ ],
            dataset_name = 'mmusculus_gene_ensembl'
    ):
        self.filters = filters
        self.attributes = attributes
        self.dataset_name = dataset_name

    def _build_xml( self, count = False ):
        """Builds an xml query to send to biomart web service"""
        if count: count_string = '1'
        else: count_string = ''
        return _query_xml_template % (
                count_string,
                self.dataset_name,
                '\n'.join( [ _attribute_xml( a ) for a in self.attributes ] ),
                '\n'.join( [ _filter_xml( n, v ) for n, v in self.filters ] )
        )

    def get_count( self ):
        handle = _query( self._build_xml( count = True ) )
        result = handle.read()
        try:
            return int( result )
        except ValueError:
            raise RuntimeError( 'Did not count from server: %s' % result )


    def __call__( self ):
        return _query( self._build_xml( count = False ) )

if '__main__' == __name__:
    Q = Query(
            filters = [
                    #( 'ensembl_gene_id', 'ENSMUSG00000029754' ),
                    #( 'ensembl_gene_id', 'ENSMUSG00000020330' ),
                    ( 'go', 'GO:0004872' ),
                    ( 'go', 'GO:0005540' ),
            ],
            attributes = [
                    'ensembl_gene_id',
                    'ensembl_transcript_id',
            ],
    )
    print Q.get_count()
