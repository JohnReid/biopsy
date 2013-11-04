#
# Copyright John Reid 2007
#


"""
Code to access KEGG database via DBGet interface. (http://www.genome.jp/kegg/docs/keggapi_manual.html)

The WSDL file for the KEGG API can be found at:
    * <URL:http://soap.genome.jp/KEGG.wsdl>

Terminology
    * 'org' is a three-letter (or four-letter) organism code used in KEGG. The list can be found at (see the description of the list_organisms method below):
          o <URL:http://www.genome.jp/kegg/catalog/org_list.html>
    * 'db' is a database name used in GenomeNet service. See the description of the list_databases method below.
    * 'entry_id' is a unique identifier of which format is the combination of the database name and the identifier of an entry joined by a colon sign as 'database:entry' (e.g. 'embl:J00231' means an EMBL entry 'J00231'). 'entry_id' includes 'genes_id', 'enzyme_id', 'compound_id', 'drug_id', 'glycan_id', 'reaction_id', 'pathway_id' and 'motif_id' described in below.
    * 'genes_id' is a gene identifier used in KEGG/GENES which consists of 'keggorg' and a gene name (e.g. 'eco:b0001' means an E. coli gene 'b0001').
    * 'enzyme_id' is an enzyme identifier consisting of database name 'ec' and an enzyme code used in KEGG/LIGAND ENZYME database. (e.g. 'ec:1.1.1.1' means an alcohol dehydrogenase enzyme)
    * 'compound_id' is a compound identifier consisting of database name 'cpd' and a compound number used in KEGG COMPOUND / LIGAND database (e.g. 'cpd:C00158' means a citric acid). Note that some compounds also have 'glycan_id' and both IDs are accepted and converted internally by the corresponding methods.
    * 'drug_id' is a drug identifier consisting of database name 'dr' and a compound number used in KEGG DRUG / LIGAND database (e.g. 'dr:D00201' means a tetracycline).
    * 'glycan_id' is a glycan identifier consisting of database name 'gl' and a glycan number used in KEGG GLYCAN database (e.g. 'gl:G00050' means a Paragloboside). Note that some glycans also have 'compound_id' and both IDs are accepted and converted internally by the corresponding methods.
    * 'reaction_id' is a reaction identifier consisting of database name 'rn' and a reaction number used in KEGG/REACTION (e.g. 'rn:R00959' is a reaction which catalyze cpd:C00103 into cpd:C00668)
    * 'pathway_id' is a pathway identifier consisting of 'path' and a pathway number used in KEGG/PATHWAY. Pathway numbers prefixed by 'map' specify the reference pathway and pathways prefixed by the 'keggorg' specify pathways specific to the organism (e.g. 'path:map00020' means a reference pathway for the cytrate cycle and 'path:eco00020' means a same pathway of which E. coli genes are marked).
    * 'motif_id' is a motif identifier consisting of motif database names ('ps' for prosite, 'bl' for blocks, 'pr' for prints, 'pd' for prodom, and 'pf' for pfam) and a motif entry name. (e.g. 'pf:DnaJ' means a Pfam database entry 'DnaJ').
    * 'ko_id' is a KO identifier consisting of 'ko' and a ko number used in KEGG/KO. KO (KEGG Orthology) is an classification of orthologous genes defined by KEGG (e.g. 'ko:K02598' means a KO group for nitrite transporter NirC genes).
    * 'ko_class_id' is a KO class identifier which is used to classify 'ko_id' hierarchically (e.g. '01110' means a 'Carbohydrate Metabolism' class).
          o <URL:http://www.genome.jp/dbget-bin/get_htext?KO>
    * 'offset' and 'limit' are both an integer and used to control the number of the results returned at once. Methods having these arguments will return first 'limit' results starting from 'offset'th.
    * 'fg_color_list' is a list of colors for the foreground (corresponding to the texts and borders of the objects on the KEGG pathway map).
    * 'bg_color_list' is a list of colors for the background (corresponding to the inside of the objects on the KEGG pathway map).

Related site:
    * <URL:http://www.genome.jp/kegg/kegg3.html>

Many of the KEGG API methods will return a set of values in a complex data structure as described below. This section summarizes all kind of these data types. Note that, the retuened values for the empty result will be

    * an empty array -- for the methods which return ArrayOf'OBJ'
    * an empty string -- for the methods which return String
    * -1 -- for the methods which return int
    * NULL -- for the methods which return any other 'OBJ'
"""

from . import lazy

wsdl = 'http://soap.genome.jp/KEGG.wsdl'
def _server():
    "Get a server handle"
    from SOAPpy import WSDL
    return WSDL.Proxy(wsdl)
server = lazy.LazyInitialiser(_server)

def _str_for_structType_list(structType_list):
    return "\n".join(("%s:%s" % (e.entry_id, e.definition)) for e in structType_list)

def print_array_of_definition(a):
    for e in a:
        print "%s %s" % (e.entry_id, e.definition)
