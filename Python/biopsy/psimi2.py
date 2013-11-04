#
# Copyright John Reid 2007
#

import biopsy
T = biopsy.transfac
from xml.etree import ElementTree as ET
import boost.graph as bgl
import itertools, cookbook, sys, copy

"""
Code to read PSI-MI XML format for protein-protein interactions.

See:
        http://www.psidev.info/index.php?q=node/60
        http://psidev.sourceforge.net/mi/rel25/doc/

To get the interaction network for IntAct try:
        networks['INTACT']
"""

class Interactor(object):
    """
    A node in an interaction network
    """

    def __init__(self, label, primary_ref, secondary_refs, aliases, taxid):
        if not primary_ref:
            #print secondary_refs
            raise ValueError('Must have primary reference for interactor')
        self.label = label.upper()
        self.primary_ref = primary_ref
        self.secondary_refs = secondary_refs
        self.aliases = [ a.upper() for a in aliases ]
        self.taxid = taxid

    def __repr__(self):
        return str(self.primary_ref)

    def __str__(self):
        return self.label

    def url(self):
        "A url for this interactor."
        try:
            return self.primary_ref.url
        except:
            for ref in self.secondary_refs:
                try:
                    return ref.url
                except:
                    pass
        return ''


class InteractionNetwork(object):
    """
    A network of interactions, typically protein-protein interactions
    """

    def __init__(self):
        self.g = bgl.Graph()
        "A boost.graph.Graph that stores the network."

        self.interactor_props = self.g.add_vertex_property('interactors', type='object')
        "A boost.graph vertex property that maps from vertices to interactors."

        self.primary_2_node = dict()
        "A dict mapping from primary references to nodes (vertices)"

        self.ref_2_nodes = cookbook.DictOfLists()
        "A dict mapping from references to lists of nodes (vertices)"

        self.name_2_nodes = cookbook.DictOfLists() # maps names to nodes
        "A dict mapping from names to lists of nodes (vertices)"

    def __deepcopy__(self, memo):
        """
        Make a deep copy my creating a new network and reinserting the interactors.
        Hence the interactors are not duplicated.
        """
        #print '__deepcopy__(%s)' % str(memo)
        result = self.__class__()
        memo[id(self)] = result
        result.__init__()
        for e in self.g.edges:
            u, v = self.g.source(e), self.g.target(e)
            i1, i2 = self.interactor_for_node(u), self.interactor_for_node(v)
            result.add_interaction(i1, i2)
        return result

    def add_interaction(self, interactor_1, interactor_2):
        "Add an interaction between 2 interactors to the network"
        self.g.add_edge(self.node_for(interactor_1), self.node_for(interactor_2))

    def interactor_for_node(self, node):
        "The interactor for this node"
        return self.interactor_props[node]

    def nodes_for_ref(self, ref):
        "The nodes that have this reference"
        if ref in self.primary_2_node:
            result = [self.primary_2_node[ref]]
        else:
            result = []
        if ref in self.ref_2_nodes:
            result.extend(self.ref_2_nodes[ref])
        return result

    def interactors_for_ref(self, ref):
        "The interactors that have this reference"
        return [self.interactor_for_node(node) for node in self.nodes_for_ref(ref)]

    def nodes_for_name(self, name):
        "The nodes that have this name"
        name = name.upper()
        return name in self.name_2_nodes and self.name_2_nodes[name] or []

    def interactors_for_name(self, name):
        "The interactors that have this name"
        return [self.interactor_for_node(node) for node in self.nodes_for_name(name)]

    def dbs_in_network(self):
        return set(r.db for r in self.ref_2_nodes.keys())

    def node_for(self, interactor):
        "The node for this interactor, adds one if necessary"
        node = self.find_node_for(interactor)
        if None == node:
            node = self.g.add_vertex()
            self.interactor_props[node] = interactor
            self.primary_2_node[interactor.primary_ref] = node

            self.ref_2_nodes[interactor.primary_ref].append(node)
            for r in interactor.secondary_refs: # add the secondary references
                self.ref_2_nodes[r].append(node)

            self.name_2_nodes[interactor.label].append(node)
            for alias in interactor.aliases: # add the aliases
                self.name_2_nodes[alias].append(node)
        return node

    def find_node_for(self, interactor):
        "The node for this interactor, or None"
        return interactor.primary_ref in self.primary_2_node and self.primary_2_node[interactor.primary_ref] or None

    def summary(self):
        "Returns a string summarising the network's size and content."
        return """# interactors: %d
# interactions: %d
DBs: %s""" % (
self.g.num_vertices(),
self.g.num_edges(),
','.join(str(db) for db in self.dbs_in_network())
        )

    def typical_primary_ref(self):
        "Returns the primary ref of the first node"
        return self.interactor_props[self.g.vertices.next()].primary_ref

    def add_labels_to_graph(self):
        """
        Add a property to the graph that contains the label for each interactor.
        Returns the graph label property.
        """
        label = self.g.add_vertex_property( 'label', 'string' )
        for v in self.g.vertices:
            label[v] = self.interactor_for_node(v).label
        return label

    def add_urls_to_graph(self):
        """
        Add a property to the graph that contains a url for each interactor.
        Returns the graph url property.
        """
        url = self.g.add_vertex_property( 'URL', 'string' )
        for v in self.g.vertices:
            url[v] = self.interactor_for_node(v).url()
        return url

    def normalise(self):
        """
        Make sure only nodes referred to are those in graph.
        """
        to_remove = []
        for ref, v in self.primary_2_node.iteritems():
            if v not in self.g.vertices:
                to_remove.append(ref)
        for ref in to_remove:
            del self.primary_2_node[ref]

        to_remove = []
        for ref, vs in self.ref_2_nodes.iteritems():
            self.ref_2_nodes[ref] = [v for v in vs if v in self.g.vertices]
            if not len(self.ref_2_nodes[ref]):
                to_remove.append(ref)
        for ref in to_remove:
            del self.ref_2_nodes[ref]

        to_remove = []
        for name, vs in self.name_2_nodes.iteritems():
            self.name_2_nodes[name] = [v for v in vs if v in self.g.vertices]
            if not len(self.name_2_nodes[name]):
                to_remove.append(name)
        for name in to_remove:
            del self.name_2_nodes[name]

    def remove_nodes_not_in_identifiers(self, identifiers):
        """
        Remove nodes in the graph that are not in the given sequence of identifiers.
        """
        to_keep = []
        for i in identifiers:
            to_keep.extend(self.nodes_for_ref(i))
        to_remove = [v for v in self.g.vertices if v not in to_keep]
        for v in to_remove:
            self.g.clear_vertex(v)
            self.g.remove_vertex(v)
        self.normalise()

def ref_for_node(node):
    "Get the reference associated with this xref node"
    db = node.get('db')
    try:
        if db == 'Entrez Protein':
            return biopsy.DbRef.parse_as(node.get('id'), biopsy.transfac.db.entrez_protein)
        if db == 'PROTEIN GI':
            return biopsy.DbRef(biopsy.db.entrez_protein, "", int(node.get('id')))
        if db == 'MGI':
            return biopsy.DbRef(biopsy.db.mgi, "", int(node.get('id')))
        return biopsy.DbRef.parse(node.get('id'))
    except:
        return None

def secondary_refs_for_xref(node, ns):
    for node in node.findall('{%s}secondaryRef' % ns):
        r = ref_for_node(node)
        if r:
            yield r

def aliases_for_interactor(names, ns):
    return [ alias.text for alias in names.findall('{%s}alias' % ns) if alias.text ]

def label_for_interactor(names, ns):
    label_node = names.find('{%s}shortLabel' % ns)
    if None == label_node:
        label_node = names.find('{%s}fullName' % ns)
    return label_node.text

def organism_for_interactor(node, ns):
    org_node = node.find('{%s}organism' % ns)
    if None == org_node:
        return ''
    else:
        return org_node.get('ncbiTaxId')

def interactor_from_node(node, ns):
    "Get the primary and secondary references associated with this node"
    xref = node.find('{%s}xref' % ns)
    names = node.find('{%s}names' % ns)
    return Interactor(
label_for_interactor(names, ns),
ref_for_node(xref.find('{%s}primaryRef' % ns)),
[ r for r in secondary_refs_for_xref(xref, ns) ],
aliases_for_interactor(names, ns),
organism_for_interactor(node, ns)
    )

def add_tree_to_network(network, tree, ns='net:sf:psidev:mi'):
    """
    Adds an ElementTree representation of a PsiMi file to a network.
    """
    root = tree.getroot()

    #sometimes entry has namespace, sometimes not
    for entry in itertools.chain(
            root.findall('{%s}entry' % ns),
            root.findall('entry'),
    ):

        # look for interactor list for entry - is indexed by id
        interactors = {}
        for interactor in entry.findall('{%s}interactorList/{%s}interactor' % (ns, ns)):
            try:
                interactors[interactor.get('id')] = interactor_from_node(interactor, ns)
            except ValueError:
                pass


        # for each interaction
        for interaction in entry.findall('{%s}interactionList/{%s}interaction' % (ns, ns)):
            interactor_refs = []
            for participants in interaction.findall('{%s}participantList' % ns):
                for participant in itertools.chain(
                        participants.findall('{%s}participant' % ns),
                        participants.findall('{%s}proteinParticipant' % ns)
                ):
                    try:
                        # can we find the interactor under the correct name
                        interactor = participant.find('{%s}interactor' % ns)
                        if None == interactor:
                            # otherwise can we find the interactor under the old name
                            interactor = participant.find('{%s}proteinInteractor' % ns)
                        if None == interactor:
                            # look up the interactor by reference to earlier in the document
                            interactorRef = participant.find('{%s}interactorRef' % ns)
                            interactor = interactors[interactorRef.text]
                        else:
                            interactor = interactor_from_node(interactor, ns)

                        # yes we did find it - check it is right tax id
                        if interactor.taxid == '' or '10090' == interactor.taxid:
                            interactor_refs.append(interactor)

                    except ValueError:
                        pass
                    except KeyError:
                        #print sys.exc_info()[1]
                        pass

            # add all combinations of interactors to network
            for i, i1 in enumerate(interactor_refs):
                for i2 in interactor_refs[i+1:]:
                    network.add_interaction(i1, i2)


_mint_files = [
        'C:/Data/Protein-Protein/MINT/Mammalia-1.psi25.xml',
        'C:/Data/Protein-Protein/MINT/Mammalia-2.psi25.xml',
        'C:/Data/Protein-Protein/MINT/Mammalia-3.psi25.xml',
        'C:/Data/Protein-Protein/MINT/Mammalia-4.psi25.xml',
        'C:/Data/Protein-Protein/MINT/Mammalia-5.psi25.xml',
        'C:/Data/Protein-Protein/MINT/Mammalia-6.psi25.xml',
        'C:/Data/Protein-Protein/MINT/Mammalia-7.psi25.xml',
        'C:/Data/Protein-Protein/MINT/Mammalia-8.psi25.xml',
        'C:/Data/Protein-Protein/MINT/Mammalia-9.psi25.xml',
        'C:/Data/Protein-Protein/MINT/Mammalia-10.psi25.xml',
        'C:/Data/Protein-Protein/MINT/Mammalia-11.psi25.xml',
]

psimi_files = {
  'BIND-downloaded'             : [ 'C:/Data/Protein-Protein/Bind/taxid10090-downloaded.psi.xml' ],
  'BIND'                                                        : [ 'C:/Data/Protein-Protein/Bind/taxid10090.1.psi.xml' ],
  'MIPS'                                                        : [ 'C:/Data/Protein-Protein/MIPS/allppis.xml' ],
  'GRID'                                                        : [ 'C:/Data/Protein-Protein/GRID/BIOGRID-ORGANISM-Mus_musculus-2.0.33.psi25.xml' ],
  'DIP'                                                         : [ 'C:/Data/Protein-Protein/DIP/Mmusc20070707.mif25' ],
  'MINT'                                                        : _mint_files,
  'INTACT'                                              : [
        'C:/Data/Protein-Protein/INTACT/mouse_small-01.xml',
        'C:/Data/Protein-Protein/INTACT/mouse_small-02.xml',
        'C:/Data/Protein-Protein/INTACT/mouse_small-03.xml',
        'C:/Data/Protein-Protein/INTACT/mouse_small-04.xml',
        'C:/Data/Protein-Protein/INTACT/mouse_small-05.xml',
        'C:/Data/Protein-Protein/INTACT/mouse_small-06.xml',
        'C:/Data/Protein-Protein/INTACT/mouse_small-07.xml',
        'C:/Data/Protein-Protein/INTACT/mouse_small-08.xml',
        ],
}

def dbs():
    "The available PPI database names."
    return psimi_files.keys()


def build_network(db_name):
    "Build an interaction network from parsed psimi trees."
    network = InteractionNetwork()
    ns = 'DIP' == db_name and 'mi' or 'net:sf:psidev:mi'
    for psimi_file in psimi_files[db_name]:
        print '%s parsing: %s' % (db_name, psimi_file)
        tree = ET.parse(open(psimi_file))
        print '%s adding: %s' % (db_name, psimi_file)
        add_tree_to_network(network, tree, ns)
        del tree
    return network


class Networks(dict):
    "Cache the interaction networks"
    def __missing__(self, k):
        self[k] = build_network(k)
        return self[k]

try:
    networks
except NameError:
    networks = Networks()


class NamesForPssm(dict):
    def __missing__(self, k):
        names = set()
        if T.trans_data.matrix == k.db or T.trans_data.site == k.db:
            for fl in k.entry.factors:
                names.update(self[fl.link])
        if T.trans_data.factor == k.db:
            for s in k.entry.synonyms:
                names.add(s)
            if k.entry.gene:
                names.add(k.entry.gene.entry.name)
            for su in k.entry.subunits:
                names.update(self[su])
        self[k] = names
        return names

names_for_pssm = NamesForPssm()

class Transfac2Network(dict):
    "Maps transfac identifiers to interactors in a PPI network."
    def __init__(self, identifier_map, network):
        self.identifier_map = identifier_map
        self.network = network

    def __missing__(self, k):
        nodes = set()
        for name in names_for_pssm[k]:
            nodes.update(self.network.interactors_for_name(name))
        for link in self.identifier_map.links(k.as_db_ref()):
            nodes.update(self.network.interactors_for_ref(link))
        self[k] = nodes
        return nodes

    def coverage(self, pssm_accs):
        """
        Analyses how many pssms map to interactors and how many interactors are
        mapped to. Returns (num_interactors, num_pssms_mapped).
        """
        interactors = set()
        pssms_mapped = 0
        for acc in pssm_accs:
            interactors.update(self[acc])
            if len(self[acc]):
                pssms_mapped += 1
        return len(interactors), pssms_mapped

    def coverage_summary(self, pssm_accs):
        "A string summarising the result of self.coverage()."
        num_interactors, num_pssms_mapped = self.coverage(pssm_accs)
        return """%d / %d pssms mapped to interactors %.1f%%
%d / %d interactors are covered %.1f%%""" % (
num_pssms_mapped,
len(pssm_accs),
100.0 * num_pssms_mapped / len(pssm_accs),
num_interactors,
self.network.g.num_vertices(),
100.0 * num_interactors / self.network.g.num_vertices(),
)


if '__main__' == __name__:
    for db in psimi_files.keys():
        print 'Getting network: %s' % db
        network = networks[db]
        print network.summary()


    for db, network in networks.iteritems():
        print '%s typical primary ref db: %s' % (db, network.typical_primary_ref().db)

        import biopsy.transfac as T
        nodes_found = cookbook.DictOfLists()
        for m in T.Matrix.all():
            if m.acc in matrix_references:
                for matrix_ref in matrix_references[m.acc]:
                    for node in network.ref_2_nodes[matrix_ref]:
                        nodes_found[m.acc].append(node)
            if m.acc in matrix_names:
                for matrix_name in matrix_names[m.acc]:
                    for node in network.nodes_for_name(matrix_name):
                        nodes_found[m.acc].append(node)
        print '%s has references for %d / %d matrices' % (db, len(nodes_found), len(T.Matrix.all()))
