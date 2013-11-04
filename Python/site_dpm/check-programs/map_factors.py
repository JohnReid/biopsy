#
# Copyright John Reid 2009
#

tp_indices = [
    1,
    6,
    14,
    17,
    24,
    28
]


def go_ids_for_genes(genes):
    import biopsy.identifiers.biomart as B
    for row in B.quick_query(
            dataset='mmusculus_gene_ensembl',
            attributes=[
                'ensembl_gene_id',
                #'go_cellular_component_id',
                'go_biological_process_id',
                #'go_molecular_function_id'
            ],
            filters=[('ensembl_gene_id', ','.join(genes))],
        ):
        yield row

def genes_for_go_id(go_id):
    import biopsy.identifiers.biomart as B
    for row in B.quick_query(
            dataset='mmusculus_gene_ensembl',
            attributes=['ensembl_gene_id'],
            filters=[('go', go_id)],
        ):
        yield row[0]


from boost.graph import Graph


class LabelledGraph(Graph):

    def __init__(self):
        Graph.__init__(self)
        self.labels = Graph.add_vertex_property(self, name='label', type='string')
        self.vertex_map = {}

    def add_labelled_vertex(self, label):
        if label in self.vertices:
            raise RuntimeError('Vertex for "%s" already in graph' % label)
        v = self.add_vertex()
        self.labels[v] = label
        self.vertex_map[label] = v
        return v

    def get_vertex(self, label):
        if label in self.vertex_map:
            return self.vertex_map[label]
        else:
            return self.add_labelled_vertex(label)


def create_graph(factors, pssm_map):
    """
    Create a bipartite graph representing which matrices map onto the factors.
    """
    import boost.graph as bgl
    g = LabelledGraph()
    vertices = {}
    for f in factors:
        for matrix, domain in pssm_map.iteritems():
            if f in domain:
                g.add_edge(g.get_vertex(matrix), g.get_vertex(f))
    return g


for tp_index in tp_indices:
    tp = transcriptional_programs[tp_index]
    print tp_index
    print tp.tp_factors
    g = create_graph(tp.tp_factors, pssm_map)
    graphviz_file = 'tp-%03d-factors.dot' % tp_index
    svg_file = 'tp-%03d-factors.svg' % tp_index
    g.write_graphviz(graphviz_file)
    os.system('dot %s -Tsvg -o %s' % (graphviz_file, svg_file))
