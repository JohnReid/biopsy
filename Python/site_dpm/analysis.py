#
# Copyright John Reid 2008, 2009
#

"""
Code to analyse sequences in preparation for the HDPM.
"""


from shared import *
import biopsy.analyse_remos.consolidate_hits as CH
import biopsy.analyse_remos.expectations as EX


from boost.graph import Graph

class LabelledGraph(Graph):
    """
    A boost.graph.Graph that indexes vertices by strings.
    """

    def __init__(self, label_property_name='label', label_property_type='string'):
        Graph.__init__(self)
        self.labels = Graph.add_vertex_property(self, name=label_property_name, type=label_property_type)
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



def find_all_aligned_sequences_in_refs(remome, refs):
    "Return a set of all aligned sequence sets that have an id in ref"
    aligned_sequences = remome.get_aligned_sequences()
    result = set()
    for aligned_seq in aligned_sequences:
        for seq_id in aligned_seq.get_sequence_ids():
            gene_id = seq_id.gene_id
            ref = biopsy.DbRef(biopsy.db.ensembl, gene_id.prefix, gene_id.num)
            #import IPython; IPython.Debugger.Pdb().set_trace()
            if ref in refs:
                result.add(aligned_seq)
                break
    logging.info('Matched %d aligned sequence sets to genes', len(result))
    return result



def get_remo_sequences(remo, centre_id, masked=True):
    """
    @return: a sequence of sequences for the remo. The first sequence is the centre sequence
    """
    sequences = biopsy.SequenceVec()
    sequences.append(remo.get_sequence_for(centre_id, masked))
    for seq_id in remo.get_sequence_ids():
        if seq_id != centre_id:
            yield remo.get_sequence_for(seq_id, masked)



def get_sequences_for_sequence_sets(remome, aligned_sequences, masked=True):
    sequence_dict = {}
    logging.info('Using %s sequences', masked and 'masked' or 'unmasked')
    for aligned_seq in aligned_sequences:
        remos = remome.get_remos_for(aligned_seq)
        gene_id = aligned_seq.centre_sequence.gene_id
        gene = biopsy.DbRef(biopsy.db.ensembl, gene_id.prefix, gene_id.num)
        if gene not in sequence_dict:
            sequence_dict[gene] = list()
        for remo in remos:
            sequences = biopsy.SequenceVec()
            sequences.extend(get_remo_sequences(remo, aligned_seq.centre_sequence, masked=masked))
            sequence_dict[gene].append(sequences)
    return sequence_dict


@output_cached_method('pssm-ensembl-map')
def get_pssm_to_ensembl_map():
    """
    @return: A map from PSSMs to sets of ENSMUSGs.
    """
    import biopsy.transfac as T
    pssm_filter = T.PssmFilter()
    pssm_map = cookbook.DictOfSets()
    for acc in biopsy.get_transfac_pssm_accessions(pssm_filter):
        for f in biopsy.get_factors_for_pssm(acc):
            f = T.Factor(f)
            if f.gene:
                for ref in f.gene.entry.db_refs:
                    if ref.db == biopsy.db.ensembl and ref.table == 'ENSMUSG':
                        pssm_map[acc].add(str(ref))
    logging.info('Found %d PSSMs that map to Ensembl mouse genes', len(pssm_map))

    return pssm_map


@global_cached_method('human-mouse-orthologs')
def get_human_mouse_orthologs():
    from biopsy.identifiers.biomart import quick_query
    logging.info('Getting human mouse orthologs from Ensembl')
    result = dict(quick_query(dataset='hsapiens_gene_ensembl', attributes=['ensembl_gene_id', 'mouse_ensembl_gene']))
    logging.info('Mapped %d human genes to mouse', len(result))
    return result


@global_cached_method('rat-mouse-orthologs')
def get_rat_mouse_orthologs():
    from biopsy.identifiers.biomart import quick_query
    logging.info('Getting rat mouse orthologs from Ensembl')
    result = dict(quick_query(dataset='rnorvegicus_gene_ensembl', attributes=['ensembl_gene_id', 'mouse_ensembl_gene']))
    logging.info('Mapped %d rat genes to mouse', len(result))
    return result


@output_cached_method('pssm-ensembl-map-min-range')
def get_pssm_to_ensembl_map_min_range():
    """
    @return: A map from PSSMs to single ENSMUSGs.

    Uses a greedy heuristic to minimise the range of the map.
    """
    import biopsy.transfac as T

    human_mouse_orthologs = get_human_mouse_orthologs()
    rat_mouse_orthologs = get_rat_mouse_orthologs()

    # Build a graph containing possible mappings from PSSMs to ENSMUSGs
    graph = LabelledGraph()
    pssms = set()
    genes = set()
    pssm_filter = T.PssmFilter()

    def add_mapping(acc, ref):
        ref = str(ref)
        pssms.add(acc)
        genes.add(ref)
        graph.add_edge(graph.get_vertex(acc), graph.get_vertex(ref))

    def deal_with_pssm(acc, use_orthologs=False):
        for f in imap(T.Factor, biopsy.get_factors_for_pssm(acc)):
            if f.gene:
                for ref in f.gene.entry.db_refs:
                    if ref.db == biopsy.db.ensembl:
                        if ref.table == 'ENSMUSG':
                            add_mapping(acc, ref)
                        elif use_orthologs:
                            if ref.table == 'ENSRNOG':
                                ortholog = rat_mouse_orthologs[str(ref)]
                                if ortholog:
                                    add_mapping(acc, ortholog)
                            elif ref.table == 'ENSG':
                                ortholog = human_mouse_orthologs[str(ref)]
                                if ortholog:
                                    add_mapping(acc, ortholog)

    for acc in biopsy.get_transfac_pssm_accessions(pssm_filter):
        deal_with_pssm(acc, use_orthologs=False)
        # if we didn't map the PSSM directly onto a mouse gene trying going through rat/human orthologs
        if acc not in pssms:
            deal_with_pssm(acc, use_orthologs=True)
    logging.info('%d PSSMs map onto to %d Ensembl mouse genes', len(pssms), len(genes))

    # Use greedy heuristic to build Many-to-1 mapping
    sorted_by_degree = [(graph.in_degree(graph.get_vertex(g)), g) for g in genes]
    sorted_by_degree.sort(reverse=True)
    pssm_map = dict()
    for degree, g in sorted_by_degree:
        for pssm_v in graph.adjacent_vertices(graph.vertex_map[g]):
            pssm = graph.labels[pssm_v]
            if pssm not in pssm_map: # if we don't already have a map for this pssm, then use this mapping
                pssm_map[pssm] = g

    del pssm_map['M00158']

    # Calculate how many Ensembl genes are in range now
    logging.info('%d PSSMs map onto to %d Ensembl mouse genes (minimum range)', len(pssms), len(set(pssm_map.values())))

    return pssm_map



@cookbook.cache_decorator.cachedmethod
def get_remome(threshold=100):
    "Get the remome of the given threshold."
    remome_file = os.path.join(biopsy.get_data_dir(), 'ReMos', '%d' % threshold, '%d.filtered' % threshold)
    logging.info('Loading remome %d from %s', threshold, remome_file)
    remome = biopsy.Remome.load(remome_file)
    logging.info('Have %d aligned sequence sets in remome', len(remome.get_aligned_sequences()))
    return remome


def get_aligned_sequences(remome_threshold=100):
    return set(a for a in get_remome(remome_threshold).get_aligned_sequences())


@output_cached_method('sequence-dict')
def get_sequence_dict(remome_threshold=100, masked=True):
    aligned_sequences = get_aligned_sequences(remome_threshold)
    return get_sequences_for_sequence_sets(get_remome(remome_threshold), aligned_sequences, masked)


@output_cached_method('remome-analysis')
def analyse_remome(remome_threshold=100, masked=True, threshold=0.01):
    return analyse_remome_sequences(get_sequence_dict(remome_threshold, masked), threshold=threshold)


def analyse_remome_sequences(sequence_dict, threshold=0.01):
    """
    Analyse the sequences.

    @return: A dict mapping gene ids to sequences of (hits, max_chain) tuples.
    """
    # get the pssm accessions we will use
    pssm_accs = biopsy.SequenceVec()
    pssm_accs.extend(get_pssm_to_ensembl_map_min_range().keys())

    result = cookbook.DictOfLists()

    for gene, sequences in sequence_dict.iteritems():
        logging.info('Analysing %d remos for %s', len(sequences), gene)
        hit_counts = None

        for seqs in sequences:
            #import IPython; IPython.Debugger.Pdb().set_trace()
            hits, max_chain, unadjusted_hits = biopsy.score_pssms_on_phylo_sequences(
              pssm_accs,
              seqs,
              threshold = threshold,
              phylo_threshold = threshold
            )
            result[gene].append((hits, max_chain))

    return result



def consolidate_gene_analysis(gene_analysis):
    """
    Take the analysis for one gene and consolidate its hits. The analysis should be an iterable over
    (hits, max_chain) tuples. The hits are filtered by the hit_threshold. The maximal chain is left
    as is.
    """
    hit_counts = None
    max_chain_counts = None
    pssm_map = get_pssm_to_ensembl_map_min_range()

    for hits, max_chain in gene_analysis:

        factor_mapped_hits = map_binders(hits, pssm_map)
        #consolidated_hits = CH.consolidate_hits(factor_mapped_hits)
        consolidated_hits = CH.maximal_chain_hits(factor_mapped_hits)
        if consolidated_hits:
            hit_counts = EX.num_hits_per_binder(consolidated_hits, hit_counts)

        if None != max_chain:
            factor_mapped_hits = map_binders(max_chain, pssm_map)
            #consolidated_hits = CH.consolidate_hits(factor_mapped_hits)
            consolidated_hits = CH.maximal_chain_hits(factor_mapped_hits)
            if consolidated_hits:
                max_chain_counts = EX.num_hits_per_binder(consolidated_hits, max_chain_counts)

    return hit_counts, max_chain_counts


@output_cached_method('consolidated-hits')
def consolidate_remome_analysis(
    remome_threshold=100,
    masked=True,
    use_max_chain=False,
    analysis_threshold=.01
):
    """
    Consolidates analysis for entire remome.
    """
    remome_analysis = analyse_remome(remome_threshold, masked, analysis_threshold)
    return dict(
      (gene, consolidate_gene_analysis(gene_analysis))
      for gene, gene_analysis in remome_analysis.iteritems()
    )



@output_cached_method('hit-counts')
def get_hit_counts_for_remome_analysis(
    remome_threshold=100,
    masked=True,
    use_max_chain=False,
    analysis_threshold=.01,
):
    """
    Get the consolidated hit counts for the remome analysis.

    @arg remome_threshold: The threshold used to build the remome.
    @arg masked: Were the sequences masked before BiFa analysis?
    @arg use_max_chain: Use the counts from the maximal chain?
    @arg hit_threshold: If not using the maximal chain filter the hits by this threshold.
    """
    consolidated_analysis = consolidate_remome_analysis(
      remome_threshold,
      masked,
      use_max_chain,
      analysis_threshold
    )
    if use_max_chain:
        return dict((gene, counts[1]) for gene, counts in consolidated_analysis.iteritems() if counts[1])
    else:
        return dict((gene, counts[0]) for gene, counts in consolidated_analysis.iteritems() if counts[0])


class AnalyseSequence(object):
    def __init__(self, pssms, threshold=.01):
        logging.info('Creating sequence analyser with threshold of %f for %d pssms', threshold, len(pssms))
        self.pssms = biopsy.SequenceVec()
        self.pssms.extend(pssms)
        self.threshold = threshold
    def __call__(self, sequence):
        import biopsy
        import biopsy.analyse_remos.consolidate_hits as CH
        hits = biopsy.score_pssms_on_sequence(self.pssms, sequence.strip('Nn'), self.threshold)
        return CH.maximal_chain_hits(hits)


def get_sequence_analyser():
    pssm_to_ensembl_map = get_pssm_to_ensembl_map_min_range()
    return AnalyseSequence(pssm_to_ensembl_map.keys(), threshold=options.analysis_threshold)


def get_analysis_union(analysis):
    "Takes analysis and returns one sequence containing all the hits."
    result = list()
    for gene, hit_seqs in analysis.iteritems():
        for hits in hit_seqs:
            result.extend(hits)
    return result


def analyse_union(analysis):
    "Take an analysis and return various statistics."
    def p_binding(hit):
        return hit.p_binding
    union = get_analysis_union(analysis)
    union.sort(key=p_binding, reverse=True)

    # find the highest scoring hit of each binder
    already_found = set()
    best_scores = list()
    for hit in union:
        if hit.binder not in already_found:
            best_scores.append((hit.binder, hit.p_binding))
            already_found.add(hit.binder)

    return union, best_scores
