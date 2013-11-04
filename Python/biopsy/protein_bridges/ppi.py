#
# Copyright John Reid 2007
#

import biopsy, biopsy.analyse_remos, biopsy.transfac, re, cookbook, csv
T = biopsy.transfac
from itertools import chain
from remos_to_proteins import get_ensembl_2_mgi, get_entrez_2_mgi
from graph_as_svg import graph_as_svg
from cookbook import lru_cache
from . import protein_bridges
protein_bridge_builder = protein_bridges.protein_bridge_builder


def fasta_file(remo):
    return remo+'.fa'

def get_ppi_network(ppi_db_filename):
    import biopsy.psimi
    return biopsy.psimi.build_graph( open(ppi_db_filename) )

def handy_chain(iterable):
    for i in iterable:
        for x in i:
            yield x

class DictOfSets(dict):
    def __missing__(self, k):
        self.__setitem__(k, set())
        return self[k]

def _hit_location_cmp(l1, l2):
    "Compares 2 hit locations, first by start position, then by end"
    if l1.start() < l2.start():
        return -1
    if l2.start() < l1.start():
        return 1
    if l1.end() < l2.end():
        return -1
    if l2.end() < l1.end():
        return 1
    return 0

def sort_hit_locations(locations):
    "Return a sorted list of the locations"
    result = [l for l in locations]
    result.sort(cmp=_hit_location_cmp)
    return result

def absorb_hit_locations( locations ):
    "Takes a sorted iterable of hit locations and merges those that overlap"
    cur_loc = None
    for l in locations:
        if None == cur_loc:
            cur_loc = l
        else:
            if cur_loc.start() > l.start():
                raise RuntimeError("absorb_hit_locations(): requires input sorted by position")
            if cur_loc.end() >= l.start():
                cur_loc = biopsy.HitLocation( cur_loc.start(), l.end() - cur_loc.start(), False )
            else:
                yield cur_loc
                cur_loc = l
    if None != cur_loc:
        yield cur_loc

class Remo(object):
    def __init__(self, name, sequence_centre_key_regex, threshold, phylo_threshold, directory='.'):
        self.name = name
        self.seqs = biopsy.parse_fasta(fasta_file(name))
        self.converted_seqs = biopsy.analyse_remos.convert_dict_seqs_for_phylo_analysis(
                self.seqs,
                centre_key_regex=sequence_centre_key_regex)
        self.threshold = threshold
        self.phylo_threshold = phylo_threshold
        self.directory = directory

    def create_analysis(
            self,
            pssm_names,
            file_tag='',
            max_threshold=0.0
    ):
        print 'Analysing %s' % self.name
        hits, max_chain, unadjusted_hits = biopsy.score_pssms_on_phylo_sequences(
                pssm_names,
                self.converted_seqs,
                threshold = self.threshold,
                phylo_threshold = self.phylo_threshold
        )
        notes = 'Species:\n\t%s\nThreshold: %f\nPhylo: %f' % (
                "\n\t".join( self.seqs.keys() ),
                self.threshold,
                self.phylo_threshold
        )
        self.write_svg(
                filename = '%s/bifa-%s%s.svg' % (self.directory, self.name, file_tag),
                max_threshold = max_threshold,
                notes = notes,
                hits = hits,
                max_chain = max_chain
        )
        return hits, max_chain

    def write_svg(self, filename, hits, max_threshold=0.0, notes="", max_chain=biopsy.HitVec()):
        build_svg_args = biopsy.BuildSvgArgs(
                min_threshold = self.threshold,
                max_threshold = max_threshold,
                file = filename,
                title = self.name,
                notes = notes,
                open_file = False
        )
        if len( hits ):
            biopsy.build_analysis_svg(
                    hits,
                    max_chain,
                    self.converted_seqs[0],
                    args = build_svg_args
            )

    def analyse(
            self,
            pssm_names,
            use_max_chain = False
    ):
        "Run the BiFa analysis on the sequences."
        full, max_chain = self.create_analysis(pssm_names)
        self.analysis = use_max_chain and max_chain or full
        self.best_hit = max(h.p_binding for h in full)

    def interactors_for_hits(self, transfac_2_network):
        "Get a map from the hits to the interactors in the network."
        result = dict()
        for hit in self.analysis:
            acc = T.TableLink(hit.binder)
            interactors = transfac_2_network[T.TableLink(hit.binder)]
            if len(interactors):
                result[hit.binder] = interactors
        return result

    def map_hits(self, transfac_2_network, tag=''):
        """
        Build a map from the hits to the interactors in the network and
        rerun analysis based on those proteins.
        """
        self.hits_to_interactors = self.interactors_for_hits(transfac_2_network)
        self.interactors = set(handy_chain(self.hits_to_interactors.values()))

        # run analysis for just those pssms that mapped into the network.
        self.mapped_pssms = biopsy.SequenceVec()
        for pssm in self.hits_to_interactors.keys():
            self.mapped_pssms.append(pssm)
        self.filtered_analysis, tmp = self.create_analysis(
                pssm_names = self.mapped_pssms,
                file_tag='-just_ppi_proteins-%s' % tag,
                max_threshold=self.best_hit
        )

    def analyse_hits_for_interactors(self, interactors, tag):
        """
        Re-analyse the remo only for those pssms that are associated with
        the given interactors.
        """
        filtered_hits = self.hits_that_share_interactors(interactors)
        pssms = biopsy.SequenceVec()
        for hit in filtered_hits:
            pssms.append(hit.binder)
        self.filtered_analysis, tmp = self.create_analysis(
                pssm_names = pssms,
                file_tag=tag,
                max_threshold=self.best_hit
        )
        return pssms

    def write_hit_locations(self):
        "Create a text file of the hit locations."
        f = open('%s/%s_protein_locations.txt' % (self.directory, self.name), 'w')
        for interactor, pssms in self.interactors_to_pssms.iteritems():
            #if 'ETS1' == str(interactor):
            #       import IPython; IPython.Debugger.Pdb().set_trace()
            f.write('%20s: ' % str(interactor))
            f.write(', '.join(
                            '%d:%d' % (l.start(), l.end())
                            for l in absorb_hit_locations(
                              sort_hit_locations(
                                      hit.location
                                            for hit
                                            in self.analysis
                                            if hit.binder in pssms
                                    )
                            )
                    )
            )
            f.write('\n')

    def write_hit_to_proteins(self):
        "Write a text file linking the hits to the interactors."
        f = open('%s/%s_hits_to_proteins.txt' % (self.directory, self.name), 'w')
        self.interactors_to_pssms = DictOfSets()
        for pssm, interactors in self.hits_to_interactors.iteritems():
            interactor_names = [str(i) for i in interactors]
            pssm_name = biopsy.transfac.TableLink(pssm).name
            f.write('%20s: ' % pssm_name)
            f.write(', '.join(interactor_names))
            f.write('\n')
            for interactor in interactors:
                self.interactors_to_pssms[interactor].add(pssm)
        f.write('\n')
        for interactor, pssms in self.interactors_to_pssms.iteritems():
            f.write('%20s' % str(interactor))
            f.write(': ')
            f.write(', '.join(biopsy.transfac.TableLink(pssm).name for pssm in pssms))
            f.write('\n')

    def hits_that_share_interactors(self, interactors):
        "Return those hits that are associated with the given interactors."
        filtered_hits = biopsy.HitVec()
        for hit in self.analysis:
            if hit.binder in self.hits_to_interactors:
                if len(self.hits_to_interactors[hit.binder].intersection(interactors)):
                    filtered_hits.append(hit)
        return filtered_hits

class RemoAnalyser(object):
    '''
    Analyses remos and builds protein bridges between them
    '''
    def threshold(self, remo):
        return remo in self.override_thresholds and self. override_thresholds[remo] or self.default_threshold

    def phylo_threshold(self, remo):
        return remo in self.override_phylo and self.override_phylo[remo] or self.default_phylo_threshold

    def __init__(self, bridge_sets, ppi_network, transfac_2_network, tag='', pssm_names=None, directory='.'):
        self.bridge_sets = bridge_sets
        self.ppi_network = ppi_network
        self.transfac_2_network = transfac_2_network
        self.directory = directory
        self.tag = tag
        if None == pssm_names:
            self.pssm_names = biopsy.get_transfac_pssm_accessions(
              biopsy.get_default_transfac_pssm_filter()
            )
        else:
            self.pssm_names = biopsy.SequenceVec()
            for name in pssm_names:
                self.pssm_names.append(name)
        self.sequence_centre_key_regex = re.compile("[Mm]ouse")
        self.default_threshold = 0.01
        self.override_thresholds = {}
        self.default_phylo_threshold = 0.001
        self.override_phylo = {}
        self.use_max_chain = False
        remo_universe = set(chain(*list(names for s, names in bridge_sets.iteritems())))
        self.remos = dict(
                (
                        name,
                        Remo(
                         name,
                         self.sequence_centre_key_regex,
                         self.threshold(name),
                         self.phylo_threshold(name),
                         directory = self.directory
                        )
                ) for name in remo_universe
        )

    def analyse(self, use_max_chain = False):
        "Analyse the remos."
        for name, remo in self.remos.iteritems():
            remo.analyse(self.pssm_names, use_max_chain=use_max_chain)
            remo.map_hits(self.transfac_2_network, tag=self.tag)
            remo.write_hit_to_proteins()
            remo.write_hit_locations()

    def file_tag(self, bridge_set, max_path_length):
        return '%s-length_%d-%s' % (bridge_set, max_path_length, self.tag)

    def build_bridges(self, max_path_length):
        bridges = dict()
        for bridge_set, bridge_remos in self.bridge_sets.iteritems():
            bridge = protein_bridge_builder(
             self.ppi_network,
             dict((remo, self.remos[remo].interactors) for remo in bridge_remos)
            )
            print '%s: Building protein bridges' % bridge_set
            bridge.build(max_path_length = max_path_length)

            # work out which hits are in the bridge
            print '%s: Finding hits that could be part of potential protein bridges' % bridge_set
            for remo_name in bridge_remos:
                self.remos[remo_name].analyse_hits_for_interactors(
                        bridge.interactors_for_remo(remo_name),
                        tag = '-for_bridge-%s' % self.file_tag(bridge_set, max_path_length),
                )

            #print 'Writing remo-remo paths'
            #bridge.write_remo_remo_paths('bridge_paths-%s.txt' % self.file_tag(bridge_set, max_path_length), max_length=max_path_length)

            print '%s: Writing graph as svg' % bridge_set
            graph_as_svg(
                    bridge.g,
                    '%s/bridges-%s' % (self.directory, self.file_tag(bridge_set, max_path_length)),
                    '-Elen=2.0',
                    '-Gscale=overlap'
            )

            print '%s: Writing graph as dot' % bridge_set
            bridge.g.write_graphviz('%s/bridges-%s.dot' % (self.directory, self.file_tag(bridge_set, max_path_length)))

            #print 'Writing graph as GraphML'
            #bridge.g.write_graphml('%s/bridges-%s.graphml' % (self.directory, self.file_tag(bridge_set, max_path_length)))

            bridges[bridge_set] = bridge
        return bridges
