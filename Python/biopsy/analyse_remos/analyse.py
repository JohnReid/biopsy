#
# Copyright John Reid 2007
#

import biopsy.transfac, biopsy.pssm_paths
import os.path, re, webbrowser
from seqs import *
from paths import *

def analyse_remos(
        remo_fasta_files,
        centre_key_regex = "[Mm]ouse",
        threshold = 0.01,
        pssm_regex = ".",
        tag = 'all',
        open_svg = False,
        include_max_chain = False
):
    pssm_filter = biopsy.transfac.PssmFilter(
            name_regex_pattern = pssm_regex
    )
    pssm_accs = biopsy.get_transfac_pssm_accessions( pssm_filter )
    hit_map = { }

    for fasta_file in remo_fasta_files:
        remo = os.path.basename( os.path.splitext( fasta_file )[0] )
        print 'Looking at %s' % os.path.basename( fasta_file )
        seqs = biopsy.parse_fasta( fasta_file )
        sequences = convert_dict_seqs_for_phylo_analysis(
                seqs,
                centre_key_regex = re.compile( centre_key_regex )
        )
        #for s in sequences:
        #       hits = biopsy.HitVec()
        #       p_bind = biopsy.score_pssm_on_sequence(
        #               pssm_name = 'M00935',
        #               threshold = threshold,
        #               sequence = s,
        #               result = hits
        #       )
        #       if 0.0 == p_bind: print s
        #       print p_bind
        #continue
        for phylo_threshold in [
                # 0.01,
                0.001,
                # 0.0001
        ]:
            print 'Phylo threshold: %f' % phylo_threshold
            hits, max_chain, unadjusted_hits = biopsy.score_pssms_on_phylo_sequences(
                    pssm_accs,
                    sequences,
                    threshold = threshold,
                    phylo_threshold = phylo_threshold
            )
            if include_max_chain: hit_map[ remo ] = (hits, max_chain)
            else: hit_map[ remo ] = hits
            notes = 'Species: %s\nMatrix filter: %s\nThreshold: %f\nPhylo: %f' % (
                    ",".join( seqs.keys() ),
                    pssm_regex,
                    threshold,
                    phylo_threshold
            )
            # if None != max_chain: print 'Max chain length: %d' % len( max_chain )
            # else: print 'No max chain'
            build_svg_args = biopsy.BuildSvgArgs(
                    min_threshold = threshold,
                    file = '%s_%s.svg' % ( remo, tag ),
                    title = remo,
                    notes = notes,
                    open_file = open_svg
            )
            if len( hits ):
                biopsy.build_analysis_svg(
                        hits,
                        max_chain,
                        sequences[0],
                        args = build_svg_args
                )
            print 'Done remo %s with phylo_threshold = %f' % ( remo, phylo_threshold )
            print
            print
    return hit_map

def build_path_graph(
        name_1,
        name_2,
        hits_1,
        hits_2,
        max_length = 6,
        open_svg = True
):
    "Builds a graph of potential protein interaction between hits of 2 remos"
    name = '%s-%s' % (name_1, name_2)
    # build graph of protein chains between B and D
    g = path_graph(
            [ h.binder for h in hits_1 ],
            [ h.binder for h in hits_2 ],
            max_length = max_length
    )
    dot_file = '%s.dot' % name
    g.write_graphviz( dot_file )
    label = 'diamonds bind %s\\ngrays bind %s' % (name_1, name_2)
    svg_file = 'bridge_%s.svg' % name
    command = 'neato %s -Tsvg -o%s -Goverlap=scale -Glabelloc=t "-Glabel=%s"' % (
            dot_file,
            svg_file,
            label
    )
    if os.system( command ):
        raise RuntimeError( 'Could not execute: %s' % command)
    if open_svg: webbrowser.open( svg_file )

def build_all_path_graphs(
        hit_map,
        max_length = 6
):
    for k1, v1 in hit_map.iteritems():
        for k2, v2 in hit_map.iteritems():
            if k2 <= k1: continue
            build_path_graph(
                    k1,
                    k2,
                    v1,
                    v2,
                    max_length = max_length,
                    open_svg = False
            )
