#
# Copyright John Reid 2008, 2009
#

"""
Code to test various things in a sandbox environment.
"""

import sys, csv, time
from gene_set_enrichment import kegg_and_inoh_gene_sets

if '__main__' == __name__:
    for i in xrange(100):
        print i
        time.sleep(1)

    raise ''



    def genes_in_term(go_analysis, term):
        return set(r.genesInTerm(go_analysis.go_data, term)[term])

    pathways = dict(kegg_and_inoh_gene_sets())

    for k, go_term, pathway_id in (
      ( 1, 'GO:0009059', 'mmu04620'),
      ( 6, 'GO:0001708', 'mmu04360'),
      (14, 'GO:0019221', 'mmu04630'),
      (24, 'GO:0009954', 'mmu04340'),
      (28, 'GO:0007178', 'mmu04350'),
    ):
        print k, go_term, pathway_id
        tp = transcriptional_programs[k]
        term_genes = genes_in_term(tp.factors_go_analysis['BP'], go_term)
        pathway = pathways[pathway_id]

        print 'In GO term and TP: %d' % len(term_genes.intersection(tp.tp_factors))
        print 'In pathway and TP: %d' % len(pathway.intersection(tp.tp_factors))
        print 'In all 3:          %d' % len(term_genes.intersection(tp.tp_factors).intersection(pathway))


    factors_30 = transcriptional_programs[30].tp_factors
    pathways_for_30 = [
      'mmu04520',
      'mmu04916',
      'mmu05213',
      'mmu05216',
      'mmu05217',
    ]
    for p1 in pathways_for_30:
        for p2 in pathways_for_30:
            if p1 == p2:
                continue
            print pathways[p1].intersection(factors_30).intersection(pathways[p2])


    raise



    def threshold_results(results, threshold):
        if not results.classic:
            return results
        return r.subset(results, [c < threshold for c in imap(float, results.classic)])



    def print_latex_table(k, size, ontology, results):
        for go_id, go_desc, significant, annotated, p_value in zip(
          r['$'](results, 'GO.ID'),
          results.Term,
          results.Significant,
          results.Annotated,
          results.classic
        ):
            print ' & '.join((str(k), str(size), go_id, ontology, go_desc, str(significant), str(annotated), '%.1f' % math.log10(float(p_value)))), ' \\\\'



    def print_go_analyses(self, f, k, threshold):
        for ontology, go_analysis in self.factors_go_analysis.iteritems():
            if None != go_analysis:
                self.print_go_analysis(f, k, 'factors', len(tp.tp_factors), ontology, threshold, go_analysis)
        for ontology, go_analysis in self.targets_go_analysis.iteritems():
            if None != go_analysis:
                self.print_go_analysis(f, k, 'targets', len(tp.tp_targets), ontology, threshold, go_analysis)

    for set_type, tp_fn in (
      ('Factors', lambda tp: (len(tp.tp_factors), tp.factors_go_analysis)),
      ('Targets', lambda tp: (len(tp.tp_targets), tp.targets_go_analysis)),
    ):
        print '%% %s' % set_type
        tp = transcriptional_programs[0]
        for k, tp in enumerate(transcriptional_programs):
            size, go_analyses = tp_fn(tp)
            for ontology, go_analysis in go_analyses.iteritems():
                if go_analysis:
                    print_latex_table(k, size, ontology, threshold_results(go_analysis.results, 1e-4))
        print


    raise



    from SOAPpy import WSDL

    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)

    results = serv.get_genes_by_pathway('path:eco00020')
    print results

    raise

    from soaplib.client import make_service_client
    from soaplib.wsgi_soap import SimpleWSGISoapApp
    from soaplib.service import soapmethod
    from soaplib.serializers.primitive import String, Integer, Array, Float

    KEGG_wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    class KeggService(SimpleWSGISoapApp):
        __tns__ = 'SOAP/KEGG'

        @soapmethod(
          String,
          _returns=String,
          _inMessage='bconvRequest',
          _outMessage='bconvResponse',
          _outVariableName='return'
        )
        def bconv(self, id):
            pass

    client = make_service_client(KEGG_wsdl, KeggService())
    print client.bconv("ncbi-gi:10047086 ncbi-gi:10047090 ncbi-geneid:14751")



    raise




    tp = transcriptional_programs[27]
    #for g in tp.tp_targets:
    #  print ensembl_names[g]

    def gene_name(gene):
        return ensembl_names.get(str(gene))

    def gene_names(genes):
        return ifilter(None, imap(gene_name, genes))

    gene_universe = set(gene_names(tp.genes))
    tp_genes = set(gene_names(tp.tp_targets))

    def load_conversions(f):
        for row in csv.reader(f, delimiter='\t'):
            gene, ensembl, unigene, entrez_gene, genbank, clone, ens_chr, start, end, strand, tmp = row
            yield gene, unigene

    tgf_beta_genes = set()
    tgf_beta_clones = set()
    for line in open('tgf-beta-targets.txt'):
        spl = line.strip().split()
        clone, symbol, desc = spl[0], spl[1], spl[2:]
        tgf_beta_genes.add(symbol)
        tgf_beta_clones.add(clone)
        #print symbol

    def unigene_for_genes(conversions, genes):
        return set(ifilter(None, imap(conversions.get, genes)))

    conversions = dict(
      chain(
        load_conversions(open('tgf-beta-conversions.txt')),
        load_conversions(open('tp-target-conversions.txt')),
        load_conversions(open('universe-conversions.txt'))
      )
    )
    tgf_beta_unigene = unigene_for_genes(conversions, tgf_beta_genes)
    tp_target_unigene = unigene_for_genes(conversions, tp_genes)
    universe_unigene = unigene_for_genes(conversions, gene_universe)

    raise


    def r_to_str(robj):
        "Returns an R object in a representation as a list of strings."
        from rpy import r
        from tempfile import mktemp
        tmpfile = mktemp()
        #logging.info('Tmpfile: %s' % tmpfile)
        try:
            r.assign('tmpobj', robj)
            r('save(tmpobj, file="%s", ascii=TRUE)' % tmpfile)
            return open(tmpfile).read()
        finally:
            if os.access(tmpfile, os.R_OK):
                os.remove(tmpfile)

    def r_from_str(s):
        "Returns an R object in a representation as a list of strings."
        from rpy import r, with_mode, NO_CONVERSION
        from tempfile import mktemp
        tmpfile = mktemp()
        #logging.info('Tmpfile: %s' % tmpfile)
        try:
            open(tmpfile, 'w').write(s)
            names = with_mode(NO_CONVERSION, lambda : r.load(file=tmpfile))()
        finally:
            if os.access(tmpfile, os.R_OK):
                os.remove(tmpfile)

    raise

    from soaplib.client import make_service_client
    client = make_service_client('http://soap.genome.jp/KEGG.wsdl', KeggService())


    def threshold_results(results, threshold):
        if not results.classic:
            return results
        return r.subset(results, [c < threshold for c in imap(float, results.classic)])

    f = open('go_mf_analysis.txt', 'w')
    #f = sys.stdout
    my_threshold = 1.e-4
    num_hits = 0
    for k in xrange(len(transcriptional_programs)):
        tp = transcriptional_programs[k]
        if tp.factors_go_analysis:
            thresholded = threshold_results(tp.factors_go_analysis.results, my_threshold)
            if thresholded.classic:
                num_hits += 1
                print >> f, '************ TP: %d **************; # factors=%d' % (k, len(tp.tp_factors))
                print >> f, str(thresholded)
        if tp.targets_go_analysis:
            thresholded = threshold_results(tp.targets_go_analysis.results, my_threshold)
            if thresholded.classic:
                num_hits += 1
                print >> f, '************ TP: %d **************; # targets=%d' % (k, len(tp.tp_targets))
                print >> f, str(thresholded)
        #print >> f

    print >> f, '# results below threshold: %d' % num_hits

    f.close()

    for k, tp in enumerate(transcriptional_programs):
        print '************ TP: %d **************; # targets=%d' % (k, len(tp.tp_targets))
        if tp.targets_go_analysis:
            print str(tp.targets_go_analysis.results)



    K=17
    theta=dpm.exp_theta()
    theta=theta[:,:K]
    bg_topic_dist=theta.sum(axis=0)
    bg_topic_dist/=bg_topic_dist.sum()
    for d in xrange(6):
        figure()
        title('Document: %d' % d)
        subplot(3,1,1)
        bar(xrange(K), bg_topic_dist)
        subplot(3,1,2)
        bar(xrange(K), theta[d])
        subplot(3,1,3)
        bar(xrange(K), theta[d]/bg_topic_dist)
