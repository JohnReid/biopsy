#
# Copyright John Reid 2011
#


"""
Analyse the UCSC promoters for binding sites.
"""

import biopsy.data.ucsc; reload(biopsy.data.ucsc)
import biopsy.data.ucsc as UCSC
import biopsy.transfac as T
import logging
from IPython.kernel import client
from itertools import imap


def parse_record_id(_id):
    fields = _id.split('_')
    assert 7 == len(fields)
    nm_id = int(fields[1])
    assert fields[2] == 'up'
    assert fields[3] == '5000'
    chr = fields[4]
    pos = int(fields[5])
    tag = fields[6]
    return nm_id, chr, pos, tag


def make_sequence_vec(l):
    import biopsy
    result = biopsy.SequenceVec()
    result.extend(l)
    return result


def get_TRANSFAC_vertebrate_pssms():
    return [str(m.acc) for m in T.Matrix.all() if m.name.startswith('V')]


def calculate_p_binding_from_hits(hits):
    from collections import defaultdict
    p_not_binding = defaultdict(lambda: 1.)
    for hit in hits:
        p_not_binding[hit.binder] *= (1.-hit.p_binding)
    return dict((binder, 1.-p) for binder, p in p_not_binding.iteritems())
    

def score_pssms_on_record(pssms, record):
    import biopsy
    promoter_hits = biopsy.score_pssms_on_sequence(make_sequence_vec(pssms), record.seq.tostring())
    promoter_p_binding = calculate_p_binding_from_hits(promoter_hits)
    return promoter_p_binding

def factor_info(f):
    return "%s;%s" % (f.name, ";".join(species.replace(',', ' -').replace(';', ' -').replace(':', ' -') for species in  f.species))
    

logging.basicConfig(level=logging.INFO)    

pssms = get_TRANSFAC_vertebrate_pssms()
logging.info('Have %d PSSMs from TRANSFAC', len(pssms))

#
# Distribute function to engines
#    
mec = client.MultiEngineClient()
mec.push_function(
    dict(
        calculate_p_binding_from_hits=calculate_p_binding_from_hits,
        make_sequence_vec=make_sequence_vec,
    )
)

#
# pass tasks to engines
#
logging.info('Passing tasks to engines')
tc = client.TaskClient()
task_ids = []
task_args = dict()
promoters = dict()
positions_analysed = set()
skipped = set()
for record in UCSC.upstream('hg19', length=5000):
    
    #
    # Parse record id
    #
    try:
        nm_id, chr, pos, tag = parse_record_id(record.id)
    except:
        logging.warning('Skipping unparseable record ID: %s', record.id)
	skipped.add(record.id)
    else:
        logging.debug('Creating task for NM_%-16d %6s %10d %s %r', nm_id, chr, pos, tag, record.seq)
        promoters[nm_id] = chr, pos, tag, record
        pos_key = chr, pos, tag

        #
        # Check we have not analysed this genomic position
        #
        if pos_key not in positions_analysed:
            
            #
            # Pass task to engine
            #
            args = (map(str, pssms), record)
            task = client.MapTask(score_pssms_on_record, args)
            task_id = tc.run(task, block=False)
            task_ids.append(task_id)
            task_args[task_id] = args
            positions_analysed.add(pos_key)
            #if len(task_ids) > 4: break # for debugging
logging.info('Skipped %d unparseable records', len(skipped))
        

#
# Write info about the PSSMs
#
logging.info('Writing PSSMs info.')
pssm_file = open('pssms.csv', 'w')
pssm_file.write("chr:pos:orientation,%s\n" % ",".join(imap(str, pssms))) # write header
for pssm in pssms:
    m = T.Matrix(pssm)
    pssm_file.write(
        "%s,%s,%s\n" % (
            pssm,
            m.name,
            ":".join(imap(factor_info, m.factors))
        )
    )
pssm_file.close()


#
# Write info about each promoter
#
logging.info('Writing transcripts info.')
transcripts_file = open('transcripts.csv', 'w')
transcripts_file.write("RefSeq,Chr,Pos,Orientation\n") # write header
for nm_id, (chr, pos, tag, record) in promoters.iteritems():
    transcripts_file.write('NM_%d,%s,%d,%s\n' % (nm_id, chr, pos, tag))
transcripts_file.close()



#
# Get results from engines
#
logging.info('Blocking on %d results...', len(task_ids))
genome_p_binding = dict()
for task_id in task_ids:
    pssms, record = task_args[task_id]
    nm_id, chr, pos, tag = parse_record_id(record.id)
    pos_key = chr, pos, tag
    logging.debug('Blocking on result for NM_%-16d %6s %10d %s %r', nm_id, chr, pos, tag, record.seq)
    promoter_p_binding = tc.get_task_result(task_id, block=True)
    genome_p_binding[pos_key] = promoter_p_binding
    
logging.info('Analysed %d promoters', len(promoters))


#
# Write results out
#
logging.info('Writing results.')
tf_nm_matrix = open('tf-nm-matrix.csv', 'w')
tf_nm_matrix.write("chr:pos:orientation,%s\n" % ",".join(imap(str, pssms))) # write header
for (chr, pos, tag), p_binding_map in genome_p_binding.iteritems():
    tf_nm_matrix.write(
        "%s:%d:%s,%s\n" % (
        chr, pos, tag, ",".join(str(p_binding_map.get(pssm, 0)) for pssm in pssms))
    )
tf_nm_matrix.close()


