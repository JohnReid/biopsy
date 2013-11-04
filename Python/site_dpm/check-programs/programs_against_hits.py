#
# Copyright John Reid 2009
#

for k, tp in enumerate(transcriptional_programs):
    logging.info('TP: %d; %d factors', k, len(tp.tp_factors))
    logging.info('TP: %d; %d targets', k, len(tp.tp_targets))
    if len(tp.tp_factors) > 1:
        for target in tp.tp_targets:
            factors_hitting_target = set(hit_counts[biopsy.transfac.DbRef.parse(target)].keys())
            logging.info(
                'TP: %d; target: %s; %3d of the %3d factors binding the target are in the program',
                k,
                target,
                len(factors_hitting_target.intersection(tp.tp_factors)),
                len(factors_hitting_target)
            )
