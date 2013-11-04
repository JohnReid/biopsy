#
# Copyright John Reid 2008
#

"""
Code to implement finding and erasing sites in sequences.
"""

def yield_sites_from_states(states, bg_states):
    """
    Yields (start, length) tuples, one per site in the state sequence.
    """
    in_site = False
    site_start = -1
    for i, q in enumerate(states):
        is_bg = q in bg_states
        # are we changing state in or out of site?
        if is_bg == in_site:
            # yes changing state...
            if is_bg:
                # changing state out of site
                site_len = (i-site_start)
                yield site_start, site_len
                in_site = False
            else:
                # changing state into site
                site_start = i
                in_site = True
    if in_site:
        site_len = len(states) - site_start
        yield site_start, site_len



def yield_sites_in_sequence(
  model,
  sequence,
  bg_states
):
    """
    Yields sites in the sequence as tuples (start, length).

    @arg bg_states: Set of states that make up the background.
    """
    log_P_star, q_star = model.viterbi(sequence)
    for site in yield_sites_from_states(q_star, bg_states):
        yield site


def remove_sites(sites, sequence, unknown_output):
    """
    Sets as unknown all those bases in site sequence.

    @arg sites: A sequence of sites represented as (start, length) tuples.
    @arg sequence: The sequence
    @arg unknown_output: The output value that represents unknown data.
    """
    for start, length in sites:
        sequence[start:start+length] = unknown_output


def remove_predicted_sites(
  model,
  sequence,
  bg_states
):
    """
    Set as unknown all those bases in sites predicted by Viterbi algorithm.
    """
    remove_sites(yield_sites_in_sequence(model, sequence, bg_states), sequence, model.M)


def sites_in_sequences(model, bg_states, sequences):
    """
    @return: (# sites, # sequences with at least one site)
    """
    # find sites...
    sites = [list(yield_sites_in_sequence(model, s, bg_states)) for s in sequences]
    return sum(len(s) for s in sites), sum(len(s) > 0 for s in sites)
