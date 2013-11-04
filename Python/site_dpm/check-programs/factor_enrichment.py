
def display_tp(dpm, k):
    phi_ratio = dpm.exp_phi() / dpm.exp_Phi()

    clf()

    subplot(311)
    plot(phi_ratio[k])
    title('phi/Phi')

    subplot(312)
    plot(dpm.E_n_kw[k])
    title('E(n_kw)')

    subplot(313)
    sorted = phi_ratio[k].copy()
    sorted.sort()
    plot(sorted[::-1])

    print len(stats.words_for_topic(k))
