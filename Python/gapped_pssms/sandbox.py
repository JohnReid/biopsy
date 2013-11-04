

from gapped_pssms.sequence import convert_fasta_sequences
from gapped_pssms.background import k_mer_log_likelihoods, forward_backward_log_likelihoods
import cPickle, numpy as N, pylab as P

def reduce_sequence(sequence, desired_length):
    if len(sequence) < desired_length:
        raise RuntimeError('Cannot reduce sequence to longer length.')
    indexes = N.arange(desired_length) * len(sequence) / desired_length
    return sequence[indexes]

print reduce_sequence(N.arange(10), 10)
print reduce_sequence(N.arange(10), 9)
print reduce_sequence(N.arange(10), 6)
print reduce_sequence(N.arange(10), 5)
print reduce_sequence(N.arange(10), 3)
print reduce_sequence(N.arange(10), 2)
print reduce_sequence(N.arange(10), 1)

raise RuntimeError('Stopping')

K = 8
model_file = 'bg-model.pickle'
fasta = '/home/john/Data/GappedPssms/apr-2009/T00759trimRM.fa'

bg_model = cPickle.load(open(model_file))
sequences = convert_fasta_sequences(fasta)
converted_seqs = [bg_model.converter.to_order_n(s) for s in sequences]


def k_mer_log_likelihoods_new(alpha, c):
    result = N.empty(len(c)-K+1)
    alpha_sum = alpha.sum(axis=1)
    for i in xrange(len(c)-K+1):
        if 0 == i:
            result[i] = alpha_sum[K] / c[:K+1].prod()
        else:
            result[i] = alpha_sum[i+K-1] / alpha_sum[i-1] / c[i:i+K].prod()
    return N.log(result)


def calculate_k_mer_scores(bg_model, converted_seqs, K):
    result = list()
    for seq in converted_seqs:
        LL, alpha, beta, c = bg_model.forward_backward(seq)
        k_mer_LLs = k_mer_log_likelihoods(K=K, LL=LL, alpha=alpha, beta=beta, c=c)
        result.append(k_mer_LLs)
    return result

for i, seq in enumerate(converted_seqs[:8]):
    LL, alpha, c = bg_model.forward(seq)
    bg_k_mer_LLs = k_mer_log_likelihoods_new(alpha, c)
    P.plot(bg_k_mer_LLs, label='Seq %d' % i)

LL, alpha, beta, c = bg_model.forward_backward(converted_seqs[2])
bg_k_mer_LLs_new = k_mer_log_likelihoods_new(alpha, c)
bg_k_mer_LLs_old = k_mer_log_likelihoods(K=K, LL=LL, alpha=alpha, beta=beta, c=c)
P.figure()
P.plot(bg_k_mer_LLs_new, label='new')
P.plot(bg_k_mer_LLs_old, label='old')
P.legend()
P.show()
