#
# Copyright John Reid 2008
#

"""
Code to test chow liu using conditioned and unconditioned data.
"""

import numpy.random as R
import chow_liu as CL
from itertools import chain, repeat


p_is_conditioned = .5
dirichlet_prior_strength = 1.
num_data = 20
K = 4 # num. bases
num_examples = 100

def create_multinomial():
    return R.dirichlet([dirichlet_prior_strength] * K)

def convert_multinomial_sample(sample):
    return list(chain(*(repeat(x, s) for x, s in enumerate(sample))))

def shuffle(l):
    R.shuffle(l)
    return l


def generate_unconditioned_example():
    return list(
      zip(
        convert_multinomial_sample(R.multinomial(num_data, create_multinomial())),
        shuffle(convert_multinomial_sample(R.multinomial(num_data, create_multinomial())))
      )
    )

def generate_conditioned_example():
    return list(
      chain(
        *(
          [(x, y) for y in convert_multinomial_sample(R.multinomial(num, create_multinomial()))]
          for x, num in enumerate(R.multinomial(num_data, create_multinomial()))
        )
      )
    )

def generate_example():
    """
    Generate one example of num_data points
    """
    is_conditioned = R.random() < p_is_conditioned
    return is_conditioned, is_conditioned and generate_conditioned_example() or generate_unconditioned_example()

conditioned_LL_ratio = 0.
unconditioned_LL_ratio = 0.
num_conditioned = 0
num_unconditioned = 0
for e in xrange(num_examples):
    is_conditioned, example = generate_example()
    #print example
    LL_ratio = CL.cross_validate(example, 2, K, num_folds=5)
    #print is_conditioned, LL_ratio
    if is_conditioned:
        conditioned_LL_ratio += LL_ratio
        num_conditioned += 1
    else:
        unconditioned_LL_ratio += LL_ratio
        num_unconditioned += 1
print conditioned_LL_ratio / num_conditioned
print unconditioned_LL_ratio / num_unconditioned
