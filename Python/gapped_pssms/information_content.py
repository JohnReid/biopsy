#
# Copyright John Reid 2009
#

"""
Code to examine information content of gapped PWMs.
"""

import numpy as N, fpconst


def inc_index(idx, max=4):
    "Increment the index."
    for i in xrange(len(idx)):
        idx[i] = (idx[i] + 1) % max
        if 0 != idx[i]:
            break
    else:
        return False
    return True



def ic(pwm, bg):
    """
    Calculate the information content of the PWM against the background in nats. Both pwm and bg are expected to be
    callable objects that return the log likelihood.
    """
    ic = 0.
    ix = N.zeros(pwm.K, dtype=int)
    while True:
        l = pwm(ix)
        if not fpconst.isNegInf(l):
            x = N.exp(l) * (l - bg(ix))
            assert not fpconst.isNaN(x)
            ic += x
        if not inc_index(ix):
            break
    return ic



def eval_pwm(pwm, X):
    "Evaluate the PWM on X (don't account for rev comps : in one direction only). Assumes pwm in log likelihood form."
    return sum(p[x] for p, x in zip(pwm, X))



class GappedPWM(object):
    "A gapped PWM."

    def __init__(self, pwm, gap_char, gap_freq):
        "Construct."
        self.gapped_pwm = pwm
        self.ungapped_pwm = pwm.copy()
        self.ungapped_pwm[gap_char:-1] = pwm[gap_char+1:]
        self.ungapped_pwm[-1,:] = .25
        self.gapped_pwm = N.log(self.gapped_pwm)
        self.gapped_pwm_r = self.gapped_pwm[::-1,::-1].copy()
        self.ungapped_pwm = N.log(self.ungapped_pwm)
        self.ungapped_pwm_r = self.ungapped_pwm[::-1,::-1].copy()
        self.gap_char = gap_char
        self.gap_freq = gap_freq
        self.K = len(pwm)

    def __call__(self, X):
        "Evaluate the gapped PWM."
        return N.log(
            .5*(
                self.gap_freq*(N.exp(eval_pwm(self.gapped_pwm, X)) + N.exp(eval_pwm(self.gapped_pwm_r, X)))
                + (1.-self.gap_freq)*(N.exp(eval_pwm(self.ungapped_pwm, X)) + N.exp(eval_pwm(self.ungapped_pwm_r, X)))
            )
        )

    def logo(self):
        "Create a logo of the gapped PWM."
        import hmm.pssm.logo as L
        transparencies = N.ones(self.K)
        transparencies[self.gap_char] = self.gap_freq
        return L.pssm_as_image(N.exp(self.gapped_pwm), size=(160*self.K,480), transparencies=transparencies)


class MultipleGappedPWM(object):
    "A gapped PWM with multiple gaps."

    def __init__(self, freqs, gaps):
        "Construct."
        freqs = (freqs.T / freqs.sum(axis=1)).T
        self.gaps = list(gaps)
        self.K = len(freqs)
        ix = N.zeros(len(self.gaps), dtype=int)
        self.pwms = []
        while True: # for each possible combination of gaps
            pwm_freqs = freqs.copy()
            missed_gaps = 0
            p = 1.
            for i in xrange(len(self.gaps)):
                gap_char, gap_freq = self.gaps[i]
                if ix[i]:
                    p *= gap_freq
                else:
                    p *= (1.-gap_freq)
                    copy_len = self.K - 1 - gap_char
                    pwm_freqs[gap_char-missed_gaps:gap_char-missed_gaps+copy_len] = freqs[gap_char+1:]
                    missed_gaps += 1
            if missed_gaps:
                pwm_freqs[-missed_gaps:] = .25
            self.pwms.append((p, PWM(pwm_freqs)))
            if not inc_index(ix, max=2):
                break


    def __call__(self, X):
        "Evaluate the gapped PWM."
        return N.log(sum(p*N.exp(pwm(X)) for p, pwm in self.pwms))

    def logo(self):
        "Create a logo of the gapped PWM."
        import hmm.pssm.logo as L
        transparencies = N.ones(self.K)
        for gap_char, gap_freq in self.gaps:
            transparencies[gap_char] = gap_freq
        return L.pssm_as_image(N.exp(self.pwms[-1][1].freqs), size=(160*self.K,480), transparencies=transparencies)

    def logos(self):
        "Create a logo for the standard PWM representing each possible combination of gaps."
        import hmm.pssm.logo as L
        return [L.pssm_as_image(N.exp(pwm.freqs), size=(160*self.K,480)) for p, pwm in self.pwms]


class PWM(object):
    "A standard PWM."

    def __init__(self, freqs):
        "Construct."
        self.freqs = freqs
        self.freqs = N.log(self.freqs)
        self.freqs_r = self.freqs[::-1,::-1].copy()
        self.K = len(freqs)

    def __call__(self, X):
        "Evaluate the gapped PWM."
        return N.log(.5*(N.exp(eval_pwm(self.freqs, X)) + N.exp(eval_pwm(self.freqs_r, X))))

    def logo(self):
        "Create a logo of the gapped PWM."
        import hmm.pssm.logo as L
        transparencies = N.ones(self.K)
        return L.pssm_as_image(N.exp(self.freqs), size=(160*self.K,480))




class BG(object):
    "Simple background model."

    def __init__(self, K):
        "Construct."
        self.p = K * N.log(.25)

    def __call__(self, X):
        return self.p




def test_index():
    idx = N.zeros(3, dtype=int)
    print idx
    while inc_index(idx):
        print idx



if '__main__' == __name__:
    def output_pwm(f, tag, gapped_pwm):
        logo = gapped_pwm.logo()
        logo.save('ic/logo-%s.png' % tag)
        logo.save('ic/logo-%s.eps' % tag)
        print >> f, ('\\includegraphics[height=100pt]{logo-%s} \\\\' % tag)
        print >> f, 'Information content = %.2f bits\n\n' % (ic(gapped_pwm, BG(gapped_pwm.K)) / N.log(2.))

    def gapped_results():
        # TRANSFAC sp1
        transfac_sp1_freqs = N.array(
            [
                [ 37      ,17      ,137      ,44 ],
                [ 26     ,4     ,204     ,1  ]   ,
                [ 0     ,3     ,230     ,2   ]  ,
                [ 0     ,0     ,234     ,1   ]  ,
                [ 40     ,159     ,5     ,31 ],
                [ 7     ,6     ,214     ,8   ]  ,
                [ 1     ,2     ,212     ,20 ],
                [ 39     ,16     ,169     ,11 ],
                [ 31     ,31     ,155     ,18 ],
                [ 31     ,119     ,40     ,45 ],
            ],
            dtype=N.float64
        )
        transfac_sp1_pwm = MultipleGappedPWM(transfac_sp1_freqs, [])

        # Sp1 frequencies and gap
        gapped_sp1_freqs = N.array(
            [
                [ 0.0933702382725,0.242248047918,0.566206988199,0.0981747256098        ],
                [ 1.16094945149e-08,0.67653521585,0.00610826470716,0.317356507834      ],
                [ 7.13315935192e-11,0.87168908789,6.22903972303e-05,0.128248621642     ],
                [ 1.43918642857e-07,0.993525315867,1.28219159736e-16,0.00647454021475  ],
                [ 1.20896878223e-08,0.999999215443,9.30869032889e-19,7.72467729618e-07 ],
                [ 0.0356652209825,0.0185904602858,1.06129495172e-11,0.945744318721     ],
                [ 0.125408435326,3.93276509946e-15,0.648809170375,0.225782394299       ],
                [ 7.98397498566e-09,0.99919057356,0.000160248916246,0.000649169539886  ],
                [ 1.88950891252e-10,0.998649985294,2.45236030926e-13,0.00135001451637  ],
                [ 1.82878927635e-15,0.944192619205,1.26450262074e-11,0.0558073807819   ],
                [ 0.17319780881,0.779436828144,6.97591450816e-16,0.0473653630463       ],
                [ 0.246824694996,0.453994213412,0.137368283321,0.16181280827           ],
                [ 0.279186787399,0.285921551901,0.285684295696,0.149207365005          ],
            ]
        )
        gapped_sp1_pwm = MultipleGappedPWM(gapped_sp1_freqs, [(5, 0.321088)])

        f = open('ic/gapped-results.tex', 'w')
        output_pwm(f, 'TRANSFAC-Sp1', transfac_sp1_pwm)
        output_pwm(f, 'Gapped-Sp1', gapped_sp1_pwm)
        f.close()
    gapped_results()



    def transfac_analysis():
        pou_gapped = MultipleGappedPWM(
            freqs=N.array(
                [
                    [ 0.076923077 , 0 , 0.923076923 , 0 ],
                    [ 0.066666667 , 0.933333333 , 0 , 0 ],
                    [ 1 , 0 , 0 , 0 ],
                    [ 0 , 0 , 0 , 1 ],
                    [ 1 , 0 , 0 , 0 ],
                    [ 0 , 0 , 0 , 1 ],
                    [ 0.75 , 0.125 , 0.125 , 0 ],
                    [ 0.571428571 , 0 , 0 , 0.428571429 ],
                    [ 0 , 0 , 0 , 1 ],
                    [ 0 , 0 , 0 , 1 ],
                    [ 1 , 0 , 0 , 0 ],
                    [ 0.1875 , 0 , 0 , 0.8125 ],
                ]
            ),
            gaps = [(5,0.1875), (7,0.875)]
        )

        pou_standard = PWM(
            freqs=N.array(
                [
                    [ 0.1875 , 0.0625 , 0.75 , 0 ],
                    [ 0.125 , 0.8125 , 0 , 0.0625 ],
                    [ 0.875 , 0 , 0 , 0.125 ],
                    [ 0.125 , 0 , 0 , 0.875 ],
                    [ 0.8125 , 0.0625 , 0 , 0.125 ],
                    [ 0.8125 , 0.0625 , 0.125 , 0 ],
                    [ 0.4375 , 0 , 0 , 0.5625 ],
                    [ 0 , 0 , 0 , 1 ],
                    [ 0.1875 , 0 , 0 , 0.8125 ],
                    [ 1 , 0 , 0 , 0 ],
                    [ 0 , 0 , 0.0625 , 0.9375 ],
                ]
            )
        )

        mef2_gapped = MultipleGappedPWM(
            freqs=N.array(
                [
                    [ 0.307692308 , 0 , 0.461538462 , 0.230769231 ],
                    [ 0 , 0.230769231 , 0.692307692 , 0.076923077 ],
                    [ 0 , 0.692307692 , 0 , 0.307692308 ],
                    [ 0 , 0 , 0 , 1 ],
                    [ 0.769230769 , 0.076923077 , 0 , 0.153846154 ],
                    [ 0 , 0 , 0 , 1 ],
                    [ 0.625 , 0 , 0 , 0.375 ],
                    [ 0.076923077 , 0 , 0 , 0.923076923 ],
                    [ 0.153846154 , 0 , 0 , 0.846153846 ],
                    [ 0 , 0 , 0 , 1 ],
                    [ 1 , 0 , 0 , 0 ],
                    [ 0.384615385 , 0.076923077 , 0.384615385 , 0.153846154 ]
                ]
            ),
            gaps = [(6,0.615384615384615),]
        )

        mef2_standard = PWM(
            freqs=N.array(
                [
                    [ 0.307692308 , 0 , 0.461538462 , 0.230769231 ],
                    [ 0 , 0.153846154 , 0.769230769 , 0.076923077 ],
                    [ 0 , 0.692307692 , 0 , 0.307692308 ],
                    [ 0 , 0.076923077 , 0 , 0.923076923 ],
                    [ 0.769230769 , 0.076923077 , 0 , 0.153846154 ],
                    [ 0 , 0 , 0 , 1 ],
                    [ 0.384615385 , 0 , 0 , 0.615384615 ],
                    [ 0.076923077 , 0 , 0 , 0.923076923 ],
                    [ 0.153846154 , 0 , 0 , 0.846153846 ],
                    [ 0.307692308 , 0 , 0 , 0.692307692 ],
                    [ 0.846153846 , 0 , 0.153846154 , 0 ],
                    [ 0.384615385 , 0.076923077 , 0.384615385 , 0.153846154 ],
                ]
            )
        )

        f = open('ic/transfac-analysis.tex', 'w')
        output_pwm(f, 'POU-gapped', pou_gapped)
        output_pwm(f, 'POU-standard', pou_standard)
        output_pwm(f, 'MEF2-gapped', mef2_gapped)
        output_pwm(f, 'MEF2-standard', mef2_standard)
        f.close()


    def output_example_pwms():
        I=.85
        O=.05
        f = open('ic/pwms.tex', 'w')
        output_pwm(
            f,
            'rev-comp',
            GappedPWM(
                pwm=N.array([
                    [O,O,O,I],
                    [O,O,O,I],
                    [O,O,O,I],
                    [O,O,O,I],
                    [I,O,O,O],
                    [I,O,O,O],
                    [I,O,O,O],
                ]),
                gap_char=3,
                gap_freq=.5
            )
        )
        output_pwm(
            f,
            'non-rev-comp',
            GappedPWM(
                pwm=N.array([
                    [O,O,O,I],
                    [O,O,O,I],
                    [O,O,O,I],
                    [O,O,O,I],
                    [O,O,O,I],
                    [O,O,O,I],
                    [O,O,O,I],
               ]),
                gap_char=3,
                gap_freq=.5
            )
        )
        output_pwm(
            f,
            'low-ent',
            GappedPWM(
                pwm=N.array([
                    [O,O,O,I],
                    [O,O,I,O],
                    [O,I,O,O],
                    [I,O,O,O],
                    [O,O,O,I],
                    [O,O,I,O],
                    [O,I,O,O],
                ]),
                gap_char=3,
                gap_freq=.5
            )
        )
        output_pwm(
            f,
            'without-gap',
            GappedPWM(
                pwm=N.array([
                    [O,O,O,I],
                    [O,O,I,O],
                    [O,I,O,O],
                    [O,O,O,I],
                    [O,O,I,O],
                    [O,I,O,O],
                ]),
                gap_char=3,
                gap_freq=1.
            )
        )
        f.close()
