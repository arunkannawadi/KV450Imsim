def match(R,Q):
    #    """ Returns idx, cuts such that
    #        R[idx] = Q[cuts]
    #    """

    import numpy as np
    from collections import Counter

    Q_argsort = np.argsort(Q)
    R_argsort = np.argsort(R)

    Q_s = Q[Q_argsort]
    Q_us, Qu_inverse = np.unique(Q_s, return_index=False, return_inverse=True)

    superset = set(R)
    subset = set(Q_us)
    assert set(Q)==set(Q_us)

    difference = subset.difference(superset)
    D = np.array(list(difference))

    RD = np.hstack([R,D])
    RD_argsort = np.argsort(RD)

    ## Track D
    D_matching_indices = np.searchsorted(RD,D,sorter=RD_argsort)
    RD_matching_indices = np.searchsorted(RD, Q_us, sorter=RD_argsort) ## Q_us = RD[RD_argsort][RD_matching_indices]
    R_matching_indices = np.array(list(set(RD_matching_indices).difference(set(D_matching_indices))))
    ## We have to do R_matching_indices in such a complicated way. Because D is guaranteed to be fully matched, but not R

    np.testing.assert_array_equal(Q_us, RD[RD_argsort][RD_matching_indices])
    np.testing.assert_array_equal(Q_us, RD[RD_argsort[RD_matching_indices]])

    ## Identify those entries in Q which don't have a match in R
    try:
        Q_us_keep_indices = np.searchsorted(Q_us, RD[RD_argsort][R_matching_indices])
    except IndexError as IE:
        print IE
        print "That's okay. It is never used again"

    cntr = Counter(Q_s)
    for k in RD[RD_argsort][D_matching_indices]:
        cntr.pop(k)

    Q_s_recon = np.repeat(cntr.keys(), repeats=cntr.values())
    Q_s_recon.sort()

    precuts = np.searchsorted(Q_s,Q_s_recon)
    addendum = np.zeros_like(precuts)
    num = 0
    for k in xrange(1,len(precuts)):
        if precuts[k]==precuts[k-1]:
            num += 1
        else:
            num = 0
        addendum[k] = num

    cuts = Q_argsort[precuts+addendum]
    idx = R_argsort[np.searchsorted(R, Q[cuts], sorter=R_argsort)]

    return idx, cuts

def test_match():
    import numpy as np

    Q = np.array([-4, -2,  0,  2,  4, -2,  8,  4, 10, 13, -2,  2, -2,  4, 13])
    R = np.array([ 1,  3,  4,  2,  5, 10,  6,  7,  9, 11, 12])

    idx, cuts = match(R,Q)

    print "Q = ", Q
    print "R = ", R

    print "Q[cuts] = ", Q[cuts]
    print "R[idx] = ", R[idx]

    print "cuts = ", cuts
    print "idx = ", idx

if __name__=='__main__':
    test_match()
