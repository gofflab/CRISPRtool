#!/usr/bin/env python

import numpy as np

weights = [
    0.000,  # position 0
    0.000,
    0.014,
    0.000,
    0.000,
    0.395,  # position 5
    0.317,
    0.000,
    0.389,
    0.079,
    0.445,  # position 10
    0.508,
    0.613,
    0.851,
    0.732,
    0.828,  # position 15
    0.615,
    0.804,
    0.685,
    0.583,
]

def term1(x, y):
    '''
    Accepts two strings with lengths equal to that of the weights table.
    '''
    assert len(x) == len(y) == len(weights)

    penalties = [(1 - w)
                 for (a, b, w) in zip(x, y, weights)
                 if (a != b)]

    return np.product(penalties)

def term2(x, y):
    '''
    Accepts two strings of equal length.
    '''
    assert len(x) == len(y)

    mismatches = [i for i, (a, b) in enumerate(zip(x, y)) if a != b]
    if len(mismatches) < 2:
        return 1.0

    pairs = zip(mismatches, mismatches[1:])
    distances = [(b - a) for (a, b) in pairs]
    mean_pairwise_distance = np.mean(distances)

    return 1.0 / (1.0 + 4.0 * (19.0 - mean_pairwise_distance) / 19.0)

def term3(x, y):
    '''
    Accepts two strings of equal length.
    '''
    assert len(x) == len(y)

    mismatches = [1 for (a, b) in zip(x, y) if a != b]
    if len(mismatches) < 2:
        return 1.0

    return 1.0 / len(mismatches) ** 2

def calculateScore(x, y):
    return 100 * term1(x, y) * term2(x, y) * term3(x, y)

