# -*- coding: utf-8 -*-
"""
mchem.profile
~~~~~~~~~~~~~

Functions for benchmarking chemical searches in MongoDB.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging
from math import ceil
import time

import numpy as np


log = logging.getLogger(__name__)


def profile_search(mols, fingerprinter, fp_collection, result_collection, threshold=0.8, count_collection=None):
    """Benchmark similarity search."""
    log.info('Benchmarking: fp: %s Threshold: %s' % (fingerprinter.name, threshold))
    result = {
        'fp': fingerprinter.name,
        'threshold': threshold,
        'total': fp_collection.count()
    }
    times = []
    for i, qmol in enumerate(mols):
        qfp = fingerprinter.generate(qmol)
        start = time.time()
        results = similarity_search(qfp, fp_collection, threshold, count_collection)
        end = time.time()
        log.debug('Query molecule %s of %s: %s results in %ss' % (i+1, len(mols), len(results), end-start))
        times.append(end - start)
    result['median_time'] = np.median(times)
    result['mean_time'] = np.mean(times)
    result_collection.insert(result)


def similarity_search(qfp, fp_collection, threshold=0.8, count_collection=None):
    """Perform a similarity search using aggregation framework.

    :param qfp: The query fingerprint
    :param threshold: The tanimoto threshold
    """
    qn = len(qfp)                           # Number of bits in query fingerprint
    qmin = int(ceil(qn * threshold))        # Minimum number of bits in results fingerprints
    qmax = int(qn / threshold)              # Maximum number of bits in results fingerprints
    ncommon = qn - qmin + 1                 # Number of fingerprint bits in which at least 1 must be in common
    if count_collection:
        reqbits = [count['_id'] for count in count_collection.find({'_id': {'$in': qfp}}).sort('count', 1).limit(ncommon)]
    else:
        reqbits = qfp[:ncommon]
    aggregate = [
        {'$match': {'count': {'$gte': qmin, '$lte': qmax}, 'bits': {'$in': reqbits}}},
        {'$project': {
            'tanimoto': {'$let': {
                'vars': {'common': {'$size': {'$setIntersection': ['$bits', qfp]}}},
                'in': {'$divide': ['$$common', {'$subtract': [{'$add': [qn, '$count']}, '$$common']}]}
            }},
        }},
        {'$match': {'tanimoto': {'$gte': threshold}}}
    ]
    response = fp_collection.aggregate(aggregate)
    return response['result']
