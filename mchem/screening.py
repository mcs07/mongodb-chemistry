# -*- coding: utf-8 -*-
"""
mchem.screening
~~~~~~~~~~~~~~~

Functions for analysing screening methods for chemical searches in MongoDB.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging
from math import ceil

import numpy as np


log = logging.getLogger(__name__)


def screen(mol, fingerprinter, fp_collection, threshold=0.8, count_collection=None, reqbits=True, counts=True):
    """Return the number of molecules remaining after screening."""
    qfp = fingerprinter.generate(mol)
    qn = len(qfp)                           # Number of bits in query fingerprint
    qmin = int(ceil(qn * threshold))        # Minimum number of bits in results fingerprints
    qmax = int(qn / threshold)              # Maximum number of bits in results fingerprints
    ncommon = qn - qmin + 1                 # Number of fingerprint bits in which at least one must be in common
    query = {}
    if reqbits:
        if count_collection:
            # Use the count_collection to specifically get the rarest required bits
            reqbits = [count['_id'] for count in count_collection.find({'_id': {'$in': qfp}}).sort('count', 1).limit(ncommon)]
        else:
            # Just randomly choose the required bits
            reqbits = qfp[:ncommon]
        query['bits'] = {'$in': reqbits}
    if counts:
        query['count'] = {'$gte': qmin, '$lte': qmax}
    remaining = fp_collection.find(query).count()
    return remaining


def test_screening(mols, fingerprinter, fp_collection, result_collection, threshold=0.8, count_collection=None, reqbits=True, counts=True):
    """Test how different screening methods can constrain the molecule collection."""
    log.info('Testing screening: fp: %s Threshold: %s' % (fingerprinter.name, threshold))
    log.info('Methods: counts: %s reqbits: %s rarest: %s' % (counts, reqbits, count_collection is not None))
    result = {
        'fp': fingerprinter.name,
        'threshold': threshold,
        'remaining': [],
        'total': fp_collection.count(),
        'reqbits': reqbits,
        'counts': counts,
        'rarest': count_collection is not None
    }
    for i, qmol in enumerate(mols):
        remain = screen(qmol, fingerprinter, fp_collection, threshold, count_collection, reqbits, counts)
        log.debug('Query molecule %s of %s: %s remaining' % (i+1, len(mols), remain))
        result['remaining'].append(remain)
    result['median_remaining'] = np.median(result['remaining'])
    result['mean_remaining'] = np.mean(result['remaining'])
    result_collection.insert(result)


