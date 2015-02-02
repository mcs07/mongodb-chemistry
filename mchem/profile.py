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
import time

import numpy as np

from .similarity import similarity_search_fp


log = logging.getLogger(__name__)


def profile_similarity(mols, fingerprinter, fp_collection, result_collection, threshold=0.8, count_collection=None):
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
        results = similarity_search_fp(qfp, fp_collection, threshold, count_collection)
        end = time.time()
        log.debug('Query molecule %s of %s: %s results in %ss' % (i+1, len(mols), len(results), end-start))
        times.append(end - start)
    result['median_time'] = np.median(times)
    result['mean_time'] = np.mean(times)
    result_collection.insert(result)
