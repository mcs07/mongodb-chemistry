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
import datetime
import logging
import random
import time
from math import ceil

import numpy as np
import pymongo

from .similarity import similarity_search

log = logging.getLogger(__name__)


def choose_random():
    """Choose 1000 random molecules from the database."""
    db = pymongo.MongoClient().mchem
    chembl_ids = [m['_id'] for m in db.chembl.find().sort('_id')]
    print(len(chembl_ids))
    random.seed(201405291515)
    rands = random.sample(chembl_ids, 1000)
    print(rands)


def profile_screening(mol, fingerprinter, fp_collection, threshold=0.8, count_collection=None):
    """Analyse the proportion of molecules that can be screened quickly."""

    qfp = fingerprinter.generate(mol)
    qn = len(qfp)                           # Number of bits in query fingerprint
    qmin = int(ceil(qn * threshold))        # Minimum number of bits in results fingerprints
    qmax = int(qn / threshold)              # Maximum number of bits in results fingerprints
    ncommon = qn - qmin + 1                 # Number of fingerprint bits in which at least one must be in common
    if count_collection:
        reqbits = [count['_id'] for count in count_collection.find({'_id': {'$in': qfp}}).sort('count', 1).limit(ncommon)]
    else:
        reqbits = qfp[:ncommon]
    total = fp_collection.count()
    remaining = fp_collection.find({'bits': {'$in': reqbits}, 'count': {'$gte': qmin, '$lte': qmax}}).count()



# Screened percentage as a function of tanimoto threshold
# TNa⩽Nb⩽Na/T constraint
# Na−TNa+1 constraint
# Rarest bits for Na−TNa+1 constraint

# Once settled on the best constraint method:
# unfolded, 2048, 1024, 512 folded
# radius 2, 3, 4

# Benchmarks:
# Re-check all, in case other factors affect performance even though screening is better

# Note the impact on results - folded has small negative effect, increased radius has effect, subjective?


#
# def profile_mongodb():
#     """"""
#
#     db = pymongo.MongoClient().chem
#     for threshold in [0.95, 0.9, 0.85, 0.8, 0.75, 0.7]:
#         report = []
#         times = []
#         counts = []
#         # Perform the queries
#         for chembl_id in chembl_ids:
#             print(chembl_id)
#             qmol = db.molecules.find_one({'chembl_id': chembl_id})
#             report.append('Query: {} - {}'.format(qmol['chembl_id'], qmol['smiles']))
#             start = time.time()
#             results = similarity4(qmol['mfp2']['bits'], threshold)
#             end = time.time()
#             counts.append(len(results))
#             report.append('Results ({})'.format(len(results)))
#             print(report[-1])
#             for r in results:
#                 report.append('{}: {}'.format(r['tanimoto'], r['chembl_id']))
#                 print(report[-1])
#             times.append(end - start)
#
#         # Produce a report of the results
#         report.append('Counts: {}'.format(' '.join(str(c) for c in counts)))
#         report.append('Mean: {}'.format(np.mean(counts)))
#         report.append('Median: {}'.format(np.median(counts)))
#         report.append('Times: {}'.format(' '.join(str(t) for t in times)))
#         report.append('Mean: {}'.format(np.mean(times)))
#         report.append('Median: {}'.format(np.median(times)))
#         report.append('95th Percentile: {}'.format(np.percentile(times, 95)))
#         report = '\n'.join(report)
#         print(report)
#         with open('benchmarks/{}-{}-{}-{}.txt'.format('unfolded', 'agg', threshold, datetime.datetime.utcnow()), 'w') as f:
#             f.write(report)

# You can profile individual queries using built in MongoDB profiling:
#db.set_profiling_level(pymongo.ALL)
#db.set_profiling_level(pymongo.OFF)
# But this isn't so useful for searches made up of multiple database queries (it will underestimate)
