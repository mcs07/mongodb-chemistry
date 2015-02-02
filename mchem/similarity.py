# -*- coding: utf-8 -*-
"""
mchem.similarity
~~~~~~~~~~~~~~~~

Functions for performing similarity search queries.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging
from math import ceil


log = logging.getLogger(__name__)


def similarity_client(mol, fingerprinter, fp_collection,  threshold=0.8, count_collection=None):
    """Perform a similarity search on the client, with initial screening to improve performance."""
    log.info('Similarity search with %s and threshold %s' % (fp_collection.name, threshold))
    qfp = fingerprinter.generate(mol)
    qn = len(qfp)                           # Number of bits in query fingerprint
    qmin = int(ceil(qn * threshold))        # Minimum number of bits in results fingerprints
    qmax = int(qn / threshold)              # Maximum number of bits in results fingerprints
    ncommon = qn - qmin + 1                 # Number of fingerprint bits in which at least one must be in common
    # Get list of bits where at least one must be in result fp. Use least popular bits if possible.
    if count_collection:
        reqbits = [count['_id'] for count in count_collection.find({'_id': {'$in': qfp}}).sort('count', 1).limit(ncommon)]
    else:
        reqbits = qfp[:ncommon]
    results = []
    for fp in fp_collection.find({'bits': {'$in': reqbits}, 'count': {'$gte': qmin, '$lte': qmax}}):
        intersection = len(set(qfp) & set(fp['bits']))
        pn = fp['count']
        tanimoto = float(intersection) / (pn + qn - intersection)
        if tanimoto >= threshold:
            results.append({'_id': fp['_id'], 'tanimoto': tanimoto})
    return results


def similarity_search(mol, fingerprinter, fp_collection, threshold=0.8, count_collection=None):
    """Perform a similarity search using aggregation framework.

    :param mol: The query molecule
    :param fingerprinter: The Fingerprinter to use
    :param fp_collection: MongoDB fingerprint collection to query
    :param threshold: The tanimoto threshold
    :param count_collection: MongoDB collection containing fingerprint bit frequencies
    """
    qfp = fingerprinter.generate(mol)
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


def similarity_suvee(mol, fingerprinter, fp_collection, threshold=0.8):
    """Perform a similarity search using aggregation using method proposed by Davy Suvee.

    Read more at http://datablend.be/?p=265 http://datablend.be/?p=256 http://datablend.be/?p=254

    :param qfp: The query fingerprint
    :param threshold: The tanimoto threshold
    """
    qfp = fingerprinter.generate(mol)
    qn = len(qfp)                           # Number of bits in query fingerprint
    qmin = qn * threshold                   # Minimum number of bits in results fingerprints
    qmax = qn / threshold                   # Maximum number of bits in results fingerprints
    aggregate = [
        {'$match': {'count': {'$gte': qmin, '$lte': qmax}}},
        {'$unwind': '$bits'},
        {'$match': {'bits': {'$in': qfp}}},
        {'$group': {
            '_id': '$_id',
            'common': {'$sum': 1},
            'count': {'$first': '$count'},
        }},
        {'$project': {'tanimoto': {'$divide': ['$common', {'$subtract': [{'$add': [qn, '$count']}, '$common']}]}}},
        {'$match': {'tanimoto': {'$gte': threshold}}}
    ]
    response = fp_collection.aggregate(aggregate)
    return response['result']
