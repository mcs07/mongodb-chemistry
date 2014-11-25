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
from rdkit import Chem

from mchem.fps import MorganFingerprinter
from mchem.similarity import similarity_search

log = logging.getLogger(__name__)


def choose_sample(db):
    """Choose 1000 random molecules from the database."""
    chembl_ids = [m['_id'] for m in db.chembl.find().sort('_id')]
    print(len(chembl_ids))
    random.seed(201405291515)
    rands = random.sample(chembl_ids, 1000)
    for rand in rands:
        print(rand)


def load_query_mols(db):
    with open('data/sample_chembl_1000.txt') as f:
        chembl_ids = f.read().strip().split('\n')
    return [Chem.Mol(mol['rdmol']) for mol in db.chembl.find({'_id': {'$in': chembl_ids}})]


def test_constraints(db, name, reqbits=True, counts=True, rarest=True):
    """Test how different constraints can screen the molecule collection."""
    qmols = load_query_mols(db)[:100]
    fingerprinter = MorganFingerprinter(radius=2, length=None)
    result_fields = {'name': name, 'fingerprint': 'morgan', 'radius': 2, 'sample': 'chembl_1000', 'total': db.chembl.m2.count()}
    for threshold in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]:
        remaining = []
        result = result_fields.copy()
        log.info('Threshold: %s' % threshold)
        result['threshold'] = threshold
        for i, qmol in enumerate(qmols):
            remain = screen(qmol, fingerprinter, db.chembl.m2, threshold, db.chembl.m2.counts, reqbits, counts, rarest)
            log.debug('Query molecule %s of %s: %s remaining' % (i+1, len(qmols), remain))
            remaining.append(remain)
        result['median_remaining'] = np.median(remaining)
        result['mean_remaining'] = np.mean(remaining)
        db.profile.constraint.insert(result)


def test_folding(db, fp_collection, count_collection, length=None):
    """Analyse the proportion of molecules that are screened for fingerprints folded to different sizes."""
    qmols = load_query_mols(db)[:100]
    fingerprinter = MorganFingerprinter(radius=2, length=length)
    result_fields = {'length': length, 'fingerprint': 'morgan', 'radius': 2, 'sample': 'chembl_1000', 'total': db.chembl.m2.count()}
    for threshold in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]:
        remaining = []
        result = result_fields.copy()
        log.info('Threshold: %s' % threshold)
        result['threshold'] = threshold
        for i, qmol in enumerate(qmols):
            remain = screen(qmol, fingerprinter, fp_collection, threshold, count_collection)
            log.debug('Query molecule %s of %s: %s remaining' % (i+1, len(qmols), remain))
            remaining.append(remain)
        result['median_remaining'] = np.median(remaining)
        result['mean_remaining'] = np.mean(remaining)
        db.profile.folding.insert(result)


def test_radius(db, fp_collection, count_collection, radius=2):
    """Analyse the proportion of molecules that are screened for fingerprints of different radius."""
    qmols = load_query_mols(db)[:100]
    fingerprinter = MorganFingerprinter(radius=radius)
    result_fields = {'radius': radius, 'fingerprint': 'morgan', 'sample': 'chembl_1000', 'total': db.chembl.m2.count()}
    for threshold in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]:
        remaining = []
        result = result_fields.copy()
        log.info('Threshold: %s' % threshold)
        result['threshold'] = threshold
        for i, qmol in enumerate(qmols):
            remain = screen(qmol, fingerprinter, fp_collection, threshold, count_collection)
            log.debug('Query molecule %s of %s: %s remaining' % (i+1, len(qmols), remain))
            remaining.append(remain)
        result['median_remaining'] = np.median(remaining)
        result['mean_remaining'] = np.mean(remaining)
        db.profile.radius.insert(result)


def screen(mol, fingerprinter, fp_collection, threshold=0.8, count_collection=None, reqbits=True, counts=True, rarest=True):
    """Return the number of molecules remaining after screening."""
    qfp = fingerprinter.generate(mol)
    qn = len(qfp)                           # Number of bits in query fingerprint
    qmin = int(ceil(qn * threshold))        # Minimum number of bits in results fingerprints
    qmax = int(qn / threshold)              # Maximum number of bits in results fingerprints
    ncommon = qn - qmin + 1                 # Number of fingerprint bits in which at least one must be in common
    query = {}
    if reqbits:
        if count_collection and rarest:
            reqbits = [count['_id'] for count in count_collection.find({'_id': {'$in': qfp}}).sort('count', 1).limit(ncommon)]
        else:
            reqbits = qfp[:ncommon]
        query['bits'] = {'$in': reqbits}
    if counts:
        query['count'] = {'$gte': qmin, '$lte': qmax}
    remaining = fp_collection.find(query).count()
    return remaining
