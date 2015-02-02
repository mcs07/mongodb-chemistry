# -*- coding: utf-8 -*-
"""
mchem.fps
~~~~~~~~~

Functions for generating fingerprints using RDKit.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from collections import defaultdict
import logging

import pymongo
from rdkit import Chem
from rdkit.Chem import AllChem


log = logging.getLogger(__name__)


def generate(mol_collection, fp_collection, fingerprinter):
    """Generate a fingerprint for all molecules in a collection.

    Example::

        generate(mols, MorganFingerprinter(radius=2))
        generate(mols, MorganFingerprinter(radius=2, length=1024))

    :param mol_collection: MongoDB database collection containing molecules.
    :param fp_collection: MongoDB database collection to store fingerprints.
    :param fingerprinter: fingerprinter instance to generate fingerprint for each molecule.
    """
    log.info('Generating %s fingerprints for %s into %s' % (fingerprinter.name, mol_collection.name, fp_collection.name))
    success, skip = 0, 0
    for molecule in mol_collection.find(timeout=False):
        log.debug('Generating %s for %s' % (fingerprinter.name, molecule['_id']))
        bits = fingerprinter.generate(Chem.Mol(molecule['rdmol']))
        fp = {
            '_id': molecule['_id'],
            'bits': bits,
            'count': len(bits)
        }
        try:
            fp_collection.insert(fp)
            log.debug('Inserted fingerprint for %s' % fp['_id'])
            success += 1
        except pymongo.errors.DuplicateKeyError:
            log.debug('Skipped %s: Fingerprint already exists' % fp['_id'])
            skip += 1
    log.info('%s successes, %s skipped' % (success, skip))
    log.info('Ensuring index on bits and counts for %s' % fp_collection.name)
    fp_collection.ensure_index('bits')
    fp_collection.ensure_index('count')


def count(fp_collection, count_collection):
    """Build collection containing total counts of all occurrences of each fingerprint bit."""
    counts = defaultdict(int)
    count_collection.drop()
    log.info('Counting fingerprint bits in %s' % count_collection.name)
    for fp in fp_collection.find(timeout=False):
        log.debug('Processing %s' % fp['_id'])
        for bit in fp['bits']:
            counts[bit] += 1
    for k, v in counts.items():
        log.debug('Saving count %s: %s' % (k, v))
        count_collection.insert({'_id': k, 'count': v})


class Fingerprinter(object):
    """Fingerprinter interface."""

    def generate(self, mol):
        """Generate this fingerprint for a molecule."""
        raise NotImplementedError('Fingerprinter subclasses must implement a generate method')

    @property
    def name(self):
        """Unique name for this fingerprint."""
        raise NotImplementedError('Fingerprinter subclasses must implement a name property')


class MorganFingerprinter(Fingerprinter):
    """Class for generating morgan fingerprints."""

    def __init__(self, radius=2, length=None):
        """Initialize with a radius and an optional length.

        :param int radius: The Morgan fingerprint radius (default: 2).
        :param length: The number of bits to optionally fold the fingerprint down to.
        """
        self.radius = radius
        self.length = length

    def generate(self, mol):
        """Generate Morgan fingerprint for a molecule.

        :param mol: The RDKit Mol to generate the fingerprint for.
        """
        if self.length:
            fp = AllChem.GetHashedMorganFingerprint(mol, radius=self.radius, nBits=self.length)
        else:
            fp = AllChem.GetMorganFingerprint(mol, radius=self.radius)
        return sorted(fp.GetNonzeroElements().keys())

    @property
    def name(self):
        """A unique identifier for this fingerprint with the current settings."""
        n = 'm%s' % self.radius
        if self.length:
            n = '%sl%s' % (n, self.length)
        return n

