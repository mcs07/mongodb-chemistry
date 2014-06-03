#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Perform chemical similarity and substructures searches in MongoDB."""

import pymongo
from rdkit import Chem
from rdkit.Chem import AllChem


db = pymongo.MongoClient().chem


def similarity1(qfp, threshold=0.8):
    """Perform a similarity search without any efficiency optimisations."""
    results = []
    for molecule in db.molecules.find():
        dbmol = Chem.Mol(molecule['rdmol'])
        mfp = list(AllChem.GetMorganFingerprintAsBitVect(dbmol, 2, nBits=2048).GetOnBits())
        intersection = len(set(qfp) & set(mfp))
        tanimoto = float(intersection) / (len(mfp) + len(qfp) - intersection)
        if tanimoto >= threshold:
            #print '{} : {} : {}'.format(tanimoto, molecule['chembl_id'], molecule['smiles'])
            results.append({'molecule': molecule, 'tanimoto': tanimoto})
    return results


def similarity2(qfp, threshold=0.8):
    """Perform a similarity search using efficiency optimisations."""
    qn = len(qfp)                           # Number of bits in query fingerprint
    qmin = int(qn * threshold)              # Minimum number of bits in results fingerprints
    qmax = int(qn / threshold)              # Maximum number of bits in results fingerprints
    ncommon = qn - qmin + 1                 # Number of fingerprint bits in which at least 1 must be in common
    # Query mfp2_counts collection for the least popular bits that are in qfp
    reqbits = [count['_id'] for count in db.mfp2_counts.find({'_id': {'$in': qfp}}).sort('count', 1).limit(ncommon)]
    #reqbits = qfp[:ncommon]
    results = []
    for molecule in db.molecules.find({'mfp2.bits': {'$in': reqbits}, 'mfp2.count': {'$gte': qmin, '$lte': qmax}}):
        intersection = len(set(qfp) & set(molecule['mfp2']['bits']))
        pn = molecule['mfp2']['count']
        tanimoto = float(intersection) / (pn + qn - intersection)
        if tanimoto >= threshold:
            #print '{} : {} : {}'.format(tanimoto, molecule['chembl_id'], molecule['smiles'])
            results.append({'molecule': molecule, 'tanimoto': tanimoto})
    return results


def similarity3(qfp, threshold=0.8):
    """Perform a similarity search using aggregation as proposed by Davy Suvee.

    :param qfp: The query fingerprint
    :param threshold: The tanimoto threshold
    """
    qn = len(qfp)                           # Number of bits in query fingerprint
    qmin = int(qn * threshold)               # Minimum number of bits in results fingerprints
    qmax = int(qn / threshold)               # Maximum number of bits in results fingerprints
    #ncommon = qn - qmin + 1                 # Number of fingerprint bits in which at least 1 must be in common
    #reqbits = [count['_id'] for count in db.mfp2_counts.find({'_id': {'$in': qfp}}).sort('count', 1).limit(ncommon)]
    #reqbits = qfp[:ncommon]
    aggregate = [
        {'$match': {'mfp2.count': {'$gte': qmin, '$lte': qmax}}},
        {'$unwind': '$mfp2.bits'},
        {'$match': {'mfp2.bits': {'$in': qfp}}},
        {'$group': {
            '_id': '$_id',
            'common': {'$sum': 1},
            'count': {'$first': '$mfp2.count'},
            'smiles': {'$first': '$smiles'},
            'chembl_id': {'$first': '$chembl_id'}
        }},
        {'$project': {
            'tanimoto': {'$divide': ['$common', {'$subtract': [{'$add': [qn, '$count']}, '$common']}]},
            'smiles': 1,
            'chembl_id': 1
        }},
        {'$match':  {'tanimoto': {'$gte': threshold}}}
    ]
    response = db.molecules.aggregate(aggregate)
    # for result in response['result']:
    #     print '%s : %s : %s' % (result['tanimoto'] * 100, result['chembl_id'], result['smiles'])
    return response['result']


def similarity4(qfp, threshold=0.8):
    """Perform a similarity search using aggregation with new operators from MongoDB 2.6.

    :param qfp: The query fingerprint
    :param threshold: The tanimoto threshold
    """
    qn = len(qfp)                           # Number of bits in query fingerprint
    qmin = int(qn * threshold)               # Minimum number of bits in results fingerprints
    qmax = int(qn / threshold)               # Maximum number of bits in results fingerprints
    ncommon = qn - qmin + 1                 # Number of fingerprint bits in which at least 1 must be in common
    reqbits = [count['_id'] for count in db.mfp2_counts.find({'_id': {'$in': qfp}}).sort('count', 1).limit(ncommon)]
    #reqbits = qfp[:ncommon]
    aggregate = [
        {'$match': {'mfp2.count': {'$gte': qmin, '$lte': qmax}, 'mfp2.bits': {'$in': reqbits}}},
        {'$project': {
            'tanimoto': {'$let': {
                'vars': {'common': {'$size': {'$setIntersection': ['$mfp2.bits', qfp]}}},
                'in': {'$divide': ['$$common', {'$subtract': [{'$add': [qn, '$mfp2.count']}, '$$common']}]}
            }},
            'smiles': 1,
            'chembl_id': 1
        }},
        {'$match':  {'tanimoto': {'$gte': threshold}}}
    ]
    response = db.molecules.aggregate(aggregate)
    # for result in response['result']:
    #     print '%s : %s : %s' % (result['tanimoto'] * 100, result['chembl_id'], result['smiles'])
    return response['result']


if __name__ == '__main__':
    # start = time.time()
    # similarity('C1=CC(=C(C=C1C(=N)N)Br)OCCCOC2=C(C=C(C=C2)C(=N)N)Br', 0.7)
    # end = time.time()
    # print 'Time: %s' % (end - start)
    # start = time.time()
    # similarity2('C1=CC(=C(C=C1C(=N)N)Br)OCCCOC2=C(C=C(C=C2)C(=N)N)Br', 0.7)
    # end = time.time()
    # print 'Time: %s' % (end - start)
    # start = time.time()
    #similarity3(qfp=[], tanimoto=0.9, fp='morganbv_fp')
    # end = time.time()
    # print 'Time: %s' % (end - start)

    #qmol = db.molecules_morganbv_fp.find_one(skip=123)
    #start = time.time()
    #similarity3(qfp=qmol['morganbv_fp'], tanimoto=0.7, fp='morganbv_fp')
    #end = time.time()
    #print 'Time: %s' % (end - start)

    qmol = db.molecules.find_one(skip=123456)
    similarity2(qmol['mfp']['bits'], 0.9)
