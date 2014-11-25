# -*- coding: utf-8 -*-
"""
mchem.plot
~~~~~~~~~~

Plot profiling results using matplotlib.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

log = logging.getLogger(__name__)


def plot_constraints(db):
    allresults = list(db.profile.constraint.find({'name': 'all'}).sort('threshold'))
    allthresholds = np.array([0] + [r['threshold'] for r in allresults])
    alldiscard = np.array([0] + [1 - r['mean_remaining'] / r['total'] for r in allresults])
    #plt.plot(allthresholds, alldiscard, 'ow', markeredgewidth=1)
    xfit = np.linspace(0, 1, 1000)
    spl = interpolate.UnivariateSpline(allthresholds, alldiscard, s=0)
    plt.plot(xfit, spl(xfit), linewidth=1.5, color='k', label='Combined')

    rareresults = list(db.profile.constraint.find({'name': 'rarest'}).sort('threshold'))
    rarethresholds = np.array([0] + [r['threshold'] for r in rareresults])
    rarediscard = np.array([0] + [1 - r['mean_remaining'] / r['total'] for r in rareresults])
    #plt.plot(rarethresholds, rarediscard, 'Dw', markeredgewidth=1)
    xfit = np.linspace(0, 1, 1000)
    spl = interpolate.UnivariateSpline(rarethresholds, rarediscard, s=0)
    plt.plot(xfit, spl(xfit), linewidth=1.5, color='r', label='Rarest')
    #
    countresults = list(db.profile.constraint.find({'name': 'counts'}).sort('threshold'))
    countthresholds = np.array([0] + [r['threshold'] for r in countresults])
    countdiscard = np.array([0] + [1 - r['mean_remaining'] / r['total'] for r in countresults])
    #plt.plot(countthresholds, countdiscard, 'sw', markeredgewidth=1)
    xfit = np.linspace(0, 1, 1000)
    spl = interpolate.UnivariateSpline(countthresholds, countdiscard, s=0)
    plt.plot(xfit, spl(xfit), linewidth=1.5, color='g', label='Counts')
    #
    reqresults = list(db.profile.constraint.find({'name': 'reqbits'}).sort('threshold'))
    reqthresholds = np.array([0] + [r['threshold'] for r in reqresults])
    reqdiscard = np.array([0] + [1 - r['mean_remaining'] / r['total'] for r in reqresults])
    #plt.plot(reqthresholds, reqdiscard, '^w', markeredgewidth=1)
    xfit = np.linspace(0, 1, 1000)
    spl = interpolate.UnivariateSpline(reqthresholds, reqdiscard, s=0)
    plt.plot(xfit, spl(xfit), linewidth=1.5, color='b', label='Required')

    plt.ylabel('Discard fraction')
    plt.xlabel('Similarity threshold')
    plt.axis([0, 1, 0, 1])
    plt.grid(True)
    #plt.show()
    #fig = plt.gcf()
    #fig.set_size_inches(10, 6)
    plt.legend(loc='upper left', numpoints=1)
    plt.savefig('img/constraints.svg', dpi=100)  # facecolor='#f2f2f2'


def plot_folding(db):
    unfresults = list(db.profile.folding.find({'length': None}).sort('threshold'))
    unfthresholds = np.array([0] + [r['threshold'] for r in unfresults])
    unfdiscard = np.array([0] + [1 - r['mean_remaining'] / r['total'] for r in unfresults])
    #plt.plot(unfthresholds, unfdiscard, 'ow', markeredgewidth=1)
    xfit = np.linspace(0, 1, 1000)
    spl = interpolate.UnivariateSpline(unfthresholds, unfdiscard, s=0)
    plt.plot(xfit, spl(xfit), linewidth=1.5, color='k', label='Unfolded')

    f2048results = list(db.profile.folding.find({'length': 2048}).sort('threshold'))
    f2048thresholds = np.array([0] + [r['threshold'] for r in f2048results])
    f2048discard = np.array([0] + [1 - r['mean_remaining'] / r['total'] for r in f2048results])
    #plt.plot(f2048thresholds, f2048discard, 'Dw', markeredgewidth=1)
    xfit = np.linspace(0, 1, 1000)
    spl = interpolate.UnivariateSpline(f2048thresholds, f2048discard, s=0)
    plt.plot(xfit, spl(xfit), linewidth=1.5, color='r', label='2048 bits')

    # f1024results = list(db.profile.folding.find({'length': 1024}).sort('threshold'))
    # f1024thresholds = np.array([0] + [r['threshold'] for r in f1024results])
    # f1024discard = np.array([0] + [1 - r['mean_remaining'] / r['total'] for r in f1024results])
    # #plt.plot(f1024thresholds, f1024discard, 'sw', markeredgewidth=1)
    # xfit = np.linspace(0, 1, 1000)
    # spl = interpolate.UnivariateSpline(f1024thresholds, f1024discard, s=0)
    # plt.plot(xfit, spl(xfit), linewidth=1.5, color='g', label='1024 bits')
    #
    f512results = list(db.profile.folding.find({'length': 512}).sort('threshold'))
    f512thresholds = np.array([0] + [r['threshold'] for r in f512results])
    f512discard = np.array([0] + [1 - r['mean_remaining'] / r['total'] for r in f512results])
    #plt.plot(f512thresholds, f512discard, '^w', markeredgewidth=1)
    xfit = np.linspace(0, 1, 1000)
    spl = interpolate.UnivariateSpline(f512thresholds, f512discard, s=0)
    plt.plot(xfit, spl(xfit), linewidth=1.5, color='b', label='512 bits')

    plt.ylabel('Discard fraction')
    plt.xlabel('Similarity threshold')
    plt.axis([0, 1, 0, 1])
    plt.grid(True)
    #plt.show()
    #fig = plt.gcf()
    #fig.set_size_inches(10, 6)
    plt.legend(loc='upper left', numpoints=1)
    plt.savefig('img/folding.svg', dpi=100)  # facecolor='#f2f2f2'


def plot_radius(db):
    r2results = list(db.profile.radius.find({'radius': 2}).sort('threshold'))
    r2thresholds = np.array([0] + [r['threshold'] for r in r2results])
    r2discard = np.array([0] + [1 - r['mean_remaining'] / r['total'] for r in r2results])
    #plt.plot(r2thresholds, r2discard, 'ow', markeredgewidth=1)
    xfit = np.linspace(0, 1, 1000)
    spl = interpolate.UnivariateSpline(r2thresholds, r2discard, s=0)
    plt.plot(xfit, spl(xfit), linewidth=1.5, color='k', label='Morgan radius 2')

    r3results = list(db.profile.radius.find({'radius': 3}).sort('threshold'))
    r3thresholds = np.array([0] + [r['threshold'] for r in r3results])
    r3discard = np.array([0] + [1 - r['mean_remaining'] / r['total'] for r in r3results])
    #plt.plot(r3thresholds, r3discard, 'Dw', markeredgewidth=1)
    xfit = np.linspace(0, 1, 1000)
    spl = interpolate.UnivariateSpline(r3thresholds, r3discard, s=0)
    plt.plot(xfit, spl(xfit), linewidth=1.5, color='r', label='Morgan radius 3')

    r4results = list(db.profile.radius.find({'radius': 4}).sort('threshold'))
    r4thresholds = np.array([0] + [r['threshold'] for r in r4results])
    r4discard = np.array([0] + [1 - r['mean_remaining'] / r['total'] for r in r4results])
    #plt.plot(r4thresholds, r4discard, '^w', markeredgewidth=1)
    xfit = np.linspace(0, 1, 1000)
    spl = interpolate.UnivariateSpline(r4thresholds, r4discard, s=0)
    plt.plot(xfit, spl(xfit), linewidth=1.5, color='b', label='Morgan radius 4')

    plt.ylabel('Discard fraction')
    plt.xlabel('Similarity threshold')
    plt.axis([0, 1, 0, 1])
    plt.grid(True)
    #plt.show()
    #fig = plt.gcf()
    #fig.set_size_inches(10, 6)
    plt.legend(loc='lower right', numpoints=1)
    plt.savefig('img/radius.svg', dpi=100)  # facecolor='#f2f2f2'


def plot_radius_hist(db):

    response = db.chembl.m2.aggregate([{'$group': {'_id': '$count', 'total': {'$sum': 1}}}])
    nos = [r['_id']-0.5 for r in response['result']]
    counts = [r['total'] for r in response['result']]
    plt.bar(nos, counts, width=1, color='b', alpha=0.5, edgecolor='b', label='Morgan radius 2')

    response = db.chembl.m3.aggregate([{'$group': {'_id': '$count', 'total': {'$sum': 1}}}])
    nos = [r['_id']-0.5 for r in response['result']]
    counts = [r['total'] for r in response['result']]
    plt.bar(nos, counts, width=1, color='r', alpha=0.5, edgecolor='r', label='Morgan radius 3')

    response = db.chembl.m4.aggregate([{'$group': {'_id': '$count', 'total': {'$sum': 1}}}])
    nos = [r['_id']-0.5 for r in response['result']]
    counts = [r['total'] for r in response['result']]
    plt.bar(nos, counts, width=1, color='k', alpha=0.5, edgecolor='k', label='Morgan radius 4')

    plt.ylabel('Number of molecules')
    plt.xlabel('Number of 1-bits in fingerprint')
    plt.axis([0, 200, 0, 45000])
    #plt.show()
    fig = plt.gcf()
    fig.set_size_inches(10, 6)
    plt.legend(loc='upper right')
    plt.savefig('img/rhist.svg', dpi=100)

    # nos = []
    # for molecule in db.molecules.find():
    #     nos.append(molecule['mfp2']['count'])
    # print len(nos)
    # nos = np.array(nos)
    # print 'Mean: {}'.format(np.mean(nos))
    # print 'Median: {}'.format(np.median(nos))
    #plt.hist(nos, bins=np.arange(min(nos)-0.5, max(nos)+0.5))
    #plt.show()




# def plot_postgres():
#     # mongodb unfolded aggregation
#     thresholds = np.array([0.95, 0.9, 0.85, 0.8, 0.75, 0.7])
#     times = np.array([6.16252422333, 19.4190740585, 60.7690811157, 166.830062866, 433.220028877, 966.069459915])
#     # mongodb 2048 aggregation
#     #thresholds2 = [0.95, 0.9, 0.8]
#     #times2 = np.array([277.276039124, 829.626083374, 3517.94052124])
#     # mongodb 512 aggregation
#     #thresholds3 = [0.95, 0.9]
#     #times3 = [1565.83452225, 4422.25348949]
#
#     # postgres unfolded
#     thresholds4 = np.array([0.95, 0.9, 0.85, 0.8, 0.75, 0.7])
#     times4 = np.array([57.6151609421, 70.2495574951, 95.799446106, 217.606544495, 509.948015213, 1000.04947186])
#     # postgres 2048
#     thresholds5 = np.array([0.95, 0.9, 0.85, 0.8, 0.75, 0.7])
#     times5 = np.array([77.0859718323, 88.0289077759, 94.8474407196, 112.540483475, 174.414634705, 416.743516922])
#     # postgres 512
#     thresholds6 = np.array([0.95, 0.9, 0.85, 0.8, 0.75, 0.7])
#     times6 = np.array([499.438047409, 670.236587524, 695.670485497, 701.292991638, 733.651995659, 763.621449471])
#
#     popt, pcov = curve_fit(fitfunc, thresholds, times)
#     plt.plot(thresholds6, times6, color='grey', marker='o', markeredgewidth=0, markerfacecolor='k', markersize=7, linewidth=1, linestyle='dashed', label='PostgreSQL 512')
#     plt.plot(thresholds5, times5, color='grey', marker='s', markeredgewidth=0, markerfacecolor='k', markersize=7, linewidth=1, linestyle='dashed', label='PostgreSQL 2048')
#     plt.plot(thresholds4, times4, color='grey', marker='^', markeredgewidth=0, markerfacecolor='k', markersize=7, linewidth=1, linestyle='dashed', label='PostgreSQL unfolded')
#     plt.plot(thresholds, times, 'ow', markeredgewidth=1, markersize=7, label='MongoDB unfolded')
#     fitpoints = np.arange(0.68, 0.97, 0.01)
#     plt.plot(fitpoints, fitfunc(fitpoints, *popt), linewidth=2)
#     plt.ylabel('Median query time (ms)')
#     plt.xlabel('Similarity threshold')
#     plt.axis([0.68, 1, 0, 1100])
#     plt.legend(loc='upper right', numpoints=1)
#     #plt.show()
#     fig = plt.gcf()
#     fig.set_size_inches(10, 6)
#     plt.savefig('img/postgres.svg', dpi=100, facecolor='#f2f2f2')
#
