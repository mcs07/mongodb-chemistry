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


def plot_screen(results, label, style):
    thresholds = np.array([0] + [r['threshold'] for r in results])
    discard = np.array([0] + [1 - r['mean_remaining'] / r['total'] for r in results])
    log.info(thresholds)
    log.info(discard)
    #plt.plot(thresholds, discard, 'ow', markeredgewidth=1)
    xfit = np.linspace(0, 1, 1000)
    spl = interpolate.UnivariateSpline(thresholds, discard, s=0)
    plt.plot(xfit, spl(xfit), style, linewidth=1.5, label=label)


def plot_screening(result_collection):
    """Plot screening ability of different screening methods."""
    #idealresults = list(result_collection.find({'ideal': True}).sort('threshold'))
    allresults = list(result_collection.find({'fp': 'm2', 'counts': True, 'reqbits': True, 'rarest': True}).sort('threshold'))
    countresults = list(result_collection.find({'counts': True, 'reqbits': False, 'rarest': False}).sort('threshold'))
    reqresults = list(result_collection.find({'counts': False, 'reqbits': True, 'rarest': False}).sort('threshold'))
    rareresults = list(result_collection.find({'counts': False, 'reqbits': True, 'rarest': True}).sort('threshold'))
    #plot_screen(idealresults, 'Ideal', 'k--')
    plot_screen(allresults, 'Combined', 'k')
    plot_screen(countresults, 'Counts', 'g')
    plot_screen(reqresults, 'Required', 'b')
    plot_screen(rareresults, 'Rarest', 'r')
    plt.ylabel('Discard fraction')
    plt.xlabel('Similarity threshold')
    plt.axis([0, 1, 0, 1])
    plt.grid(True)
    #plt.show()
    #fig = plt.gcf()
    #fig.set_size_inches(10, 6)
    plt.legend(loc='upper left', numpoints=1)
    plt.savefig('img/screening.svg', dpi=100)  # facecolor='#f2f2f2'


def plot_folding(result_collection):
    unfresults = list(result_collection.find({'fp': 'm2', 'counts': True, 'reqbits': True, 'rarest': True}).sort('threshold'))
    plot_screen(unfresults, 'Unfolded', 'k')
    f2048results = list(result_collection.find({'fp': 'm2l2048', 'counts': True, 'reqbits': True, 'rarest': True}).sort('threshold'))
    plot_screen(f2048results, '2048 bits', 'r')
    f1024results = list(result_collection.find({'fp': 'm2l1024', 'counts': True, 'reqbits': True, 'rarest': True}).sort('threshold'))
    plot_screen(f1024results, '1024 bits', 'g')
    f512results = list(result_collection.find({'fp': 'm2l512', 'counts': True, 'reqbits': True, 'rarest': True}).sort('threshold'))
    plot_screen(f512results, '512 bits', 'b')
    plt.ylabel('Discard fraction')
    plt.xlabel('Similarity threshold')
    plt.axis([0, 1, 0, 1])
    plt.grid(True)
    #plt.show()
    #fig = plt.gcf()
    #fig.set_size_inches(10, 6)
    plt.legend(loc='upper left', numpoints=1)
    plt.savefig('img/folding.svg', dpi=100)  # facecolor='#f2f2f2'


def plot_radius(result_collection):
    r2results = list(result_collection.find({'fp': 'm2', 'counts': True, 'reqbits': True, 'rarest': True}).sort('threshold'))
    plot_screen(r2results, 'Morgan radius 2', 'k')
    r3results = list(result_collection.find({'fp': 'm3', 'counts': True, 'reqbits': True, 'rarest': True}).sort('threshold'))
    plot_screen(r3results, 'Morgan radius 3', 'r')
    r4results = list(result_collection.find({'fp': 'm4', 'counts': True, 'reqbits': True, 'rarest': True}).sort('threshold'))
    plot_screen(r4results, 'Morgan radius 4', 'b')
    plt.ylabel('Discard fraction')
    plt.xlabel('Similarity threshold')
    plt.axis([0, 1, 0, 1])
    plt.grid(True)
    #plt.show()
    #fig = plt.gcf()
    #fig.set_size_inches(10, 6)
    plt.legend(loc='lower right', numpoints=1)
    plt.savefig('img/radius.svg', dpi=100)  # facecolor='#f2f2f2'


# def plot_radius_hist(db):
#
#     response = db.chembl.m2.aggregate([{'$group': {'_id': '$count', 'total': {'$sum': 1}}}])
#     nos = [r['_id']-0.5 for r in response['result']]
#     counts = [r['total'] for r in response['result']]
#     plt.bar(nos, counts, width=1, color='b', alpha=0.5, edgecolor='b', label='Morgan radius 2')
#
#     response = db.chembl.m3.aggregate([{'$group': {'_id': '$count', 'total': {'$sum': 1}}}])
#     nos = [r['_id']-0.5 for r in response['result']]
#     counts = [r['total'] for r in response['result']]
#     plt.bar(nos, counts, width=1, color='r', alpha=0.5, edgecolor='r', label='Morgan radius 3')
#
#     response = db.chembl.m4.aggregate([{'$group': {'_id': '$count', 'total': {'$sum': 1}}}])
#     nos = [r['_id']-0.5 for r in response['result']]
#     counts = [r['total'] for r in response['result']]
#     plt.bar(nos, counts, width=1, color='k', alpha=0.5, edgecolor='k', label='Morgan radius 4')
#
#     plt.ylabel('Number of molecules')
#     plt.xlabel('Number of 1-bits in fingerprint')
#     plt.axis([0, 200, 0, 45000])
#     #plt.show()
#     fig = plt.gcf()
#     fig.set_size_inches(10, 6)
#     plt.legend(loc='upper right')
#     plt.savefig('img/rhist.svg', dpi=100)
#
#     # nos = []
#     # for molecule in db.molecules.find():
#     #     nos.append(molecule['mfp2']['count'])
#     # print len(nos)
#     # nos = np.array(nos)
#     # print 'Mean: {}'.format(np.mean(nos))
#     # print 'Median: {}'.format(np.median(nos))
#     #plt.hist(nos, bins=np.arange(min(nos)-0.5, max(nos)+0.5))
#     #plt.show()


# def plot_benchmarks():
#     # unfolded
#     thresholds = np.array([0.95, 0.9, 0.85, 0.8, 0.75, 0.7])
#     times = np.array([6.16252422333, 19.4190740585, 60.7690811157, 166.830062866, 433.220028877, 966.069459915])
#     # 2048
#     thresholds2 = [0.95, 0.9, 0.8]
#     times2 = np.array([277.276039124, 829.626083374, 3517.94052124])
#     # 512
#     thresholds3 = [0.95, 0.9]
#     times3 = [1565.83452225, 4422.25348949]
#
#
#     popt, pcov = curve_fit(fitfunc, thresholds, times)
#     #plt.plot(thresholds3, times3, color='grey', marker='s', markeredgewidth=0, markerfacecolor='k', markersize=7, linewidth=1, linestyle='dashed', label='512 bit')
#     #plt.plot(thresholds2, times2, color='grey', marker='o', markeredgewidth=0, markerfacecolor='k', markersize=7, linewidth=1, linestyle='dashed', label='2048 bit')
#     plt.plot(thresholds, times, 'ow', markeredgewidth=1, label='unfolded')
#     fitpoints = np.arange(0.68, 0.97, 0.01)
#     plt.plot(fitpoints, fitfunc(fitpoints, *popt), linewidth=2)
#     plt.ylabel('Median query time (ms)')
#     plt.xlabel('Similarity threshold')
#     plt.axis([0.68, 1, 0, 1100])
#     #plt.show()
#     fig = plt.gcf()
#     fig.set_size_inches(10, 6)
#     #plt.legend(loc='upper left', numpoints=1)
#     plt.savefig('img/bench-fold.svg', dpi=100)


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
