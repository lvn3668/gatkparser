# Author: Lalitha Viswanathan
# Calculate Depth of coverage GATK Utilities
# Affiliation: Stanford Health Care
from itertools import chain
from typing import Union

import numpy as np


########################################################
def depthofcoveragesamplesummaryandsampleintervalsummary(filename: str) -> dict[str, dict]:
    """

    :rtype: object
    """
    resultsdata: dict[str, Union[list, list[str]]] = {
        'header': list[str],
        'rows': list[str]
    }
    depthofcoverage: dict[str, dict] = {
        'coverage10x': {},
        'meancoverage': {},
        'granularQ1': {},
        'mediancoverage': {},
        'granularQ3': {},
        'coverage10xbins': {},
        'meancoveragebins': {},
        'mediancoveragebins': {}}
    with open(filename) as f:
        line = f.readline()
        resultsdata['header'] = line.strip().split()
        while line:
            resultsdata['rows'].append(line.rstrip().split())
            try:
                line = f.readline().rstrip()
                if line:
                    vals: list[str] = line.strip().split('\t')
                    depthofcoverage['coverage10x'][vals[0]] = float(vals[8])
                    depthofcoverage['meancoverage'][vals[0]] = float(vals[4])
                    if '>' not in vals[5]:
                        depthofcoverage['granularQ1'][vals[0]] = float(vals[5])
                    if '>' not in vals[6]:
                        depthofcoverage['mediancoverage'][vals[0]] = float(vals[6])
                    if '>' not in vals[7]:
                        depthofcoverage['granularQ3'][vals[0]] = float(vals[7])
            except StopIteration:
                break
    ########################################################
    # Binning for 10xcoverage bins
    unique, counts = np.unique(depthofcoverage['coverage10x'].values(), return_counts=True)
    unique_and_counts = dict(zip(unique, counts))
    data = list(chain.from_iterable([k] * v for k, v in unique_and_counts.items()))
    bins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130,
            140, 150, 160, 170, 180, 190, 200]
    hist = np.histogram(data, bins=bins)[0]
    for r in zip(bins, bins[1:], hist):
        depthofcoverage['coverage10xbins'][str(r[0]) + "-" + str(r[1])] = r[2]
    ########################################################
    # Binning for mean coverage
    unique, counts = np.unique(depthofcoverage['meancoverage'].values(), return_counts=True)
    unique_and_counts = dict(zip(unique, counts))
    data = list(chain.from_iterable([k] * v for k, v in unique_and_counts.items()))
    bins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130,
            140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250]
    hist = np.histogram(data, bins=bins)[0]
    for r in zip(bins, bins[1:], hist):
        depthofcoverage['meancoveragebins'][str(r[0]) + "-" + str(r[1])] = r[2]
    ########################################################
    # Binning for median coverage
    unique, counts = np.unique(depthofcoverage['mediancoverage'].values(), return_counts=True)
    unique_and_counts = dict(zip(unique, counts))
    data = list(chain.from_iterable([k] * v for k, v in unique_and_counts.items()))
    bins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130,
            140, 150, 160, 170,
            180, 190, 200, 210, 220, 230, 240, 250]
    hist = np.histogram(data, bins=bins)[0]
    for r in zip(bins, bins[1:], hist):
        depthofcoverage['mediancoveragebins'][str(r[0]) + "-" + str(r[1])] = r[2]

    if depthofcoverage['mediancoveragebins'] and depthofcoverage['meancoveragebins'] \
            and depthofcoverage['coverage10xbins']:
        return depthofcoverage
    else:
        raise Exception('Depth of Coverage sample interval summary not parsed correctly')

########################################################
