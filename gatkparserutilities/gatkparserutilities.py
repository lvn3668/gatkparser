# Author: Lalitha Viswanathan
# GATK Parser Utilities
# Affiliation: Stanford Health Care
########################################################
import json
from typing import Any

import depthofcoverage as depthcoverage


########################################################
# check if GATK O/P file is a valid JSON
def is_json(myjson: json) -> bool:
    try:
        json_object: json = json.loads(myjson)
    except ValueError as e:
        print('invalid json: %s' % e)
        return False
    return True


########################################################
# Parse summary statistics at (mbq = 0)
def depthofcoverage0(filename: str) -> dict[str, dict[str, dict]]:
    results = {'depthofcoverage0': depthcoverage.depthofcoveragesamplesummaryandsampleintervalsummary(filename)}
    if (results):
        return results
    else:
        raise Exception('Summary statistics at mbq0 not parsed correctly')


########################################################
# mean base quality at that position = 10)
def depthofcoverage10(filename: str) -> dict[str, dict[str, dict]]:
    results = {'depthofcoverage10': depthcoverage.depthofcoveragesamplesummaryandsampleintervalsummary(filename)}
    if results:
        return results
    else:
        raise Exception('Depth of coverage at mbq10 not parsed correctly')


########################################################
# mean base quality = 20)
def depthofcoverage20(filename: str) -> dict[str, dict[str, dict]]:
    depthofcoverage: dict[str, dict[str, dict]] = dict[str, dict[str, dict]]
    depthofcoverage = {
        'depthofcoverage20': depthcoverage.depthofcoveragesamplesummaryandsampleintervalsummary(filename)}
    if depthofcoverage:
        return depthofcoverage
    else:
        raise Exception('Depth of coverage at mbq20 not parsed correctly')


########################################################

# 30% of the reads span the interval
def depthofcoverage30(filename) -> dict[str, dict[str, dict]]:
    results: dict[str, dict[str, dict]] = dict[str, dict[str, dict]]
    results = {'depthofcoverage30': depthcoverage.depthofcoveragesamplesummaryandsampleintervalsummary(filename)}
    if results:
        return results
    else:
        raise Exception('Depth of coverage at mbq30 not parsed correctly')
########################################################
