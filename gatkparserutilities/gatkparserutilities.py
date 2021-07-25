# Author: Lalitha Viswanathan
# GATK Parser Utilities
########################################################
import json

from gatkparser.depthofcoverage import depthofcoveragesamplesummaryandsampleintervalsummary


def is_json(myjson: json):
    try:
        json_object: json = json.loads(myjson)
    except ValueError as e:
        print('invalid json: %s' % e)
        return False
    return True


def depthofcoverage0(filename: str):
    results = {'depthofcoverage0': depthofcoveragesamplesummaryandsampleintervalsummary(filename)}
    return results


########################################################
def depthofcoverage10(filename: str):
    results = {'depthofcoverage10': depthofcoveragesamplesummaryandsampleintervalsummary(filename)}
    return results


########################################################
def depthofcoverage20(filename: str):
    depthofcoverage: dict[str, dict[str, dict]] = dict[str, dict[str, dict]]
    depthofcoverage = {'depthofcoverage20': depthofcoveragesamplesummaryandsampleintervalsummary(filename)}
    if depthofcoverage:
        return depthofcoverage
    else:
        raise Exception('Depth of coverage mbq20 sample interval summary not parsed correctly')

    ########################################################


def depthofcoverage30(filename):
    results: dict[str, dict[str, dict]] = dict[str, dict[str, dict]]
    results = {'depthofcoverage30': depthofcoveragesamplesummaryandsampleintervalsummary(filename)}
    return results

########################################################
