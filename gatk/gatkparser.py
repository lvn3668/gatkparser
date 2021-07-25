#!/usr/bin/env python3
# Author: Lalitha Viswanathan
# GATK Results Parser
#
import re, json
from argparse import ArgumentParser
from typing import Union, Optional, Any

import numpy as np
from itertools import chain
# import gatkparserutilities as gatkp

########################################################
from gatkparserutilities.gatkparserutilities import depthofcoverage30, depthofcoverage0, \
    depthofcoverage20, \
    depthofcoverage10


def generic_parser(filename, tablename) -> object:
    """

    :param filename:
    :param tablename:
    :return:
    :rtype: object
    """
    results = None

    with open(filename) as f:
        for line in f:
            if line.startswith('#:GATKTable:%s' % tablename):

                # Capture the table name and description from the second line.
                m = re.match(r'^#:GATKTable:(.*):(.*)$', line)
                if not m:
                    continue
                assert tablename == m.groups()[0]
                tablename = m.groups()[0]
                tabledescription: str = m.groups()[1]

                # Capture the column header info
                nextrow: str = f.readline().rstrip()
                tableheader: list[str] = re.split('\s+', nextrow)
                tableheader.pop(0)  # Drop the table name

                # Capture data from all rows that start with the current tablename
                rows = []
                nextrow: str = f.readline().rstrip()
                while nextrow.startswith(tablename):
                    rowdata: list[str] = re.split('\s+', nextrow)
                    rowdata.pop(0)  # Drop the table name from the row
                    rows.append(rowdata)  # Save data to list of rows
                    nextrow: str = f.readline().rstrip()

                results = {
                    'description': tabledescription,
                    'header': tableheader,
                    'rows': rows
                }
    return results


########################################################

########################################################
def depthofcoverageparser(gatkoutputfilename, parsertype):
    """

    :param gatkoutputfilename:
    :param parsertype:
    :return:
    """
    parsefxn: str = None
    if parsertype == 'depthofcoveragembq20_sample_summary':
        parsefxn = 'depthofcoverage20'
    if parsertype == 'depthofcoveragembq20_sample_interval_summary':
        parsefxn = 'depthofcoverage20'
    depthofcoverage: Any = parsefxn(gatkoutputfilename)
    if depthofcoverage:
        # ['coverage10x', 'granularQ1', 'mediancoverage', 'coverage10xbins',
        # 'meancoverage', 'mediancoveragebins', 'granularQ3', 'meancoveragebins']
        del depthofcoverage['depthofcoverage20']['coverage10x']
        del depthofcoverage['depthofcoverage20']['granularQ1']
        del depthofcoverage['depthofcoverage20']['mediancoverage']
        del depthofcoverage['depthofcoverage20']['meancoverage']
        del depthofcoverage['depthofcoverage20']['granularQ3']
        return json.dumps(depthofcoverage)
    else:
        raise Exception('Depth of Coverage sample interval summary not parsed correctly')


########################################################
if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('gatkoutputfilename')
    parser.add_argument('sampleid')
    parser.add_argument('runid')
    parser.add_argument('parsertype')
    args = parser.parse_args()
    parsefxn = None

    if args.parsertype == 'depthofcoveragembq0_sample_summary':
        parsefxn = depthofcoverage0
    if args.parsertype == 'depthofcoveragembq0_sample_interval_summary':
        parsefxn = depthofcoverage0
    if args.parsertype == 'depthofcoveragembq10_sample_summary':
        parsefxn = depthofcoverage10
    if args.parsertype == 'depthofcoveragembq10_sample_interval_summary':
        parsefxn = depthofcoverage10
    if args.parsertype == 'depthofcoveragembq20_sample_summary':
        parsefxn = depthofcoverage20
    if args.parsertype == 'depthofcoveragembq20_sample_interval_summary':
        parsefxn = depthofcoverage20
    if args.parsertype == 'depthofcoveragembq30_sample_summary':
        parsefxn = depthofcoverage30
    if args.parsertype == 'depthofcoveragembq30_sample_interval_summary':
        parsefxn = depthofcoverage30
    if args.parsertype == 'varianteval':
        parsefxn = varianteval
    depthofcoverage: dict = parsefxn(args.gatkoutputfilename)
    # print "####", depthofcoverage['depthofcoverage20'].keys()
    print(" 10x coverage bins")
    # pprint(depthofcoverage['depthofcoverage20']['coverage10xbins'])
    print(" Mean coverage bins")
    # pprint(depthofcoverage['depthofcoverage20']['meancoveragebins'])
    print(" Median coverage bins")
    if depthofcoverage:
        jsonfilename: str = args.sampleid + "_" + args.runid + "_" + args.parsertype + "_subset_json.out"
        open(jsonfilename, 'w').close()
        with open(jsonfilename, 'w') as f:
            f.write('%s	%s	%s	%s\n' % (args.sampleid, args.runid, 'GATK', json.dumps(depthofcoverage)))

#############################################################################
