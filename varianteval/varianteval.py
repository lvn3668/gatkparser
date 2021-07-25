# Author Lalitha Viswanathan
# Utility for variant eval (GATK parsing results)
from typing import Optional, Union
import re
from gatk.gatkparser import generic_parser


def varianteval_CountVariants(filename: str) -> dict:
    """

    :rtype: object
    :param filename:
    :return:
    """
    tablename: str = 'CountVariants'
    resultsname: str = 'varianteval_CountVariants'

    fullresults: Optional[dict[str, Union[str, list[str], list[list[str]]]]] = generic_parser(filename, tablename)

    columns = []
    # Get a subset of the results; only SNP and Indel rows, and select columns
    col_varianttype: int = fullresults['header'].index('VariantType')
    columns.append(col_varianttype)
    columns.append(fullresults['header'].index('nProcessedLoci'))
    columns.append(fullresults['header'].index('nVariantLoci'))
    columns.append(fullresults['header'].index('nSNPs'))
    columns.append(fullresults['header'].index('nMNPs'))
    columns.append(fullresults['header'].index('nInsertions'))
    columns.append(fullresults['header'].index('nDeletions'))
    columns.append(fullresults['header'].index('nComplex'))
    columns.append(fullresults['header'].index('nNoCalls'))
    columns.append(fullresults['header'].index('nHets'))
    columns.append(fullresults['header'].index('nHomVar'))
    columns.append(fullresults['header'].index('hetHomRatio'))

    rownums: list[int] = []
    for i in range(0, len(fullresults['rows'])):
        if fullresults['rows'][i][col_varianttype] == 'SNP':
            rownums.append(i)
        if fullresults['rows'][i][col_varianttype] == 'INDEL':
            rownums.append(i)

    header: list[str] = []
    rows: list[list[str]] = []
    for col in columns:
        header.append(fullresults['header'][col])
    for rownum in rownums:
        row: list[str] = []
        for col in columns:
            row.append(fullresults['rows'][rownum][col])
        rows.append(row)

    results = {resultsname: {}}
    results[resultsname]['header'] = header
    results[resultsname]['rows'] = rows

    return results


########################################################
def varianteval_dbsnp_CountVariants(filename: str) -> dict:
    """

    :param filename:
    :return:
    """
    tablename: str = 'CountVariants'
    resultsname: str = 'varianteval_dbsnp_CountVariants'

    fullresults: str = generic_parser(filename, tablename)
    columns: list[int] = []
    # Get a subset of the results; only SNP and Indel rows, and select columns
    col_novelty: int = fullresults['header'].index('Novelty')
    columns.append(col_novelty)
    columns.append(fullresults['header'].index('nHets'))
    columns.append(fullresults['header'].index('nHomVar'))
    columns.append(fullresults['header'].index('hetHomRatio'))

    rownums: list[int] = []
    for i in range(0, len(fullresults['rows'])):
        if fullresults['rows'][i][col_novelty] == 'all':
            rownums.append(i)

    header: list[str] = []
    rows: list[list[str]] = []
    for col in columns:
        header.append(fullresults['header'][col])
    for rownum in rownums:
        row = []
        for col in columns:
            row.append(fullresults['rows'][rownum][col])
        rows.append(row)

    results: dict[str, dict] = {resultsname: {}}
    results[resultsname]['header'] = header
    results[resultsname]['rows'] = rows

    return results


########################################################
def varianteval_TiTvVariantEvaluator(filename: str) -> dict:
    """

    :param filename:
    :return:
    """
    tablename: str = 'TiTvVariantEvaluator'
    resultsname: str = 'varianteval_TiTvVariantEvaluator'

    results: dict[str, object] = {resultsname: generic_parser(filename, tablename)}

    # Now we want to pull out some key results that can be displayed without the whole table.
    # Find the row with SNP's.

    col_varianttype: int = results[resultsname]['header'].index('VariantType')
    col_nTi: int = results[resultsname]['header'].index('nTi')
    col_nTv: int = results[resultsname]['header'].index('nTv')
    col_TiTvRatio: int = results[resultsname]['header'].index('tiTvRatio')

    for row in results[resultsname]['rows']:
        if row[col_varianttype] == 'SNP':
            break

    results[resultsname]['variantType'] = row[col_varianttype]
    results[resultsname]['nTi'] = row[col_nTi]
    results[resultsname]['nTv'] = row[col_nTv]
    results[resultsname]['TiTvRatio'] = row[col_TiTvRatio]

    return results


########################################################
def varianteval_dbsnp_TiTvVariantEvaluator(filename: str) -> dict:
    """

    :param filename:
    :return:
    """
    tablename: str = 'TiTvVariantEvaluator'
    resultsname: str = 'varianteval_dbsnp_TiTvVariantEvaluator'

    fullresults: Optional[dict[str, Union[str, list[str], list[list[str]]]]] = generic_parser(filename, tablename)

    # Now we want to pull out some key results that can be displayed without the whole table.

    columns: list[int] = []
    col_novelty: int = fullresults['header'].index('Novelty')
    columns.append(col_novelty)
    columns.append(fullresults['header'].index('nTi'))
    columns.append(fullresults['header'].index('nTv'))
    columns.append(fullresults['header'].index('tiTvRatio'))

    rownums: list[int] = []
    for i in range(len(fullresults['rows'])):
        if fullresults['rows'][i][col_novelty] == 'all':
            rownums.append(i)

    header: list[str] = []
    rows: list[list[str]] = []
    for col in columns:
        header.append(fullresults['header'][col])

    for rownum in rownums:
        row: list[str] = []
        for col in columns:
            row.append(fullresults['rows'][rownum][col])
        rows.append(row)

    results = {resultsname: {}}
    results[resultsname]['header'] = header
    results[resultsname]['rows'] = rows

    return results


########################################################
def varianteval(filename: str) -> dict:
    """

    :param filename:
    :return:
    """
    resultstables: dict[str, dict[str, Union[str, list[str], list[list[str]]]]] = {}

    with open(filename) as f:
        for line in f:
            if line.startswith('#:GATKTable:'):
                # New table found
                # Ignore the first line with format information.
                # Capture the table name and description from the second line.
                nextrow: str = f.readline()
                m = re.match(r'^#:GATKTable:(.*):(.*)$', nextrow)
                if not m:
                    continue
                tablename = m.groups()[0]
                tabledescription = m.groups()[1]

                # Capture the column header info
                nextrow: str = f.readline().rstrip()
                tableheader: list[str] = re.split('\s+', nextrow)
                tableheader.pop(0)  # Drop the table name

                # Capture data from all rows that start with the current tablename
                rows: list[list[str]] = []
                nextrow: str = f.readline().rstrip()
                while nextrow.startswith(tablename):
                    rowdata: list[str] = re.split('\s+', nextrow)
                    rowdata.pop(0)  # Drop the table name from the row
                    rows.append(rowdata)  # Save data to list of rows
                    nextrow: str = f.readline().rstrip()

                resultstables[tablename] = {
                    'description': tabledescription,
                    'header': tableheader,
                    'rows': rows
                }
    results = {'varianteval': resultstables}
    return results

