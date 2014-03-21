import os
import time
import sys
import math
sys.path.append("./py_modules/")

#import psyco
#psyco.full()

import col
import col_util
import stats


################################################################
## definition of string constants
################################################################
OUTPUT_DIR  = './output/'
REPORT_EXT  = '_v.1.0_report.txt'
OUT_XML_EXT = '_est.xml'



def ln(rawString):
    
    mantissa, exponent = rawString.split('e')

    return math.log( float(mantissa) ) + float(exponent) * math.log(10)


def likelihoodRatioTest(rawString_n, rawString_d, df):

    test_stat = -2 * (ln(rawString_n) - ln(rawString_d))

    print round(test_stat, 2)


    return stats.chisqprob(test_stat, df)


def calculateLRT(tupleList):

    for tuple in tupleList:
        modelName, n, d, df = tuple
        print modelName, likelihoodRatioTest(n, d, df)


def countFreeParameters(resultsEntry):

    df = 0
    for name, value, comment in resultsEntry['parameters']:
        if name.find('freq') != -1:
            df += len(value) - 1 
        elif comment != "final":
            df += 1

    return df


def calcLRT(LRT_FILE, results):

    file = open(LRT_FILE, 'w')

    for i in range(1, len(results)):
        df_im1         = countFreeParameters(results[i-1])
        df_i           = countFreeParameters(results[i])
        likelihood_im1 = results[i-1]['likelihood']
        likelihood_i   = results[i]['likelihood']

        test_stat      = -2 * (ln(likelihood_i) - ln(likelihood_im1))
        delta_df       = df_im1 - df_i
        test_prob      = stats.chisqprob(test_stat, delta_df)

        # print entries for latex table
        file.write(results[i]['model_name'] + " " + "&" + " " + "some description" + " " + "&" + " " + str(df_i) + " " + "&" + " " +
                   str(delta_df) + " " + "&" + " " + likelihood_i + " " + "&" + " " + str(round(test_stat, 2)) + " " + "&" + " " +
                   str(round(test_prob, 4)) + " " + "\\\\" + "\n")

    print "LTR saved to file '" + LRT_FILE + "'"


def createXMLfile(templateFileName, outFileName, taxonList, treeString):
    """ Look for INSERT_TAG in templateFile, and inserst taxa element and tree element. """

    INSERT_TAG = 'ZZZZZ'


    # make taxaList
    taxaList = []
    taxaList.append('      <taxa id="TA_all">\n')
    for taxon in taxonList:
        taxaList.append('        <taxon id="' + taxon.strip() + '"/>\n')
    taxaList.append('      </taxa>\n')

    # make taxaList
    treeList = []    
    treeList.append('          <tree id="T_all">\n')
    treeList.append('            <newick>\n')
    treeList.append(treeString)
    treeList.append('            </newick>\n')
    treeList.append('          </tree>\n')

    # make xmlList
    inFile = open(templateFileName, 'r')
    xmlList = []
    for line in inFile.readlines():
        xmlList.append(line)
        if (line.find(INSERT_TAG) != -1):
            xmlList += taxaList
            xmlList.append('\n')
            xmlList += treeList
            
    # make out file   
    outFile = open(outFileName, 'w')
    for line in xmlList:
        outFile.write(line)
    outFile.close()


def extractTaxons(colSeq):
    """ Return list of taxon names """
    
    columnIndexes      = col_util.symbolsColumnIndexes(colSeq)
    (names, sequences) = col_util.namesAndSequences(colSeq, columnIndexes)
    return names


def getTaxonListAndTreeString(baseName, testDataFileName):
    """ Return taxon list and tree string touple. Expects tree file with tree suffix """
    TREE_PREFIX = './tree/'
    TREE_SUFFIX = '.tree'

    treeString = open(TREE_PREFIX + baseName + TREE_SUFFIX).readlines()[0]

    colSeqList = col.ColumnSeqList()
    colSeqList.read(testDataFileName)
    for entry in colSeqList:
        taxonList     = extractTaxons(entry)
        return (taxonList, treeString)
    
    
def InsertAtTag(list, insertList, tag):
    """ Insert list into insertList after line(s) with tag """
    for i in range(len(list)):
        if list[i].find(tag) != -1:
            list[i:i] = insertList
            return list


def replaceInFile(inFile, outFile, original, replacement):
    """ Open inFile, relace all instances of original with replacement and save to outFile """

    lines  = open(inFile, 'r').readlines()
    output = open(outFile, 'w')

    for line in lines:
        output.write( line.replace(original, replacement) )

    output.close()


def saveEquiFreq(xmlFile, equiFreqFile):
    tag = 'SM_loop.freq'
    lines   = open(xmlFile, 'r').readlines()
    freqLines = []

    for i in range(len(lines)):
        if lines[i].find(tag) != -1:
            freqLines = lines[i:i+6]

    #check
    if not freqLines:
        raise Exception, "no freqlines found"

    for i in range(4):
        freqLines[i+1] = freqLines[i+1].replace('/>', 'final = "" />')

    outFile = open(equiFreqFile, 'w')
    for line in freqLines:
        outFile.write(line)


def insertEquiFreq(xmlFile, equiFreqFile):
    # inser equiFreqs
    xmlLines   = open(xmlFile).readlines()
    freqLines  = open(equiFreqFile).readlines()
    InsertAtTag(xmlLines, freqLines, 'UUUUU')
    outFile = open(xmlFile, 'w')

    for line in xmlLines:
        outFile.write(line)

def extractLikelihood(reportFile):
    file = open(reportFile, 'r')

    for line in file:
        if line.find('Total inside prob') != -1:
            return line.split()[-1]

    raise Exception, "No likelihood found in report file:'" + reportFile + "'."


def getNameAndValue(xmlLine):
    """ Parse parameter line in xml format and return name and value attributes """
    
    nameBegin = xmlLine.index('name')
    name      = xmlLine[nameBegin:].split('"')[1]

    valueBegin = xmlLine.index('value')
    value      = xmlLine[valueBegin:].split('"')[1]

    comment    = ''
    if xmlLine.find('final') != -1:
        comment = 'final'

    return (name, value, comment)


def getParameters(xmlFile):
    """ Return list of all shared parameters """
    results  = []
    xmlLines = open(xmlFile, 'r').readlines()

    SHARED_PARAM_TAG      = 'Definition of shared parameters'
    END_COMPUND_PARAM_TAG = '</compoundParameter>'
    PARAM_TAG             = '<parameter'

    foundSharedParam        = False
    i = 0
    while not foundSharedParam and i < len(xmlLines):
        if xmlLines[i].find(SHARED_PARAM_TAG) != -1:
            foundSharedParam = True
        i += 1

    if not foundSharedParam:
        print "getParameters: warning: no shared parameters found in '" + xmlFile + "'."
        return results

    foundEndCompoundParam   = False
    while not foundEndCompoundParam and i < len(xmlLines):
        if xmlLines[i].find(PARAM_TAG) != -1:
            results.append( getNameAndValue(xmlLines[i]) )
        if xmlLines[i].find(END_COMPUND_PARAM_TAG) != -1:
            foundEndCompoundParam = True
        i += 1

    EQUI_FREQ_TAG = '<frequencies id='
    i = 0
    while i < len(xmlLines):
        if xmlLines[i].find(EQUI_FREQ_TAG) != -1:
            results.append( getEquiFreq(xmlLines, i) )
        i += 1
    
    return results


def getEquiFreq(xmlLines, startIndex):
    idBegin = xmlLines[startIndex].index('id')
    id      = xmlLines[startIndex][idBegin:].split('"')[1]

    v = []
    i = startIndex
    while xmlLines[i].find('</frequencies>') == -1:
        if xmlLines[i].find('parameter') != -1:
            name, value, comment = getNameAndValue(xmlLines[i])
            v.append(value)
        i += 1

    return (id, v, '')


def collectResults(baseName, modelName):

    REPORT_FILE = OUTPUT_DIR + baseName + '_' + modelName + REPORT_EXT
    XML_FILE    = OUTPUT_DIR + baseName + '_' + modelName + OUT_XML_EXT

    results = {}
    results['model_name']   = modelName
    results['seq_name']     = baseName
    results['likelihood']   = extractLikelihood(REPORT_FILE)
    results['parameters']   = getParameters(XML_FILE)

    return results


def genericRunModel(baseName, modelName, taxonList, treeString):
    
    TEMPLATE_NAME = 'models/' + modelName + '_template.xml'
    OUT_FILE_NAME = OUTPUT_DIR + baseName + '_' + modelName + '.xml'

    createXMLfile(TEMPLATE_NAME, OUT_FILE_NAME, taxonList, treeString)
    replaceInFile(OUT_FILE_NAME, OUT_FILE_NAME, 'likelihood_test_scfg', baseName + '_' + modelName)
    os.system('./bin/CORS_SCFG '+ OUT_FILE_NAME)

    return collectResults(baseName, modelName)


def saveResults(resultFile, results):
    #write results to file
    file = open(resultFile, 'w')
    for dict_entry in results:
        file.write( "Results for '" + dict_entry['model_name'] + "' using sequence '" + dict_entry['seq_name'] + "'\n")
        file.write( 'likelihood' + ':\t' + dict_entry['likelihood'] + '\n')
        for name, value, comment in dict_entry['parameters']:
            file.write(name + ':\t' + str(value) + '\t\t' + comment + '\n')
        file.write('\n')

    # print to screen
    os.system('cat ' + resultFile)
    
  
def runTestModelHierachy(testData, baseName, modelList):

    TEST_DATA = testData
    BASE_NAME = baseName

    RESULT_FILE = BASE_NAME + '_results.txt'
    ANNO_DATA   = './data/anno.col'
    
    taxonList, treeString = getTaxonListAndTreeString(BASE_NAME, TEST_DATA)
    os.system( 'rm ' + ANNO_DATA + '; cp ' + TEST_DATA + ' ' + ANNO_DATA )
    
    results = []

    for model in modelList:
        results.append( genericRunModel(BASE_NAME, model, taxonList, treeString) )
        saveResults(RESULT_FILE, results)

    LRT_FILE = BASE_NAME + '_LRT.txt'
    calcLRT(LRT_FILE, results)


def main():

    modelList = ['model_full', 'model_start', 'model_1', 'model_2',
                 'model_3', 'model_4', 'model_5', 'model_6',
                 'model_7', 'model_8']


    runTestModelHierachy('HCV_all_type_1a_and_1b_TRAIN_all.col', 'HCV_all_type_1a_and_1b', modelList)


# Run main
if __name__ == '__main__':
    main()
