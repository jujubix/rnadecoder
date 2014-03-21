#! /usr/bin/env python
import os
import copy
import sys

####################################################################################
# path to col and col_util
####################################################################################
sys.path.append("./python_scripts/py_modules/")

import col
import col_util
import testModelHierarchy

####################################################################################
# common string constants
####################################################################################
TRUE_PAIRED_POSITIONS  = 'pairedPositions'
ORIGINAL_POSITIONS     = 'positions'
CODON_MASK             = 'labels type="codonmask"'
PAIRING_MASK           = 'labels type="pairingmask"'
TEST_COL               = './test.col'
TEST_CYK                = './output/test_cyk.col'
POST_PROB              = 'posteriorProbabilities'


def getPerformance(colFileName):
    list       = col.ColumnSeqList()
    list.read(colFileName)

    perf = []
    for entry in list:
         perf.append(col_util.RNAssPerformanceMeasures(entry) )

    return perf


def makeSCFGXMLFile(scfgTemplate, colFileName, scfgName, baseName):
    taxonList, treeString = testModelHierarchy.getTaxonListAndTreeString(baseName, colFileName)
    testModelHierarchy.createXMLfile(scfgTemplate, scfgName, taxonList, treeString)


def makeSingleSeqSCFGXMLFile(scfgTemplate, scfgName, baseName, name):
    treeString = '(' + name + ');'
    taxonList  = [name]
    testModelHierarchy.createXMLfile(scfgTemplate, scfgName, taxonList, treeString)
    

def ComparativePred(scfgTemplate, colFileName, baseName):

    scfgName = baseName + '_tmp.xml'
    makeSCFGXMLFile(scfgTemplate, colFileName, scfgName, baseName)
    os.system('cp -f ' + colFileName + ' ' + TEST_COL)
    os.system("RNA-decoder " + scfgName)


def getSymbolColumnIndex(entry, name):

    for i in range(entry.size() ):
        if entry.columnTag(i).find('symbol') != -1:
            if entry.columnTag(i).find(name) != -1:            
                return i
    else:
        raise Exception, "symbol column with name '" + name + "' not found."



def makePairingMask(colSeq, fixFraction):
    ppList           = colSeq.column(POST_PROB)
    pairingMask      = copy.deepcopy(colSeq.column(PAIRING_MASK) )

    rankedpp         = [x for x in enumerate(ppList)]
    rankedpp.sort(lambda x,y: cmp(x[1], y[1]) )

    L                = len(rankedpp)
    #indexList = [x[0] for x in rankedpp[L - int(fixFraction * L):] ]
    indexList = [x[0] for x in rankedpp if float(x[1]) > fixFraction ]
   
    for i in range(len(pairingMask)):
        if i not in indexList:
            pairingMask[i] = '*'

    return pairingMask



def makePairingMask_2(colSeq, fixFraction):
    ppList           = colSeq.column(POST_PROB)
    pairingMask      = copy.deepcopy(colSeq.column(PAIRING_MASK) )

    rankedpp         = [x for x in enumerate(ppList)]
    rankedpp.sort(lambda x,y: cmp(x[1], y[1]) )

    L                = len(rankedpp)
    #indexList = [x[0] for x in rankedpp[L - int(fixFraction * L):] ]
    indexList = [x[0] for x in rankedpp if float(x[1]) > fixFraction ]
    indexList = [x for x in indexList if not pairingMask[x] == '.']
   
    for i in range(len(pairingMask)):
        if i not in indexList:
            pairingMask[i] = '*'

    print pairingMask

    return pairingMask


def makeSingleSeqColFile(colFileName, fixFraction, name, prefix = './tmp/'):

    list       = col.ColumnSeqList()
    list.read(colFileName)
    
    newList  = col.ColumnSeqList()
    for entry in list:
        newColSeq = col.ColumnSeq(name + '_' + entry.feature('ENTRY') )

        newPairingMask = makePairingMask_2(entry, fixFraction)
        newColSeq.setColumn(PAIRING_MASK, newPairingMask)
        newColSeq.setColumn(CODON_MASK, entry.column(CODON_MASK) )

        newColSeq.setColumn(TRUE_PAIRED_POSITIONS, entry.column(TRUE_PAIRED_POSITIONS) )
        newColSeq.setColumn(ORIGINAL_POSITIONS, entry.column(ORIGINAL_POSITIONS) )

        newColSeq.setColumn(POST_PROB, entry.column(POST_PROB) )

        i      = getSymbolColumnIndex(entry, name)
        tag    = entry.columnTag(i)
        column = entry.columnByIndex(i)

        newColSeq.setColumn(tag, column)

        newList.setEntry(newColSeq)

    fileName = prefix + name + '.col'
    newList.write(fileName)
    return fileName

    
def singleSeqPred(scfgTemplate, singleSeqColFileName, baseName, name, prefix = './tmp/'):

    scfgName = prefix + baseName + '_tmp.xml'
    makeSingleSeqSCFGXMLFile(scfgTemplate, scfgName, baseName, name)
    os.system('cp -f ' + singleSeqColFileName + ' ' + TEST_COL)
    os.system("RNA-decoder " + scfgName)


def twoStepPrediction(scfgTemplate_cmp, scfgTemplate_single, colFileName, baseName, nameList = None, fixFraction = 0.3):

    ComparativePred(scfgTemplate_cmp, colFileName, baseName)

    for name in nameList:
        singleSeqColFileName = makeSingleSeqColFile(TEST_CYK, fixFraction, name, './output/')

        os.system('cp -f ' + TEST_CYK + ' ./output/test_cmp_pd.col' )
        
        singleSeqPred(scfgTemplate_single, singleSeqColFileName, baseName, name, './output/')


def main():
    if len(sys.argv) < 7 or len(sys.argv) > 7:
        raise Exception, "usage: " + sys.argv[0] + " <scfg_cmp_template.xml> <scfg_single_template.xml> <colFileName> <baseName> <name> <trustFraction>"

    scfgTemplate_cmp    = sys.argv[1]
    scfgTemplate_single = sys.argv[2]
    colFileName         = sys.argv[3]
    baseName            = sys.argv[4]
    nameList            = []
    nameList.append(sys.argv[5])
    fixFraction         = float(sys.argv[6])

    twoStepPrediction(scfgTemplate_cmp, scfgTemplate_single, colFileName, baseName, nameList, fixFraction)


#run main
if __name__ == '__main__':
    main()



    
# def singleSeqPred(scfgTemplate, singleSeqColFileName, baseName, name, prefix = './tmp/'):

#     scfgName = prefix + baseName + '_tmp.xml'
#     makeSingleSeqSCFGXMLFile(scfgTemplate, scfgName, baseName, name)
#     os.system('cp -f ' + singleSeqColFileName + ' ' + TEST_COL)
#     os.system("RNA-decoder " + scfgName)


# def twoStepPrediction(scfgTemplate, colFileName, baseName, nameList = None, fixFraction = 0.3):

#     ComparativePred(scfgTemplate, colFileName, baseName)

#     print 'comparative performance sn_s, sp_s, sn_p, sp_p:'
#     print getPerformance(TEST_CYK)

#     for name in nameList:
#         singleSeqColFileName = makeSingleSeqColFile(TEST_CYK, fixFraction, name, './output/')
#         singleSeqPred(scfgTemplate, singleSeqColFileName, baseName, name, './output/')

#         print 'single performance sn_s, sp_s, sn_p, sp_p:'
#         print getPerformance(TEST_CYK)


# def main():
#     if len(sys.argv) < 6 or len(sys.argv) > 6:
#         raise Exception, "usage: " + sys.argv[0] + " <scfg_template.xml> <colFileName> <baseName> <name> <trustFraction>"

#     scfgTemplate = sys.argv[1]
#     colFileName  = sys.argv[2]
#     baseName     = sys.argv[3]
#     nameList     = []
#     nameList.append(sys.argv[4])
#     fixFraction = float(sys.argv[5])

#     twoStepPrediction(scfgTemplate, colFileName, baseName, nameList, fixFraction)
