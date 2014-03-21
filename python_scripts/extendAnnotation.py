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

####################################################################################
# common string constants
####################################################################################
PAIRED_POSITION = 'pairedPosition'
PAIRING_MASK    = 'labels type="pairingmask"'
POST_PROB       = 'posteriorProbabilities'


def removeRNAPair(col, pairCol, j, j_pair):
    pairCol[j]      = '.' 
    pairCol[j_pair] = '.'

    return


def getName(tag):
    begin = tag.find('"') + 1
    end   = tag.find('"', begin)
    return tag[begin:end]
    

def separateSymbolSequences(colSeq):

    cols         = ['pairedPositions', 'positions', 'labels type="pairingmask"', 'pairedPosition', 'posteriorProbabilities']
    
    colSeqList   = col.ColumnSeqList()    

    for i in range( colSeq.columnCount() ):
        if colSeq.columnTag(i).find('symbol') != -1:
            name      = getName(colSeq.columnTag(i))
            subColSeq = col.ColumnSeq(colSeq.feature('ENTRY') + '_' + name)
            subColSeq.setColumn(colSeq.columnTag(i), copy.deepcopy(colSeq.columnByIndex(i)))
            for colName in cols:
                subColSeq.setColumn(colName, copy.deepcopy(colSeq.column(colName)))

            colSeqList.setEntry(subColSeq)

    return colSeqList


def replacePP(colSeq, invert = False):
    pairCol         = colSeq.column(PAIRED_POSITION)
    pairingMask     = colSeq.column(PAIRING_MASK)
    postProb        = colSeq.column(POST_PROB)
    
    for i in range(len(pairCol)):
        if pairingMask[i] == '.' and pairCol[i] != '.':
            pair_i = int(pairCol[i]) - 1  # counting starts at one
            if i > pair_i:
                if invert:
                    postProb[i]      = str(1-float(postProb[i]))
                    postProb[pair_i] = str(1-float(postProb[pair_i]))
                else:                    
                    postProb[i]      = '-1'
                    postProb[pair_i] = '-1'
        if pairingMask[i] != '.' and pairCol[i] == '.':
            if invert:
                postProb[i]      = str(1-float(postProb[i]))
            else:
                postProb[i]      = '-1'

    return

def rectifyPairingMask(colSeq):
    pairCol         = colSeq.column(PAIRED_POSITION)
    pairingMask     = colSeq.column(PAIRING_MASK)

    for i in range(len(pairCol)):
        if pairCol[i] != '.':
            pair_i = int(pairCol[i]) - 1  # counting starts at one
            if i < pair_i:
                pairingMask[i]      = '('
                pairingMask[pair_i] = ')'
        else:
            pairingMask[i]      = '.'

    return


def extendRNAStem(col, pairCol, j, j_pair, min_loop):
    L = len(col)

    k = 1
    while True:
        if j < j_pair and (j - k < 0 or j_pair + k >= L):                              # exceeding seq length
            return
        if abs((j - k) - (j_pair + k)) - 1 < min_loop:        # min loop length
            return
        if (pairCol[j - k] != '.') or (pairCol[j_pair + k] != '.'):
            return 
        if not col_util.RNAPair(col[j - k], col[j_pair + k]):
            return

        pairCol[j - k]      = str(j_pair + k + 1)
        pairCol[j_pair + k] = str(j - k + 1)

        k += 1

    return


def extendAllRNAStems(colSeq, min_loop):
    col = []
    for i in range( colSeq.columnCount() ):
        if colSeq.columnTag(i).find('symbol') != -1:
            col = colSeq.columnByIndex(i)

    pairCol         = colSeq.column(PAIRED_POSITION)
    pairingMask     = colSeq.column(PAIRING_MASK)

    for j in range(len(col)):
        if pairCol[j] != '.':
            j_pair = int(pairCol[j])-1 # counting starts at one
            if not col_util.RNAPair(col[j], col[j_pair]):
                removeRNAPair(col, pairCol, j, j_pair)                  
            else:
                extendRNAStem(col, pairCol, j, j_pair, min_loop)

    return


def extendRNAStemsInColseq(colSeq, min_loop, invertPP):
    colSeqList      = separateSymbolSequences(colSeq)

    for entry in colSeqList:
        extendAllRNAStems(entry, min_loop)
        replacePP(entry, invertPP)
        rectifyPairingMask(entry)
        col_util.addPosition(entry,1)
        
    return colSeqList


def extendRNAStemsInFile(in_file, out_file, min_loop = 4, invertPP=False):
    colSeqList  = col.ColumnSeqList()
    colSeqList.read(in_file)

    split_list = col.ColumnSeqList()

    for entry in colSeqList:
        split_list.add(extendRNAStemsInColseq(entry, min_loop, invertPP))

    split_list.write(out_file)


def main():
    if len(sys.argv) != 3:
        raise Exception, "usage: " + sys.argv[0] + " <infile.col> <outfile.col>"

    inFileName  = sys.argv[1]
    outFileName = sys.argv[2]

    MIN_LOOP_LENGTH = 4
    INVERT_PP       = False
    extendRNAStemsInFile(inFileName, outFileName, MIN_LOOP_LENGTH, INVERT_PP)

    print 'finished!'


#run main
main()
