# made by Jakob skou Pedersen - Nov 1 2003

import col
import copy
import time


def namesAndSequences_old(colSeq):
    """Return tuple consisting of list of names and list of sequences"""
    symbols   = colSeq.column("symbols")
    N         = len(symbols[0])
    L         = len(symbols)
    sequences = []
    names     = []

    for i in range(N):
        sequences.append([])
        taxon = "taxon" + str(i+1)
        names.append(colSeq.feature(taxon))

    for i in range(N):
        for j in range(L):
            sequences[i].append(symbols[j][i])

    return (names, sequences)



def namesAndSequences(colSeq, columnIndexes):
    """Return tuple consisting of list of names and list of sequences"""

    names     = []
    sequences = []
    for index in columnIndexes:
        # get name
        tag   = colSeq.columnTag(index)
        begin = tag.find('"') + 1
        end   = tag.find('"', begin)
        name  = tag[begin:end]
        names.append(name)
        # get sequence
        sequences.append( colSeq.columnByIndex(index) )

    return (names, sequences)


def writePhylipFile(fileName, names, sequences):
    """ helper function for outputting sequences in phylip format """
   
    N         = len(names)
    L         = len(sequences[0])

    # write phy file
    file = open(fileName, 'w')
    file.write(str(N) + " " + str(L) + "\n")
    for i in range(N):
        file.write(names[i] + "  \t")
        for j in range(L):
            file.write(sequences[i][j])
        file.write("\n")


def colSeq2phy_old(colSeq, fileName):
    """Convert colSeq to phylip sequential format and save to 'fileName' """

    symbols = colSeq.column("symbols")
    (names, sequences) = namesAndSequences_old(colSeq)
    # write phy file
    writePhylipFile(fileName, names, sequences)

    #     # write phy file
    #     file = open(fileName, 'w')
    #     file.write(str(N) + " " + str(L) + "\n")
    #     for i in range(N):
    #         file.write(names[i] + "  \t")
    #         for j in range(L):
    #             file.write(sequences[i][j])
    #         file.write("\n")


def colSeq2phy(colSeq, fileName):
    """ Convert colSeq to phylip sequential format and save to 'fileName' """

    columnIndexes = symbolsColumnIndexes(colSeq)
    (names, sequences) = namesAndSequences(colSeq, columnIndexes)

    # check
    if not columnIndexes:
        raise Exception, "no 'symbols' column in colSeq"

    # write phy file
    writePhylipFile(fileName, names, sequences)


def colSeq2strict_phy(colSeq, fileName):
    """ Convert colSeq to phylip sequential format and save to 'fileName'. Sequence i is  renamed to seq_i. """

    columnIndexes = symbolsColumnIndexes(colSeq)
    (names, sequences) = namesAndSequences(colSeq, columnIndexes)

    # check
    if not columnIndexes:
        raise Exception, "no 'symbols' column in colSeq"

    for i in range(len(names)):
        name     = names[i][:10]#'seq_' + str(i)
        names[i] = name.ljust(10)

   
    N         = len(names)
    L         = len(sequences[0])

    # write phy file
    file = open(fileName, 'w')
    file.write(str(N) + " " + str(L) + "\n")
    for i in range(N):
        file.write(names[i])
        for j in range(L):
            file.write(sequences[i][j])
        file.write("\n")


def symbolsColumnIndexes(colSeq):
    """ Return list of all columns in colSeq which has 'symbols' in their tag """

    symbolsColumnIndexes= []
    for i in range(colSeq.columnCount() ):
        if colSeq.columnTag(i).find('symbols') != -1:
            symbolsColumnIndexes.append(i)

    return symbolsColumnIndexes


def makeSymbolsList_old(sequences):
    N = len(sequences)
    L = len(sequences[0])

    symbols = []
    for i in range(L):
        symbol = ''
        for j in range(N):
            symbol += sequences[j][i]
        symbols.append(symbol)

    return symbols


def namesAndSeqs2col_old(id, names, sequences):
    """ Take name and sequence list and return a col seq object """

    colSeq = col.ColumnSeq(id)
    n = 1
    for name in names:
        colSeq.setFeature('taxon' + str(n), name)
        n += 1
    
    colSeq.setColumn('symbols', makeSymbolsList_old(sequences) )

    return colSeq


def readAln(fileName):
    sequences = []
    names     = []

    file = open(fileName, 'r')
    for line in file:
        name, sequence = line.split(':', 1)
        names.append(name.strip())
        sequences.append(sequence.strip())


    return (names, sequences)


def aln2col(fileName):

    baseName, suffix = fileName.strip('./').split('.')
    names, sequences = readAln(fileName)

    return namesAndSeqs2col_old(baseName, names, sequences)


def strList2FloatList(strList):
    return [float(x) for x in strList]


def phy2col(inFileName):
    """Convert phylip sequential format to colSeq  """

    file  = open(inFileName, 'r')
    lines = file.readlines()

    #parse header
    header = lines[0].split() 
    N      = int(header[0])
    L      = int(header[1])

    # parse body
    names     = []
    sequences = []

    for line in lines[1:]:

        if len(line) != 0:
            line_list = line.split()
            names.append(line_list[0])
            list = []
            for char in line_list[1]:
                list += char
            sequences.append(list)

    #check
    if len(names) != N:
        raise col.Error, "Phylip file '" + inFileName + "' claims to have " + str(N) + " names but has " + str(len(names)) + " ." 

    # create ColumnSeq object
    colSeq = col.ColumnSeq()

    for i in range(len(names)):
        colSeq.setColumn('symbols taxa="' + names[i] + '"', sequences[i])

    entryName = ""
    if inFileName[-4:] == '.phy':
        entryName = inFileName[:-4]
    else:
        entryName = inFileName
    colSeq.setFeature('ENTRY', entryName)

    return colSeq


def pip2col(inFileName):
    """Convert pip format to colSeq  """

    file  = open(inFileName, 'r')
    lines = file.readlines()

    #parse header
    header = lines[0].split() 
    N      = int(header[0])
    L      = int(header[1])

    #debug 
    print str(N)
    print str(L)

    # parse body
    names     = []
    sequences = []

    for line in lines[1:1+N]:
        names.append( line.strip() )

    for line in lines[1+N:1+2*N]:
        list = []
        for char in line:
            list += char
        sequences.append(list)

    #check
    if len(names) != N:
        raise col.Error, "Phlip file '" + inFileName + "' claims to have " + str(N) + " names but has " + str(len(names)) + " ." 

    #check
    if len(sequences) != N:
        raise col.Error, "Phlip file '" + inFileName + "' claims to have " + str(N) + " sequences but has " + str(len(sequences)) + " ." 

    # create ColumnSeq object
    colSeq = col.ColumnSeq()

    for i in range(len(names)):
        colSeq.setColumn('symbols taxa="' + names[i] + '"', sequences[i])

    entryName = ""
    if inFileName[-4:] == '.txt':
        entryName = inFileName[:-4]
    else:
        entryName = inFileName
    colSeq.setFeature('ENTRY', entryName)

    return colSeq


def initColumnSeqWithFeatures(colSeq):
    """ return new ColumnSeqList which have all features of colSeqList """
    # initialize new ColumnSeq
    newColSeq = col.ColumnSeq( colSeq.feature('ENTRY') )
    for i in range(colSeq.featureCount() ):
        tag = colSeq.featureTag(i)
        newColSeq.setFeature(tag, colSeq.feature(tag) )

    return newColSeq


def reverseComplement(entry):
    """ Reverse complement all symbol collumns  """

    metaNuc       = ['a','c','g','t','A','C','G','T','m','M','r','R','w','W','s','S','y','Y','k','K','b','B','d','D','h','H','v','V','n','N','x','X','?','~','o','O','-']
    compl_metaNuc = ['t','g','c','a','T','G','C','A','k','K','y','Y','w','W','s','S','r','R','m','M','v','V','h','H','d','D','b','B','n','N','x','X','?','~','o','O','-']

    for i in range( entry.columnCount() ):
        if entry.columnTag(i).find('symbols') != -1:
            column = []
            for symbol in entry.columnByIndex(i):

                #complementing
                index = metaNuc.index(symbol)
                if (index == -1):
                    raise col.Error, "Unknown symbol '" + symbol + "'."

                column += compl_metaNuc[index]
                
            column.reverse()
            entry.setColumn(entry.columnTag(i), column)


def slice(entry, begin, end):
    """ Extracts the subsection [begin;end) """

    #check
    if begin < 0 or end > entry.size():
        raise col.Error, "Bounds out of range."

    colSeq = initColumnSeqWithFeatures(entry)
    for i in range( entry.columnCount() ):
        colSeq.setColumn( entry.columnTag(i), entry.columnByIndex(i)[begin:end] )
        
    return colSeq


def split(entry, size, offSet = 1):
    """ Split entry into sub sections of size 'size'. Return
    ColumnSeqList with the new entries names 'original_entry_name' +
    '_<from>_<to>'"""
            
    entry_name = entry.feature('ENTRY')
    splitList  = col.ColumnSeqList()
    
    end = 0
    begin = 0
    L = entry.size()
    while end < L:
        end += size
        if end > L:
            end = L
        subEntry     = slice(entry, begin, end)
        subEntryName = entry_name + '_'  + str(begin + offSet) + '_' + str(end + offSet)
        subEntry.setFeature('ENTRY', subEntryName)
        splitList.setEntry(subEntry)
        begin += size

    return splitList


def listSplit(entry, sectionList, offSet = 1):

    entry_name = entry.feature('ENTRY')
    splitList  = col.ColumnSeqList()

    L = entry.size()

    for touple in sectionList:
        begin, end = touple

        # check
        if begin < 0 or begin > end or end > L:
            raise Exception, "begin '" + str(begin) + "' or end '" + end + "' out of range "
        
        subEntry     = slice(entry, begin, end)
        subEntryName = entry_name + '_'  + str(begin + offSet) + '_' + str(end + offSet)
        subEntry.setFeature('ENTRY', subEntryName)
        splitList.setEntry(subEntry)

    return splitList


def grep(entry, colName, value):
    """ Extracts all rows i for which entry.column(colName)[i] equals value """

    colSeq = initColumnSeqWithFeatures(entry)
    column = entry.column(colName)

    grepMask = []
    for symbol in column:
        if symbol == value:
            grepMask.append(True)
        else:
            grepMask.append(False)

    for i in range( entry.columnCount() ):
        colSeq.setColumn( entry.columnTag(i), maskFilterList(entry.columnByIndex(i), grepMask) )
        
    return colSeq



def cat(colSeqList):

    """ Return a copy of the first entry in list with columns being
    concatenations of columns from all entries in list. All entries
    must have identical columns """

    firstEntry = colSeqList.entryByIndex(0)
    colTagList = []
    colList = {}
    for i in range( firstEntry.columnCount() ):
        tag = firstEntry.columnTag(i)
        colTagList.append(tag)
        colList[tag] = []

    # check and column cat
    for entry in colSeqList:
        for tag in colTagList:
            if not entry.hasColumn(tag):
                raise Exception, "colSeq entry '" + entry.feature('ENTRY') + "' miss column with tag '" + tag + "'."
            colList[tag] += entry.column(tag)

    # build return object
    catColSeq = initColumnSeqWithFeatures(colSeqList.entryByIndex(0) )
    for tag in colTagList:
        catColSeq.setColumn(tag, colList[tag])

    return catColSeq


def maskFilterList(list, boolMask):
    """ Create copy of list including all entries for which boolmask equals true """

    #check
    if (len(list) != len(boolMask) ):
        raise Exception, "list and boolMask differ in size"

    newList = []
    for i in range(len(list)):
        if boolMask[i]:
            newList.append(list[i])

    return newList


def makeBoolMask(colSeq, colName, symbol):
    """ Return a boolMask (a list of true and false values), with
    value true at positions where column with colName equals symbol
    """

    boolMask = []
    for s in colSeq.column(colName):
        if s == symbol:
            boolMask.append(True)
        else:
            boolMask.append(False)

    return boolMask


def makeBoolMaskByValue(colSeq, colName, minValue):
    """ Return a boolMask (a list of true and false values), with
    value true at positions where value of column with colName is at
    least minValue. Ignore '.' values (given true value) """

    boolMask = []
    for s in colSeq.column(colName):
        if s == '.':
            boolMask.append(True)
        elif float(s) >= minValue:
            boolMask.append(True)
        else:
            boolMask.append(False)

    return boolMask


def addPosition(entry, begin=0):
    """ Add a column containing seq position entry. Default to start counting at 0 """

    #check
    if begin < 0:
        raise col.Error, "Bounds out of range. begin='" + str(begin) +"'."

    column = []
    for i in range(entry.size()):
        column.append( str(begin + i) )
            
    entry.setColumn('position', column)



def removeGapRows(entry, colTag):
    """ Add a column containing seq position entry. Default to start counting at 0 """

    cols   = []
    noGap  = []

    found = None
    # find col index
    for i in range( entry.columnCount() ):
        if entry.columnTag(i) == colTag:
            noGap = entry.columnByIndex(i)
            cols.append(entry.columnByIndex(i) )
            found = 1
        else:
            cols.append(entry.columnByIndex(i) )

    if not found:
        Exception, "colTag '" + colTag + "' not found in entry '" + entry.feature('ENTRY')

    # create gap mask
    gapMask = copy.deepcopy(noGap)
    L = 0
    for i in range(len(noGap)):
        if noGap[i] == '-' or noGap[i] == '~':
            gapMask[i] = False
        else:
            gapMask[i] = True
            L += 1

    # modify columns
    for k in range( entry.columnCount() ):
        col    = entry.columnByIndex(k)
        newCol = L*['*']
        j = 0
        for i in range(len(gapMask)):
            if gapMask[i]:
                newCol[j] = col[i]
                j += 1
        entry.setColumn(entry.columnTag(k), newCol)


def RNAssPerformanceMeasures(colSeq):
    """ type s require only some pairing, type p require paring to be correct """
    PRED_PAIRED_POSITIONS  = 'pairedPosition'
    TRUE_PAIRED_POSITIONS  = 'pairedPositions'
    ORIGINAL_POSITIONS     = 'positions'

    truePairList    = shiftList( colSeq.column(TRUE_PAIRED_POSITIONS), int( colSeq.column(ORIGINAL_POSITIONS)[0] ) - 1)
    predPairList    = colSeq.column(PRED_PAIRED_POSITIONS)

    TN_s = 0.0
    FN_s = 0.0
    TP_s = 0.0
    FP_s = 0.0

    TN_p = 0.0
    FN_p = 0.0
    TP_p = 0.0
    FP_p = 0.0

    for i in range(len(predPairList) ):
        if truePairList[i] == '.' and predPairList[i] == '.' :
            TN_s += 1
            TN_p += 1

        elif truePairList[i] != '.' and predPairList[i] == '.' :
            FN_s += 1
            if i > int(truePairList[i]):
                FN_p += 1

        elif truePairList[i] != '.' and predPairList[i] != '.' :
            TP_s += 1
            if i > int(truePairList[i]):
                if int(truePairList[i]) == int(predPairList[i]):
                    TP_p += 1

        elif truePairList[i] == '.' and predPairList[i] != '.' :
            FP_s += 1
            if i > int(predPairList[i]):
                FP_p += 1
        else:
            print "this can not happen!"
            

    sn_s = 0
    if TP_s != 0:
        sn_s = TP_s/(TP_s+FN_s)

    sp_s = 0
    if TP_s != 0:
        sp_s = TP_s/(TP_s+FP_s)



    sn_p = 0
    if TP_p != 0:
        sn_p = TP_p/(TP_p+FN_p)

    sp_p = 0
    if TP_p != 0:
        sp_p = TP_p/(TP_p+FP_p)

    print "TN_s FN_s TP_s FP_s: ", TN_s, FN_s, TP_s, FP_s
    print "TN_p FN_p TP_p FP_p: ", TN_p, FN_p, TP_p, FP_p
    return (sn_s, sp_s, sn_p, sp_p)  # (sn, sp)


def RNAPair(nuc_r, nuc_l):
    RNA_PAIRS = ['AT', 'TA', 'GC', 'CG', 'GT', 'TG']

    pair = nuc_r + nuc_l
    pair = pair.upper().replace('U','T')

    return pair in RNA_PAIRS


def shiftList(list, shift):
    """ Return copy of list with shift subtracted values except '.' """ 
    shiftList = []

    for value in list:
        if value == '.':
            shiftList.append(value)
        else:
            shiftList.append( str( int(value) - shift) )

    return shiftList


def addPairingFractionColumn(colSeq):
    """ Add column measuring fraction of symbols which pair correctly
    according to a 'pairedPositions' and a position column. """

    PAIRING_FRACTION = 'pairingFraction'
    PAIRED_POSITIONS = 'pairedPositions'
    POSITIONS        = 'positions'

    pairCol          = shiftList( colSeq.column(PAIRED_POSITIONS), int( colSeq.column(POSITIONS)[0] ) )
    pairCount        = colSeq.size() * [0.0]
    L                = len(pairCol)

    for j in range(L):
        symbolCount = 0
        for i in range( colSeq.columnCount() ):
            if colSeq.columnTag(i).find('symbol') != -1:
                col = colSeq.columnByIndex(i)
                if pairCol[j] != '.':
                    pairedPos = int(pairCol[j])
                    if col[j] != '~' and col[pairedPos] != '~':
                        symbolCount += 1
                        if RNAPair(col[j], col[pairedPos]):
                            pairCount[j] += 1
        # normalize
        if (symbolCount):
            pairCount[j] = pairCount[j]/symbolCount
        else:
            pairCount[j] = 0.0

    pairFractionList = []
    for i in range( colSeq.size() ):
        if pairCol[i] != '.':
            pairFractionList.append( str( round(pairCount[i], 5) ) )
        else:
            pairFractionList.append( str('.') )
            
    colSeq.setColumn(PAIRING_FRACTION, pairFractionList)


       
