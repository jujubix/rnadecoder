# made by Jakob skou Pedersen - Nov 1 2003

import copy

class ColumnSeq:
    """Class containing biological sequence information. Modelled over\
    the col format: 'http://www.bioinf.au.dk/colformat/index.html' """

    # member functions

    def __init__(self, entryName = 'unspecified'):
        # data members
        self.featureTags_ = []   #list : enumerate feature tags (remember order)
        self.features_    = {}   #list : map featureTags to features
        self.columnTags_  = []   #list : enumerate columnTags
        self.columns_     = []   #list : same enumeration as columnTags

        self.setFeature('ENTRY', entryName)


    def setColumn(self, tag, column):
        """Set or reset column"""

        if (self.columnTags_.count(tag) == 0):
            self.columnTags_.append(tag)
            self.columns_.append(column)
        else:
            self.columns_[ self.columnTags_.index(tag) ] = column


    def setFeature(self, tag, feature):
        if (self.featureTags_.count(tag) == 0):
            self.featureTags_.append(tag)
            self.features_[tag] = feature
        else:
            self.features_[tag] = feature

        
    def column(self, tag):
        return self.columns_[self.columnTags_.index(tag)]


    def columnByIndex(self, i):
        return self.columns_[i]


    def columnTag(self, i):
        return self.columnTags_[i]


    def hasColumn(self, tag):
        return self.columnTags_.count(tag)


    def delColumn(self, tag):
        index = self.columnTags_.index(tag)
        del self.columnTags_[index]
        del self.columns_[index]


    def feature(self, tag):
        return self.features_[tag]


    def featureTag(self, i):
        return self.featureTags_[i]


    def featureCount(self):
        return len(self.featureTags_)


    def hasFeature(self, tag):
        return self.features_.has_key(tag)


    def delFeature(self, tag):
        index = self.featureTags_.index(tag)
        del self.featureTags_[index]
        del self.features_[index]


    def columnCount(self):
        return len(self.columns_)


    def size(self):
        if self.columnTags_:
            return len(self.columns_[0])
        else:
            return 0


class ColumnSeqList:
    "Class containing list of ColumnSeqs"

    # member functions

    def __init__(self):
        # data members
        self.commonFeatures_  = {} #dict : map featureTags to features
        self.entryNames_      = [] #list : enumerate entry names
        self.entries_         = [] #list : store entries in same order as entryNames

    def read(self, fileName):
        """Read and parse file in col format into ColumnSeqList object"""
        self.erase()

        lines = open(fileName).readlines()

        N = 0;
        for line in lines:
            if line.find("**********") != -1:
                N += 1
        for i in range(len(lines)):
            lines[i] = lines[i].strip() # remove the '\n'

        i = 0;
        # parse header
        while (lines[i].find("==========") == -1):

            while lines[i].isspace() :
                i += 1

            if (lines[i][0] != ';') or (not lines[i].count('\t') ) :
                raise Error, "ColumnSeqList file '" + fileName + "' has corrupt format in line " + str(i+1)

            (tag, feature) = lines[i].lstrip('; \t').split('\t',1)
            self.setCommonFeature(tag, feature)

            i += 1

        
        # parse entries
        for entry in range(N):

            featureTags  = []   #list : enumerate feature tags (remember order)
            features     = {}   #dict : map featureTags to features
            columnTags   = []   #list : enumerate columnTags
            columns      = []   #list : follow enumeration in columnTags

            # parse entry header
            i += 1

            while lines[i].find("----------") == -1:

                if lines[i].isspace() or (len(lines[i]) == 0):
                    i += 1
                    
                else:
                    if lines[i][0] != ';':
                        raise Error, "ColumnSeqList file '" + fileName + "' has corrupt format in line " + str(i+1)
                    tag     = ''
                    feature = ''
                    if lines[i].count('\t'):
                        (tag, feature) = lines[i].lstrip('; \t').split('\t',1)
                    else:
                        (tag, feature) = lines[i].lstrip('; \t').split('  ',1)  # acommodate Bjarne Knudsen's seperation
                        (tag, feature) = ( tag.strip(), feature.strip() )
                
                    if tag[:3] == "COL":
                        # assume columnTag appear in correct order
                        columnTags.append(feature)
                        columns.append([])
                    else:
                        featureTags.append(tag)
                        features[tag] = feature

                    i += 1
            
            # parse columns
            i += 1
            while lines[i].find("**********") == -1:

                if lines[i].isspace() or ( len(lines[i]) == 0):
                    i += 1
                else:
                    values = lines[i].strip().split()
                    for j in range(len(values)):
                        columns[j].append(values[j])
                    i += 1

            # create ColumnSeq object
            colSeq = ColumnSeq()

            for tag in featureTags:
                colSeq.setFeature(tag, features[tag])
            for j in range(len( columnTags ) ):
                colSeq.setColumn(columnTags[j], columns[j])
            if ( not colSeq.hasFeature("ENTRY") ):
                colSeq.setFeature("ENTRY", "unspecified_" + str(entry) )
            if ( not colSeq.hasFeature("TYPE") ):
                colSeq.setFeature("TYPE", "unspecified")

            self.setEntry(colSeq)


    def write(self, fileName):
        """Write self in col format into file"""
        file = open(fileName, 'w')

        #write header
        for tag in self.commonFeatures_:
            file.write("; " + tag + "\t" + self.commonFeatures_[tag] + "\n")
        file.write("; " + "==========\n")

        #write entries
        for entry in self:
            
            for tag in entry.featureTags_:
                file.write("; " + tag + "\t" + entry.features_[tag] + "\n")
            for j in range(len( entry.columnTags_ ) ):
                file.write("; " + "COL " + str(j) + "\t" + entry.columnTags_[j] + "\n")
            file.write("; " + "----------" + "\n")

            colCount = len(entry.columns_)
            if colCount:
                for k in range(len(entry.columns_[0] ) ):
                    for col in range(colCount):
                        file.write(entry.columns_[col][k])

                        if col == colCount - 1:
                            file.write("\n")
                        else:
                            file.write("\t")
            file.write("; " + "**********\n")

        
    def setCommonFeature(self, tag, feature):
        self.commonFeatures_[tag] = feature


    def setEntry(self, columnSeq):
        name = columnSeq.feature("ENTRY")

        if (self.entryNames_.count(name) == 0):
            self.entryNames_.append(name)
            self.entries_.append(columnSeq)
        else:
            index = self.entryNames_.index(name)
            self.entries_[index] = columnSeq


    def entry(self, name):
        index = self.entryNames_.index(name)
        return self.entries_[index]


    def entryByIndex(self, i):
        return self.entries_[i]


    def size(self):
        return len(self.entries_)


    def add(self, colSeqList):
        """ add all entries in colSeqList to self. Will overwrite entries of same name """
        for entry in colSeqList:
            newEntry = copy.deepcopy(entry)
            self.setEntry(newEntry)

        return self

    def erase(self):
        # data members
        self.commonFeatures_  = {} 
        self.entryNames_      = [] 
        self.entries_         = [] 

    def __getitem__(self, i):
        return self.entries_[i] 

    def __iter__(self):
        """Iterate over entries"""
        iter = ColumnSeqListIter()
        iter.data  = self.entries_
        iter.index = -1
        return iter



class Error(Exception):
    """Exception type used in this module"""
    pass


class ColumnSeqListIter:
    def next(self):
        self.index += 1
        if self.index == len(self.data):
            raise StopIteration
        return self.data[self.index]
        

# # test
# colList = ColumnSeqList()
# colList.read('musacasa_label.col')

# print colList.entryNames_
# for entry in colList:
#     print entry.feature("ENTRY")

# colList.write('tmp.col')

# c = ColumnSeq()
