import string
import Orange
from AZutilities import getCinfonyDesc
from AZutilities import dataUtilities

data = dataUtilities.DataTable("Assay1350AZOdesc.txt")

# Read nonIID cmpds
fid = open("nonIIDlist.txt")
IDlist = []
for line in fid:
    ID = string.strip(line)
    IDlist.append(ID)
fid.close()
print "Number of non IID cmpds; ", len(IDlist)

# Select nonIID cmpds from full data set
nonIID = []
IID = []
for ex in data:
    if ex["ID"].value in IDlist:
        nonIID.append(ex)
    else:
       IID.append(ex)
print "# IID ", len(IID)
print "# nonIID ", len(nonIID)

nonIIDdata = dataUtilities.DataTable(nonIID)
IIDdata = dataUtilities.DataTable(IID)

# Partition the nonIID cmpds into a test and an external test set
indices2 = Orange.data.sample.SubsetIndices2(p0=0.5)
ind = indices2(nonIIDdata)
nonIIDtest = nonIIDdata.select(ind, 0)
nonIIDextTest = nonIIDdata.select(ind, 1)

nonIIDtest.save("nonIIDtestAZOdesc.txt")
nonIIDextTest.save("nonIIDextTestAZOdesc.txt")

# Select about 100 cmpds from the IID set as an IID external test set
indices2 = Orange.data.sample.SubsetIndices2(p0=0.07)
ind = indices2(IIDdata)
IIDset = IIDdata.select(ind, 1)
IIDextTest = IIDdata.select(ind, 0)

IIDset.save("IIDsetAZOdesc.txt")
IIDextTest.save("IIDextTestAZOdesc.txt")






