import string

def appendCluster(fileName, IIDlist):
    fid = open(fileName)
    for line in fid:
        lineList = string.split(line, "\t")
        ID = lineList[3]
        if ID not in IIDlist:
            IIDlist.append(ID)
        else:
            print "Duplicate in cluster"
    fid.close()
    return IIDlist
    
IIDlist = []
#idxList = ["8", "19", "20", "28", "31", "32"]
idxList = ["2", "1"]
for idx in idxList:
    fileName = "Clusters/clusterData_"+idx+".txt"
    IIDlist = appendCluster(fileName, IIDlist)

fid = open("nonIIDlist.txt", "w")
for elem in IIDlist:
    fid.write(elem+"\n")
fid.close()
