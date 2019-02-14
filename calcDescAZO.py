import string
from AZutilities import getCinfonyDesc
from AZutilities import dataUtilities


# Parameters
assay = "Assay1350"


def parseToFile(assay, data):
    descList = getCinfonyDesc.getAvailableDescs("rdkPhysChem")
    print descList    
    data = getCinfonyDesc.getCinfonyDescResults(data, descList)
    data.save(assay+"AZOdesc.txt")


data = dataUtilities.DataTable(assay+".txt")
print "Lenght of molDict ", str(len(data))

# Calculate descriptors and write to file
parseToFile(assay, data)
