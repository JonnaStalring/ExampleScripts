import string
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem
import numpy
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina
from rdkit.Chem import rdMolDescriptors
#from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw

descList = []

def getFps(smilesList):

    molList = []
    for smiles in smilesList:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            molList.append(mol)
        else:
            print smiles
    fps = [FingerprintMols.FingerprintMol(x) for x in molList] # Topological
    #fps = [AllChem.GetMorganFingerprint(x, 2) for x in molList]
    #print "Length of data and fp ", len(data), len(fps)
    return fps, molList


def GenerateLowerTriangularDistanceMatrix(MolsFingerprints):
    """Generate a lower triangular distance matrix without the diagonal."""

    DistanceMatrix = []
    NumFPs = len(MolsFingerprints)
    for Index1 in range(0, NumFPs):
        for Index2 in range(0, Index1):
            Distance =  1 - DataStructs.FingerprintSimilarity(MolsFingerprints[Index1], MolsFingerprints[Index2],)
            DistanceMatrix.append(Distance)

    return DistanceMatrix


fid = open("Assay1350AZOdesc.txt")
smilesList = []
IDList = []
responseList = []
dataList = []
header = fid.readline()
for line in fid:
    lineList = string.split(line, "\t")
    smilesList.append(lineList[2])
    IDList.append(lineList[3])
    dataList.append(line)
    responseList.append(string.strip(lineList[186]))
fid.close()

# Calculate fingerprints for train and test sets
fps, molList = getFps(smilesList)

print "Clustering ...."
DistanceMatrix = GenerateLowerTriangularDistanceMatrix(fps)
ClusteredMolIndices = Butina.ClusterData(DistanceMatrix, len(fps), 0.25, isDistData = True)
print "End of clustering"
print "Number of clusters ", len(ClusteredMolIndices)

# Get fp clusters
idx = 0
for Cluster in ClusteredMolIndices:
    idx = idx + 1
    if len(Cluster) > 10:
        print "***************** In cluster ", idx
        print len(Cluster)
        molCluster = []
        IDCluster = []
        responseCluster = []
        for molIdx in Cluster:
            molCluster.append(molList[molIdx])
            IDCluster.append(IDList[molIdx])
            responseCluster.append(responseList[molIdx])
        #MCS = rdFMCS.FindMCS(molCluster)
        #print MCS.smartsString
        LC = HC = 0
        for elem in responseCluster:
            if elem == "Inactive":
                LC = LC + 1
            else:
                HC = HC + 1
        print "Inactive  within cluster ", LC
        print "Active  within cluster ", HC

        fid = open("Clusters/clusterData_"+str(idx)+".txt", "w")
        fid.write(header)
        for ID in IDCluster:
            cmpdIdx = IDList.index(ID)
            fid.write(dataList[cmpdIdx])
        fid.close()
        
        #MCSmol = Chem.MolFromSmarts(MCS.smartsString)
        #Draw.MolToFile(MCSmol, "ClusterMCS/MCS"+str(len(Cluster))+"_"+str(LC)+"_"+str(HC)+".png", size=(300, 300))


