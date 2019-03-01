import Orange
from AZutilities import evalUtilities
from AZutilities import dataUtilities
from AZutilities import AZOrangePredictor
from trainingMethods import AZorngRF
import orngTest
import orngStat
import orange
import string
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem
import numpy
from rdkit.Chem import AllChem
import commands



def getDist(queryFp, fps):

    distList = []
    for fp in fps:
        dist = DataStructs.FingerprintSimilarity(queryFp,fp) # Tanimoto
        #dist = DataStructs.DiceSimilarity(queryFp,fp) # Dice
        distList.append(dist)
    distList.sort()
    #medDist = numpy.median(distList[len(distList)-9:len(distList)])
    medDist = numpy.median(distList[len(distList)-4:len(distList)])
    NN1Dist = distList[len(distList)-1:len(distList)][0]
    return medDist, NN1Dist


def getFps(data):

    molList = []
    for ex in data:
        mol = Chem.MolFromSmiles(ex["Smiles"].value)
        if mol:
            molList.append(mol)
        else:
            print ex["Smiles"].value
            print ex["Leonumber"].value
    fps = [FingerprintMols.FingerprintMol(x) for x in molList] # Topological
    #fps = [AllChem.GetMorganFingerprint(x, 2) for x in molList]
    #print "Length of data and fp ", len(data), len(fps)
    return fps

THRS = 0.75

model = AZorngRF.RFread("OI_RFmodel")
predictor = AZOrangePredictor.AZOrangePredictor("OI_RFmodel")

train = dataUtilities.DataTable("BioActivityAZOdesc.txt")


# Calculate fingerprints for train and test sets
fps = getFps(train)


#smiles = test[idx]["Smiles"].value
smiles = "CC(C)n1c(/C=C/[C@H](O)C[C@H](O)CC(=O)O)c(-c2ccc(F)cc2)c2ccccc21"
smiles = "Cc1cc(=Nc2cc(CN3CCOCC3)c3nc(C)c(Cc4ccc(Cl)cc4F)n3n2)[nH][nH]1"  # Train set
#smiles = "Cc1nc2c(CN3CCOCC3)cc(NC3=CC(C)NN3)nn2c1Cc1ccc(Cl)cc1F"  # From Drawing - Wrong no tautomer
smiles = "Cc1cc(Nc2cc(CN3CCOCC3)c3nc(C)c(Cc4ccc(Cl)cc4F)n3n2)[nH]n1" #From drawing of Galilei structure
#smiles = "Cc1cc(=Nc2cc(CN3CCOCC3)c3nc(C)c(Cc4ccc(Cl)cc4F)n3n2)[nH][nH]1" # Canonicalized from drawing in Galilei
cmd = "env -i HOME='$HOME' bash -l -c './cleanSmiles.sh "+'"'+smiles+'"'+"'"
print cmd
status, cleanSmiles = commands.getstatusoutput(cmd)
print cleanSmiles
predictor.getDescriptors(cleanSmiles)
prediction, prob = predictor.predict(True)
# Normalize prob
prob = 200*abs(prob)
if prob > 100:
    prob = 0.98
#print "****************** Prediction with AZOrangePredictor ******************"
print prediction
print prob

# Create an Orange data set to calculate the fp of the smiles
features = [Orange.data.variable.String("Smiles")]
domain = Orange.data.Domain(features)
smiData = Orange.data.Table(domain)
smiData.append([cleanSmiles])
fpsSmiles = getFps(smiData) 
distSmi, distNN1 = getDist(fpsSmiles[0], fps)
print "Dist from smiles ", distSmi
print "Dist from smiles ", distNN1


























