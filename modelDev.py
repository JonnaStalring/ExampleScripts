from AZutilities import evalUtilities
from AZutilities import dataUtilities
from trainingMethods import AZorngRF
from trainingMethods import AZorngCvBoost
from trainingMethods import AZorngCvSVM
import orngTest
import orngStat
import orange
import string

data = dataUtilities.DataTable("IIDsetAZOdesc.txt")
test = dataUtilities.DataTable("nonIIDtestAZOdesc.txt")

descList = ["ID", "Smiles", "Conc", "Effect", "Conc_1",  "Effect_1", "ID_1", "origSmiles_1", "BioActivity_1"]
data = dataUtilities.attributeDeselectionData(data, descList)

# Deselect descriptors with no variance
descList = ["rdk.fr_dihydropyridine", "rdk.fr_nitroso", "rdk.fr_benzodiazepine", "rdk.fr_thiocyan", "rdk.VSA_EState4" ,"rdk.VSA_EState6" \
,"rdk.VSA_EState7" ,"rdk.VSA_EState1" ,"rdk.VSA_EState2" ,"rdk.VSA_EState3" ,"rdk.SlogP_VSA9" ,"rdk.SMR_VSA8" ,"rdk.fr_diazo" \
,"rdk.fr_prisulfonamd" ,"rdk.fr_isocyan" ,"rdk.fr_azide" ,"rdk.fr_isothiocyan"]
data = dataUtilities.attributeDeselectionData(data, descList)
print "Length domain ", len(data.domain)

learner = AZorngCvSVM.CvSVMLearner(C=32, gamma=0.03125)
#learner = AZorngRF.RFLearner()
#learner = AZorngRF.RFLearner(stratify = "Yes") # No effect
#learner = AZorngCvBoost.CvBoostLearner()
#learner.stratify = "Yes" # No effect
#learner.priors = {"Active":0.80, "Inactive":0.20}


# Test set accuracy
model = learner(data)
res = orngTest.testOnData([model], data)
CM = evalUtilities.ConfMat(res)[0]
CA = round(orngStat.CA(res)[0],3)
MCC = round(evalUtilities.calcMCC(CM),3)
# TH, FL, FH, TL
resList = [str(CM[0][0]), str(CM[0][1]), str(CM[1][0]), str(CM[1][1]), str(CA), str(MCC)]
wrtStr = string.join(resList, "\t")
print "nonIID test set results"
print wrtStr


# CV accuracy
res = orngTest.crossValidation([learner], data, strat=orange.MakeRandomIndices.StratifiedIfPossible, folds = 10)
CM = evalUtilities.ConfMat(res)[0]
CA = round(orngStat.CA(res)[0],3)
MCC = round(evalUtilities.calcMCC(CM),3)
# TH, FL, FH, TL
resList = [str(CM[0][0]), str(CM[0][1]), str(CM[1][0]), str(CM[1][1]), str(CA), str(MCC)]
wrtStr = string.join(resList, "\t")
print "CV results"
print wrtStr

# Test set accuracy
model = learner(data)
res = orngTest.testOnData([model], test)
CM = evalUtilities.ConfMat(res)[0]
CA = round(orngStat.CA(res)[0],3)
MCC = round(evalUtilities.calcMCC(CM),3)
# TH, FL, FH, TL
resList = [str(CM[0][0]), str(CM[0][1]), str(CM[1][0]), str(CM[1][1]), str(CA), str(MCC)]
wrtStr = string.join(resList, "\t")
print "nonIID test set results"
print wrtStr




