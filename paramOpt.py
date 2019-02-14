import AZLearnersParamsConfig
from AZutilities import dataUtilities
from AZutilities import paramOptUtilities
from trainingMethods import AZorngCvSVM

data = dataUtilities.DataTable("IIDsetAZOdesc.txt")

trainFile = "/ChemistryData/jgw/projects/OI/Assay1350/trainFile.txt"

descList = ["ID", "Smiles", "Conc", "Effect", "Conc_1",  "Effect_1", "ID_1", "origSmiles_1", "BioActivity_1"]
data = dataUtilities.attributeDeselectionData(data, descList)

# Deselect descriptors with no variance
descList = ["rdk.fr_dihydropyridine", "rdk.fr_nitroso", "rdk.fr_benzodiazepine", "rdk.fr_thiocyan", "rdk.VSA_EState4" ,"rdk.VSA_EState6" \
,"rdk.VSA_EState7" ,"rdk.VSA_EState1" ,"rdk.VSA_EState2" ,"rdk.VSA_EState3" ,"rdk.SlogP_VSA9" ,"rdk.SMR_VSA8" ,"rdk.fr_diazo" \
,"rdk.fr_prisulfonamd" ,"rdk.fr_isocyan" ,"rdk.fr_azide" ,"rdk.fr_isothiocyan"]
data = dataUtilities.attributeDeselectionData(data, descList)
print "Length domain ", len(data.domain)

print "Running param opt"
data.save(trainFile)
learner = AZorngCvSVM.CvSVMLearner()
evalM = "AZutilities.evalUtilities.CA"
fMin = False
optimizer = paramOptUtilities.Appspack()
learnerName = "CvSVMLearner"
runPath = "/ChemistryData/jgw/projects/OI/Assay1350/ParamOpt"
pars = AZLearnersParamsConfig.API(learnerName)
tunedPars = optimizer(learner=learner,\
                    dataSet=trainFile,\
                    evaluateMethod = evalM,\
                    useParameters = pars.getParametersDict(),\
                    findMin=fMin,\
                    useStd = False,\
                    runPath = runPath,\
                    verbose = 0)
verbTunedPars = optimizer.getTunedParameters()
print "Returned: ", tunedPars
print "====================== optimization Done ==========================="
print "Learner optimized flag = ", learner.optimized
print "Tuned parameters = ", tunedPars[1]
print "Best optimization result = ", tunedPars[0]
print "Best result index from intRes file:",verbTunedPars["ResIdx"]

