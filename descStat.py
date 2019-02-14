import Orange
from AZutilities import dataUtilities

data = dataUtilities.DataTable("Assay1350AZOdesc.txt")

bas = Orange.statistics.basic.Domain(data)

print "%20s %5s %5s %5s" % ("feature", "min", "max", "avg")
for a in bas:
    if a:
        print "%20s %5.3f %5.3f %5.3f" % (a.variable.name, a.min, a.max, a.avg)
        #if a.min == a.max:
        #    print a.variable.name

#for ex in data:
#    print ex["rdk.fr_dihydropyridine"]
