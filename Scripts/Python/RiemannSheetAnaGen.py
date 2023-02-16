import sys
from matplotlib.colors import LinearSegmentedColormap
from plotstyle import SetupPlotStyle
SetupPlotStyle()
from scipy.interpolate import griddata
import sys
import RiemannSheetAna_py as rs

print ("Number of arguments:", len(sys.argv), "arguments")
print ("sys.argv[0]:", sys.argv[0])
print ("sys.argv[1]:", sys.argv[1])
print ("Argument List:", str(sys.argv))

theRiemannSheetAna = []

theRiemannSheetAna.append(rs.RiemannSheetAna_py(list(sys.argv)))

sheetSize = theRiemannSheetAna[0].sheetSize()

print ("sheet size:", sheetSize)

def binary(d):
    numbers = []
    i = 1
    for x in range(d):
        numbers.append(i)
        i *= 2
        if i > d:
            i /= 2
            i = int(i)
            break
    numbers.reverse()
    theList = []
    theListstr = ""
    i = 0
    while len(theList) < len(numbers):
        if numbers[i] <= d:
            theList.append("+")
            d -= numbers[i]
        else:
            theList.append("-")
        i += 1
    for x in theList:
        theListstr += str(x)
    return theListstr

sheetList = []

for i in range(0, sheetSize*sheetSize):
    currentSheet=binary(i)
    while len(currentSheet) < sheetSize:
        currentSheet = "-"+ currentSheet
    sheetList.append(currentSheet)
    
sheetList.reverse()
for sheetIdx in range(0,sheetSize*sheetSize):
    print("sheetList[",sheetIdx,"]= ", sheetList[sheetIdx])

for i in range(0, sheetSize*sheetSize):
    if i > 0:
        theRiemannSheetAna.append(rs.RiemannSheetAna_py(list(sys.argv)))
    theRiemannSheetAna[i].SetSheet(sheetList[i])
