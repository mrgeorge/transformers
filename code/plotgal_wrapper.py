#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import numpy as np
import plotgalSMvR
import re # to remove decimal from string

imDir="../images/"
imListFile=imDir+"acs_mem.list"
plotDir="../plots/"

zStep=0.1
zBins=np.arange(0.,1.,zStep)

for zz in zBins:
    zMin=zz
    zMax=zMin+zStep
    plotFile=plotDir+"candy_z"+re.sub("\D","",str(zMin))+"_"+re.sub("\D","",str(zMax))+".pdf"
    print plotFile
    plotgalSMvR.main(imDir, imListFile, plotFile, zMin, zMax)

