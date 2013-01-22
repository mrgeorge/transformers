#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import numpy as np
import plotgalSMvR
import re # to remove decimal from string

imDir="../images/twenty_arcsec/"
imListFile=imDir+"acs_mem.list"
plotDir="../plots/"

zStep=0.1
zBins=np.arange(0.1,1.,zStep)
minMh=12.0
maxMh=15.0

#morphTypes=["all","spheroidal","bulge+disk","latedisk","xdisk","irr"]
morphTypes=["all"]
colorType="all"
morphLegend=False
zType="zp"

for zz in zBins:
    for morph in morphTypes:
        minZ=zz
        maxZ=minZ+zStep
        plotFile=plotDir+"candy_z"+re.sub("\D","",str(minZ))+"_"+re.sub("\D","",str(maxZ))+"_m"+re.sub("\D","",str(minMh))+"_"+re.sub("\D","",str(maxMh))+"_"+morph+".pdf"
        print plotFile
        plotgalSMvR.main(imDir, imListFile, plotFile, minZ, maxZ, zStep, minMh, maxMh, colorType, morph, morphLegend, zType)
