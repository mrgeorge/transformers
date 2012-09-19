#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

from plot_colormorph import *
from esutil import cosmology

cosmo=cosmology.Cosmo(h=0.72,omega_m=0.258,omega_l=0.742) # WMAP5

# TO DO:

# create histogram plots showing color bimodality in different mass, z, morph selections
# overlay dashed line showing 3.5 boundary
# label different curves, or use different panels?
# use N for y-axis and a simple bin width in color, so that relative sizes of different curves can clearly show differences in population sizes

# panel ideas (first just make as single plots):
#  3 mass bins at low z
#  3 mass bins at high z
#  3 morph types at low z (all M?)
# etc

def setupHistPlot(xtitle=r"NUV - r$^+$"):

    # use helvetica and latex
    plt.clf()
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
    plt.rc('text', usetex=True)
    plt.rc('axes', linewidth=1.5)
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    plt.xlabel(xtitle,fontsize='medium')
    plt.ylabel("N",fontsize='medium')

def histPlotLegend():
    plt.legend(prop={'size':12},loc='upper right',frameon=False)

def oplotZTitle(zRange,xPos,yPos,fontsize):
    zStr=r"z = [{}-{})".format(zRange[0],zRange[1])
    plt.gca().text(xPos,yPos,zStr,verticalalignment="top",fontsize=fontsize)

def oplotSMTitle(smRange,xPos,yPos,fontsize):
    smStr=r"log(M$_{\star}$/M$_{\odot}$)=["+str(smRange[0])+", "+str(smRange[1])+")"
    plt.gca().text(xPos,yPos,smStr,verticalalignment="top",fontsize=fontsize)

def addVolumeAxis(zRange):
    area=1.64 * (np.pi/180.)**2
    skyFrac=area/(4.*np.pi)
    vol=cosmo.V(zRange[0],zRange[1])*skyFrac # Mpc^3

    yMax=2.

    vFactor=vol/(10**4)

    ax=plt.gca()

    # copy N axis over to the right size
    ax2=ax.twinx()
    ax2.set_ylabel("N")
    ax2.set_ylim(np.array([0,yMax])*vFactor)

    # now change the left axis to show number density
    ax.set_ylabel(r"N / V$_{\rm comoving}$ [10$^{-4}$ Mpc$^{-3}$]")
    ax.set_ylim(np.array([0,yMax])*vFactor)
    tickLocs=np.linspace(0,yMax,num=5)
    ax.yaxis.set_ticks(tickLocs*vFactor)
    ax.yaxis.set_ticklabels(tickLocs)

def makeHistPlot(arrs, labels, colors, linestyles, zRange, smRange):

    cRange=np.array([-1.,6.5])
    cBin=0.5
    nCBins=(cRange[1]-cRange[0])/cBin

    lw=2

    setupHistPlot()
    plt.xlim(cRange)
    plt.hist(arrs,histtype='step',range=cRange,bins=nCBins,label=labels,color=colors,lw=lw)
    histPlotLegend()

    addVolumeAxis(zRange)

    xlim=plt.gca().get_xlim()
    ylim=plt.gca().get_ylim()
    xPos=xlim[0]+0.05*(xlim[1]-xlim[0])
    yPos=ylim[0]+0.95*(ylim[1]-ylim[0])
    fontsize=14
    oplotZTitle(zRange,xPos,yPos,fontsize)
    oplotSMTitle(smRange,xPos,yPos-0.07*(ylim[1]-ylim[0]),fontsize)


def selectHistData(data,morph,zMin=0,zMax=1,smMin=2,smMax=15,morphType="All",zType="zp"):

    if(zType == "zb"): # zbest, use speczs when available
        zphot="ZPHOT"
    else: # zp for photoz only
        zphot="PHOTOZ_NON_COMB"

    morphDict={"Spheroidal":[1],"Bulge+Disk":[2],"Late Disk":[3],"All":[0,1,2,3]}

    sel=((data[zphot] >= zMin) &
         (data[zphot] < zMax) &
         (data['KEVIN_MSTAR'] >= smMin) &
         (data['KEVIN_MSTAR'] < smMax) &
         (map(lambda x: x in morphDict[morphType],morph))
        )

    return data[sel]['MNUV_MR']

def defineHistCategories(option,data,morph,zType):

    smBins=np.array([9.8,10.3,10.7,11.5])
    zBins=np.array([0.2,0.5,0.8,1.0])
    morphTypes=np.array(["Late Disk","Bulge+Disk","Spheroidal"])

    colors=["blue","green","red"]
    #    linestyles=["solid","dashed","dotted"]
    linestyles=["-","--",":"]
   
    if(option == "mass_lowz"):
        zRange=zBins[0:2]
        smRange=np.array([smBins[0],smBins[-1]])
        arrs=[selectHistData(data,morph,zMin=zBins[0],zMax=zBins[1],smMin=smBins[ii],smMax=smBins[ii+1],morphType="All",zType=zType) for ii in range(smBins.size-1)]
        labels=["[{}-{})".format(smBins[ii],smBins[ii+1]) for ii in range(smBins.size-1)]
        colors=colors
        linestyles=linestyles
    elif(option == "mass_highz"):
        zRange=zBins[-2:]
        smRange=np.array([smBins[0],smBins[-1]])
        arrs=[selectHistData(data,morph,zMin=zBins[-2],zMax=zBins[-1],smMin=smBins[ii],smMax=smBins[ii+1],morphType="All",zType=zType) for ii in np.arange(smBins.size-2)+1]
        labels=["[{}-{})".format(smBins[ii],smBins[ii+1]) for ii in np.arange(smBins.size-2)+1]
        colors=colors[1:]
        linestyles=linestyles[1:]
    elif(option == "morph_lowz_highm"):
        zRange=zBins[0:2]
        smRange=smBins[-2:]
        arrs=[selectHistData(data,morph,zMin=zBins[0],zMax=zBins[1],smMin=smBins[-2],smMax=smBins[-1],morphType=morphTypes[ii],zType=zType) for ii in range(morphTypes.size)]
        labels=[morphTypes[ii] for ii in range(morphTypes.size)]
        colors=colors
        linestyles=linestyles
    elif(option == "morph_highz_highm"):
        zRange=zBins[-2:]
        smRange=smBins[-2:]
        arrs=[selectHistData(data,morph,zMin=zBins[-2],zMax=zBins[-1],smMin=smBins[-2],smMax=smBins[-1],morphType=morphTypes[ii],zType=zType) for ii in range(morphTypes.size)]
        labels=[morphTypes[ii] for ii in range(morphTypes.size)]
        colors=colors
        linestyles=linestyles

    return (arrs,np.array(labels),np.array(colors),np.array(linestyles),zRange,smRange)


# MAIN - if called from command line
if __name__ == '__main__':

    minz=0.2
    maxz=1.0
    minmh=13.0
    maxmh=14.0
    ztype="zp" # use zb for zbest (i.e. specz if available), else zp for photoz only

    plotDir="../plots/"

    # read group and galaxy catalogs
    acsFile="../../code/lensing18_20110914.fits"
    groupFile="../../code/group5_20110914.fits"
    acs,group=readCatalogs(acsFile,groupFile)
    acs,group=cleanCatalogs(acs,group,minz,maxz,ztype)

    # assign halo mass and single morph class for each galaxy
    halomass=assignHaloMass(acs,group,ztype)
    morph=assignMorph(acs['ZEST_TYPE'],acs['ZEST_BULGE'])

    options=np.array(["mass_lowz","mass_highz","morph_lowz_highm","morph_highz_highm"])

    for opt in options:
        print opt
        arrs,labels,colors,linestyles,zRange,smRange=defineHistCategories(opt,acs,morph,ztype)
        print "making"
        makeHistPlot(arrs, labels, colors, linestyles, zRange, smRange)
        print "saving"
        plt.savefig("../plots/{}.pdf".format(opt))

