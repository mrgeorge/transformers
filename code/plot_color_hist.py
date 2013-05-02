#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

from plot_colormorph import *
from esutil import cosmology

cosmo=cosmology.Cosmo(h=0.72,omega_m=0.258,omega_l=0.742) # WMAP5

# create histogram plots showing color bimodality in different mass, z, morph selections
# overlay dashed line showing 3.5 boundary

def histPlotLegend():
    xpad=0.1
    #    plt.legend(prop={'size':12},numpoints=1,loc='center right',bbox_to_anchor=(-1-xpad,0.5),frameon=False)
    plt.legend(prop={'size':12},numpoints=1,loc='center right',bbox_to_anchor=(-xpad,0.5),frameon=False)
    
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

def makeHistPlot(ax, arrs, labels, colors, hatchstyles, fill, alpha, zRange, smRange, cRange, cBin, logOpt):

    cSplit=3.5
    nCBins=(cRange[1]-cRange[0])/cBin

    lw=2

    #    setupHistPlot()
    #    plt.xlim(cRange)
    [ax.hist(arrs[ii],histtype='step',hatch=hatchstyles[ii],range=cRange,bins=nCBins,label=labels[ii],color=colors[ii],facecolor=colors[ii],lw=lw,log=logOpt,fill=fill[ii],alpha=alpha[ii]) for ii in range(len(arrs))]

    ax.plot((cSplit,cSplit),ax.get_ylim(),color="gray",linestyle=':',lw=lw)
    
    #    histPlotLegend()

    # addVolumeAxis(zRange)

    #    xlim=plt.gca().get_xlim()
    #    ylim=plt.gca().get_ylim()
    #    xPos=xlim[0]+0.05*(xlim[1]-xlim[0])
    #    yPos=ylim[0]+0.95*(ylim[1]-ylim[0])
    #    fontsize=14
    #    oplotZTitle(zRange,xPos,yPos,fontsize)
    #    oplotSMTitle(smRange,xPos,yPos-0.07*(ylim[1]-ylim[0]),fontsize)

def makeScatterPlot(ax,xarrs,yarrs,labels,colors,pointstyles,zRange,smRange):

    [ax.scatter(xarrs[ii],yarrs[ii],color=colors[ii],label=labels[ii],marker=pointstyles[ii],s=1,edgecolor='none') for ii in range(len(xarrs))]

def selectHistData(data,morph,zRange=np.array([0,1]),smRange=np.array([2,15]),morphType="All",zType="zp"):

    if(zType == "zb"): # zbest, use speczs when available
        zphot="ZPHOT"
    else: # zp for photoz only
        zphot="PHOTOZ_NON_COMB"

    morphDict={"Spheroidal":[1],"Bulge+Disk":[2],"Late Disk":[3],"All":[0,1,2,3]}

    sel=((data[zphot] >= zRange[0]) &
         (data[zphot] < zRange[1]) &
         (data['KEVIN_MSTAR'] >= smRange[0]) &
         (data['KEVIN_MSTAR'] < smRange[1]) &
         (map(lambda x: x in morphDict[morphType],morph))
        )

    return sel

def makeHistMultiPlot(histPlotFile,acs,morph,zbins,smbins,labels,colors,hatchstyles,fill,alpha,logOpt):

    cRange=np.array([-1.,6.5])
    cBin=0.5

    if(logOpt=="log"):
        ylim=np.array([9,990])
        logBool=True
    else:
        ylim=np.array([0,499])
        logBool=False

    figsize=((smbins.size-1)*2.5+1,(zbins.size-1)*1.6+1) # default (8,6)
    axarr=setupPlotArray(zbins.size-1,smbins.size-1,xtitle=r"NUV-r$^+$",ytitle="N",xlim=cRange,ylim=ylim,figsize=figsize)

    # loop over sm bins and z bins selecting data and plotting histograms in each panel
    for sm in range(smbins.size-1):
         for zz in range(zbins.size-1):
              if(complete[sm,zz]==1):
                  zRange=zbins[zz:zz+2]
                  smRange=smbins[sm:sm+2]
                  
                  # select arrays
                  sels=[selectHistData(acs,morph,zRange=zRange,smRange=smRange,morphType=morphTypes[ii],zType=ztype) for ii in range(morphTypes.size)]
                  arrs=[acs[sel]['MNUV_MR'] for sel in sels]

                  # plot histograms
                  makeHistPlot(axarr[zz,sm],arrs,labels,colors,hatchstyles,fill,alpha,zRange,smRange,cRange,cBin,logBool)

                  if(zz==0):
                       setSMTitle(axarr[zz,sm],smbins,sm)
                  if(sm==smbins.size-2):
                       setZTitle(axarr[zz,sm],zbins,zz,12,logBool)

              else:
                   hidePanel(axarr,zz,sm)

    # add legend to first panel
    plt.sca(axarr[-1,1]) # change focus to panel adjacent to bottom left to get labels
    histPlotLegend()
                  
    # save figure
    plt.savefig(histPlotFile)

def makeScatterMultiPlot(scatterPlotFile,acs,morph,zbins,smbins,labels,colors,pointstyles):
    
    xlim=np.array([-0.3,1.49])
    ylim=np.array([1.,5.9])

    axarr=setupPlotArray(zbins.size-1,smbins.size-1,xtitle=r"NUV-R",ytitle="R-J",xlim=xlim,ylim=ylim)

    # loop over sm bins and z bins selecting data and plotting histograms in each panel
    for sm in range(smbins.size-1):
         for zz in range(zbins.size-1):
              if(complete[sm,zz]==1):
                  zRange=zbins[zz:zz+2]
                  smRange=smbins[sm:sm+2]
                  
                  # select arrays
                  sels=[selectHistData(acs,morph,zRange=zRange,smRange=smRange,morphType=morphTypes[ii],zType=ztype) for ii in range(morphTypes.size)]
                  xarrs=[(acs[sel]['MR']-acs[sel]['MJ']) for sel in sels]
                  yarrs=[(acs[sel]['MNUV']-acs[sel]['MR']) for sel in sels]

                  # plot color-color scatter diagram
                  makeScatterPlot(axarr[zz,sm],xarrs,yarrs,labels,colors,pointstyles,zRange,smRange)

                  if(zz==0):
                       setSMTitle(axarr[zz,sm],smbins,sm)
                  if(sm==smbins.size-2):
                       setZTitle(axarr[zz,sm],zbins,zz,12,False)

              else:
                   hidePanel(axarr,zz,sm)

    # add legend to first panel
    histPlotLegend()
                  
    # save figure
    plt.savefig(scatterPlotFile)




# MAIN - if called from command line
if __name__ == '__main__':

    minz=0.2
    maxz=1.0
    minmh=13.0
    maxmh=14.0
    ztype="zp" # use zb for zbest (i.e. specz if available), else zp for photoz only

    plotDir="../plots/"
    plotFmt="pdf" # pdf or eps

    # read group and galaxy catalogs
    acsFile="../../code/lensing18_20110914.fits"
    groupFile="../../code/group5_20110914.fits"
    acs,group=readCatalogs(acsFile,groupFile)
    acs,group=cleanCatalogs(acs,group,minz,maxz,ztype)

    # assign halo mass and single morph class for each galaxy
    halomass=assignHaloMass(acs,group,ztype)
    morph=assignMorph(acs,"zest")

    #    smbins=np.array([9.8,10.3,10.7,11.5])
    #    zbins=np.array([0.2,0.5,0.8,1.0])
    #    complete=np.array([[1,1,0],[1,1,1],[1,1,1]]) # update by hand with smbins, zbins
    smbins=np.array([9.8,10.3,10.8])
    zbins=np.array([0.2,0.5,0.8,1.0])
    complete=np.array([[1,1,0],[1,1,1]]) # update by hand with smbins, zbins

    morphTypes=np.array(["Late Disk","Bulge+Disk","Spheroidal"])

    labels=[morphTypes[ii] for ii in range(morphTypes.size)]
    #    colors=["blue","green","red"]
    colors=["indigo","goldenrod","gray"]
#    linestyles=["-","--",":"] # ["solid","dashed","dotted"]
    hatchstyles=['/','\\','']
    pointstyles=['o','o','o']
    fill=np.array([False,False,True])
    alpha=np.array([1,1,0.3])

    # options for multi-panel hist plot
    logOpt="log"
    histPlotFile=plotDir+"color_hist_{}.{}".format(logOpt,plotFmt)

    makeHistMultiPlot(histPlotFile,acs,morph,zbins,smbins,labels,colors,hatchstyles,fill,alpha,logOpt)

    # options for multi-panel color-color scatter plot
    scatterPlotFile=plotDir+"color_scatter.{}".format(plotFmt)

#    makeScatterMultiPlot(scatterPlotFile,acs,morph,zbins,smbins,labels,colors,pointstyles)
