#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import fitsio
import numpy as np
import matplotlib
matplotlib.use('Agg') # must appear before importing pyplot to get plots w/o GUI
import matplotlib.pyplot as plt
import string

def cleanCatalogs(acs,group,minz,maxz,ztype):

    if(ztype=="zb"):
         zphot="ZPHOT"
    else:
         zphot="PHOTOZ_NON_COMB"

    sel=((acs['KEVIN_MSTAR'] > 2) &
         (acs[zphot] >= minz) &
         (acs[zphot] < maxz) &
         (acs['MAG_AUTO'] < 24.2) &
         (acs['GOOD_ZPHOT_LENS'] == 1) &
         (acs['MU_CLASS'] == 1) &
         (acs['TYPE'] == 0)) # added 9/26/12 to remove masked regions/X-ray sources/star SED types from Ilbert's catalog
    acs=acs[sel]

    sel=((group['FLAG_INCLUDE'] == 1) &
         (group['ZPHOT'] > minz) &
         (group['ZPHOT'] <= maxz))
    group=group[sel]

    return (acs,group)

def assignHaloMass(acs,group,ztype):
    halomass=np.zeros(acs.size)

    if(ztype=="zb"):
         groupID="GROUP_ID_BEST_SPECZ"
    else:
         groupID="GROUP_ID_BEST"
    for gg in range(group.size):
        mem=(acs[groupID] == group['ID'][gg])
        halomass[mem]=group['LENSING_M200'][gg]
        
    return halomass

def assignColour(acs,colourType):
# spelling note: for this array use colour to avoid confusion with e.g. plotting color
# oiQ = extinction-corrected rest-frame template color (Olivier Ilbert quenching def)
# oiC = evolving color cut, no dust correction (Olivier Ilbert color def)
# mrg = NUV-R, R-J color-color cut from quench_mrg catalog (based on KB color-color cuts)

    colour=np.zeros(acs.size)

    if(colourType == "oiQ"):
        blue=acs['MNUV_MR'] < 3.5
        red=acs['MNUV_MR'] >= 3.5
        colour[blue]=0
        colour[red]=1

    elif(colourType == "oiC"):
        blue=(acs['MNUV']-acs['MR'] < 0.5*acs['KEVIN_MSTAR'] - 0.8*acs['PHOTOZ_NON_COMB'] - 0.5)
        red=(acs['MNUV']-acs['MR'] >= 0.5*acs['KEVIN_MSTAR'] - 0.8*acs['PHOTOZ_NON_COMB'] - 0.5)
        bad=((acs['MNUV'] < -30) | (acs['MR'] < -30)) # small number flagged as -99, so not a good color
        colour[blue]=0
        colour[red]=1
        colour[bad]=-1

    elif(colourType == "mrg"):
        blue=(acs['QUENCH_MRG'] == 0)
        red=(acs['QUENCH_MRG'] == 1)
        bad=(acs['QUENCH_MRG'] < 0)
        colour[blue]=0
        colour[red]=1
        colour[bad]=-1

    elif(colourType == "kbQ"):
        blue=(acs['KEVIN_QUENCH_FLAG'] == 0)
        red=(acs['KEVIN_QUENCH_FLAG'] == 1)
        bad=(acs['KEVIN_QUENCH_FLAG'] < 0)
        colour[blue]=0
        colour[red]=1
        colour[bad]=-1

    return colour

def assignMorph(acs,morphType):

    if(morphType == "zest"):
         morph=np.zeros(acs.size)

         early=((acs['ZEST_TYPE'] == 1) |
               ((acs['ZEST_TYPE'] == 2) & (acs['ZEST_BULGE'] == 0)))
         edisk=((acs['ZEST_TYPE'] == 2) &
                (acs['ZEST_BULGE'] == 1))
         ldisk=((acs['ZEST_TYPE'] == 2) &
               ((acs['ZEST_BULGE'] == 2) | (acs['ZEST_BULGE'] ==3)))

         morph[early]=1
         morph[edisk]=2
         morph[ldisk]=3

    elif(morphType == "tasca1"):
         morph=np.zeros(acs.size)

         ell=(acs['MORPH_TASCA1'] == 1)
         spiral=(acs['MORPH_TASCA1'] == 2)
         irr=(acs['MORPH_TASCA1'] == 3)

         morph[ell]=1
         morph[spiral]=3

    elif(morphType == "tasca2"):
         morph=np.zeros(acs.size)

         ell=(acs['MORPH_TASCA2'] == 1)
         spiral=(acs['MORPH_TASCA2'] == 2)
         irr=(acs['MORPH_TASCA2'] == 3)

         morph[ell]=1
         morph[spiral]=3
    elif(morphType == "tasca3"):
         morph=np.zeros(acs.size)

         ell=(acs['MORPH_TASCA3'] == 1)
         spiral=(acs['MORPH_TASCA3'] == 2)
         irr=(acs['MORPH_TASCA3'] == 3)

         morph[ell]=1
         morph[spiral]=3

    elif(morphType == "cassata"):
         morph=np.zeros(acs.size)

         ell=(acs['MORPH_CASSATA'] == 1)
         spiral=(acs['MORPH_CASSATA'] == 2)
         irr=(acs['MORPH_CASSATA'] == 3)

         morph[ell]=1
         morph[spiral]=3

    elif(morphType == "morph2005"):
         morph=np.zeros(acs.size)

         ell=((acs['MORPH2005_GINI'] > 0.43) & (acs['MORPH2005_CONC'] >= 0.3))
         
         # no separate spiral/irr classification, so let's call all others spirals
         morph[:]=3
         
         morph[ell]=1
         
    # irreg and unclassified remain 0
    return morph

def census(acs,group,halomass,colour,morph,smbins,zbins,cbins,mbins,rbins,minmh,maxmh,ztype,smtype,centype):

    nSMbins=smbins.size-1
    nzbins=zbins.size-1
    ncbins=cbins.size-1
    nmbins=mbins.size-1
    nrbins=rbins.size-1

    sat=np.zeros((nSMbins,nzbins,ncbins,nmbins,nrbins))
    cen=np.zeros((nSMbins,nzbins,ncbins,nmbins))
    field=np.zeros_like(cen)

    if(ztype=="zb"):
         mmggScale="MMGG_SCALE_SPECZ"
         groupFlag="GROUP_FLAG_BEST_SPECZ"
         zphot="ZPHOT"
         pMem="P_MEM_BEST_SPECZ"
         distCen="DIST_BCG_R200_SPECZ"
    else: # zp
         mmggScale="MMGG_SCALE"
         groupFlag="GROUP_FLAG_BEST"
         zphot="PHOTOZ_NON_COMB"
         pMem="P_MEM_BEST"
         distCen="DIST_BCG_R200"

    if(centype=="cl"):
        distCen="DIST_CL_R200" # no specz version so just use this for either ztype

    if(smtype=="kb"):
         stellarMass="KEVIN_MSTAR"
    elif(smtype=="oi"):
         stellarMass="ILBERT_MASS_MED"

    for sm in range(nSMbins):
        for zz in range(nzbins):
            for cc in range(ncbins):
                for mm in range(nmbins):
                    censel=((halomass >= minmh) &
                            (halomass < maxmh) &
                            (acs[mmggScale] == 1) &
                            (acs[groupFlag] == 1) &
                            (acs[stellarMass] >= smbins[sm]) &
                            (acs[stellarMass] < smbins[sm+1]) &
                            (acs[zphot] >= zbins[zz]) &
                            (acs[zphot] < zbins[zz+1]) &
                            (colour >= cbins[cc]) &
                            (colour < cbins[cc+1]) &
                            (morph >= mbins[mm]) &
                            (morph < mbins[mm+1]))
                    fieldsel=((acs[pMem] <= 0) &
                              (acs[stellarMass] >= smbins[sm]) &
                              (acs[stellarMass] < smbins[sm+1]) &
                              (acs[zphot] >= zbins[zz]) &
                              (acs[zphot] < zbins[zz+1]) &
                              (colour >= cbins[cc]) &
                              (colour < cbins[cc+1]) &
                              (morph >= mbins[mm]) &
                              (morph < mbins[mm+1]))

                    cen[sm][zz][cc][mm]=len(censel.nonzero()[0])
                    field[sm][zz][cc][mm]=len(fieldsel.nonzero()[0])

                    for rr in range(nrbins):
                        satsel=((halomass >= minmh) &
                                (halomass < maxmh) &
                                (acs[mmggScale] == 0) &
                                (acs[pMem] >= 0.5) &
                                (acs[groupFlag] == 1) &
                                (acs[stellarMass] >= smbins[sm]) &
                                (acs[stellarMass] < smbins[sm+1]) &
                                (acs[zphot] >= zbins[zz]) &
                                (acs[zphot] < zbins[zz+1]) &
                                (colour >= cbins[cc]) &
                                (colour < cbins[cc+1]) &
                                (morph >= mbins[mm]) &
                                (morph < mbins[mm+1]) &
                                (acs[distCen] >= rbins[rr]) &
                                (acs[distCen] < rbins[rr+1]))

                        sat[sm][zz][cc][mm][rr]=len(satsel.nonzero()[0])

    return (cen,sat,field)

def sliceArr(arr,sm=-1,zz=-1,cc=-1,mm=-1,rr=-1):
# return an array of reduced dimensions
# positive inputs slice to include only that bin
# -1 leaves the dimension intact
# -2 sums over the dimension

    ndim=arr.ndim
    nslice=(sm>=0)+(zz>=0)+(cc>=0)+(mm>=0)+(rr>=0)
    nsum=(sm==-2)+(zz==-2)+(cc==-2)+(mm==-2)+(rr==-2)

    if((rr != -1) & (ndim < 5)):
        print "Error in sliceArr: can't slice or sum on r"

    if(nsum == arr.ndim):
        arr=arr.sum()
        return arr
    
    if(rr == -2):
        arr=arr.sum(4)
    if(mm == -2):
        arr=arr.sum(3)
    if(cc == -2):
        arr=arr.sum(2)
    if(zz == -2):
        arr=arr.sum(1)
    if(sm == -2):
        arr=arr.sum(0)

    if(ndim == 4):
        ind="{},{},{},{}".format(sm,zz,cc,mm).replace("-1",":").replace("-2,","").replace(",-2","")
    elif(ndim == 5):
        ind="{},{},{},{},{}".format(sm,zz,cc,mm,rr).replace("-1",":").replace("-2,","").replace(",-2","")

    command="arr=arr[{}]".format(ind)

    exec(command)

    if(ndim-nslice-nsum == 1):
        arr=arr.reshape(arr.size)

    return arr

def bootfrac(frac,denominator,nSample=500):
    numerator=frac*denominator

    if(numerator > denominator):
         print "Error in bootfrac: numerator > denominator"

    if(denominator == 0):
         #         print "Warning in bootfrac: denominator == 0"
         return -1

    if(denominator < nSample ** (1./denominator)):
         #        print "Warning in bootfrac: too few elements"
         return -1
        
    arr=np.concatenate((np.ones(int(numerator)),np.zeros(int(denominator)-int(numerator))))
    sel=np.random.random_integers(0,denominator-1,(nSample,denominator))
    avgArr=arr[sel].sum(1)/denominator
    bootMean=avgArr.mean()
    bootErr=avgArr.std()

    return bootErr

def getCMFracRad(arr,sm,zz,nsplit):
# return a set of slices of the population census that give relevant population fractions
# for nsplit=6, return the fractions for each of the color+morph cuts (vs. R)
# for nsplit=2, return the red fraction and the early type fraction (vs. R)

    tot=sliceArr(arr,sm=sm,zz=zz,cc=-2,mm=-2,rr=-1)

    if(nsplit==6):
         rearly=sliceArr(arr,sm=sm,zz=zz,cc=2,mm=1,rr=-1)
         redisk=sliceArr(arr,sm=sm,zz=zz,cc=2,mm=2,rr=-1)
         rldisk=sliceArr(arr,sm=sm,zz=zz,cc=2,mm=3,rr=-1)
         bearly=sliceArr(arr,sm=sm,zz=zz,cc=1,mm=1,rr=-1)
         bedisk=sliceArr(arr,sm=sm,zz=zz,cc=1,mm=2,rr=-1)
         bldisk=sliceArr(arr,sm=sm,zz=zz,cc=1,mm=3,rr=-1)

         return (1.*rearly/tot,1.*redisk/tot,1.*rldisk/tot,1.*bearly/tot,1.*bedisk/tot,1.*bldisk/tot,tot)

    elif(nsplit==2):
         red=sliceArr(arr,sm=sm,zz=zz,cc=2,mm=-2,rr=-1)
         early=sliceArr(arr,sm=sm,zz=zz,cc=-2,mm=1,rr=-1)

         return (1.*red/tot,1.*early/tot,tot)

    else: # nsplit != 2 or 6
         print "Error in getCMFracRad: nsplit should = 2 or 6"
         return -1

def getCMErr(fracs):

     #    vboot=np.vectorize(bootfrac)

    tot=fracs[-1]
    #    errs=np.zeros((len(fracs)-1,tot.size))
    #    for ii in range(len(fracs)-1):
    #         errs[ii,:]=vboot(fracs[ii],tot)
    #         print tot,fracs[ii],errs[ii,:]

    if(tot.size == 1):
         errs=np.array([bootfrac(fracs[ii],tot) for ii in range(len(fracs)-1)])
    else:
         errs=np.array([[bootfrac(fracs[ii][jj],tot[jj]) for jj in range(tot.size)] for ii in range(len(fracs)-1)])

    return errs

def getCMFracZ(arr,sm,rr,nsplit):
# return a set of slices of the population census that give relevant population fractions
# for nsplit=6, return the fractions for each of the color+morph cuts (vs. z)
# for nsplit=2, return the red fraction and the early type fraction (vs. z)

# rr=0 is taken to mean centrals, rr=-1 for field, else satellites (with rr=index+1 in rbins)
# for indexing arrays, we define rrInd based on the above

    if(rr==0): # centrals
         rrInd=-1 
    elif(rr==-1): # field
         rrInd=-1
    else: # satellites
         rrInd=rr-1

    tot=sliceArr(arr,sm=sm,zz=-1,cc=-2,mm=-2,rr=rrInd)

    if(nsplit==6):
         rearly=sliceArr(arr,sm=sm,zz=-1,cc=2,mm=1,rr=rrInd)
         redisk=sliceArr(arr,sm=sm,zz=-1,cc=2,mm=2,rr=rrInd)
         rldisk=sliceArr(arr,sm=sm,zz=-1,cc=2,mm=3,rr=rrInd)
         bearly=sliceArr(arr,sm=sm,zz=-1,cc=1,mm=1,rr=rrInd)
         bedisk=sliceArr(arr,sm=sm,zz=-1,cc=1,mm=2,rr=rrInd)
         bldisk=sliceArr(arr,sm=sm,zz=-1,cc=1,mm=3,rr=rrInd)

         return (1.*rearly/tot,1.*redisk/tot,1.*rldisk/tot,1.*bearly/tot,1.*bedisk/tot,1.*bldisk/tot,tot)

    elif(nsplit==2):
         red=sliceArr(arr,sm=sm,zz=-1,cc=2,mm=-2,rr=rrInd)
         early=sliceArr(arr,sm=sm,zz=-1,cc=-2,mm=1,rr=rrInd)

         return (1.*red/tot,1.*early/tot,tot)

    else: # nsplit != 2 or 6
         print "Error in getCMFracZ: nsplit should = 2 or 6"
         return -1

def makeBlueSphPlot(zbins,smVals,field):
    plt.clf()
    plt.xlim((9,11.7))
    plt.ylim((0,1))

    colors=np.array(['blue','green','yellow','red'])

    for zz in range(zbins.size-1):
         blueSph=sliceArr(field,sm=-1,zz=zz,cc=1,mm=1,rr=-1)
         allSph=sliceArr(field,sm=-1,zz=zz,cc=-2,mm=1,rr=-1)
         frac=1.*blueSph/allSph

         plt.plot(smVals,frac,color=colors[zz])


def setupPlotArray(nrows,ncols,xtitle=r'Distance from Group Center [R/R$_{200{\rm c}}$]',ytitle=r'Fraction $|_{\rm M_{\star},z}$',figsize=(8,6),xlim=(-0.1,1.19),ylim=(-0.03,0.68)):
    # use helvetica and latex
    plt.clf()
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
    plt.rc('text', usetex=True)
    plt.rc('axes', linewidth=1.5)
    plt.rc('xtick',labelsize=12)
    plt.rc('ytick',labelsize=12)

    # make figure with separate panels for SM and z (or SM and R) bins with shared axes
    fig,axarr=plt.subplots(nrows,ncols,sharex=True,sharey=True,figsize=figsize)
    fig.subplots_adjust(hspace=0,wspace=0)

#    plt.xlabel(r'R/R$_{200{\rm c}}$',fontsize='medium')
#    plt.ylabel(r'Fraction $|_{\rm M_{\star},z}$',fontsize='medium')
    fig.text(0.5,0.02,xtitle,horizontalalignment='center',rotation='horizontal',fontsize='medium')
    fig.text(0.03,0.5,ytitle,verticalalignment='center',rotation='vertical',fontsize='medium')
    plt.xlim(xlim)
    plt.ylim(ylim)

    # hide x ticks for top plots and y ticks for right plots
    plt.setp([[a.get_xticklabels() for a in axarr[b,:]] for b in range(0,nrows-1)], visible=False)
    plt.setp([[a.get_yticklabels() for a in axarr[:,b]] for b in range(1,ncols)], visible=False)

    return axarr

def defineStyles(labelflag):
     if(labelflag == 1):
          labels=np.array(['Red Spheroidal','Red Bulge+Disk','Red Late Disk','Blue Spheroidal','Blue Bulge+Disk','Blue Late Disk'])
     else:
          labels=np.repeat(None,6)
          
     offsets=np.linspace(-1,1,num=6)*0.025
     colors=np.array(['red','orangered','salmon','darkblue','dodgerblue','skyblue'])
     markers=np.tile(['o','p','<'],2)
     linestyles=np.tile(['-','--',':'],2)

     return (labels,offsets,colors,markers,linestyles)

def oplotRad(ax,arr,sm,zz,rad,errflag,labelflag):

     nsplit=6 # how many color morph categories are there
     yvals=getCMFracRad(arr,sm,zz,nsplit)
     if(errflag == 1):
          lw=2
          ms=8
          errs=getCMErr(yvals)
     else:
          lw=0
          ms=6
          errs=np.repeat(None,nsplit)

     labels,offsets,colors,markers,linestyles=defineStyles(labelflag)

     if((np.min(errs) >= 0.) | (errs[0]==None)): # check that there are enough to get error bars
          for ii in range(nsplit):
               ax.errorbar(rad+offsets[ii],yvals[ii],yerr=errs[ii],color=colors[ii],marker=markers[ii],ms=ms,ls=linestyles[ii],lw=lw,label=labels[ii])

     else:
          print "no errs for errflag={},sm={}, zz={}".format(errflag,sm,zz)

def oplotZ(ax,arr,sm,rr,zVals,complete,errflag,labelflag):

     nsplit=6 # how many color morph categories are there
     yvals=getCMFracZ(arr,sm,rr,nsplit)
     if(errflag == 1):
          lw=2
          ms=8
          errs=getCMErr(yvals)
     else:
          lw=0
          ms=6
          errs=np.repeat(None,nsplit)

     labels,offsets,colors,markers,linestyles=defineStyles(labelflag)

     # clean points with too few objects or incomplete SM/z
     yvals=cleanPlot(yvals,errs,complete)

     for ii in range(nsplit):
          ax.errorbar(zVals+offsets[ii],yvals[ii],yerr=errs[ii],color=colors[ii],marker=markers[ii],ms=ms,ls=linestyles[ii],lw=lw,label=labels[ii])

def cleanPlot(yvals,errs,complete):
     nsplit=len(yvals)-1
     nbins=yvals[0].size

     minPop=5

     for jj in range(nbins):
          if((yvals[-1][jj] < minPop) | (complete[jj] == 0)):
               for ii in range(nsplit):
                    yvals[ii][jj]=np.nan

     return yvals

def setSMTitle(ax,smbins,sm):
     smstr=r"log(M$_{\star}$/M$_{\odot}$)=["+str(smbins[sm])+", "+str(smbins[sm+1])+")"
     ax.set_title(smstr,fontsize=12,fontweight='heavy')

def setZTitle(ax,zbins,zz,fontsize):
     zstr=r"z=["+str(zbins[zz])+", "+str(zbins[zz+1])+")"
     xpos=1.05*(ax.get_xlim()[1]-ax.get_xlim()[0])+ax.get_xlim()[0]
     ypos=np.mean(ax.get_ylim())
     ax.text(xpos,ypos,zstr,fontsize=fontsize,rotation=270,fontweight='heavy',verticalalignment='center')

def setRadTitle(ax,rbins,rr,fontsize):
     if(rr==0):
          rstr=r"Centrals"
     elif((rr==rbins.size) | (rr==-1)):
          rstr=r"Field"
     else:
          rstr=r"R/R$_{200{\rm c}}$=["+str(rbins[rr-1])+", "+str(rbins[rr])+")"
     xpos=1.05*(ax.get_xlim()[1]-ax.get_xlim()[0])+ax.get_xlim()[0]
     ypos=np.mean(ax.get_ylim())
     ax.text(xpos,ypos,rstr,fontsize=fontsize,rotation=270,fontweight='heavy',verticalalignment='center')

def plotRadLegend():
#     plt.legend(('Red Spheroidal','Red Bulge+Disk','Red Disk','Blue Spheroidal','Blue Bulge+Disk','Blue Disk'))
     xpad=0.1
     plt.legend(prop={'size':12},numpoints=1,loc='center right', bbox_to_anchor=(-1-xpad,0.5),frameon=False)

def plotBamfordLegend():
     # get list of legend entries and reverse them
     fontsize=9
     handles,labels=plt.gca().get_legend_handles_labels()
     plt.legend(handles[::-1],labels[::-1],title=r'log(M$_{\star}$/M$_{\odot}$) =',prop={'size':fontsize},numpoints=1,loc='lower left',markerscale=0.7,frameon=False)
     plt.gca().get_legend().get_title().set_fontsize(fontsize)

def hidePanel(axarr,zz,sm):
     axarr[zz,sm].set_frame_on(False)
     axarr[zz,sm].get_xaxis().set_visible(False)
     axarr[zz,sm].get_yaxis().set_visible(False)

     # put the axis labels on neighboring panels back
     plt.setp(axarr[zz-1,sm].get_xticklabels(), visible=True)
     plt.setp(axarr[zz,sm+1].get_yticklabels(), visible=True)

def makeRadPlot(plotFile,zbins,smbins,satRad,complete,sat,satCorr,field,cen):
# plot fractions of color/morph types vs R/R200 in separate panels for z and SM bins

    axarr=setupPlotArray(zbins.size-1,smbins.size-1)
    cenRad=0.
    fieldRad=1.1

    for sm in range(smbins.size-1):
         for zz in range(zbins.size-1):
              if(complete[sm,zz]==1):
                   oplotRad(axarr[zz,sm],sat,sm,zz,satRad,1,1)
                   oplotRad(axarr[zz,sm],satCorr,sm,zz,satRad,0,0)
                   oplotRad(axarr[zz,sm],field,sm,zz,fieldRad,1,0)
                   oplotRad(axarr[zz,sm],cen,sm,zz,cenRad,1,0)

                   if(zz==0):
                        setSMTitle(axarr[zz,sm],smbins,sm)
                   if(sm==smbins.size-2):
                        setZTitle(axarr[zz,sm],zbins,zz,12)

              else:
                   hidePanel(axarr,zz,sm)

    plotRadLegend()
    plt.savefig(plotFile)

def makeZPlot(plotFile,smbins,rbins,zVals,complete,sat,satCorr,field,cen):
# plot fractions of color/morph types vs z in separate panels for R and SM bins
     
    axarr=setupPlotArray(rbins.size+1,smbins.size-1,xtitle='Redshift',ytitle=r'Fraction $|_{\rm M_{\star},R/R_{200{\rm c}}}$',figsize=(8,8),xlim=(0.15,1.08),ylim=(-0.03,0.75))

    for sm in range(smbins.size-1):
         oplotZ(axarr[0,sm],cen,sm,0,zVals,complete[sm,:],1,0)
         setSMTitle(axarr[0,sm],smbins,sm)
         for rr in range(1,rbins.size):
              oplotZ(axarr[rr,sm],sat,sm,rr,zVals,complete[sm,:],1,1)
              oplotZ(axarr[rr,sm],satCorr,sm,rr,zVals,complete[sm,:],0,0)
         oplotZ(axarr[-1,sm],field,sm,-1,zVals,complete[sm,:],1,0)


    for rr in range(rbins.size-1+2):
         setRadTitle(axarr[rr,-1],rbins,rr,12)

    plt.savefig(plotFile)

     
def oplotBrokenAxis(axarr,xBreak,xWidth):
     # from http://stackoverflow.com/questions/5656798/python-matplotlib-is-there-a-way-to-make-a-discontinuous-axis

     ds=0.015 # how big to make the diagonal lines in axes coordinates
     # arguments to pass plot, just so we don't keep repeating them
     kwargs = dict(color='black', clip_on=False)

     for ax in axarr.flat: 
          yMin,yMax=ax.axes.get_ylim()
          ax.plot((xBreak-ds,xBreak+ds),(yMin-ds,yMin+ds), **kwargs) # left diagonal
          ax.plot((xBreak+xWidth-ds,xBreak+xWidth+ds),(yMin-ds,yMin+ds), **kwargs) # right diagonal
          ax.plot((xBreak-ds,xBreak+ds),(yMax-ds,yMax+ds), **kwargs) # left diagonal
          ax.plot((xBreak+xWidth-ds,xBreak+xWidth+ds),(yMax-ds,yMax+ds), **kwargs) # right diagonal

def setupBamfordPlot(zbins):
    # use helvetica and latex
    plt.clf()
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
    plt.rc('text', usetex=True)
    plt.rc('axes', linewidth=1.5)
    plt.rc('xtick',labelsize=12)
    plt.rc('ytick',labelsize=12)

    # make figure with separate panels for SM and z bins with shared axes
    fig,axarr=plt.subplots(zbins.size-1,2,sharex=True,sharey=True)
    fig.subplots_adjust(hspace=0,wspace=0)

    fig.text(0.5,0.02,r'Distance from Group Center [R/R$_{200{\rm c}}$]',horizontalalignment='center',rotation='horizontal',fontsize='medium')
    fig.text(0.03,0.5,r'Fraction $|_{\rm M_{\star},z}$',verticalalignment='center',rotation='vertical',fontsize='medium')
    plt.xlim((-0.03,1.19))
    plt.ylim((-0.03,1.03))

    # hide x ticks for top plots and y ticks for right plots
    plt.setp([[a.get_xticklabels() for a in axarr[b,:]] for b in range(0,zbins.size-2)], visible=False)
    plt.setp([a.get_yticklabels() for a in axarr[:,1]], visible=False)

    # overlay "broken axis symbol"
#    xBreak=1.05
#    xWidth=0.03
#    oplotBrokenAxis(axarr,xBreak,xWidth)

    return axarr

def oplotBamford(axs,complete,smbins,arr,zz,rad,errflag,labelflag):

     nsplit=2 # red, elliptical
     nsmbins=smbins.size-1
     offsets=np.linspace(-1,1,num=nsmbins)*0.025
     colors=np.array(['red','green'])
     markers=np.array(['s','o'])
     linestyles=np.array([':','--','-'])
     alphas=np.array([0.5,0.75,1.0])

     for sm in range(nsmbins):
          if(complete[sm,zz]==1):
               yvals=getCMFracRad(arr,sm,zz,nsplit)

               if(errflag == 1):
                    lw=2
                    ms=8
                    errs=getCMErr(yvals)
               else:
                    lw=0
                    ms=6
                    errs=np.repeat(None,nsplit)

               if(labelflag == 1):
                    label="[{}-{})".format(smbins[sm],smbins[sm+1])
               else:
                    label=None

               if((np.min(errs) >= 0.) | (errs[0]==None)): # check that there are enough to get error bars
                    for ii in range(nsplit):
                         axs[ii].errorbar(rad+offsets[sm],yvals[ii],yerr=errs[ii],color=colors[ii],marker=markers[ii],ms=ms,ls=linestyles[sm],lw=lw,alpha=alphas[sm],label=label)

def setBamfordTitle(axs):
     axs[0].set_title('Red Fraction',fontsize='medium',fontweight='heavy')
     axs[1].set_title('Spheroidal Fraction',fontsize='medium',fontweight='heavy')

def makeBamfordPlot(plotFile,zbins,smbins,satRad,complete,sat,satCorr,field,cen):
# plot red frac and early type frac vs R/R200 with separate curves for different SM and separate panels for different z

    axarr=setupBamfordPlot(zbins)
    cenRad=0.
    fieldRad=1.1
    
    for zz in range(zbins.size-1):
         oplotBamford(axarr[zz,:],complete,smbins,sat,zz,satRad,1,1)
         oplotBamford(axarr[zz,:],complete,smbins,satCorr,zz,satRad,0,0)
         oplotBamford(axarr[zz,:],complete,smbins,field,zz,fieldRad,1,0)
         oplotBamford(axarr[zz,:],complete,smbins,cen,zz,cenRad,1,0)

         if(zz==0):
              setBamfordTitle(axarr[zz,:])
         setZTitle(axarr[zz,1],zbins,zz,'medium')

    plt.sca(axarr[0,0])
    plotBamfordLegend()     

    plt.savefig(plotFile)

def readCatalogs(acsFile,groupFile):
    acs=fitsio.read(acsFile,ext=1)
    group=fitsio.read(groupFile,ext=1)
    return (acs,group)

def getContamination(rbins):
# calculate purity vs mag and group-centric radius
# then correct for membership contamination using the typical magnitude in each SM/z/c/m bin
# ported from IDL version memstat_separability.pro

     magbins=np.array([17,20.,20.5,21.0,21.5,22.0,22.5,23.0,23.5,24.0,24.2])
     nmagbins=magbins.size-1
     nrbins=rbins.size-1

     magVals=np.zeros((nmagbins,nrbins))

     totPhotNonMockMem=np.zeros((nmagbins,nrbins))
     totPhotMemMockNon=np.zeros((nmagbins,nrbins))
     totPhotMemMockMem=np.zeros((nmagbins,nrbins))

     mockDir="../../mocks/"
     nMocks=10
     for nm in range(nMocks):

          # open catalogs
          acsFile="{}acs_mock4.{}.fits".format(mockDir,str(nm))
          groupFile="{}group_mock2.{}.fits".format(mockDir,str(nm))
          acs=fitsio.read(acsFile,ext=1)
          group=fitsio.read(groupFile,ext=1)

          # exclude galaxies fainter than 24.2 and stellar mass < 8.5 and z < 0.2
          sel=((acs['MAG_AUTO'] < 24.2) &
               (acs['KEVIN_MSTAR'] > 8.5) &
               (acs['Z'] > 0.2))
          acs=acs[sel]

          nGroups=group.size
          nAcs=acs.size

          # keep track of galaxies that belong in good halos that are NOT in the
          # catalog due to edge effects, etc.
          haloInd=np.zeros(nAcs)
          groupInd=np.zeros(nAcs)
          for ii in range(nAcs):
               sel=(group['ID'] == acs[ii]['HALOID']).nonzero()[0]
               if(len(sel) > 0):
                    haloInd[ii]=sel[0]
               else:
                    haloInd[ii]=-1

               sel=(group['ID'] == acs[ii]['GROUP_ID_BEST'])
               if(len(sel) > 0):
                    groupInd[ii]=sel[0]
               else:
                    groupInd[ii]=-1

          bad=(haloInd == -1)
          good=(haloInd >= 0)
          
          d_halo_r200=np.zeros(nAcs)
          d_halo_r200[bad]=-1.
          d_halo_r200[good]=[np.sqrt((acs[ii]['ALPHA_J2000']-group[haloInd[ii]]['ALPHA_J2000'])**2 + (acs[ii]['DELTA_J2000']-group[haloInd[ii]]['DELTA_J2000'])**2)*3600./group[haloInd[ii]]['LENSING_R200_AS'] for ii in good.nonzero()[0]]

          memThresh=0.5
          for ii in range(nmagbins):
               for jj in range(nrbins):
                    photNonMockMem=((acs['P_MEM_BEST'] < memThresh) &
                                    (acs['MAG_AUTO'] > magbins[ii]) &
                                    (acs['MAG_AUTO'] <= magbins[ii+1]) &
                                    (d_halo_r200 >= rbins[jj]) &
                                    (d_halo_r200 < rbins[jj+1]) &
                                    (haloInd >= 0))
                    photMemMockNon=((acs['P_MEM_BEST'] >= memThresh) &
                                    (acs['MAG_AUTO'] > magbins[ii]) &
                                    (acs['MAG_AUTO'] <= magbins[ii+1]) &
                                    (acs['DIST_BCG_R200'] >= rbins[jj]) &
                                    (acs['DIST_BCG_R200'] < rbins[jj+1]) &
                                    (haloInd == -1))
                    photMemMockMem=((acs['P_MEM_BEST'] >= memThresh) &
                                    (acs['MAG_AUTO'] > magbins[ii]) &
                                    (acs['MAG_AUTO'] <= magbins[ii+1]) &
                                    (acs['DIST_BCG_R200'] >= rbins[jj]) &
                                    (acs['DIST_BCG_R200'] < rbins[jj+1]) &
                                    (d_halo_r200 <= 1) &
                                    (haloInd >= 0))

                    totPhotNonMockMem[ii,jj]+=len(photNonMockMem.nonzero()[0])
                    totPhotMemMockNon[ii,jj]+=len(photMemMockNon.nonzero()[0])
                    totPhotMemMockMem[ii,jj]+=len(photMemMockMem.nonzero()[0])

                    magVals[ii,jj]+=np.mean(acs[photMemMockNon | photMemMockMem]['MAG_AUTO'])

     #end loop over mock lightcones

     magVals=magVals/nMocks

     totPhotMem=totPhotMemMockNon + totPhotMemMockMem
     totMockMem=totPhotNonMockMem + totPhotMemMockMem

     purity=1.*totPhotMemMockMem/totPhotMem
     completeness=1.*totPhotMemMockMem/totMockMem

     return (purity,completeness,magVals)

def getMags(acs,colour,morph,smbins,zbins,cbins,mbins,ztype):

     nSMbins=smbins.size-1
     nzbins=zbins.size-1
     ncbins=cbins.size-1
     nmbins=mbins.size-1

     mags=np.zeros((nSMbins,nzbins,ncbins,nmbins))

     if(ztype=="zb"):
         zphot="ZPHOT"
     else:
         zphot="PHOTOZ_NON_COMB"

     for sm in range(nSMbins):
          for zz in range(nzbins):
               for cc in range(ncbins):
                    for mm in range(nmbins):
                         sel=((acs['KEVIN_MSTAR'] >= smbins[sm]) &
                                   (acs['KEVIN_MSTAR'] < smbins[sm+1]) &
                                   (acs[zphot] >= zbins[zz]) &
                                   (acs[zphot] < zbins[zz+1]) &
                                   (colour >= cbins[cc]) &
                                   (colour < cbins[cc+1]) &
                                   (morph >= mbins[mm]) &
                                   (morph < mbins[mm+1]))
                         mags[sm][zz][cc][mm]=np.median(acs[sel]['MAG_AUTO'])

     return mags

def contaminationCorrection(sat,field,acs,colour,morph,smbins,zbins,cbins,mbins,rbins,ztype):

     satCorr=sat.copy() # make a copy to avoid changing sat

     mags=getMags(acs,colour,morph,smbins,zbins,cbins,mbins,ztype)
     purity,completeness,magVals=getContamination(rbins)

     nSMbins=smbins.size-1
     nzbins=zbins.size-1
     ncbins=cbins.size-1
     nmbins=mbins.size-1
     nrbins=rbins.size-1

     for sm in range(nSMbins):
          for zz in range(nzbins):
               fieldtot=sliceArr(field,sm=sm,zz=zz,cc=-2,mm=-2,rr=-1)
               for cc in range(ncbins):
                    for mm in range(nmbins):
                         fieldpop=sliceArr(field,sm=sm,zz=zz,cc=cc,mm=mm,rr=-1)
                         fieldfrac=1.*fieldpop/fieldtot
                         for rr in range(nrbins):
                              if((mags[sm][zz][cc][mm] > np.max(magVals[:,rr])) |
                                 (mags[sm][zz][cc][mm] < np.min(magVals[:,rr]))):
                                   print mags[sm][zz][cc][mm],sm,zz,cc,mm
                              p=np.interp(mags[sm][zz][cc][mm],magVals[:,rr],purity[:,rr])
                              factor=1.-(1-p)*fieldfrac
                              satCorr[sm][zz][cc][mm][rr] *= factor

     return satCorr
                         
def printCensusTable(censusTableFile,cen,sat,field,zbins,smbins):

     out=open(censusTableFile,'w')

     out.write("% written by printCensusTable in plot_colormorph.py\n")
     out.write("% remove incomplete sections by hand\n")

     nzbins=zbins.size-1
     nsmbins=smbins.size-1

     for zz in range(nzbins):
          out.write("\sidehead{$z={}-{}$}\n".format(zbins[zz],zbins[zz+1]))
          out.write("Centrals & "+ string.join([str(int(cen[sm][zz][0][0])) for sm in range(nsmbins)], " & ") + " \\\\ \n")
          out.write("Satellites & "+ string.join([str(int(sat[sm][zz][0][0][0])) for sm in range(nsmbins)], " & ") + " \\\\ \n")
          out.write("Field & "+ string.join([str(int(field[sm][zz][0][0])) for sm in range(nsmbins)], " & ") + " \\\\ \n")




# MAIN - if called from command line
if __name__ == '__main__':

    minz=0.2
    maxz=1.0
    minmh=13.0
    maxmh=14.0
    ztype="zp" # use zb for zbest (i.e. specz if available), else zp for photoz only
    smtype="kb" # for census selecting in SM range, use kb or oi
    colourType="oiQ" # oiQ, oiC, mrg, kbQ
    morphType="zest" # zest, tasca1, tasca2, tasca3, cassata, morph2005
    centype="cf" # mmgg or cl to test effect of miscentering on radial trend

    plotDir="../plots/"

    # read group and galaxy catalogs
    acsFile="../../code/lensing18_20110914_morph.fits"
    groupFile="../../code/group5_20110914.fits"
    acs,group=readCatalogs(acsFile,groupFile)
    acs,group=cleanCatalogs(acs,group,minz,maxz,ztype)

    # assign halo mass, colour class, and single morph class for each galaxy
    halomass=assignHaloMass(acs,group,ztype)
    colour=assignColour(acs,colourType)
    morph=assignMorph(acs,morphType)

    # setup bins in SM, z, color, morphology, and group-centric radius
    smbins=np.array([9.8,10.3,10.7,11.5])
    zbins=np.array([0.2,0.8,1.0])
    complete=np.array([[1,0],[1,1],[1,1]]) # update by hand with smbins, zbins
    #    zbins=np.array([0.2,0.5,0.8,1.0])
    #   complete=np.array([[1,1,0],[1,1,1],[1,1,1]]) # update by hand with smbins, zbins
    #smbins=np.array([9.6,9.8,10.0,10.2,10.4,10.6,10.8,11.0,11.2,11.4,11.6])
    #zbins=np.array([0.2,0.4,0.6,0.8,1.0])

    cbins=np.array([-2.5,-0.5,0.5,1.5]) # -2=missing, -1=bad, 0=blue,1=red
    mbins=np.array([-0.5,0.5,1.5,2.5,3.5]) # 0=missing/bad, 1=spheroidal, 2=bulge+disk, 3=late disk
    rbins=np.array([0.01,0.33,0.66,1.0])

    satRad=[np.median((rbins[rr],rbins[rr+1])) for rr in range(rbins.size-1)]
    zVals=[np.median((zbins[zz],zbins[zz+1])) for zz in range(zbins.size-1)]
    smVals=[np.median((smbins[sm],smbins[sm+1])) for sm in range(smbins.size-1)]

    # put galaxies in grid of bins
    cen,sat,field=census(acs,group,halomass,colour,morph,smbins,zbins,cbins,mbins,rbins,minmh,maxmh,ztype,smtype,centype)

    # apply corrections to the satellite population for contamination from field galaxies
    satCorr=contaminationCorrection(sat,field,acs,colour,morph,smbins,zbins,cbins,mbins,rbins,ztype)

    # make plot of fraction of color/morph types vs R/R200 with separate panels for each SM, z bin.
    radPlotFile=plotDir+"satrad_{}_{}_{}_{}_mh{}-{}.pdf".format(colourType,morphType,ztype,centype,minmh,maxmh)
    makeRadPlot(radPlotFile,zbins,smbins,satRad,complete,sat,satCorr,field,cen)

    bamfordPlotFile=plotDir+"bamford_{}_{}_{}_{}_mh{}-{}.pdf".format(colourType,morphType,ztype,centype,minmh,maxmh)
    makeBamfordPlot(bamfordPlotFile,zbins,smbins,satRad,complete,sat,satCorr,field,cen)


    # now make different bins for the plot of z-dependence
    
    minmh=13.0
    maxmh=14.0

    smbins=np.array([9.8,10.3,10.7,11.5])
    zbins=np.array([0.2,0.5,0.8,1.0])
    complete=np.array([[1,1,0],[1,1,1],[1,1,1]]) # update by hand with smbins, zbins

    cbins=np.array([-2.5,-0.5,0.5,1.5]) # -2=missing, -1=bad, 0=blue,1=red
    mbins=np.array([-0.5,0.5,1.5,2.5,3.5])
    rbins=np.array([0.01,0.5,1.0])

    satRad=[np.median((rbins[rr],rbins[rr+1])) for rr in range(rbins.size-1)]
    zVals=[np.median((zbins[zz],zbins[zz+1])) for zz in range(zbins.size-1)]

    # put galaxies in grid of bins
    cen,sat,field=census(acs,group,halomass,colour,morph,smbins,zbins,cbins,mbins,rbins,minmh,maxmh,ztype,smtype,centype)

    # apply corrections to the satellite population for contamination from field galaxies
    satCorr=contaminationCorrection(sat,field,acs,colour,morph,smbins,zbins,cbins,mbins,rbins,ztype)


    zPlotFile=plotDir+"satz_{}_{}_{}_{}_mh{}-{}.pdf".format(colourType,morphType,ztype,centype,minmh,maxmh)
    makeZPlot(zPlotFile,smbins,rbins,zVals,complete,sat,satCorr,field,cen)



    # a different binning for census table
    minmh=13.0
    maxmh=14.0

    smbins=np.array([9.8,10.3,10.7,11.5,15.0]) # last bin is to get all centrals >11.5
    zbins=np.array([0.2,0.5,0.8,1.0])

    cbins=np.array([-2.5,1.5]) # all
    mbins=np.array([-0.5,3.5]) # all
    rbins=np.array([0.01,1.0]) # all

    #cen,sat,field=census(acs,group,halomass,colour,morph,smbins,zbins,cbins,mbins,rbins,minmh,maxmh,ztype,smtype)

    #censusTableFile=plotDir+"census_{}_mh{}-{}.tex".format(ztype,minmh,maxmh)
    #printCensusTable(censusTableFile,cen,sat,field,zbins,smbins)
