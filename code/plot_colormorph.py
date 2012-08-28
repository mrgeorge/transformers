#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import fitsio
import numpy as np
import matplotlib
matplotlib.use('Agg') # must appear before importing pyplot to get plots w/o GUI
import matplotlib.pyplot as plt

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
         (acs['MU_CLASS'] == 1))
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

def assignMorph(zestType,zestBulge):

    morph=np.zeros(zestType.size)

    early=((zestType == 1) |
         ((zestType == 2) & (zestBulge == 0)))
    edisk=((zestType == 2) &
           (zestBulge == 1))
    ldisk=((zestType == 2) &
           ((zestBulge == 2) | (zestBulge ==3)))

    morph[early]=1
    morph[edisk]=2
    morph[ldisk]=3

    # irreg and unclassified remain 0
    return morph

def census(acs,group,halomass,morph,smbins,zbins,cbins,mbins,rbins,minmh,maxmh,ztype):

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
         distBCG="DIST_BCG_R200_SPECZ"
    else:
         mmggScale="MMGG_SCALE"
         groupFlag="GROUP_FLAG_BEST"
         zphot="PHOTOZ_NON_COMB"
         pMem="P_MEM_BEST"
         distBCG="DIST_BCG_R200"

    for sm in range(nSMbins):
        for zz in range(nzbins):
            for cc in range(ncbins):
                for mm in range(nmbins):
                    censel=((halomass >= minmh) &
                            (halomass < maxmh) &
                            (acs[mmggScale] == 1) &
                            (acs[groupFlag] == 1) &
                            (acs['KEVIN_MSTAR'] >= smbins[sm]) &
                            (acs['KEVIN_MSTAR'] < smbins[sm+1]) &
                            (acs[zphot] >= zbins[zz]) &
                            (acs[zphot] < zbins[zz+1]) &
                            (acs['MNUV_MR'] >= cbins[cc]) &
                            (acs['MNUV_MR'] < cbins[cc+1]) &
                            (morph >= mbins[mm]) &
                            (morph < mbins[mm+1]))
                    fieldsel=((acs[pMem] <= 0) &
                              (acs['KEVIN_MSTAR'] >= smbins[sm]) &
                              (acs['KEVIN_MSTAR'] < smbins[sm+1]) &
                              (acs[zphot] >= zbins[zz]) &
                              (acs[zphot] < zbins[zz+1]) &
                              (acs['MNUV_MR'] >= cbins[cc]) &
                              (acs['MNUV_MR'] < cbins[cc+1]) &
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
                                (acs['KEVIN_MSTAR'] >= smbins[sm]) &
                                (acs['KEVIN_MSTAR'] < smbins[sm+1]) &
                                (acs[zphot] >= zbins[zz]) &
                                (acs[zphot] < zbins[zz+1]) &
                                (acs['MNUV_MR'] >= cbins[cc]) &
                                (acs['MNUV_MR'] < cbins[cc+1]) &
                                (morph >= mbins[mm]) &
                                (morph < mbins[mm+1]) &
                                (acs[distBCG] >= rbins[rr]) &
                                (acs[distBCG] < rbins[rr+1]))

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

def getCMFrac(arr,sm,zz,nsplit):
# return a set of slices of the population census that give relevant population fractions
# for nsplit=6, return the fractions for each of the color+morph cuts (vs. R)
# for nsplit=2, return the red fraction and the early type fraction (vs. R)

    tot=sliceArr(arr,sm=sm,zz=zz,cc=-2,mm=-2,rr=-1)

    if(nsplit==6):
         rearly=sliceArr(arr,sm=sm,zz=zz,cc=1,mm=1,rr=-1)
         redisk=sliceArr(arr,sm=sm,zz=zz,cc=1,mm=2,rr=-1)
         rldisk=sliceArr(arr,sm=sm,zz=zz,cc=1,mm=3,rr=-1)
         bearly=sliceArr(arr,sm=sm,zz=zz,cc=0,mm=1,rr=-1)
         bedisk=sliceArr(arr,sm=sm,zz=zz,cc=0,mm=2,rr=-1)
         bldisk=sliceArr(arr,sm=sm,zz=zz,cc=0,mm=3,rr=-1)

         return (1.*rearly/tot,1.*redisk/tot,1.*rldisk/tot,1.*bearly/tot,1.*bedisk/tot,1.*bldisk/tot,tot)

    elif(nsplit==2):
         red=sliceArr(arr,sm=sm,zz=zz,cc=1,mm=-2,rr=-1)
         early=sliceArr(arr,sm=sm,zz=zz,cc=-2,mm=1,rr=-1)

         return (1.*red/tot,1.*early/tot,tot)

    else: # nsplit != 2 or 6
         print "Error in getCMFrac: nsplit should = 2 or 6"
         return -1

def getCMErr(fracs):

    vboot=np.vectorize(bootfrac)

    tot=fracs[-1]
    errs=np.zeros((len(fracs)-1,tot.size))
    for ii in range(len(fracs)-1):
         errs[ii,:]=vboot(fracs[ii],tot)

    return errs

def setupRadPlot(zbins,smbins):
    # use helvetica and latex
    plt.clf()
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
    plt.rc('text', usetex=True)
    plt.rc('axes', linewidth=1.5)
    plt.rc('xtick',labelsize=12)
    plt.rc('ytick',labelsize=12)

    # make figure with separate panels for SM and z bins with shared axes
    fig,axarr=plt.subplots(zbins.size-1,smbins.size-1,sharex=True,sharey=True)
    fig.subplots_adjust(hspace=0,wspace=0)

#    plt.xlabel(r'R/R$_{200{\rm c}}$',fontsize='medium')
#    plt.ylabel(r'Fraction $|_{\rm M_{\star},z}$',fontsize='medium')
    fig.text(0.5,0.02,r'R/R$_{200{\rm c}}$',horizontalalignment='center',rotation='horizontal',fontsize='medium')
    fig.text(0.03,0.5,r'Fraction $|_{\rm M_{\star},z}$',verticalalignment='center',rotation='vertical',fontsize='medium')
    plt.xlim((-0.1,1.1))
    plt.ylim((-0.03,0.68))

    # hide x ticks for top plots and y ticks for right plots
    plt.setp([[a.get_xticklabels() for a in axarr[b,:]] for b in range(0,zbins.size-2)], visible=False)
    plt.setp([[a.get_yticklabels() for a in axarr[:,b]] for b in range(1,smbins.size-1)], visible=False)

    return axarr

def oplotRad(ax,arr,sm,zz,rad,errflag,labelflag):

     nsplit=6
     yvals=getCMFrac(arr,sm,zz,nsplit)
     if(errflag == 1):
          lw=2
          ms=8
          errs=getCMErr(yvals)
     else:
          lw=0
          ms=6
          errs=np.repeat(None,nsplit)

     if(labelflag == 1):
          labels=np.array(['Red Spheroidal','Red Bulge+Disk','Red Disk','Blue Spheroidal','Blue Bulge+Disk','Blue Disk'])
     else:
          labels=np.repeat(None,nsplit)
          
     offsets=np.linspace(-1,1,num=nsplit)*0.025
     colors=np.array(['red','orangered','salmon','darkblue','dodgerblue','skyblue'])
     markers=np.tile(['o','p','<'],2)
     linestyles=np.tile(['-','--',':'],2)

     if((np.min(errs) >= 0.) | (errs[0]==None)): # check that there are enough to get error bars

          for ii in range(nsplit):
               ax.errorbar(rad+offsets[ii],yvals[ii],yerr=errs[ii],color=colors[ii],marker=markers[ii],ms=ms,ls=linestyles[ii],lw=lw,label=labels[ii])

     else:
          print "no errs for errflag={},sm={}, zz={}".format(errflag,sm,zz)

def setSMTitle(ax,smbins,sm):
     smstr=r"log(M$_{\star}$/M$_{\odot}$)=["+str(smbins[sm])+", "+str(smbins[sm+1])+")"
     ax.set_title(smstr,fontsize=12,fontweight='heavy')

def setZTitle(ax,zbins,zz,fontsize):
     zstr=r"z=["+str(zbins[zz])+", "+str(zbins[zz+1])+")"
     xpos=1.05*(ax.get_xlim()[1]-ax.get_xlim()[0])+ax.get_xlim()[0]
     ypos=np.mean(ax.get_ylim())
     ax.text(xpos,ypos,zstr,fontsize=fontsize,rotation=270,fontweight='heavy',verticalalignment='center')

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

def makeRadPlot(plotFile,zbins,smbins,complete,sat,satCorr,field,cen):
# plot fractions of color/morph types vs R/R200 in separate panels for z and SM bins

    axarr=setupRadPlot(zbins,smbins)
    for sm in range(smbins.size-1):
         for zz in range(zbins.size-1):
              if(complete[sm,zz]==1):
                   oplotRad(axarr[zz,sm],sat,sm,zz,satRad,1,1)
                   oplotRad(axarr[zz,sm],satCorr,sm,zz,satRad,0,0)
                   oplotRad(axarr[zz,sm],field,sm,zz,1.,1,0)
                   oplotRad(axarr[zz,sm],cen,sm,zz,0.,1,0)

                   if(zz==0):
                        setSMTitle(axarr[zz,sm],smbins,sm)
                   if(sm==smbins.size-2):
                        setZTitle(axarr[zz,sm],zbins,zz,12)

              else:
                   hidePanel(axarr,zz,sm)

    plotRadLegend()
    plt.savefig(plotFile)

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

    fig.text(0.5,0.02,r'R/R$_{200{\rm c}}$',horizontalalignment='center',rotation='horizontal',fontsize='medium')
    fig.text(0.03,0.5,r'Fraction $|_{\rm M_{\star},z}$',verticalalignment='center',rotation='vertical',fontsize='medium')
    plt.xlim((-0.1,1.1))
    plt.ylim((-0.03,1.03))

    # hide x ticks for top plots and y ticks for right plots
    plt.setp([[a.get_xticklabels() for a in axarr[b,:]] for b in range(0,zbins.size-2)], visible=False)
    plt.setp([a.get_yticklabels() for a in axarr[:,1]], visible=False)

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
               yvals=getCMFrac(arr,sm,zz,nsplit)

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

def makeBamfordPlot(plotFile,zbins,smbins,complete,sat,satCorr,field,cen):
# plot red frac and early type frac vs R/R200 with separate curves for different SM and separate panels for different z

    axarr=setupBamfordPlot(zbins)
    for zz in range(zbins.size-1):
         oplotBamford(axarr[zz,:],complete,smbins,sat,zz,satRad,1,1)
         oplotBamford(axarr[zz,:],complete,smbins,satCorr,zz,satRad,0,0)
         oplotBamford(axarr[zz,:],complete,smbins,field,zz,1.,1,0)
         oplotBamford(axarr[zz,:],complete,smbins,cen,zz,0.,1,0)

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

def getMags(acs,morph,smbins,zbins,cbins,mbins,ztype):

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
                                   (acs['MNUV_MR'] >= cbins[cc]) &
                                   (acs['MNUV_MR'] < cbins[cc+1]) &
                                   (morph >= mbins[mm]) &
                                   (morph < mbins[mm+1]))
                         mags[sm][zz][cc][mm]=np.median(acs[sel]['MAG_AUTO'])

     return mags

def contaminationCorrection(sat,field,acs,morph,smbins,zbins,cbins,mbins,rbins,ztype):

     satCorr=sat.copy() # make a copy to avoid changing sat

     mags=getMags(acs,morph,smbins,zbins,cbins,mbins,ztype)
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
                         
# MAIN - if called from command line
if __name__ == '__main__':

    minz=0.2
    maxz=1.0
    minmh=13.0
    maxmh=14.0
    ztype="zb" # use zb for zbest (i.e. specz if available), else zp for photoz only

    plotDir="../plots/"

    # read group and galaxy catalogs
    acsFile="../../code/lensing18_20110914.fits"
    groupFile="../../code/group5_20110914.fits"
    acs,group=readCatalogs(acsFile,groupFile)
    acs,group=cleanCatalogs(acs,group,minz,maxz,ztype)

    # assign halo mass and single morph class for each galaxy
    halomass=assignHaloMass(acs,group,ztype)
    morph=assignMorph(acs['ZEST_TYPE'],acs['ZEST_BULGE'])

    # setup bins in SM, z, color, morphology, and group-centric radius
    smbins=np.array([9.8,10.3,10.7,11.5])
    zbins=np.array([0.2,0.8,1.0])
    complete=np.array([[1,0],[1,1],[1,1]]) # update by hand with smbins, zbins
#    smbins=np.array([10.3,12.0])
#    zbins=np.array([0.2,1.0])
#    complete=np.array([[1]])

    cbins=np.array([-1.,3.5,7.])
    mbins=np.array([-0.5,0.5,1.5,2.5,3.5])
    rbins=np.array([0.01,0.33,0.66,1.0])

    satRad=rbins[:-1]+0.5*(rbins[1]-rbins[0])

    # put galaxies in grid of bins
    cen,sat,field=census(acs,group,halomass,morph,smbins,zbins,cbins,mbins,rbins,minmh,maxmh,ztype)

    # apply corrections to the satellite population for contamination from field galaxies
    satCorr=contaminationCorrection(sat,field,acs,morph,smbins,zbins,cbins,mbins,rbins,ztype)

    # make plot of fraction of color/morph types vs R/R200 with separate panels for each SM, z bin.
    radPlotFile=plotDir+"satrad_{}_mh{}-{}.pdf".format(ztype,minmh,maxmh)
    makeRadPlot(radPlotFile,zbins,smbins,complete,sat,satCorr,field,cen)

    bamfordPlotFile=plotDir+"bamford_{}_mh{}-{}.pdf".format(ztype,minmh,maxmh)
    makeBamfordPlot(bamfordPlotFile,zbins,smbins,complete,sat,satCorr,field,cen)
    
