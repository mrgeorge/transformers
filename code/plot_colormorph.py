#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import fitsio
import numpy as np
import matplotlib
matplotlib.use('Agg') # must appear before importing pyplot to get plots w/o GUI
import matplotlib.pyplot as plt

def cleanCatalogs(acs,group,minz,maxz):

    sel=((acs['KEVIN_MSTAR'] > 2) &
         (acs['ZPHOT'] >= minz) &
         (acs['ZPHOT'] < maxz) &
         (acs['MAG_AUTO'] < 24.2) &
         (acs['GOOD_ZPHOT_LENS'] == 1) &
         (acs['MU_CLASS'] == 1))
    acs=acs[sel]

    sel=((group['FLAG_INCLUDE'] == 1) &
         (group['ZPHOT'] > minz) &
         (group['ZPHOT'] <= maxz))
    group=group[sel]

    return (acs,group)

def assignHaloMass(acs,group):
    halomass=np.zeros(acs.size)
    for gg in range(group.size):
        mem=(acs['GROUP_ID_BEST_SPECZ'] == group['ID'][gg])
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

def census(acs,group,halomass,morph,smbins,zbins,cbins,mbins,rbins,minmh,maxmh):

    nSMbins=smbins.size-1
    nzbins=zbins.size-1
    ncbins=cbins.size-1
    nmbins=mbins.size-1
    nrbins=rbins.size-1

    sat=np.zeros((nSMbins,nzbins,ncbins,nmbins,nrbins))
    cen=np.zeros((nSMbins,nzbins,ncbins,nmbins))
    field=np.zeros_like(cen)

    for sm in range(nSMbins):
        for zz in range(nzbins):
            for cc in range(ncbins):
                for mm in range(nmbins):
                    censel=((halomass >= minmh) &
                            (halomass < maxmh) &
                            (acs['MMGG_SCALE_SPECZ'] == 1) &
                            (acs['GROUP_FLAG_BEST_SPECZ'] == 1) &
                            (acs['KEVIN_MSTAR'] >= smbins[sm]) &
                            (acs['KEVIN_MSTAR'] < smbins[sm+1]) &
                            (acs['ZPHOT'] >= zbins[zz]) &
                            (acs['ZPHOT'] < zbins[zz+1]) &
                            (acs['MNUV_MR'] >= cbins[cc]) &
                            (acs['MNUV_MR'] < cbins[cc+1]) &
                            (morph >= mbins[mm]) &
                            (morph < mbins[mm+1]))
                    fieldsel=((acs['P_MEM_BEST_SPECZ'] <= 0) &
                              (acs['KEVIN_MSTAR'] >= smbins[sm]) &
                              (acs['KEVIN_MSTAR'] < smbins[sm+1]) &
                              (acs['ZPHOT'] >= zbins[zz]) &
                              (acs['ZPHOT'] < zbins[zz+1]) &
                              (acs['MNUV_MR'] >= cbins[cc]) &
                              (acs['MNUV_MR'] < cbins[cc+1]) &
                              (morph >= mbins[mm]) &
                              (morph < mbins[mm+1]))

                    cen[sm][zz][cc][mm]=len(censel.nonzero()[0])
                    field[sm][zz][cc][mm]=len(fieldsel.nonzero()[0])

                    for rr in range(nrbins):
                        satsel=((halomass >= minmh) &
                                (halomass < maxmh) &
                                (acs['MMGG_SCALE_SPECZ'] == 0) &
                                (acs['P_MEM_BEST_SPECZ'] >= 0.5) &
                                (acs['GROUP_FLAG_BEST_SPECZ'] == 1) &
                                (acs['KEVIN_MSTAR'] >= smbins[sm]) &
                                (acs['KEVIN_MSTAR'] < smbins[sm+1]) &
                                (acs['ZPHOT'] >= zbins[zz]) &
                                (acs['ZPHOT'] < zbins[zz+1]) &
                                (acs['MNUV_MR'] >= cbins[cc]) &
                                (acs['MNUV_MR'] < cbins[cc+1]) &
                                (morph >= mbins[mm]) &
                                (morph < mbins[mm+1]) &
                                (acs['DIST_BCG_R200_SPECZ'] >= rbins[rr]) &
                                (acs['DIST_BCG_R200_SPECZ'] < rbins[rr+1]))

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

def bootfrac(numerator,denominator,nSample=500):
    if(numerator > denominator):
        print "Error in bootfrac: numerator > denominator"

    if(denominator == 0):
         #         print "Warning in bootfrac: denominator == 0"
         return -1

    if(denominator < nSample ** (1./denominator)):
        print "Warning in bootfrac: too few elements"
        return -1
        
    arr=np.concatenate((np.ones(numerator),np.zeros(denominator-numerator)))
    sel=np.random.random_integers(0,denominator-1,(nSample,denominator))
    avgArr=arr[sel].sum(1)/denominator
    bootMean=avgArr.mean()
    bootErr=avgArr.std()

    return bootErr

def getCMPop(arr,sm,zz):
    tot=sliceArr(arr,sm=sm,zz=zz,cc=-2,mm=-2,rr=-1)
    rother=sliceArr(arr,sm=sm,zz=zz,cc=1,mm=0,rr=-1)
    rearly=sliceArr(arr,sm=sm,zz=zz,cc=1,mm=1,rr=-1)
    redisk=sliceArr(arr,sm=sm,zz=zz,cc=1,mm=2,rr=-1)
    rldisk=sliceArr(arr,sm=sm,zz=zz,cc=1,mm=3,rr=-1)
    bother=sliceArr(arr,sm=sm,zz=zz,cc=0,mm=0,rr=-1)
    bearly=sliceArr(arr,sm=sm,zz=zz,cc=0,mm=1,rr=-1)
    bedisk=sliceArr(arr,sm=sm,zz=zz,cc=0,mm=2,rr=-1)
    bldisk=sliceArr(arr,sm=sm,zz=zz,cc=0,mm=3,rr=-1)

    return (tot,rother,rearly,redisk,rldisk,bother,bearly,bedisk,bldisk)

def getCMErr(tot,rother,rearly,redisk,rldisk,bother,bearly,bedisk,bldisk):

    vboot=np.vectorize(bootfrac)

    rothererr=vboot(rother,tot)
    rearlyerr=vboot(rearly,tot)
    rediskerr=vboot(redisk,tot)
    rldiskerr=vboot(rldisk,tot)
    bothererr=vboot(bother,tot)
    bearlyerr=vboot(bearly,tot)
    bediskerr=vboot(bedisk,tot)
    bldiskerr=vboot(bldisk,tot)

    return (rothererr,rearlyerr,rediskerr,rldiskerr,bothererr,bearlyerr,bediskerr,bldiskerr)

def plotSatRad(sat,sm,zz,rbins):

    tot,rother,rearly,redisk,rldisk,bother,bearly,bedisk,bldisk=getCMPop(sat,sm,zz)
    rothererr,rearlyerr,rediskerr,rldiskerr,bothererr,bearlyerr,bediskerr,bldiskerr=getCMErr(tot,rother,rearly,redisk,rldisk,bother,bearly,bedisk,bldisk)

    rad=rbins[:-1]+0.5*(rbins[1]-rbins[0])

    # use helvetica and latex
    plt.clf()
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
    plt.rc('text', usetex=True)
    plt.rc('axes',linewidth=1.5)

    plt.xlabel(r'R/R$_{200{\rm c}}$',fontsize='medium')
    plt.ylabel(r'Fraction $|_{\rm M_{\star},z}$',fontsize='medium')
    plt.xlim((-0.05,1.05))
    plt.ylim((-0.01,0.8))

    ms=10
    lw=3
    print tot,rearly/tot,redisk/tot,bldisk/tot
    plt.errorbar(rad,rearly/tot,yerr=rearlyerr,color='red',marker='o',ms=ms,ls='-',lw=lw,label='Red Spheroidal')
    plt.errorbar(rad,redisk/tot,yerr=rediskerr,color='orangered',marker='p',ms=ms,ls='--',lw=lw,label='Red Bulge+Disk')
    plt.errorbar(rad,rldisk/tot,yerr=rldiskerr,color='salmon',marker='<',ms=ms,ls=':',lw=lw,label='Red Disk')
    plt.errorbar(rad,bearly/tot,yerr=bearlyerr,color='darkblue',marker='o',ms=ms,ls='-',lw=lw,label='Blue Spheroidal')
    plt.errorbar(rad,bedisk/tot,yerr=bediskerr,color='dodgerblue',marker='p',ms=ms,ls='--',lw=lw,label='Blue Bulge+Disk')
    plt.errorbar(rad,bldisk/tot,yerr=bldiskerr,color='skyblue',marker='<',ms=ms,ls=':',lw=lw,label='Blue Disk')

def oplotSatCorr(satCorr,sm,zz,rbins):
    tot,rother,rearly,redisk,rldisk,bother,bearly,bedisk,bldisk=getCMPop(satCorr,sm,zz)
    rad=rbins[:-1]+0.5*(rbins[1]-rbins[0])
    ms=100
    print tot,rearly/tot,redisk/tot,bldisk/tot
    plt.scatter(rad,rearly/tot,color='red',marker='o',s=ms)
    plt.scatter(rad,redisk/tot,color='orangered',marker='p',s=ms)
    plt.scatter(rad,rldisk/tot,color='salmon',marker='<',s=ms)
    plt.scatter(rad,bearly/tot,color='darkblue',marker='o',s=ms)
    plt.scatter(rad,bedisk/tot,color='dodgerblue',marker='p',s=ms)
    plt.scatter(rad,bldisk/tot,color='skyblue',marker='<',s=ms)

def oplotField(field,sm,zz):

     tot,rother,rearly,redisk,rldisk,bother,bearly,bedisk,bldisk=getCMPop(field,sm,zz)
     rothererr,rearlyerr,rediskerr,rldiskerr,bothererr,bearlyerr,bediskerr,bldiskerr=getCMErr(tot,rother,rearly,redisk,rldisk,bother,bearly,bedisk,bldisk)

     fieldRad=1.0

     ms=10
     plt.errorbar(fieldRad,rearly/tot,yerr=rearlyerr,color='red',marker='o',ms=ms)
     plt.errorbar(fieldRad,redisk/tot,yerr=rediskerr,color='orangered',marker='p',ms=ms)
     plt.errorbar(fieldRad,rldisk/tot,yerr=rldiskerr,color='salmon',marker='<',ms=ms)
     plt.errorbar(fieldRad,bearly/tot,yerr=bearlyerr,color='darkblue',marker='o',ms=ms)
     plt.errorbar(fieldRad,bedisk/tot,yerr=bediskerr,color='dodgerblue',marker='p',ms=ms)
     plt.errorbar(fieldRad,bldisk/tot,yerr=bldiskerr,color='skyblue',marker='<',ms=ms)

def oplotCen(cen,sm,zz):

     tot,rother,rearly,redisk,rldisk,bother,bearly,bedisk,bldisk=getCMPop(cen,sm,zz)
     rothererr,rearlyerr,rediskerr,rldiskerr,bothererr,bearlyerr,bediskerr,bldiskerr=getCMErr(tot,rother,rearly,redisk,rldisk,bother,bearly,bedisk,bldisk)

     if(np.min(rearlyerr) > 0.): # check that there are enough centrals to get error bars
          cenRad=0.0
          ms=10
          plt.errorbar(cenRad,rearly/tot,yerr=rearlyerr,color='red',marker='o',ms=ms)
          plt.errorbar(cenRad,redisk/tot,yerr=rediskerr,color='orangered',marker='p',ms=ms)
          plt.errorbar(cenRad,rldisk/tot,yerr=rldiskerr,color='salmon',marker='<',ms=ms)
          plt.errorbar(cenRad,bearly/tot,yerr=bearlyerr,color='darkblue',marker='o',ms=ms)
          plt.errorbar(cenRad,bedisk/tot,yerr=bediskerr,color='dodgerblue',marker='p',ms=ms)
          plt.errorbar(cenRad,bldisk/tot,yerr=bldiskerr,color='skyblue',marker='<',ms=ms)

def setTitle(smbins,zbins,sm,zz):
     smstr=r"log(M$_{\star}$/M$_{\odot}$)=["+str(smbins[sm])+", "+str(smbins[sm+1])+")"
     zstr=r"z=["+str(zbins[zz])+", "+str(zbins[zz+1])+")"
     plt.title(r"{}; {}".format(smstr,zstr))

def plotLegend():
#     plt.legend(('Red Spheroidal','Red Bulge+Disk','Red Disk','Blue Spheroidal','Blue Bulge+Disk','Blue Disk'))
     plt.legend(prop={'size':12},numpoints=1,loc=9)

def readCatalogs(acsFile,groupFile):
    acs=fitsio.read(acsFile,ext=1)
    group=fitsio.read(groupFile,ext=1)
    return (acs,group)

def getContamination(rbins):
# calculate purity vs mag and group-centric radius
# then correct for membership contamination using the typical magnitude in each SM/z/c/m bin
# ported from IDL version memstat_separability.pro

     magbins=np.array([17,20,21,22.5,24.2])
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

def getMags(acs,morph,smbins,zbins,cbins,mbins):

     nSMbins=smbins.size-1
     nzbins=zbins.size-1
     ncbins=cbins.size-1
     nmbins=mbins.size-1

     mags=np.zeros((nSMbins,nzbins,ncbins,nmbins))

     for sm in range(nSMbins):
          for zz in range(nzbins):
               for cc in range(ncbins):
                    for mm in range(nmbins):
                         sel=((acs['KEVIN_MSTAR'] >= smbins[sm]) &
                                   (acs['KEVIN_MSTAR'] < smbins[sm+1]) &
                                   (acs['ZPHOT'] >= zbins[zz]) &
                                   (acs['ZPHOT'] < zbins[zz+1]) &
                                   (acs['MNUV_MR'] >= cbins[cc]) &
                                   (acs['MNUV_MR'] < cbins[cc+1]) &
                                   (morph >= mbins[mm]) &
                                   (morph < mbins[mm+1]))
                         mags[sm][zz][cc][mm]=np.median(acs[sel]['MAG_AUTO'])

     return mags

def contaminationCorrection(sat,field,acs,morph,smbins,zbins,cbins,mbins,rbins):

     satCorr=sat.copy() # make a copy to avoid changing sat

     mags=getMags(acs,morph,smbins,zbins,cbins,mbins)
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
# TO DO - check that magVals extend beyond mags so that edges are properly interpolated
                              p=np.interp(mags[sm][zz][cc][mm],magVals[:,rr],purity[:,rr])
                              factor=1.-(1-p)*fieldfrac
                              satCorr[sm][zz][cc][mm][rr] *= factor

     return satCorr
                         
# MAIN - if called from command line
if __name__ == '__main__':

    minz=0.2
    maxz=1.0
    minmh=13.6
    maxmh=14.0

    plotDir="../plots/"

    acsFile="../../code/lensing18_20110914.fits"
    groupFile="../../code/group5_20110914.fits"
    acs,group=readCatalogs(acsFile,groupFile)
    acs,group=cleanCatalogs(acs,group,minz,maxz)

    halomass=assignHaloMass(acs,group)
    morph=assignMorph(acs['ZEST_TYPE'],acs['ZEST_BULGE'])

    smbins=np.array([9.8,10.3,10.7,11.5])
    zbins=np.array([0.2,0.8,1.0])
    complete=np.array([[1,0],[1,1],[1,1]]) # update by hand with smbins, zbins
#    smbins=np.array([10.3,12.0])
#    zbins=np.array([0.2,1.0])
#    complete=np.array([[1]])

    cbins=np.array([-1.,3.5,7.])
    mbins=np.array([-0.5,0.5,1.5,2.5,3.5])
    rbins=np.array([0.01,0.33,0.66,1.0])


    cen,sat,field=census(acs,group,halomass,morph,smbins,zbins,cbins,mbins,rbins,minmh,maxmh)

    satCorr=contaminationCorrection(sat,field,acs,morph,smbins,zbins,cbins,mbins,rbins)

    for sm in range(smbins.size-1):
         for zz in range(zbins.size-1):
              if(complete[sm,zz]==1):
                   plotFile=plotDir+"satrad_sm{}-{}_z{}-{}.pdf".format(smbins[sm],smbins[sm+1],zbins[zz],zbins[zz+1])

                   plotSatRad(sat,sm,zz,rbins)
                   oplotSatCorr(satCorr,sm,zz,rbins)
                   setTitle(smbins,zbins,sm,zz)
                   oplotField(field,sm,zz)
                   oplotCen(cen,sm,zz)
                   plotLegend()
                   plt.savefig(plotFile)

