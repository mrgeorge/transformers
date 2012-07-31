#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import fitsio
import numpy as np

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

    if(rr >= 0 & arr.ndim < 5):
        print "Error in sliceArr: can't cut on r"

    if(sm >= 0):
        arr=arr[sm]
        if(zz >= 0):
            arr=arr[zz]
            if(cc >= 0):
                arr=arr[cc]
                if(mm >= 0):
                    arr=arr[mm]
                    if(rr >= 0):
                        arr=arr[rr]
                # mm == -1
                elif(rr >= 0):
                    arr=arr[:][rr]
            # cc == -1
            elif(mm >= 0):
                arr=arr[:][mm]
                if(rr >= 0):
                    arr=arr[:][rr]
            # cc,mm == -1
            elif(rr >= 0):
                arr=arr[:][:][rr]
        # zz == -1
        elif(cc >= 0):
            arr=arr[:][cc]
            if(mm >= 0):
                arr=arr[:][mm]
                if(rr >= 0):
                    arr=arr[:][rr]
            # zz,mm == -1
            elif(rr >= 0):
                arr=arr[:][:][rr]
        # zz,cc == -1
        elif(mm >= 0):
            arr=arr[:][:][mm]
            if(rr >= 0):
                arr=arr[:][:][rr]
        # zz,cc,mm == -1
        elif(rr >= 0):
            arr=arr[:][:][:][rr]
    # sm == -1
    elif(zz >= 0):
        arr=arr[:][zz]
        if(cc >= 0):
            arr=arr[:][cc]
            if(mm >= 0):
                arr=arr[:][mm]
                if(rr >= 0):
                    arr=arr[:][rr]
            # sm,mm == -1
            elif(rr >= 0):
                arr=arr[:][:][rr]
        # sm,cc == -1
        elif(mm >= 0):
            arr=arr[:][:][mm]
            if(rr >= 0):
                arr=arr[:][:][rr]
        # sm,cc,mm == -1
        elif(rr >= 0):
            arr=arr[:][:][:][rr]
    # sm,zz == -1
    elif(cc >= 0):
        arr=arr[:][:][cc]
        if(mm >= 0):
            arr=arr[:][:][mm]
            if(rr >= 0):
                arr=arr[:][:][rr]
        # sm,zz,mm == -1
        elif(rr >= 0):
            arr=arr[:][:][:][rr]
    # sm,zz,cc == -1
    elif(mm >= 0):
        arr=arr[:][:][:][mm]
        if(rr >= 0):
            arr=arr[:][:][:][rr]
    # sm,zz,cc,mm == -1
    elif(rr >= 0):
        arr=arr[:][:][:][:][rr]

    return arr


# MAIN - if called from command line
if __name__ == '__main__':

    minz=0.2
    maxz=1.0
    minmh=13.6
    maxmh=14.0

    acs=fitsio.read("~/data/cosmos/code/lensing18_20110914.fits",ext=1)
    group=fitsio.read("~/data/cosmos/code/group5_20110914.fits",ext=1)
    acs,group=cleanCatalogs(acs,group,minz,maxz)

    halomass=assignHaloMass(acs,group)
    morph=assignMorph(acs['ZEST_TYPE'],acs['ZEST_BULGE'])

    smbins=np.array([9.8,10.3,10.7,11.5])
    zbins=np.array([0.2,0.8,1.0])
    cbins=np.array([-1.,3.5,7.])
    mbins=np.array([0.5,1.5,2.5,3.5])
    rbins=np.array([0.01,0.33,0.66,1.0])

    complete=np.array([[1,0],[1,1],[1,1]]) # update by hand with smbins, zbins

    cen,sat,field=census(acs,group,halomass,morph,smbins,zbins,cbins,mbins,rbins,minmh,maxmh)
    
    
    
