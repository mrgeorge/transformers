#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import fitsio
import matplotlib
matplotlib.use('Agg') # must appear before importing pyplot to get plots w/o GUI
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from scipy import ndimage
from esutil import cosmology
import mrg

cosmo=cosmology.Cosmo(h=0.72,omega_m=0.258,omega_l=0.742) # WMAP5

def readData(imListFile):
    # read data file in IRSA format with image filenames
    names=(('id','i4'),
           ('flag','i2'),
           ('ra','f4'),
           ('dec','f4'),
           ('zphot','f4'),
           ('zbest','f4'),
           ('pmem','f4'),
           ('pmemspec','f4'),
           ('sm','f4'),
           ('mhalo','f4'),
           ('r','f4'),
           ('ztype','i2'),
           ('zbulge','i2'),
           ('nomatch','i2'),
           ('rgflag','i2'),
           ('rgsize','f4'),
           ('rgsersic','f4'),
           ('rgba','f4'),
           ('rgpa','f4'),
           ('color','f4'),
           ('ssfr','f4'),
           ('ebv','f4'),
           ('theta','f4'),
           ('filename','S128')
           )
    data=np.loadtxt(imListFile,skiprows=1,dtype={'names':[x[0] for x in names], 'formats':[x[1] for x in names]})

    return data

def selectData(data, minZ, maxZ, minMh, maxMh, colorSel, morph,zType):

    minP=0.5 # membership probability cut

    if(zType=="zb"):
       zTag="zbest"
       pTag="pmemspec"
    else:
       zTag="zphot"
       pTag="pmem"

    sel=((data['flag'] == 1) &    # only include flag=1
         (data[zTag] > minZ) &    # make a tighter z cut
         (data[zTag] < maxZ) &
         (data[pTag] > minP) &
         (data['mhalo'] > minMh) &
         (data['mhalo'] < maxMh))
    data=data[sel]

    if(morph == "all"):
       sel=range(data.size)
    elif(morph == "spheroidal"):
       sel=((data['ztype'] == 1) |
            ((data['ztype'] == 2) & (data['zbulge'] == 0)))
    elif(morph == "bulge+disk"):
       sel=((data['ztype'] == 2) &
            (data['zbulge'] == 1))
    elif(morph == "latedisk"):
       sel=((data['ztype'] == 2) &
            ((data['zbulge'] == 2) | (data['zbulge'] == 3)))
    elif(morph == "xdisk"):
       sel=((data['ztype'] == 2) &
            (data['zbulge'] == 3))
    elif(morph == "irr"):
       sel=(data['ztype'] == 3)

    data=data[sel]


    if(colorSel =="all"):
       sel=range(data.size)
    elif(colorSel == "red"):
       sel=(data['color'] > 3.5)
    elif(colorSel == "green"):
       sel=((data['color'] > 1.2) &
            (data['color'] < 3.5))
    elif(colorSel == "blue"):
       sel=(data['color'] < 1.2)
    elif(colorSel == "notRed"):
       sel=(data['color'] < 3.5)

    data=data[sel]

    return data

def scatterCentrals(data, minR):
    # assign small negative R values for centrals for display purposes
    cen=data['r'] == 0.
    nCen=len(cen.nonzero()[0])
    offset=np.linspace(0.2*minR,0.8*minR,num=nCen)
    data['r'][cen] += offset

    return data

def setMorphTitle(plt, morph):
    # print morphological category plotted
    if(morph == "all"):
       morphStr="All Morphologies"
    elif(morph == "spheroidal"):
       morphStr="Spheroidal"
    elif(morph == "bulge+disk"):
       morphStr="Bulge+Disk"
    elif(morph == "latedisk"):
       morphStr="Late Disk"
    elif(morph == "xdisk"):
       morphStr="Extreme Disk"
    elif(morph == "irr"):
       morphStr="Irregular"
       
    plt.title(morphStr)

def oplotScaleBar(plt, xBar, yBar, zMed, imSize, imSizePlot, barFontSize):
    # add scale to show galaxy sizes
    dAng=cosmo.Da(0,zMed)
    angStr=str(round(imSize,1))+r'"'
    kpcStr=str(int(round(dAng*(imSize/3600.)*(np.pi/180.)*1000)))+r' kpc'
    barStr=r'{} $\approx$ {}'.format(kpcStr,angStr)

    plt.plot((xBar,xBar-imSizePlot),(yBar,yBar),linestyle='solid',color='black',linewidth=3)
    plt.text(xBar-0.5*imSizePlot,yBar+0.015,barStr,size=barFontSize,horizontalalignment='center',verticalalignment='bottom')

def oplotZRange(plt, xText, yText, minZ, maxZ, zFontSize):
    # add text to show z range
    zStr=r''+str(minZ)+' $<$ z $<$ '+str(maxZ)
    plt.text(xText,yText,zStr,horizontalalignment='right',verticalalignment='top',fontsize=zFontSize)

def getSMLimit(maxZ):
    # return the stellar mass limit at maxZ
    smLimitFile="../../auxfiles/sm_limits.dat"
    z_lim, sm_lim=np.loadtxt(smLimitFile,unpack=True,comments="#")
    thisSMLimit=np.interp(maxZ,z_lim,sm_lim)

    return thisSMLimit

def oplotSMLimit(plt, smLimit, minSM, maxSM):
    # add shaded region to show SM incompleteness
    plt.fill_between((0,1),(smLimit-minSM)/(maxSM-minSM),0,color='gray',alpha=0.5)

def oplotCentralRegion(plt, minR, maxR, minSM, maxSM, smLimit):
    # add shaded region on left to show central area
    plt.fill_betweenx(((smLimit-minSM)/(maxSM-minSM),1),0,(0-minR)/(maxR-minR),color='gray',alpha=0.1)

def oplotMorphText(plt,yMorph,morphFontSize):
    xmin,xmax=plt.gca().get_xlim()
    xleft=xmin+0.03*(xmax-xmin)
    xmid=xmin+0.5*(xmax-xmin)
    xright=xmin+0.97*(xmax-xmin)
    xMorph=np.array([xleft,xmid,xright])
    halignment=np.array(['left','center','right'])
    morphText=np.array(['Spheroidal','Bulge+Disk','Late Disk'])
    for ss in range(morphText.size):
       plt.text(xMorph[ss],yMorph,morphText[ss],horizontalalignment=halignment[ss],verticalalignment='top',fontsize=morphFontSize)


def setupPlot(minZ, zMed, maxZ, imSize, imSizePlot, minR, maxR, minSM, maxSM, morph, thisPanel, nPanels, legendFlag):
    # setup plot of SM vs R with colorbar

    if(nPanels==1):
       lmarg=0.03
       rmarg=0.08
       bmarg=0.1
       tmarg=0.02
       nCols=1
       nRows=1
       rect=np.array([lmarg,bmarg,(1.-rmarg-lmarg),(1.-tmarg-bmarg)])

       tickLabelSize=14
       figSize=(8,6)
       zFontSize='medium'
       barFontSize='x-small'
       morphFontSize='medium'
    else:
       lmarg=0.1
       rmarg=0.01
       bmarg=0.2
       tmarg=0.01
       nCols=np.min([3,nPanels])
       nRows=int(np.ceil(float(nPanels)/nCols))
       thisCol=int(np.remainder(thisPanel,nCols))
       thisRow=int(np.floor(float(thisPanel)/nCols))
       panelWidth=(1.-rmarg-lmarg)/nCols
       panelHeight=(1.-tmarg-bmarg)/nRows
       rect=np.array([lmarg+thisCol*panelWidth,bmarg+(nRows-thisRow-1)*panelHeight,panelWidth,panelHeight])

       tickLabelSize=14
       figSize=(8,9)
       zFontSize=10
       barFontSize=7
       morphFontSize=9

    if(thisPanel==0):
           # use helvetica and latex
           plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
           plt.rc('text', usetex=True)
           plt.rc('axes',linewidth=1.5)
           plt.rc('xtick',labelsize=tickLabelSize)
           plt.rc('ytick',labelsize=tickLabelSize)
           xStr=r'Distance from Group Center [R/R$_{200{\rm c}}$]'
           yStr=r'Stellar Mass [log(M$_{\star}/$M$_{\odot}$)]'
                  
           # start the plot with axes
           fig=plt.figure(1,figsize=figSize)

    # there will be a main figure with the galaxies and a second for the colorbar
    # a set of axes will be created for each panel and for the colorbar
    if(nPanels==1): 
           ax1=plt.axes(rect)

           plt.xlabel(xStr,fontsize='medium')
           plt.ylabel(yStr,fontsize='medium')

    else:
           ax1=plt.axes(rect)

           if(thisPanel==0):
              fig.text(lmarg+0.5*(1.-lmarg-rmarg),0.75*bmarg,xStr,fontsize='medium',horizontalalignment='center',verticalalignment='center',rotation='horizontal')
              fig.text(0.22*lmarg,bmarg+0.5*(1.-bmarg-tmarg),yStr,fontsize='medium',horizontalalignment='center',verticalalignment='center',rotation='vertical')
           
           if(thisCol != 0):
                  plt.setp(ax1.get_yticklabels(),visible=False)
           if((thisRow != nRows) & (thisPanel+nCols < nPanels)):
                  plt.setp(ax1.get_xticklabels(),visible=False)

    if(legendFlag):
       # overlay morph category text
       yMorph=0.97
       oplotMorphText(plt,yMorph,morphFontSize)
       plt.setp(ax1.get_xticklines(),visible=False)
       plt.setp(ax1.get_yticklines(),visible=False)
       plt.minorticks_off()
       axThick=5
       axColor='black'
       ax1.axhline(linewidth=axThick,color=axColor)
       ax1.axvline(linewidth=axThick,color=axColor)
       ax1.axhline(y=1,linewidth=axThick,color=axColor)
       ax1.axvline(x=1,linewidth=axThick,color=axColor)
    else:
       # z range text (upper right)
       xText=0.95
       yText=0.95
       oplotZRange(plt, xText, yText, minZ, maxZ, zFontSize)

       # scale bar (below z range)
       xBar=xText
       yBar=0.82
       oplotScaleBar(plt, xBar, yBar, zMed, imSize, imSizePlot, barFontSize)
    
    # print morph type (title)
#    setMorphTitle(plt, morph)

       # bottom shaded region for SM limit
       smLimit=getSMLimit(maxZ)
       oplotSMLimit(plt, smLimit, minSM, maxSM)

       # left shaded region for centrals
       oplotCentralRegion(plt, minR, maxR, minSM, maxSM, smLimit)

       plt.minorticks_on()

    # axes
    plt.xlim((0,1))
    plt.ylim((0,1))

    xticks=np.linspace(0,1,num=6)
    plt.xticks((xticks-minR)/(maxR-minR),xticks)

    yticks=np.linspace(10,12,num=3)
    plt.yticks((yticks-minSM)/(maxSM-minSM),yticks)

    marg=np.array([lmarg,rmarg,bmarg,tmarg])
    return (plt,ax1,marg)

def prepareImage(img, imSize, fullImSize, angle, floor):
    # clean up the image for plotting: crop, filter, and rotate

    # select inner part of image
    xmid=img.shape[0]/2
    ymid=img.shape[1]/2
    xhalf=0.5*(imSize/fullImSize)*img.shape[0]
    yhalf=0.5*(imSize/fullImSize)*img.shape[1]
    img=img[xmid-xhalf:xmid+xhalf,ymid-yhalf:ymid+yhalf] # cutout image

    # filter the image (gaussian is too blurry, median is ok)
    img=ndimage.median_filter(img,3)

    # rotate the image to align to central
    img=ndimage.rotate(img, angle, reshape=False)

    # filter out other sources in the image
    mask=img > floor
    label_im, nb_labels = ndimage.label(mask)
    if(nb_labels == 0):
        imgBelowFloor=True
        noCentralSource=True
    else:
        imgBelowFloor=False
        mainLabel=label_im[img.shape[0]/2,img.shape[1]/2] # get object number for central source
        if(mainLabel == 0):
            noCentralSource=True
        else:
            noCentralSource=False
            mask2=label_im != mainLabel # select sources other than the main one
            img[mask2]=floor # and set them to the floor value so they aren't plotted

            mask3=img < floor
            img[mask3]=floor

    return (img, imgBelowFloor, noCentralSource)

def threeColorCmap(c0, rgb, c1, pivot, nShades):
    # make a color map that moves from c0 to rgb to c1, switching at pivot
    
    # each row is x, y0, y1
    cdict={'red':  ((  0.0,  c0[0],  c0[0]), # flux=0 -> white
                    (pivot, rgb[0], rgb[0]), # pivot flux -> max color
                    (  1.0,  c1[0],  c1[0])),# flux=1 -> black
           'green':((  0.0,  c0[1],  c0[1]), # flux=0 -> white
                    (pivot, rgb[1], rgb[1]), # pivot flux -> max color
                    (  1.0,  c1[1],  c1[1])),# flux=1 -> black
           'blue': ((  0.0,  c0[2],  c0[2]), # flux=0 -> white
                    (pivot, rgb[2], rgb[2]), # pivot flux -> max color
                    (  1.0,  c1[2],  c1[2]))}# flux=1 -> black

    my_cmap=matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,N=nShades)
    return my_cmap

def normalizeImage(img, floor, cmap, cmin, cmax):
   # convert MxN image to RGBA (MxNx4 RGB+alpha transparency) image
   # renormalized with log scaling
        
   peak=np.max(img)
   normImg=cmap(cmin+(cmax-cmin)*np.log(img/floor)/np.log(peak/floor))

   # set alpha transparency channel (img~floor -> 0, else ~1)
   normImg[:,:,3]=(np.log(img/floor)/np.log(peak/floor))**0.1
   # trim edges by making them transparent
   eps=0.0
   low=img <= (1.+eps)*floor
   normImg[low,3]=0.

   return normImg

def oplotGalaxy(ax1, img, rScale, smScale, floor, imSizePlot, c0, rgb, c1, pivot, cmin, cmax, nShades):
    # plot the image of the galaxy with an appropriate color map
    
    # define color map using RGB color and white/gray gradient
    my_cmap=threeColorCmap(c0,rgb,c1,pivot,nShades)

    # normalize and log-scale image, and get RGBA array
    normImg=normalizeImage(img, floor, my_cmap, cmin, cmax)

    ax1.imshow(normImg,origin='lower',extent=[rScale-0.5*imSizePlot,rScale+0.5*imSizePlot,smScale-0.5*imSizePlot,smScale+0.5*imSizePlot],interpolation='nearest')

def plotGalaxies(data,imDir,fullImSize,imSize,imSizePlot,ax1,minR,maxR,minSM,maxSM,minColor,maxColor,c0,c1,cmap,pivot,cmin,cmax,nShades,legendFlag):
# take list of galaxies to be plotted

       # record stats of plotted galaxies
       nCen=0 # number of centrals plotted
       nShown=0 # total number of galaxies plotted
       nBelowFloor=0 # number not shown because whole image is below floor
       nNoSource=0 # number not shown because a source wasn't found in the middle of image

       for ii in range(data.shape[0]):

           # read and manipulate thumbnail image
           img=fitsio.read(imDir+data['filename'][ii])

           # try to calculate a surface brightness floor based on R_e (not currently used)
#           colMat=np.resize(range(img.shape[0]),(img.shape[0],img.shape[0]))
#           rowMat=colMat.T
#           maxPos=np.unravel_index(np.argmax(img),img.shape)
#           pixDistMat=np.sqrt((rowMat-maxPos[0])**2 + (colMat-maxPos[1])**2)
#           reDeg=mrg.rp2deg(data['rgsize'][ii],cosmo.Da(0,data['z'][ii])*1000.)
#           rePix=reDeg*3600.*img.shape[0]/fullImSize
#           sel=((pixDistMat > 0.4*rePix) & (pixDistMat < 0.6*rePix))
#           floor=np.median(img[sel])
#           print "Re (pixels), Npix, Floor",rePix,len(sel.nonzero()[0]),floor

           floor=0.01
#           floor=0.01*((1.+zMed)/(1.+1.))**(-4)

           img,imgBelowFloor,noSource=prepareImage(img, imSize, fullImSize, data['theta'][ii], floor)

           if(imgBelowFloor):
               nBelowFloor+=1
           elif(noSource):
               nNoSource+=1
           else: # image passes cuts
               nShown+=1
               if(data['r'][ii] <= 0):
                   nCen+=1
                
               rScale=(data['r'][ii]-minR)/(maxR-minR)
               smScale=(data['sm'][ii]-minSM)/(maxSM-minSM)
               colorScale=(data['color'][ii]-minColor)/(maxColor-minColor)

               # define RGB color for this galaxy
               rgb=(cmap(colorScale))[0:3]

               # plot the galaxy
               if(legendFlag):
                  oplotGalaxy(ax1, img, data['r'][ii], data['sm'][ii], floor, imSizePlot, c0, rgb, c1, pivot, cmin, cmax, nShades)
               else:
                  oplotGalaxy(ax1, img, rScale, smScale, floor, imSizePlot, c0, rgb, c1, pivot, cmin, cmax, nShades)
  
       # print plot stats
       print "nShown, nCen, nBelowFloor, nNoSource", nShown, nCen, nBelowFloor, nNoSource

def plotMorphLegend(data,imDir,fullImSize,imSize,imSizePlot,ax1,minR,maxR,minSM,maxSM,minColor,maxColor,c0,c1,cmap,pivot,cmin,cmax,nShades,zType):

    # set up color and morph bins
    colorBins=np.linspace(minColor,maxColor,num=7)
    morphNames=np.array(["spheroidal","bulge+disk","latedisk"])

    # set up 3 columns
    xBuffer=0.09
    yBuffer=0.04
    xLeft,xRight=ax1.get_xlim()
    xWidth=(xRight-xLeft-6*xBuffer)/3.
    xVals=np.array([xLeft+xBuffer+0.5*xWidth,xLeft+3*xBuffer+1.5*xWidth,xLeft+5*xBuffer+2.5*xWidth])

    # trim low mass galaxies
    smCut=10.
    sel=(data['sm'] > smCut)
    data=data[sel]

    for mm in range(morphNames.size):
       for cc in range(colorBins.size-1):
           # select a few random galaxies within data
           nGal=3
           selData=selectData(data, -1, 10, 0, 20, morphNames[mm], "all", zType) # data have already been selected in z and Mh range so just set wide limits
           sel=((selData['color'] > colorBins[cc]) & (selData['color'] <= colorBins[cc+1]))
           selData=selData[sel]

           if(selData.size < nGal):
              print "Not enough galaxies for morphology legend in type=",morphNames[mm],colorBins[cc]
              nGal=selData.size
           if(nGal > 0):
              # plot selected galaxies as normal, just changing positions within axes
              selData[:nGal]['r']=xVals[mm] + np.linspace(-0.5*xWidth,0.5*xWidth,num=nGal)
              selData[:nGal]['sm']=(selData[:nGal]['color']-minColor)/(maxColor-minColor) * (1.-2.*yBuffer)+yBuffer
              plotGalaxies(selData[:nGal],imDir,fullImSize,imSize,imSizePlot,ax1,minR,maxR,minSM,maxSM,minColor,maxColor,c0,c1,cmap,pivot,cmin,cmax,nShades,True)

def oplot2dColorbar(plt, nPanels, marg, nColors, nShades, c0, c1, cmap, pivot, cmin, cmax, minColor, maxColor):
    # add a two-dimensional colorbar

    cbar2d_width=0.6 # sets how big the colorbar looks

    # set up subplot
    if(nPanels==1): # vertical colorbar on the right
       rect=np.array([1.-2*marg[1],marg[2],marg[1],1.-marg[3]-marg[2]])
       ax2=plt.axes(rect)

       # set axes
       ax2.yaxis.tick_right()
       ax2.yaxis.set_label_position("right")
       plt.xlabel(r'$\mu$',fontsize='medium')
       plt.ylabel(r'NUV$ - $R',fontsize='medium')

       plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
       plt.tick_params(axis='y',which='major',left='on')

       plt.xlim((0,cbar2d_width))
       plt.ylim((minColor,maxColor))
    else: # horizontal colorbar on the bottom
       rect=np.array([marg[0]+0.2*(1.-marg[1]-marg[0]),0.1*marg[2],0.6*(1.-marg[1]-marg[0]),0.6*(marg[2])])
       ax2=plt.axes(rect)
    
       # set axes
       plt.xlabel(r'NUV$ - $R',fontsize='medium')
       plt.ylabel(r'$\mu$',fontsize='medium')

       plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
       plt.tick_params(axis='x',which='major',bottom='on')

       plt.xlim((minColor,maxColor))
       plt.ylim((0,cbar2d_width))

    
    # pixel grid for colorbar image
    cbar2d_im=np.linspace(0,nColors*nShades-1,num=nColors*nShades).reshape(nColors,nShades)/(nColors*nShades)

    # color dictionary and index arrays
    # first dim is 3 entries for each row in the colorbar image
    # second dim is x, y0, y1 for each of those 3 coordinates
    # third dim is for r,g,b
    cbar2d_cdict_rgb=np.zeros((3*nColors,3,3))
    ones=np.arange(0,nColors)*3 # left side, white
    twos=ones+1 # pivot x, color
    threes=twos+1 # right side, gray

    # assign RGB values to each color
    rgb_2d=(cmap(1.*np.arange(nColors)/nColors))[:,0:3]

    # fill in color dictionary
    for cc in range(3):
       cbar2d_cdict_rgb[ones,0,cc]=np.linspace(0,nColors-1,num=nColors)/nColors # left side
       cbar2d_cdict_rgb[ones,1,cc]=c0[cc] + cmin/pivot*(rgb_2d[:,cc]-c0[cc]) # white end
       cbar2d_cdict_rgb[ones,2,cc]=c0[cc] + cmin/pivot*(rgb_2d[:,cc]-c0[cc])
       cbar2d_cdict_rgb[twos,0,cc]=cbar2d_cdict_rgb[ones,0,cc]+pivot/nColors # middle pivot
       cbar2d_cdict_rgb[twos,1,cc]=rgb_2d[:,cc] # color
       cbar2d_cdict_rgb[twos,2,cc]=rgb_2d[:,cc]
       cbar2d_cdict_rgb[threes,0,cc]=cbar2d_cdict_rgb[ones,0,cc]+float(nShades-1)/(nShades*nColors) # right side
       cbar2d_cdict_rgb[threes,1,cc]=rgb_2d[:,cc] + (cmax-pivot)/(1.-pivot)*(c1[cc]-rgb_2d[:,cc]) # gray end
       cbar2d_cdict_rgb[threes,2,cc]=rgb_2d[:,cc] + (cmax-pivot)/(1.-pivot)*(c1[cc]-rgb_2d[:,cc])
       cbar2d_cdict_rgb[-1,0,cc]=1. # set last xval=1 by hand

    cbar2d_cdict={'red':cbar2d_cdict_rgb[:,:,0],
                  'green':cbar2d_cdict_rgb[:,:,1],
                  'blue':cbar2d_cdict_rgb[:,:,2]}

    # define color for each of nColors and nShades
    my_cmap_2d=matplotlib.colors.LinearSegmentedColormap('my_colormap_2d',cbar2d_cdict,N=nColors*nShades)

    # plot the colorbar
    if(nPanels==1):
       ax2.imshow(cbar2d_im,origin='lower',extent=[0,cbar2d_width,minColor,maxColor],cmap=my_cmap_2d,interpolation='nearest')
    else:
       ax2.imshow(np.transpose(cbar2d_im),origin='lower',extent=[minColor,maxColor,0,cbar2d_width],cmap=my_cmap_2d,interpolation='nearest')

def main(imDir, imListFile, plotFile, minZ, maxZ, zBin, minMh, maxMh, colorSel, morph, morphLegend, zType):

    # set plot ranges
    minR=-0.15
    maxR=1.05
    minSM=9.3
    maxSM=12.1
    minColor=0. # min M(NUV)-M(R) for colorbar
    maxColor=6. # max M(NUV)-M(R) for colorbar



    # color params
    pivot=0.6    # (0.-1.) sets flux for max color (lower->white, higher->gray)
    grayness=0.5 # (0.-1.) sets how dark the centers of galaxies look (0.->black)
    nColors=256
    nShades=256

    # define red to blue colormap and white/gray colors
    cmap=matplotlib.cm.jet
    c0=np.array((1.,1.,1.)) # RGB for white, left side of colorbar
    c1=grayness*np.array((1.,1.,1.)) # RGB for gray or black, right side of colorbar
    # cmin, cmax set range of shades from white to gray used
    cmin=0.25 # 0 -> pure white, higher to leave some color
    cmax=1.0 # 1 -> pure gray/black, lower to leave some color

    # read data
    fullData=readData(imListFile)

    # determine number of panels for redshift bins and loop over them
    zVals=np.arange(minZ,maxZ,zBin)
    nPanels=zVals.size + morphLegend

    # where should the morphology legend go (if morphLegend==True)
    morphPanel=4
    morphZ=np.array([0.4,0.6])

    for zz in range(nPanels):
       if(morphLegend):
          if(zz < morphPanel):
              thisMinZ=minZ + zz*zBin
              thisMaxZ=thisMinZ + zBin
              legendFlag=False
          elif(zz == morphPanel):
              thisMinZ=morphZ[0]
              thisMaxZ=morphZ[1]
              legendFlag=True
          elif(zz > morphPanel):
              thisMinZ=minZ + (zz-1)*zBin
              thisMaxZ=thisMinZ + zBin
              legendFlag=False
       else:
          thisMinZ=minZ + zz*zBin
          thisMaxZ=thisMinZ + zBin
          legendFlag=False

       # trim full data to redshift range
       data=selectData(fullData, thisMinZ, thisMaxZ, minMh, maxMh, colorSel, morph, zType) # cut on flag=1, z-range, halo mass, color, morph
       if(data.size == 0):
           print "no data in range"
           return
       data=scatterCentrals(data, minR) # add offsets to centrals for display

       if(zType=="zb"):
           zTag="zbest"
       else:
           zTag="zphot"
       zMed=np.median(data[zTag])

       # get image properties
       img=fitsio.read(imDir+data['filename'][0]) # read first one to get params
       fullImSize=20. # fits file is 20" on a side
       imSize=mrg.rp2deg(30.,cosmo.Da(0,zMed)*1000)*3600 # plot a square of this size on a side (30 kpc in arcsec)
       imSizePlot=0.2 # fraction of axis range

       # set up the plot
       plt,ax1,marg=setupPlot(thisMinZ, zMed, thisMaxZ, imSize, imSizePlot, minR, maxR, minSM, maxSM, morph, zz, nPanels, legendFlag)

       if(legendFlag):
           plotMorphLegend(data,imDir,fullImSize,imSize,imSizePlot,ax1,minR,maxR,minSM,maxSM,minColor,maxColor,c0,c1,cmap,pivot,cmin,cmax,nShades,zType)
       else:
           plotGalaxies(data,imDir,fullImSize,imSize,imSizePlot,ax1,minR,maxR,minSM,maxSM,minColor,maxColor,c0,c1,cmap,pivot,cmin,cmax,nShades,legendFlag)


    # end loop over redshift bins / plot panels

    # plot the 2d colorbar
    oplot2dColorbar(plt, nPanels, marg, nColors, nShades, c0, c1, cmap, pivot, cmin, cmax, minColor, maxColor)
    
    # save figure
    plt.savefig(plotFile,dpi=300)


# MAIN - if plotgalSMvR.py is called from command line
if __name__ == '__main__':
    import sys
    if(len(sys.argv)!=13):
        print "Calling sequence: plotgalSMvR.py imDir imListFile plotFile minZ maxZ zBin minMh maxMh colorSel morph morphLegend zType"
        exit

    imDir=sys.argv[1] # str
    imListFile=sys.argv[2] # str
    plotFile=sys.argv[3] # str
    minZ=float(sys.argv[4])
    maxZ=float(sys.argv[5])
    zBin=float(sys.argv[6])
    minMh=float(sys.argv[7])
    maxMh=float(sys.argv[8])
    colorSel=sys.argv[9] # str
    morph=sys.argv[10] # str
    morphLegend=(sys.argv[11] in ["True","true","T","t","TRUE"]) # bool
    zType=sys.argv[12] # str zb for zbest, zp for photoz only

    main(imDir, imListFile, plotFile, minZ, maxZ, zBin, minMh, maxMh, colorSel, morph, morphLegend, zType)
