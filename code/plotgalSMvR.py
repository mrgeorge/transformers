#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import sys
import fitsio
import matplotlib
matplotlib.use('Agg') # must appear before importing pyplot to get plots w/o GUI
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import numpy as np
from scipy import ndimage
from esutil import cosmology


def readData(imListFile):
    # read data file in IRSA format with image filenames
    names=('id','flag','ra','dec','z','sm','mhalo','r','type','bulge','color','theta','filename')
    formats=('i4','i2','f4','f4','f4','f4','f4','f4','i2','i2','f4','f4','S128')
    data=np.loadtxt(imListFile,skiprows=1,dtype={'names':names, 'formats':formats})

    return data

def selectData(data, zMin, zMax, minMh, maxMh):

    sel=((data['flag'] == 1) &    # only include flag=1
         (data['z'] > zMin) &    # make a tighter z cut
         (data['z'] < zMax) &
         (data['mhalo'] > minMh) &
         (data['mhalo'] < maxMh))
    data=data[sel]

    return data

def scatterCentrals(data, minR):
    # assign small negative R values for centrals for display purposes
    cen=data['r'] == 0.
    nCen=len(cen.nonzero()[0])
    offset=np.linspace(0.2*minR,0.8*minR,num=nCen)
    data['r'][cen] += offset

    return data

def oplotScaleBar(plt, xBar, yBar, zMed, imSize, imSizePlot):
    # add scale to show galaxy sizes
    cosmo=cosmology.Cosmo(h=0.72,omega_m=0.258,omega_l=0.742) # WMAP5
    dAng=cosmo.Da(0,zMed)
    barStr=str(imSize)+r'" $\approx$ '+str(int(round(dAng*(imSize/3600.)*(np.pi/180.)*1000)))+' kpc'

    plt.plot((xBar,xBar-imSizePlot),(yBar,yBar),linestyle='solid',color='black',linewidth=3)
    plt.text(xBar-0.5*imSizePlot,yBar+0.015,barStr,size='x-small',horizontalalignment='center',verticalalignment='bottom')

def oplotZRange(plt, xText, yText, zMin, zMax):
    # add text to show z range
    zStr=r''+str(zMin)+' $<$ z $<$ '+str(zMax)
    plt.text(xText,yText,zStr,horizontalalignment='right',verticalalignment='top')

def getSMLimit(zMax):
    # return the stellar mass limit at zMax
    smLimitFile="../../auxfiles/sm_limits.dat"
    z_lim, sm_lim=np.loadtxt(smLimitFile,unpack=True,comments="#")
    thisSMLimit=np.interp(zMax,z_lim,sm_lim)

    return thisSMLimit

def oplotSMLimit(plt, smLimit, minSM, maxSM):
    # add shaded region to show SM incompleteness
    plt.fill_between((0,1),(smLimit-minSM)/(maxSM-minSM),0,color='gray',alpha=0.5)

def oplotCentralRegion(plt, minR, maxR, minSM, maxSM, smLimit):
    # add shaded region on left to show central area
    plt.fill_betweenx(((smLimit-minSM)/(maxSM-minSM),1),0,(0-minR)/(maxR-minR),color='gray',alpha=0.1)

def setupPlot(zMin, zMed, zMax, imSize, imSizePlot, minR, maxR, minSM, maxSM):
    # setup plot of SM vs R with colorbar
    
    # use helvetica and latex
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
    plt.rc('text', usetex=True)
    plt.rc('axes',linewidth=1.5)

    # start the plot with axes
    plt.figure(1)

    # main gridspec to divide figure into a different cells
    # there will be a main figure with the galaxies and a second for the colorbar
    nCells=8 # main will take up a fraction (nCells-1)/nCells of the full width
    gs=gridspec.GridSpec(1,nCells)
    ax1=plt.subplot(gs[:,:-1])

    plt.xlabel(r'R/R$_{200{\rm c}}$',fontsize='medium')
    plt.ylabel(r'log(M$_{\star}$)',fontsize='medium')

    # z range text (upper right)
    xText=0.95
    yText=0.95
    oplotZRange(plt, xText, yText, zMin, zMax)

    # scale bar (below z range)
    xBar=xText
    yBar=0.82
    oplotScaleBar(plt, xBar, yBar, zMed, imSize, imSizePlot)
    
    # bottom shaded region for SM limit
    smLimit=getSMLimit(zMax)
    oplotSMLimit(plt, smLimit, minSM, maxSM)

    # left shaded region for centrals
    oplotCentralRegion(plt, minR, maxR, minSM, maxSM, smLimit)

    # axes
    plt.xlim((0,1))
    plt.ylim((0,1))

    xticks=np.linspace(0,1,num=6)
    plt.xticks((xticks-minR)/(maxR-minR),xticks)

    yticks=np.linspace(10,12,num=3)
    plt.yticks((yticks-minSM)/(maxSM-minSM),yticks)

    plt.minorticks_on()

    return (plt,ax1,gs)

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

def oplotGalaxy(ax1, img, rScale, smScale, floor, imSizePlot, c0, rgb, c1, pivot, nShades):
    # plot the image of the galaxy with an appropriate color map
    
    # define color map using RGB color and white/gray gradient
    my_cmap=threeColorCmap(c0,rgb,c1,pivot,nShades)

    ax1.imshow(np.log((img-floor)/(np.max(img)-floor)),origin='lower',extent=[rScale-0.5*imSizePlot,rScale+0.5*imSizePlot,smScale-0.5*imSizePlot,smScale+0.5*imSizePlot],cmap=my_cmap,alpha=1.)

def oplot2dColorbar(plt, gs, nColors, nShades, c0, c1, cmap, pivot, minColor, maxColor, cbar2d_width):
    # add a two-dimensional colorbar

    # set up subplot
    ax2=plt.subplot(gs[:,-1])

    # pixel grid for colorbar image
    cbar2d_im=np.linspace(0,nColors*nShades-1,num=nColors*nShades).reshape(nColors,nShades)/(nColors*nShades)

    # color dictionary and index arrays
    cbar2d_cdict_r=np.zeros((3*nColors,3))
    cbar2d_cdict_g=np.zeros((3*nColors,3))
    cbar2d_cdict_b=np.zeros((3*nColors,3))
    ones=np.arange(0,nColors)*3 # left side, white
    twos=ones+1 # pivot x, color
    threes=twos+1 # right side, gray

    # assign RGB values to each color
    rgb_2d=np.zeros((nColors,3))
    for cc in range(nColors):
        rgb_2d[cc,:]=matplotlib.colors.colorConverter.to_rgb(cmap(float(cc)/nColors))

    # fill in color dictionary
    # red
    cbar2d_cdict_r[ones,0]=np.linspace(0,nColors-1,num=nColors)/nColors # left side
    cbar2d_cdict_r[ones,1]=c0[0] # white
    cbar2d_cdict_r[ones,2]=c0[0]
    cbar2d_cdict_r[twos,0]=cbar2d_cdict_r[ones,0]+pivot/nColors # middle pivot
    cbar2d_cdict_r[twos,1]=rgb_2d[:,0] # color
    cbar2d_cdict_r[twos,2]=rgb_2d[:,0]
    cbar2d_cdict_r[threes,0]=cbar2d_cdict_r[ones,0]+float(nShades-1)/(nShades*nColors) # right side
    cbar2d_cdict_r[threes,1]=c1[0] # gray
    cbar2d_cdict_r[threes,2]=c1[0]
    cbar2d_cdict_r[-1,0]=1. # set last xval=1 by hand
    # green
    cbar2d_cdict_g[ones,0]=np.linspace(0,nColors-1,num=nColors)/nColors # left side
    cbar2d_cdict_g[ones,1]=c0[1] # white
    cbar2d_cdict_g[ones,2]=c0[1]
    cbar2d_cdict_g[twos,0]=cbar2d_cdict_g[ones,0]+pivot/nColors # middle pivot
    cbar2d_cdict_g[twos,1]=rgb_2d[:,1] # color
    cbar2d_cdict_g[twos,2]=rgb_2d[:,1]
    cbar2d_cdict_g[threes,0]=cbar2d_cdict_g[ones,0]+float(nShades-1)/(nShades*nColors) # right side
    cbar2d_cdict_g[threes,1]=c1[1] # gray
    cbar2d_cdict_g[threes,2]=c1[1]
    cbar2d_cdict_g[-1,0]=1. # set last xval=1 by hand
    # blue
    cbar2d_cdict_b[ones,0]=np.linspace(0,nColors-1,num=nColors)/nColors # left side
    cbar2d_cdict_b[ones,1]=c0[2] # white
    cbar2d_cdict_b[ones,2]=c0[2]
    cbar2d_cdict_b[twos,0]=cbar2d_cdict_b[ones,0]+pivot/nColors # middle pivot
    cbar2d_cdict_b[twos,1]=rgb_2d[:,2] # color
    cbar2d_cdict_b[twos,2]=rgb_2d[:,2]
    cbar2d_cdict_b[threes,0]=cbar2d_cdict_b[ones,0]+float(nShades-1)/(nShades*nColors) # right side
    cbar2d_cdict_b[threes,1]=c1[2] # gray
    cbar2d_cdict_b[threes,2]=c1[2]
    cbar2d_cdict_b[-1,0]=1. # set last xval=1 by hand

    cbar2d_cdict={'red':cbar2d_cdict_r,
                  'green':cbar2d_cdict_g,
                  'blue':cbar2d_cdict_b}

    # define color for each of nColors and nShades
    my_cmap_2d=matplotlib.colors.LinearSegmentedColormap('my_colormap_2d',cbar2d_cdict,N=nColors*nShades)

    # plot the colorbar
    ax2.imshow(cbar2d_im,origin='lower',extent=[0,cbar2d_width,minColor,maxColor],cmap=my_cmap_2d,interpolation='nearest')

    # set axes
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")

    plt.xlabel(r'$\mu$',fontsize='medium')
    plt.ylabel(r'NUV$ - $R',fontsize='medium')

    plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
    plt.tick_params(axis='y',which='major',left='on')

    plt.xlim((0,cbar2d_width))
    plt.ylim((minColor,maxColor))

    
def main(imDir, imListFile, plotFile, zMin, zMax):

    # set plot ranges
    minR=-0.15
    maxR=1.05
    minSM=9.3
    maxSM=12.1
    minColor=0.
    maxColor=6.
    minMh=13.5
    maxMh=14.0

    # read data
    data=readData(imListFile)
    data=selectData(data, zMin, zMax, minMh, maxMh) # cut on flag=1, z-range, halo mass
    if(data.size == 0):
        print "no data in range"
        return
    data=scatterCentrals(data, minR) # add offsets to centrals for display
    zMed=np.median(data['z'])

    # get image properties
    img=fitsio.read(imDir+data['filename'][0]) # read first one to get params
    fullImSize=10. # fits file is 10" on a side
    imSize=5 # plot a square of this size on a side (arcsec)
    imSizePlot=0.2 # fraction of axis range

    # other aesthetics for plotting
    floor=0.01 # noise level, fluxes lower than this are dropped
    cbar2d_width=0.6 # sets how big the colorbar looks

    # color params
    pivot=0.7    # (0.-1.) sets flux for max color (lower->white, higher->gray)
    grayness=0.5 # (0.-1.) sets how dark the centers of galaxies look (0.->black)
    nColors=256
    nShades=256

    # define red to blue colormap and white/gray colors
    cmap=matplotlib.cm.jet
    c0=np.array((1.,1.,1.)) # RGB for white
    c1=grayness*np.array((1.,1.,1.)) # RGB for gray or black

    # setup the plot
    plt,ax1,gs=setupPlot(zMin, zMed, zMax, imSize, imSizePlot, minR, maxR, minSM, maxSM)

    # record stats of plotted galaxies
    nCen=0 # number of centrals plotted
    nShown=0 # total number of galaxies plotted
    nBelowFloor=0 # number not shown because whole image is below floor
    nNoSource=0 # number not shown because a source wasn't found in the middle of image

    for ii in range(data.shape[0]):

        # read and manipulate thumbnail image
        img=fitsio.read(imDir+data['filename'][ii])
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
            rgb=matplotlib.colors.colorConverter.to_rgb(cmap(colorScale))

            # plot the galaxy
            oplotGalaxy(ax1, img, rScale, smScale, floor, imSizePlot, c0, rgb, c1, pivot, nShades)
        #end for loop over galaxies
    # plot the 2d colorbar
    oplot2dColorbar(plt, gs, nColors, nShades, c0, c1, cmap, pivot, minColor, maxColor, cbar2d_width)

    # print plot stats
    print "nShown, nCen, nBelowFloor, nNoSource", nShown, nCen, nBelowFloor, nNoSource
    
    # save figure
    plt.savefig(plotFile)


# MAIN - if plotgalSMvR.py is called from command line
if __name__ == '__main__':
    if(len(sys.argv)!=6):
        print "Calling sequence: plotgalSMvR.py imDir imListFile plotFile zMin zMax"
        exit

    imDir=sys.argv[1]
    imListFile=sys.argv[2]
    plotFile=sys.argv[3]
    zMin=float(sys.argv[4])
    zMax=float(sys.argv[5])

    main(imDir, imListFile, plotFile, zMin, zMax)
    
