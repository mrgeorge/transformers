#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

# Make a schematic/cartoon plot of blue late disks and red bulge+disks for
# an ensemble cluster, similar to plotgalArt but to illustrate the main trend
# from the transformers paper

from plotgalSMvR import *
import scipy.stats # for ranking data


imDir="../images/twenty_arcsec/"
imListFile=imDir+"acs_mem.list"

minZ=0.2
maxZ=0.5
zMed=np.median([minZ,maxZ])
minSM=9.8
maxSM=10.3
minMh=13.0
maxMh=14.0
zType="zp"

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

# get data
fullData=readData(imListFile)

dataCen=selectData(fullData,minZ,maxZ,minMh,maxMh,"all","all",zType)
sel=(dataCen["r"] == 0)
dataCen=dataCen[sel]
dataCen=dataCen[np.random.randint(dataCen.size)]

sel=((fullData['sm'] > minSM) & (fullData['sm'] < maxSM))
fullData=fullData[sel]

colorSel="blue"
morph="latedisk"
dataBlue=selectData(fullData, minZ, maxZ, minMh, maxMh, colorSel,morph, zType)

colorSel="red"
morph="bulge+disk"
dataRed=selectData(fullData, minZ, maxZ, minMh, maxMh, colorSel,morph, zType)


# setup plot
fig=plt.figure(1)
plt.clf()
plt.xlim((-0.1,1.1))
plt.ylim((-0.1,1.1))
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
plt.rc('text', usetex=True)
plt.rc('axes',linewidth=1.5)
ax1=plt.gca()
ax1.set_frame_on(False)
ax1.get_xaxis().set_visible(False)
ax1.get_yaxis().set_visible(False)
ax2=fig.add_axes(ax1)

img=fitsio.read(imDir+dataBlue['filename'][0]) # read first one to get params
fullImSize=20. # fits file is 20" on a side
imSize=mrg.rp2deg(30.,cosmo.Da(0,zMed)*1000)*3600 # plot a square of this size on a side (30 kpc in arcsec)
imSizePlot=0.2 # fraction of axis range


# convert coordinates
    
# for a polar plot showing distance from group center and azimuthal angle
dataBlueX=dataBlue['r']*np.cos(dataBlue['theta']*(np.pi/180.))
dataBlueY=dataBlue['r']*np.sin(dataBlue['theta']*(np.pi/180.))
dataRedX=dataRed['r']*np.cos(dataRed['theta']*(np.pi/180.))
dataRedY=dataRed['r']*np.sin(dataRed['theta']*(np.pi/180.))

# for an ordered grid to separate galaxies
#nx=np.ceil(np.sqrt(data.size))
#datax=(np.tile(range(int(nx)),nx)*2/nx-1)[0:data.size]
#datay=(np.repeat(range(int(nx)),nx)*2/nx-1)[0:data.size]
    
# fudge x and y coords with these values since plotGalaxies expects SM and R
dataBlue['r']=dataBlueX
dataBlue['sm']=dataBlueY
dataRed['r']=dataRedX
dataRed['sm']=dataRedY
dataCen['sm']=0.


lw=2
circle1=plt.Circle((0.5,0.5),1./2,color="gray",alpha=0.09,lw=lw)
circle2=plt.Circle((0.5,0.5),2./3/2,color="gray",alpha=0.12,lw=lw)
circle3=plt.Circle((0.5,0.5),1./3/2,fill=True,color="gray",alpha=0.15,lw=lw)
ax1.add_artist(circle1)
ax1.add_artist(circle2)
ax1.add_artist(circle3)
plt.text(0.5,1.01,r"R$_{200c}$",color="gray",horizontalalignment="center",size="xx-small",alpha=0.5)

plotGalaxies(dataBlue,imDir,fullImSize,imSize,imSizePlot,ax2,-1.0,1.0,-1.0,1.0,minColor,maxColor,c0,c1,cmap,pivot,cmin,cmax,nShades,False)
plotGalaxies(dataRed,imDir,fullImSize,imSize,imSizePlot,ax2,-1.0,1.0,-1.0,1.0,minColor,maxColor,c0,c1,cmap,pivot,cmin,cmax,nShades,False)
#plotGalaxies(dataCen,imDir,fullImSize,imSize,imSizePlot,ax2,-1.0,1.0,-1.0,1.0,minColor,maxColor,c0,c1,cmap,pivot,cmin,cmax,nShades,False)
    
for oo in fig.findobj():
    oo.set_clip_on(False)

figSize=(10,10)
plt.savefig("../plots/stack_cartoon.pdf",bbox_inches="tight",dpi=300,figsize=figSize)

