#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

from plotgalSMvR import *
import scipy.stats # for ranking data


imDir="../images/twenty_arcsec/"
imListFile=imDir+"acs_mem.list"

minZ=0.2
maxZ=0.5
zMed=np.median([minZ,maxZ])
minSM=9.8
maxSM=10.4
minMh=12.0
maxMh=15.0
colorSel="blue"
morph="all"
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
data=selectData(fullData, minZ, maxZ, minMh, maxMh, colorSel,morph, zType)
sel=((data['sm'] > minSM) & (data['sm'] < maxSM))
data=data[sel]


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

img=fitsio.read(imDir+data['filename'][0]) # read first one to get params
fullImSize=20. # fits file is 20" on a side
imSize=mrg.rp2deg(30.,cosmo.Da(0,zMed)*1000)*3600 # plot a square of this size on a side (30 kpc in arcsec)
imSizePlot=0.2 # fraction of axis range


# convert coordinates
    
# for a polar plot showing distance from group center and azimuthal angle
datax=data['r']*np.cos(data['theta']*(np.pi/180.))
datay=data['r']*np.sin(data['theta']*(np.pi/180.))

# for an ordered grid to separate galaxies
#nx=np.ceil(np.sqrt(data.size))
#datax=(np.tile(range(int(nx)),nx)*2/nx-1)[0:data.size]
#datay=(np.repeat(range(int(nx)),nx)*2/nx-1)[0:data.size]
    
# fudge x and y coords with these values since plotGalaxies expects SM and R
data['r']=datax
data['sm']=datay

plotGalaxies(data,imDir,fullImSize,imSize,imSizePlot,ax1,-1.0,1.0,-1.0,1.0,minColor,maxColor,c0,c1,cmap,pivot,cmin,cmax,nShades,False)
    
for oo in fig.findobj():
    oo.set_clip_on(False)

figSize=(10,10)
plt.savefig("../plots/art_blue.pdf",bbox_inches="tight",dpi=300,figsize=figSize)

