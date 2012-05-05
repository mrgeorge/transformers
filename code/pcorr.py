#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import plotgalSMvR

def pcorr(x,y,z,method):
    # Compute partial correlation coefficient between x and y given z
    # x and y are n-element arrays
    # z is an m x n element array
    # method is 'p' for Pearson or 's' for Spearman rank
    # This is a port of the "var-covar" method pcor.mat written in R at
    # http://www.yilab.gatech.edu/pcor.R
    # See http://www.yilab.gatech.edu/pcor.html for more info

    # Note: I've followed some of their weird naming conventions
    #  e.g., Sxx is the covariance between x and y
    # Also note that cov returns different things in R and in numpy

    nData=len(x)
    nControl=(np.shape(z))[0] # number of control variables

    # The only difference between Pearson and Spearman is that
    # we rank the variables first for Spearman. In the large-N
    # limit they should also have the same null distribution (for pvalue)
    if(method=='p'): # pearson
          xc=x
          yc=y
          zc=z
    elif(method=='s'): # spearman rank
          xc=stats.rankdata(x)
          yc=stats.rankdata(y)
          zc=np.zeros_like(z)
          for ii in range(nControl):
               zc[ii]=stats.rankdata(z[ii])
    else:
         print "Error in pcorr: must specify method"
         exit

    Sxx=np.cov(xc,yc)

    Sxz=np.zeros(2*nControl).reshape(2,nControl)
    Szz=np.zeros(nControl**2).reshape(nControl,nControl)
    for ii in range(nControl):
         Sxz[0,ii]=(np.cov(xc,zc[ii]))[0,1]
         Sxz[1,ii]=(np.cov(yc,zc[ii]))[0,1]
         for jj in range(nControl):
              Szz[ii,jj]=(np.cov(zc[ii],zc[jj]))[0,1]
              
    # Check that Szz is positive definite before inverting
    if(np.min(stats.stats.linalg.eigvals(Szz)) < 0):
         print "Error in pcorr: Szz is not positive definite"
         exit
         
    SzzInv=np.linalg.inv(Szz)
    Sxxz=Sxx - np.dot(np.dot(Sxz,SzzInv),Sxz.T)

    coeff=Sxxz[0,1]/(np.sqrt(Sxxz[0,0])*np.sqrt(Sxxz[1,1]))

    pvalue, significance=pcorr_pvalue(coeff, nData, nControl)

    return (coeff, pvalue, significance)

def pcorr_pvalue(coeff, n, gn):
    # Calculate 2-tailed P-value for partial correlation coefficient
    # (this is just statistics-speak for saying "calculate the area under a gaussian with mean=0, std=1, in both tails beyond a threshold")
    # n = number of data points
    # gn = number of "given" variables that are conditioned on
    # See http://www.yilab.gatech.edu/pcor.html for more info
     
    statistic=coeff * np.sqrt((n-2.-gn)/(1.-coeff**2))
    pvalue=2.*stats.norm.cdf(-np.abs(statistic))

    return (pvalue, statistic)

def selectData(data, minZ, maxZ, minMh, maxMh):

    sel=((data['flag'] == 1) &    # only include flag=1
         (data['z'] > minZ) &    # make a tighter z cut
         (data['z'] < maxZ) &
         (data['mhalo'] > minMh) &
         (data['mhalo'] < maxMh) &
         (data['nomatch'] == 0) &
         (data['rgflag'] == 0) &
         (data['rgba'] < 1)
         )
    data=data[sel]

    return data

def getOrientation(theta, rgpa):
    # get angle (0 to 90) between major axis of disk and direction to central
    # theta is clockwise angle between central, satellite, and east (-180 to +180)
    # rgpa is counter-clockwise angle between north and major axis of satellite (-90 to +90)
    # returns orientation, which is angle between major axis and direction to central (0 to 90)

    # there's probably a better way to do this but it's easy to picture the quadrant cases with satellite at center, N up, E left
    sel3n=((theta >= 0) & 
           (theta < 90) &
           (rgpa < 0))
    sel3p=((theta >= 0) & 
           (theta < 90) &
           (rgpa >= 0))
    sel4n=((theta >= 90) & 
           (theta < 180) &
           (rgpa < 0))
    sel4p=((theta >= 90) & 
           (theta < 180) &
           (rgpa >= 0))
    sel2n=((theta >= -90) & 
           (theta < 0) &
           (rgpa < 0))
    sel2p=((theta >= -90) & 
           (theta < 0) &
           (rgpa >= 0))
    sel1n=((theta >= -180) & 
           (theta < -90) &
           (rgpa < 0))
    sel1p=((theta >= -180) & 
           (theta < -90) &
           (rgpa >= 0))

    orientation=np.zeros_like(theta)
    orientation[sel1n]=(np.abs(-(90. + theta) + rgpa))[sel1n]
    orientation[sel1p]=(-(90. + theta) + rgpa)[sel1p]
    orientation[sel2n]=((90. + theta) - rgpa)[sel2n]
    orientation[sel2p]=(np.abs((90. - rgpa) + theta))[sel2p]
    orientation[sel3n]=(np.abs(theta - (90. + rgpa)))[sel3n]
    orientation[sel3p]=(theta + (90. - rgpa))[sel3p]
    orientation[sel4n]=((theta - 90.) - rgpa)[sel4n]
    orientation[sel4p]=(np.abs((theta - 90.) - rgpa))[sel4p]

    gt90=(orientation > 90.)
    orientation[gt90]=180. - orientation[gt90]

    bad=(rgpa < -90) # flag value is -99
    orientation[bad]=-99

    return orientation


def main():
    imDir="../images/"
    imListFile=imDir+"acs_mem.list"
    plotDir="../plots/"

    minZ=0.1
    maxZ=1.0
    minMh=13.5
    maxMh=14.0

    data=plotgalSMvR.readData(imListFile)
    data=selectData(data, minZ, maxZ, minMh, maxMh)

    varStr=['z','mhalo','sm','r','rgsize','rgba','rgsersic','color','ebv']
    variables=[data[x] for x in varStr]

    #    variables=(data['z'],
               # data['mhalo'],
               # data['sm'],
               #           data['r'],
               # data['rgsize'],
               #data['rgba'],
               # data['rgsersic'],
               # data['color'],
               # data['ebv'],
               #getOrientation(data['theta'],data['rgpa'])
               #)

    nVars=len(variables)
    corrMatrix=np.zeros(nVars**2).reshape(nVars,nVars)
    pMatrix=np.zeros_like(corrMatrix)
    sigMatrix=np.zeros_like(corrMatrix)
    method='s' # spearman rank
    for vv in range(nVars):
         for uu in range(nVars):
              controls=~((np.arange(nVars)==vv) | (np.arange(nVars)==uu)) # all vars except uu or vv
              x=variables[vv]
              y=variables[uu]
              z=[variables[ww] for ww in controls.nonzero()[0]]
              pcoeff, pval, sig=pcorr(x,y,z,method)
              corrMatrix[vv,uu]=pcoeff
              pMatrix[vv,uu]=pval
              sigMatrix[vv,uu]=sig

    plt.imshow(corrMatrix,interpolation='nearest',cmap=plt.cm.binary,origin='lower')
    xlocs, xlabels=plt.xticks(range(nVars),varStr)
    ylocs, ylabels=plt.yticks(range(nVars),varStr)
    plt.setp(xlabels,'rotation','vertical')
    plt.colorbar()
    plt.show()

    print corrMatrix
    print pMatrix
    print sigMatrix

    print "end"

# MAIN - if called from command line
if __name__ == '__main__':

    main()
    
