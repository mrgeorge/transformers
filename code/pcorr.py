#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import matplotlib
matplotlib.use('Agg') # must appear before importing pyplot to get plots w/o GUI
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import plotgalSMvR
import matplotlib.gridspec as gridspec

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
         (data['rgba'] < 1) &
         (data['rgsize'] < 100)
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

def plotScatterMatrix(variables, varStr, varTicks, scatterPlotFile):

    # use helvetica and latex
    plt.clf()
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':8})
    plt.rc('text', usetex=True)
    plt.rc('axes',linewidth=1.5)

    nVars=len(variables)
    fig=plt.figure(1)
    gs=gridspec.GridSpec(nVars,nVars)
    gs.update(wspace=0,hspace=0)

    for vv in range(nVars):
         for uu in range(nVars):
              x=variables[vv]
              y=variables[uu]

              ax=fig.add_subplot(gs[(nVars-1-uu)*nVars+vv])
              ax.scatter(x,y,s=0.1)

              if(uu == 0):
                  plt.xlabel(varStr[vv])
                  ax.xaxis.set_ticks(varTicks[vv])
              else:
                  ax.set_xticklabels([])
              if(vv == 0):
                  plt.ylabel(varStr[uu])
                  ax.yaxis.set_ticks(varTicks[uu])
              else:
                  ax.set_yticklabels([])

    print 'saving scatter plot'
    plt.savefig(scatterPlotFile,bbox_inches='tight')
    print 'done saving scatter plot'

def pcorrMatrix(variables, method):

    nVars=len(variables)
    corrMatrix=np.zeros(nVars**2).reshape(nVars,nVars)
    pMatrix=np.zeros_like(corrMatrix)
    sigMatrix=np.zeros_like(corrMatrix)
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

    return (corrMatrix,pMatrix,sigMatrix)

def oplotText(fig, corrMatrix, sigMatrix):

    nVars=len(corrMatrix)
    x=np.indices((nVars,nVars))[0]
    y=x.transpose()
    yoff=0.1
    for ii in range(nVars**2):
        if(corrMatrix.flatten()[ii] > 0.999):
            color='white'
            sigStr=''
        else:
            color='black'
            sigStr='%4.1f' % np.abs(sigMatrix.flatten()[ii])
            
        fig.text(x.flatten()[ii],y.flatten()[ii]+yoff,'%5.2f' % corrMatrix.flatten()[ii],verticalalignment='bottom', horizontalalignment='center',size='xx-small',color=color)
        fig.text(x.flatten()[ii],y.flatten()[ii]-yoff,sigStr,verticalalignment='top', horizontalalignment='center',size='xx-small',color=color)
        

def plotCorrMatrix(variables, varStr, corrMatrix, pMatrix, sigMatrix, corrPlotFile, pPlotFile, sigPlotFile):

    # use helvetica and latex
    plt.clf()
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
    plt.rc('text', usetex=True)
    plt.rc('axes',linewidth=1.5)

    nVars=len(variables)
    my_cmap=matplotlib.colors.LinearSegmentedColormap.from_list('BlWBk',('blue','white','black'))
    
    plt.imshow(corrMatrix,interpolation='nearest',origin='lower',cmap=my_cmap,vmin=-1,vmax=1)
    oplotText(plt, corrMatrix, sigMatrix)
    xlocs, xlabels=plt.xticks(range(nVars),varStr)
    ylocs, ylabels=plt.yticks(range(nVars),varStr)
    plt.setp(xlabels,'rotation','vertical')
    plt.tick_params(which='both',length=0)
    plt.colorbar()
    
    plt.savefig(corrPlotFile,bbox_inches='tight')

def getVariables(varNames, data):
    
    indTheta=(varNames == 'theta')
    indSize=(varNames == 'rgsize')
    if(len(indTheta.nonzero()[0]) == 0):
        variables=[data[x] for x in varNames]
    else:
        theta=getOrientation(data['theta'],data['rgpa'])
        if(varNames[0] == 'theta'):
            variables=[theta]
        else:
            variables=[data[varNames[0]]]
        for ii in range(1,len(varNames)):
            if(varNames[ii] == 'theta'):
                variables.append(theta)
            else:
                variables.append(data[varNames[ii]])

    if(len(indSize.nonzero()[0] > 0)):
        variables[indSize.nonzero()[0]]=np.log10(variables[indSize.nonzero()[0]])

    return variables

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

    allVars={'z':       ('z', (0,0.5,1)),
             'sm':      (r'log(M$_{\star}$/M$_{\odot}$)', (9.5,10.5,11.5)),
             'mhalo':   (r'log(M$_{\rm 200c}$/M$_{\odot}$)', (13.5, 13.7, 13.9)),
             'r':       (r'R/R$_{\rm 200c}$', (0,0.5,1.)),
             'rgsize':  (r'log(R$_e$/kpc)', (0, 1, 2)),
             'rgsersic':(r'n$_S$', (0, 4, 8)),
             'color':   ('NUV-R', (0,3,6)),
             'ebv':     ('E(B-V)', (0,0.3,0.6)),
             'theta':   (r'$\theta$ (deg)', (0,45,90))
             }

    choiceVars=np.array(['z','sm','mhalo','r','rgsize','rgsersic','color'])
    variables=getVariables(choiceVars, data)
    varStr=[allVars[x][0] for x in choiceVars]
    varTicks=[allVars[x][1] for x in choiceVars]

    method='s' # spearman rank

    scatterPlotFile=plotDir+"scattermatrix.pdf"
    plotScatterMatrix(variables, varStr, varTicks, scatterPlotFile)

    corrPlotFile=plotDir+"corrmatrix.pdf"
    pPlotFile=plotDir+"pmatrix.pdf"
    sigPlotFile=plotDir+"sigmatrix.pdf"
    corrMatrix, pMatrix, sigMatrix = pcorrMatrix(variables, method)

    plotCorrMatrix(variables, varStr, corrMatrix, pMatrix, sigMatrix, corrPlotFile, pPlotFile, sigPlotFile)


# MAIN - if called from command line
if __name__ == '__main__':

    main()
    
