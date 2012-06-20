#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import matplotlib
matplotlib.use('Agg') # must appear before importing pyplot to get plots w/o GUI
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import plotgalSMvR
import matplotlib.gridspec as gridspec
import re # to remove decimal from string

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
     
    if(np.abs(coeff) > 0.9999):
          statistic=999. * coeff
          pvalue=0.
    else:
          statistic=coeff * np.sqrt((n-2.-gn)/(1.-coeff**2))
          pvalue=2.*stats.norm.cdf(-np.abs(statistic))

    return (pvalue, statistic)

def setupBins(limOption):

    if(limOption == 'all'):
          minZ=[0.1]
          maxZ=[1.0]
          minMh=[12.0]
          maxMh=[15.0]
          minSM=[8.0]
          maxSM=[13.0]
          zTicks=[(0., 0.5, 1.0)]
          smTicks=[(9.5, 10.5, 11.5)]
          mhTicks=[(13.0, 13.5, 14.0)]
    elif(limOption == 'complete'):
          minZ=[0.1]
          maxZ=[1.0]
          minMh=[13.5]
          maxMh=[15.0]
          minSM=[10.3]
          maxSM=[13.0]
          zTicks=[(0., 0.5, 1.0)]
          smTicks=[(10.4, 11.0, 11.6)]
          mhTicks=[(13.6, 13.8, 14.0)]
    elif(limOption == 'thresh'):
          minZ=[0.2,0.5,0.8]
          maxZ=[0.5,0.8,1.0]
          minMh=[13.3,13.5,13.6]
          maxMh=[15.0,15.0,15.0]
          minSM=[9.4,9.8,10.3]
          maxSM=[13.0,13.0,13.0]
          zTicks=[(0.2, 0.3, 0.5), (0.5, 0.6, 0.7), (0.8, 0.9, 1.0)]
          smTicks=[(9.5, 10.5, 11.5),(10, 10.8, 11.6),(10.4, 11.0, 11.6)]
          mhTicks=[(13.3, 13.6, 13.9),(13.5, 13.7, 13.9),(13.6, 13.8, 14.0)]
    elif(limOption == 'bins'):
          minZ=np.repeat([0.2,0.5,0.8],2)
          maxZ=np.repeat([0.5,0.8,1.0],2)
          minMh=np.repeat([13.3,13.5,13.6],2)
          maxMh=np.repeat([14.0,14.0,14.0],2)
          minSM=[9.4,10.5,9.8,10.5,10.3,10.7]
          maxSM=[10.5,13.0,10.5,13.0,10.7,13.0]
          zTicks=[(0.2, 0.3, 0.5), (0.2, 0.3, 0.5),
                  (0.5, 0.6, 0.7), (0.5, 0.6, 0.7),
                  (0.8, 0.9, 1.0), (0.8, 0.9, 1.0)]
          smTicks=[(9.5, 10.5, 11.5),(9.5, 10.5, 11.5),
                   (10, 10.8, 11.6),(10, 10.8, 11.6),
                   (10.4, 11.0, 11.6),(10.4, 11.0, 11.6)]
          mhTicks=[(13.3, 13.6, 13.9),(13.3, 13.6, 13.9),
                   (13.5, 13.7, 13.9),(13.5, 13.7, 13.9),
                   (13.6, 13.8, 14.0),(13.6, 13.8, 14.0)]

    return (minZ, maxZ, minMh, maxMh, minSM, maxSM, zTicks, mhTicks, smTicks)

def selectData(data, csOption='censat',
               minZ=0.1, maxZ=1.0,
               minMh=12.0, maxMh=15.0,
               minSM=8.0, maxSM=13.0):

    print 'Ngal before cuts:', data.size
      
    sel=((data['flag'] == 1) &    # only include flag=1
         (data['z'] > minZ) &    # make a tighter z cut
         (data['z'] < maxZ) &
         (data['mhalo'] > minMh) &
         (data['mhalo'] < maxMh) &
         (data['sm'] > minSM) &
         (data['sm'] < maxSM) &
         (data['nomatch'] == 0) &
         (data['rgflag'] == 0) &
         (data['rgba'] < 1) & # there are a few to toss with b/a > 1 ?
         (data['rgsize'] < 100) # there are a few with extremely high values to toss ?
         )
    data=data[sel]

    if(csOption == 'cen'):
          sel=(data['r'] == 0.)
          data=data[sel]
    elif(csOption == 'sat'):
          sel=(data['r'] > 0.)
          data=data[sel]

    print 'Ngal after cuts:', data.size
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

def plotScatterMatrix(variables, varStr, varTicks, corrMatrix, sigMatrix, scatterPlotFile, csOption, minZ, maxZ, minMh, maxMh, minSM, maxSM):
    # make multi-panel plot including scatter plot of each pair of variables
    # with the background colors to show the partial correlation matrix

    my_cmap=getCorrCMap()
    # kludge to get mappable object for colorbar
    img=plt.imshow(corrMatrix,interpolation='nearest',origin='lower',cmap=my_cmap,vmin=-1,vmax=1)
    plt.clf()

    # use helvetica and latex
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':10})
    plt.rc('text', usetex=True)
    plt.rc('axes',linewidth=1.)

    nVars=len(variables)
    fig=plt.figure(1)
    gs=gridspec.GridSpec(nVars,nVars)
    gs.update(wspace=0,hspace=0,right=0.93)

    gscolor=gridspec.GridSpec(1,1)
    gscolor.update(left=0.95,right=0.98)

    # plot range and size
    xstart=0
    ystart=0
    xend=0.93
    yend=1.
    xlen=(xend-xstart)/nVars
    ylen=(xend-xstart)/nVars

    fig.set_size_inches(7,7)

    # add title text
    matrixLegend('scat',fig, variables[0].size, csOption, minZ, maxZ, minMh, maxMh, minSM, maxSM)
    
    for xx in range(nVars):
         for yy in range(nVars):
              x=variables[xx]
              y=variables[yy]

              # create subplot with background color = correlation coefficient
              bgcolor=matplotlib.colors.colorConverter.to_rgb(my_cmap((corrMatrix[xx,yy]+1.)/2.))
              ax=fig.add_subplot(gs[(nVars-1-yy)*nVars+xx], axisbg=bgcolor)

              if(xx == yy): # plot histogram along diagonal
                    # dummy plot to get axis range
                    ax.scatter(x,y)
                    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=3,prune='lower'))
                    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                    xlim=ax.get_xlim()
                    #xmaj=ax.xaxis.get_major_locator()
                    #xmin=ax.xaxis.get_minor_locator()
                    xmaj=ax.xaxis.get_ticklocs()
                    
                    # clear scatter plot and make histogram
                    ax.cla()
                    n, bins, patches=ax.hist(x, bins=np.sqrt(len(x)), range=xlim, color='white', histtype='step', normed=1)

                    # set x axis to look like scatter plot
                    ax.set_xticks(xmaj)
                    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                    ax.set_xlim(xlim)
                    
                    # set y axis to look like x (it's arbitrarily normalized for the histogram)
                    ymin=0.
                    ymax=1.2*np.max(n)
                    ylim=(ymin,ymax)
                    yscale=(ymax-ymin)/(xlim[1]-xlim[0])
                    
                    ax.set_yticks(ymin + (xmaj-xlim[0])*yscale)
                    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                    ax.set_yticklabels(xmaj)
                    ax.set_ylim(ylim)

              else: # scatter plot on off-diagonals
                    ax.scatter(x,y,s=1,color='black',edgecolors='none')

                    # setup axes
                    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3,prune='lower'))
                    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3,prune='lower'))
                    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())

              # ax.xaxis.set_ticks(varTicks[xx])
              # ax.yaxis.set_ticks(varTicks[yy])

              # add axis labels along edges, else don't label
              if(yy == 0):
                  plt.xlabel(varStr[xx])
              else:
                  ax.set_xticklabels([])
              if(xx == 0):
                  plt.ylabel(varStr[yy])
              else:
                  ax.set_yticklabels([])

              # add correlation and significance to each panel
              oplotScatterText(ax,corrMatrix[xx,yy],sigMatrix[xx,yy])

    # add color bar
    gs.update(hspace=0)
    axcolor=fig.add_subplot(gscolor[0,0])
    cbar=plt.colorbar(img, cax=axcolor)
    cbar.ax.set_ylabel(r'$\rho_{xy|\{V\}}$')

    print 'saving scatter plot'
    plt.savefig(scatterPlotFile,bbox_inches='tight',pad_inches=0.25)
    print 'done saving scatter plot'

def oplotScatterText(ax, corr, sig):
    # add correlation coefficient and significance to a panel in the scatter plot

    # where to put the text
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    x=0.06 * (x1-x0) + x0
    y=0.94 * (y1-y0) + y0

    if(corr > 0.9999):
        color='white'
        sigStr=''
    else:
        color='black'
        sigStr=', '+'%4.1f' % np.abs(sig)
            
    corrStr='%5.2f' % corr
    ax.text(x,y,corrStr+sigStr,verticalalignment='top', horizontalalignment='left',size='x-small',color=color)
        
def pcorrMatrix(variables, controlVars, method):
    # compute the partial correlation matrix for a set of variables
    # variables is list (nVar long) of arrays (each nGal long)
    # controlVars is an ndarray giving the indices of which variables should be held fixed (when not being correlated)

    nVars=len(variables)
    corrMatrix=np.zeros(nVars**2).reshape(nVars,nVars)
    pMatrix=np.zeros_like(corrMatrix)
    sigMatrix=np.zeros_like(corrMatrix)
    for xx in range(nVars):
        for yy in range(nVars):
            controlInd= ~((controlVars==xx) | (controlVars==yy)) # all of controlVars except yy and xx
            x=variables[xx]
            y=variables[yy]
            z=[variables[ww] for ww in controlInd.nonzero()[0]]
            pcoeff, pval, sig=pcorr(x,y,z,method)
            corrMatrix[xx,yy]=pcoeff
            pMatrix[xx,yy]=pval
            sigMatrix[xx,yy]=sig

    return (corrMatrix,pMatrix,sigMatrix)

def oplotText(fig, corrMatrix, sigMatrix):
    # print correlation and significance in each panel of correlation matrix plot

    nVars=len(corrMatrix)
    x=np.indices((nVars,nVars))[0]
    y=x.transpose()
    yoff=0.1
    for ii in range(nVars**2):
        if(corrMatrix.flatten()[ii] > 0.9999):
            color='white'
            sigStr=''
        else:
            color='black'
            sigStr='%4.1f' % np.abs(sigMatrix.flatten()[ii])
            
        fig.text(x.flatten()[ii],y.flatten()[ii]+yoff,'%5.2f' % corrMatrix.flatten()[ii],verticalalignment='bottom', horizontalalignment='center',size='xx-small',color=color)
        fig.text(x.flatten()[ii],y.flatten()[ii]-yoff,sigStr,verticalalignment='top', horizontalalignment='center',size='xx-small',color=color)
        
def matrixLegend(type, fig, nGal, csOption, minZ, maxZ, minMh, maxMh, minSM, maxSM):
    # print plot info as title or legend (e.g., ranges for z, SM, Mh; Ngal; Cen vs Sat)

    if(csOption == 'cen'):
          csStr='Centrals'
    elif(csOption == 'sat'):
          csStr='Satellites'
    if(csOption == 'censat'):
          csStr='Centrals + Satellites'
    nStr='N = '+str(nGal)
    zStr=str(minZ)+'$<$z$<$ '+str(maxZ)
    if(maxSM < 13.):
          smStr=str(minSM)+r'$<$log(M$_{\star}$/M$_{\odot}$)$<$'+str(maxSM)
    else:
          smStr=r'log(M$_{\star}$/M$_{\odot}$)$>$'+str(minSM)
    if(maxMh < 15.):
          mhStr=str(minMh)+r'$<$log(M$_{\rm 200c}$/M$_{\odot}$)$<$'+str(maxMh)
    else:
          mhStr=r'log(M$_{\rm 200c}$/M$_{\odot}$)$>$'+str(minMh)


    if(type == 'corr'):
          fullStr=csStr+'\n'+nStr+'\n'+zStr+'\n'+smStr+'\n'+mhStr
          fig.figtext(-0.05,-0.17,fullStr,size='x-small',ha='left',va='bottom')
    elif(type == 'scat'):
          #          fig.figtext(0.95,0.05,fullStr,size='small',ha='left',va='bottom')
          fullStr=csStr+', '+nStr+', '+zStr+', '+smStr+', '+mhStr
          fig.suptitle(fullStr, x=0.5, y=0.93, fontsize=12)

def getCorrCMap():
    # define blue-white-black color map for correlation matrix
    my_cmap=matplotlib.colors.LinearSegmentedColormap.from_list('BlWBk',('blue','white','black'))
    return my_cmap

def plotCorrMatrix(variables, varStr, corrMatrix, sigMatrix, corrPlotFile, csOption, minZ, maxZ, minMh, maxMh, minSM, maxSM):
    # plot partial correlation matrix with colored boxes
      
    # use helvetica and latex
    plt.clf()
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
    plt.rc('text', usetex=True)
    plt.rc('axes',linewidth=1.2)

    nVars=len(variables)
    my_cmap=getCorrCMap()
    
    plt.imshow(corrMatrix,interpolation='nearest',origin='lower',cmap=my_cmap,vmin=-1,vmax=1)
    oplotText(plt, corrMatrix, sigMatrix)
    xlocs, xlabels=plt.xticks(range(nVars),varStr)
    ylocs, ylabels=plt.yticks(range(nVars),varStr)
    plt.setp(xlabels,'rotation','vertical')
    plt.tick_params(which='both',length=0)
    plt.colorbar()

    matrixLegend('core',plt, variables[0].size, csOption, minZ, maxZ, minMh, maxMh, minSM, maxSM)

    plt.savefig(corrPlotFile,bbox_inches='tight')

def getVariables(varNames, data):
    # convert input data to variable vectors
    # most of the work here is to compute theta if needed and put it in order
    
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

def main(csOption='censat', limOption='all'):
    imDir="../images/"
    imListFile=imDir+"acs_mem.list"
    plotDir="../plots/"

    method='s' # spearman rank

    minZ, maxZ, minMh, maxMh, minSM, maxSM, zTicks, mhTicks, smTicks = setupBins(limOption)
    nBins=len(minZ)
    if(csOption == 'cen'):
        choiceVars=np.array(['z','sm','mhalo','rgsize','rgsersic','color'])
        controlVars=np.array([0,1,2]) # indices of choiceVars that will be held fixed
    else:
        choiceVars=np.array(['z','sm','mhalo','r','rgsize','rgsersic','color'])
        controlVars=np.array([0,1,2,3]) # indices of choiceVars that will be held fixed

    dataAll=plotgalSMvR.readData(imListFile)

    for ii in range(nBins):
        data=selectData(dataAll, csOption,
                        minZ=minZ[ii], maxZ=maxZ[ii],
                        minMh=minMh[ii], maxMh=maxMh[ii],
                        minSM=minSM[ii], maxSM=maxSM[ii]
                       )
        if(data.size < 10):
              print 'Not enough galaxies, skipping'
              continue

        allVars={'z':       ('z', zTicks[ii]),
                 'sm':      (r'log(M$_{\star}$/M$_{\odot}$)', smTicks[ii]),
                 'mhalo':   (r'log(M$_{\rm 200c}$/M$_{\odot}$)', mhTicks[ii]),
                 'r':       (r'R/R$_{\rm 200c}$', (0,0.5,1.)),
                 'rgsize':  (r'log(R$_e$/kpc)', (0, 1, 2)),
                 'rgsersic':(r'n$_S$', (0, 4, 8)),
                 'color':   ('NUV-R', (0,3,6)),
                 'ebv':     ('E(B-V)', (0,0.3,0.6)),
                 'theta':   (r'$\theta$ (deg)', (0,45,90))
                 }

        variables=getVariables(choiceVars, data)
        varStr=[allVars[x][0] for x in choiceVars]
        varTicks=[allVars[x][1] for x in choiceVars]

        binStr="_z"+re.sub("\D","",str(minZ[ii]))+"_"+re.sub("\D","",str(maxZ[ii]))+"_mh"+re.sub("\D","",str(minMh[ii]))+"_"+re.sub("\D","",str(maxMh[ii]))+"_sm"+re.sub("\D","",str(minSM[ii]))+"_"+re.sub("\D","",str(maxSM[ii]))
        scatterPlotFile=plotDir+"scat_"+csOption+binStr+".pdf"
        corrPlotFile=plotDir+"corr_"+csOption+binStr+".pdf"

        corrMatrix, pMatrix, sigMatrix = pcorrMatrix(variables, controlVars, method)

        plotScatterMatrix(variables, varStr, varTicks, corrMatrix, sigMatrix, scatterPlotFile, csOption, minZ[ii], maxZ[ii], minMh[ii], maxMh[ii], minSM[ii], maxSM[ii])

        # plotCorrMatrix(variables, varStr, corrMatrix, sigMatrix, corrPlotFile, csOption, minZ[ii], maxZ[ii], minMh[ii], maxMh[ii], minSM[ii], maxSM[ii])


# MAIN - if called from command line
if __name__ == '__main__':
    import sys

    if(len(sys.argv) == 1):
          main()
    elif(len(sys.argv) == 2):
          main(csOption=sys.argv[1])
    elif(len(sys.argv) == 3):
          main(csOption=sys.argv[1], limOption=sys.argv[2])
    else:
          print "Calling sequence: pcorr.py csOption limOption"
    
