#! /usr/bin/env python

# IMPORT PYTHON MODULES
import sys
import string
import scipy.stats as stats
from os import path
import numpy as np
from itertools import combinations

# define class of gaussian_kde with user-defined covariance matrix
class gaussian_kde_set_covariance(stats.gaussian_kde):
	def __init__(self, dataset, covariance):
		self.covariance = covariance
		stats.gaussian_kde.__init__(self, dataset)
	def _compute_covariance(self):
		self.inv_cov = np.linalg.inv(self.covariance)
		self._norm_factor = np.sqrt(np.linalg.det(2*np.pi*self.covariance)) * self.n

def EstimateRel(data,pdfoutname,parNames,parSpace=['snr_sum','snr_max','n_pix'],logPars=[1,1,1],autoKernel=True,negPerBin=50,skellamTol=-0.2,kernel=[0.15,0.05,0.1],usecov=0,doscatter=1,docontour=1,doskellam=1,dostats=0,saverel=1,threshold=0.99,Nmin=0,dV=0.2,fMin=0,verb=0,makePlot=False):

	# matplotlib stuff
	if makePlot:
		import matplotlib
		# the following line is necessary to run SoFiA remotely
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt

	########################################
	### BUILD ARRAY OF SOURCE PARAMETERS ###
	########################################

	idCOL=parNames.index('id')
        ftotCOL=parNames.index('snr_sum')
        fmaxCOL=parNames.index('snr_max')
        fminCOL=parNames.index('snr_min')

        # get columns of requested parameters
	parCol = []
	for ii in range(0,len(parSpace)): parCol.append(parNames.index(parSpace[ii]))
	
	# get position and number of positive and negative sources
	pos=data[:,ftotCOL]>0
	neg=data[:,ftotCOL]<=0
	Npos=pos.sum()
	Nneg=neg.sum()
	if not Npos:
		sys.stderr.write("ERROR: no positive sources found; cannot proceed.\n")
		sys.exit(1)
	elif not Nneg:
		sys.stderr.write("ERROR: no negative sources found; cannot proceed.\n")
		sys.exit(1)

	# get array of relevant source parameters (and take log of them is requested)
	ids=data[:,idCOL]
        print '# Working in parameter space [',
	pars = np.empty((data.shape[0],0))
	for ii in range(0,len(parSpace)):
		print parSpace[ii],
		if parSpace[ii] == 'snr_max':
			parsTmp = data[:,fmaxCOL]*pos-data[:,fminCOL]*neg
			if logPars[ii]: parsTmp = np.log10(parsTmp)
			pars = np.concatenate((pars,parsTmp.reshape(-1,1)),axis=1)
		elif parSpace[ii] == 'snr_sum':
			parsTmp = abs(data[:,parCol[ii]].reshape(-1,1))
			if logPars[ii]: parsTmp = np.log10(parsTmp)
			pars = np.concatenate((pars,parsTmp),axis=1)
		else:
			parsTmp = data[:,parCol[ii]].reshape(-1,1)
			if logPars[ii]: parsTmp = np.log10(parsTmp)
			pars = np.concatenate((pars,parsTmp),axis=1)

        print ']'
	pars=np.transpose(pars)


	#################################################################
	### SET PARAMETERS TO WORK WITH AND GRIDDING/PLOTTNG FOR EACH ###
	#################################################################


	# axis labels when plotting
	labs = []
	for ii in range(0,len(parSpace)):
		labs.append('')
		if logPars[ii]: labs[ii] += 'log '
		labs[ii] += parSpace[ii]


	# axes limits when plotting
	pmin,pmax=pars.min(axis=1),pars.max(axis=1)
	pmin,pmax=pmin-0.1*(pmax-pmin),pmax+0.1*(pmax-pmin)
	lims = [[pmin[i],pmax[i]] for i in range(0,len(parSpace))]

	# grid on which to evaluate Np and Nn in order to plot contours
	grid = [[pmin[i],pmax[i],0.02*(pmax[i]-pmin[i])] for i in range(0,len(parSpace))]

	# calculate the number of rows and columns in figure
	projections = [subset for subset in combinations(range(0,len(parSpace)),2)]
	nr=int(np.floor(np.sqrt(len(projections))))
	nc=int(np.ceil(float(len(projections))/nr))


	###############################################
	### SET SMOOTHING KERNEL IN PARAMETER SPACE ###
	###############################################

	# If autoKernel is True the initial kernel is taken as a scaled version of the covariance matrix of the negative sources.
	# The kernel size along each axis is such that the number of sources per kernel width (sigma**2) is equal to 'negPerBin'.
        # Optionally, the user can decide to use only the diagonal terms of the covariance matrix.
        # The kernel is then grown until convergence is reached on the Skellam plot.
	# If autoKernel is False, use the kernel given by 'kernel' parameter (argument of EstimateRel); this is sigma, and is squared
        #    to be consistent with the auto kernel above.
	if autoKernel:
		#kernel=(pars[:,neg].max(axis=1)-pars[:,neg].min(axis=1))/(float(neg.sum())/negPerBin) # kernel from pars range
                parcov=np.sqrt(np.cov(pars[:,neg]))
                kernel=parcov/(float(Nneg)/negPerBin)     # kernel from covariance matrix
                if not usecov: kernel=np.diag(np.diag(kernel)) # kernel from diagonal of covariance matrix
	        grow_kernel=1 # set to 1 to start the kernel growing loop below; if autoKernel=0 will do just one pass
	        deltOLD=-1e+9 # used to stop kernel growth if P-N stops moving closer to zero [NOT USED CURRENTLY]
                deltplot=[]
                kernelIter=0.
	else: kernel,grow_kernel=np.identity(len(kernel))*np.array(kernel)**2,0

        if grow_kernel:
	        print '# Starting with (co-)variance kernel:'
                print kernel
                print '# Growing kernel...'
                sys.stdout.flush()
        else:
	        print '# Using user-defined variance kernel:'
                print kernel
                sys.stdout.flush()

	while grow_kernel:
                kernel*=float(negPerBin+kernelIter)/(negPerBin+kernelIter-1) 

		################################
		### EVALUATE N-d RELIABILITY ###
		################################

		if verb: print '  estimate normalised positive and negative density fields ...'
		
		# Calculate density fields Np and Nn using 'kernel' covariance
		#tmpCov = []
		#for ii in range(0,len(parSpace)):
		#	tmp = []
		#	for jj in range(0,ii):
		#		tmp.append(0)
		#	tmp.append(kernel[ii]**2)
		#	for jj in range(ii+1,len(parSpace)):
		#		tmp.append(0)
		#	tmpCov.append(tmp)
		
		#setcov = np.array(tuple(tmpCov))

                #print setcov
                #setcov=(np.sqrt(np.cov(pars))/(float(neg.sum())/negPerBin))**2
                #print setcov
                #exit()
		
		#if verb:
		#	print '  using diagonal kernel with sigma:',
		#	for kk in kernel: print '%.4f'%kk,
		#	print

		Np=gaussian_kde_set_covariance(pars[:,pos],kernel)
		Nn=gaussian_kde_set_covariance(pars[:,neg],kernel)

		#############################
		### PRINT STATS TO SCREEN ###
		#############################

		if docontour or dostats or doskellam:
			# volume within which to calculate the 
                        dV=(2*np.sqrt(np.linalg.eig(kernel)[0])).prod()

			# calculate the reliability at the location of positive sources
			if verb: print '  from the density fields, calculating the reliability at the location of positive sources ...'
			Nps=Np(pars[:,pos])*Npos
			Nns=Nn(pars[:,pos])*Nneg
			Rs=(Nps-Nns)/Nps

			# The reliability must be <=1. If not, something is wrong.
			if Rs.max()>1:
				sys.stderr.write("ERROR: maximum reliability larger than 1 -- something is wrong.\n")
				sys.exit(1)

			# calculate the number of positive and negative sources at the location of negative sources
			nNps=Np(pars[:,neg])*Npos
			nNns=Nn(pars[:,neg])*Nneg

			# I have verified that the integral NpI is equal to 0.85*Nps for a variety of
			# reasonable kernels, so I use 0.85*Nps as a proxy for NpI (same for Nns)
			if verb:
				# calculate the reliability at the location of negative sources
				nRs=(nNps-nNns)/nNps
				# calculate formal uncertainty of R at the location of positive and negative sources
				dRs=np.sqrt(Nps*Nns*(Nps+Nns))/Nps**2/np.sqrt(dV*0.85)
				dnRs=np.sqrt(nNps*nNns*(nNps+nNns))/nNps**2/np.sqrt(dV*0.85)
				print '  multiplying by 0.85*dV to get a proxy of integral of density fields within a +/-1 sigma(kernel) box ...'
				print '  source density at the location of positive sources (per +/-1 sigma(kernel) volume element):'
				print '             positive: %3.1f - %3.1f'%(0.85*dV*Nps.min(),0.85*dV*Nps.max())
				print '             negative: %3.1f - %3.1f'%(0.85*dV*Nns.min(),0.85*dV*Nns.max())
				print '  positive + negative: %3.1f - %3.1f'%(0.85*dV*(Nps+Nns).min(),0.85*dV*(Nps+Nns).max())

				print '  median error on R at location of negative sources: %.2f'%np.median(dnRs)
				print '  R<0 at the location of %3i/%3i negative sources'%((nRs<0).sum(),nRs.shape[0])
				print '  R<0 at the location of %3i/%3i negative sources within 1-sigma error bar'%(((nRs+1*dnRs)<0).sum(),nRs.shape[0])

				print '  median error on R at location of positive sources: %.2f'%np.median(dRs)
				print '  R<0 at the location of %3i/%3i positive sources'%((Rs<0).sum(),Rs.shape[0])
				print '  R<0 at the location of %3i/%3i positive sources within 1-sigma error bar'%(((Rs+1*dRs)<0).sum(),Rs.shape[0])


			# find pseudoreliable sources
			# (taking maximum(Rs,0) in order to include objects with Rs<0 if threshold==0)
			pseudoreliable=np.maximum(Rs,0)>=threshold
			# these are called pseudoreliable because some objets may be discarded later based on additional criteria below

			# find reliable sources
			# (taking maximum(Rs,0) in order to include objects with Rs<0 if threshold==0)
			# Nmin is by default zero so the line below normally selects (np.maximum(Rs,0)>=threshold)*(ftot[pos].reshape(-1,)>fMin)
			if 'snr_sum' in parSpace:
				ftot = pars[parSpace.index('snr_sum')]
			else:
				ftot = np.log10(abs(data[:,ftotCOL]).reshape(-1,1))
			reliable=(np.maximum(Rs,0)>=threshold)*((Nps+Nns)*0.85*dV>Nmin)*(ftot[pos].reshape(-1,)>fMin)
		
			# calculate quantities needed for comparison to Skellam distribution
			# multiplying by 0.85*dV to get a proxy of integral of density fields within a +/-1 sigma(kernel) box
			delt=np.sort((nNps-nNns)/np.sqrt(nNps+nNns)*np.sqrt(0.85*dV))

			if verb:
				print '  negative sources found:'
				print '    %20s: %4i'%('total',Nneg)
				print '  positive sources found:'
				print '    %20s: %4i'%('total',Npos),
				print
				print '                  R>%.2f: %4i'%(threshold,pseudoreliable.sum()),
				print
				print '     R>%.2f, N(3sig)>%3i: %4i'%(threshold,Nmin,reliable.sum()),
				print

		else: grow_kernel=0

                deltplot.append([negPerBin+kernelIter,delt[delt.shape[0]/2]])

		if delt[delt.shape[0]/2]>skellamTol or negPerBin+kernelIter==Nneg:
                        grow_kernel=0
	                print '# Found good kernel:'
                        print kernel
                        sys.stdout.flush()
                #elif: delt[delt.shape[0]/2]<deltOLD: grow_kernel=0
		else:
			deltOLD=delt[delt.shape[0]/2]
			#print deltOLD,'<',skellamTol; sys.stdout.flush()
			kernelIter+=1
        deltplot=np.array(deltplot)

	
	####################
	### SKELLAM PLOT ###
	####################

	if doskellam and makePlot:
		fig0=plt.figure()
		plt.plot(np.arange(-10,10,0.01),stats.norm().cdf(np.arange(-10,10,0.01)),'k-')
		plt.plot(np.arange(-10,10,0.01),stats.norm(scale=0.4).cdf(np.arange(-10,10,0.01)),'k:')
		plt.plot(delt,np.arange(1,delt.shape[0]+1,1,dtype=float)/delt.shape[0],'r-',drawstyle='steps-post')
		plt.xlim(-3,3)
		plt.ylim(0,1)
		plt.xlabel('(P-N)/sqrt(N+P)')
		plt.ylabel('cumulative distribution')
		plt.legend(('Gaussian (sigma=1)','Gaussian (sigma=0.4)','negative sources'),loc='upper left',prop={'size':15})
		plt.plot([0,0],[0,1],'k--')
		#plt.title('sigma(kde) = %.3f, %.3f, %.3f'%(kernel[0],kernel[1],kernel[2]),fontsize=20)
		fig0.savefig('%s_skel.pdf'%pdfoutname,rasterized=True)

                fig3=plt.figure()
                plt.plot(deltplot[:,0],deltplot[:,1],'ko-')
		plt.xlabel('negPerBin')
		plt.ylabel('delta')
                plt.axhline(y=skellamTol,linestyle='--',color='r')
		fig3.savefig('%s_delt.pdf'%pdfoutname,rasterized=True)


	############################
	### SCATTER PLOT SOURCES ###
	############################

	if doscatter and makePlot:
		if verb: print '  plotting sources ...'
		fig1=plt.figure(figsize=(18,4.5*nr))
		plt.subplots_adjust(left=0.06,bottom=0.15/nr,right=0.97,top=1-0.08/nr,wspace=0.35,hspace=0.25)

		n_p=0
		for jj in projections:
			if verb: print '    projection %i/%i'%(projections.index(jj)+1,len(projections))
			n_p,p1,p2=n_p+1,jj[0],jj[1]
			plt.subplot(nr,nc,n_p)
			plt.scatter(pars[p1,pos],pars[p2,pos],marker='o',c='b',s=10,edgecolor='face',alpha=0.5)
			plt.scatter(pars[p1,neg],pars[p2,neg],marker='o',c='r',s=10,edgecolor='face',alpha=0.5)
			plt.xlim(lims[p1][0],lims[p1][1])
			plt.ylim(lims[p2][0],lims[p2][1])
			plt.xlabel(labs[p1])
			plt.ylabel(labs[p2])
		fig1.savefig('%s_scat.pdf'%pdfoutname,rasterized=True)

	#####################
	### PLOT CONTOURS ###
	#####################

	if docontour and makePlot:
		levs=10**np.arange(-1.5,2,0.5)

		if verb: print '  plotting contours ...'
		fig2=plt.figure(figsize=(18,4.5*nr))
		plt.subplots_adjust(left=0.06,bottom=0.15/nr,right=0.97,top=1-0.08/nr,wspace=0.35,hspace=0.25)
		n_p=0
		for jj in projections:
			if verb: print '    projection %i/%i'%(projections.index(jj)+1,len(projections))
			n_p,p1,p2=n_p+1,jj[0],jj[1]
			g1,g2=grid[p1],grid[p2]
			x1=np.arange(g1[0],g1[1],g1[2])
			x2=np.arange(g2[0],g2[1],g2[2])
			pshape=(x2.shape[0],x1.shape[0])

			# get array of source parameters on current projection
			parsp=np.concatenate((pars[p1:p1+1],pars[p2:p2+1]),axis=0)

			# derive Np and Nn density fields on the current projection

			# Np and Nn calculated with *authomatic* covariance
			#Np=stats.kde.gaussian_kde(parsp[:,pos])
			#Nn=stats.kde.gaussian_kde(parsp[:,neg])

			# Np and Nn calculated with *input* covariance
                        setcov=kernel[p1:p2+1:p2-p1,p1:p2+1:p2-p1]
			Np=gaussian_kde_set_covariance(parsp[:,pos],setcov)
			Nn=gaussian_kde_set_covariance(parsp[:,neg],setcov)

			# evaluate density  fields on grid on current projection
			g=np.transpose(np.transpose(np.mgrid[slice(g1[0],g1[1],g1[2]),slice(g2[0],g2[1],g2[2])]).reshape(-1,2))
			Np=Np(g)
			Nn=Nn(g)
			Np=Np/Np.sum()*Npos
			Nn=Nn/Nn.sum()*Nneg
			Np.resize(pshape)
			Nn.resize(pshape)
			plt.subplot(nr,nc,n_p)
			plt.contour(x1,x2,Np,origin='lower',colors='b',levels=levs,zorder=2)
			plt.contour(x1,x2,Nn,origin='lower',colors='r',levels=levs,zorder=1)

			if reliable.sum(): plt.scatter(pars[p1,pos][reliable],pars[p2,pos][reliable],marker='o',s=10,edgecolor='k',facecolor='k',zorder=4)
			if (pseudoreliable*(reliable==False)).sum(): plt.scatter(pars[p1,pos][pseudoreliable*(reliable==False)],pars[p2,pos][pseudoreliable*(reliable==False)],marker='x',s=40,edgecolor='0.5',zorder=3)

			plt.xlim(lims[p1][0],lims[p1][1])
			plt.ylim(lims[p2][0],lims[p2][1])
			plt.xlabel(labs[p1])
			plt.ylabel(labs[p2])
		fig2.savefig('%s_cont.pdf'%pdfoutname,rasterized=True)

	###############################
	### ADD Np and Nn TO TABLES ###
	###############################

	# this allows me not to calculate R everytime I want to do
	# some plot analysis, but just read it from the file
	if saverel:
		if not (docontour or dostats):
			Nps=Np(pars[:,pos])*Npos
			Nns=Nn(pars[:,pos])*Nneg
		Np=np.zeros((data.shape[0],))
		Np[pos]=Nps
		Nn=np.zeros((data.shape[0],))
		Nn[pos]=Nns
		R=-np.ones((data.shape[0],)) # R will be -1 for negative sources
		# set R to zero for positive sources if R<0 because of Nn>Np
		R[pos]=np.maximum(0,(Np[pos]-Nn[pos])/Np[pos])
		data=np.concatenate((data,Np.reshape(-1,1),Nn.reshape(-1,1),R.reshape(-1,1)),axis=1)

	data=[list(jj) for jj in list(data)]
	return data,ids[pos][reliable].astype(int)
