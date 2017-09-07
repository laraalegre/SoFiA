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

def EstimateRel(data, pdfoutname, parNames, parSpace=['snr_sum', 'snr_max', 'n_pix'], logPars=[1, 1, 1], autoKernel=True, scaleKernel=1, negPerBin=1, skellamTol=-0.5, kernel=[0.15, 0.05, 0.1], usecov=False, doscatter=1, docontour=1, doskellam=1, dostats=0, saverel=1, threshold=0.99, fMin=0, verb=0, makePlot=False):
	# matplotlib stuff
	if makePlot:
		import matplotlib
		# the following line is necessary to run SoFiA remotely
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt
	
	########################################
	### BUILD ARRAY OF SOURCE PARAMETERS ###
	########################################
	
	idCOL   = parNames.index('id')
	ftotCOL = parNames.index('snr_sum')
	fmaxCOL = parNames.index('snr_max')
	fminCOL = parNames.index('snr_min')
	
	# get columns of requested parameters
	parCol = []
	for ii in range(len(parSpace)): parCol.append(parNames.index(parSpace[ii]))
	
	# get position and number of positive and negative sources
	pos  = data[:,ftotCOL] >  0
	neg  = data[:,ftotCOL] <= 0
	Npos = pos.sum()
	Nneg = neg.sum()
	
	if not Npos:
		sys.stderr.write("ERROR: no positive sources found; cannot proceed.\n")
		sys.exit(1)
	elif not Nneg:
		sys.stderr.write("ERROR: no negative sources found; cannot proceed.\n")
		sys.exit(1)
	
	# get array of relevant source parameters (and take log of them is requested)
	ids = data[:,idCOL]
	print '# Working in parameter space [',
	pars = np.empty((data.shape[0], 0))
	for ii in range(len(parSpace)):
		print parSpace[ii],
		if parSpace[ii] == 'snr_max':
			parsTmp = data[:,fmaxCOL] * pos - data[:,fminCOL] * neg
			if logPars[ii]: parsTmp = np.log10(parsTmp)
			pars = np.concatenate((pars, parsTmp.reshape(-1, 1)), axis=1)
		elif parSpace[ii] == 'snr_sum':
			parsTmp = abs(data[:,parCol[ii]].reshape(-1, 1))
			if logPars[ii]: parsTmp = np.log10(parsTmp)
			pars = np.concatenate((pars, parsTmp), axis=1)
		else:
			parsTmp = data[:,parCol[ii]].reshape(-1, 1)
			if logPars[ii]: parsTmp = np.log10(parsTmp)
			pars = np.concatenate((pars, parsTmp), axis=1)
	
	print ']'
	pars = np.transpose(pars)
	
	
	#################################################################
	### SET PARAMETERS TO WORK WITH AND GRIDDING/PLOTTNG FOR EACH ###
	#################################################################
	
	# axis labels when plotting
	labs = []
	for ii in range(len(parSpace)):
		labs.append('')
		if logPars[ii]: labs[ii] += 'log '
		labs[ii] += parSpace[ii]
	
	
	# axes limits when plotting
	pmin, pmax = pars.min(axis=1), pars.max(axis=1)
	pmin, pmax = pmin - 0.1 * (pmax - pmin), pmax + 0.1 * (pmax - pmin)
	lims = [[pmin[i], pmax[i]] for i in range(len(parSpace))]
	
	# grid on which to evaluate Np and Nn in order to plot contours
	grid = [[pmin[i], pmax[i], 0.02 * (pmax[i] - pmin[i])] for i in range(len(parSpace))]
	
	# calculate the number of rows and columns in figure
	projections = [subset for subset in combinations(range(len(parSpace)), 2)]
	nr = int(np.floor(np.sqrt(len(projections))))
	nc = int(np.ceil(float(len(projections)) / nr))
	
	
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
		# set the kernel shape to that of the variance or covariance matrix
		kernel = np.cov(pars[:,neg])
		kernelType = 'covariance'
		if np.isnan(kernel).sum():
			sys.stderr.write("ERROR: The reliability cannot be calculated because the smoothing kernel\n")
			sys.stderr.write("       derived from %i negative sources contains NaNs.\n"%pars[:,neg].shape[1])
			sys.stderr.write("       A good kernel is required to calculate the density field of positive\n")
			sys.stderr.write("       and negative sources in parameter space.\n")
			sys.stderr.write("       Try increase the number of negative sources by changing the source.\n")
			sys.stderr.write("       finding and/or filtering settings.\n")
			raise SystemExit(1)
		if not usecov:
			kernel = np.diag(np.diag(kernel))
			kernelType = 'variance'
		kernelIter = 0.0
		deltplot = []
		
		# scale the kernel size as requested by the user or by the auto-scale algorithm (scaleKernel>0)
		if scaleKernel:
			# scale kernel size as requested by the user
			# note that the scale factor is squared because users are asked to give a factor to apply to sqrt(kernel)  
			kernel *= scaleKernel**2
			print '# Using a kernel with the shape of the %s and size scaled by a factor %.2f.' % (kernelType, scaleKernel)
			print '# The sqrt(kernel) size is:'
			print np.sqrt(np.abs(kernel))
		else:
			# scale kernel size to start the kernel-growing loop
			# the scale factor for sqrt(kernel) is elevated to the power of 1./len(parCol)
			kernel *= ((negPerBin + kernelIter) / Nneg)**(2.0 / len(parCol))
			print '# Will find the best kernel as a scaled version of the %s:' % kernelType
			print '# Starting from the kernel with sqrt(kernel) size:'
			print np.sqrt(np.abs(kernel))
			print '# Growing kernel...'
			sys.stdout.flush()
		
		#deltOLD=-1e+9 # used to stop kernel growth if P-N stops moving closer to zero [NOT USED CURRENTLY]
		if doskellam and makePlot: fig0 = plt.figure()
	else:
		print '# Using user-defined variance kernel with sqrt(kernel) size:' # Note that the user must give sigma, which then gets squared
		print np.array(kernel)
		sys.stdout.flush()
		kernel = np.identity(len(kernel)) * np.array(kernel)**2
	
	grow_kernel = 1 # set to 1 to start the kernel growing loop below;
	                # this loop will estimate the reliability, check whether the kernel is large enough, and if not pick a larger kernel;
	                # if autoKernel=0 or scaleKernel=0 will do just one pass (i.e., will not grow the kernel)
	while grow_kernel:
		################################
		### EVALUATE N-d RELIABILITY ###
		################################
		
		if verb: print '#  estimate normalised positive and negative density fields ...'
		
		Np = gaussian_kde_set_covariance(pars[:,pos], kernel)
		Nn = gaussian_kde_set_covariance(pars[:,neg], kernel)
		
		# calculate the number of positive and negative sources at the location of positive sources
		Nps = Np(pars[:,pos]) * Npos
		Nns = Nn(pars[:,pos]) * Nneg
		
		# calculate the number of positive and negative sources at the location of negative sources
		nNps = Np(pars[:,neg]) * Npos
		nNns = Nn(pars[:,neg]) * Nneg
		
		# calculate the reliability at the location of positive sources
		Rs = (Nps - Nns) / Nps
		
		# The reliability must be <=1. If not, something is wrong.
		if Rs.max() > 1:
			sys.stderr.write("ERROR: maximum reliability larger than 1 -- something is wrong.\n")
			sys.exit(1)
		
		# find pseudoreliable sources
		# (taking maximum(Rs,0) in order to include objects with Rs<0 if threshold==0; Rs may be <0 becauase of insufficient statistics)
		pseudoreliable = np.maximum(Rs, 0) >= threshold
		# these are called pseudoreliable because some objets may be discarded later based on additional criteria below

		# find reliable sources
		# (taking maximum(Rs,0) in order to include objects with Rs<0 if threshold==0; Rs may be <0 becauase of insufficient statistics)
		#reliable=(np.maximum(Rs,0)>=threshold)*(data[pos,ftotCOL].reshape(-1,)>fMin)*(data[pos,fmaxCOL].reshape(-1,)>4)
		reliable = (np.maximum(Rs, 0) >= threshold) * (data[pos, ftotCOL].reshape(-1,) > fMin)
		
		if autoKernel:
			# calculate quantities needed for comparison to Skellam distribution
			delt = (nNps - nNns) / np.sqrt(nNps + nNns)
			deltstd = delt.std()
			deltmed = np.median(delt)
			deltmin = delt.min()
			deltmax = delt.max()
			
			if deltmed / deltstd > -100 and doskellam and makePlot:
				plt.hist(delt / deltstd, bins=np.arange(deltmin / deltstd, max(5.1, deltmax / deltstd), 0.01), cumulative=True, histtype='step', color=(min(1, float(negPerBin + kernelIter) / Nneg), 0,0), normed=True)
				deltplot.append([((negPerBin + kernelIter) / Nneg)**(1.0 / len(parCol)), deltmed / deltstd])
			
			print ' iteration, median, width, median/width = %3i, %9.2e, %9.2e, %9.2e' % (kernelIter, deltmed, deltstd, deltmed / deltstd)
			sys.stdout.flush()
			
			if scaleKernel: grow_kernel = 0
			elif deltmed / deltstd > skellamTol or negPerBin + kernelIter >= Nneg:
				grow_kernel = 0
				print '# Found good kernel after %i kernel growth iterations. The sqrt(kernel) size is:' % kernelIter
				print np.sqrt(np.abs(kernel))
				sys.stdout.flush()
			elif deltmed / deltstd < 5 * skellamTol:
				kernel *= (float(negPerBin + kernelIter + 20) / (negPerBin + kernelIter))**(2.0 / len(parCol)) 
				kernelIter += 20
			elif deltmed / deltstd < 2 * skellamTol:
				kernel *= (float(negPerBin + kernelIter + 10) / (negPerBin + kernelIter))**(2.0 / len(parCol))
				kernelIter += 10
			elif deltmed / deltstd < 1.5 * skellamTol:
				kernel *= (float(negPerBin + kernelIter + 3) / (negPerBin + kernelIter))**(2.0 / len(parCol))
				kernelIter += 3
			else:
				kernel *= (float(negPerBin + kernelIter + 1) / (negPerBin + kernelIter))**(2.0 / len(parCol))
				kernelIter += 1
		else:
			grow_kernel = 0
	
	
	####################
	### SKELLAM PLOT ###
	####################
	
	if autoKernel and doskellam and makePlot:
		plt.plot(np.arange(-10, 10, 0.01), stats.norm().cdf(np.arange(-10, 10, 0.01)), 'k-')
		plt.plot(np.arange(-10, 10, 0.01), stats.norm(scale=0.4).cdf(np.arange(-10, 10, 0.01)), 'k:')
		plt.legend(('Gaussian (sigma=1)', 'Gaussian (sigma=0.4)'), loc='lower right', prop={'size':13})
		plt.hist(delt / deltstd, bins=np.arange(deltmin / deltstd, max(5.1, deltmax / deltstd), 0.01), cumulative=True, histtype='step', color='r', normed=True)
		plt.xlim(-5, 5)
		plt.ylim(0, 1)
		plt.xlabel('(P-N)/sqrt(N+P)')
		plt.ylabel('cumulative distribution')
		plt.plot([0, 0], [0, 1], 'k--')
		fig0.savefig('%s_skel.pdf' % pdfoutname, rasterized=True)
		
		if not scaleKernel:
			fig3 = plt.figure()
			deltplot = np.array(deltplot)
			plt.plot(deltplot[:,0], deltplot[:,1], 'ko-')
			plt.xlabel('kernel size (1D-sigma, aribtrary units)')
			plt.ylabel('median/std of (P-N)/sqrt(P+N)')
			plt.axhline(y=skellamTol, linestyle='--', color='r')
			fig3.savefig('%s_delt.pdf' % pdfoutname, rasterized=True)
	
	
	############################
	### SCATTER PLOT SOURCES ###
	############################
	
	specialids = []
	
	if doscatter and makePlot:
		if verb: print '  plotting sources ...'
		fig1 = plt.figure(figsize=(18, 4.5 * nr))
		plt.subplots_adjust(left=0.06, bottom=0.15/nr, right = 0.97, top=1-0.08/nr, wspace=0.35, hspace=0.25)
		
		n_p = 0
		for jj in projections:
			if verb: print '    projection %i/%i' % (projections.index(jj) + 1, len(projections))
			n_p, p1, p2 = n_p + 1, jj[0], jj[1]
			plt.subplot(nr, nc, n_p)
			plt.scatter(pars[p1,pos], pars[p2,pos], marker='o', c='b', s=10, edgecolor='face', alpha=0.5)
			plt.scatter(pars[p1,neg], pars[p2,neg], marker='o', c='r', s=10, edgecolor='face', alpha=0.5)
			for si in specialids: plt.plot(pars[p1,ids==si], pars[p2,ids==si], 'kd', zorder=10000, ms=7, mfc='none', mew=2)
			plt.xlim(lims[p1][0], lims[p1][1])
			plt.ylim(lims[p2][0], lims[p2][1])
			plt.xlabel(labs[p1])
			plt.ylabel(labs[p2])
		fig1.savefig('%s_scat.pdf' % pdfoutname, rasterized=True)
	
	
	#####################
	### PLOT CONTOURS ###
	#####################
	
	if docontour and makePlot:
		levs = 10**np.arange(-1.5, 2, 0.5)
		
		if verb: print '  plotting contours ...'
		fig2 = plt.figure(figsize=(18, 4.5 * nr))
		plt.subplots_adjust(left=0.06, bottom=0.15/nr, right=0.97, top=1-0.08/nr, wspace=0.35, hspace=0.25)
		n_p = 0
		for jj in projections:
			if verb: print '    projection %i/%i' % (projections.index(jj) + 1, len(projections))
			n_p, p1, p2 = n_p + 1, jj[0], jj[1]
			g1, g2 = grid[p1], grid[p2]
			x1 = np.arange(g1[0], g1[1], g1[2])
			x2 = np.arange(g2[0], g2[1], g2[2])
			pshape = (x2.shape[0], x1.shape[0])
			
			# get array of source parameters on current projection
			parsp = np.concatenate((pars[p1:p1+1], pars[p2:p2+1]), axis=0)
			
			# derive Np and Nn density fields on the current projection
			setcov = kernel[p1:p2+1:p2-p1,p1:p2+1:p2-p1]
			Np = gaussian_kde_set_covariance(parsp[:,pos], setcov)
			Nn = gaussian_kde_set_covariance(parsp[:,neg], setcov)
			
			# evaluate density  fields on grid on current projection
			g = np.transpose(np.transpose(np.mgrid[slice(g1[0], g1[1], g1[2]), slice(g2[0], g2[1], g2[2])]).reshape(-1, 2))
			Np = Np(g)
			Nn = Nn(g)
			Np = Np / Np.sum() * Npos
			Nn = Nn / Nn.sum() * Nneg
			Np.resize(pshape)
			Nn.resize(pshape)
			plt.subplot(nr, nc, n_p)
			plt.contour(x1, x2, Np, origin='lower', colors='b', levels=levs, zorder=2)
			plt.contour(x1, x2, Nn, origin='lower', colors='r', levels=levs, zorder=1)
			
			if reliable.sum(): plt.scatter(pars[p1,pos][reliable], pars[p2,pos][reliable], marker='o', s=10, edgecolor='k', facecolor='k', zorder=4)
			if (pseudoreliable * (reliable == False)).sum(): plt.scatter(pars[p1,pos][pseudoreliable * (reliable == False)], pars[p2,pos][pseudoreliable * (reliable == False)], marker='x', s=40, edgecolor='0.5', facecolor='0.5', zorder=3)
			for si in specialids: plt.plot(pars[p1,ids==si], pars[p2,ids==si], 'kd', zorder=10000, ms=7, mfc='none', mew=2)
			plt.xlim(lims[p1][0], lims[p1][1])
			plt.ylim(lims[p2][0], lims[p2][1])
			plt.xlabel(labs[p1])
			plt.ylabel(labs[p2])
		fig2.savefig('%s_cont.pdf' % pdfoutname, rasterized=True)
	
	
	#################################
	### ADD Np, Nn AND R TO TABLE ###
	#################################
	
	# this allows me not to calculate R everytime I want to do
	# some plot analysis, but just read it from the file
	if saverel:
		if not (docontour or dostats):
			Nps = Np(pars[:,pos]) * Npos
			Nns = Nn(pars[:,pos]) * Nneg
		Np = np.zeros((data.shape[0],))
		Np[pos] = Nps
		Nn = np.zeros((data.shape[0],))
		Nn[pos] = Nns
		R = -np.ones((data.shape[0],)) # R will be -1 for negative sources
		# set R to zero for positive sources if R<0 because of Nn>Np
		R[pos] = np.maximum(0, (Np[pos] - Nn[pos]) / Np[pos])
		data = np.concatenate((data, Np.reshape(-1, 1), Nn.reshape(-1, 1), R.reshape(-1, 1)), axis=1)
	
	data = [list(jj) for jj in list(data)]
	return data, ids[pos][reliable].astype(int)
