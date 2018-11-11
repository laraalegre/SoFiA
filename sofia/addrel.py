#! /usr/bin/env python

import numpy as np
import scipy.stats as stats
from itertools import combinations
from sofia import error as err

# =======================================================================
# CLASS: Define class of gaussian_kde with user-defined covariance matrix
# =======================================================================

class gaussian_kde_set_covariance(stats.gaussian_kde):
	def __init__(self, dataset, covariance):
		self.covariance = covariance
		stats.gaussian_kde.__init__(self, dataset)
	def _compute_covariance(self):
		self.inv_cov = np.linalg.inv(self.covariance)
		self._norm_factor = np.sqrt(np.linalg.det(2 * np.pi * self.covariance)) * self.n


# ================================================
# FUNCTION: Estimate the reliability of detections
# ================================================

def EstimateRel(data, pdfoutname, parNames, parSpace=["snr_sum", "snr_max", "n_pix"], logPars=[1, 1, 1], autoKernel=True, scaleKernel=1, negPerBin=1, skellamTol=-0.5, kernel=[0.15, 0.05, 0.1], usecov=False, doscatter=1, docontour=1, doskellam=1, dostats=0, saverel=1, threshold=0.99, fMin=0, verb=0, makePlot=False):

	# Set negPerBin to be >=1
	negPerBin=max(1.,negPerBin)
	
	# Always work on logarithmic parameter values
	if 0 in logPars: err.warning("  Setting all reliability.logPars entries to 1")
	logPars=[1 for pp in parSpace]
	
	# Import Matplotlib if diagnostic plots requested
	if makePlot:
		import matplotlib
		# The following line is necessary to run SoFiA remotely
		matplotlib.use("Agg")
		import matplotlib.pyplot as plt
	
	# --------------------------------
	# Build array of source parameters
	# --------------------------------
	
	idCOL   = parNames.index("id")
	ftotCOL = parNames.index("snr_sum")
	fmaxCOL = parNames.index("snr_max")
	fminCOL = parNames.index("snr_min")
	
	# Get columns of requested parameters
	parCol = []
	for ii in range(len(parSpace)): parCol.append(parNames.index(parSpace[ii]))
	
	# Get position and number of positive and negative sources
	pos  = data[:, ftotCOL] >  0
	neg  = data[:, ftotCOL] <= 0
	Npos = pos.sum()
	Nneg = neg.sum()
	
	err.ensure(Npos, "No positive sources found; cannot proceed.")
	err.ensure(Nneg, "No negative sources found; cannot proceed.")
	
	# Get array of relevant source parameters (and take log of them if requested)
	ids = data[:,idCOL]
	pars = np.empty((data.shape[0], 0))
	
	for ii in range(len(parSpace)):
		if parSpace[ii] == "snr_max":
			parsTmp = data[:,fmaxCOL] * pos - data[:,fminCOL] * neg
			if logPars[ii]: parsTmp = np.log10(parsTmp)
			pars = np.concatenate((pars, parsTmp.reshape(-1, 1)), axis=1)
		elif parSpace[ii] == "snr_sum" or parSpace[ii] == "snr_mean":
			parsTmp = abs(data[:,parCol[ii]].reshape(-1, 1))
			if logPars[ii]: parsTmp = np.log10(parsTmp)
			pars = np.concatenate((pars, parsTmp), axis=1)
		else:
			parsTmp = data[:,parCol[ii]].reshape(-1, 1)
			if logPars[ii]: parsTmp = np.log10(parsTmp)
			pars = np.concatenate((pars, parsTmp), axis=1)
	
	err.message("  Working in parameter space " + str(parSpace))
	pars = np.transpose(pars)
	
	
	# ----------------------------------------------------------
	# Set parameters to work with and gridding/plotting for each
	# ----------------------------------------------------------
	
	# Axis labels when plotting
	labs = []
	for ii in range(len(parSpace)):
		labs.append("")
		if logPars[ii]: labs[ii] += "log "
		labs[ii] += parSpace[ii]
	
	# Axis limits when plotting
	pmin, pmax = pars.min(axis=1), pars.max(axis=1)
	pmin, pmax = pmin - 0.1 * (pmax - pmin), pmax + 0.1 * (pmax - pmin)
	lims = [[pmin[i], pmax[i]] for i in range(len(parSpace))]
	
	# Grid on which to evaluate Np and Nn in order to plot contours
	grid = [[pmin[i], pmax[i], 0.02 * (pmax[i] - pmin[i])] for i in range(len(parSpace))]
	
	# Calculate the number of rows and columns in figure
	projections = [subset for subset in combinations(range(len(parSpace)), 2)]
	nr = int(np.floor(np.sqrt(len(projections))))
	nc = int(np.ceil(float(len(projections)) / nr))
	
	
	# ---------------------------------------
	# Set smoothing kernel in parameter space
	# ---------------------------------------
	
	# If autoKernel is True, then the initial kernel is taken as a scaled version of the covariance matrix
	# of the negative sources. The kernel size along each axis is such that the number of sources per kernel
	# width (sigma**2) is equal to "negPerBin". Optionally, the user can decide to use only the diagonal
	# terms of the covariance matrix. The kernel is then grown until convergence is reached on the Skellam
	# plot. If autoKernel is False, then use the kernel given by "kernel" parameter (argument of EstimateRel);
	# this is sigma, and is squared to be consistent with the auto kernel above.
	
	if autoKernel:
		# Set the kernel shape to that of the variance or covariance matrix
		kernel = np.cov(pars[:, neg])
		kernelType = "covariance"
		# Check if kernel matrix can be inverted
		try:
			np.linalg.inv(kernel)
		except:
			err.error(
				"The reliability cannot be calculated because the smoothing kernel\n"
				"derived from " + str(pars[:,neg].shape[1]) + " negative sources cannot be inverted.\n"
				"This is likely due to an insufficient number of negative sources.\n"
				"Try to increase the number of negative sources by changing the\n"
				"source finding and/or filtering settings.", fatal=True, frame=True)
		
		if np.isnan(kernel).sum():
			err.error(
				"The reliability cannot be calculated because the smoothing kernel\n"
				"derived from " + str(pars[:,neg].shape[1]) + " negative sources contains NaNs.\n"
				"A good kernel is required to calculate the density field of positive\n"
				"and negative sources in parameter space.\n"
				"Try to increase the number of negative sources by changing the\n"
				"source finding and/or filtering settings.", fatal=True, frame=True)
		
		if not usecov:
			kernel = np.diag(np.diag(kernel))
			kernelType = "variance"
		
		kernelIter = 0.0
		deltplot = []
		
		# Scale the kernel size as requested by the user or by the auto-scale algorithm (scaleKernel > 0)
		if scaleKernel:
			# Scale kernel size as requested by the user
			# Note that the scale factor is squared because users are asked to give a factor to apply to sqrt(kernel)
			kernel *= scaleKernel**2
			err.message("  Using kernel with shape of %s and size scaled by factor %.2f." % (kernelType, scaleKernel))
			err.message("  The sqrt(kernel) size is:")
			err.message(str(np.sqrt(np.abs(kernel))))
		else:
			# Scale kernel size to start the kernel-growing loop
			# The scale factor for sqrt(kernel) is elevated to the power of 1.0 / len(parCol)
			kernel *= ((negPerBin + kernelIter) / Nneg)**(2.0 / len(parCol))
			err.message("  Will find the best kernel as a scaled version of the %s:" % kernelType)
			err.message("  Starting from the kernel with sqrt(kernel) size:")
			err.message("  " + str(np.sqrt(np.abs(kernel))))
			err.message("  Growing kernel...")
		
		#deltOLD=-1e+9 # Used to stop kernel growth if P-N stops moving closer to zero [NOT USED CURRENTLY]
		if doskellam and makePlot: fig0 = plt.figure()
	else:
		# Note that the user must give sigma, which then gets squared
		err.message("  Using user-defined variance kernel with sqrt(kernel) size:")
		err.message(str(np.array(kernel)))
		kernel = np.identity(len(kernel)) * np.array(kernel)**2
	
	# Set grow_kernel to 1 to start the kernel growing loop below.
	grow_kernel = 1
	
	# This loop will estimate the reliability, check whether the kernel is large enough,
	# and if not pick a larger kernel. If autoKernel = 0 or scaleKernel = 0, we will do
	# just one pass (i.e., we will not grow the kernel).
	while grow_kernel:
		# ------------------------
		# Evaluate N-d reliability
		# ------------------------
		
		if verb: err.message("   estimate normalised positive and negative density fields ...")
		
		Np = gaussian_kde_set_covariance(pars[:,pos], kernel)
		Nn = gaussian_kde_set_covariance(pars[:,neg], kernel)
		
		# Calculate the number of positive and negative sources at the location of positive sources
		Nps = Np(pars[:,pos]) * Npos
		Nns = Nn(pars[:,pos]) * Nneg
		
		# Calculate the number of positive and negative sources at the location of negative sources
		nNps = Np(pars[:,neg]) * Npos
		nNns = Nn(pars[:,neg]) * Nneg
		
		# Calculate the reliability at the location of positive sources
		Rs = (Nps - Nns) / Nps
		
		# The reliability must be <= 1. If not, something is wrong.
		err.ensure(Rs.max() <= 1, "Maximum reliability greater than 1; something is wrong.\nPlease ensure that enough negative sources are detected\nand decrease your source finding threshold if necessary.", frame=True)
		
		# Find pseudo-reliable sources (taking maximum(Rs, 0) in order to include objects with Rs < 0
		# if threshold == 0; Rs may be < 0 because of insufficient statistics)
		# These are called pseudo-reliable because some objects may be discarded later based on additional criteria below
		pseudoreliable = np.maximum(Rs, 0) >= threshold

		# Find reliable sources (taking maximum(Rs, 0) in order to include objects with Rs < 0 if
		# threshold == 0; Rs may be < 0 because of insufficient statistics)
		#reliable=(np.maximum(Rs, 0)>=threshold) * (data[pos, ftotCOL].reshape(-1,) > fMin) * (data[pos, fmaxCOL].reshape(-1,) > 4)
		reliable = (np.maximum(Rs, 0) >= threshold) * ((data[pos, ftotCOL] / np.sqrt(data[pos, parNames.index("n_pix")])).reshape(-1,) > fMin)
		
		if autoKernel:
			# Calculate quantities needed for comparison to Skellam distribution
			delt = (nNps - nNns) / np.sqrt(nNps + nNns)
			deltstd = delt.std()
			deltmed = np.median(delt)
			deltmin = delt.min()
			deltmax = delt.max()
			
			if deltmed / deltstd > -100 and doskellam and makePlot:
				plt.hist(delt / deltstd, bins=np.arange(deltmin / deltstd, max(5.1, deltmax / deltstd), 0.01), cumulative=True, histtype="step", color=(min(1, float(negPerBin + kernelIter) / Nneg), 0,0), normed=True)
				deltplot.append([((negPerBin + kernelIter) / Nneg)**(1.0 / len(parCol)), deltmed / deltstd])
			
			err.message("  iteration, median, width, median/width = %3i, %9.2e, %9.2e, %9.2e" % (kernelIter, deltmed, deltstd, deltmed / deltstd))
			
			if scaleKernel: grow_kernel = 0
			elif deltmed / deltstd > skellamTol or negPerBin + kernelIter >= Nneg:
				grow_kernel = 0
				err.message("  Found good kernel after %i kernel growth iterations. The sqrt(kernel) size is:" % kernelIter)
				err.message(np.sqrt(np.abs(kernel)))
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
	
	
	# ------------
	# Skellam plot
	# ------------
	
	if autoKernel and deltmed / deltstd > -100 and doskellam and makePlot:
		plt.plot(np.arange(-10, 10, 0.01), stats.norm().cdf(np.arange(-10, 10, 0.01)), "k-")
		plt.plot(np.arange(-10, 10, 0.01), stats.norm(scale=0.4).cdf(np.arange(-10, 10, 0.01)), "k:")
		plt.legend(("Gaussian (sigma=1)", "Gaussian (sigma=0.4)"), loc="lower right", prop={"size":13})
		plt.hist(delt / deltstd, bins=np.arange(deltmin / deltstd, max(5.1, deltmax / deltstd), 0.01), cumulative=True, histtype="step", color="r", normed=True)
		plt.xlim(-5, 5)
		plt.ylim(0, 1)
		plt.xlabel("(P-N)/sqrt(N+P)")
		plt.ylabel("cumulative distribution")
		plt.plot([0, 0], [0, 1], "k--")
		fig0.savefig("%s_rel_skellam.pdf" % pdfoutname, rasterized=True)
		
		if not scaleKernel:
			fig3 = plt.figure()
			deltplot = np.array(deltplot)
			plt.plot(deltplot[:,0], deltplot[:,1], "ko-")
			plt.xlabel("kernel size (1D-sigma, aribtrary units)")
			plt.ylabel("median/std of (P-N)/sqrt(P+N)")
			plt.axhline(y=skellamTol, linestyle="--", color="r")
			fig3.savefig("%s_rel_skellam-delta.pdf" % pdfoutname, rasterized=True)
	
	
	# -----------------------
	# Scatter plot of sources
	# -----------------------
	
	specialids = []
	
	if doscatter and makePlot:
		if verb: err.message("  plotting sources ...")
		fig1 = plt.figure(figsize=(18, 4.5 * nr))
		plt.subplots_adjust(left=0.06, bottom=0.15/nr, right = 0.97, top=1-0.08/nr, wspace=0.35, hspace=0.25)
		
		n_p = 0
		for jj in projections:
			if verb: err.message("    projection %i/%i" % (projections.index(jj) + 1, len(projections)))
			n_p, p1, p2 = n_p + 1, jj[0], jj[1]
			plt.subplot(nr, nc, n_p)
			plt.scatter(pars[p1,pos], pars[p2,pos], marker="o", c="b", s=10, edgecolor="face", alpha=0.5)
			plt.scatter(pars[p1,neg], pars[p2,neg], marker="o", c="r", s=10, edgecolor="face", alpha=0.5)
			for si in specialids: plt.plot(pars[p1, ids==si], pars[p2, ids==si], "kd", zorder=10000, ms=7, mfc="none", mew=2)
			plt.xlim(lims[p1][0], lims[p1][1])
			plt.ylim(lims[p2][0], lims[p2][1])
			plt.xlabel(labs[p1])
			plt.ylabel(labs[p2])
			plt.grid(color='k',linestyle='-',linewidth=0.2)
		fig1.savefig("%s_rel_scatter.pdf" % pdfoutname, rasterized=True)
	
	
	# -------------
	# Plot contours
	# -------------
	
	if docontour and makePlot:
		levs = 10**np.arange(-1.5, 2, 0.5)
		
		if verb: err.message("  plotting contours ...")
		fig2 = plt.figure(figsize=(18, 4.5 * nr))
		plt.subplots_adjust(left=0.06, bottom=0.15/nr, right=0.97, top=1-0.08/nr, wspace=0.35, hspace=0.25)
		n_p = 0
		for jj in projections:
			if verb: err.message("    projection %i/%i" % (projections.index(jj) + 1, len(projections)))
			n_p, p1, p2 = n_p + 1, jj[0], jj[1]
			g1, g2 = grid[p1], grid[p2]
			x1 = np.arange(g1[0], g1[1], g1[2])
			x2 = np.arange(g2[0], g2[1], g2[2])
			pshape = (x2.shape[0], x1.shape[0])
			
			# Get array of source parameters on current projection
			parsp = np.concatenate((pars[p1:p1+1], pars[p2:p2+1]), axis=0)
			
			# Derive Np and Nn density fields on the current projection
			setcov = kernel[p1:p2+1:p2-p1,p1:p2+1:p2-p1]
			try:
				Np = gaussian_kde_set_covariance(parsp[:,pos], setcov)
				Nn = gaussian_kde_set_covariance(parsp[:,neg], setcov)
			except:
				err.error(
					"Reliability  determination  failed  because of issues  with the\n"
					"smoothing kernel.  This is likely due to an insufficient number\n"
					"of negative detections. Please review your filtering and source\n"
					"finding settings to ensure that a sufficient number of negative\n"
					"detections is found.", fatal=True, frame=True)
			
			# Evaluate density fields on grid on current projection
			g = np.transpose(np.transpose(np.mgrid[slice(g1[0], g1[1], g1[2]), slice(g2[0], g2[1], g2[2])]).reshape(-1, 2))
			Np = Np(g)
			Nn = Nn(g)
			Np = Np / Np.sum() * Npos
			Nn = Nn / Nn.sum() * Nneg
			Np.resize(pshape)
			Nn.resize(pshape)
			plt.subplot(nr, nc, n_p)
			plt.contour(x1, x2, Np, origin="lower", colors="b", levels=levs, zorder=2)
			plt.contour(x1, x2, Nn, origin="lower", colors="r", levels=levs, zorder=1)

			# Plot Integrated SNR threshold
			if fMin>0:
				if (parSpace[jj[0]],parSpace[jj[1]])==("snr_sum","snr_mean"):
					xArray=np.arange(lims[p1][0],lims[p1][1]+(lims[p1][1]-lims[p1][0])/100,(lims[p1][1]-lims[p1][0])/100)
					plt.plot(xArray,np.log10(fMin)*2-xArray,'k:')
				if (parSpace[jj[0]],parSpace[jj[1]])==("snr_mean","snr_sum"):
					yArray=np.arange(lims[p2][0],lims[p2][1]+(lims[p2][1]-lims[p2][0])/100,(lims[p2][1]-lims[p2][0])/100)
					plt.plot(np.log10(fMin)*2-yArray,yArray,'k:')

			if reliable.sum(): plt.scatter(pars[p1,pos][reliable], pars[p2,pos][reliable], marker="o", s=10, edgecolor="k", facecolor="k", zorder=4)
			if (pseudoreliable * (reliable == False)).sum(): plt.scatter(pars[p1,pos][pseudoreliable * (reliable == False)], pars[p2,pos][pseudoreliable * (reliable == False)], marker="x", s=40, edgecolor="0.5", facecolor="0.5", zorder=3)
			for si in specialids: plt.plot(pars[p1,ids==si], pars[p2,ids==si], "kd", zorder=10000, ms=7, mfc="none", mew=2)
			plt.xlim(lims[p1][0], lims[p1][1])
			plt.ylim(lims[p2][0], lims[p2][1])
			plt.xlabel(labs[p1])
			plt.ylabel(labs[p2])
			plt.grid(color='k',linestyle='-',linewidth=0.2)
		fig2.savefig("%s_rel_contour.pdf" % pdfoutname, rasterized=True)
	
	
	# -------------------------
	# Add Np, Nn and R to table
	# -------------------------
	
	# This allows me not to calculate R every time I want to do some plot analysis,
	# but just read it from the file
	if saverel:
		if not (docontour or dostats):
			Nps = Np(pars[:,pos]) * Npos
			Nns = Nn(pars[:,pos]) * Nneg
		Np = np.zeros((data.shape[0],))
		Np[pos] = Nps
		Nn = np.zeros((data.shape[0],))
		Nn[pos] = Nns
		R = -np.ones((data.shape[0],)) # R will be -1 for negative sources
		# Set R to zero for positive sources if R < 0 because of Nn > Np
		R[pos] = np.maximum(0, (Np[pos] - Nn[pos]) / Np[pos])
		data = np.concatenate((data, Np.reshape(-1, 1), Nn.reshape(-1, 1), R.reshape(-1, 1)), axis=1)
	
	data = [list(jj) for jj in list(data)]
	return data, ids[pos][reliable].astype(int)
