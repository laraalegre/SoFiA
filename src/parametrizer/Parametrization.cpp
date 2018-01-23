#include <iostream>
#include <cmath>
#include <limits>
#include <new>

#include "helperFunctions.h"
#include "DataCube.h"
#include "Parametrization.h"
#include "BFfit.h"
#include "Measurement.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SOURCE_LOOP_START for(long x = subRegionX1; x <= subRegionX2; ++x) { for(long y = subRegionY1; y <= subRegionY2; ++y) { for(long z = subRegionZ1; z <= subRegionZ2; ++z) { double fluxValue = static_cast<double>(dataCube->getData(x, y, z)); if(static_cast<unsigned short>(maskCube->getData(x, y, z)) == source->getSourceID())
#define SOURCE_LOOP_END } } }



// ----------- //
// Constructor //
// ----------- //

Parametrization::Parametrization()
{
	// Initialisation of data:
	dataCube = 0;
	maskCube = 0;
	source   = 0;
	dataSize = 0;
	
	// Initialisation of all parameters:
	noiseSubCube         = 0.0;
	noiseSpectrumMax     = 0.0;
	centroidX            = 0.0;
	centroidY            = 0.0;
	centroidZ            = 0.0;
	errCentroidX         = 0.0;
	errCentroidY         = 0.0;
	errCentroidZ         = 0.0;
	lineWidthW20         = 0.0;
	lineWidthW50         = 0.0;
	errlineWidthW20      = 0.0;
	errlineWidthW50      = 0.0;
	lineWidthWm50        = 0.0;
	meanFluxWm50         = 0.0;
	peakFlux             = 0.0;
	peakFluxSpec         = 0.0;
	totalFlux            = 0.0;
	intSNR               = 0.0;
	ellMaj               = 0.0;
	ellMin               = 0.0;
	ellPA                = 0.0;
	ell3sMaj             = 0.0;
	ell3sMin             = 0.0;
	ell3sPA              = 0.0;
	kinematicPA          = 0.0;
	flagKinematicPA      = false;
	flagWarp             = false;
	
	busyFitSuccess       = 0;
	busyFunctionChi2     = 0.0;
	busyFunctionCentroid = 0.0;
	busyFunctionW20      = 0.0;
	busyFunctionW50      = 0.0;
	busyFunctionFpeak    = 0.0;
	busyFunctionFint     = 0.0;
	
	for(size_t i = 0; i < 2 * BUSYFIT_FREE_PARAM + 1; i++)
	{
		busyFitParameters[i]    = 0.0;
		//busyFitUncertainties[i] = 0.0;
	}
	
	return;
}



// -------------------------------------------------- //
// Interface function to call parameterisation module //
// -------------------------------------------------- //

int Parametrization::parametrize(DataCube<float> *d, DataCube<short> *m, Source *s, bool doBF)
{
	doBusyFunction = doBF;
	
	if(loadData(d, m, s) != 0)
	{
		std::cerr << "Error (Parametrization): No data found; source parametrisation failed.\n";
		return 1;
	}
	
	if(createIntegratedSpectrum() != 0)
	{
		std::cerr << "Error (Parametrization): Failed to create integrated spectrum.\n";
		return 1;
	}
	
	if(measureCentroid() != 0)
	{
		std::cerr << "Warning (Parametrization): Failed to measure source centroid.\n";
	}
	
	if(measureFlux() != 0)
	{
		std::cerr << "Warning (Parametrization): Source flux measurement failed.\n";
	}
	
	if(measureLineWidth() != 0)
	{
		std::cerr << "Warning (Parametrization): Failed to measure source line width.\n";
	}
	
	if(fitEllipse() != 0)
	{
		std::cerr << "Warning (Parametrization): Ellipse fit failed.\n";
	}
	
	if(kinematicMajorAxis() != 0)
	{
		std::cerr << "Warning (Parametrization): Measurement of kinematic PA failed.\n";
	}
	
	if(doBusyFunction)
	{
		if(fitBusyFunction() != 0)
		{
			std::cerr << "Warning (Parametrization): Failed to fit Busy Function.\n";
		}
	}
	
	if(writeParameters() != 0)
	{
		std::cerr << "Error (Parametrization): Failed to write parameters to source.\n";
		return 1;
	}
	
	return 0;
}



// ----------------------------------------------- //
// Function to load data and determine local noise //
// ----------------------------------------------- //

int Parametrization::loadData(DataCube<float> *d, DataCube<short> *m, Source *s)
{
	dataCube = 0;
	maskCube = 0;
	source = 0;
	
	//data.clear();            // Clear all previously defined data.
	
	if(d == 0 or m == 0 or s == 0)
	{
		std::cerr << "Error (Parametrization): Cannot load data; invalid pointer provided.\n";
		return 1;
	}
	
	if(!d->isDefined() or !m->isDefined() or !s->isDefined())
	{
		std::cerr << "Error (Parametrization): Cannot load data; source or data cube undefined.\n";
		return 1;
	}
	
	if(d->getSize(0) != m->getSize(0) or d->getSize(1) != m->getSize(1) or d->getSize(2) != m->getSize(2))
	{
		std::cerr << "Error (Parametrization): Mask and data cube have different sizes.\n";
		return 1;
	}
	
	double posX = s->getParameter("x");
	double posY = s->getParameter("y");
	double posZ = s->getParameter("z");
	
	if(posX < 0.0 or posY < 0.0 or posZ < 0.0 or posX >= static_cast<double>(d->getSize(0)) or posY >= static_cast<double>(d->getSize(1)) or posZ >= static_cast<double>(d->getSize(2)))
	{
		std::cerr << "Error (Parametrization): Source position outside cube range.\n";
		return 1;
	}
	
	dataCube = d;
	maskCube = m;
	source   = s;
	
	// Define sub-region to operate on:
	
	if(source->parameterDefined("x_min") and source->parameterDefined("x_max"))
	{
		searchRadiusX = static_cast<long>(source->getParameter("x_max") - source->getParameter("x_min"));
	}
	else
	{
		searchRadiusX = PARAMETRIZATION_DEFAULT_SPATIAL_RADIUS;
		std::cerr << "Warning (MaskOptimization): No bounding box defined; using default search radius\n";
		std::cerr << "                            in the spatial domain instead.\n";
	}
	
	if(source->parameterDefined("y_min") and source->parameterDefined("y_max"))
	{
		searchRadiusY = static_cast<long>(source->getParameter("y_max") - source->getParameter("y_min"));
	}
	else
	{
		searchRadiusY = PARAMETRIZATION_DEFAULT_SPATIAL_RADIUS;
		std::cerr << "Warning (MaskOptimization): No bounding box defined; using default search radius\n";
		std::cerr << "                            in the spatial domain instead.\n";
	}
	
	if(source->parameterDefined("z_min") and source->parameterDefined("z_max"))
	{
		searchRadiusZ = static_cast<long>(0.6 * (source->getParameter("z_max") - source->getParameter("z_min")));
	}
	else
	{
		searchRadiusZ = PARAMETRIZATION_DEFAULT_SPECTRAL_RADIUS;
		std::cerr << "Warning (MaskOptimization): No bounding box defined; using default search radius\n";
		std::cerr << "                            in the spectral domain instead.\n";
	}
	
	subRegionX1 = static_cast<long>(posX) - searchRadiusX;
	subRegionX2 = static_cast<long>(posX) + searchRadiusX;
	subRegionY1 = static_cast<long>(posY) - searchRadiusY;
	subRegionY2 = static_cast<long>(posY) + searchRadiusY;
	subRegionZ1 = static_cast<long>(posZ) - searchRadiusZ;
	subRegionZ2 = static_cast<long>(posZ) + searchRadiusZ;
	
	if(subRegionX1 < 0L) subRegionX1 = 0L;
	if(subRegionY1 < 0L) subRegionY1 = 0L;
	if(subRegionZ1 < 0L) subRegionZ1 = 0L;
	if(subRegionX2 >= dataCube->getSize(0)) subRegionX2 = dataCube->getSize(0) - 1L;
	if(subRegionY2 >= dataCube->getSize(1)) subRegionY2 = dataCube->getSize(1) - 1L;
	if(subRegionZ2 >= dataCube->getSize(2)) subRegionZ2 = dataCube->getSize(2) - 1L;
	
	// Extract all pixels belonging to the source and calculate local noise:
	dataSize = 0;
	noiseSubCube = 0.0;
	std::vector<double> rmsMad;
	
	SOURCE_LOOP_START ++dataSize;
	else if(maskCube->getData(x, y, z) == 0 and not std::isnan(fluxValue)) rmsMad.push_back(static_cast<double>(fluxValue));
	SOURCE_LOOP_END
	
	if(dataSize == 0)
	{
		std::cerr << "Error (Parametrization): No data found for source " << source->getSourceID() << ".\n";
		return 1;
	}
	
	if(rmsMad.size() == 0)
	{
		std::cerr << "Warning (Parametrization): Noise calculation failed for source " << source->getSourceID() << ".\n";
	}
	else
	{
		// Calculate the rms via the median absolute deviation (MAD):
		double m = median(rmsMad);
		for(size_t i = 0; i < rmsMad.size(); ++i) rmsMad[i] = fabs(rmsMad[i] - m);
		noiseSubCube = 1.4826 * median(rmsMad);
		// NOTE: The factor of 1.4826 above is the approximate theoretical conversion factor
		//       between the MAD and the standard deviation under the assumption that the 
		//       data samples follow a normal distribution.
	}
	
	return 0;
}



// ---------------- //
// Measure centroid //
// ---------------- //

int Parametrization::measureCentroid()
{
	if(dataSize == 0)
	{
		std::cerr << "Error (Parametrization): No valid data found.\n";
		return 1;
	}
	
	double sum   = 0.0;
	centroidX    = 0.0;
	centroidY    = 0.0;
	centroidZ    = 0.0;
	errCentroidX = 0.0;
	errCentroidY = 0.0;
	errCentroidZ = 0.0;
	
	// Centroid:
	SOURCE_LOOP_START
	{
		if(fluxValue > 0.0)                // NOTE: Only positive pixels considered here!
		{
			centroidX += fluxValue * static_cast<double>(x);
			centroidY += fluxValue * static_cast<double>(y);
			centroidZ += fluxValue * static_cast<double>(z);
			sum       += fluxValue;
		}
	} SOURCE_LOOP_END
	
	centroidX /= sum;
	centroidY /= sum;
	centroidZ /= sum;
	
	// Uncertainties:
	SOURCE_LOOP_START
	{
		if(fluxValue > 0.0)                // NOTE: Only positive pixels considered here!
		{
			errCentroidX += (static_cast<double>(x) - centroidX) * (static_cast<double>(x) - centroidX);
			errCentroidY += (static_cast<double>(y) - centroidY) * (static_cast<double>(y) - centroidY);
			errCentroidZ += (static_cast<double>(z) - centroidZ) * (static_cast<double>(z) - centroidZ);
		}
	} SOURCE_LOOP_END
	
	errCentroidX = sqrt(errCentroidX) * noiseSubCube / sum;
	errCentroidY = sqrt(errCentroidY) * noiseSubCube / sum;
	errCentroidZ = sqrt(errCentroidZ) * noiseSubCube / sum;
	
	return 0;
}



// ----------------------------------------------------- //
// Measure peak flux, integrated flux and integrated S/N //
// ----------------------------------------------------- //

int Parametrization::measureFlux()
{
	if(dataSize == 0)
	{
		std::cerr << "Error (Parametrization): No valid data found.\n";
		return 1;
	}
	
	totalFlux = 0.0;
	peakFlux  = -std::numeric_limits<double>::max();
	
	// Sum over all pixels (including negative ones):
	SOURCE_LOOP_START
	{
		totalFlux += fluxValue;
		if(peakFlux < fluxValue) peakFlux = fluxValue;
	} SOURCE_LOOP_END
	
	// Calculate integrated SNR:
	intSNR = totalFlux / (noiseSubCube * sqrt(static_cast<double>(dataSize)));
	// WARNING The integrated SNR would need to be divided by the beam solid angle,
	// WARNING but this will need to be done outside this module as no header
	// WARNING information known at this stage. The current value would only be
	// WARNING correct if the noise in adjacent pixels were uncorrelated.
	
	return 0;
}



// --------------------- //
// Fit ellipse to source //
// --------------------- //

int Parametrization::fitEllipse()
{
	// (1) Fitting ellipse to intensity-weighted image:
	
	if(dataSize == 0)
	{
		std::cerr << "Error (Parametrization): No valid data found.\n";
		return 1;
	}
	
	if(totalFlux <= 0.0)
	{
		std::cerr << "Error (Parametrization): Cannot fit ellipse, source flux <= 0.\n";
		return 1;
	}
	
	double momX  = 0.0;
	double momY  = 0.0;
	double momXY = 0.0;
	double sum   = 0.0;
	
	SOURCE_LOOP_START
	{
		if(fluxValue > 0.0)            // NOTE: Only positive pixels considered here!
		{
			momX  += static_cast<double>((x - source->getParameter("x")) * (x - source->getParameter("x"))) * fluxValue;
			momY  += static_cast<double>((y - source->getParameter("y")) * (y - source->getParameter("y"))) * fluxValue;
			momXY += static_cast<double>((x - source->getParameter("x")) * (y - source->getParameter("y"))) * fluxValue;
			sum += fluxValue;
		}
	} SOURCE_LOOP_END
	
	momX  /= sum;
	momY  /= sum;
	momXY /= sum;
	
	ellPA  = 0.5 * atan2(2.0 * momXY, momX - momY);
	ellMaj = sqrt(2.0 * (momX + momY + sqrt((momX - momY) * (momX - momY) + 4.0 * momXY * momXY)));
	ellMin = sqrt(2.0 * (momX + momY - sqrt((momX - momY) * (momX - momY) + 4.0 * momXY * momXY)));
	
	// WARNING: Converting PA from rad to deg:
	ellPA *= 180.0 / M_PI;
	
	// WARNING: Subtracting 90° from PA here, because astronomers like to have 0° pointing up.
	//          This means that PA will no longer have the mathematically correct orientation!
	ellPA -= 90.0;
	if(ellPA < -90.0) ellPA += 180;
	// NOTE:    PA should now be between -90° and +90°.
	//          Note that the PA is relative to the pixel grid, not the coordinate system!
	
	
	
	// (2) Fitting ellipse to moment map of uniformly weighted pixels above 3 sigma:
	
	size_t sizeX = subRegionX2 - subRegionX1 + 1;
	size_t sizeY = subRegionY2 - subRegionY1 + 1;
	
	double *momentMap;
	long   *maskMap;
	
	try
	{
		momentMap = new double[sizeX * sizeY];
		maskMap   = new long[sizeX * sizeY];
	}
	catch(std::bad_alloc &badAlloc)
	{
		std::cerr << "Error (Parametrization): Memory allocation failed: " << badAlloc.what() << '\n';
		return 1;
	}
	
	for(size_t i = 0; i < sizeX * sizeY; ++i)
	{
		momentMap[i] = 0.0;
		maskMap[i] = 0;
	}
	
	SOURCE_LOOP_START
	{
		momentMap[x - subRegionX1 + sizeX * (y - subRegionY1)] += fluxValue;
		maskMap[x - subRegionX1 + sizeX * (y - subRegionY1)] += 1;
	} SOURCE_LOOP_END
	
	for(size_t x = 0; x < sizeX; ++x)
	{
		for(size_t y = 0; y < sizeY; ++y)
		{
			// Mask all valid pixels above 3 × RMS:
			if(maskMap[x + sizeX * y] > 0)
			{
				if(momentMap[x + sizeX * y] / (sqrt(static_cast<double>(maskMap[x + sizeX * y])) * noiseSubCube) > 3.0) maskMap[x + sizeX * y] = 1;
				else maskMap[x + sizeX * y] = 0;
			}
		}
	}
	
	delete[] momentMap;      // No longer needed, so let's release the memory again.
	
	momX  = 0.0;
	momY  = 0.0;
	momXY = 0.0;
	sum   = 0.0;
	
	// Determine ellipse parameters on mask (where all pixels > 3 sigma have been
	// set to 1, all other pixels to 0):
	for(size_t x = 0; x < sizeX; ++x)
	{
		for(size_t y = 0; y < sizeY; ++y)
		{
			momX  += static_cast<double>((x - source->getParameter("x") + subRegionX1) * (x - source->getParameter("x") + subRegionX1) * maskMap[x + sizeX * y]);
			momY  += static_cast<double>((y - source->getParameter("y") + subRegionY1) * (y - source->getParameter("y") + subRegionY1) * maskMap[x + sizeX * y]);
			momXY += static_cast<double>((x - source->getParameter("x") + subRegionX1) * (y - source->getParameter("y") + subRegionY1) * maskMap[x + sizeX * y]);
			sum += static_cast<double>(maskMap[x + sizeX * y]);
		}
	}
	
	delete[] maskMap;        // No longer needed, so let's release the memory again.
	
	momX  /= sum;
	momY  /= sum;
	momXY /= sum;
	
	ell3sPA  = 0.5 * atan2(2.0 * momXY, momX - momY);
	ell3sMaj = sqrt(2.0 * (momX + momY + sqrt((momX - momY) * (momX - momY) + 4.0 * momXY * momXY)));
	ell3sMin = sqrt(2.0 * (momX + momY - sqrt((momX - momY) * (momX - momY) + 4.0 * momXY * momXY)));
	
	// WARNING: Converting PA from rad to deg:
	ell3sPA *= 180.0 / M_PI;
	
	// WARNING: Subtracting 90° from PA here, because astronomers like to have 0° pointing up.
	//          This means that PA will no longer have the mathematically correct orientation!
	ell3sPA -= 90.0;
	if(ell3sPA < -90.0) ell3sPA += 180;
	// NOTE:    PA should now be between -90° and +90°.
	//          Note that the PA is relative to the pixel grid, not the coordinate system!
	
	return 0;
}



// ---------------------------- //
// Measure kinematic major axis //
// ---------------------------- //

int Parametrization::kinematicMajorAxis()
{
	if(dataSize == 0)
	{
		std::cerr << "Error (Parametrization): No valid data found.\n";
		return 1;
	}
	
	// Determine the flux-weighted centroid in each channel:
	size_t size = subRegionZ2 - subRegionZ1 + 1;
	std::vector <double> cenX(size, 0.0);
	std::vector <double> cenY(size, 0.0);
	std::vector <double> sum(size, 0.0);
	size_t firstPoint = size;
	size_t lastPoint  = 0;
	
	SOURCE_LOOP_START
	{
		if(fluxValue > 3.0 * noiseSubCube)       // NOTE: Only values > 3 sigma considered here!
		{
			cenX[z - subRegionZ1] += fluxValue * static_cast<double>(x);
			cenY[z - subRegionZ1] += fluxValue * static_cast<double>(y);
			sum [z - subRegionZ1] += fluxValue;
		}
	} SOURCE_LOOP_END
	
	unsigned int counter = 0;
	
	for(size_t i = 0; i < size; ++i)
	{
		if(sum[i] > 0.0)
		{
			cenX[i] /= sum[i];
			cenY[i] /= sum[i];
			++counter;
			
			if(i < firstPoint) firstPoint = i;
			if(i > lastPoint)  lastPoint  = i;
		}
	}
	
	if(counter < 2)
	{
		std::cerr << "Warning (Parametrization): Cannot determine kinematic major axis. Source too faint.\n";
		flagKinematicPA = true;
		return 1;
	}
	else if(counter == 2)
	{
		std::cerr << "Warning (Parametrization): Kinematic major axis derived from just 2 data points.\n";
		flagKinematicPA = true;
	}
	
	// Fit a straight line to the set of centroids (using user-defined weights):
	double sumW   = 0.0;
	double sumX   = 0.0;
	double sumY   = 0.0;
	double sumXX  = 0.0;
	double sumYY  = 0.0;
	double sumXY  = 0.0;
	double weight = 1.0;
	
	// Using linear regression:
	/*for(size_t i = 0; i < size; ++i)
	{
		if(sum[i] > 0.0)
		{
			// Set the desired weights here (defaults to 1 otherwise):
			//weight = sum[i];
			
			sumW  += weight;
			sumX  += weight * cenX[i];
			sumY  += weight * cenY[i];
			sumXX += weight * cenX[i] * cenX[i];
			sumXY += weight * cenX[i] * cenY[i];
		}
	}
	
	double slope = (sumW * sumXY  - sumX * sumY)  / (sumW * sumXX - sumX * sumX);
	//double inter = (sumXX * sumY - sumX * sumXY) / (sum * sumXX - sumX * sumX);*/
	
	// Using Deming (orthogonal) regression:
	for(size_t i = 0; i < size; ++i)
	{
		if(sum[i] > 0.0)
		{
			// Set the desired weights here (defaults to 1 otherwise):
			weight = sum[i] * sum[i];
			
			sumW += weight;
			sumX += weight * cenX[i];
			sumY += weight * cenY[i];
		}
	}
	
	sumX /= sumW;
	sumY /= sumW;
	
	for(size_t i = 0; i < size; ++i)
	{
		if(sum[i] > 0.0)
		{
			// Set the desired weights here (defaults to 1 otherwise):
			weight = sum[i] * sum[i];
			
			sumXX += weight * (cenX[i] - sumX) * (cenX[i] - sumX);
			sumYY += weight * (cenY[i] - sumY) * (cenY[i] - sumY);
			sumXY += weight * (cenX[i] - sumX) * (cenY[i] - sumY);
		}
	}
	
	double slope = (sumYY - sumXX + sqrt((sumYY - sumXX) * (sumYY - sumXX) + 4.0 * sumXY * sumXY)) / (2.0 * sumXY);
	//double inter = sumY - slope * sumX;
	
	// Calculate position angle of kinematic major axis:
	kinematicPA = atan(slope);
	
	// Check orientation of approaching/receding side of disc:
	double fullAngle = atan2(cenY[lastPoint] - cenY[firstPoint], cenX[lastPoint] - cenX[firstPoint]);
	
	// Correct for full angle and astronomers' favourite definition of PA:
	double difference = fabs(atan2(sin(fullAngle) * cos(kinematicPA) - cos(fullAngle) * sin(kinematicPA), cos(fullAngle) * cos(kinematicPA) + sin(fullAngle) * sin(kinematicPA)));
	if(difference > M_PI / 2.0)
	{
		kinematicPA += M_PI;
		difference  -= M_PI;
	}
	
	if(fabs(difference) > M_PI / 6.0) flagWarp = true;     // WARNING: Warping angle of 30° hard-coded here!
	
	kinematicPA = 180.0 * kinematicPA / M_PI - 90.0;
	while(kinematicPA <    0.0) kinematicPA += 360.0;
	while(kinematicPA >= 360.0) kinematicPA -= 360.0;
	// NOTE: PA should now be between 0° (pointing up) and 360° and refer to the side of the disc
	//       that occupies the upper end of the channel range covered by the source.
	//       PAs are relative to the pixel grid, not the WCS!
	
	return 0;
}



// -------------------------- //
// Create integrated spectrum //
// -------------------------- //

int Parametrization::createIntegratedSpectrum()
{
	if(dataSize == 0)
	{
		std::cerr << "Error (Parametrization): No valid data found.\n";
		return 1;
	}
	
	spectrum.clear();
	noiseSpectrum.clear();
	
	for(long i = 0; i <= subRegionZ2 - subRegionZ1; ++i)
	{
		spectrum.push_back(0.0);
		noiseSpectrum.push_back(0.0);
	}
	
	std::vector<size_t> counter(spectrum.size(), 0);
	
	// Extract spectrum...
	SOURCE_LOOP_START
	{
		spectrum[z - subRegionZ1] += fluxValue;
		counter [z - subRegionZ1] += 1;
	} SOURCE_LOOP_END
	
	// ...and determine noise per channel:
	// WARNING: This still needs some consideration. Is is wise to provide the rms per
	// WARNING: channel as the uncertainty? Or would one rather use the S/N? If the rms
	// WARNING: is zero (because there are no data) the BF fitting will just produce NaNs.
	// WARNING: Even worse: the noise will not scale with sqrt(N), because pixels are
	// WARNING: spatially correlated!
	peakFluxSpec = 0.0;
	noiseSpectrumMax = 0.0;
	for(size_t i = 0; i < noiseSpectrum.size(); ++i)
	{
		if(counter[i] > 0) noiseSpectrum[i] = sqrt(static_cast<double>(counter[i])) * noiseSubCube;
		else noiseSpectrum[i] = std::numeric_limits<double>::infinity();
		
		// Record peak flux density and maximum noise level:
		if(std::isfinite(spectrum[i]) and peakFluxSpec < spectrum[i]) peakFluxSpec = spectrum[i];
		if(std::isfinite(noiseSpectrum[i]) and noiseSpectrumMax < noiseSpectrum[i]) noiseSpectrumMax = noiseSpectrum[i];
	}
	
	return 0;
}



// ------------------ //
// Measure line width //
// ------------------ //

int Parametrization::measureLineWidth()
{
	if(dataSize == 0 or spectrum.empty())
	{
		std::cerr << "Error (Parametrization): No valid data found.\n";
		return 1;
	}
	
	// Determine maximum:
	double specMax = 0.0;    // WARNING: This assumes that sources are always positive!
	
	for(size_t i = 0; i < spectrum.size(); ++i)
	{
		if(spectrum[i] > specMax) specMax = spectrum[i];
	}
	
	
	
	// (1) Determine w50:
	
	size_t i = 0;
	
	while(i < spectrum.size() and spectrum[i] < specMax / 2.0) ++i;
	
	if(i >= spectrum.size())
	{
		std::cerr << "Error (Parametrization): Calculation of W50 failed (1).\n";
		lineWidthW50 = 0.0;
		return 1;
	}
	
	lineWidthW50 = static_cast<double>(i);
	if(i > 0) lineWidthW50 -= (spectrum[i] - specMax / 2.0) / (spectrum[i] - spectrum[i - 1]);     // Interpolate if not at edge.
	
	i = spectrum.size() - 1;
	
	while(i >= 0 and spectrum[i] < specMax / 2.0) --i;
	
	if(i < 0)
	{
		std::cerr << "Error (Parametrization): Calculation of W50 failed (2).\n";
		lineWidthW50 = 0.0;
		return 1;
	}
	
	lineWidthW50 = static_cast<double>(i) - lineWidthW50;
	if(i < spectrum.size() - 1) lineWidthW50 += (spectrum[i] - specMax / 2.0) / (spectrum[i] - spectrum[i + 1]);  // Interpolate if not at edge.
	
	if(lineWidthW50 <= 0.0)
	{
		std::cerr << "Error (Parametrization): Calculation of W50 failed (3).\n";
		lineWidthW50 = 0.0;
		return 1;
	}
	
	// Uncertainty
	errlineWidthW50 = noiseSpectrumMax * lineWidthW50 / peakFluxSpec;
	// WARNING: There are an awful lot of assumptions going into this uncertainty estimate,
	//          including the assumptions of a Gaussian line profile and a constant RMS
	//          across the spectrum! Furthermore, the RMS is estimated as sqrt(N) times 
	//          the local RMS in the cube. Therefore, this should only be considered as a 
	//          rough estimate of the uncertainty (the same holds for w20)!
	
	
	
	// (2) Determine w20:
	
	i = 0;
	
	while(i < spectrum.size() and spectrum[i] < specMax / 5.0) ++i;
	
	if(i >= spectrum.size())
	{
		std::cerr << "Error (Parametrization): Calculation of W20 failed (1).\n";
		lineWidthW20 = 0.0;
		return 1;
	}
	
	lineWidthW20 = static_cast<double>(i);
	if(i > 0) lineWidthW20 -= (spectrum[i] - specMax / 5.0) / (spectrum[i] - spectrum[i - 1]);     // Interpolate if not at edge.
	
	i = spectrum.size() - 1;
	
	while(i >= 0 and spectrum[i] < specMax / 5.0) --i;
	
	if(i < 0)
	{
		std::cerr << "Error (Parametrization): Calculation of W20 failed (2).\n";
		lineWidthW20 = 0.0;
		return 1;
	}
	
	lineWidthW20 = static_cast<double>(i) - lineWidthW20;
	if(i < spectrum.size() - 1) lineWidthW20 += (spectrum[i] - specMax / 5.0) / (spectrum[i] - spectrum[i + 1]);  // Interpolate if not at edge.
	
	if(lineWidthW20 <= 0.0)
	{
		std::cerr << "Error (Parametrization): Calculation of W20 failed (3).\n";
		lineWidthW20 = 0.0;
		return 1;
	}
	
	// Uncertainty
	errlineWidthW20 = noiseSpectrumMax * lineWidthW20 / peakFluxSpec;
	// WARNING: See warnings in w50 section above.
	
	
	
	// (3) Determine Wm50:
	
	double sum = 0.0;
	double bound90l = 0.0;
	double bound90r = 0.0;
	lineWidthWm50 = 0.0;
	meanFluxWm50 = 0.0;
	i = 0;
	
	while(sum < 0.05 * totalFlux and i < spectrum.size())
	{
		sum += spectrum[i];
		++i;
	}
	
	if(i >= spectrum.size())
	{
		std::cerr << "Error (Parametrization): Calculation of Wm50 failed (1).\n";
		return 1;
	}
	
	i--;
	
	bound90l = static_cast<double>(i);
	if(i > 0) bound90l -= (spectrum[i] - 0.05 * totalFlux) / (spectrum[i] - spectrum[i - 1]);      // Interpolate if not at edge.
	
	i = spectrum.size() - 1;
	sum = 0.0;
	
	while(sum < 0.05 * totalFlux and i >= 0)
	{
		sum += spectrum[i];
		--i;
	}
	
	if(i < 0)
	{
		std::cerr << "Error (Parametrization): Calculation of Wm50 failed (2).\n";
		return 1;
	}
	
	++i;
	
	bound90r = static_cast<double>(i);
	if(i < spectrum.size() - 1) bound90r += (spectrum[i] - 0.05 * totalFlux) / (spectrum[i] - spectrum[i + 1]);   // Interpolate if not at edge.
	
	meanFluxWm50 = 0.9 * totalFlux / (bound90r - bound90l);
	
	if(meanFluxWm50 <= 0)
	{
		std::cerr << "Error (Parametrization): Calculation of Wm50 failed (3).\n";
		meanFluxWm50 = 0.0;
		return 1;
	}
	
	i = 0;
	
	while(spectrum[i] < 0.5 * meanFluxWm50 and i < spectrum.size()) ++i;
	
	if(i >= spectrum.size())
	{
		std::cerr << "Error (Parametrization): Calculation of Wm50 failed (4).\n";
		meanFluxWm50 = 0.0;
		return 1;
	}
	
	lineWidthWm50 = static_cast<double>(i);
	if(i > 0) lineWidthWm50 -= (spectrum[i] - 0.5 * meanFluxWm50) / (spectrum[i] - spectrum[i - 1]);      // Interpolate if not at edge.
	
	i = spectrum.size() - 1;
	
	while(spectrum[i] < 0.5 * meanFluxWm50 and i >= 0) --i;
	
	if(i < 0)
	{
		std::cerr << "Error (Parametrization): Calculation of Wm50 failed (5).\n";
		lineWidthWm50 = 0.0;
		meanFluxWm50 = 0.0;
		return 1;
	}
	
	lineWidthWm50 = static_cast<double>(i) - lineWidthWm50;
	if(i < spectrum.size() - 1) lineWidthWm50 += (spectrum[i] - 0.5 * meanFluxWm50) / (spectrum[i] - spectrum[i + 1]);   // Interpolate if not at edge.
	
	if(lineWidthWm50 <= 0)
	{
		std::cerr << "Error (Parametrization): Calculation of Wm50 failed (6).\n";
		lineWidthWm50 = 0.0;
		meanFluxWm50 = 0.0;
		return 1;
	}
	
	return 0;
}



// ----------------- //
// Fit Busy Function //
// ----------------- //

int Parametrization::fitBusyFunction()
{
	if(dataSize == 0 or spectrum.empty())
	{
		std::cerr << "Error (Parametrization): No valid data found.\n";
		return 1;
	}
	
	// Create spectral axis:
	std::vector<double> channels;
	for(size_t i = 0; i < spectrum.size(); ++i) channels.push_back(static_cast<double>(i));
	
	// Create covariance matrix:
	double **fitCov;
	fitCov = new double*[BUSYFIT_FREE_PARAM];
	for(int i = 0; i < BUSYFIT_FREE_PARAM; ++i) fitCov[i] = new double[BUSYFIT_FREE_PARAM];
	
	// Define a few parameters:
	int bestNoP = 8;
	int nSeeds  = 1000;
	int iterMax = 30;
	int verbose = 0;
	
	// Carry out the fitting:
	busyFitSuccess = FitBusyFunc(spectrum.size(), &channels[0], &spectrum[0], &noiseSpectrum[0], &busyFitParameters[0], fitCov, bestNoP, nSeeds, iterMax, verbose);
	
	// Repeat to refine fit:
	busyFitSuccess = FitBusyFunc(spectrum.size(), &channels[0], &spectrum[0], &noiseSpectrum[0], &busyFitParameters[0], fitCov, bestNoP, -1, iterMax, verbose);
	
	// Determine observational parameters:
	double **tmpPar = new double *[1];
	double **tmpObsVals = new double *[1];
	tmpPar[0] = new double[8];
	tmpObsVals[0] = new double[7];
	for(size_t i = 0; i < 8; ++i) tmpPar[0][i] = busyFitParameters[2*i];
	CalcObsParams(1, tmpPar, spectrum.size(), &channels[0], tmpObsVals, verbose);
	busyFunctionFint     = tmpObsVals[0][0];
	busyFunctionFpeak    = tmpObsVals[0][1];
	busyFunctionW50      = tmpObsVals[0][3];
	busyFunctionW20      = tmpObsVals[0][5];
	busyFunctionCentroid = 0.0; // ALERT: Russell's code does not calculate centroid!
	
	// Delete temporary variables again:
	delete[] tmpPar[0];
	delete[] tmpObsVals[0];
	delete[] tmpPar;
	delete[] tmpObsVals;
	
	// Delete the covariance matrix again
	// (we don't really need it anyway):
	for(int i = 0; i < 8; ++i) delete[] fitCov[i];
	delete[] fitCov;
	
	// Extract chi^2:
	busyFunctionChi2 = busyFitParameters[2 * BUSYFIT_FREE_PARAM];
	
	//BusyFit busyFit;
	//busyFit.setup(spectrum.size(), &spectrum[0], &noiseSpectrum[0], 2, true, false);
	
	//busyFitSuccess = busyFit.fit();
	
	//busyFit.getResult(&busyFitParameters[0], &busyFitUncertainties[0], busyFunctionChi2);
	//busyFit.getParameters(busyFunctionCentroid, busyFunctionW50, busyFunctionW20, busyFunctionFpeak, busyFunctionFint);
	
	// Correct spectral parameters for the shift caused by operating on a sub-cube:
	busyFitParameters[2] += subRegionZ1;
	busyFitParameters[4] += subRegionZ1;
	busyFitParameters[6] += subRegionZ1;
	//busyFunctionCentroid += subRegionZ1;
	
	return 0;
}

/*int Parametrization::fitBusyFunction()
{
	if(dataSize == 0 or spectrum.empty())
	{
		std::cerr << "Error (Parametrization): No valid data found.\n";
		return 1;
	}
	
	BusyFit busyFit;
	busyFit.setup(spectrum.size(), &spectrum[0], &noiseSpectrum[0], 2, true, false);
	
	busyFitSuccess = busyFit.fit();
	
	busyFit.getResult(&busyFitParameters[0], &busyFitUncertainties[0], busyFunctionChi2);
	busyFit.getParameters(busyFunctionCentroid, busyFunctionW50, busyFunctionW20, busyFunctionFpeak, busyFunctionFint);
	
	// Correct spectral parameters for the shift caused by operating on a sub-cube:
	busyFitParameters[4] += subRegionZ1;
	busyFitParameters[5] += subRegionZ1;
	busyFunctionCentroid += subRegionZ1;
	
	return 0;
}*/



// ------------------------ //
// Assign results to source //
// ------------------------ //

int Parametrization::writeParameters()
{
	source->setParameter("id",        source->getSourceID());
	source->setParameter("x",         centroidX);
	source->setParameter("y",         centroidY);
	source->setParameter("z",         centroidZ);
	source->setParameter("err_x",     errCentroidX);
	source->setParameter("err_y",     errCentroidY);
	source->setParameter("err_z",     errCentroidZ);
	source->setParameter("w50",       lineWidthW50);
	source->setParameter("w20",       lineWidthW20);
	source->setParameter("err_w50",   errlineWidthW50);
	source->setParameter("err_w20",   errlineWidthW20);
	source->setParameter("wm50",      lineWidthWm50);
	source->setParameter("f_wm50",    meanFluxWm50);
	source->setParameter("f_peak",    peakFlux);
	source->setParameter("f_int",     totalFlux);
	source->setParameter("snr_int",   intSNR);
	source->setParameter("ell_maj",   ellMaj);
	source->setParameter("ell_min",   ellMin);
	source->setParameter("ell_pa",    ellPA);
	source->setParameter("ell3s_maj", ell3sMaj);
	source->setParameter("ell3s_min", ell3sMin);
	source->setParameter("ell3s_pa",  ell3sPA);
	source->setParameter("kin_pa",    kinematicPA);
	//source->setParameter("flag_kin",  flagKinematicPA);
	//source->setParameter("flag_warp", flagWarp);
	
	source->setParameter("rms",       noiseSubCube);
	
	if(doBusyFunction)
	{
		source->setParameter("bf_flag",   busyFitSuccess);
		source->setParameter("bf_chi2",   busyFunctionChi2);
		source->setParameter("bf_a",      busyFitParameters[0]);
		source->setParameter("bf_b1",     busyFitParameters[1]);
		source->setParameter("bf_b2",     busyFitParameters[3]);
		source->setParameter("bf_c",      busyFitParameters[5]);
		source->setParameter("bf_xe",     (busyFitParameters[4] + busyFitParameters[1]) / 2.0);
		source->setParameter("bf_xp",     busyFitParameters[6]);
		source->setParameter("bf_w",      busyFitParameters[4] - busyFitParameters[2]);
		source->setParameter("bf_z",      busyFunctionCentroid);
		source->setParameter("bf_w20",    busyFunctionW20);
		source->setParameter("bf_w50",    busyFunctionW50);
		source->setParameter("bf_f_peak", busyFunctionFpeak);
		source->setParameter("bf_f_int",  busyFunctionFint);
	}
	
	return 0;
}
