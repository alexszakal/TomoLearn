#pragma once

void launchRayDrivenProjectionKernel(const double* phantomData, std::array<int, 2> numberOfPixels, std::array<double, 2> pixSizes,
		                   int numAngles, const double* anglesData, int pixNum, double detWidth,
						   double* sinogramData);

void launchRayDrivenBackprojectionKernel(const double* sinogram, int numAngles, const double* anglesData, int pixNum, double detWidth,
		                                 double* detPixCenters,
										 double* backProjection, const std::array<int,2>& numberOfPixels, const std::array<double,2>& resolution);
