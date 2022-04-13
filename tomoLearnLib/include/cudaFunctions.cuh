#pragma once

void launchRayDrivenKernel(const double* phantomData, std::array<int, 2> numberOfPixels, std::array<double, 2> pixSizes,
		                   int numAngles, const double* anglesData, int pixNum, double detWidth,
						   double* sinogramData);
