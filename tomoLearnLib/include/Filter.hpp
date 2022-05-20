#pragma once
#include <Eigen/Dense>

enum class FilterType{RamLak,
                      SheppLogan,
                      Cosine,
                      Hann,
                      Hamming};
/**
 * @brief Filter class implements different filters that are used for the filtered backprojection algorithm
 */
class Filter{
public:
	Filter(FilterType filterType, double cutOffIn = 1.0);

	void operator()(Eigen::MatrixXcd& fftOfSinogram);
private:
	FilterType filterType;  /** Type of filter */
	double cutOff;  /** Cut-off frequency of the filter */

	void RamLakFilter    (Eigen::MatrixXcd& freqFilter);
	void SheppLoganFilter(Eigen::MatrixXcd& freqFilter);
	void CosineFilter    (Eigen::MatrixXcd& freqFilter);
	void HannFilter   (Eigen::MatrixXcd& freqFilter);
	void HammingFilter   (Eigen::MatrixXcd& freqFilter);
};
