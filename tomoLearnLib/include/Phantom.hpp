/*
 * Phantom.hpp
 *
 *  Created on: 2021. febr. 16.
 *      Author: szakal
 */

#ifndef PHANTOM_HPP_
#define PHANTOM_HPP_

#include <CImg.h>

#include <Object2D.hpp>
#include <string>
#include <thread>

class Phantom : public Object2D{
public:
	Phantom();
	Phantom(const std::string& label, const std::string& imageFilePath,
			const std::array<double, 2>& objPixSizes,
			bool convertFromHUtoLA=false);
	Phantom(const std::string& label, const Eigen::MatrixXd inData,
			const std::array<double, 2> objPixSizes = {0.1, 0.1});
	Phantom(std::string label, const Object2D& dataPar);
	Phantom(std::string label,
			const std::array<int,2>& numberOfPixels,
			const std::array<double,2>& objPixSizes,
			const std::vector<double>& rhos,
			const std::vector<double>& alphas,
			const std::vector< std::array<double,2> >& centers,
			const std::vector< std::array<double,2> >& axes );

	Phantom operator*(double coeff) const;
	Phantom operator+(double addVal) const;
	Phantom operator-(double subVal) const;

	std::string getLabel() const;

	std::array<Phantom,2> calculateQuadRegTerms() const;
	std::array<Phantom,2> calculateHuberRegTerms(double delta) const;
	std::array<Phantom,2> calculateGibbsRegTerms(double delta) const;

	friend Phantom operator/(const Phantom& lhs, const Phantom& rhs);
	friend Phantom operator*(const Phantom& lhs, const Phantom& rhs);
	friend Phantom operator+(const Phantom& lhs, const Phantom& rhs);

private:
	std::string label;
};




#endif /* PHANTOM_HPP_ */
