/*
 * Phantom.hpp
 *
 *  Created on: 2021. febr. 16.
 *      Author: szakal
 */

#ifndef PHANTOM_HPP_
#define PHANTOM_HPP_

#include <CImg.h>

#include <TomoLearnS/Object2D.hpp>
#include <string>
#include <thread>

class Phantom : public Object2D{
public:
	Phantom(const std::string& label, const std::string& imageFilePath,
			const std::array<double, 2>& objPixSizes);

private:
	std::string label;


};




#endif /* PHANTOM_HPP_ */
