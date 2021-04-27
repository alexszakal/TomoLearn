#include <TomoLearnS/Filter.hpp>

#define EIGEN_FFTW_DEFAULT
#include <unsupported/Eigen/FFT>

#include <matplotlibcpp/matplotlibcpp_old.h>

#include <iostream>
#include <cmath>

Filter::Filter(FilterType filterType, double cutOffIn):filterType{filterType}{
	if( cutOffIn>0.0 && cutOffIn<=1.0){
		cutOff=cutOffIn;
	}
	else{
		std::cout << "WARNING: Incorrect cutOff value! Valid range: 0 < CutOff <=1.0 \n cutOff was set to 1.0 !!!!!!";
		cutOff=1.0;
	}
}

//#define FILTERING_DEBUG
void Filter::operator()(Eigen::MatrixXcd& fftOfSinogram){

	int pixNumPadded = fftOfSinogram.rows();

	//Construct the filter
	Eigen::MatrixXcd freqFilter = Eigen::MatrixXcd::Zero(pixNumPadded,1);
	Eigen::MatrixXd filter = Eigen::MatrixXd::Zero(pixNumPadded,1);

	filter(0)=1.0/4;
	for(int i=0; i<pixNumPadded/4; ++i){
		int idx= i*2+1;
		filter(idx) = -1.0/std::pow(idx*M_PI, 2);
		filter(pixNumPadded-idx) = filter(idx);
	}
	Eigen::FFT<double> fft2;
	freqFilter.col(0)=fft2.fwd(filter.col(0));

	switch (filterType){
		case FilterType::RamLak:
			RamLakFilter(freqFilter);
			break;
		case FilterType::SheppLogan:
			SheppLoganFilter(freqFilter);
			break;
		case FilterType::Cosine:
			CosineFilter(freqFilter);
			break;
		case FilterType::Hann:
			HannFilter(freqFilter);
			break;
		case FilterType::Hamming:
			HammingFilter(freqFilter);
			break;
	}

#ifdef FILTERING_DEBUG  //Show the filter
	matplotlibcpp::figure(3);
	Eigen::MatrixXd absVector = freqFilter.imag().array().pow(2) + freqFilter.real().array().pow(2) ;
	absVector = absVector.array().pow(0.5);
	std::cout << "\n H(0)= " << absVector(0);
	matplotlibcpp::plot(std::vector<float> (&absVector(0), absVector.data()+absVector.cols()*absVector.rows()) );
	matplotlibcpp::show();
#endif

	//Multiply with filter
	for(int i=0; i<fftOfSinogram.cols(); ++i){
		fftOfSinogram.col(i) = fftOfSinogram.col(i).array() * freqFilter.array();
	}
}

void Filter::RamLakFilter(Eigen::MatrixXcd& freqFilter){
	int pixNumPadded=freqFilter.rows();

	//Low-pass filter to suppress the noise
	int maxFreq=pixNumPadded*cutOff/2;
	for(int i=maxFreq+1; i<=pixNumPadded/2; ++i ){
		freqFilter(i)=0;
		freqFilter(pixNumPadded-i)=0;
	}
}

void Filter::SheppLoganFilter(Eigen::MatrixXcd& freqFilter){
	int pixNumPadded=freqFilter.rows();

	int maxFreq=pixNumPadded*cutOff/2;
	//freqFilter(0) *=1; //Shepp-Logan is 1 at x=0
	for (int i=1; i<=pixNumPadded/2; ++i){
		if(i>maxFreq){
			freqFilter(i)=0;
			freqFilter(pixNumPadded-i)=0;
		}
		else{
			double x = M_PI *i / (2 * maxFreq) ;
			freqFilter(i) *= sin( x ) / x;
			freqFilter(pixNumPadded-i)=freqFilter(i);
		}
	}
}

void Filter::CosineFilter(Eigen::MatrixXcd& freqFilter){
	int pixNumPadded=freqFilter.rows();

	int maxFreq=pixNumPadded*cutOff/2;
	//freqFilter(0) *=1; //Cosine is 1 at x=0
	for (int i=1; i<=pixNumPadded/2; ++i){
		if(i>maxFreq){
			freqFilter(i)=0;
			freqFilter(pixNumPadded-i)=0;
		}
		else{
			double x = M_PI *i / (2 * maxFreq) ;
			freqFilter(i) *= cos( x );
			freqFilter(pixNumPadded-i)=freqFilter(i);
		}
	}
	return;
}


void Filter::HannFilter(Eigen::MatrixXcd& freqFilter){
	int pixNumPadded=freqFilter.rows();

	int maxFreq=pixNumPadded*cutOff/2;
	//freqFilter(0) *=1; //Hann is 1 at x=0
	for (int i=1; i<=pixNumPadded/2; ++i){
		if(i>maxFreq){
			freqFilter(i)=0;
			freqFilter(pixNumPadded-i)=0;
		}
		else{
			double x = M_PI *i / ( maxFreq) ;
			freqFilter(i) *= (1+cos( x ))*0.5;
			freqFilter(pixNumPadded-i)=freqFilter(i);
		}
	}
	return;
}

void Filter::HammingFilter(Eigen::MatrixXcd& freqFilter){
	int pixNumPadded=freqFilter.rows();
	double alpha = 25.0/46;

	int maxFreq=pixNumPadded*cutOff/2;
	//freqFilter(0) *=1; //HAmming is 1 at x=0
	for (int i=1; i<=pixNumPadded/2; ++i){
		if(i>maxFreq){
			freqFilter(i)=0;
			freqFilter(pixNumPadded-i)=0;
		}
		else{
			double x = M_PI *i / (1 * maxFreq) ;
			freqFilter(i) *= alpha+(1-alpha)*cos( x );
			freqFilter(pixNumPadded-i)=freqFilter(i);
		}
	}
	return;
}

