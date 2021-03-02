#include <TomoLearnS/Filter.hpp>

#include <unsupported/Eigen/FFT>

#include <iostream>

Filter::Filter(FilterType filterType, double cutOffIn):filterType{filterType}{
	if( cutOffIn>0.0 && cutOffIn<=1.0){
		cutOff=cutOffIn;
	}
	else{
		std::cout << "WARNING: Incorrect cutOff value! Valid range: 0 < CutOff <=1.0 \n cutOff was set to 1.0 !!!!!!";
		cutOff=1.0;
	}
}

void Filter::operator()(Eigen::MatrixXcd& fftOfSinogram){

	switch (filterType){
		case FilterType::RamLak:
			RamLakFilter(fftOfSinogram);
			break;
		case FilterType::SheppLogan:
			//SheppLoganFilter(fftOfSinogram);
			break;
		case FilterType::Cosine:
			//CosineFilter(fftOfSinogram);
			break;
		case FilterType::Hamming:
			//HammingFilter(fftOfSinogram);
			break;
		case FilterType::Hanning:
			//HanningFilter(fftOfSinogram);
			break;
	}
}

void Filter::RamLakFilter(Eigen::MatrixXcd& fftOfSinogram){
	int pixNumPadded = fftOfSinogram.rows();

	//Construct the filter
	Eigen::MatrixXcd freqFilter = Eigen::MatrixXd::Zero(pixNumPadded,1);
	Eigen::MatrixXd filter = Eigen::MatrixXd::Zero(pixNumPadded,1);

	filter(0)=1/(4*0.107422*0.107422);
	for(int i=0; i<pixNumPadded/4; ++i){
		int idx= i*2+1;
		filter(idx) = -1/std::pow(idx*M_PI*0.107422, 2);
		filter(pixNumPadded-idx) = filter(idx);
	}
	Eigen::FFT<double> fft2;
	freqFilter.col(0)=fft2.fwd(filter.col(0));

	//Low-pass filter to suppress the noise
	/*int maxFreq=1024;
	for(int i=maxFreq+1; i<=pixNumPadded/2; ++i ){
		freqFilter(i)=0;
		freqFilter(pixNumPadded-i)=0;
	}*/

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

void SheppLoganFilter(Eigen::MatrixXcd& fftOfSinogram){return;}
void CosineFilter(Eigen::MatrixXcd& fftOfSinogram){return;}
void HammingFilter(Eigen::MatrixXcd& fftOfSinogram){return;}
void HanningFilter(Eigen::MatrixXcd& fftOfSinogram){return;}

