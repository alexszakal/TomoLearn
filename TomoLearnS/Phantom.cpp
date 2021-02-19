#include <TomoLearnS/Phantom.hpp>
#include <iostream>

//Phantom::Phantom():Object2D(){
//
//}

Phantom::Phantom(const std::string& label,
		         const std::string& imageFilePath,
				 const std::array<double, 2>& objPixSizes):
				                                          Object2D{imageFilePath, objPixSizes},
														  label{label}{

}

Phantom::~Phantom(){
	if(!cimg_window.is_closed()){
	    cimg_window.close();
	    cimg_window.flush();
	}
	if(displayThread.joinable())
		displayThread.join();
}

void Phantom::display(){
	if(!(cimg_window.is_closed())){
		return;
	}
	if(displayThread.joinable()){
		displayThread.join();
	}

	std::array<int, 2> numberOfPixels = getNumberOfPixels();
	std::array<double, 2> objPixSizes = getPixSizes();

	std::string title{"TEST"};

	std::stringstream ss;
	ss << title << " Label: " << label << " " << numberOfPixels[0] << " x " << numberOfPixels[1] << " points; "
			                                     << objPixSizes[0] << " x " << objPixSizes[1] << " pixSize"
									             << numberOfPixels[0]*objPixSizes[0] << " x "
												 << numberOfPixels[1]*objPixSizes[1] << " w x h";
	std::string longTitle = ss.str();

	cimg_image.assign(numberOfPixels[0], numberOfPixels[1], 1, 1);

	const Eigen::MatrixXd& objData = getDataAsEigenMatrixRef();

	double maxInt = objData.maxCoeff();
	double minInt = objData.minCoeff();
	double normFactor = 65530/(maxInt-minInt);
	for(int i=0; i<numberOfPixels[0]; ++i){
		for(int j=0; j<numberOfPixels[1]; ++j){
			cimg_image(i,j) = static_cast<uint16_t>((objData(i,j)-minInt)*normFactor);
		}
	}

	//cimg_image.display(cimg_window, true,0, true);
	displayThread = std::thread([this](){cimg_image.display(cimg_window, true,0, true);});
	while( cimg_window.is_closed() ){
		cimg_window.wait(10); //Wait 10 milliseconds
	}
	cimg_window.set_title("%s", longTitle.c_str());
}




