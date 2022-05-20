#include <iostream>
#include <cmath>
#include <chrono>
#include <algorithm>

#include <CImg.h>

#ifdef Success       //Because otherwise Eigen not compile (EIgen <-> CImg interference)
  #undef Success
#endif

#include <Object2D.hpp>
#include <Gen1CT.hpp>

#include <config.h>

void testRadonTransform(const std::string& phantomName, projectorType projectAlgo);

int main(){

#if ENABLE_CUDA
	std::cout << "\n \n CUDA enabled!!!!" ;
#else
	std::cout << "\n \n CUDA disabled!!!" ;
#endif

	testRadonTransform("modSL", projectorType::rayDriven );

	return 0;
}

/***
 * Compare the numerical and analytic Radon transform of an ellipse
 * @param phantomName Name of the phantom
 * @param projectAlgo Type of projector used
 */
void testRadonTransform(const std::string& phantomName, projectorType projectAlgo){

	std::cout << "Parallel beam projection simulation" << std::endl;

	double PhantomWidthInMM = 102.4;

	std::vector<double> rhos;
	std::vector<double> alphas;
	std::vector<std::array<double,2>> centers;
    std::vector<std::array<double,2>> axes;

	//Parameters of the Shepp-Logan phantom
    if(phantomName == "SL"){
    	rhos = { 1.0, -0.98,
		        -0.02, -0.02,
		         0.01, 0.01,
				 0.01, 0.01,
				 0.01, 0.01 };

    	alphas = { 0.0/180.0*M_PI, 0.0/180.0*M_PI,
	              -18.0/180.0*M_PI, 18.0/180.0*M_PI,
				   0.0/180.0*M_PI, 0.0/180.0*M_PI,
				   0.0/180.0*M_PI, 0.0/180.0*M_PI,
				   0.0/180.0*M_PI, 0.0/180.0*M_PI };

    	centers = { {0.0*PhantomWidthInMM/2, 0.0*PhantomWidthInMM/2}, {0.0*PhantomWidthInMM/2, -0.0184*PhantomWidthInMM/2},
	                {0.22*PhantomWidthInMM/2, 0.0*PhantomWidthInMM/2}, {-0.22*PhantomWidthInMM/2, 0.0*PhantomWidthInMM/2},
	                {0.0*PhantomWidthInMM/2, 0.35*PhantomWidthInMM/2}, {0.0*PhantomWidthInMM/2, 0.1*PhantomWidthInMM/2},
					{0.0*PhantomWidthInMM/2, -0.1*PhantomWidthInMM/2}, {-0.08*PhantomWidthInMM/2, -0.605*PhantomWidthInMM/2},
					{0*PhantomWidthInMM/2, -0.606*PhantomWidthInMM/2}, {0.06*PhantomWidthInMM/2, -0.605*PhantomWidthInMM/2} };
    	axes = { {0.69*PhantomWidthInMM/2, 0.92*PhantomWidthInMM/2}, {0.6624*PhantomWidthInMM/2, 0.874*PhantomWidthInMM/2},
			     {0.11*PhantomWidthInMM/2, 0.31*PhantomWidthInMM/2}, {0.16*PhantomWidthInMM/2, 0.41*PhantomWidthInMM/2},
				 {0.21*PhantomWidthInMM/2, 0.25*PhantomWidthInMM/2}, {0.046*PhantomWidthInMM/2, 0.046*PhantomWidthInMM/2},
				 {0.046*PhantomWidthInMM/2, 0.046*PhantomWidthInMM/2}, {0.046*PhantomWidthInMM/2, 0.023*PhantomWidthInMM/2},
				 {0.023*PhantomWidthInMM/2, 0.023*PhantomWidthInMM/2}, {0.023*PhantomWidthInMM/2, 0.046*PhantomWidthInMM/2} };
    } else if(phantomName=="modSL"){
    	rhos = { 2.0, -0.8,
    		    -0.2, -0.2,
    		     0.1, 0.1,
    			 0.1, 0.1,
    			 0.1, 0.1 };

        alphas = { 0.0/180.0*M_PI, 0.0/180.0*M_PI,
    	         -18.0/180.0*M_PI, 18.0/180.0*M_PI,
    	     	   0.0/180.0*M_PI, 0.0/180.0*M_PI,
    			   0.0/180.0*M_PI, 0.0/180.0*M_PI,
    			   0.0/180.0*M_PI, 0.0/180.0*M_PI };

        centers = { {0.0*PhantomWidthInMM/2, 0.0*PhantomWidthInMM/2}, {0.0*PhantomWidthInMM/2, -0.0184*PhantomWidthInMM/2},
    	            {0.22*PhantomWidthInMM/2, 0.0*PhantomWidthInMM/2}, {-0.22*PhantomWidthInMM/2, 0.0*PhantomWidthInMM/2},
    	            {0.0*PhantomWidthInMM/2, 0.35*PhantomWidthInMM/2}, {0.0*PhantomWidthInMM/2, 0.1*PhantomWidthInMM/2},
    	         	{0.0*PhantomWidthInMM/2, -0.1*PhantomWidthInMM/2}, {-0.08*PhantomWidthInMM/2, -0.605*PhantomWidthInMM/2},
    				{0*PhantomWidthInMM/2, -0.606*PhantomWidthInMM/2}, {0.06*PhantomWidthInMM/2, -0.605*PhantomWidthInMM/2} };
        axes = { {0.69*PhantomWidthInMM/2, 0.92*PhantomWidthInMM/2}, {0.6624*PhantomWidthInMM/2, 0.874*PhantomWidthInMM/2},
    			 {0.11*PhantomWidthInMM/2, 0.31*PhantomWidthInMM/2}, {0.16*PhantomWidthInMM/2, 0.41*PhantomWidthInMM/2},
    			 {0.21*PhantomWidthInMM/2, 0.25*PhantomWidthInMM/2}, {0.046*PhantomWidthInMM/2, 0.046*PhantomWidthInMM/2},
    			 {0.046*PhantomWidthInMM/2, 0.046*PhantomWidthInMM/2}, {0.046*PhantomWidthInMM/2, 0.023*PhantomWidthInMM/2},
    			 {0.023*PhantomWidthInMM/2, 0.023*PhantomWidthInMM/2}, {0.023*PhantomWidthInMM/2, 0.046*PhantomWidthInMM/2} };
    }
    else{
    	std::cout << "\nphantomName parameter not recognized. Possible values: \"SL\" or \"modSL\" ";
    	std::cout << "\nAborting testRadonTransform function";
    	return;
    }

    //Generate phantom
	Phantom ellipsePhantom("ellipsePhantom",
			               {1024, 1024},
						   {0.1, 0.1},
						   rhos, //rhos
						   alphas, //alphas
						   centers,  //centers
						   axes //axes
                           );
	ellipsePhantom.display();

	Phantom centerDotPhantom( "centerDotPhantom",
	            {1024, 1024},
				   {0.1, 0.1},
				   std::vector<double>{1000}, //rhos
				   std::vector<double>{0.0}, //alphas
				   std::vector<std::array<double,2>> {{0,0}},  //centers
				   std::vector<std::array<double,2>> {{10,10}} //axes
	            );
	//centerDotPhantom.display();

	//generate Radon Transform Numerically
	int detWidthInMM { 110 };
	int detPixNum { 512 };
	Gen1CT ct(detWidthInMM, detPixNum);
	ct.setI0(0.0);

	const int numProjections{180};
	Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(numProjections, 0.0/180.0 * M_PI,
			179.0 / 180 * M_PI);

	ct.addPhantom(ellipsePhantom);
	ct.addPhantom(centerDotPhantom);

	ct.measure("ellipsePhantom", angles, "Sinogram", projectAlgo);

	ct.displayMeasurement("Sinogram");

	CTScan numericalSinogram = ct.getMeasurement("Sinogram");

	CTScan analyticSinogram("analyticCenterDotSinogram",
							detWidthInMM,
							detPixNum,
							angles,
							rhos, //rhos
							alphas, //alphas
							centers,  //centers
							axes, //axes
							0.0); //I0=0 -> do not draw Poisson statistics

	analyticSinogram.display("AnalyticResult");

	const Eigen::MatrixXd& numRes = numericalSinogram.getDataAsEigenMatrixRef();
	const Eigen::MatrixXd& anRes = analyticSinogram.getDataAsEigenMatrixRef();

	const Eigen::MatrixXd relativeError((numRes.array()-anRes.array())/((numRes.array() + anRes.array()+0.000001) * 0.5)*100); // @suppress("Invalid arguments")

	CTScan metric("metric", relativeError.cwiseAbs(), detWidthInMM, angles, 0.0);
	metric.display();

	std::cout << "\nDifference was normalized with the average of the corresponding pixel values.";
	std::cout << "\nMaximal relative error: " << relativeError.cwiseAbs().maxCoeff() << "%";
	std::cout << "\nAverage error: " << relativeError.cwiseAbs().mean() << "%";

	metric.display();

	std::cout<<"\n Press ENTER to continue";
	std::cin.get();
}


