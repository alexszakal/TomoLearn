#include <iostream>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <string>

#include <CImg.h>

#ifdef Success       //Because otherwise Eigen not compile (EIgen <-> CImg interference)
  #undef Success
#endif

#include <Object2D.hpp>
#include <Gen1CT.hpp>

#include <matplotlibcpp.h>

#include <config.h>

void testHaoGaoTransform_CPU( const std::string& phantomName, bool useGPU );

//TODO: A ct.compareRowPhantomAndReconst() Mukodjon. HA fajlbol olvasunk, akkor 1000-et ki kell vonni, mert akkor kapjuk meg HU unitban!

//DEBUG: A szurt szinogram eltunik amikor a visszaallitas megjelenik
//TODO: Valahogy a szurt szinogramokat is el kell menteni (lehetne egy map, ahol a key a filter osztaly?? )

//TODO: Visszavetitest felgyorsitani
//TODO: Gyorsabb elorevetites a cache jobb hasznalataval

//TTOK sanitizers:
//-Dokumentacioba Properties -> C/C++Build -> CMake4eclipse -> Symbols -> Cmake cache entries -> ENABLE_SANITIZER_ADDRESS:BOOL:ON
//- Az address symbolizer is kell ahhoz hogy kodsorokat irjon ki: LAunch config. -> Environmentbe:
                                                           // -> ASAN_SYMBOLIZER_PATH=/usr/lib/llvm-6.0/bin/llvm-symbolizer

//Parallel geometry
int main(){

	//testHaoGaoTransform_CPU( "SD" );
	testHaoGaoTransform_CPU( "modSL", true );

	std::cin.ignore();

	return 0;
}

/***
 * Compare the numerical and analytic X-ray transform of an ellipse
 * the ray-driven method developed by Hao Gao is used
 * @param phantomName "SL" Shepp-Logan or "modSL" modified Shepp-Logan phantoms are available
 * @param useGPU Use the GPU acceleration
 */
void testHaoGaoTransform_CPU(const std::string& phantomName, bool useGPU){

	std::cout << "Parallel beam projection simulation using the method proposed by Hao Gao" << std::endl;

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
    	std::cout << "\nphantomName is not \"SL\" or \"modSL\". Using one of the pixelated phantoms ";
    }

	//generate Radon Transform Numerically
	int detWidthInMM { 110 };
	int detPixNum { 512 };
	Gen1CT ct(detWidthInMM, detPixNum);
	ct.setI0(0.0);

	const int numProjections{180};
	Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(numProjections, 0.0/180.0 * M_PI,
			179.0 / 180 * M_PI);

    //Add phantoms
    if( (phantomName == "SL") or (phantomName == "modSL")){
    	Phantom ellipsePhantom(phantomName,
    			               {1024, 1024},
						       {0.1, 0.1},
						       rhos, //rhos
						       alphas, //alphas
						       centers,  //centers
						       axes //axes
                               );
    	ct.addPhantom(ellipsePhantom);
    }
	ct.addPhantom("SD", "Phantoms/SingleDot.png"); //Single dot Phantom

	ct.displayPhantom(phantomName);

	if(useGPU){
#if ENABLE_CUDA
		ct.measure(phantomName, angles, "HaoGaoSinogram", projectorType::rayDriven_GPU);
#else
		std::cout << "\nCUDA is not allowed, fallback to CPU projector!!";
		ct.measure(phantomName, angles, "HaoGaoSinogram", projectorType::rayDriven);
#endif
	}
	else{
		ct.measure(phantomName, angles, "HaoGaoSinogram", projectorType::rayDriven);
	}

	ct.displayMeasurement("HaoGaoSinogram");

	CTScan numericalSinogram = ct.getMeasurement("HaoGaoSinogram");

	if( (phantomName == "SL") or (phantomName == "modSL")){

		CTScan analyticSinogram("analyticSinogram",
				                detWidthInMM,
				                detPixNum,
				                angles,
				                rhos, //rhos
				                alphas, //alphas
				                centers,  //centers
				                axes, //axes
                                0.0); //I0
		analyticSinogram.display("AnalyticResult");

		const Eigen::MatrixXd& numRes = numericalSinogram.getDataAsEigenMatrixRef();
		const Eigen::MatrixXd& anRes = analyticSinogram.getDataAsEigenMatrixRef();

		const Eigen::MatrixXd relativeError((numRes.array()-anRes.array())/((numRes.array() + anRes.array()+0.000001) * 0.5)*100); // @suppress("Invalid arguments")

		CTScan metric("metric", relativeError.cwiseAbs(), detWidthInMM, angles, 0.0);
		metric.display();

		std::cout << "\nDifference was normalized with the average of the corresponding pixel values.";
		std::cout << "\nMaximal relative error: " << relativeError.cwiseAbs().maxCoeff() << "%";
		std::cout << "\nAverage error: " << relativeError.cwiseAbs().mean() << "%";
	}

	//int tmpi;
	//std::cin>>tmpi;
}
