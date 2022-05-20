#include <vector>
#include <Reconst.hpp>

/***
 * Returns the convergence curve
 * @return Returns the convergence curve or an empty vector if it is not available
 */
std::vector<double> Reconst::getConvergenceCurve() const{
	return convergenceCurve;
}
