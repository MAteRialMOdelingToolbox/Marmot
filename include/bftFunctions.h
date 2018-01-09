#pragma once 
#include <string>

//external!
extern "C" bool warningToMSG(const std::string& message); // must return a 'false';
extern "C" bool notificationToMSG(const std::string& message); // must return a 'true';

namespace bft{
    namespace Functions
    {
        double linearInterpolation(double x, double x0, double x1, double y0, double y1);
        double exp(double x);
        double getExponentPowerTen(const double x);
    	double radToDeg(const double alpha);
		double degToRad(const double alpha);
    }
}

