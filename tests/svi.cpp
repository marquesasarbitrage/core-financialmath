#include <iostream>
#include <cassert>
#include <map>
#include "../include/core-financialmath/svi.hpp"

bool isClose(double a, double b, double tol = 1e-3) {

    return std::abs(a - b) <= tol;
}

void testRawSVIConstructor1() {

    SVI rawSVI(std::make_tuple(-0.0410,0.1331,0.3586,1.0,0.3060), 1.0);

    assert(isClose(rawSVI.getTotalVariance(-1.5), 0.1642151669, 1e-8));
    assert(isClose(rawSVI.getTotalVariance(.4), 0.09390017924, 1e-8));
    assert(isClose(rawSVI.getImpliedVariance(-.8), 0.1155180767, 1e-8));
    assert(isClose(rawSVI.getImpliedVariance(.8), 0.1224671653, 1e-8));
    assert(isClose(rawSVI.getImpliedVolatility(-.2), 0.2978374318, 1e-8));
    assert(isClose(rawSVI.getImpliedVolatility(.7), 0.3369683615, 1e-8));

    assert(isClose(rawSVI.getTotalVarianceATM(), 0.08579391231, 1e-8));
    assert(isClose(rawSVI.getImpliedVarianceATM(), 0.08579391231, 1e-8));
    assert(isClose(rawSVI.getImpliedVolatilityATM(), 0.2929059786, 1e-8));

    assert(isClose(rawSVI.getRiskNeutralDensity(-1.5), 0.4212048496, 1e-8));
    assert(isClose(rawSVI.getRiskNeutralDensity(0.0), 1.055453702, 1e-8));
    assert(isClose(rawSVI.getRiskNeutralDensity(0.9), 0.4417043688, 1e-8));

    assert(isClose(rawSVI.getLocalVariance(-1.5), 0.291713995, 1e-1));
    assert(isClose(rawSVI.getLocalVariance(0.0), 0.08128628681, 1e-1));
    assert(isClose(rawSVI.getLocalVariance(0.9), 0.2198565415, 1e-1));
    assert(isClose(rawSVI.getLocalVariance(1.5), 0.5020725112, 1e-1));

    assert(isClose(rawSVI.getLocalVolatility(-1.5), 0.5401055406, 1e-1));
    assert(isClose(rawSVI.getLocalVolatility(0.0), 0.2851075005, 1e-1));
    assert(isClose(rawSVI.getLocalVolatility(0.9), 0.4688886238, 1e-1));
    assert(isClose(rawSVI.getLocalVolatility(1.5), 0.7085707525, 1e-1));

    assert(rawSVI.getTimeToMaturity()==1.0);

    assert(isClose(rawSVI.getCallSlope(),0.5934621096, 1e-8));
    assert(isClose(rawSVI.getPutSlope(),0.315361948, 1e-8));
    assert(isClose(rawSVI.getMinimumVariance(),0.08571539426, 1e-8));
    assert(isClose(rawSVI.getSkewATM(),-0.007168934643, 1e-8));
    //assert(rawSVI.getRawParameters() == std::make_tuple(-0.0410,0.1331,0.3586,0.4153,0.3060));

    assert(rawSVI.isFreeOfButterflyArbitrage()==true);


}

void testRawSVIConstructor2() {

    SVI rawSVI(std::make_tuple(-0.0410,0.1331,0.3586,1.0,0.3060), .8);

    assert(isClose(rawSVI.getTotalVariance(-1.5), 0.1642151669, 1e-8));
    assert(isClose(rawSVI.getTotalVariance(.4), 0.09390017924, 1e-8));
    assert(isClose(rawSVI.getImpliedVariance(-.8), 0.1443975959, 1e-8));
    assert(isClose(rawSVI.getImpliedVariance(.8), 0.1530839567, 1e-8));
    assert(isClose(rawSVI.getImpliedVolatility(-.2), 0.3329923719, 1e-8));
    assert(isClose(rawSVI.getImpliedVolatility(.7), 0.3767420813, 1e-8));

    assert(isClose(rawSVI.getTotalVarianceATM(), 0.08579391231, 1e-8));
    assert(isClose(rawSVI.getImpliedVarianceATM(), 0.1072423904, 1e-8));
    assert(isClose(rawSVI.getImpliedVolatilityATM(), 0.3274788396, 1e-8));

    assert(isClose(rawSVI.getRiskNeutralDensity(-1.5), 0.4212048496, 1e-8));
    assert(isClose(rawSVI.getRiskNeutralDensity(0.0), 1.055453702, 1e-8));
    assert(isClose(rawSVI.getRiskNeutralDensity(0.9), 0.4417043688, 1e-8));

    assert(isClose(rawSVI.getLocalVariance(0.0), 0.1016078585, 1e-1));
    assert(isClose(rawSVI.getLocalVariance(0.2), 0.1049961905, 1e-1));
    assert(isClose(rawSVI.getLocalVariance(0.4), 0.1237613799, 1e-1));
    assert(isClose(rawSVI.getLocalVariance(0.9), 0.2748225508, 1e-1));

    assert(isClose(rawSVI.getLocalVolatility(0.0), 0.3187598759, 1e-1));
    assert(isClose(rawSVI.getLocalVolatility(0.2), 0.3240311567, 1e-1));
    assert(isClose(rawSVI.getLocalVolatility(0.4), 0.3517973563, 1e-1));
    assert(isClose(rawSVI.getLocalVolatility(0.9), 0.5242352056, 1e-1));

    assert(rawSVI.getTimeToMaturity()==.8);

    assert(isClose(rawSVI.getCallSlope(),0.5934621096, 1e-8));
    assert(isClose(rawSVI.getPutSlope(),0.315361948, 1e-8));
    assert(isClose(rawSVI.getMinimumVariance(),0.1071442428, 1e-8));
    assert(isClose(rawSVI.getSkewATM(),-0.007168934643, 1e-8));

    assert(rawSVI.isFreeOfButterflyArbitrage()==true);

}

void testRawSVIAxelVogt() {

    SVI rawSVI(std::make_tuple(-0.0410,0.1331,0.3586,0.4153,0.3060), 1.0);

    assert(isClose(rawSVI.getCallSlope(),1.316798219, 1e-8));
    assert(isClose(rawSVI.getPutSlope(),0.6997381041, 1e-8));
    assert(isClose(rawSVI.getMinimumVariance(),0.01162490324, 1e-8));
    assert(isClose(rawSVI.getSkewATM(),-0.1752111408, 1e-8));

    assert(rawSVI.isFreeOfButterflyArbitrage()==false);

    assert(isClose(rawSVI.getTotalVariance(-1.5), 0.1367819611, 1e-8));
    assert(isClose(rawSVI.getTotalVariance(.4), 0.01623656962, 1e-8));
    assert(isClose(rawSVI.getImpliedVariance(-.8), 0.0756291293, 1e-8));
    assert(isClose(rawSVI.getImpliedVariance(.8), 0.05764411607, 1e-8));
    assert(isClose(rawSVI.getImpliedVolatility(-.2), 0.1699864855, 1e-8));
    assert(isClose(rawSVI.getImpliedVolatility(.7), 0.210857945, 1e-8));

    assert(isClose(rawSVI.getTotalVarianceATM(), 0.01742625256, 1e-8));
    assert(isClose(rawSVI.getImpliedVarianceATM(), 0.01742625256, 1e-8));
    assert(isClose(rawSVI.getImpliedVolatilityATM(), 0.1320085321, 1e-8));

    assert(isClose(rawSVI.getRiskNeutralDensity(-1.5), 0.2478306953, 1e-8));
    assert(isClose(rawSVI.getRiskNeutralDensity(0.0), 1.038649731, 1e-8));
    assert(isClose(rawSVI.getRiskNeutralDensity(0.9), -0.03268513071, 1e-8));

    assert(isClose(rawSVI.getLocalVariance(0.0), 0.01677779527, 1e-2));
    assert(rawSVI.getLocalVariance(0.7)<0.0);
    assert(rawSVI.getLocalVariance(0.9)<0.0);

    assert(std::isnan(rawSVI.getLocalVolatility(0.7)));
    assert(std::isnan(rawSVI.getLocalVolatility(0.9)));

    assert(rawSVI.getTimeToMaturity()==1.0);

}

void testJWSVIConstructor() {

    SVI jwSVI(0.01742625,-0.1752111,0.0116249,0.8564763,0.6997381,1.0); 

    assert(isClose(jwSVI.getTotalVariance(-1.5), 0.1380576036, 1e-6));
    assert(isClose(jwSVI.getTotalVariance(.4), 0.01514407736, 1e-6));

    assert(isClose(jwSVI.getImpliedVariance(-.8), 0.07639878745, 1e-6));
    assert(isClose(jwSVI.getImpliedVariance(.8), 0.04372812104, 1e-6));

    assert(isClose(jwSVI.getImpliedVolatility(-.2), 0.1702976916, 1e-6));
    assert(isClose(jwSVI.getImpliedVolatility(.7), 0.186893459, 1e-6));

    assert(isClose(jwSVI.getTotalVarianceATM(), 0.01742633055, 1e-6));
    assert(isClose(jwSVI.getImpliedVarianceATM(), 0.01742633055, 1e-6));
    assert(isClose(jwSVI.getImpliedVolatilityATM(), 0.1320088276, 1e-6));

    assert(isClose(jwSVI.getRiskNeutralDensity(-1.5), 0.2492476732, 1e-6));
    assert(isClose(jwSVI.getRiskNeutralDensity(0.0), 1.041529529, 1e-6));
    assert(isClose(jwSVI.getRiskNeutralDensity(0.9), 0.01074175286, 1e-6));

    assert(jwSVI.getTimeToMaturity()==1.0);

    assert(isClose(jwSVI.getCallSlope(),0.8564763, 1e-8));
    assert(isClose(jwSVI.getPutSlope(),0.6997381, 1e-8));
    assert(isClose(jwSVI.getMinimumVariance(),0.0116249, 1e-8));
    assert(isClose(jwSVI.getSkewATM(),-0.1752111, 1e-8));

    assert(jwSVI.isFreeOfButterflyArbitrage()==true);

}

int main() {

    testRawSVIConstructor1();
    testRawSVIConstructor2();
    testRawSVIAxelVogt();
    testJWSVIConstructor();
    std::cout << "All tests for the SVI object have been passed successfully!" << std::endl; 
    return 0;
}