#include <iostream>
#include <cassert>
#include <map>
#include "../include/core-financialmath/baroneadesiwhaley.hpp"

bool isClose(double a, double b, double tol = 1e-3) {
    return std::abs(a - b) <= tol;
}

void testTemplateBAW(double mu, double r,double sigma, double T, bool isCall, std::vector<double> baw_prices)
{
    std::map<double,double> baw_prices_map =  { 
        {80, baw_prices[0]},
        {90, baw_prices[1]},
        {100, baw_prices[2]},
        {110, baw_prices[3]},
        {120, baw_prices[4]}
    };
    std::vector<double> spot_list = {80, 90,100,110,120};
    double price, F;
    for (double s: spot_list){
        price = BaroneAdesiWhaley::getPrice(s,mu,100,T,sigma,r,isCall);
        assert(isClose(price, baw_prices_map.at(s), 1e-2));
    }
};

void testBaroneAdesiWhalley() {

    // Call and mu = 0.04
    std::vector<double> baw_prices = {.05,.85,4.44,11.66,20.90};
    testTemplateBAW(.04,.08,.2,.25,true,baw_prices);
    baw_prices =  {.05,.84,4.40,11.55,20.69};
    testTemplateBAW(.04,.12,.2,.25,true,baw_prices);
    baw_prices =  {1.29,3.82,8.35,14.8,22.72};
    testTemplateBAW(.04,.08,.4,.25,true,baw_prices);
    baw_prices =  {.41,2.18,6.5,13.42,22.06};
    testTemplateBAW(.04,.08,.2,.5,true,baw_prices);

    // Call and mu = 0.00
    baw_prices =  {.04,.70,3.93,10.81,20.02};
    testTemplateBAW(.00,.08,.2,.25,true,baw_prices);
    baw_prices =  {.04,.70,3.9,10.75,20.0};
    testTemplateBAW(.00,.12,.2,.25,true,baw_prices);
    baw_prices =  {1.17,3.53,7.84,14.08,21.86};
    testTemplateBAW(.00,.08,.4,.25,true,baw_prices);
    baw_prices =  {.30,1.72,5.48,11.9,20.34};
    testTemplateBAW(.00,.08,.2,.5,true,baw_prices);

    // Call and T = 3.0
    baw_prices =  {4.2,7.54,12.03,17.64,24.30};
    testTemplateBAW(.00,.08,.2,3.0,true,baw_prices);
    baw_prices =  {6.97,11.62,17.4,24.09,31.49};
    testTemplateBAW(.04,.08,.2,3.0,true,baw_prices);
    baw_prices =  {2.52,4.97,8.67,13.88,20.88};
    testTemplateBAW(-0.04,.08,.2,3.0,true,baw_prices);

    // Call and mu = -0.04
    baw_prices =  {.03,.59,3.52,10.31,20.0};
    testTemplateBAW(-.04,.08,.2,.25,true,baw_prices);
    baw_prices =  {.03,.59,3.51,10.29,20.0};
    testTemplateBAW(-.04,.12,.2,.25,true,baw_prices);
    baw_prices =  {1.07,3.28,7.41,13.5,21.23};
    testTemplateBAW(-.04,.08,.4,.25,true,baw_prices);
    baw_prices =  {.23,1.39,4.72,10.96,20.0};
    testTemplateBAW(-.04,.08,.2,.5,true,baw_prices);

    // Put and mu = 0.04
    baw_prices =  {20.0,10.58,3.93,0.94,0.15};
    testTemplateBAW(.00,.08,.2,.25,false,baw_prices);
    baw_prices =  {20.0,10.53,3.9,0.93,0.15};
    testTemplateBAW(.00,.12,.2,.25,false,baw_prices);
    baw_prices =  {20.93,13.39,7.84,4.23,2.12};
    testTemplateBAW(.00,.08,.4,.25,false,baw_prices);
    baw_prices =  {20.04,11.48,5.48,2.15,0.7};
    testTemplateBAW(.00,.08,.2,.5,false,baw_prices);

    // Put and mu = -0.04
    baw_prices =  {20.42,11.25,4.4,1.12,0.18};
    testTemplateBAW(-.04,.08,.2,.25,false,baw_prices);
    baw_prices =  {20.25,11.15,4.35,1.11,0.18};
    testTemplateBAW(-.04,.12,.2,.25,false,baw_prices);
    baw_prices =  {21.46,13.93,8.27,4.52,2.3};
    testTemplateBAW(-.04,.08,.4,.25,false,baw_prices);
    baw_prices =  {20.98,12.64,6.37,2.65,0.92};
    testTemplateBAW(-.04,.08,.2,.5,false,baw_prices);

    // Put and mu = 0.0
    baw_prices =  {20.0,10.58,3.93,0.94,0.15};
    testTemplateBAW(.0,.08,.2,.25,false,baw_prices);
    baw_prices =  {20.0,10.53,3.9,0.93,0.15};
    testTemplateBAW(0.0,.12,.2,.25,false,baw_prices);
    baw_prices =  {20.93,13.39,7.84,4.23,2.12};
    testTemplateBAW(.0,.08,.4,.25,false,baw_prices);
    baw_prices =  {20.04,11.48,5.48,2.15,0.7};
    testTemplateBAW(.0,.08,.2,.5,false,baw_prices);

    // Put and T=3.0
    baw_prices =  {26.25,20.64,15.99,12.22,9.23};
    testTemplateBAW(-0.04,.08,.2,3.0,false,baw_prices);
    baw_prices =  {22.4,16.5,12.03,8.69,6.22};
    testTemplateBAW(.0,.08,.2,3.0,false,baw_prices);
    baw_prices =  {20.33,13.56,9.11,6.12,4.12};
    testTemplateBAW(.04,.08,.2,3.0,false,baw_prices);
    baw_prices =  {20.0,11.63,6.96,4.26,2.64};
    testTemplateBAW(.08,.08,.2,3.0,false,baw_prices);

}

int main() {

    testBaroneAdesiWhalley();
    std::cout << "All tests of the Barone-Adesi Whaley object have been passed successfully" << std::endl;
    return 0;
}