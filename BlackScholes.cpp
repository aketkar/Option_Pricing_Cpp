//Author: Aishwarya Ketkar
//Date: 02/15/2015
#include <iostream> 
#include <cmath>
#include <stdio.h>
#include <vector>
using namespace std;
//Calculates the normal distribution CDF
double NormalCdfCalc (double x) 
{
    double a1 =  0.2548;
    double a2 = -0.2844;
    double a3 =  1.4214;
    double a4 = -1.4531;
    double a5 =  1.0614;
    double p  =  0.3275;
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
 
    return 0.5*(1.0 + sign*y);
}
//Pricing the Call option
void BlackScholesCall(double StockPriceAt0, double TimeToMaturity, double StrikePrice, 
                        double RiskFreeRate, double Sigma){
    double PriceOfCall;
    double d1, d2;
    d1= (log(StockPriceAt0/StrikePrice)+(RiskFreeRate+(Sigma*Sigma/2)*TimeToMaturity))/
        pow(TimeToMaturity,StrikePrice/ 2);
    d2= d1- pow(TimeToMaturity,StrikePrice/ 2);
    PriceOfCall= StockPriceAt0*NormalCdfCalc(d1) - (NormalCdfCalc(d2)*StrikePrice*exp(-RiskFreeRate*TimeToMaturity));
    cout << PriceOfCall << " is the call price" <<endl;
}
//Pricing the Put Option
void BlackScholesPut(double StockPriceAt0, double TimeToMaturity, double StrikePrice, 
                        double RiskFreeRate, double Sigma){
    double PriceOfPut;
    double d1, d2;
    d1= (log(StockPriceAt0/StrikePrice)+(RiskFreeRate+(Sigma*Sigma/2)*TimeToMaturity))/
        pow(TimeToMaturity,StrikePrice/ 2);
    d2= d1- pow(TimeToMaturity,StrikePrice/ 2);
    PriceOfPut= NormalCdfCalc(-d2)*StrikePrice*exp(-RiskFreeRate*TimeToMaturity) - StockPriceAt0*NormalCdfCalc(-d1);
    cout << PriceOfPut << " is the put price" <<endl;
}

//Driver function
int main() {
    double S0= 18;
    double T= 4;
    double K= 1.10;
    double r= 0.05;
    double sigma= 2.43;
    BlackScholesCall(S0, T, K, r, sigma);
    BlackScholesPut(S0, T, K, r, sigma);
return 0;
}
