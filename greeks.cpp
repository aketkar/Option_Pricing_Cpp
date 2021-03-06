//Author: Aishwarya Ketkar
//Date: 02/15/2015
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <vector>
using namespace std;
double NormalCdfCalc(double x)
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
void CallDelta(double StockPriceAt0, double Sigma, double TimeToMaturity, double StrikePrice, double RiskFreeRate,
    double AnnualDivYield) {
    double d1;
    d1= (log(StockPriceAt0/StrikePrice)+(RiskFreeRate- AnnualDivYield + (Sigma*Sigma/2))*TimeToMaturity)/
        pow(TimeToMaturity,StrikePrice/ 2);
    double DeltaCall= exp(-AnnualDivYield*TimeToMaturity)*NormalCdfCalc(d1);
    cout << "Delta of this call option is " <<DeltaCall <<endl;
}
void PutDelta(double StockPriceAt0, double Sigma, double TimeToMaturity, double StrikePrice, double RiskFreeRate,
    double AnnualDivYield) {
    double d1;
    d1= (log(StockPriceAt0/StrikePrice)+(RiskFreeRate- AnnualDivYield + (Sigma*Sigma/2))*TimeToMaturity)/
        pow(TimeToMaturity,StrikePrice/ 2);
    double DeltaPut= -exp(-AnnualDivYield*TimeToMaturity)*NormalCdfCalc(-d1);
    cout << "Delta of this put option is " <<DeltaPut <<endl;
}
void StockGamma(double StockPriceAt0, double Sigma, double TimeToMaturity, double StrikePrice, double RiskFreeRate,
    double AnnualDivYield) {
    double d1;
    d1= (log(StockPriceAt0/StrikePrice)+(RiskFreeRate- AnnualDivYield + (Sigma*Sigma/2))*TimeToMaturity)/
        pow(TimeToMaturity,StrikePrice/ 2);
    double gamma= exp(-AnnualDivYield*TimeToMaturity)*NormalCdfCalc(d1)/StockPriceAt0*Sigma*pow(TimeToMaturity, 0.5);
    cout << "Gamma of this option is " <<gamma <<endl;
}
void ThetaCall(double StockPriceAt0, double Sigma, double TimeToMaturity, double StrikePrice, double RiskFreeRate,
    double AnnualDivYield) {
    double d1= (log(StockPriceAt0/StrikePrice)+(RiskFreeRate- AnnualDivYield + (Sigma*Sigma/2))*TimeToMaturity)/
        pow(TimeToMaturity,StrikePrice/ 2);
    double d2= d1- pow(TimeToMaturity,StrikePrice/ 2);
    double ThetaCall= -exp(-AnnualDivYield*TimeToMaturity)*StockPriceAt0*NormalCdfCalc(d1)*Sigma/ 2*sqrt(TimeToMaturity)-
                    RiskFreeRate*StrikePrice*exp(-RiskFreeRate*TimeToMaturity)*NormalCdfCalc(d2)+ 
                    AnnualDivYield*StockPriceAt0*exp(-AnnualDivYield*TimeToMaturity)*NormalCdfCalc(d1);
    cout << "Theta of this call option is "<<ThetaCall <<endl;
}
void ThetaPut(double StockPriceAt0, double Sigma, double TimeToMaturity, double StrikePrice, double RiskFreeRate,
    double AnnualDivYield){
    double d1= (log(StockPriceAt0/StrikePrice)+(RiskFreeRate- AnnualDivYield + (Sigma*Sigma/2))*TimeToMaturity)/
        pow(TimeToMaturity,StrikePrice/ 2);
    double d2= d1- pow(TimeToMaturity,StrikePrice/ 2);
    double ThetaPut= -exp(-AnnualDivYield*TimeToMaturity)*StockPriceAt0*NormalCdfCalc(d1)*Sigma/ 2*sqrt(TimeToMaturity)+
                    RiskFreeRate*StrikePrice*exp(-RiskFreeRate*TimeToMaturity)*NormalCdfCalc(-d2)- 
                    AnnualDivYield*StockPriceAt0*exp(-AnnualDivYield*TimeToMaturity)*NormalCdfCalc(-d1);
    cout << "Theta of this put option is "<<ThetaPut<<endl;;
}
void vega(double StockPriceAt0, double Sigma, double TimeToMaturity, double StrikePrice, double RiskFreeRate,
    double AnnualDivYield) {
    double d1;
    d1= (log(StockPriceAt0/StrikePrice)+(RiskFreeRate- AnnualDivYield + (Sigma*Sigma/2))*TimeToMaturity)/
        pow(TimeToMaturity,StrikePrice/ 2);
    double vega= StockPriceAt0*exp(-AnnualDivYield*TimeToMaturity)*NormalCdfCalc(d1)*sqrt(TimeToMaturity);
    cout << "Vega of this option is " <<vega <<endl;
}
void RhoCall(double StockPriceAt0, double Sigma, double TimeToMaturity, double StrikePrice, double RiskFreeRate,
    double AnnualDivYield) {
    double d1= (log(StockPriceAt0/StrikePrice)+(RiskFreeRate- AnnualDivYield + (Sigma*Sigma/2))*TimeToMaturity)/
        pow(TimeToMaturity,StrikePrice/ 2);
    double d2= d1- pow(TimeToMaturity,StrikePrice/ 2);
    double RhoC= StrikePrice*TimeToMaturity*exp(-RiskFreeRate*TimeToMaturity)*NormalCdfCalc(d2);
    cout << "The Rho for this call option is: "<<RhoC <<endl;
}
void RhoPut(double StockPriceAt0, double Sigma, double TimeToMaturity, double StrikePrice, double RiskFreeRate,
    double AnnualDivYield) {
    double d1= (log(StockPriceAt0/StrikePrice)+(RiskFreeRate- AnnualDivYield + (Sigma*Sigma/2))*TimeToMaturity)/
        pow(TimeToMaturity,StrikePrice/ 2);
    double d2= d1- pow(TimeToMaturity,StrikePrice/ 2);
    double RhoP= -StrikePrice*TimeToMaturity*exp(-RiskFreeRate*TimeToMaturity)*NormalCdfCalc(-d2);
    cout << "The Rho for this call option is: "<<RhoP <<endl;
}
