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
double BlackScholesCall(double StockPriceAt0, double TimeToMaturity, double StrikePrice, 
                        double RiskFreeRate, double Sigma){
    double PriceOfCall;
    double d1, d2;
    d1= (log(StockPriceAt0/StrikePrice)+(RiskFreeRate+(Sigma*Sigma/2)*TimeToMaturity))/
        pow(TimeToMaturity,StrikePrice/ 2);
    d2= d1- pow(TimeToMaturity,StrikePrice/ 2);
    PriceOfCall= StockPriceAt0*NormalCdfCalc(d1) - (NormalCdfCalc(d2)*StrikePrice*exp(-RiskFreeRate*TimeToMaturity));
    return PriceOfCall;
}
double BlackScholesPut(double StockPriceAt0, double TimeToMaturity, double StrikePrice, 
                        double RiskFreeRate, double Sigma){
    double PriceOfPut;
    double d1, d2;
    d1= (log(StockPriceAt0/StrikePrice)+(RiskFreeRate+(Sigma*Sigma/2)*TimeToMaturity))/
        pow(TimeToMaturity,StrikePrice/ 2);
    d2= d1- pow(TimeToMaturity,StrikePrice/ 2);
    PriceOfPut= NormalCdfCalc(-d2)*StrikePrice*exp(-RiskFreeRate*TimeToMaturity) - StockPriceAt0*NormalCdfCalc(-d1);
    cout << "Put price based on Black Scholes model:" << PriceOfPut <<endl;
    return PriceOfPut;
}
double PutCallParity(double StockPriceAt0, double TimeToMaturity, double StrikePrice, 
                        double RiskFreeRate, double Sigma){
    double PriceOfPut;
    double call= BlackScholesCall(StockPriceAt0, TimeToMaturity, StrikePrice, RiskFreeRate, Sigma);
    PriceOfPut= call - StockPriceAt0 + StrikePrice*exp(-RiskFreeRate*TimeToMaturity);
    cout << "Put price based on the put call parity model:"<<PriceOfPut <<endl;
    return PriceOfPut;
}

    
int main() {
    double S0= 18;
    double T= 4;
    double K= 1.10;
    double r= 0.05;
    double sigma= 2.43;
    double callPrice= BlackScholesCall(S0, T, K, r, sigma);
    cout << "Call price based on BlackScholes model: " << callPrice << endl;
    double PutPrice= BlackScholesPut(S0, T, K, r, sigma);
    double Parity= PutCallParity(S0, T, K, r, sigma);
    double epsilon= 0.00000001;
    if (fabs(PutPrice- Parity)< epsilon) {
        cout << "There! Put-Call parity indeed holds!" <<endl;
    }
    else {
        cout << "Not equal!" <<endl;
    }
return 0;
}