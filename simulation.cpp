//Author: Aishwarya Ketkar
//MonteCarlo Pricing engine to price plain Vanilla call and put options
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <random>
using namespace std;
double GenerateRandomNormals(double mu, double sigma) {
	random_device rd;
    mt19937 gen(rd());
    normal_distribution<> distribution(mu, sigma);
    double newno= (distribution(gen));
 return newno;
}
double MonteCarloEuroCallPrice(double Sims, double S0, double K, double r, double T, double Sigma) {
	double StockSpot= S0*exp((r-0.5*Sigma*Sigma)*T);
	double TotalPayoffSum= 0.0;
	for (int i=0; i<=Sims; i++) {
		double GaussVar= GenerateRandomNormals(r, Sigma);
		double St= StockSpot- exp(sqrt(Sigma*Sigma*T)*GaussVar);
		TotalPayoffSum += max(St-K, 0.0);	
	}
double newPayoff= exp(-r*T)*(TotalPayoffSum/Sims);
return newPayoff;
}
double MonteCarloEuroPutPrice(double Sims, double S0, double K, double r, double T, double Sigma) {
	double StockSpot= S0*exp((r-0.5*Sigma*Sigma)*T);
	double TotalPayoffSum= 0.0;
	for (int i=0; i<=Sims; i++) {
		double GaussVar= GenerateRandomNormals(r, Sigma);
		double St= StockSpot- exp(sqrt(Sigma*Sigma*T)*GaussVar);
		TotalPayoffSum += max(K-St, 0.0);	
	}
double newPayoff= exp(-r*T)*(TotalPayoffSum/Sims);
return newPayoff;
}
