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
int main () {
	double S0 = 21;
	double T= 1; //time steps
	double Sims= 100000; //number of simulations
	double K= 20;
	double r= 0.05;
	double sigma= 0.2;
	double CallPrice= MonteCarloEuroCallPrice(Sims, S0, K, r, T, sigma);
	cout << "The call price is: "<< CallPrice <<endl;
	double PutPrice= MonteCarloEuroPutPrice(Sims, S0, K, r, T, sigma);
	cout << "The put price is: "<< PutPrice <<endl;
return 0;
	
}