#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>
#include <vector>
using namespace std;
int main()
{
    random_device rd;
    mt19937 generator(rd());
 
    normal_distribution<> distributionParams(0,1.2);
 
    std::vector<double> v;
    for(int n=0; n<10000; ++n) {
        double newno= (distributionParams(generator));
         //gen is the engine name
        v.push_back(newno);
    }
    cout <<v.at(2) <<endl;
    /*for(auto p : hist) {
        std::cout << std::fixed << std::setprecision(1) << std::setw(2)
                  << p.first << ' ' << std::string(p.second/200, '*') << '\n';
    }*/
}
