# include <array>

# include <vector>

#include <iostream>
#include <cmath>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

int main() 
{
    plt::plot({1,3,2,4});
    plt::show();
    std::array<std::array<double,2>,2> k;
    std::vector<std::vector<double> > K(10, std::vector<double>(10,0.0));
/*
    k[0][0] = 1.0;
    k[0][1] = 1.0;
    k[1][0] = 1.0;
    k[1][1] = 1.0;

    std::vector<double> R(5, 0.0 );
    return R(1,3);
    */
    std::cout<< pow(5,2)  <<std::endl;
    return pow(5,2);
}
