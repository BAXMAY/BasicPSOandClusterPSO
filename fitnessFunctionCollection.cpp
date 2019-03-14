#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <random>
#include <vector>
#include <ctime>

#define Nvariables  30
#define Rand()      ((double)rand()/RAND_MAX);

using namespace std;

double sphere(vector<double> x) 
{
    double fit = 0.0;

    for(int k = 0 ; k < Nvariables ; k++ ) {
        fit +=  x[k] * x[k];
    }

    return fit;
}

double schwefel(vector<double> x)
{
    double fit = 0, temp = 0;

    for(int i = 0 ; i < Nvariables ; ++i ) 
    {
        for(int j = 0; j < i; ++j)
        {
            temp += x[j];
        }      
        fit +=  temp * temp;
    }

    return fit;
}

double rosenbrock(vector<double> x)
{
    double fit = 0;
    
    for(int i = 0 ; i < Nvariables - 1 ; ++i ) 
    {
        fit +=  100 * pow(x[i] * x[i] - x[i + 1], 2) + pow(x[i] - 1.0, 2);
    }

    return fit;
}

double quarticNoise(vector<double> x)
{
    double fit = 0;
    
    for(int i = 0 ; i < Nvariables ; ++i ) 
    {
        fit +=  (1.0 + i) * pow(x[i], 4);
    }

    return fit + Rand();
}

double ackley(vector<double> x)
{
    double fit = 0, sum1 = 0, sum2 = 0;
    
    for(int i = 0 ; i < Nvariables ; ++i ) 
    {
        sum1 += x[i] * x[i];
        sum2 += cos(2 * M_PI * x[i]);    
    }

    sum1 = -0.2 * sqrt(sum1 / Nvariables);
    sum2 = sum2 / Nvariables;
    fit = -20.0 * exp(sum1) - exp(sum2) + 20.0 + M_E;

    return fit;
}

double griewank(vector<double> x)
{
    double fit = 0, sum = 0, mul = 1;
    
    for(int i = 0 ; i < Nvariables ; ++i ) 
    {
        sum += x[i] * x[i];
        mul *= cos(x[i] / (i + 1));
    }

    fit = sum / 4000 - mul + 1;

    return fit;
}

double rastrigin(vector<double> x)
{
    double fit = 0, temp = 0;
    
    for(int i = 0 ; i < Nvariables ; ++i ) 
    {
        fit +=  x[i] * x[i] - 10.0 * cos(2.0 * M_PI * x[i]) + 10.0;
    }

    return fit;
}

double penalized(vector<double> x)
{
    vector<double> y;
    y.resize(Nvariables);
    double fit = 0;
    double sumu = 0;

    for(int i = 0 ; i < Nvariables ; ++i ) 
    {
        y[i] = 1.0 + (x[i] + 1.0)/4.0;
    }
    
    for (int i = 0; i < Nvariables; i++)
    {
        if (x[i] > 10) sumu += 100 * pow(x[i] - 10, 4);
        else if (x[i] < -10) sumu += 100 * pow(-x[i] - 10, 4);
        else sumu += 0;
    }

    for (int i = 0; i < Nvariables - 1; i++)
    {
        fit += (y[i] - 1.0) * (y[i] - 1.0) * (1.0 + 10.0 * pow(sin(M_PI * y[i + 1]), 2.0));
    }

    fit += 10.0 * pow(sin(M_PI * y[0]), 2.0);    // term 1
    fit += (y[Nvariables - 1] - 1.0) * (y[Nvariables - 1] - 1.0); // term 3
    fit = fit * M_PI / Nvariables;
    fit += sumu; // term 4

    return fit;
}