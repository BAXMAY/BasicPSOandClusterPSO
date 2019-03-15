#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <random>
#include <vector>
#include <ctime>

using namespace std;

#define Nparticles  30
#define Nvariables  30
#define T_MAX       1000
#define NFC_MAX     1000000
#define W_0         0.72894
#define W_T         0.4
#define MAX_V       2.0
#define c1          1.49618
#define c2          1.49618

#define Rand()      ((double)rand()/RAND_MAX);

int probNo;

double sphere(vector<double> x) 
{
    // -100 + (100 - (-100)) * Rand();
    double fit = 0.0;

    for(int k = 0 ; k < Nvariables ; k++ ) {
        fit +=  x[k] * x[k];
    }

    return fit;
}

double schwefel(vector<double> x)
{
    // -5.12 + (5.12 - (-5.12)) * Rand();
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
    // -30 + (30 - (-30)) * Rand();
    double fit = 0;
    
    for(int i = 0 ; i < Nvariables - 1 ; ++i ) 
    {
        fit +=  100 * pow(x[i] * x[i] - x[i + 1], 2) + pow(x[i] - 1.0, 2);
    }

    return fit;
}

double quarticNoise(vector<double> x)
{
    // -1.28 + (1.28 - (-1.28)) * Rand();
    double fit = 0;
    
    for(int i = 0 ; i < Nvariables ; ++i ) 
    {
        fit +=  (1.0 + i) * pow(x[i], 4);
    }

    return fit + Rand();
}

double ackley(vector<double> x)
{
    // -32 + (32 - (-32)) * Rand();
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
    // -600 + (600 - (-600)) * Rand();
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
    // -5.12 + (5.12 - (-5.12)) * Rand();
    double fit = 0, temp = 0;
    
    for(int i = 0 ; i < Nvariables ; ++i ) 
    {
        fit +=  x[i] * x[i] - 10.0 * cos(2.0 * M_PI * x[i]) + 10.0;
    }

    return fit;
}

double penalized(vector<double> x)
{
    // -50 + (50 - (-50)) * Rand();
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

struct Particle
{
    vector<double> x;
    vector<double> v;
    vector<double> xBest;
    double fitness;
    double pBest;
};

class Swarm
{
  public:
    vector<Particle> P;
    vector<Particle> centroids;
    int gBestIndex; 
    double gBestValue;
    double w;
    int nfc;
    double maxV[Nvariables];
    void initialize();
    void evolution();
    double randx();
    double fitness(vector<double> x);
    void evaluate(int index);
    void evaluateSwarm();
    void updateBest(int index);
    void calculateVMax();
    void particleMovement();
    void print();
    Swarm();
};

Swarm::Swarm()
{
    P.resize(Nparticles);
    for(size_t i = 0; i < Nparticles; ++i)
    {
        P[i].x.resize(Nvariables);
        P[i].v.resize(Nvariables); 
        P[i].xBest.resize(Nvariables);         
    }

    centroids.reserve(5);
    for(size_t i = 0; i < 5; ++i)
    {
        P[i].x.resize(Nvariables);
        P[i].v.resize(Nvariables); 
        P[i].xBest.resize(Nvariables);         
    }
}

void Swarm::evolution() 
{
    double dw = (W_0 - W_T) / (NFC_MAX / Nparticles);
    w = W_0;
    
    initialize();

    nfc = 0;
    while(nfc < NFC_MAX) {
        calculateVMax();
        particleMovement();  
        evaluateSwarm();
        w -= dw;     
    }

    ofstream file;
    file.open("BPSO.txt", ios::out | ios::app);
    file << "===========================================" << "\n\n";
    file.close();
}

double Swarm::randx()
{
    switch (::probNo)
    {
    case 0:
        return -100 + (100 - (-100)) * Rand();
        break;
    case 1:
        return -5.12 + (5.12 - (-5.12)) * Rand();
        break;
    case 2:
        return -30 + (30 - (-30)) * Rand();
        break;
    case 3:
        return -1.28 + (1.28 - (-1.28)) * Rand();
        break;
    case 4:
        return -32 + (32 - (-32)) * Rand();
        break;
    case 5:
        return -600 + (600 - (-600)) * Rand();
        break;
    case 6:
        return -5.12 + (5.12 - (-5.12)) * Rand();
        break;
    case 7:
        return -50 + (50 - (-50)) * Rand();
        break;
    default:
        break;
    }
}

void Swarm::initialize() 
{
    gBestValue = numeric_limits<double>::max();

    for(int i = 0; i < Nparticles; ++i) 
    {
        for(int j = 0; j < Nvariables; ++j) 
        {
            P[i].x[j] = randx();
            P[i].v[j] = 0.0;
        }
        evaluate(i);
        updateBest(i);

        if (P[i].fitness < gBestValue) 
        {
            gBestValue = P[i].fitness;
            gBestIndex = i;
        }
    }

    printf("0 : ");
    printf(" = %g\n", gBestValue);

    ofstream file;
    file.open("BPSO.txt", ios::out | ios::app);
    file << "0," << gBestValue<< "\n";
    file.close();
}

void Swarm::evaluate(int index)
{   
    P[index].fitness = fitness(P[index].x);
}

double Swarm::fitness(vector<double> x)
{
    switch (::probNo)
    {
    case 0:
        return sphere(x);
        break;
    case 1:
        return schwefel(x);
        break;
    case 2:
        return rosenbrock(x);
        break;
    case 3:
        return quarticNoise(x);
        break;
    case 4:
        return ackley(x);
        break;
    case 5:
        return griewank(x);
        break;
    case 6:
        return rastrigin(x);
        break;
    case 7:
        return penalized(x);
        break;
    default:
        break;
    }
}

void Swarm::calculateVMax() {
    double xmin[Nparticles], xmax[Nparticles];

    for (int d = 0; d < Nvariables; d++) {
        xmin[d] = xmax[d] = P[0].x[d];
        for (int n = 1; n < Nparticles; n++) {
            double pos = P[n].x[d];
            if (pos < xmin[d])
                xmin[d] = pos;
            if (pos > xmax[d])
                xmax[d] = pos;
        }
        maxV[d] = xmax[d] - xmin[d];
    }
}

void Swarm::particleMovement() {
    int n, d;

    for (n = 0; n < Nparticles ; n++) {
        Particle par = P[n];
        Particle bPar = P[gBestIndex];
        // update velocities
        for(d = 0; d < Nvariables ; d++ ) {
            double r1 = Rand();
            double r2 = Rand();
            P[n].v[d] = w * par.v[d] + c1 * r1 * (P[n].xBest[d] - P[n].x[d]) + c2 * r2 * (P[gBestIndex].x[d] - P[n].x[d]);
            // check v with its dimensional maxV
            if ( P[n].v[d] > maxV[d] ) P[n].v[d] = maxV[d];
            else if ( P[n].v[d] < -maxV[d] ) P[n].v[d] = -maxV[d];
        }
        // update positions
        for (int d = 0; d < Nvariables; d++)
        {
            P[n].x[d] += P[n].v[d];
        }
    }
}

void Swarm::evaluateSwarm() {
    for(int i = 0; i < Nparticles; i++) {
        evaluate(i);    
        nfc++;

        if (nfc % 5000 == 0) {
            Particle best = P[gBestIndex];
            printf("%d : ", nfc);
            printf(" = %g\n", best.pBest);

            ofstream file;
            file.open("BPSO.txt", ios::out | ios::app);
            file << nfc << "," << best.pBest << "\n";
            file.close();
        }
    }

    for(int n = 0; n < Nparticles; n++) {
        if (P[n].fitness < P[n].pBest ) {
            P[n].pBest =  P[n].fitness;
            P[n].xBest =  P[n].x;

            if (P[n].fitness < gBestValue) {
                gBestIndex = n;
                gBestValue = P[n].fitness;
            }
        }
    }
}

void Swarm::updateBest(int index)
{
    P[index].xBest = P[index].x;
    P[index].pBest = P[index].fitness;
}

void Swarm::print()
{
    //cout << gBestValue << endl;
    for(size_t i = 0; i < Nparticles; ++i)
    {
        for(size_t j = 0; j < Nvariables; ++j)
        {
            cout << "P" << i << " X" << j << " : " << P[i].x[j] << endl;
        }
        cout << "///////////////////////////////////////\n" << endl;
    }

    for(size_t i = 0; i < 5; ++i)
    {
        for(size_t j = 0; j < Nvariables; ++j)
        {
            cout << "Cluster" << i << " Var" << j << " : " << centroids[i].x[j] << endl;
        }
        cout << "Cluster" << i << " fitness" << " : " << centroids[i].fitness << endl;
        cout << "///////////////////////////////////////" << endl;
    }
    
}     

int main(int argc, const char *argv[])
{
    char probName[8][14] = {"Sphere", "Schwefel", "Rosenbrock", "QuarticNoise", "Ackley", "Griewank", "Rastrigin", "Penalized"};
    for(::probNo = 0; ::probNo < 8; ::probNo++)
    {
        ofstream file;
        file.open("BPSO.txt", ios::out | ios::app);
        file << probName[::probNo] << "\n";
        file << "#####################################" << "\n\n";
        file.close();
        cout << probName[::probNo] << "\n\n";
        for (int i = 0; i < 10; i++)
        {
            srand(time(0));
            Swarm sw;
            sw.evolution();
        }
    }
    
    //sw.print();
    return 0;
}
