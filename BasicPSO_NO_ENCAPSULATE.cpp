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
#define W_0         0.9
#define W_T         0.4
#define MAX_V       2.0
#define c1          2.0
#define c2          2.0

#define Rand()      ((double)rand()/RAND_MAX);

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
}

void Swarm::initialize() 
{
    gBestValue = numeric_limits<double>::max();

    for(int i = 0; i < Nparticles; ++i) 
    {
        for(int j = 0; j < Nvariables; ++j) 
        {
            P[i].x[j] = -5.12 + 10.24 * Rand();
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
}

void Swarm::evaluate(int index)
{   
    P[index].fitness = fitness(P[index].x);
}

double Swarm::fitness(vector<double> x)
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
        for (d = 0; d < Nvariables ; d++) {
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
    srand(time(0));
    Swarm sw;
    sw.evolution();
    //sw.print();
    return 0;
}
