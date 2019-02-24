#include <iostream>
#include <cstdlib>
#include <fstream>
#include <limits>
using namespace std;

#define Nparticles 30
#define T_MAX      1000
#define NFC_MAX    1000000
#define W_0        0.9
#define W_T        0.4
#define MAX_V      2.0
#define c1         2.0      
#define c2         2.0
#define Nvariables 30

#define Rand()     ((double)rand()/RAND_MAX);

class Particle {
    public:
        double *x;
        double *v;
        double fitness;
        double pBest;
        double *xBest;
};

class Swarm {
    private:
        Particle *P;
        int gBest; //index
        double gBestValue;
    public:
        Swarm();
        void setGBestIndex(int index);
        void setGBestValue(double gBestValue);
        void setParticleXV(int i, int j, double x, double v);
        void setParticleFitness(int i, double fit);
        void setParticlePBest(int i, double value, double *x);
        Particle getParticleValue(int i);
        int getGbest();
        double getGbestValue();
        ~Swarm();
};

Swarm::Swarm() {
    gBest = 0;

    P = (Particle *)malloc(sizeof(Particle)*Nparticles);
    if (P==NULL) {
        fprintf(stderr, "Cannot allocate memory for %d Particles\n", Nparticles);
        exit (1);
    } 

    for(int i = 0; i < Nparticles; i++) {
        P[i].x = (double *)malloc(sizeof(double)*Nvariables);
        if (P[i].x==NULL) {
            fprintf(stderr, "Cannot allocate memory for %d x\n", Nvariables);
            exit (1);
        } 
        P[i].v = (double *)malloc(sizeof(double)*Nvariables);
        if (P[i].v==NULL) {
            fprintf(stderr, "Cannot allocate memory for %d v\n", Nvariables);
            exit (1);
        } 
        P[i].xBest = (double *)malloc(sizeof(double)*Nvariables);
        if (P[i].xBest==NULL) {
            fprintf(stderr, "Cannot allocate memory for %d xBest\n", Nvariables);
            exit (1);
        } 
    }
    
}

void Swarm::setGBestIndex(int index) {
    this->gBest = index;
}

void Swarm::setGBestValue(double gBestValue) {
    this->gBestValue = gBestValue;
}

void Swarm::setParticleXV(int i, int j, double x, double v) {
    if (x != NULL) P[i].x[j] = x;
    if (v != NULL) P[i].v[j] = v;
}

void Swarm::setParticleFitness(int i, double fit) {
    P[i].fitness = fit;
}

void Swarm::setParticlePBest(int i, double value, double *x) {
    P[i].pBest = value;
    for(int j = 0; j < Nvariables; j++)
    {
        P[i].xBest[j] = P[i].x[j];
    }    
}

Particle Swarm::getParticleValue(int i) {
    return P[i];
}

int Swarm::getGbest() {
    return gBest;
}

double Swarm::getGbestValue() {
    return gBestValue;
}

Swarm::~Swarm() {

}

class PSO {
    private:
        int nfc;
        double w;
        double *maxV; 
    public:
        Swarm swarm;
        PSO();
        void initialize();
        void evolution();
        void updateBest(int i);
        void calculateVMax();
        void particleMovement();
        void evaluate(int i);
        void evaluateSwarm();
        ~PSO();
};

PSO::PSO() {
    maxV = (double *)malloc(sizeof(double)*Nvariables);
    if (maxV==NULL) {
        fprintf(stderr, "Cannot allocate memory for %d maxV\n", Nvariables);
        exit (1);
    } 
}

void PSO::evolution() {
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

void PSO::initialize() {
    swarm.setGBestValue(numeric_limits<double>::max());

    for(int i = 0; i < Nparticles; i++) {
        for(int j = 0; j < Nvariables; j++) {
            double x = -5.12 + 10.24 * Rand();
            //double x = Rand();
            double v = 0.0;
            swarm.setParticleXV(i, j, x, v);
        }
        evaluate(i);
        updateBest(i);

        double fitness = swarm.getParticleValue(i).fitness;
        double gbest = swarm.getGbestValue();

        if (fitness < gbest) {
            swarm.setGBestValue(fitness);
            swarm.setGBestIndex(i);
        }

    }

    Particle best = swarm.getParticleValue(swarm.getGbest());
    printf("0 : ");
    for (int j = 0; j < Nvariables; j++)
        printf("%g ", best.xBest[j]);
    printf(" = %e\n", best.pBest);
}

void PSO::evaluate(int i) {
    int index = i;
    // double fitness = 0.0;

    // for(int k = 0 ; k < Nvariables ; k++ ) {
    //     int x = swarm.getParticleValue(index).x[k];
    //     fitness +=  x * x;
    // }

    double fitness = 0, temp = 0;
    for(int k = 0 ; k < Nvariables ; k++ ) {
        for(int l = 0; l < k; l++)
        {
            double x = swarm.getParticleValue(index).x[l];
            temp += x;
        }
        
        fitness +=  temp * temp;
    }

    
    swarm.setParticleFitness(index, fitness);
}

void PSO::evaluateSwarm() {
    for(int i = 0; i < Nparticles; i++) {
        evaluate(i);    
        nfc++;

        if (nfc % 5000 == 0) {
            Particle best = swarm.getParticleValue(swarm.getGbest());
            printf("%d : ", nfc);
            // for (int j = 0; j < Nvariables; j++)
            //     printf("%g ", best.xBest[j]);
            printf(" = %g\n", best.pBest);

            ofstream file;
            file.open("BPSO.txt", ios::out | ios::app);
            file << nfc << "," << best.pBest << "\n";
            file.close();
            //cout << "PSO SPHERE nfc" << nfc << " \tbestfit " << best.pBest << "\n";
        }
    }

    for(int n = 0; n < Nparticles; n++) {
        Particle par = swarm.getParticleValue(n);
        if (par.fitness < par.pBest ) {
            swarm.setParticlePBest(n, par.fitness, par.x);

            if (par.fitness < swarm.getGbestValue() ) {
                swarm.setGBestIndex(n);
                swarm.setGBestValue(par.fitness);
            }
        }
    }
}

void PSO::updateBest(int i) {
    Particle par = swarm.getParticleValue(i);
    swarm.setParticlePBest(i, par.fitness, par.x);
}

void PSO::calculateVMax() {
    double xmin[Nparticles], xmax[Nparticles];

    for (int d = 0; d < Nvariables; d++) {
        xmin[d] = xmax[d] = swarm.getParticleValue(0).x[d];
        for (int n = 1; n < Nparticles; n++) {
            double pos = swarm.getParticleValue(n).x[d];
            if (pos < xmin[d])
                xmin[d] = pos;
            if (pos > xmax[d])
                xmax[d] = pos;
        }
        maxV[d] = xmax[d] - xmin[d];
    }
}

void PSO::particleMovement() {
    int n, d;

    for (n = 0; n < Nparticles ; n++) {
        Particle par = swarm.getParticleValue(n);
        Particle bPar = swarm.getParticleValue(swarm.getGbest());
        // update velocities
        for(d = 0; d < Nvariables ; d++ ) {
            double r1 = Rand();
            double r2 = Rand();
            double newV = w * par.v[d] + c1 * r1 * (par.xBest[d] - par.x[d]) + c2 * r2 * (bPar.x[d] - par.x[d]);
            swarm.setParticleXV(n, d, NULL, newV);
            // check v with its dimensional maxV
            if ( swarm.getParticleValue(n).v[d] > maxV[d] ) swarm.setParticleXV(n, d, NULL, maxV[d]);
            else if ( swarm.getParticleValue(n).v[d] < -maxV[d] ) swarm.setParticleXV(n, d, NULL, -maxV[d]);
        }
        // update positions
        Particle newPar = swarm.getParticleValue(n);
        for (d = 0; d < Nvariables ; d++) {
            //newPar.x[d] = newPar.x[d] + newPar.v[d];
            swarm.setParticleXV(n, d, newPar.x[d] + newPar.v[d], NULL);
        }
    }
}

PSO::~PSO() {

}

int main() {

    PSO pso;
    pso.evolution();

    return 0;
}