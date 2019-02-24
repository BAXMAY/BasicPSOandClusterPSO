#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <random>
#include <vector>

using namespace std;

#define Nparticles  30
#define T_MAX       1000
#define NFC_MAX     1000000
#define W_0         0.9
#define W_T         0.4
#define MAX_V       2.0
#define c1          2.0
#define c2          2.0
#define Nvariables  30
#define Niterations 300

#define Rand()      ((double)rand()/RAND_MAX);

struct Particle
{
    vector<double> x;
    vector<double> v;
    vector<double> xBest;
    double fitness;
    double pBest;
};

vector<Particle> k_means(const vector<Particle> &particles, size_t k);

class Swarm
{
  public:
    vector<Particle> P;
    vector<Particle> centroids;
    int gBestIndex; 
    double gBestValue;
    double w;
    void initialize();
    void evaluate(int index);
    void updateBest(int index);
    void print();
    Swarm(/* args */);
};

Swarm::Swarm(/* args */)
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
    
    double dw = (W_0 - W_T) / (NFC_MAX / Nparticles);
    w = W_0;
    
    initialize();

    centroids = k_means(P, 5);
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
}

void Swarm::evaluate(int index)
{   
    double fitness = 0, temp = 0;
    for(int i = 0 ; i < Nvariables ; ++i ) 
    {
        for(int j = 0; j < i; ++j)
        {
            temp += P[index].x[j];
        }      
        fitness +=  temp * temp;
    }

    P[index].fitness = fitness;
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

double square(double val)
{
    return val * val;
}

double eucildean_distance(Particle first, Particle second) {
    double distance = 0.0;

    for(size_t i = 0; i < first.x.size(); ++i)
    {
        distance += square(first.x[i] - second.x[i]);
    }
    distance += square(first.fitness - second.fitness);
    
    return sqrt(distance);
}

vector<Particle> k_means(const vector<Particle> &particles, size_t k)
{
    static random_device seed;
    static mt19937 random_number_generator(seed());
    uniform_int_distribution<size_t> indices(0, particles.size() - 1);

    vector<Particle> means(k); // cluster centroid particle
    for(auto& cluster : means)
    {
        cluster = particles[indices(random_number_generator)];
    }

    // Compute E(t)
    // double currentE;
    // double nextE = NULL;    
    
    vector<size_t> assignments(particles.size()); // which cluster that particle belong to
    //do
    for(size_t iteration = 0; iteration < Niterations; ++iteration) 
    {
        //if (nextE != NULL) 
        //{
        //    currentE = nextE;
        //}

        // Find assignments.
        for (size_t particle = 0; particle < particles.size(); ++particle)
        {
            auto best_distance = std::numeric_limits<double>::max();
            size_t best_cluster = 0;
            for (size_t cluster = 0; cluster < k; ++cluster)
            {
                const double distance = eucildean_distance(particles[particle], means[cluster]); //distance
                if (distance < best_distance)
                {
                    best_distance = distance;
                    best_cluster = cluster;
                }
            }
            assignments[particle] = best_cluster;
        }

        // Sum up and count particles for each cluster.
        vector<Particle> new_means(k);
        for (size_t i = 0; i < k; ++i)
        {
            new_means[i].x.resize(Nvariables);
            new_means[i].v.resize(Nvariables);
            new_means[i].xBest.resize(Nvariables);
        }
        vector<size_t> counts(k, 0); // number of particle in cluster, use for find new mean
        for (size_t particle = 0; particle < particles.size(); ++particle)
        {
            const auto cluster = assignments[particle];
            for(size_t var = 0; var < particles[particle].x.size(); ++var)
            {
                new_means[cluster].x[var] += particles[particle].x[var];
            }
            counts[cluster] += 1;
        }

        // Divide sums by counts to get new centroids.
        for (size_t cluster = 0; cluster < k; ++cluster)
        {
            double fitness = 0, sum = 0;
            // Turn 0/0 into 0/1 to avoid zero division.
            const auto count = max<size_t>(1, counts[cluster]);
            for(size_t var = 0; var < new_means[cluster].x.size(); ++var)
            {
                means[cluster].x[var] += new_means[cluster].x[var] / count;
                // Compute New Centroid fitness
                for (int j = 0; j < var; ++j)
                {
                    sum += means[cluster].x[j];
                }
                fitness += sum * sum;
            }
            means[cluster].fitness = fitness;
        }
    }
        // Compute E(t+1)
        //nextE = quality(particles, means);
        
    //} while (currentE != nextE);
    
    return means;
}

int main(int argc, const char *argv[])
{
    Swarm sw;
    
    sw.print();
    return 0;
}