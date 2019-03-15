#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <random>
#include <vector>
#include <ctime>

using namespace std;

#define Niterations 5
#define Nparticles  30
#define Nvariables  30
#define T_MAX       1000
#define NFC_MAX     1000000
#define W_0         0.72894
#define W_T         0.4
#define MAX_V       2.0
#define c1          1.49618
#define c2          1.49618
#define K           5

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

int randomParticle()
{
    return rand() % Nparticles;
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

vector<size_t> k_means(const vector<Particle> &particles, size_t k);

class Swarm
{
  public:
    vector<Particle> P;
    vector<Particle> candidates;
    vector<size_t> cluster;
    int nfc;
    int gBestIndex; 
    int dominants[K];
    double w;
    double gBestValue;    
    double maxV[Nvariables];    
    void initialize();
    void evolution();
    double randx();
    double fitness(vector<double> x);
    void evaluate(int index);
    void evaluateSwarm(int clusIndex);
    void updateBest(int index);
    void updateGBest();
    void calculateVMax();
    void particleMovement(int clusIndex);
    void chooseDominant(int clusIndex);
    void spsa(int clusIndex);
    void print();
    vector<size_t> k_means(const vector<Particle> &particles, size_t k);
    Swarm();
};

Swarm::Swarm()
{
    candidates.resize(K);
    P.resize(Nparticles);
    for(size_t i = 0; i < Nparticles; ++i)
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
    // cluster = k_means(P, K);
    // for(size_t i = 0; i < Nparticles; i++)
    //     {
    //         if(nfc == 0)
    //         cout << cluster[i] << endl;
    //     }
    while(nfc < NFC_MAX) {
        cluster = k_means(P, K);
        // for(size_t i = 0; i < Nparticles; i++)
        // {
        //     if(nfc == 0)
        //     cout << cluster[i] << endl;
        // }
        
        for(size_t c = 0; c < K; ++c)
        {
            chooseDominant(c);
            spsa(c);
            for(int i = 0; i < 50; ++i)
            {
                calculateVMax();
                particleMovement(c);  
                evaluateSwarm(c);
            }             
        }    
        updateGBest();
        w -= dw;
    }

    ofstream file;
    file.open("ClusterPSO.txt", ios::out | ios::app);
    file << "===========================================" << "\n\n";
    file.close();
}

void Swarm::chooseDominant(int clusIndex)
{
    int min = 0;
    for(int i = 1; i < Nparticles; ++i)
    {
        if(min == 0 && cluster[min] != clusIndex && cluster[i] == clusIndex) min = i;
        else if(P[i].fitness < P[min].fitness && cluster[i] == clusIndex) min = i;
    }
    dominants[clusIndex] = min;
}

void Swarm::spsa(int clusIndex)
{
    int domIndex = dominants[clusIndex];
    
    //double a = 1, c = 1, A = 1, alpha = 1, gamma = 1/6;
    double a = 1, c = 1, A = 1, alpha = 0.602, gamma = 0.101;
    vector<double> ghat;
    // Execute SPSA ALGORITHM
    for(int k = 1; k <= 50; ++k)
    {
        // bernoulli deltaK
        std::default_random_engine generator;
        std::bernoulli_distribution distribution(0.5);
        vector<double> ghat; // deltaK
        ghat.resize(Nvariables);
        for (int i = 0; i < Nvariables; ++i) {
            if (distribution(generator)) ghat[i] = 1;
            else ghat[i] = -1;
        } 
        
        // ak and ck
        double ak = a / pow(A + k, alpha);
        double ck = c / pow(k, gamma);
        
        // evaluate positive
        vector<double> pos;
        pos.resize(Nvariables);
        for(int i = 0; i < Nvariables; ++i)
        {
            pos[i] = P[domIndex].x[i] + (ck * ghat[i]);
        }
        double positive = fitness(pos);

        // evaluate negative
        vector<double> neg;
        neg.resize(Nvariables);
        for(int i = 0; i < Nvariables; ++i)
        {
            neg[i] = P[domIndex].x[i] - (ck * ghat[i]);
        }
        double negative = fitness(neg);

        // calculate ghat
        double coff = (positive - negative) / (2 * ck);
        for(int i = 0; i < Nvariables; ++i)
        {
            ghat[i] = coff / ghat[i];
        }
        
        // P[domIndex].x -= ghat; 
        for(int i = 0; i < Nvariables; ++i)
        {
            P[domIndex].x[i] -= ak * ghat[i];
        }
        
        // evaluate new x and record pbest
        evaluate(domIndex);
        if (P[domIndex].fitness < P[domIndex].pBest)
        {
            P[domIndex].pBest = P[domIndex].fitness;
            P[domIndex].xBest = P[domIndex].x;
            dominants[clusIndex] = domIndex;
            candidates[clusIndex] = P[domIndex];
        }
    }
    //nfc++;
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
    //cout << gBestValue << endl;

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

    for(int i = 0; i < K; ++i)
    {
        candidates[i] = P[gBestIndex];
    }

    printf("0 : ");
    printf(" = %g\n", gBestValue);

    ofstream file;
    file.open("ClusterPSO.txt", ios::out | ios::app);
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

void Swarm::calculateVMax() 
{
    double xmin[Nparticles], xmax[Nparticles];

    for (int d = 0; d < Nvariables; d++) 
    {
        xmin[d] = xmax[d] = P[0].x[d];
        for (int n = 1; n < Nparticles; n++) 
        {
            double pos = P[n].x[d];
            if (pos < xmin[d]) xmin[d] = pos;
            if (pos > xmax[d]) xmax[d] = pos;
        }
        maxV[d] = xmax[d] - xmin[d];
    }
}

void Swarm::particleMovement(int clusIndex)
{
    for (int n = 0; n < Nparticles; n++)
    {
        // update velocities
        if (cluster[n] == clusIndex && n != dominants[clusIndex])
        {
            //cout << n << endl;
            for (int d = 0; d < Nvariables; d++)
            {
                double r1 = Rand();
                double r2 = Rand();
                double sk = Rand();
                P[n].v[d] = w * P[n].v[d] 
                        + c1 * r1 * (P[n].xBest[d] - P[n].x[d]) 
                        + c2 * r2 * (P[gBestIndex].x[d] - P[n].x[d]); 
                        + sk * (P[n].x[d] - P[dominants[clusIndex]].x[d]);
                // check v with its dimensional maxV
                if (P[n].v[d] > maxV[d])
                    P[n].v[d] = maxV[d];
                else if (P[n].v[d] < -maxV[d])
                    P[n].v[d] = -maxV[d];
            }
            // update positions
            //cout << n << " : ";
            for (int d = 0; d < Nvariables; d++)
            {
                P[n].x[d] += P[n].v[d];
                //cout << P[n].x[d] << " ";
            }
            //cout << "\n";
        }
    }
    //cout << "\n";
}

void Swarm::evaluateSwarm(int clusIndex) 
{
    for(int i = 0; i < Nparticles; i++) {
        if(cluster[i] == clusIndex) // i != dominants[clusIndex]
        {
            evaluate(i);    
            nfc++;

            if (nfc % 5000 == 0) 
            {
                printf("%d : ", nfc);
                printf(" = %g\n", gBestValue);

                ofstream file;
                file.open("ClusterPSO.txt", ios::out | ios::app);
                file << nfc << "," << gBestValue << "\n";
                file.close();
            }
        }       
    }

    for(int n = 0; n < Nparticles; n++) {
        if(cluster[n] == clusIndex) // n != dominants[clusIndex]
        {
            if (P[n].fitness < P[n].pBest)
            {
                P[n].pBest = P[n].fitness;
                P[n].xBest = P[n].x;

                // if (P[n].fitness < gBestValue)
                // {
                //     gBestIndex = n;
                //     gBestValue = P[n].fitness;
                // }

                if (P[n].fitness < P[dominants[clusIndex]].fitness && P[n].fitness < candidates[clusIndex].fitness)
                {
                    dominants[clusIndex] = n;
                    candidates[clusIndex] = P[n];
                }
            }
        }
    }
}

void Swarm::updateBest(int index)
{
    P[index].xBest = P[index].x;
    P[index].pBest = P[index].fitness;
}

void Swarm::updateGBest()
{
    gBestValue = candidates[0].pBest;
    gBestIndex = dominants[0];
    for(int i = 1; i < K; ++i)
    {
        if(candidates[i].pBest < gBestValue)
        {
            gBestValue = candidates[i].pBest;
            gBestIndex = dominants[i];
        }
    }
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
}

vector<size_t> Swarm::k_means(const vector<Particle> &particles, size_t k)
{
    static random_device seed;
    static mt19937 random_number_generator(seed());
    uniform_int_distribution<size_t> indices(0, Nparticles - 1);

    vector<Particle> means(k); // cluster centroid particle
    for(auto& cluster : means)
    {
        int center = (int) Nparticles * Rand();
        cluster = particles[center];
    }

    // Compute E(t)
    // double currentE;
    // double nextE = NULL;    
    
    vector<size_t> assignments(particles.size()); // which cluster that particle belong to
    //do
    for(size_t iteration = 0; iteration < Niterations; ++iteration) 
    {
        // currentE = nextE;

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
            // Turn 0/0 into 0/1 to avoid zero division.
            const auto count = max<size_t>(1, counts[cluster]);
            for(size_t var = 0; var < new_means[cluster].x.size(); ++var)
            {
                means[cluster].x[var] += new_means[cluster].x[var] / count;
            }
            // Compute New Centroid fitness
            means[cluster].fitness = fitness(means[cluster].x);
        }
    }

    //     // Compute E(t+1)
    //     nextE = 0;
    //     double sum[k];
    //     for (size_t particle = 0; particle < particles.size(); ++particle)
    //     {
    //         const auto cluster = assignments[particle];
    //         for (size_t var = 0; var < particles[particle].x.size(); ++var)
    //         {
    //             sum[cluster] += square(particles[particle].x[var] - means[cluster].x[var]);
    //         }
    //     }
    //     for (size_t cluster = 0; cluster < k; cluster++)
    //     {
    //         nextE += sum[cluster];
    //     }
    //     cout << nextE << endl;

    // } while (currentE != nextE);

    return assignments;
}

int main(int argc, const char *argv[])
{
    char probName[8][14] = {"Sphere", "Schwefel", "Rosenbrock", "QuarticNoise", "Ackley", "Griewank", "Rastrigin", "Penalized"};
    for(::probNo = 0; ::probNo < 8; ::probNo++)
    {
        ofstream file;
        file.open("ClusterPSO.txt", ios::out | ios::app);
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
