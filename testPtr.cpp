#include <iostream>
#include <cstdlib>
using namespace std;

class Particle {
    public:
        double *x;
        double *v;
        double fitness;
        double pBest;
        double *xBest;
        Particle() {
            x = (double *)malloc(sizeof(double) * 30);
            if (x == NULL)
            {
                fprintf(stderr, "Cannot allocate memory for %d x\n", 30);
                exit(1);
            }
            v = (double *)malloc(sizeof(double) * 30);
            if (v == NULL)
            {
                fprintf(stderr, "Cannot allocate memory for %d v\n", 30);
                exit(1);
            }
            xBest = (double *)malloc(sizeof(double) * 30);
            if (xBest == NULL)
            {
                fprintf(stderr, "Cannot allocate memory for %d xBest\n", 30);
                exit(1);
            }
        }
};

int main() {
    Particle p1, p2;
    p1.x[0] = 2;
    p2.x[0] = 3;
    cout << p1.x << endl;
    cout << p1.x[0] << endl;
    cout << p2.x << endl;
    cout << p2.x[0] << endl;
    *(p1.x + 1) = 9;
    cout << p1.x + 1 << endl;
    cout << p1.x[1] << endl;
    cout << p2.x << endl;
    cout << p2.x[0] << endl;

    return 0;
}