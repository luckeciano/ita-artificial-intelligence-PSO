#include <iostream>
#include <ctime>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;



#define BEST_VALUE_EIGHT_QUEENS 55
int computeCostFunction (vector<int> queens) {
    int res = 0;
    for (int i = 0; i < queens.size(); i++) {
        for (int j = i + 1; j < queens.size(); j++) {
            if (queens[i] != queens[j] && abs(i - j) != abs(queens[i] - queens[j]))
                res++;
        }
    }
    return res;
}

class Particle {
public:
    Particle(vector<int> p, vector<double> v) {
        actualParticle = p;
        prevBestParticle = p;
        velocity = v;
        prevBestParticleFitness = -1;
    }
    vector<int> getActualParticle() {
        return actualParticle;
    }
    vector<int> getPrevBestParticle() {
        return prevBestParticle;
    }
    int getPrevPartFitnessFunc() {
        if (prevBestParticleFitness != -1) return prevBestParticleFitness;
        return computeCostFunction(prevBestParticle);
    }
    void setParticle(int pos, int value) {
        actualParticle[pos] = value;
    }
    void setBestParticle (vector<int> p, int newFitValue) {
        prevBestParticle = p;
        prevBestParticleFitness = newFitValue;
    }
    void setVelocity (int pos, double value) {
        velocity[pos] = value;
    }
    double getVelocity(int pos) {
        return velocity[pos];
    }

private:
    vector<int> actualParticle;
    vector<int> prevBestParticle;
    vector<double> velocity;
    int prevBestParticleFitness;
};
vector<Particle> randomInitialization(int N, int nParticles, int vMax) {
    vector<Particle> result;
    for (int i = 0; i < nParticles; i++) {
        vector<int> newParticle;
        vector<double> newVelocity;
        for (int j = 0; j < N; j++) {
            newParticle.push_back(rand()%N);
            newVelocity.push_back(-vMax + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2*vMax))));

        }
        Particle p(newParticle, newVelocity);
        result.push_back(p);
    }
    return result;
}

vector<int> computeBestParticle(vector<Particle> p) {
    int bestFunc = -1;
    vector<int> res;
    for (int i = 0; i < p.size(); i++) {
        int aux = computeCostFunction(p[i].getActualParticle());
        if (bestFunc < aux) {
            aux = bestFunc;
            res = p[i].getActualParticle();
        }
    }
    return res;
}
void printParticle (vector<int> p) {
    for (int i = 0; i < p.size(); i++)
        cout << p[i] << " ";
    cout << endl;
}


vector<int> PSO (int N, double w, double vMax, double phiOne, double phiTwo, double C1, double C2, int nParticles) {
    vector<Particle> particles = randomInitialization(N, nParticles, vMax);;
    vector<int> bestParticle = computeBestParticle(particles);
    printParticle(bestParticle);
    int bestParticleFitnessValue = computeCostFunction(bestParticle);
    while (bestParticleFitnessValue != BEST_VALUE_EIGHT_QUEENS) {
        for (int i = 0; i < nParticles; i++) {
            int fitnessFunc = computeCostFunction(particles[i].getActualParticle());
            if (fitnessFunc > particles[i].getPrevPartFitnessFunc()) {
                particles[i].setBestParticle(particles[i].getActualParticle(), fitnessFunc);
            }
            if (fitnessFunc > bestParticleFitnessValue) {
                bestParticle = particles[i].getActualParticle();
                bestParticleFitnessValue = fitnessFunc;
            }
            for (int j = 0; j < N; j++) {
                particles[i].setVelocity(j, w * particles[i].getVelocity(j) +
                                         C1 * phiOne * (particles[i].getPrevBestParticle()[j] - particles[i].getActualParticle()[j]) +
                                         C2 * phiTwo * (bestParticle[j] - particles[i].getActualParticle()[j]));
                particles[i].setParticle(j, int(fabs(particles[i].getActualParticle()[j] + particles[i].getVelocity(j)))%N);
            }

        }
    }
    return bestParticle;

}

int main() {
    const clock_t begin_time = clock();
    //PARAMETERS
    int N = 11;
    double w = 1;
    double vMax = 3;
    double phiOne = static_cast <float> (rand()) /( static_cast <float> (RAND_MAX));
    double phiTwo = static_cast <float> (rand()) /( static_cast <float> (RAND_MAX));
    double C1 = 5;
    double C2 = 2;
    int particles = 20;
    vector<int> result = PSO(N, w, vMax, phiOne, phiTwo, C1, C2, particles);
    const clock_t end_time = clock();
    cout << "Final Result From PSO:" << endl;
    for (int i = 0; i < result.size(); i++) {
        cout << result[i] << " ";
    }
    cout << "Time: " << double((end_time - begin_time)/double(CLOCKS_PER_SEC)) << " seconds." << endl;
    system("Pause");
    return 0;
}
