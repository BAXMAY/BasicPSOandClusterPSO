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