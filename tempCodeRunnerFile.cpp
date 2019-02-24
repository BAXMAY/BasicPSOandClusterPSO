for(size_t iteration = 0; iteration < Niterations; ++iteration) 
    // {
    //     //if (nextE != NULL) 
    //     //{
    //     //    currentE = nextE;
    //     //}

    //     // Find assignments.
    //     for (size_t particle = 0; particle < particles.size(); ++particle)
    //     {
    //         auto best_distance = std::numeric_limits<double>::max();
    //         size_t best_cluster = 0;
    //         for (size_t cluster = 0; cluster < k; ++cluster)
    //         {
    //             const double distance = eucildean_distance(particles[particle], means[cluster]); //distance
    //             if (distance < best_distance)
    //             {
    //                 best_distance = distance;
    //                 best_cluster = cluster;
    //             }
    //         }
    //         assignments[particle] = best_cluster;
    //     }

    //     // Sum up and count particles for each cluster.
    //     vector<Particle> new_means(k);
    //     vector<size_t> counts(k, 0); // number of particle in cluster, use for find new mean
    //     for (size_t particle = 0; particle < particles.size(); ++particle)
    //     {
    //         const auto cluster = assignments[particle];
    //         for(size_t var = 0; var < particles[particle].x.size(); ++var)
    //         {
    //             new_means[cluster].x[var] += particles[particle].x[var];
    //         }
    //         counts[cluster] += 1;
    //     }

    //     // Divide sums by counts to get new centroids.
    //     for (size_t cluster = 0; cluster < k; ++cluster)
    //     {
    //         double fitness = 0, sum = 0;
    //         // Turn 0/0 into 0/1 to avoid zero division.
    //         const auto count = max<size_t>(1, counts[cluster]);
    //         for(size_t var = 0; var < new_means[cluster].x.size(); ++var)
    //         {
    //             means[cluster].x[var] += new_means[cluster].x[var] / count;
    //             // Compute New Centroid fitness
    //             for (int j = 0; j < var; ++j)
    //             {
    //                 sum += means[cluster].x[j];
    //             }
    //             fitness += sum * sum;
    //         }
    //         means[cluster].fitness = fitness;
    //     }
    // }