#include <cstdio>
#include <cstdlib>
#include <cmath>

#ifndef __OPENMP
#include <omp.h>
inline int t_num() {return omp_get_thread_num();}
#else
inline int t_num() {return 0;}
#endif

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Parameter setting the transition in the scaling function from linear to log
const int lam=10;
const double flam=static_cast<double>(lam);

// Maximum number of bins
const unsigned int n=1000;

// Scaling function used for binning data
inline unsigned int l_scale(unsigned int x) {
    return x<=lam?x:static_cast<unsigned int>(flam*(1.+log(x/flam))+0.5);
}

// Inverse scaling function
inline double il_scale(double x) {
    return x<=flam?x:flam*exp(x/flam-1);
}

// Bin size in the scaling function
inline unsigned int l_size(unsigned int x) {
    return x<lam?1:int(il_scale(x+0.5))-int(il_scale(x-0.5));
}

// Truncated power law distribution
inline double t_power(double a,double b,gsl_rng *rng) {
    return b*pow(gsl_rng_uniform_pos(rng),1/(1-a));
}

int main() {
    unsigned long b[n];
    for(unsigned int l=0;l<n;l++) b[l]=0;

    // Generate samples
#pragma omp parallel
    {
        unsigned int j,k;

        // Allocate the GSL random number generator
        gsl_rng *rng;
        rng=gsl_rng_alloc(gsl_rng_taus2);
        gsl_rng_set(rng,t_num()+1);
    
        // Allocate local bins
        unsigned long c[n];
        for(j=0;j<n;j++) c[j]=0;

#pragma omp for
        for(unsigned long i=0;i<100000000;i++) {
            double v=t_power(1.2,0.0001,rng)*t_power(2.16,0.0111,rng);
            if(v>1000000) continue;
            k=gsl_ran_poisson(rng,v);

            j=l_scale(k);
            if(j>=n) {
                fprintf(stderr,"Count %u and bin number %u out of range\n",k,j);
            } else c[j]++;
        }

        // Add local bins to global bins
        for(j=0;j<n;j++) {
#pragma omp atomic
            b[j]+=c[j];
        }
        gsl_rng_free(rng);
    }

    unsigned int bmax=n;
    while(b[--bmax]==0) {
        if(bmax==0) {fputs("No counts\n",stderr);return 1;}
    }

    for(unsigned int l=0;l<=bmax;l++) printf("%g %lu %g\n",il_scale(l),b[l],b[l]/double(l_size(l)));
}
