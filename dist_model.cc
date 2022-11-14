#include <cstdio>
#include <cstdlib>
#include <cmath>

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
    unsigned int i,j,k,b[n],bmax=0;

    // Allocate the GSL random number generator
    gsl_rng *rng;
    rng=gsl_rng_alloc(gsl_rng_taus2);

    // Clear the bins for frequency distribution
    for(j=0;j<n;j++) b[j]=0;

    // Generate samples
    for(i=0;i<100000000;i++) {
        k=gsl_ran_poisson(rng,t_power(2.16,0.0111,rng));
        
        j=l_scale(k);
        if(j>=n) {
            fprintf(stderr,"Count %u and bin number %u out of range\n",k,j);
        } else {
            if(j>bmax) bmax=j;
            b[j]++;
        }
    }

    for(j=0;j<=bmax;j++) printf("%g %u %g\n",il_scale(j),b[j],b[j]/double(l_size(j)));
}
