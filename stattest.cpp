#include "stattest.h"

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <utility>

#include <cmath>

using std::cout;
using std::lgamma;
using std::vector;
using std::pair;

void pdp(vector<pair<double,size_t>> &v){
    for (const auto &c : v){
        std::cout << c.first << ", " << c.second << "\n";
    }
    std::cout << "\n";
}
template<typename D>
void pdp(vector<D> &v){
    for (D c : v){
        std::cout << c << " ";
    }
    std::cout << "\n";
}

double lbinom( int n, int k){
    if ( n == 0 || k == 0){
        return 0;
    }
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1); 
}

double lhyper_geom_pmf( int x, int n, int m, int N){
    return lbinom( n, x) + lbinom(N-n, m-x) - lbinom( N, m );
}

double hyper_geom_pmf( int x, int n, int m, int N){
    return exp(lhyper_geom_pmf(x,n,m,N));
}

//Upper
double hyper_geom_cdf( int x, int n, int m, int N){
    double cdf = 0;
    for( int i = x; i <= n; ++i){
        cdf += hyper_geom_pmf(i,n,m,N);
    }
    return cdf;
}


enum class pvalue_correction{
    BH,BY,
};

class hypothesis_testing{
    public:
    vector<double>  corr_pvals;
    vector<bool>    null_rejected;
};

// Used statsmodel as a reference
// https://github.com/statsmodels/statsmodels/blob/main/statsmodels/stats/multitest.py
//
hypothesis_testing benjamini_hochberg_fdr_correction(const vector<double> &pvalues, double alpha){

    hypothesis_testing test;
    vector<pair<double,size_t>> paired_pvalues;
    for(int i = 0; i < pvalues.size(); ++i){
        paired_pvalues.push_back(std::make_pair(pvalues[i],i));
    }
    
    //Sort by the p-values
    std::sort(paired_pvalues.begin(),paired_pvalues.end(),
            [] (const pair<double,int> &a, const pair<double,int> &b)
                { return a.first > b.first;});
  
    int k;

    double min_p = 1;//paired_pvalues[0].first / ecdf_vector[0];

    for ( k = 0; k < paired_pvalues.size(); ++k){
        double q = alpha * (k+1) /  paired_pvalues.size();
        //if(q > paired_pvalues[k].first){
        double new_p =  paired_pvalues[k].first / (k+1) * paired_pvalues.size();
        if ( new_p < min_p){
            min_p = new_p;
        }
        paired_pvalues[k].first = min_p;

    }
    //Sort p-values back to original order
    std::sort(paired_pvalues.begin(),paired_pvalues.end(),
            [] (const pair<double,int> &a, const pair<double,int> &b)
                { return a.second < b.second;});
   

    for( const auto &pr : paired_pvalues){
        test.corr_pvals.push_back(pr.first);
        test.null_rejected.push_back( pr.first < alpha);
    }

    return test;
}

hypothesis_testing benjamini_yekutieli_fdr_correction(const vector<double> &pvalues, double alpha){
    hypothesis_testing test;
    return test;
}
hypothesis_testing correct_pvalues(const vector<double> &pvalues, double alpha, pvalue_correction method){

    switch(method){
        case pvalue_correction::BH:
            return benjamini_hochberg_fdr_correction(pvalues, alpha);
        case pvalue_correction::BY:
            return benjamini_yekutieli_fdr_correction(pvalues, alpha);
        default:
            throw std::invalid_argument("Correction Method Not implemented!");
    }
}

int main(int argc, char **argv){

    std::vector<double> pvals{{0.074,0.205,0.041,0.039,0.001,0.042,0.060,0.008,0.05,0.049,0.025}};

    auto cp = correct_pvalues( pvals, 0.05, pvalue_correction::BH);
    pdp(pvals);

    pdp(cp.corr_pvals);
    pdp(cp.null_rejected);
    return 0;
}
