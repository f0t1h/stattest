#include "stattest.h"

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <utility>
#include <tuple>
#include <cmath>

using std::cout;
using std::lgamma;
using std::vector;
using std::pair;
using std::tuple;
using std::get;

void pdp(vector<pair<double,size_t>> &v){
    for (const auto &c : v){
        cout << c.first << ", " << c.second << "\n";
    }
    cout << "\n";
}
template<typename D>
void pdp(vector<D> &v){
    for (D c : v){
        cout << c << " ";
    }
    cout << "\n";
}

double approx_digamma(double n){
    return std::log(n) - 1.0 / (2 * n);
}

// I can convert this to a look-up table instead.
double approx_harmonic(int n){
    return approx_digamma(n + 1) + 0.5772156649;
}

double lbinom( int n, int k){
    if ( n == 0 || k == 0){
        return 0;
    }
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1); 
}

double lhyper_geom_pmf( int x, int n, int m, int N){
    return lbinom(n, x) + lbinom(N-n, m-x) - lbinom(N, m);
}

double hyper_geom_pmf( int x, int n, int m, int N){
    return exp(lhyper_geom_pmf(x,n,m,N));
}

//Upper
//Computes upper tail of fisher_exact test
double hyper_geom_cdf( int x, int n, int m, int N){
    double cdf = 0;
    for( int i = x; i <= n; ++i){
        cdf += hyper_geom_pmf(i,n,m,N);
    }
    return cdf;
}

enum class pvalue_corrector{
    BENJAMINI_HOCHBERG,
    BENJAMINI_YEKUTIELI,
    BONFERRONI,
    HOLM_BONFERRONI,
    SIDAK,
};

class hypothesis_testing{
    public:
    vector<double>  corr_pvals;
    vector<bool>    null_rejected;
};

// Used statsmodel as a reference
// https://github.com/statsmodels/statsmodels/blob/main/statsmodels/stats/multitest.py
//
hypothesis_testing fdr_correction(const vector<double> &pvalues, double alpha, pvalue_corrector method){

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


    double new_p = 1;
    double m = paired_pvalues.size();
    for ( k = 0; k < paired_pvalues.size(); ++k){
//        double q = alpha * (k+1) /  paired_pvalues.size();
        double p = paired_pvalues[k].first;
        switch(method){
            case pvalue_corrector::BENJAMINI_HOCHBERG:
                p = p * m / (1.0 * k+1);
                break;
            case pvalue_corrector::BENJAMINI_YEKUTIELI:
                p = p * approx_harmonic(k+1) * m / (1.0 * k+1);
                break;
        }
        paired_pvalues[k].first = p;
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


hypothesis_testing fixed_alpha_correction(const vector<double> &pvalues, double alpha, pvalue_corrector method){
    hypothesis_testing test;
    double corrected_alpha = alpha;
    switch(method){
        case pvalue_corrector::BONFERRONI:
            corrected_alpha = alpha / pvalues.size();
            break;
        case pvalue_corrector::SIDAK:
            corrected_alpha = 1 - pow((1 - alpha),1.0/pvalues.size());
            break;
        default:
            throw std::invalid_argument("Method not defined");
    }

    for( double p : pvalues){
        double corrected_p = p;

        switch(method){
            case pvalue_corrector::BONFERRONI:
                p = p * pvalues.size();
                break;
            case pvalue_corrector::SIDAK:
                p = 1 - pow((1.0 - p),pvalues.size());
                alpha = 1 - pow((1 - alpha),1.0/pvalues.size());
                break;
            default:
                throw std::invalid_argument("Method not defined");
        }

        test.corr_pvals.push_back(corrected_p);
        test.null_rejected.push_back( p < corrected_alpha);       
    }
    return test;
}
hypothesis_testing holm_bonferroni_method(const vector<double> &pvalues, double alpha){
    hypothesis_testing test;
    
    vector<tuple<double,size_t,bool>> paired_pvalues;
    for(int i = 0; i < pvalues.size(); ++i){
        paired_pvalues.push_back(std::make_tuple(pvalues[i],i,false));
    }
    
    //Sort by the p-values
    std::sort(paired_pvalues.begin(),paired_pvalues.end(),
            [] (const tuple<double,int,bool> &a, const tuple<double,int,bool> &b)
                { return get<0>(a) < get<0>(b);});
 
    int cnt = 0;
    size_t m = paired_pvalues.size();
    for( auto &pr : paired_pvalues){
        double alpha_prime = alpha / (m - cnt);
        cnt+=1;
        if (get<0>(pr) < alpha_prime){
            get<2>(pr) = true;
        }
    }
    //Sort p-values back to original order
    std::sort(paired_pvalues.begin(),paired_pvalues.end(),
            [] (const tuple<double,int,bool> &a, const tuple<double,int,bool> &b)
                { return get<1>(a) < get<1>(b);});

    for( auto &pr : paired_pvalues){
        test.corr_pvals.push_back(get<0>(pr));
        test.null_rejected.push_back(get<2>(pr));       
    }
    return test;
}

hypothesis_testing multiple_test(const vector<double> &pvalues, double alpha, pvalue_corrector method){
    switch(method){
        case pvalue_corrector::BENJAMINI_HOCHBERG:
        case pvalue_corrector::BENJAMINI_YEKUTIELI:
            return fdr_correction(pvalues, alpha, method);
        case pvalue_corrector::BONFERRONI:
        case pvalue_corrector::SIDAK:
            return fixed_alpha_correction(pvalues, alpha, method);
        case pvalue_corrector::HOLM_BONFERRONI:
            return holm_bonferroni_method(pvalues, alpha);
        default:
            throw std::invalid_argument("Correction Method Not implemented!");
    }
}

int main(int argc, char **argv){

    std::vector<double> pvals{{0.074,0.205,0.041,0.039,0.001,0.042,0.060,0.004,0.05,0.049,0.025}};


    pdp(pvals);
    auto cp = multiple_test( pvals, 0.05, pvalue_corrector::BENJAMINI_HOCHBERG);
    pdp(cp.corr_pvals);
    pdp(cp.null_rejected);
    cp = multiple_test( pvals, 0.05, pvalue_corrector::BENJAMINI_YEKUTIELI);
    pdp(cp.corr_pvals);
    pdp(cp.null_rejected);
    cout << hyper_geom_cdf( 30, 50, 60, 100 ) << "|\n";
    //for(int i = 1; i< 50; ++i){
    //    cout << "H_i" << i << " = " <<approx_harmonic(i) << "\n";
    //}
    return 0;
}
