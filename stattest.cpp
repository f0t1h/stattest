
#include "stattest.h"

#include <iostream>
using std::cout;

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
int main(int argc, char **argv){
    cout << argc << "\t" << argv[0] << "\n";
    std::vector<double> pvals{{0.074,0.205,0.041,0.039,0.001,0.042,0.060,0.004,0.05,0.049,0.025}};


    pdp(pvals);
    auto cp = multiple_test( pvals, 0.05, pvalue_corrector::BENJAMINI_HOCHBERG);
    pdp(cp.corr_pvals);
    pdp(cp.null_rejected);
    cp = multiple_test( pvals, 0.05, pvalue_corrector::BENJAMINI_YEKUTIELI);
    pdp(cp.corr_pvals);
    pdp(cp.null_rejected);
    cout << hyper_geom_cdf( 30, 50, 60, 100 ) << "|\n";

    return 0;
}
