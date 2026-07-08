// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// Encodes the full state (l, eHQ, iHQ, eHG, iHG, eP, iP, n) as a string
inline std::string state_key(int l, int eHQ, int iHQ, int eHG, int iHG, int eP, int iP, int n) {
    return std::to_string(l)   + "_" + std::to_string(eHQ) + "_" +
           std::to_string(iHQ) + "_" + std::to_string(eHG) + "_" +
           std::to_string(iHG) + "_" + std::to_string(eP)  + "_" +
           std::to_string(iP)  + "_" + std::to_string(n);
}

// [[Rcpp::export]]
NumericVector fast_trial_cpp_P(int nmaxmax, double gamma, double rC,
                             String cohort, double beta, double zeta, double eps) {

    // Population setup (cohorting types)
    int NHQ = 0, NHG = 0, NP = 12;
    if      (cohort == "S") { NHQ = 2;  NHG = 11; }
    else if (cohort == "M") { NHQ = 6;  NHG = 7;  }
    else if (cohort == "N") { NHQ = 13; NHG = 0;  }

    // Parameters  
    double deltaP = 1.0 / 6.0;
    // Hypoexponential rates for isolation stages (index 0..4)
    std::vector<double> Hyp = {0.0, 6.69470955813761, 0.555473683554152,
                                0.55547824700609, 0.555478430962926};

    double rI = 0.95;

    // Beta_O[l] for l = 0,...,K
    std::vector<double> Beta_O(Hyp.size(), 0.0);
    if (NHQ > 0) {
        double bo = (1.0 - rI) * beta * 2.0 / NHQ; 
        for (int i = 1; i <= 4; i++) Beta_O[i] = bo; //Fixed infection rate (but can vary with removal stage)
    }

    double Beta_H   = eps * (2.0 / 3.0) * beta;   // cross-cohort HCW-to-HCW (we set eps=1 or 1-rC in the paper)
    double Beta_HQ  = (2.0 / 3.0) * beta;   // within HQ cohort
    double Beta_HG  = (2.0 / 3.0) * beta;   // within HG cohort
    double Beta_HQP = (1.0 - rC) * beta;    // between HQ and patient
    double Beta_HGP = beta;                  // between HG and patient
    double Beta_P = 0;                       // within patient group (but we assume separate rooms, hence Beta_P = 0)

    double zeta_HQ = zeta, zeta_HG = zeta, zeta_P = zeta;
    double gamma_HQ = gamma, gamma_HG = gamma, gamma_P = gamma;

    //  Initial state
    int l_init   = Hyp.size() - 1; // l_init = 4
    int eHQ_init = 0, iHQ_init = 0;
    int eHG_init = 0, iHG_init = 0;
    int eP_init  = 0, iP_init  = 0;

    int q_init = eHQ_init + iHQ_init;  // total HQ infected at start = 0
    int g_init = eHG_init + iHG_init;  // total HG infected at start = 0


    // Calculating outbreak probability weights
    int K = l_init;
    NumericVector outbreak_weights(K + 1); // Weights for l = 1..K
    double total_po_sum = 0.0;
    for (int l = 1; l <= K; l++) {
    double prod = 1.0;
    // Probability of NO exposure in stages higher than l
    for (int q = l + 1; q <= K; q++) {
        double rate_HQ = NHQ * Beta_O[q];
        prod *= Hyp[q] / (Hyp[q] + rate_HQ);
    }

    // Probability of exposure happening IN stage l
    double rate_HQ_l = NHQ * Beta_O[l];
    prod *= rate_HQ_l / (rate_HQ_l + Hyp[l]);
    
    outbreak_weights[l] = prod;
    total_po_sum += prod;
    }
     // Normalize weights so they sum to 1 (Conditioning on an outbreak)
    for (int l = 1; l <= K; l++) {
    outbreak_weights[l] /= total_po_sum;
    }
    
    std::unordered_map<std::string, double> VALUE_STORE;
    VALUE_STORE.reserve(1 << 20); // reserve 2^20 entries
    
    NumericVector p_set;
    double P_val  = 0.0;
    

    // =========================================================
    // MAIN LOOP
    //   for nmax in 0:nmaxmax
    //     for n in 0:nmax
    //       for l in 0:l_init
    //         for h in (NHQ+NHG):h_init   [for level h_init = q_init+g_init]
    //           qlow(h) = max(0, h - NHG)
    //           qLB(h) = max(qlow(h), q_init)
    //           qupp(h) = min(NHQ, h)
    //           for q in qupp(h):qLB(h) 
    //              for iHQ in q:iHQ_init
    //                  for iHG in h-q:iHG_init
    //                      ePmax = min (eP_init+nmax-n, NP)
    //                      for eP in 0:ePmax
    //                          ipmax = min(iP_init+ePmax-eP, NP-eP)
    //                          for iP in 0:ipmax
    //                                solve linear equation shown in paper

for (int nmax = 0; nmax <= nmaxmax; nmax++) {
        VALUE_STORE.clear();

        for (int n = 0; n <= nmax; n++) {
            for (int l_ = 0; l_ <= l_init; l_++) {
                for (int h = NHQ+NHG; h >= q_init+g_init; h--){
                    int qlow = std::max(0, h - NHG);
                    int qLB = std::max(qlow, q_init);
                    int qupp = std::min(NHQ, h);
                    for (int q = qupp; q >= qLB; q--) {
                        for (int iHQ = q; iHQ >= iHQ_init; iHQ--) {
                            for (int iHG = h - q; iHG >= iHG_init; iHG--) {
                            int ePmax = std::min(eP_init + nmax - n, NP);
                            for (int eP = 0; eP <= ePmax; eP++) {
                                int iPmax = std::min(iP_init + ePmax - eP, NP - eP);
                                    for (int iP = 0; iP <= iPmax; iP++){
                                      int eHQ = q-iHQ;
                                      int eHG = (h-q) - iHG;
                                      double Gamma = gamma_HQ * iHQ + gamma_HG * iHG + gamma_P * iP;
                                      double lambdaHQ = (NHQ - eHQ - iHQ) * (Beta_O[l_] + iHQ * Beta_HQ + iHG * Beta_H + iP * Beta_HQP);
                                      double lambdaHG = (NHG - eHG - iHG) * (iHG * Beta_HG + iHQ * Beta_H + iP * Beta_HGP);
                                      double lambdaP = (NP - eP - iP) * (iHQ * Beta_HQP + iHG * Beta_HGP);
                                      double etaHQ = eHQ * zeta_HQ;
                                      double etaHG = eHG * zeta_HG;
                                      double etaP = eP * zeta_P;
                                      double delta_iso = Hyp[l_];
                                      double mu_IP = iP * deltaP;
                                      double mu_EP = eP * deltaP;

                                      double Theta = Gamma + lambdaHQ + lambdaHG + lambdaP + etaHQ + etaHG + etaP + delta_iso + mu_IP + mu_EP;
 


                                        if (n==0 && l_==0 && h==0 && q==0 && iHQ==0 && iHG==0 && eP==0 && iP==0) {
                                            VALUE_STORE[state_key(l_, eHQ, iHQ, eHG, iHG, eP, iP, n)] = 1.0;
                                        } else {
                                            VALUE_STORE[state_key(l_, eHQ, iHQ, eHG, iHG, eP, iP, n)] = 0.0;
                                        }
                                        
                                        if (Theta > 0){

                                        if (n==0){
                                            VALUE_STORE[state_key(l_, eHQ, iHQ, eHG, iHG, eP, iP, n)] = 
                                            (Gamma +
                                            lambdaHQ*VALUE_STORE[state_key(l_, eHQ+1, iHQ, eHG, iHG, eP, iP, n)] +
                                            lambdaHG*VALUE_STORE[state_key(l_, eHQ, iHQ, eHG+1, iHG, eP, iP, n)] +
                                            delta_iso*VALUE_STORE[state_key(l_-1, eHQ, iHQ, eHG, iHG, eP, iP, n)] +
                                            etaHQ*VALUE_STORE[state_key(l_, eHQ-1, iHQ+1, eHG, iHG, eP, iP, n)] +
                                            etaHG*VALUE_STORE[state_key(l_, eHQ, iHQ, eHG-1, iHG+1, eP, iP, n)] +
                                            etaP*VALUE_STORE[state_key(l_, eHQ, iHQ, eHG, iHG, eP-1, iP+1, n)] +
                                            mu_IP*VALUE_STORE[state_key(l_, eHQ, iHQ, eHG, iHG, eP, iP-1, n)] +
                                            mu_EP*VALUE_STORE[state_key(l_, eHQ, iHQ, eHG,  iHG, eP-1, iP, n)])/Theta;
                                        } else {
                                            VALUE_STORE[state_key(l_, eHQ, iHQ, eHG, iHG, eP, iP, n)] =
                                            (lambdaHQ*VALUE_STORE[state_key(l_, eHQ+1, iHQ, eHG, iHG, eP, iP, n)] +
                                            lambdaHG*VALUE_STORE[state_key(l_, eHQ, iHQ, eHG+1, iHG, eP, iP, n)] +
                                            lambdaP*VALUE_STORE[state_key(l_, eHQ, iHQ, eHG, iHG, eP+1, iP, n-1)] +
                                            delta_iso*VALUE_STORE[state_key(l_-1, eHQ, iHQ, eHG, iHG, eP, iP, n)] +
                                            etaHQ*VALUE_STORE[state_key(l_, eHQ-1, iHQ+1, eHG, iHG, eP, iP, n)] +
                                            etaHG*VALUE_STORE[state_key(l_, eHQ, iHQ, eHG-1, iHG+1, eP, iP, n)] +
                                            etaP*VALUE_STORE[state_key(l_, eHQ, iHQ, eHG, iHG, eP-1, iP+1, n)] +
                                            mu_IP*VALUE_STORE[state_key(l_, eHQ, iHQ, eHG, iHG, eP, iP-1, n)] +
                                            mu_EP*VALUE_STORE[state_key(l_, eHQ, iHQ, eHG,  iHG, eP-1, iP, n)])/Theta;
                                        }} else { VALUE_STORE[state_key(l_, eHQ, iHQ, eHG, iHG, eP, iP, n)] =  0.0;
                                          }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

std::string init_key = state_key(l_init, eHQ_init, iHQ_init,
                                         eHG_init, iHG_init, eP_init, iP_init, nmax);
        double pv = 0.0;
        auto it = VALUE_STORE.find(init_key);
        if (it != VALUE_STORE.end()) pv = it->second;

        Rcpp::Rcout << "Finished - " << nmax << "/" << nmaxmax
                    << " - Store size: " << VALUE_STORE.size()
                    << " - pv: " << pv << "\n";

        if (nmax > 0 && p_set.size() > 0) {
            P_val += pv / (1.0 - p_set[0]);
        }
        // mean_ += nmax * pv;
        //p_set.push_back(pv);
        double weighted_n_sum = 0.0;
        for (int l = 1; l <= K; l++) {
            std::string k = state_key(l, 1, 0, 0, 0, 0, 0, nmax);
            auto it = VALUE_STORE.find(k);
            if (it != VALUE_STORE.end()) {
                weighted_n_sum += outbreak_weights[l] * it->second;
            }
        }  // end for l
        p_set.push_back(weighted_n_sum);

    } // end nmax
    return p_set;
}   

