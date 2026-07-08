#include <Rcpp.h>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// Helper to create the string key fast
inline std::string state_key(int l, int eHQ, int iHQ, int eHG, int iHG, int eP, int iP, int n) {
    return std::to_string(l) + "_" + std::to_string(eHQ) + "_" + 
           std::to_string(iHQ) + "_" + std::to_string(eHG) + "_" + 
           std::to_string(iHG) + "_" + std::to_string(eP) + "_" + 
           std::to_string(iP) + "_" + std::to_string(n);
}

// [[Rcpp::export]]
NumericVector fast_trial_cpp(int nmaxmax, double gamma, double rC, String cohort, double beta, double zeta, double eps) {
    
    // --- Setup Variables ---
    int NHQ = 0, NHG = 0, NP = 12;
    if (cohort == "S") { NHQ = 2; NHG = 11; }
    else if (cohort == "M") { NHQ = 6; NHG = 7; }
    else if (cohort == "N") { NHQ = 13; NHG = 0; }
    
    // Rates and Parameters
    double deltaP = 1.0 / 6.0;
    std::vector<double> Hyp = {0.0, 6.69470955813761, 0.555473683554152, 0.55547824700609, 0.555478430962926};
    double rI = 0.95;
    
    std::vector<double> Beta_O(5, 0.0);
    if (NHQ > 0) {
        for(int i = 1; i <= 4; i++) {
            Beta_O[i] = (1.0 - rI) * beta * 2.0 / NHQ;
        }
    }
    
    double Beta_H   = eps * (2.0 / 3.0) * beta;
    double Beta_HQ  = (2.0 / 3.0) * beta;
    double Beta_HG  = (2.0 / 3.0) * beta;
    double Beta_HQP = (1.0 - rC) * beta;
    double Beta_HGP = beta;
    double Beta_P   = 0.0;
    
    double zeta_HQ = zeta, zeta_HG = zeta, zeta_P = zeta;
    double gamma_HQ = gamma, gamma_HG = gamma, gamma_P = gamma;

    // Initial states
    int l_init = 4;
    int eHQ_init = 0, iHQ_init = 0;
    int eHG_init = 0, iHG_init = 0;
    int eP_init = 0, iP_init = 0;
    
    int h_init = iHQ_init + iHG_init;
    int q_init = iHQ_init;
    int p_init = iP_init;
    
    std::unordered_map<std::string, double> VALUE_STORE;
    
    // Quick helper to fetch from map safely and default to 0 if missing
    auto get_val = [&](int l, int ehq, int ihq, int ehg, int ihg, int ep, int ip, int nn) {
        std::string k = state_key(l, ehq, ihq, ehg, ihg, ep, ip, nn);
        auto it = VALUE_STORE.find(k);
        return (it != VALUE_STORE.end()) ? it->second : 0.0;
    };
    
    double mean_ = 0;
    double P_val = 0;
    NumericVector p_set;
    
    
    for (int nmax = 0; nmax <= nmaxmax; nmax++) {
        for (int n = 0; n <= nmax; n++) {
            for (int l_ = 0; l_ <= l_init; l_++) {
                
                int hmax = std::min(h_init + nmax - n, NHQ + NHG);
                int hmin = h_init;
                
                for (int h = hmax; h >= hmin; h--) {
                    int qlow = std::max(0, h - NHG);
                    int qupp = std::min(h, NHQ);
                    
                    for (int q = std::min(q_init + h - h_init, qupp); q >= std::max(q_init, qlow); q--) {
                        for (int iHQ = q; iHQ >= iHQ_init; iHQ--) { 
                            for (int iHG = (h - q); iHG >= iHG_init; iHG--) {
                                
                                int ePmax = std::min(eP_init + nmax - n - (h - h_init), NP);
                                
                                for (int eP = 0; eP <= ePmax; eP++) {
                                    int iPmax = std::min(iP_init + ePmax - eP, NP - eP);
                                    
                                    for (int iP = 0; iP <= iPmax; iP++) {
                                        
                                        int eHQ = q - iHQ;
                                        int eHG = h - q - iHG;
                                        
                                        std::string key_x = state_key(l_, eHQ, iHQ, eHG, iHG, eP, iP, n);
                                        
                                        // If it exists, skip
                                        if (VALUE_STORE.count(key_x) > 0) continue;
                                        
                                        // Rates
                                        double mu_IP = iP * deltaP;
                                        double mu_EP = eP * deltaP;
                                        double delta_iso = Hyp[l_]; 
                                        
                                        double lambda_HQ = (NHQ - eHQ - iHQ) * (Beta_O[l_] + iHQ*Beta_HQ + iHG*Beta_H + iP*Beta_HQP);
                                        double lambda_HG = (NHG - eHG - iHG) * (iHG*Beta_HG + iHQ*Beta_H + iP*Beta_HGP);
                                        double lambda_P  = (NP - eP - iP) * (iHQ*Beta_HQP + iHG*Beta_HGP + iP*Beta_P);
                                        
                                        double eta_HQ = eHQ * zeta_HQ;
                                        double eta_HG = eHG * zeta_HG;
                                        double eta_P  = eP * zeta_P;
                                        
                                        double gamma_rate = iHQ*gamma_HQ + iHG*gamma_HG + iP*gamma_P;
                                        
                                        double Theta = mu_IP + mu_EP + delta_iso + lambda_HQ + lambda_HG + lambda_P + eta_HQ + eta_HG + eta_P + gamma_rate;
                                        
                                        double value = 0.0;
                                        
                                        
                                        if (n == 0) {
                                            if (l_ + iHQ + iHG + iP + eHQ + eHG + eP == 0) {
                                                value = 1.0;
                                            } else {
                                                double term_gamma  = gamma_rate;
                                                double term_deltaL = (l_ > 0)  ? get_val(l_-1, eHQ, iHQ, eHG, iHG, eP, iP, n) * delta_iso : 0.0;
                                                double term_muE    = (eP > 0)  ? get_val(l_, eHQ, iHQ, eHG, iHG, eP-1, iP, n) * mu_EP : 0.0;
                                                double term_muI    = (iP > 0)  ? get_val(l_, eHQ, iHQ, eHG, iHG, eP, iP-1, n) * mu_IP : 0.0;
                                                double term_etaHQ  = (eHQ > 0) ? get_val(l_, eHQ-1, iHQ+1, eHG, iHG, eP, iP, n) * eta_HQ : 0.0;
                                                double term_etaHG  = (eHG > 0) ? get_val(l_, eHQ, iHQ, eHG-1, iHG+1, eP, iP, n) * eta_HG : 0.0;
                                                double term_etaP   = (eP > 0)  ? get_val(l_, eHQ, iHQ, eHG, iHG, eP-1, iP+1, n) * eta_P : 0.0;
                                                
                                                if (Theta > 0.0) {
                                                    double value = (1.0 / Theta) * (term_gamma + term_deltaL + term_muE + term_muI + term_etaHQ + term_etaHG + term_etaP);
                                                    VALUE_STORE[key_x] = value;
                                                }
                                            }
                                        } 
                                        else if (n > 0) {
                                            if (l_ + h + iP + eP == 0) {
                                                value = 0.0;
                                            } else {
                                                double term_deltaL = (l_ > 0)  ? get_val(l_-1, eHQ, iHQ, eHG, iHG, eP, iP, n) * delta_iso : 0.0;
                                                double term_muE    = (eP > 0)  ? get_val(l_, eHQ, iHQ, eHG, iHG, eP-1, iP, n) * mu_EP : 0.0;
                                                double term_muI    = (iP > 0)  ? get_val(l_, eHQ, iHQ, eHG, iHG, eP, iP-1, n) * mu_IP : 0.0;
                                                double term_etaHQ  = (eHQ > 0) ? get_val(l_, eHQ-1, iHQ+1, eHG, iHG, eP, iP, n) * eta_HQ : 0.0;
                                                double term_etaHG  = (eHG > 0) ? get_val(l_, eHQ, iHQ, eHG-1, iHG+1, eP, iP, n) * eta_HG : 0.0;
                                                double term_etaP   = (eP > 0)  ? get_val(l_, eHQ, iHQ, eHG, iHG, eP-1, iP+1, n) * eta_P : 0.0;
                                                
                                                double term_lambdaHQ = (NHQ - eHQ - iHQ > 0) ? get_val(l_, eHQ+1, iHQ, eHG, iHG, eP, iP, n-1) * lambda_HQ : 0.0;
                                                double term_lambdaHG = (NHG - eHG - iHG > 0) ? get_val(l_, eHQ, iHQ, eHG+1, iHG, eP, iP, n-1) * lambda_HG : 0.0;
                                                double term_lambdaP  = (NP  - eP  - iP  > 0) ? get_val(l_, eHQ, iHQ, eHG, iHG, eP+1, iP, n-1) * lambda_P : 0.0;
                                                
                                                if (Theta > 0.0) {
                                                    double value = (1.0 / Theta) * (term_deltaL + term_muE + term_muI + term_etaHQ + term_etaHG + term_etaP + term_lambdaHQ + term_lambdaHG + term_lambdaP);
                                                    VALUE_STORE[key_x] = value;
                                                }
                                            }
                                        }
                                        
                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // Final P extraction per nmax
        std::string init_key = state_key(l_init, eHQ_init, iHQ_init, eHG_init, iHG_init, eP_init, iP_init, nmax);
        double p = VALUE_STORE.count(init_key) ? VALUE_STORE[init_key] : 0.0;
        
        if (nmax > 0 && p_set.size() > 0) {
            P_val += p / (1.0 - p_set[0]);
        }
        mean_ += nmax * p;
        p_set.push_back(p);
        
        // if (P_val > 0.9999) break;
    }
    
    return p_set;
}
