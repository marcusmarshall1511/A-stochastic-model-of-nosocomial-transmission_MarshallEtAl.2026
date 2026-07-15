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
NumericVector fast_trial_cpp_HQ(int nmaxmax, double gamma, double rC,
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

    //  Initial state (mirrors R x_init) 
    int l_init   = Hyp.size() - 1; // l_init = 4
    int eHQ_init = 0, iHQ_init = 0;
    int eHG_init = 0, iHG_init = 0;
    int eP_init  = 0, iP_init  = 0;

    int q_init = eHQ_init + iHQ_init;  // total HQ infected at start = 0
    int g_init = eHG_init + iHG_init;  // total HG infected at start = 0

    //  Patient state index helpers  
    // Maps (eP, iP) -> in order: p=0,1,...,NP; eP=0..p
    // for each p, there are p+1 states (eP=0..p, iP=p-eP)
    int dimB = 0;
    for (int p = 0; p <= NP; p++) dimB += (p + 1); //dimB gives size of sub-sub-level and is used to build matrices/vectors in loop

    // Patient tuple [eP][iP] -> 0-based index
    std::vector<std::vector<int>> pat_idx(NP + 1, std::vector<int>(NP + 1, -1));
    std::vector<std::pair<int,int>> pat_rev(dimB); // 0-based index -> (eP, iP)
    {
        int idx = 0;
        for (int p = 0; p <= NP; p++) {
            for (int eP = 0; eP <= p; eP++) {
                int iP = p - eP;
                pat_idx[eP][iP] = idx;
                pat_rev[idx] = {eP, iP};
                idx++;
            }
        }
    }

    //  VALUE_STORE maps full state key -> alpha value 
    std::unordered_map<std::string, double> VALUE_STORE;
    VALUE_STORE.reserve(1 << 20); // reserve 2^20 entries

    auto get_val = [&](int l, int ehq, int ihq, int ehg, int ihg,
                       int ep, int ip, int nn) -> double {
        auto it = VALUE_STORE.find(state_key(l, ehq, ihq, ehg, ihg, ep, ip, nn));
        return (it != VALUE_STORE.end()) ? it->second : 0.0;
    };

    //  Outputs
    NumericVector p_set;
    double P_val  = 0.0;
    // double mean_  = 0.0; 

    // =========================================================
    // MAIN LOOP
    //   for nmax in 0:nmaxmax
    //     for n in 0:nmax
    //       for l_ in 0:l_init
    //         q  = q_init + nmax - n   [single value for each nmax,n]

    //    CONSTRUCT MATRIX AND BOUNDARY VECTOR

    //           for iHQ in iHQ_init:q        [i.e. iHQ from 0 up to q]
    //             for g in NHG:g_init        [NHG to g_init]
    //               for iHG in g:iHG_init    [g to iHG_init]
    //                 -> solve subsubmatrix system
    // =========================================================
    for (int nmax = 0; nmax <= nmaxmax; nmax++) {

        VALUE_STORE.clear();  // Each nmax, clear VALUE_STORE

        for (int n = 0; n <= nmax; n++) { // increasing n
            for (int l_ = 0; l_ <= l_init; l_++) { // decreasing l
                int q = q_init + nmax - n;
                if (q > NHQ) continue; // ensures exposed+infectious dont exceed NHQ

                // iHQ ranges from iHQ_init (0) up to q
                for (int iHQ = q; iHQ >= iHQ_init; iHQ--) { //given q, start with no exposed HQ (all infectious)
                    int eHQ = q - iHQ;

                    // g ranges from NHG to g_init
                    for (int g = NHG; g >= g_init; g--) { 

                        // iHG ranges from g to iHG_init
                        for (int iHG = g; iHG >= iHG_init; iHG--) { // given g, start with no exposed HG (all infectious)
                            int eHG = g - iHG;

                            // Build and solve (I - B) alpha = b
                            // over patient states (eP, iP) with eP+iP <= NP



                            // Initialise matrix and boundary vector
                            std::vector<double> B(dimB * dimB, 0.0);
                            std::vector<double> b_vec(dimB, 0.0);


                            // Calculate rates
                            for (int idx = 0; idx < dimB; idx++) {
                                auto [eP, iP] = pat_rev[idx];
                                int p = eP + iP;

                                double delta_iso = Hyp[l_];

                                double mu_IP = iP * deltaP;
                                double mu_EP = eP * deltaP;

                                double SHQ = NHQ - eHQ - iHQ;
                                double SHG = NHG - eHG - iHG;
                                double SP  = NP  - eP  - iP;

                                double lambda_HQ = SHQ * (Beta_O[l_] + iHQ * Beta_HQ
                                                         + iHG * Beta_H + iP * Beta_HQP);
                                double lambda_HG = SHG * (iHG * Beta_HG + iHQ * Beta_H
                                                         + iP * Beta_HGP);
                                double lambda_P  = SP  * (iHQ * Beta_HQP + iHG * Beta_HGP + iP*Beta_P);  // Beta_P = 0 so patient-to-patient infection not considered

                                double eta_HQ = eHQ * zeta_HQ;
                                double eta_HG = eHG * zeta_HG;
                                double eta_P  = eP  * zeta_P;

                                double gamma_rate = iHQ * gamma_HQ + iHG * gamma_HG
                                                  + iP  * gamma_P;

                                double Theta = mu_IP + mu_EP + delta_iso
                                             + lambda_HQ + lambda_HG + lambda_P
                                             + eta_HQ + eta_HG + eta_P + gamma_rate;

                    
                                // Boundary / absorbing state
                                bool is_zero_state = (l_ + q + g + p == 0); //if zero state, this is TRUE

                                if (is_zero_state) {
                                    b_vec[idx] = (n == 0) ? 1.0 : 0.0; //TRUE returns 1, FALSE returns 0
                                } else if (Theta <= 0.0) {
                                    b_vec[idx] = 0.0; // if no transistions are possible (absorbing state)
                                } else {
                                    double rhs = 0.0; //initialise entry for boundary vector

                                    if (n == 0) {
                                        rhs += gamma_rate; //detection
                                    } else {
                                        if (SHQ > 0) {
                                            rhs += lambda_HQ * get_val(l_, eHQ + 1, iHQ,
                                                                        eHG, iHG, eP, iP, n - 1); //HQ exposure
                                        }  
                                    }
                                    if (SHG > 0) {
                                            rhs += lambda_HG * get_val(l_, eHQ, iHQ,
                                                                        eHG + 1, iHG, eP, iP, n); //HG exposure
                                        }
                                    if (l_ > 0) {
                                        rhs += delta_iso * get_val(l_ - 1, eHQ, iHQ,
                                                                   eHG, iHG, eP, iP, n); // progression to removal of isolated patient
                                    }
                                    if (eHQ > 0) {
                                        rhs += eta_HQ * get_val(l_, eHQ - 1, iHQ + 1,
                                                                eHG, iHG, eP, iP, n); // maturation of exposed HQ to infectious
                                    }
                                    if (eHG > 0) {
                                        rhs += eta_HG * get_val(l_, eHQ, iHQ,
                                                                eHG - 1, iHG + 1, eP, iP, n); // maturation of exposed HG to infectious
                                    }

                                    b_vec[idx] = rhs / Theta;

                                    //  Build matrix for in-subsublevel transitions (patients)

                                    if (eP > 0) {
                                        int idx_y_zeta = pat_idx[eP - 1][iP + 1];
                                        B[idx * dimB + idx_y_zeta] = eta_P / Theta; //// maturation of exposed P to infectious

                                        int idx_y_muE = pat_idx[eP - 1][iP];
                                        B[idx * dimB + idx_y_muE] += mu_EP / Theta; // removal of exposed P
                                    }
                                    if (iP > 0) {
                                        int idx_y = pat_idx[eP][iP - 1];
                                        B[idx * dimB + idx_y] += mu_IP / Theta; // removal of infectious P
                                    }
                                    if (SP > 0 && eP + 1 + iP <= NP) {
                                        int idx_y = pat_idx[eP + 1][iP];
                                        B[idx * dimB + idx_y] += lambda_P / Theta; // P exposure
                                    }
                                }
                            } 
                            // Matrix and Boundary vector built
                            //  Solve (I - B) alpha = b for this subsublevel
                            
                            arma::mat Bmat(B.data(), dimB, dimB, false, true); 
                            Bmat = Bmat.t();   // .t() = transpose; undoes the row/col mismatch
                            arma::vec bvec(b_vec.data(), dimB, false, true);
                            arma::mat I = arma::eye(dimB, dimB);

                            arma::vec alpha = arma::solve(I - Bmat, bvec);



                            //  Store solved alpha values into VALUE_STORE 
                            for (int idx = 0; idx < dimB; idx++) {
                                int eP = pat_rev[idx].first;
                                int iP = pat_rev[idx].second;
                                VALUE_STORE[state_key(l_, eHQ, iHQ, eHG, iHG, eP, iP, n)] = alpha(idx);
                            }

                        } 
                    } 
                } 
            }
        } 

        //  Extract probability for this nmax 
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
        p_set.push_back(pv);

        // Optional early exit (mirrors the commented-out R break)
        // if (P_val > 0.9999) break;

    } // end nmax
    
    double p0 = p_set[0];
    
    for (int i = 1; i < p_set.size(); i++) {
      p_set[i] /= (1.0 - p0);
    }
    
    p_set[0] = 0.0;
    
    
    return p_set;
}
