//
//  cesna.h
//  cesna
//
//  Created by 金宇超 on 14-12-27.
//  Copyright (c) 2014年 ccjinyc. All rights reserved.
//

#ifndef __cesna__cesna__
#define __cesna__cesna__

#include <stdio.h>
#include "graph.h"
#include "debug.h"
#include <random>
#include <vector>
#include <set>
#ifdef __GNUC__
#include <unordered_map>
using namespace std;
#else
#include <unordered_map>
using namespace std;
#endif

#define TOINFER
#define OBSERVED
#define GIVEN

class cesna {
    // Graph, as A or G, include N and E
    OBSERVED ugraph* _g;
    // X_uk, node attributes, X[node_id][attr_id] = 0 or 1
    // indicate if the node has the attr
    OBSERVED unordered_map<int,vector<int>> X;
    // IV, par 1
    // F_uc, communities memberships, F[node_id][community_id] = float, N*C
    TOINFER unordered_map<int,vector<float>> F;
    // helper vector, save the sum of each community
    vector<float> SumFV;
    // IV, par 1
    // W, logistic weight parameters arg, W[k][community_id] = float, K*(C+1)
    TOINFER vector<vector<float>> W;
    
    // P,Q 中间变量
    // Node ids
    vector<int> nids;
    
    std::random_device _rd;
    std::mt19937 _mt; // random number generator
    int n_communities; // number of communities
    int n_attributes; // number of attributes, K
    float LassoCoef; // to be set
    double MinValW;
    double MaxValW;
    double MinVal;
    double MaxVal;
public:
    // lambda, regularization hyperparameter, in EQ4,Page4
    GIVEN float lambda;
    
    cesna(ugraph* g, unordered_map<int,vector<int>> X): _g(g), X(X),_mt(_rd()) {
        for (ugraph::nodeI ni = g->getNodeItBegin(); ni != g->getNodeItEnd(); ni++) {
            nids.push_back(ni->first);
            F[ni->first] = vector<float>();
            //D({ std::cout << ni->first << std::endl; })
            n_attributes = (int)X[ni->first].size();
        }
        D( {
         std::cout << "GRAPH HAVE " << nids.size() << " NODES" << std::endl;
        } )
    }
    
    int estimateCommuNumber(); // estimate C number
    void setCommunityNumber(int c);
    void calculate(double StepAlpha, double StepBeta);
    
    double dot(vector<float> v1, vector<float> v2) {
        double r = 0.0;
        if( v1.size() == v2.size() ) {
            for (int i = 0; i < v1.size(); i++) {
                r += v1[i]*v2[i];
            }
        }
        return r;
    }

	double inline PredictAttrK(vector<float>& FU, vector<float>& WK) {
		double DP = 0.0;
        for (int c = 0; c < n_communities; c++) {
			DP += FU[c] * WK[c];
		}
		DP += WK[WK.size()-1];
		return 1.0 / ( 1.0 + exp(-DP));;
	}

	template <typename T>
	int Sign(const T& Val) {
		return Val<0?-1:(Val>0?1:0);
	}

	double inline Norm2(const vector<float>& UV) {
		double N = 0.0;
        for (int i = 0; i < UV.size(); i++) {
			N += UV[i]*UV[i];
		}
		return N;
	}

	void  GradientForWK(vector<float>& GradV, const int K) {
		//GradV.Gen(NumComs + 1);
		GradV.push_back(1);

		for (int u = 0; u < F.size(); u++) {
			//if (HOKIDSV[u].IsKey(K)) { continue; }
			double Pred = PredictAttrK(F[nids[u]], W[K]);					//Q(u,k) get in
            for (int c = 0; c < n_communities; c++) {
				GradV[c] += (X[nids[u]][K] - Pred) * F[nids[u]][c];
			}
			GradV[n_communities] += (X[nids[u]][K] - Pred);
		}
    
		for (int c = 0; c < GradV.size() - 1; c++) {
			GradV[c] -= LassoCoef * Sign(W[K][c]);
		}	  
	}

	double inline GetAttr(const int& NID, const int& K) {
        return X[NID][K];
	}

    double LikelihoodForRow(const int UID, vector<float>& FU) {
        double L = 0.0;
        
        ugraph::node* NI = _g->getNode(UID);
        for (int e = 0; e < NI->getDeg(); e++) {
            int v = NI->getNeighbors()[e];
            if (v == UID) { continue; }
            L += log (1.0 - exp(-dot(FU, F[v])) + 1.0 * dot(FU, F[v]) );
        }
        for (int HI = 0; HI < FU.size(); HI++) {
            L -= 1.0 * (SumFV[HI] - F[UID][HI] * FU[HI]);
        }
        //add regularization
//        if (RegCoef > 0.0) { //L1
//            L -= RegCoef * Sum(FU);
//        }
//        if (RegCoef < 0.0) { //L2
//            L += RegCoef * Norm2(FU);
//        }
        L *= (1.0 - 0.5);
        // add attribute part
        for (int k = 0; k < n_attributes; k++) {
            L += 0.5 * LikelihoodAttrKForRow(UID, k, FU, W[k]);
        }
        return L;
    }
    
	double LikelihoodAttrKForRow(const int UID, int K, vector<float>& FU, vector<float>& WK) {
		double Prob = PredictAttrK(FU, WK);
		double L = 0.0;
		if (GetAttr(UID, K)) { 
			L = Prob == 0.0? -100.0: log(Prob);
		} else {
			L = Prob == 1.0? -100.0: log(1.0 - Prob);
		}
		return L;
	}

	double LikelihoodForWK(int K, vector<float>& WK) {
		double L = 0.0;
		for (int u = 0; u < F.size(); u++) {
			//if (HOKIDSV[u].IsKey(K)) { continue; }
			L += LikelihoodAttrKForRow(nids[u], K, F[nids[u]], WK);
		}
		for (int c = 0; c < WK.size() - 1; c++) {
			L -= LassoCoef * fabs(WK[c]);
		} 
		return L;
	}

    double GetStepSizeByLineSearch(const int UID, vector<float>& DeltaV, vector<float>& GradV, const double& Alpha, const double& Beta, const int MaxIter = 10) {
        double StepSize = 1.0;
        double InitLikelihood = LikelihoodForRow(UID, F[UID]);
        vector<float> NewVarV(DeltaV.size());
        for(int iter = 0; iter < MaxIter; iter++) {
            for (int i = 0; i < DeltaV.size(); i++){
                int CID = i;
                double NewVal = F[UID][CID] + StepSize * DeltaV[CID];
                if (NewVal < MinVal) { NewVal = MinVal; }
                if (NewVal > MaxVal) { NewVal = MaxVal; }
                NewVarV[i] = NewVal;
            }
            if (LikelihoodForRow(UID, NewVarV) < InitLikelihood + Alpha * StepSize * dot(GradV, DeltaV)) {
                StepSize *= Beta;
            } else {
                break;
            }
            if (iter == MaxIter - 1) { 
                StepSize = 0.0;
                break;
            }
        }
        return StepSize;
    }
    
	double GetStepSizeByLineSearchForWK(int K, vector<float>& DeltaV, vector<float>& GradV, double& Alpha, double& Beta, int MaxIter = 10) {
		double StepSize = 1.0;
		double InitLikelihood = LikelihoodForWK(K, W[K]);
		vector<float> NewVarV(DeltaV.size());
		//IAssert(DeltaV.size() == n_communities + 1);
		for(int iter = 0; iter < MaxIter; iter++) {
			for (int c = 0; c < DeltaV.size(); c++) {
				double NewVal = W[K][c] + StepSize * DeltaV[c];
				if (NewVal < MinValW) {
					NewVal = MinValW;
				}
				if (NewVal > MaxValW) {
					NewVal = MaxValW;
				}
				NewVarV[c] = NewVal;
			}
			if (LikelihoodForWK(K, NewVarV) < InitLikelihood + Alpha * StepSize * dot(GradV, DeltaV)) {
				StepSize *= Beta;
			} else {
				break;
			}
			if (iter == MaxIter - 1) { 
				StepSize = 0.0;
				break;
			}
		}
		return StepSize;
	}
};

#endif /* defined(__cesna__cesna__) */
