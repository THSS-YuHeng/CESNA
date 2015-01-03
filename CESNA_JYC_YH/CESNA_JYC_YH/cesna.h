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
    OBSERVED vector<vector<float>> X;
    // IV, par 1
    // F_uc, communities memberships, F[node_id][community_id] = float, N*C
    TOINFER vector<vector<float>> F;
    // helper vector, save the sum of each community
    vector<float> SumFV;
    // IV, par 1
    // W, logistic weight parameters arg, W[k][community_id] = float, K*(C+1)
    TOINFER vector<vector<float>> W;
    
    // P,Q 中间变量
    
    std::random_device _rd;
    std::mt19937 _mt; // random number generator
    int n_communities; // number of communities
    int n_attributes; // number of attributes, K
    
    double MinValW;
    double MaxValW;
public:
    // lambda, regularization hyperparameter, in EQ4,Page4
    GIVEN float lambda;
    
    cesna(ugraph* g): _g(g),_mt(_rd()) {}
    int estimateCommuNumber(); // estimate C number
    void setCommunityNumber(int c);
    void calculate();
    
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
		for (vector<float> FI = FU.BegI(); FI < FU.EndI(); FI++) {
			DP += FI.GetDat() * WK[FI.GetKey()];
		}
		DP += WK.Last();
		return 1.0 / ( 1.0 + exp(-DP));;
	}

	template <typename T>
	int Sign(const T& Val) {
		return Val<0?-1:(Val>0?1:0);
	}

	double inline Norm2(const TIntFltH& UV) {
		double N = 0.0;
		for (TIntFltH::TIter HI = UV.BegI(); HI < UV.EndI(); HI++) {
			N += HI.GetDat() * HI.GetDat();
		}
		return N;
	}

	double inline GetCom(const int& NID, const int& CID) {
		if (F[NID].IsKey(CID)) {
			return F[NID].GetDat(CID);
		} else {
			return 0.0;
		}
	}

	double inline GetAttr(const int& NID, const int& K) {
		if (X[NID].IsKey(K)) {
			return 1.0;
		} else {
			return 0.0;
		}
	}

	void  GradientForWK(vector<float>& GradV, const int K) {
		//GradV.Gen(NumComs + 1);
		for (int u = 0; u < F.size; u++) {
			//if (HOKIDSV[u].IsKey(K)) { continue; }
			double Pred = PredictAttrK(F[u], W[K]);					//Q(u,k) get in 
			for (TIntFltH::TIter CI = F[u].BegI(); CI < F[u].EndI(); CI++) {
				GradV[CI.GetKey()] += (GetAttr(u, K) - Pred) * GetCom(u, CI.GetKey());
			}
			GradV[n_communities] += (GetAttr(u, K) - Pred);
		}
    
		for (int c = 0; c < GradV.Len() - 1; c++) {
			GradV[c] -= LassoCoef * Sign(W[K][c]);
		}	  
	}
};

#endif /* defined(__cesna__cesna__) */
