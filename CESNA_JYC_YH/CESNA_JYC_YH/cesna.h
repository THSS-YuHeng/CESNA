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
};

#endif /* defined(__cesna__cesna__) */
