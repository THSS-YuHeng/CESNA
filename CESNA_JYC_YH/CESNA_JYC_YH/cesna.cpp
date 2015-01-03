//
//  cesna.cpp
//  cesna
//
//  Created by 金宇超 on 14-12-27.
//  Copyright (c) 2014年 ccjinyc. All rights reserved.
//

#include "cesna.h"
#include <algorithm> // for std::shuffle
#define NegWgt (1.0)
int cesna::estimateCommuNumber() {
    return 10;
}

void cesna::calculate(double StepAlpha, double StepBeta) {
    // F_est, W_est = argmax log(P(G,X|F,W)), max likelihood
    // log(P(G,X|F,W)) = l~G + L~X
    // l~G = log(P(G|F)), l~X = log(P(X|F,W))
    // l~G = sigma([u,v in E] log(1 - exp(-FuFvT)) )- sigma([u,v not in E] FuFvT )
    // l~X = sigma( Xuk*logQuk + (1-Xuk)log(1-Quk) )
    // Fu => vector of node u
    // Quk eq2, [page 3], 1 / (1 + exp(-sigma[c] Wkc*Fuc))
    // invoke l1~regularization on W,
    // ----> F_est, W_est = argmax log(P(G,X|F,W)) - labmda|W|1
    // --------------------------------------------------------
    // block coordinate ascent approach.
    // update Fu for each node u, by fixing both W and Fv for all other node v
    // then update W while fixing community memberships F
    // 从而将非凸的优化问题转化为凸优化
    // --------------------------------------------------------
    // estimate
	n_communities = -1;
    if ( n_communities == -1 ) {
        n_communities = estimateCommuNumber();
    }
    SumFV = std::vector<float>(n_communities);
    for (int c = 0; c < n_communities; c++) {
        SumFV[c] = 0.0;
    }
    // init F, W
    for (int u = 0; u < _g->graphNodeMapSize(); u++) {
        F.push_back(vector<float>());
        for (int c = 0; c < n_communities; c++) {
            F[u].push_back(0.0);
        }
    }
    for (int k = 0; k < n_attributes; k++) {
        W.push_back(vector<float>());
        for (int c = 0; c < n_communities; c++) {
            W[c].push_back(0.0);
        }
    }
    // grad ascent
    SETUPSTAMP
    int iter = 0;
    int maxiter = 5; // to be set
    vector<int> shuffleU; // for shuffle
    for (int i = 0; i < _g->graphNodeMapSize(); i++) {
        shuffleU.push_back(i);
    }
    while (iter < maxiter) {
        std::shuffle(shuffleU.begin(), shuffleU.end(), _mt);
        // for every node u
        // Updating community memberships. P4,EQ(5),EQ(6)
        for (int i = 0; i < shuffleU.size(); i++) {
            int uid = shuffleU[i];
            // part1 dlg/dfu =
            //          sigma[v in Nei(u)]( Fvc * ( exp(-FuFv)) / ( 1 - exp(-FuFv)) )
            //         -sigma[v not in Nei(u)]( FuFv )
            // grad for node u
            ugraph::node* n = _g->getNode(uid);
            int deg = n->getDeg();
            vector<float> gradV(n_communities);
            vector<float> gradU(n_communities);
            vector<float> predV(deg); // exp(-FuFv)
            float val = 0.0;
            // 求gradv
            for (int ni = 0; ni < deg; ni++) {
                int nid = n->getNeighbors()[ni];
                // prediction 函数存疑, 1/1-pnocom是什么意思, pnocom = 1/N
                // 为什么要做这个修正，论文中没体现出来
                // 会不会是什么特别的方法
                // predV -> Puv? P3右半边上面的定义
                predV[ni] = exp(-dot(F[uid], F[nid]));//TODO Prediction(uid, nid);
            }
            for (int c = 0; c < n_communities; c++) {
                double val = 0.0;
                for (int ni = 0; ni < deg; ni++) {
                    int nid = n->getNeighbors()[ni];
                    // NegWgt = 1.0?
                    val += predV[ni] * F[nid][c] / (1.0 - predV[ni]) + NegWgt * F[nid][c];
                }
                // 计算v not in Nei的部分
                val -= NegWgt * (SumFV[c] - F[uid][c]);
                gradV[c] = val;
                SumFV[c] = 0.0;
                for (int sumi = 0; sumi < _g->graphNodeMapSize(); sumi++) {
                    SumFV[c] += F[sumi][c];
                }
            }
            // L1 L2 regularization
            // add attribute part
            // part2 dlx/dfu = sigma[k]( Xuk - Quk ) * Wkc
            vector<float> attrV(n_attributes);
            for (int k = 0; k < n_attributes; k++) {
                attrV[k] = 0.0;//TODO PredictAttrK(F[UID], W[k]);
            }
            for (int c = 0; c < gradV.size(); c++) {
                for (int k = 0; k < n_attributes; k++) {
                    gradV[c] += 1.0 * (X[uid][k] - attrV[k]) * W[c][k];
                }
            }
            for (int c = 0; c < gradV.size(); c++) {
                if (F[uid][c] == 0.0 && gradV[c] < 0.0) { continue; }
                if (fabs(gradV[c]) < 0.0001) { continue; }
                gradU[c] = gradV[c];
            }
            for (int c = 0; c < gradU.size(); c++) {
                if (gradU[c] >= 10) { gradU[c] = 10; }
                if (gradU[c] <= -10) { gradU[c] = -10; }
            }
            // norm2 GradV < 1e-4
            // step size
            float LearnRate = 0.0;//TODO GetStepSizeByLineSearch
            // 更新Fuc
            // P4,Eq(6), Fnewuc = max ( 0, Folduc + alpha(dlg/dfu+ dlx/dfu))
            for (int ci = 0; ci < gradV.size(); ci++) {
                double Change = LearnRate * gradV[ci];
                double NewFuc = F[uid][ci] + Change;
                if (NewFuc <= 0.0) {
                    F[uid][ci] = 0.0;
                } else {
                    F[uid][ci] = NewFuc;
                }
            }
        }
        
        for (int k = 0; k < n_attributes; k++) {
            vector<float> gradWV(n_communities);
            // calc W
            // grad W
            // TODO GradientForWK
			GradientForWK(gradWV, k);
            // TODO Norm2 gradWV
			if (Norm2(gradWV) < 1e-4) { continue; }
            // step size
            double learnRate = 0.0; // TODO GetStepSizeByLineSearchForWK
			learnRate = GetStepSizeByLineSearchForWK(k, gradWV, gradWV, StepAlpha, StepBeta);
            // 更新 Wkc
            if (learnRate == 0.0) { continue; }
            for (int c = 0; c < gradWV.size(); c++){
                W[k][c] += learnRate * gradWV[c];
                if (W[k][c] < MinValW) { W[k][c] = MinValW; }
                if (W[k][c] > MaxValW) { W[k][c] = MaxValW; }
            }
        }
        std::cout << "iter: " << iter << " " << GETSTAMP << std::endl;
        iter++;
    }
    
    // get com
    // dump com
}