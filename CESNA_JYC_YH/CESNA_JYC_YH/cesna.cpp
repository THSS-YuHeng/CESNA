//
//  cesna.cpp
//  cesna
//
//  Created by 金宇超 on 14-12-27.
//  Copyright (c) 2014年 ccjinyc. All rights reserved.
//

#include "cesna.h"

void cesna::calculate() {
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
}