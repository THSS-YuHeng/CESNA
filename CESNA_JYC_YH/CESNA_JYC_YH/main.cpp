//
//  main.cpp
//  cesna
//
//  Created by 金宇超 on 14-12-27.
//  Copyright (c) 2014年 ccjinyc. All rights reserved.
//

#include <iostream>
#include "debug.h"
#include "graph.h"
#include "cesna.h"

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    SETUPSTAMP
    STAMP("test stamp", {});
    PROFILE("test profile",
        {
            for(int i = 0; i < 1000000; i++) {
            }
        }
    )
    STAMP("test stamp", {});

    ugraph ug;
    for (int i = 0; i < 5; i++) {
        ugraph::node n(i);
        ug.addNode(n);
    }
    
	double StepAlpha = 0.1, StepBeta = 0.1;
    cesna c(&ug);
    c.calculate(StepAlpha, StepBeta);
    
    return 0;
}
