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
    ug.addNode(ugraph::node());
    std::cout << ug.isNode(0) << std::endl;
    ugraph::node* un = ug.getNode(0); un->addNeighbor(1);
    std::cout << un->getDeg() << std::endl;
    ugraph::node* un2 = ug.getNode(0); un2->addNeighbor(2);
    un->addNeighbor(3);
    std::cout << un->getDeg() << std::endl;
    return 0;
}
