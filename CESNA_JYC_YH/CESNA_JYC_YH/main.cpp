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
#include <fstream>

#include <unordered_map>

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
    std::unordered_map<int,int> ind_id_map;
    std::unordered_map<int, vector<int>> citemap;
    std::vector<int> idv;
    char buffer[2000];
    ifstream ifs("/Users/jin-yc10/Development/data_homework/Hw4/618506165_6_Project4/cora/cora.cites.txt");
    if (! ifs.is_open())
    { cout << "Error opening file"; exit (1); }
    int citen = 0;
    int id, cite;
    while (ifs >> id >> cite)
    {
        if( id == 0 && cite == 0){std::cout << "all zero" << endl; break;}
        std::unordered_map<int,int>::iterator fi = ind_id_map.find(id);
        if(fi == ind_id_map.end()) {
            //std::cout << id << std::endl;
            ind_id_map[id] = (int)idv.size();
            idv.push_back(id);
            citemap[id] = vector<int>();
            ugraph::node n(id);
            ug.addNode(n);
        }
        std::unordered_map<int,int>::iterator ci = ind_id_map.find(cite);
        if(ci == ind_id_map.end()) {
            //std::cout << cite << std::endl;
            ind_id_map[cite] = (int)idv.size();
            idv.push_back(cite);
            ugraph::node n(cite);
            ug.addNode(n);
        }
//        std::cout << citen << " " << id << " " << cite << endl;
        citemap[id].push_back(cite);
        citen++;
    }
    int eid = 0;
    for (std::unordered_map<int, vector<int>>::iterator i = citemap.begin();
         i != citemap.end(); i++) {
        vector<int> cites = i->second;
//        std::cout << ind_id_map[i->first] << " :";
        for (int x: cites) {
//            std::cout << ind_id_map[x] << " ";
            ug.getNode(i->first)->addNeighbor(x);
            ugraph::edge e1(eid, i->first, x);
            eid++;
            ug.addEdge(e1);
            ug.getNode(x)->addNeighbor(i->first);
            ugraph::edge e2(eid, x, i->first);
            eid++;
            ug.addEdge(e2);
        }
//        std::cout << std::endl;
    }
//    std::cout << idv.size() << std::endl;
//    std::cout << ug.graphNodeMapSize() << std::endl;
    for (ugraph::nodeI ni = ug.getNodeItBegin(); ni != ug.getNodeItEnd(); ni++) {
        //std::cout << ni->first << ": ";
        std::unordered_map<int,int>::iterator fi = ind_id_map.find(ni->first);
        ugraph::node n = ni->second;
//        for (int nei: n.getNeighbors()) {
//            std::cout << nei << " ";
//        }
//        std::cout << std::endl;
    }
//    std::cout << ug.graphEdgeMapSize() << std::endl;
//    std::cout << eid << " " << citen << std::endl;

    unordered_map<int,vector<int>> X;
#ifdef __APPLE__
    ifstream ifattrs("/Users/jin-yc10/Development/data_homework/Hw4/618506165_6_Project4/cora/cora.content.txt");
#else
    ifstream ifattrs("/Users/jin-yc10/Development/data_homework/Hw4/618506165_6_Project4/cora/cora.content.txt");
#endif
    if (! ifattrs.is_open())
    { cout << "Error opening file"; exit (1); }
    int attrid;
    while ( ifattrs >> attrid ) {
        X[attrid] = vector<int>();
        for (int i = 0; i < 1433; i++) {
            int at;
            ifattrs >> at;
            X[attrid].push_back(at);
        }
        char cat[30];
        ifattrs >> cat;
//        std::cout << attrid << " " << cat << " " << X[attrid].size() << endl;
    }
    
	double StepAlpha = 0.1, StepBeta = 0.1;
    cesna c(&ug, X);
    c.calculate(StepAlpha, StepBeta);
    
    return 0;
}
