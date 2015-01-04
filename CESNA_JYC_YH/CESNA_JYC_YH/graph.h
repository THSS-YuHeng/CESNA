//
//  graph.h
//  cesna
//
//  Created by 金宇超 on 14-12-27.
//  Copyright (c) 2014年 ccjinyc. All rights reserved.
//

#ifndef __cesna__graph__
#define __cesna__graph__

#include <stdio.h>
#include <vector>
#ifdef __GNUC__
#include <unordered_map>
using namespace std;
#else
#include <unordered_map>
using namespace std;
#endif

class ugraph;

// undirected graph
// node, edge iterator
// judge if node, edge is in graph


class ugraph {
public:
    class node {
        int                 nodeId;
        vector<int>         neighborNodeIds; // node ids
        vector<int>         attrIds;
    public:
        node(): nodeId(-1) {}
        node(int id): nodeId(id) {}
        void addNeighbor( int neighNodeId ) { neighborNodeIds.push_back(neighNodeId); }
        int getId() const   { return nodeId; }
        int getDeg() const  { return (int)neighborNodeIds.size(); }
        vector<int> getNeighbors() const { return neighborNodeIds; }
        bool isNeighbor(int vid) { return std::find(neighborNodeIds.begin(),                                                    neighborNodeIds.end(), vid) == neighborNodeIds.end(); }
    };
    typedef unordered_map<int, node>::iterator nodeI;
    
    class edge {
        int edgeId;
        int srcNodeId, desNodeId;
    public:
        edge():edgeId(-1),srcNodeId(-1),desNodeId(-1){}
        edge(int eid, int src, int des){ edgeId = eid; srcNodeId = src; desNodeId = des; }
        int getId() const   { return edgeId; }
        int getSrcId()const { return srcNodeId; }
        int getDesId()const { return desNodeId; }
    };
    typedef unordered_map<int, edge>::iterator edgeI;
private:
    unordered_map< int, node >   nodeIdMap;
    
public:
    void addNode(node n)    { nodeIdMap[n.getId()] = n; }
    node* getNode(int nid)  { return &nodeIdMap[nid]; }
    bool isNode(int nid)    { return !(nodeIdMap.find(nid) == getNodeItEnd()); }
    void delNode(int nid)   { nodeIdMap.erase(nid); }
    nodeI getNodeItBegin()  { return nodeIdMap.begin(); }
    nodeI getNodeItEnd()    { return nodeIdMap.end(); }
    
private:
    unordered_map< int, edge >   edgeIdMap;
public:
    void addEdge(edge e)    { edgeIdMap[e.getId()] = e; }
    edge* getEdge(int eid)  { return &edgeIdMap[eid]; }
    bool isEdge(int eid)    { return !(edgeIdMap.find(eid) == getEdgeItEnd()); }
    void delEdge(int eid)   { edgeIdMap.erase(eid); }
    edgeI getEdgeItBegin()  { return edgeIdMap.begin(); }
    edgeI getEdgeItEnd()    { return edgeIdMap.end(); }
    
    int graphNodeMapSize()  { return (int)nodeIdMap.size();}
    int graphEdgeMapSize()  { return (int)edgeIdMap.size();}
    
};

class graphAdapter {
    virtual ugraph load() = 0;
};

class dblp_article : graphAdapter{
    ugraph load(){
        return ugraph();
    }
};

#endif /* defined(__cesna__graph__) */
