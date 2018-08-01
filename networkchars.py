# -*- coding: utf-8 -*-
"""
author:vishal
"""
import snap
import numpy as np
import matplotlib.pyplot as plt
# Setup
erdosRenyi = None
smallWorld = None
collabNet = None

# Problem 1.1
def genErdosRenyi(N=5242, E=14484):
    Graph = snap.TUNGraph.New()
    # Adding the nodes
    for i in range(N):
        Graph.AddNode(i)
    # Adding the edges
    while Graph.GetEdges() < E:
        nodes = np.random.choice(N,2,replace = False)
        if not Graph.IsEdge(nodes[0],nodes[1]):
            Graph.AddEdge(nodes[0],nodes[1])
    return Graph

def genCircle(N=5242):
    Graph = snap.TUNGraph.New()
# Add nodes
    for i in range(N):
        Graph.AddNode(i)
# Add edges
    for i in range(N-1):
        Graph.AddEdge(i,i+1) 
    Graph.AddEdge(N-1,0)
    return Graph

def connectNbrOfNbr(Graph, N=5242):
# Add edges
    for i in range(N-2):
        Graph.AddEdge(i,i+2)

    Graph.AddEdge(N-3,N-1)
    Graph.AddEdge(N-2,0)
    Graph.AddEdge(N-1,1)
    return Graph

def connectRandomNodes(Graph, M=4000):
# Add edges
    target_numedges = Graph.GetEdges() + M
    num_nodes = Graph.GetNodes()
    while Graph.GetEdges() < target_numedges:
        nodes = np.random.choice(num_nodes,2,replace = False)
        if not Graph.IsEdge(nodes[0],nodes[1]):
            Graph.AddEdge(nodes[0],nodes[1])
    return Graph
def genSmallWorld(N=5242, E=14484):
    Graph = genCircle(N)
    Graph = connectNbrOfNbr(Graph, N)
    Graph = connectRandomNodes(Graph, 4000)
    return Graph
def loadCollabNet(path):
    Graph = snap.LoadEdgeList(snap.PUNGraph, path, 0, 1, '\t')
    snap.DelSelfEdges(Graph)
    return Graph
def getDataPointsToPlot(Graph):
    odegree_list = []
    for NI in Graph.Nodes():
        odegree_list.append(NI.GetOutDeg())
    from collections import Counter
    cnt = Counter(odegree_list)
    X_cnt = cnt.keys()
    Y_cnt = cnt.values()
    X, Y = [], []
    for i in range(len(odegree_list)):
        X.append(i)
        if i in X_cnt:
            idx = X_cnt.index(i)
            Y.append(Y_cnt[idx])
        else:
            Y.append(0)
    return X, Y
def Q1_1():
    global erdosRenyi, smallWorld, collabNet
    erdosRenyi = genErdosRenyi(5242, 14484)
    smallWorld = genSmallWorld(5242, 14484)
    collabNet = loadCollabNet("ca-GrQc.txt")
    x_erdosRenyi, y_erdosRenyi = getDataPointsToPlot(erdosRenyi)
    plt.loglog(x_erdosRenyi, y_erdosRenyi, color = 'y', label = 'Erdos Renyi Network')
    x_smallWorld, y_smallWorld = getDataPointsToPlot(smallWorld)
    plt.loglog(x_smallWorld, y_smallWorld, linestyle = 'dashed', color = 'b', label = 'Small World Network')
    x_collabNet, y_collabNet = getDataPointsToPlot(collabNet)
    plt.loglog(x_collabNet, y_collabNet, linestyle = 'dotted', color = 'r', label = 'Collaboration Network')
    plt.xlabel('Node Degree (log)')
    plt.ylabel('Node Proportion with the Given Degree (log)')
    plt.title('Degree Distribution of Erdos Renyi, Small World, and Collaboration Networks')
    plt.legend()
    plt.show()
# Execute code for Q1.1
Q1_1()
# Problem 1.2
# Find max degree of all 3 graphs for plotting (add 2 for padding)
maxdeg = max([erdosRenyi.GetNI((snap.GetMxDegNId(erdosRenyi))).GetDeg(),
                smallWorld.GetNI((snap.GetMxDegNId(smallWorld))).GetDeg(),
                collabNet.GetNI((snap.GetMxDegNId(collabNet))).GetDeg()]) + 2
# Erdos Renyi
def calcQk(Graph, maxDeg=maxdeg):
    q_k = np.zeros(maxDeg)
    ks, ns = getDataPointsToPlot(Graph)
    total_nodes = np.sum(ns)
    ks ,ns = ks[:maxDeg], ns[:maxDeg]
    ps = [x/float(total_nodes) for x in ns]
    denom_k = np.sum([x*y for x,y in zip(ks,ps)])
    numerator = [x*y for x,y in zip(ks[1:],ps[1:])] # k+1*Pk+1
    numerator.append(0)
    q_k = numerator/denom_k
    return q_k

def calcExpectedDegree(Graph):
    ed = 0.0
    ks, ns = getDataPointsToPlot(Graph)
    total_nodes = Graph.GetNodes()
    ps = [x/float(total_nodes) for x in ns]
    ed = np.sum([x*y for x,y in zip(ks,ps)])
    return ed

def calcExpectedExcessDegree(Graph, qk):
    eed = 0.0
    eed = np.sum([x*y for x,y in zip(qk,range(len(qk)))])
    return eed

def Q1_2_a():
    print 'Output for Q1_2_a :'
    print ' '
    qk_erdosRenyi = calcQk(erdosRenyi, maxdeg)
    qk_smallWorld = calcQk(smallWorld, maxdeg)
    qk_collabNet = calcQk(collabNet, maxdeg)

    plt.loglog(range(maxdeg), qk_erdosRenyi, color = 'y', label = 'Erdos Renyi Network')
    plt.loglog(range(maxdeg), qk_smallWorld, linestyle = 'dashed', color = 'r', label = 'Small World Network')
    plt.loglog(range(maxdeg), qk_collabNet, linestyle = 'dotted', color = 'b', label = 'Collaboration Network')

    plt.xlabel('k Degree')
    plt.ylabel('Excess Degree Distribution')
    plt.title('Excess Degree Distribution of Erdos Renyi, Small World, and Collaboration Networks')
    plt.legend()
    plt.show()

    # Calculate Expected Degree
    ed_erdosRenyi = calcExpectedDegree(erdosRenyi)
    ed_smallWorld = calcExpectedDegree(smallWorld)
    ed_collabNet = calcExpectedDegree(collabNet)
    print 'Expected Degree for Erdos Renyi: %f' % ed_erdosRenyi
    print 'Expected Degree for Small World: %f' % ed_smallWorld
    print 'Expected Degree for Collaboration Network: %f' % ed_collabNet
    print ''
    print ''

    # Calculate Expected Excess Degree
    eed_erdosRenyi = calcExpectedExcessDegree(erdosRenyi, qk_erdosRenyi)
    eed_smallWorld = calcExpectedExcessDegree(smallWorld, qk_smallWorld)
    eed_collabNet = calcExpectedExcessDegree(collabNet, qk_collabNet)
    print 'Expected Excess Degree for Erdos Renyi: %f' % (eed_erdosRenyi)
    print 'Expected Excess Degree for Small World: %f' % (eed_smallWorld)
    print 'Expected Excess Degree for Collaboration Network: %f' % (eed_collabNet)
    print ''
    print ''

# Execute code for Q1.2a
Q1_2_a()

# Problem 1.3 - Clustering Coefficient 

def calcClusteringCoefficient(Graph):
    C = 0.0
    c_arrays =[]
    for NI in Graph.Nodes():
        k = NI.GetOutDeg()
        if k>=2:
            k_denom = k*(k-1)/2
            neibor_nodes =[]
            connected_neibor_edges = 0
            for Id in NI.GetOutEdges():
                neibor_nodes.append(Id)
            for id_n1 in range(len(neibor_nodes)):
                for id_n2 in range(id_n1+1,len(neibor_nodes)):
                    n1 = neibor_nodes[id_n1]
                    n2 = neibor_nodes[id_n2]
                    connected = Graph.IsEdge(n1,n2)
                    if connected:
                        connected_neibor_edges +=1
            c = connected_neibor_edges/float(k_denom)
            c_arrays.append(c)
        else:
            c_arrays.append(0)
    C = np.mean(c_arrays)
    return C

def Q1_3():
    C_erdosRenyi = calcClusteringCoefficient(erdosRenyi)
    C_smallWorld = calcClusteringCoefficient(smallWorld)
    C_collabNet = calcClusteringCoefficient(collabNet) 
    print 'Output for Q1_3:'
    print ' '
    print('Clustering Coefficient for Erdos Renyi Network: %f' % C_erdosRenyi)
    print('Clustering Coefficient for Small World Network: %f' % C_smallWorld)
    print('Clustering Coefficient for Collaboration Network: %f' % C_collabNet)
# Execute code for Q1.3
Q1_3()



