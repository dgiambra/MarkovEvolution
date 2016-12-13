# -*- coding: utf-8 -*-
def evolution(gen_list,generation):
    gen=[]
    for i in gen_list:
        a = random.randint(0, 3)
        if a == 0 & generation == 1:
            a == 2
        if a==0 & generation != 1:
            if len(gen_list)==1:
                a = 2
        for x in range(a):
                species = list(i)
                for z in range(len(species)):
                    a = random.randint(0, 4)
                    if a == 4:
                        ##print(species[z])
                        species[z] = random.choice(amino_acids)
                gen.append("".join(species))
                
    return gen

def overall(gen_1, totalgenerations, gen = 1):
    tree = [gen_1]
    for i in range(totalgenerations):
        tree.append(evolution(tree[i],gen))
        gen += 1
    return tree

def weighter(tree):
    x = 0
    dict1 = {}
    with pymp.Parallel(4) as p:
        for i in tree:
            for sequence in i:
                dict1[x] = []
                for h in tree:
                    for sequence2 in h:
                        dict1[x].append(Bio.pairwise2.align.globalxx(sequence, sequence2)[0][4])
                x +=1
    return dict1
                
                

def grapher(N, r, T, gens):
    '''This function is a Markoc Chain Monte Carlo Simulator.
    Parameters
    ----------
     gens : int
            number of generations
        N : int
            number of iterations
        r : int
            adjustable parameter
        T : int
            adjustable Parameters
    Returns
    -------
    graphs : list
             most likely graphs from each iteration
    '''
    import networkx as nx
    from random import choice, random
    import math
    amino_acids = ['G','C','A','T',"-"]
    luca = "GCATTACTAATATATTAGTAAATCAGAGTAGTA"
    gen_1 = [luca]
    tree = overall(gen_1,gens)
    dict1 = weighter(tree)
    x = 0 
    G = nx.Graph()
    with pymp.Parallel(4) as p:
        for i in tree:
                for sequence in i:
                    G.add_node(x)
                    x+=1
    with pymp.Parallel(4) as p:
        for i in G.nodes():
            for x in G.nodes():
                G.add_edge(i,x,weight = dict1[i][x])
    G=nx.minimum_spanning_tree(G)

    graphs = []
    for i in list(range(N)):
        H = nx.Graph()
        I = nx.Graph()
        J= nx.Graph()
        with pymp.Parallel(4) as p:
            for edge in nx.edges(G):
                H.add_edge(edge[0], edge[1], weight = G[edge[0]][edge[1]]['weight'])
                I.add_edge(edge[0], edge[1], weight = G[edge[0]][edge[1]]['weight'])

        
        random_edge = choice(H.edges())
        random_edge_2 = choice(H.edges())
        if random_edge == random_edge_2:
            random_edge_2 = choice(G.edges())
        start_node = random_edge[0]
        start_node_2 = random_edge_2[0]
        end_node = random_edge[1]
        end_node_2 = random_edge_2[1]
        H.remove_edge(random_edge[0],random_edge[1])
        H.remove_edge(random_edge_2[0],random_edge_2[1])
        H.add_edge(start_node, end_node_2, weight = dict1[start_node][end_node_2])
        H.add_edge(start_node_2, end_node, weight = dict1[start_node_2][end_node])
        for edge in nx.edges(H):
            J.add_edge(edge[0], edge[1], weight = H[edge[0]][edge[1]]['weight'])
        
        
        
        
        theta_g = theta(G, r)
        theta_h = theta(H, r)
        fx = math.exp(-(theta_g-theta_h)/T)
        b_g = 0
        b_h = 0
        for i in G.edges():
            I.remove_edge(i[0],i[1])
            if not nx.is_connected(I):
                b_g +=1
                I.add_edge(i[0], i[1])
            else:
                I.add_edge(i[0], i[1])
        for i in H.edges():
            J.remove_edge(i[0],i[1])
            if not nx.is_connected(J):
                b_h +=1
                J.add_edge(i[0], i[1])
            else:
                J.add_edge(i[0], i[1])
        q_h_g = 1/(len(G.nodes())*(len(G.nodes())-1)/2-b_h)
        q_g_h = 1/(len(G.nodes())*(len(G.nodes())-1)/2-b_g)
        # print lines for classical debugging
        # print(1/(len(x)*(len(x)-1)/2-b_h))
        # print(len(x))
        # print(b_h)
        # print(q_h_g)
        aij = fx*q_h_g/q_g_h
        # print(nx.info(G))
        # print(nx.info(H))
        #print(aij)
        # print(fx)
        if aij > random():
            graphs.append(H)
            G = H
        else:
            graphs.append(G)
    return graphs

def topOnePercent(graphs):
    graph_dict = {}
    for i in graphs:
        if i in graph_dict:
            graph_dict[i] += 1
        else: 
            graph_dict[i] = 1
    num = int(len(graph_dict)/100)
    import operator
    sorted_x = sorted(graph_dict.items(), key=operator.itemgetter(1))
    # print(len(graph_dict))
    # print(num)
    top_one_percent = sorted_x[-num:]
    # print(top_one_percent)
    return top_one_percent

def expectedDegree(graphs, node):
    from statistics import mean
    degreeList = []
    for i in graphs:
        degreeList.append(i.degree(node))
    return mean(degreeList)

def expectedNumberofEdges(graphs):
    from statistics import mean
    edges = []
    for i in graphs:
        edges.append(i.number_of_edges())
    return mean(edges)

def expectedShortestPathLength(graphs, node_a, node_b):
    from statistics import mean
    import networkx as nx
    pathList=[]
    for i in graphs:
        pathList.append(nx.dijkstra_path_length(i, node_a, node_b))
    return mean(pathList)
    
    
def theta(G, r):
    '''This function calculates the relative probability parameter for a graph.
       Called in grapher function.
    Parameters
    ----------
        G : graph
            graph to be analyzed
        r : int
            adjustable parameter
    Returns
    -------
    theta : float
            parameter for probability calculation
    '''
    import networkx as nx
    theta = 0
    for item in G.edges(data='weight'):
        y = item[2]
        theta += r*y
    theta += len(G.edges())
    return theta
