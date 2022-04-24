import networkx as nx

name_list = ['com-dblp', 'com-lj', 'com-orkut', 'com-youtube']

for name in name_list:
    G = nx.Graph()
    max_core = 0
    d_sum = 0
    d_avg = 0.0
    d_max = 0
    with open("origin/{}.ungraph.txt".format(name), 'r') as file:
        id_dict = {}
        count_id = 0
        count = 0
        for line in file.readlines():
            if count > 3:
                if int(line.split()[0]) not in id_dict.keys():
                    start = count_id
                    id_dict[int(line.split()[0])] = count_id
                    count_id += 1
                else:
                    start = id_dict[int(line.split()[0])]
                if int(line.split()[1]) not in id_dict.keys():
                    end = count_id
                    id_dict[int(line.split()[1])] = count_id
                    count_id += 1
                else:
                    end = id_dict[int(line.split()[1])]
                G.add_edge(start, end)
            count = count + 1
    pr = nx.pagerank(G, alpha=0.85)
    core = nx.core_number(G)
    nodes_list = list(G.nodes)
    nodes_list = [int(x) for x in nodes_list]
    nodes_list.sort()
    dict_deg = dict(G.degree())
    for i in range(len(G.nodes)):
        if max_core < core[nodes_list[i]]:
            max_core = core[nodes_list[i]]
        if d_max < dict_deg[i]:
            d_max = dict_deg[i]
        d_sum += dict_deg[i]
    d_avg = d_sum * 1.0 / len(nodes_list)
    with open("preprocess/{}.txt".format(name), 'w') as file:
        edges_list = list(G.edges)
        file.write(str(len(G.nodes())) + '  ' +
                   str(len(G.edges())) + '  ' +
                   str(d_max) + '  ' +
                   str(d_avg) + '  ' +
                   str(max_core) + '\n')
        for i in range(len(G.nodes())):
            file.write(str(nodes_list[i]) + '  ' +
                       '{:.10f}'.format(pr[nodes_list[i]]) + '  ' +
                       str(core[nodes_list[i]]) + '\n')
        for i in range(len(G.edges())):
            file.write(str(edges_list[i][0]) + '  ' +
                       str(edges_list[i][1]) + '\n')


text_list = ['Email-Enron']

for text in text_list:
    G = nx.Graph()
    max_core = 0
    d_sum = 0
    d_avg = 0.0
    d_max = 0
    with open("origin/{}.txt".format(text), 'r') as file:
        id_dict = {}
        count_id = 0
        count = 0
        for line in file.readlines():
            if count > 3:
                if int(line.split()[0]) not in id_dict.keys():
                    start = count_id
                    id_dict[int(line.split()[0])] = count_id
                    count_id += 1
                else:
                    start = id_dict[int(line.split()[0])]
                if int(line.split()[1]) not in id_dict.keys():
                    end = count_id
                    id_dict[int(line.split()[1])] = count_id
                    count_id += 1
                else:
                    end = id_dict[int(line.split()[1])]
                G.add_edge(start, end)
            count = count + 1
    pr = nx.pagerank(G, alpha=0.85)
    core = nx.pagerank(G)
    nodes_list = list(G.nodes)
    nodes_list = [int(x) for x in nodes_list]
    nodes_list.sort()
    dict_deg = dict(G.degree())
    for i in range(len(G.nodes)):
        if max_core < core[nodes_list[i]]:
            max_core = core[nodes_list[i]]
        if d_max < dict_deg[i]:
            d_max = dict_deg[i]
        d_sum += dict_deg[i]
    d_avg = d_sum * 1.0 / len(nodes_list)
    with open("preprocess/{}.txt".format(text.split('-')[0]), 'w') as file:
        edges_list = list(G.edges)
        file.write(str(len(G.nodes())) + '  ' +
                   str(len(G.edges())) + '  ' +
                   str(d_max) + '  ' +
                   str(d_avg) + '  ' +
                   str(max_core) + '\n')
        for i in range(len(G.nodes())):
            file.write(str(nodes_list[i]) + '  ' +
                       '{:.10f}'.format(pr[nodes_list[i]]) + '  ' +
                       str(core[nodes_list[i]]) + '\n')
        for i in range(len(G.edges())):
            file.write(str(edges_list[i][0]) + '  ' + str(edges_list[i][1]) + '\n')
