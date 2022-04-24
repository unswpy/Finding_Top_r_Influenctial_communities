from xml.dom.minidom import parse

import networkx as nx
import xml.dom.minidom

DOMTree = xml.dom.minidom.parse("origin/dblp.xml")
collection = DOMTree.documentElement
print("======read data======")
articles = collection.getElementsByTagName("article")

count = 0
edge_count = 0
name_to_id_dict = {}  # author name to author id
id_to_name_dict = {}  # author id to author name
edge_dict = {}  # edge
with open("origin/dblp_case.txt", "w") as file:
    pass
for article in articles:
    author = article.getElementsByTagName('author')
    if len(author) > 1:
        num_list = []
        for num in range(len(author)):
            name = author[num].childNodes[0].data
            if name in name_to_id_dict.keys():
                num_list.append(name_to_id_dict[name])
            else:
                name_to_id_dict[name] = count
                id_to_name_dict[count] = name
                count = count + 1
                num_list.append(name_to_id_dict[name])
        for i in range(len(num_list)):
            for j in range(len(num_list)):
                if num_list[i] in edge_dict.keys():
                    edge_dict[num_list[i]].add(num_list[j])
                else:
                    edge_dict[num_list[i]] = set()
                    edge_dict[num_list[i]].add(num_list[j])
                if num_list[j] in edge_dict.keys():
                    edge_dict[num_list[j]].add(num_list[i])
                else:
                    edge_dict[num_list[j]] = set()
                    edge_dict[num_list[j]].add(num_list[i])
    # print(len(author))
    # print(author[0].childNodes[0].data)

with open("origin/dblp_case.txt", "a") as file:
    for start, end_set in edge_dict.items():
        for end in iter(end_set):
            if start < end:
                file.write(str(start) + '\t' + str(end) + '\n')
                edge_count = edge_count + 1

with open("origin/dblp_case.txt", 'r+') as file:
    content = file.read()
    file.seek(0, 0)
    file.write('#\n' + '#\n' + '# Nodes: ' + str(count) + ' Edges: ' + str(edge_count) + '\n' + '#\n' + content)

with open("origin/dblp_name_list.txt", "w") as file:
    for i in range(count):
        file.write(str(i) + '   ' + str(id_to_name_dict[i]) + '\n')

