import networkx as nx
import pandas as pd
import numpy as np

def read_sequence(file):
    sequence = []
    with open(file,'r') as rf:
        while True:
            line = rf.readline()[:-1]
            if not line:
                break
            itemset = line.split(' ')
            sequence.append(itemset)
    return sequence

def read_block(file):
    blocks = {}
    blocks_count = {}
    with open(file, 'r') as rf:
        while True:
            line = rf.readline()[:-2]
            if not line:
                break
            itemset = line.split(' ')[1:]
            for i in itemset:
                if i.startswith('-'):
                    item = i[1:].split('_')
                    block = item[0]
                    block_object = item[1]
                else:
                    item = i.split('_')
                    block = item[0]
                    block_object = item[1]
                if block not in blocks_count.keys():
                    blocks_count[block] = 1
                else:
                    blocks_count[block] += 1
                if block not in blocks.keys():
                    blocks[block] = []
                    blocks[block].append([blocks_count[block],block_object])
                else:
                    blocks[block].append([blocks_count[block], block_object])
    return blocks


class GeneContentBuilder:

    def __init__(self,block_sequence,species,homo_threshold = 2):
        self.__non_dup_sequence = block_sequence
        self.__non_dup_sequence = self.__process_duplication()
        self.__non_dup_homo_sequence = self.__process_homology_species(species,homo_threshold)
        self.__block_dictionary = self.__build_block_dictionary()
        self.__ancestor_sequence = self.__get_final_block_seqence()

    def __process_duplication(self):
        print('process dup')
        filter_dup_block_sequence = {}
        for i in self.__non_dup_sequence.keys():
            filter_dup_block_sequence[i] = []
            for j in self.__non_dup_sequence[i].keys():
                filter_gene_familly = []
                gene_dup_count = {}
                for k in self.__non_dup_sequence[i][j][0]:
                    if k not in gene_dup_count.keys():
                        gene_dup_count[k] = 1
                    else:
                        gene_dup_count[k] += 1
                for k in gene_dup_count.keys():
                    if gene_dup_count[k] > 1:
                        filter_gene_familly.append(k)
                new_sequence = []
                new_sequence_name = []
                for k in range(len(self.__non_dup_sequence[i][j][0])):
                    if self.__non_dup_sequence[i][j][0][k] not in filter_gene_familly:
                        new_sequence.append(self.__non_dup_sequence[i][j][0][k])
                        new_sequence_name.append(self.__non_dup_sequence[i][j][1][k])
                filter_dup_block_sequence[i].append([new_sequence, new_sequence_name])
        return filter_dup_block_sequence

    def __process_homology_species(self,species,homo_threshold):
        print('process homo')
        filter_homo_sequence = {}
        for i in self.__non_dup_sequence.keys():
            gene_homo = {}
            count = 1
            for j in self.__non_dup_sequence[i]:
                for k in range(len(j[0])):
                    if j[0][k] not in gene_homo:
                        gene_homo[j[0][k]] = []
                        gene_homo[j[0][k]].append(species[count])
                    else:
                        if species[count] not in gene_homo[j[0][k]]:
                            gene_homo[j[0][k]].append(species[count])
                count += 1
            save_genes = []
            for j in gene_homo.keys():
                if len(gene_homo[j]) == homo_threshold:
                    save_genes.append(j)
            print(i+': '+ str(len(gene_homo)) + ' ' + str(len(save_genes)))
            filter_homo_sequence[i] = []
            for j in self.__non_dup_sequence[i]:
                new_sequence = []
                new_sequence_name = []
                for k in range(len(j[0])):
                    if j[0][k] in save_genes:
                        new_sequence.append(j[0][k])
                        new_sequence_name.append(j[1][k])
                filter_homo_sequence[i].append([new_sequence, new_sequence_name])
        return filter_homo_sequence

    def __build_block_dictionary(self):
        block_dictionary = {}
        for i in self.__non_dup_homo_sequence.keys():
            block_dictionary[i] = {}
            for j in self.__non_dup_homo_sequence[i]:
                sequence = j[0]
                name = j[1]
                for k in range(len(sequence)):
                    if sequence[k] not in block_dictionary[i].keys():
                        block_dictionary[i][sequence[k]] = [name[k]]
                    else:
                        block_dictionary[i][sequence[k]].append(name[k])
        return block_dictionary

    def __build_edges(self,sequence_list):
        gene_count = {}
        gene_rank = {}
        for j in sequence_list:
            for k in range(len(j)):
                if j[k] not in gene_count.keys():
                    gene_count[j[k]] = 1
                    gene_rank[j[k]] = k
                else:
                    gene_count[j[k]] += 1
                    gene_rank[j[k]] += k
        for j in gene_rank.keys():
            gene_rank[j] = gene_rank[j] / gene_count[j]
        gene_rank['E'] = 0
        for j in sequence_list:
            gene_rank['E'] += len(j)
        gene_rank['E'] = gene_rank['E'] / len(sequence_list)
        gene_count['E'] = len(sequence_list)
        edges_dir = {}
        for j in sequence_list:
            sequence = j
            if len(sequence) == 0:
                continue
            for k in range(len(sequence)):
                if k == 0:
                    continue
                else:
                    edge = (sequence[k - 1], sequence[k])
                    key = edge[0] + '@' + edge[1]
                    if key not in edges_dir.keys():
                        edges_dir[key] = 1
                    else:
                        edges_dir[key] += 1
            edge = (sequence[len(sequence) - 1], 'E')
            key = edge[0] + '@' + edge[1]
            if key not in edges_dir.keys():
                edges_dir[key] = 1
            else:
                edges_dir[key] += 1
        edges = []
        for j in edges_dir.keys():
            node = j.split('@')
            edges.append((node[0], node[1], edges_dir[j]))
        return edges,gene_count,gene_rank

    def __find_topological(self,edges,gene_count, gene_rank):
        DG = nx.DiGraph()
        DG.add_weighted_edges_from(edges)
        topological = []
        candidate = {}
        while True:
            while True:
                ok = 0
                degree = list(DG.in_degree)
                for i in degree:
                    if i[1] == 0:
                        ok = 1
                        topological.append(i[0])
                        for j in DG[i[0]]:
                            if j in candidate.keys():
                                candidate[j]['weight'] = candidate[j]['weight'] + DG[i[0]][j]['weight']
                            else:
                                candidate[j] = DG[i[0]][j]
                        DG.remove_node(i[0])
                if ok == 0:
                    break
            if len(degree) == 0:
                break

            if len(candidate.keys()) == 0:
                min_node = ''
                min_weight = 100000
                for i in gene_rank.keys():
                    if gene_rank[i] < min_weight:
                        min_node = i
                        min_weight = gene_rank[i]
                topological.append(min_node)
                for i in DG[min_node]:
                    if i in candidate.keys():
                        candidate[i]['weight'] = candidate[i]['weight'] + DG[min_node][i]['weight']
                    else:
                        candidate[i] = DG[min_node][i]
                DG.remove_node(min_node)
            else:
                new_candidate = {}
                for i in candidate.keys():
                    if i not in topological:
                        new_candidate[i] = candidate[i]
                candidate = new_candidate
                # print(degree)
                max_weight = -1
                max_weight_node = []
                for i in candidate.keys():
                    if candidate[i]['weight'] / gene_count[i] > max_weight:
                        max_weight = candidate[i]['weight'] / gene_count[i]
                        max_weight_node = [i]
                    elif candidate[i]['weight'] / gene_count[i] == max_weight:
                        max_weight_node.append(i)
                    else:
                        continue
                # print(max_weight_node)
                # print(max_weight)
                for i in max_weight_node:
                    topological.append(i)
                    next = DG[i]
                    for j in next.keys():
                        if j not in topological:
                            if j in candidate.keys():
                                candidate[j]['weight'] = candidate[j]['weight'] + next[j]['weight']
                            else:
                                candidate[j] = next[j]
                    # print(i)
                    DG.remove_node(i)
                # print(self.__topological)
                # print(candidate)
                #
                # print('______')
        return topological

    def __get_final_block_seqence(self):
        print('build ancestor block sequence')
        ancestor_sequence = {}
        for i in self.__non_dup_homo_sequence.keys():
            sequence = []
            for j in self.__non_dup_homo_sequence[i]:
                sequence.append(j[0])
            edges, gene_count, gene_rank = self.__build_edges(sequence)
            topological = self.__find_topological(edges, gene_count, gene_rank)
            new_sequence = []
            new_sequence_name = []
            for j in range(len(topological)):
                if topological[j] == 'E':
                    continue
                new_sequence.append(topological[j])
                new_sequence_name.append(self.__block_dictionary[i][topological[j]][-1])
            ancestor_sequence[i] = [new_sequence, new_sequence_name]
        return ancestor_sequence

    def out_synteny_file(self,outsynteny,outsynteny_name,support,block_suffix=''):
        out_synteny_file = open(outsynteny, 'w')
        out_synteny_name_file = open(outsynteny_name, 'w')
        for i in self.__ancestor_sequence.keys():
            sequence = self.__ancestor_sequence[i][0]
            sequence_name = self.__ancestor_sequence[i][1]
            out_synteny_file.write(i+block_suffix+':'+str(support)+' ')
            out_synteny_name_file.write(i + block_suffix + ':' + str(support) + ' ')
            for j in sequence:
                out_synteny_file.write(j + ' ')
            for j in sequence_name:
                out_synteny_name_file.write(j + ' ')
            out_synteny_file.write('\n')
            out_synteny_name_file.write('\n')
        out_synteny_name_file.close()
        out_synteny_file.close()

    def out_block_dictionary(self,out_block_gene_name):
        out_block_gene_name = open(out_block_gene_name,'w')
        for i in self.__block_dictionary.keys():
            for j in self.__block_dictionary[i].keys():
                out_block_gene_name.write(i+'\t'+j+'\t')
                for k in self.__block_dictionary[i][j]:
                    out_block_gene_name.write(k+' ')
                out_block_gene_name.write('\n')
        out_block_gene_name.close()


    def get_ancestor_sequence(self):
        return self.__ancestor_sequence








