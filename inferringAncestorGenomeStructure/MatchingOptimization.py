import gurobipy as gp
from gurobipy import *
import pandas as pd
import numpy as np

class MatchingOptimization:

    def __init__(self,file_list,matching_dim1 = 4,matching_dim2 = 2,relation1 = 1,relation2 = 2,self_match = False):
        self.__endpoint_lists = []
        self.__matching_dim1 = matching_dim1
        self.__matching_dim2 = matching_dim2
        self.__relation1 = relation1
        self.__relation2 = relation2
        self.__relabel_block_sequences = []
        self.__compress_adjacency_matrixs = []
        self.__self_match = self_match
        for i in file_list:
            adjacency_list, relabel_block_order = self.__assumed_block_label(i)
            self.__relabel_block_sequences.append(relabel_block_order)
            compress_adjacency_matrix, endpoint_list = self.__build_assumed_matrix(adjacency_list)
            self.__compress_adjacency_matrixs.append(compress_adjacency_matrix)
            self.__endpoint_lists.append(endpoint_list)

        candidate_compress_adjacency_matrix = self.__compress_adjacency_matrixs[0]
        guided_compress_adjacency_matrix = self.__compress_adjacency_matrixs[1]
        candidate_adjacency_index = self.__endpoint_lists[0]
        self.__match_pairs = []

        for i in candidate_compress_adjacency_matrix:
            match_pair = []
            adj1 = i[-1]
            compare_key1 = ''
            if adj1[0] != '$':
                item = adj1[0].split('@')
                compare_key1 += item[0]
            if adj1[1] != '$':
                item = adj1[1].split('@')
                compare_key1 += item[0]

            for j in guided_compress_adjacency_matrix:
                adj2 = j[-1]
                compare_key2 = ''
                if adj2[0] != '$':
                    item = adj2[0].split('@')
                    compare_key2 += item[0]
                if adj2[1] != '$':
                    item = adj2[1].split('@')
                    compare_key2 += item[0]
                if compare_key1 == compare_key2:
                    match_pair.append(j)
            self.__match_pairs.append([i, match_pair])
        self.__k = int((len(candidate_adjacency_index) - 1) / (self.__matching_dim1*2))

    def optimization(self):
        try:
            self.__m = gp.Model()
            match_matrix = self.__m.addVars(self.__k,
                                            self.__matching_dim1,
                                            self.__matching_dim2,
                                            vtype=GRB.BINARY,
                                            name="matching_matrix")
            self.__m.update()
            self.__m.setObjective(gp.quicksum(
                (i[0][2] * match_matrix[int(i[0][0] / (self.__matching_dim1*2)),
                                        i[0][0] % self.__matching_dim1,
                                        j[0] % self.__matching_dim2] + (1 - i[0][2])) *
                (i[0][3] * match_matrix[int(i[0][1] / (self.__matching_dim1*2)),
                                        i[0][1] % self.__matching_dim1,
                                        j[1] % self.__matching_dim2] + 1 - i[0][3])
                for i in self.__match_pairs for j in i[1]
            ), GRB.MAXIMIZE)

            self.__m.addConstrs((
                gp.quicksum(match_matrix[i, j, l] for l in range(self.__matching_dim2)) == self.__relation1
                for i in range(self.__k)
                for j in range(self.__matching_dim1)), name='row_unique'
            )
            self.__m.addConstrs((
                gp.quicksum(match_matrix[i, l, j] for l in range(self.__matching_dim1)) == self.__relation2
                for i in range(self.__k)
                for j in range(self.__matching_dim2)), name='col_unique'
            )
            if self.__self_match:

                self.__m.addConstrs((
                    match_matrix[i,j,j] == 0
                    for i in range(self.__k)
                    for j in range(self.__matching_dim1)), name='diagonal'
                )

                self.__m.addConstrs((
                    match_matrix[i, j, k] == match_matrix[i, k, j]
                    for i in range(self.__k)
                    for j in range(self.__matching_dim1)
                    for k in range(self.__matching_dim2)), name='symmetry'
                )

            self.__m.optimize()
            print('Obj: %g' % self.__m.objVal)
        except gp.GurobiError as e:
            print('Error code ' + str(e.errno) + ': ' + str(e))
        except AttributeError:
            print('Encountered an attribute error')

    def matching_relation(self):
        result = []
        for v in self.__m.getVars():
            result.append(v.x)
        result = np.asarray(result)
        result = result.reshape((self.__k, self.__matching_dim1, self.__matching_dim2))
        self.__match_relations = {}
        for i in range(len(result)):
            column = []
            index = []
            for j in range(self.__matching_dim2):
                item = self.__endpoint_lists[1][i * self.__matching_dim2*2 + 1].split('@')
                column.append(item[0][:-1] + '@' + str(j + 1))
            for j in range(self.__matching_dim1):
                item = self.__endpoint_lists[1][i * self.__matching_dim2 * 2 + 1].split('@')
                index.append(item[0][:-1] + '@' + str(j + 1))
            match = pd.DataFrame(result[i], columns=column, index=index)
            match_relation = match.to_dict()
            for j in match_relation.keys():
                for k in match_relation[j].keys():
                    if match_relation[j][k] == 1:
                        self.__match_relations[k] = j

    def output_matching_relation(self,outfile):
        outfile = open(outfile, 'w')
        for i in self.__match_relations.keys():
            key1 = i.split('@')
            key2 = self.__match_relations[i].split('@')
            # block dim1 dim2
            outfile.write(key1[0] + ' ' + key1[1] + ' ' + key2[1] + '\n')
        outfile.close()

    def output_new_sequence(self,candidate_file,guided_file):
        transform_block_orders = []
        for i in self.__relabel_block_sequences[0]:
            transform_block_order = []
            chr_type = i[0]
            for j in i[1:]:
                if j.startswith('-'):
                    # key dim1 value dim2
                    item1 = j[1:].split('_')
                    key = item1[0] + '@' + item1[1]
                    item2 = self.__match_relations[key].split('@')
                    transform_block_order.append('-' + item2[0]+'_'+item2[1])
                else:
                    item1 = j.split('_')
                    key = item1[0] + '@' + item1[1]
                    item2 = self.__match_relations[key].split('@')
                    transform_block_order.append(item2[0] + '_' + item2[1])
            transform_block_orders.append([chr_type]+transform_block_order)
        candidate_file = open(candidate_file, 'w')
        guided_file = open(guided_file, 'w')
        for i in transform_block_orders:
            line = ''
            for j in i:
                line += j + ' '
            line += '\n'
            candidate_file.write(line)
        candidate_file.close()
        for i in self.__relabel_block_sequences[1]:
            line = ''
            for j in i:
                line += j + ' '
            line += '\n'
            guided_file.write(line)
        guided_file.close()

    def __assumed_block_label(self, file):
        adjacency_list = []
        block_objects = {}
        relabel_block_order = []
        with open(file) as df:
            while True:
                line = df.readline()[:-2]
                if not line:
                    break
                itemset = line.split(' ')
                chr_type = itemset[0]
                new_block_order = []
                for i in itemset[1:]:
                    block = ''
                    if i.startswith('-'):
                        block = i[1:]
                        if block not in block_objects.keys():
                            block_objects[block] = 1
                            new_block = '-' + block + '_' + str(block_objects[block])
                            block_objects[block] += 1
                        else:
                            new_block = '-' + block + '_' + str(block_objects[block])
                            block_objects[block] += 1
                    else:
                        block = i
                        if block not in block_objects.keys():
                            block_objects[block] = 1
                            new_block = block + '_' + str(block_objects[block])
                            block_objects[block] += 1
                        else:
                            new_block = block + '_' + str(block_objects[block])
                            block_objects[block] += 1
                    new_block_order.append(new_block)
                last = ''
                start = ''
                for i in range(len(new_block_order)):
                    if i == 0:
                        if chr_type == 's':
                            if new_block_order[i].startswith('-'):
                                block = new_block_order[i][1:].split('_')[0]
                                copy_number = new_block_order[i][1:].split('_')[1]
                                adjacency_list.append(['$', block + 'b' + '@' + copy_number])
                                last = block + 'a' + '@' + copy_number
                            else:
                                block = new_block_order[i].split('_')[0]
                                copy_number = new_block_order[i].split('_')[1]
                                adjacency_list.append(['$', block + 'a' + '@' + copy_number])
                                last = block + 'b' + '@' + copy_number
                        else:
                            if new_block_order[i].startswith('-'):
                                block = new_block_order[i][1:].split('_')[0]
                                copy_number = new_block_order[i][1:].split('_')[1]
                                last = block + 'a' + '@' + copy_number
                                start = block + 'b' + '@' + copy_number
                            else:
                                block = new_block_order[i].split('_')[0]
                                copy_number = new_block_order[i].split('_')[1]
                                last = block + 'b' + '@' + copy_number
                                start = block + 'a' + '@' + copy_number
                    else:
                        if new_block_order[i].startswith('-'):
                            block = new_block_order[i][1:].split('_')[0]
                            copy_number = new_block_order[i][1:].split('_')[1]
                            adjacency_list.append([last, block + 'b' + '@' + copy_number])
                            last = block + 'a' + '@' + copy_number
                        else:
                            block = new_block_order[i].split('_')[0]
                            copy_number = new_block_order[i].split('_')[1]
                            adjacency_list.append([last, block + 'a' + '@' + copy_number])
                            last = block + 'b' + '@' + copy_number
                if chr_type == 's':
                    adjacency_list.append([last, '$'])
                else:
                    adjacency_list.append([last,start])
                relabel_block_order.append([chr_type] + new_block_order)
        return adjacency_list, relabel_block_order

    def __build_assumed_matrix(self,adjacency_list):
        endpoint_list = []
        for i in adjacency_list:
            for j in i:
                if j not in endpoint_list:
                    endpoint_list.append(j)
        endpoint_list = sorted(endpoint_list)
        adjacency_matrix = {}
        for i in endpoint_list:
            adjacency_matrix[i] = {}
            for j in endpoint_list:
                adjacency_matrix[i][j] = 0
        for i in adjacency_list:
            adjacency_matrix[i[0]][i[1]] += 1
            adjacency_matrix[i[1]][i[0]] += 1
        adjacency_matrix = pd.DataFrame(adjacency_matrix)
        adjacency_matrix = np.asarray(adjacency_matrix)
        compress_adjacency_matrix = []
        for i in range(len(adjacency_matrix)):
            for j in range(len(adjacency_matrix[i])):
                if adjacency_matrix[i][j] == 1:
                    adjacency = [endpoint_list[i], endpoint_list[j]]
                    adjacency = sorted(adjacency)
                    if adjacency[0] == endpoint_list[i] and adjacency[1] == endpoint_list[j]:
                        if i == 0 and j == 0:
                            compress_adjacency_matrix.append([i, j, 0, 0, adjacency])
                        if i == 0 and j != 0:
                            compress_adjacency_matrix.append([i, j - 1, 0, 1, adjacency])
                        if i != 0 and j == 0:
                            compress_adjacency_matrix.append([i - 1, j, 1, 0, adjacency])
                        if i != 0 and j != 0:
                            compress_adjacency_matrix.append([i - 1, j - 1, 1, 1, adjacency])
                    if adjacency[0] == endpoint_list[j] and adjacency[1] == endpoint_list[i]:
                        if i == 0 and j == 0:
                            compress_adjacency_matrix.append([j, i, 0, 0, adjacency])
                        if i == 0 and j != 0:
                            compress_adjacency_matrix.append([j - 1, i, 1, 0, adjacency])
                        if i != 0 and j == 0:
                            compress_adjacency_matrix.append([j, i - 1, 0, 1, adjacency])
                        if i != 0 and j != 0:
                            compress_adjacency_matrix.append([j - 1, i - 1, 1, 1, adjacency])
        return compress_adjacency_matrix, endpoint_list



