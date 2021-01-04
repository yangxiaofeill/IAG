import gurobipy as gp
from gurobipy import *
import pandas as pd
import numpy as np
import copy
from dataSturcture.adjacencyMatrix import AdjacencyMatrix
from util.transformToAdjacency import transformToAdjacency

class GGAP:
    def __init__(self, adjacency_file,species_list,duptype=2):
        self.species_list = species_list
        self.adjacency_file = adjacency_file
        self.__duptype = duptype
        # 观测的值向量，
        # 邻接的对应名字，
        # 邻接的分区range信息，
        # 邻接的对称信息，
        # 祖先邻接的候选集合
        # 对角信息，
        # 矩阵的列名信息
        self.__observation_adjacency_vectors_value, \
        self.__adjacency_name, \
        self.__vector_range_value, \
        self.__vector_symmetry_value, \
        self.__candidate_adjacency, \
        self.__diagonal_value, \
        self.__matrix_items = self.__building_compress_calculation_matrix()
        self.__variable_number = len(self.__observation_adjacency_vectors_value[0])
        self.__dim = len(self.__vector_range_value)

    def evaluation_ancestor_adjacencies(self,outdir):
        species_adjacencies = {}
        for i in range(len(self.species_list)):
            species_adjacencies[self.species_list[i]] = {}
            for j in range(len(self.__adjacency_name)):
                if self.__observation_adjacency_vectors_value[i][j] != 0:
                    key = self.__adjacency_name[j].split('@')
                    key = sorted(key)
                    key = key[0] + '@' + key[1]
                    if key not in species_adjacencies[self.species_list[i]].keys():
                        species_adjacencies[self.species_list[i]][key] = self.__observation_adjacency_vectors_value[i][j]
                    else:
                        species_adjacencies[self.species_list[i]][key] += \
                        self.__observation_adjacency_vectors_value[i][j]

        np_adjacency_matrix = np.asarray(self.__adjacency_matrix)
        ancestor_adjacencies = {}
        for i in range(len(self.__matrix_items)):
            for j in range(len(self.__matrix_items)):
                if np_adjacency_matrix[i][j] != 0:
                    key = sorted([self.__matrix_items[i],self.__matrix_items[j]])
                    key = key[0] + '@' + key[1]
                    if key not in ancestor_adjacencies.keys():
                        ancestor_adjacencies[key] = int(np_adjacency_matrix[i][j])
                    else:
                        ancestor_adjacencies[key] += int(np_adjacency_matrix[i][j])

        evaluation_list = {}
        for i in ancestor_adjacencies.keys():
            evaluation_list[i] = {}
            for j in species_adjacencies.keys():
                if i in species_adjacencies[j].keys():
                    evaluation_list[i][j] = int(species_adjacencies[j][i]/2)
                else:
                    evaluation_list[i][j] = 0
        core_adjacencies = []
        homo_adjacencies = []
        possible_adjacencies = []
        infer_adjacencies = []
        for i in evaluation_list.keys():
            count = 0
            species = ''
            for j in evaluation_list[i].keys():
                if evaluation_list[i][j] != 0:
                    count += evaluation_list[i][j]
                    species+=j+'_'
            species = species[:-1]
            if count == 0:
                infer_adjacencies.append([i,count,species])
            if count == 1:
                possible_adjacencies.append([i,count,species])
            if count > 1 and count < self.__duptype + 1:
                homo_adjacencies.append([i,count,species])
            if count == self.__duptype + 1:
                core_adjacencies.append([i,count,species])
        evaluation_list = pd.DataFrame(evaluation_list)

        infer_adjacencies = pd.DataFrame(infer_adjacencies, columns=['adjacency', 'number', 'species'])
        infer_adjacencies.to_csv(outdir + 'infer_adjacencies.xls', sep='\t', index=None)
        core_adjacencies = pd.DataFrame(core_adjacencies,columns=['adjacency','number','species'])
        core_adjacencies.to_csv(outdir+'core_adjacencies.xls',sep='\t',index=None)
        homo_adjacencies = pd.DataFrame(homo_adjacencies, columns=['adjacency', 'number','species'])
        homo_adjacencies.to_csv(outdir+'homo_adjacencies.xls',sep='\t',index=None)
        possible_adjacencies = pd.DataFrame(possible_adjacencies, columns=['adjacency', 'number','species'])
        possible_adjacencies.to_csv(outdir+'possible_adjacencies.xls',sep='\t',index=None)
        evaluation_list.to_csv(outdir+'adjacencies_evaluation.xls',sep='\t')


    def __building_compress_calculation_matrix(self):
        adjacency_relations = []
        matrix_items = []
        with open(self.adjacency_file) as af:
            while True:
                adjacency = []
                line = af.readline()[:-2]
                if not line:
                    break
                itemset = line.split(' ')
                index = 0
                while index < len(itemset):
                    adjacency.append([itemset[index], itemset[index + 1]])
                    if itemset[index] not in matrix_items:
                        matrix_items.append(itemset[index])
                    if itemset[index + 1] not in matrix_items:
                        matrix_items.append(itemset[index + 1])
                    index += 2
                adjacency_relations.append(adjacency)
        matrix_items = sorted(matrix_items)
        # 先构建大矩阵，
        # 然后生成值向量
        # 表示每个行每个值包括的邻接对的向量，
        # 同时记录对称位置向量，
        # 记录每个项起止位置的向量
        adjacency_matrix = {}
        for i in matrix_items:
            adjacency_matrix[i] = {}
            for j in matrix_items:
                adjacency_matrix[i][j] = 0
        observation_adjacency_matrixs = []
        for i in adjacency_relations:
            observation_adjacency_matrix = copy.deepcopy(adjacency_matrix)
            for j in i:
                observation_adjacency_matrix[j[0]][j[1]] += 1
                observation_adjacency_matrix[j[1]][j[0]] += 1
            observation_adjacency_matrixs.append(observation_adjacency_matrix)
        # 遍历两个生成完的矩阵，构建并集向量
        candidate_observation_adjacency_matrix = observation_adjacency_matrixs[0]
        guided_observation_adjacency_matrix = observation_adjacency_matrixs[1]

        # 每个端点的邻接表格
        adjacency_table = {}
        for i in matrix_items:
            if i == '$':
                union_items = []
                for j in matrix_items:
                    union_items.append(j)
            else:
                union_items = []
                union_items.append('$')
                for j in candidate_observation_adjacency_matrix[i].keys():
                    if candidate_observation_adjacency_matrix[i][j] != 0:
                        if j not in union_items:
                            union_items.append(j)
                for j in guided_observation_adjacency_matrix[i].keys():
                    if guided_observation_adjacency_matrix[i][j] != 0:
                        if j not in union_items:
                            union_items.append(j)

            adjacency_table[i] = {}
            for j in union_items:
                adjacency_table[i][j] = 0
        # 分区
        range_vector = {}
        adjacency_vector = {}
        start = 0
        count = 0
        for i in adjacency_table.keys():
            range_vector[i] = []
            range_vector[i].append(start)
            end = start
            for j in adjacency_table[i].keys():
                key = i + '@' + j
                adjacency_vector[key] = count
                end += 1
                count += 1
            range_vector[i].append(end)
            start = end
        # 对称vector，表示对称点的位置
        symmetry_vector = {}
        for i in adjacency_vector.keys():
            key = i.split('@')
            key_sym = key[1] + '@' + key[0]
            symmetry_vector[i] = adjacency_vector[key_sym]

        # 清空adjacency_vector
        empty_adjacency_vector = {}
        for i in adjacency_vector:
            empty_adjacency_vector[i] = 0
        observation_adjacency_vectors = []
        for i in observation_adjacency_matrixs:
            observation_adjacency_vector = copy.deepcopy(empty_adjacency_vector)
            for j in i.keys():
                for k in i[j].keys():
                    if i[j][k] != 0:
                        key = j + '@' + k
                        observation_adjacency_vector[key] += i[j][k]
            observation_adjacency_vectors.append(observation_adjacency_vector)
        # 数值向量化
        adjacency_name = list(observation_adjacency_vectors[0].keys())
        observation_adjacency_vectors_value = []
        for i in observation_adjacency_vectors:
            observation_adjacency_vector_value = []
            for j in i.keys():
                observation_adjacency_vector_value.append(i[j])
            observation_adjacency_vectors_value.append(observation_adjacency_vector_value)

        candidate_adjacency = []
        vector_range_value = []
        for i in range_vector.keys():
            nonZero = []
            ranges = range_vector[i]
            vector_range_value.append(ranges)
            sequence = observation_adjacency_vectors_value[0][ranges[0]:ranges[1]]
            for j in range(len(sequence)):
                if sequence[j] != 0:
                    nonZero.append(j + ranges[0])
            if ranges[0] not in nonZero:
                nonZero.append(ranges[0])
            candidate_adjacency.append(nonZero)
        vector_symmetry_value = []
        for i in symmetry_vector.keys():
            vector_symmetry_value.append(symmetry_vector[i])
        # 检查对称等于自己的是错误
        diagonal_value = []
        for i in range(len(vector_symmetry_value)):
            if vector_symmetry_value[i] == i:
                diagonal_value.append(i)
        return observation_adjacency_vectors_value, \
               adjacency_name, \
               vector_range_value, \
               vector_symmetry_value, \
               candidate_adjacency, \
               diagonal_value, \
               matrix_items

    def optimization(self):
        try:
            self.__alpha = self.__duptype - 1
            self.__m = gp.Model()
            # 定义整数，需要添加范围约束
            ancestor = self.__m.addVars(self.__variable_number,
                                        vtype=GRB.BINARY, name="ancestor")
            self.__m.update()
            self.__m.setObjective(gp.quicksum(
                self.__alpha *
                (self.__duptype * ancestor[i + self.__vector_range_value[0][1]] -
                 self.__observation_adjacency_vectors_value[0][i + self.__vector_range_value[0][1]]) *
                (self.__duptype * ancestor[i + self.__vector_range_value[0][1]] -
                 self.__observation_adjacency_vectors_value[0][i + self.__vector_range_value[0][1]]) +
                (ancestor[i + self.__vector_range_value[0][1]] -
                 self.__observation_adjacency_vectors_value[1][i + self.__vector_range_value[0][1]]) *
                (ancestor[i + self.__vector_range_value[0][1]] -
                 self.__observation_adjacency_vectors_value[1][i + self.__vector_range_value[0][1]])
                for i in range(self.__variable_number - self.__vector_range_value[0][1])
            ), GRB.MINIMIZE)

            self.__m.addConstrs((
                ancestor[i] - ancestor[self.__vector_symmetry_value[i]] == 0
                for i in range(self.__variable_number)), name='symmetry'
            )

            self.__m.addConstrs((
                gp.quicksum(ancestor[j] for j in self.__candidate_adjacency[i + 1]) == 1
                for i in range(self.__dim - 1)), name='dup_minimun'
            )

            self.__m.addConstrs((
                gp.quicksum(ancestor[j + self.__vector_range_value[i + 1][0]]
                            for j in range(self.__vector_range_value[i + 1][1] - self.__vector_range_value[i + 1][0])) == 1
                for i in range(self.__dim - 1)), name='row_unique'
            )
            self.__m.addConstrs(
                (ancestor[i] == 0 for i in self.__diagonal_value),
                name='Diagonal'
            )
            self.__m.optimize()
            print('Obj: %g' % self.__m.objVal)
        except gp.GurobiError as e:
            print('Error code ' + str(e.errno) + ': ' + str(e))

        except AttributeError:
            print('Encountered an attribute error')

    def ancestor_adjacency_matrix(self):
        result = []
        for v in self.__m.getVars():
            result.append(v.x)
        self.__adjacency_matrix = {}
        for i in self.__matrix_items:
            self.__adjacency_matrix[i] = {}
            for j in self.__matrix_items:
                self.__adjacency_matrix[i][j] = 0
        for i in range(len(self.__adjacency_name)):
            if result[i] != 0:
                keys = self.__adjacency_name[i].split('@')
                self.__adjacency_matrix[keys[0]][keys[1]] = result[i]
        self.__adjacency_matrix = pd.DataFrame(self.__adjacency_matrix)
        adjacency_matrix = AdjacencyMatrix()
        adjacency_matrix.readFromMatrix(self.__adjacency_matrix)
        return adjacency_matrix


# def main():
#     dir = 'D:/InferAncestorGenome/' \
#           'simulatiedData/simpleSimulation/' \
#           'MultiGGHP/divergence_change/threeStrategyTest2/'
#     sp1 = 'species.sequence.7'
#     sp2 = 'species.sequence.9'
#     outdir = 'D:/InferAncestorGenome/simulatiedData/' \
#              'simpleSimulation/MultiGGHP/' \
#              'divergence_change/threeStrategyTest2/1/'
#     adj_file = outdir + '9_7.adj'
#     filelist = [outdir + sp2 + '.matching', outdir + sp1 + '.matching']
#     transformToAdjacency(filelist, adj_file)
#     output_sequence_file = outdir + '8.sequence'
#     output_matrix_file = outdir + '8.matrix.xls'
#     species_list = ['9', '7']
#     ggap = GGAP(adj_file, species_list, duptype=2)
#     ggap.optimization()
#     adjacency_matrix = ggap.ancestor_adjacency_matrix()
#     evaluation_dir = outdir
#     ggap.evaluation_ancestor_adjacencies(evaluation_dir)
#     adjacency_matrix.output(output_matrix_file)
#     adjacency_matrix.assemble()
#     adjacency_matrix.out_assembly(output_sequence_file)

#
# if __name__ == '__main__':
#     main()