# 1、构建结构
from inferringAncestorGenomeStructure.MatchingOptimization import MatchingOptimization
from inferringAncestorGenomeStructure.GGAP import GGAP
from util.transformToAdjacency import transformToAdjacency
from inferringAncestorGenomeStructure.GMP import GMP
from inferringGeneContent.GeneContentBuilder import *


in_workdir = 'D:/InferAncestorGenome/DCW_HN1_YMR_newpair/IAG/inputFiles/'
out_workdir = 'D:/InferAncestorGenome/DCW_HN1_YMR_newpair/IAG/outputFiles/'
# 1. DGGHP
file = ['PT.filter.block','PS.filter.block']
outcandidatefile = out_workdir + '1_PT_pre/PT.filter.block.matching'
outguidedfile = out_workdir + '1_PT_pre/PS.filter.block.matching'
outmatching = out_workdir + '1_PT_pre/matching_PT_PS.txt'
filelist = []
for i in file:
    filelist.append(in_workdir + i)
mo = MatchingOptimization(filelist,matching_dim1=4,matching_dim2=2,
                              relation1=1,relation2=2)
mo.optimization()
mo.matching_relation()
mo.output_matching_relation(outmatching)
mo.output_new_sequence(outcandidatefile,outguidedfile)

adj_file = out_workdir + '1_PT_pre/PT_PS.adj'
filelist = [outcandidatefile,outguidedfile]
transformToAdjacency(filelist,adj_file)
output_sequence_file = out_workdir + '1_PT_pre/PT_preWGD_ancestor.block'

output_matrix_file = out_workdir + '1_PT_pre/PT_preWGD_ancestor.matrix.xls'
species_list = ['PT', 'PS']
ggap = GGAP(adj_file,species_list,duptype=2)
ggap.optimization()
adjacency_matrix = ggap.ancestor_adjacency_matrix()
ggap.evaluation_ancestor_adjacencies(out_workdir + '1_PT_pre/')
adjacency_matrix.output(output_matrix_file)
adjacency_matrix.assemble()
adjacency_matrix.out_assembly(output_sequence_file,remove_bar=False)

# 2. DGMP
def doublePR(infile,outfile):
    outfile = open(outfile,'w')
    sequence = []
    with open(infile,'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            sequence.append(line)
    for i in sequence:
        outfile.write(i)
    for i in sequence:
        outfile.write(i)

doublePR(in_workdir + 'PR.filter.block',out_workdir + '/2_PS_post/PR.filter.block.double')

outcandidatefile = out_workdir + '/2_PS_post/PR.filter.block.double.matching'
outguidedfile = out_workdir + '/2_PS_post/PS.filter.block.matching'
outmatching = out_workdir + '/2_PS_post/matching_PRdouble_PS.txt'

filelist = [out_workdir + '/2_PS_post/PR.filter.block.double',in_workdir + 'PS.filter.block']

mo = MatchingOptimization(filelist,matching_dim1=2,matching_dim2=2,
                              relation1=1,relation2=1)
mo.optimization()
mo.matching_relation()
mo.output_matching_relation(outmatching)
mo.output_new_sequence(outcandidatefile,outguidedfile)

adj_file =  out_workdir + '/2_PS_post/prePT_PS_PRdouble.adj'
filelist = [out_workdir + '1_PT_pre/PT_preWGD_ancestor.block',
            out_workdir + '/2_PS_post/PS.filter.block.matching',
            out_workdir + '/2_PS_post/PR.filter.block.double.matching']

transformToAdjacency(filelist, adj_file)
output_sequence_file = out_workdir + '/2_PS_post/PS_postWGD_ancestor.block'
output_matrix_file = out_workdir + '/2_PS_post/PS_postWGD_ancestor.matrix.xls'
species_list = ['PT_preWGD_ancestor', 'PS','PR']
gmp = GMP(adj_file,species_list)
gmp.optimization()
adjacency_matrix = gmp.ancestor_adjacency_matrix()
gmp.evaluation_ancestor_adjacencies(out_workdir + '/2_PS_post/')

adjacency_matrix.assemble()
output_sequence_file_with_bar = out_workdir + '/2_PS_post/PS_postWGD_ancestor.block.withbar'
adjacency_matrix.out_assembly(output_sequence_file_with_bar,remove_bar=False)
adjacency_matrix.out_assembly(output_sequence_file,remove_bar=True)

# 3. GGHP
adj_file = out_workdir + '/3_PS_pre/PS_postWGD_PR.adj'
filelist = [out_workdir + '/2_PS_post/PS_postWGD_ancestor.block',
            in_workdir + 'PR.filter.block']

transformToAdjacency(filelist,adj_file)
output_sequence_file =  out_workdir + '/3_PS_pre/PS_preWGD_ancestor.block'

species_list = ['PS_postWGD_ancestor','PR']
ggap = GGAP(adj_file,species_list,duptype=2)
ggap.optimization()
adjacency_matrix = ggap.ancestor_adjacency_matrix()
ggap.evaluation_ancestor_adjacencies(out_workdir + '/3_PS_pre/')
adjacency_matrix.output(output_matrix_file)
adjacency_matrix.assemble()
adjacency_matrix.out_assembly(output_sequence_file,remove_bar=False)

# 2、构建gene content
def builderGeneContent():
    multi_synteny_file = in_workdir + 'multiblock_synteny.txt'
    multi_synteny_name_file = in_workdir + 'multiblock_synteny.name.txt'
    block_sequence = {}
    with open(multi_synteny_file,'r') as msf:
        while True:
            line = msf.readline()[:-2]
            if not line:
                break
            itemset = line.split(' ')
            header = itemset[0].split(':')
            genes = itemset[1:]
            if header[0] not in block_sequence.keys():
                block_sequence[header[0]] = {}
            block_sequence[header[0]][int(header[1])] = [genes]
    with open(multi_synteny_name_file,'r') as msf:
        while True:
            line = msf.readline()[:-2]
            if not line:
                break
            itemset = line.split(' ')
            header = itemset[0].split(':')
            genenames = itemset[1:]
            block_sequence[header[0]][int(header[1])].append(genenames)

    dcw_block_file = out_workdir + '1_PT_pre/PT.filter.block.matching'
    hn1_block_file = out_workdir + '1_PT_pre/PS.filter.block.matching'
    dcw_block = read_block(dcw_block_file)
    hn1_block = read_block(hn1_block_file)

    # 按照匹配分堆，将尾号一致的分到一个集合
    dan_copy1_sequence = {}
    dan_copy2_sequence = {}
    for i in dcw_block.keys():
        dan_copy1_sequence[i] = {}
        dan_copy2_sequence[i] = {}
        for j in dcw_block[i]:
            if j[1] == '1':
                dan_copy1_sequence[i][j[0]] = block_sequence[i][j[0]]
            else:
                dan_copy2_sequence[i][j[0]] = block_sequence[i][j[0]]

    for i in hn1_block.keys():
        for j in hn1_block[i]:
            if j[1] == '1':
                dan_copy1_sequence[i][j[0] + 4] = block_sequence[i][j[0] + 4]
            else:
                dan_copy2_sequence[i][j[0] + 4] = block_sequence[i][j[0] + 4]

    species = {1: 'PT', 2: 'PT', 3: 'PS'}
    homo_threshold = 2
    dan_copy1_gene_content = GeneContentBuilder(dan_copy1_sequence,species,homo_threshold)
    dan_copy1_synteny_file = out_workdir + '4_block_gene_content/dan_copy1_synteny_ps.txt'
    dan_copy1_synteny_name_file = out_workdir + '4_block_gene_content/dan_copy1_synteny.name_ps.txt'
    dan_copy1_block_gene_name = out_workdir + '4_block_gene_content/dan_copy1_synteny.genename_ps.txt'
    support = 3
    block_suffix = '_1'
    dan_copy1_gene_content.out_synteny_file(dan_copy1_synteny_file,dan_copy1_synteny_name_file,support,block_suffix)
    dan_copy1_gene_content.out_block_dictionary(dan_copy1_block_gene_name)
    dan_copy2_gene_content = GeneContentBuilder(dan_copy2_sequence,species,homo_threshold)
    dan_copy2_synteny_file = out_workdir + '4_block_gene_content/dan_copy2_synteny_ps.txt'
    dan_copy2_synteny_name_file = out_workdir + '4_block_gene_content/dan_copy2_synteny.name_ps.txt'
    dan_copy2_block_gene_name = out_workdir + '4_block_gene_content/dan_copy2_synteny.genename_ps.txt'
    support = 3
    block_suffix = '_2'
    dan_copy2_gene_content.out_synteny_file(dan_copy2_synteny_file, dan_copy2_synteny_name_file, support, block_suffix)
    dan_copy2_gene_content.out_block_dictionary(dan_copy2_block_gene_name)
    print('---------------han-------------')
    # 6个重建一个
    han_sequence = {}
    for i in block_sequence.keys():
        han_sequence[i] = {}
        for j in range(6):
            han_sequence[i][j + 1] = block_sequence[i][j + 1]
    # 保守要求在两个物种中存在
    species = {1: 'PT', 2: 'PT', 3: 'PT', 4: 'PT', 5: 'PS', 6: 'PS'}
    homo_threshold = 2
    han_gene_content = GeneContentBuilder(han_sequence, species, homo_threshold)
    han_synteny_file = out_workdir + '4_block_gene_content/han_synteny_ps.txt'
    han_synteny_name_file = out_workdir + '4_block_gene_content/han_synteny.name_ps.txt'
    han_synteny_gene_name_file = out_workdir + '4_block_gene_content/han_synteny.genename_ps.txt'
    support = 6
    block_suffix = ''
    han_gene_content.out_synteny_file(han_synteny_file, han_synteny_name_file, support, block_suffix)
    han_gene_content.out_block_dictionary(han_synteny_gene_name_file)

builderGeneContent()