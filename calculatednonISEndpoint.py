from inferringAncestorGenomeStructure.MatchingOptimization import MatchingOptimization
import os
def transformToAdjacency_for_distance(file):
        adjacency_list = []
        with open(file) as df:
            while True:
                line = df.readline()[:-2]
                if not line:
                    break
                item = line.split(' ')
                chr_type = item[0]
                last = ''
                start = ''
                sequence = item[1:]
                for j in range(len(sequence)):
                    if j == 0:
                        if chr_type == 's':
                            if sequence[j].startswith('-'):
                                adjacency_list.append(['$', sequence[j][1:] + 'b'])
                                last = sequence[j][1:] + 'a'
                            else:
                                adjacency_list.append(['$', sequence[j] + 'a'])
                                last = sequence[j] + 'b'
                        else:
                            if sequence[j].startswith('-'):
                                last = sequence[j][1:] + 'a'
                                start = sequence[j][1:] + 'b'
                            else:
                                last = sequence[j] + 'b'
                                start = sequence[j] + 'a'
                    else:
                        if sequence[j].startswith('-'):
                            adjacency_list.append([last, sequence[j][1:] + 'b'])
                            last = sequence[j][1:] + 'a'
                        else:
                            adjacency_list.append([last, sequence[j] + 'a'])
                            last = sequence[j] + 'b'
                if chr_type == 's':
                    adjacency_list.append([last, '$'])
                else:
                    adjacency_list.append([last, start])

        new_adjacency_list = []
        for j in adjacency_list:
            new_adjacency_list.append(j[0]+'@'+j[1])
        return new_adjacency_list

def calculateBPReuseRateForPrePT(filelist,ploidy):
    mo = MatchingOptimization(filelist,matching_dim1=4, matching_dim2=2,
                                   relation1=1, relation2=2)
    mo.optimization()
    mo.matching_relation()
    mo.output_new_sequence(filelist[0]+'.relabelforbp',filelist[1]+'.relabelforbp')
    filelist = [filelist[0]+'.relabelforbp',
                filelist[1]+'.relabelforbp']

    endpoints = {}
    for i in filelist:
        adj = transformToAdjacency_for_distance(i)
        for j in adj:
            endpoint1 = j.split('@')[0]
            endpoint2 = j.split('@')[1]
            if endpoint1 not in endpoints.keys():
                endpoints[endpoint1] = [endpoint2]
            else:
                if endpoint2 not in endpoints[endpoint1]:
                    endpoints[endpoint1].append(endpoint2)

            if endpoint2 not in endpoints.keys():
                endpoints[endpoint2] = [endpoint1]
            else:
                if endpoint1 not in endpoints[endpoint2]:
                    endpoints[endpoint2].append(endpoint1)
    rate = 0

    for i in endpoints.keys():
        if i == '$':
            continue
        if len(endpoints[i]) == ploidy:
            rate += 1
    rate = round(rate / (len(endpoints.keys()) - 1), 4)
    os.remove(filelist[0])
    os.remove(filelist[1])
    return rate

def calculateBPReuseRateForPostPs(filelist,ploidy):
    mo = MatchingOptimization([filelist[1],filelist[0]],
                                   matching_dim1=2, matching_dim2=2,
                                   relation1=1, relation2=1)
    mo.optimization()
    mo.matching_relation()

    mo.output_new_sequence(filelist[1] + '.relabelforbp',filelist[0] + '.relabelforbp')

    mo = MatchingOptimization([filelist[2],filelist[0]],
                                   matching_dim1=2, matching_dim2=2,
                                   relation1=1, relation2=1)
    mo.optimization()
    mo.matching_relation()

    mo.output_new_sequence(filelist[2] + '.relabelforbp',filelist[0] + '.relabelforbp')
    filelist = [filelist[0] + '.relabelforbp',
                filelist[1] + '.relabelforbp',
                filelist[2] + '.relabelforbp']
    endpoints = {}
    for i in filelist:
        adj = transformToAdjacency_for_distance(i)
        for j in adj:
            endpoint1 = j.split('@')[0]
            endpoint2 = j.split('@')[1]
            if endpoint1 not in endpoints.keys():
                endpoints[endpoint1] = [endpoint2]
            else:
                if endpoint2 not in endpoints[endpoint1]:
                    endpoints[endpoint1].append(endpoint2)

            if endpoint2 not in endpoints.keys():
                endpoints[endpoint2] = [endpoint1]
            else:
                if endpoint1 not in endpoints[endpoint2]:
                    endpoints[endpoint2].append(endpoint1)
    rate = 0

    for i in endpoints.keys():
        if i == '$':
            continue
        if len(endpoints[i]) == ploidy:
            rate += 1
    rate = round(rate / (len(endpoints.keys()) - 1), 4)
    os.remove(filelist[0])
    os.remove(filelist[1])
    os.remove(filelist[2])
    return rate

def calculateBPReuseRateForPrePS(filelist,ploidy):
    endpoints = {}
    for i in filelist:
        adj = transformToAdjacency_for_distance(i)
        for j in adj:
            endpoint1 = j.split('@')[0]
            endpoint2 = j.split('@')[1]
            if endpoint1 not in endpoints.keys():
                endpoints[endpoint1] = [endpoint2]
            else:
                if endpoint2 not in endpoints[endpoint1]:
                    endpoints[endpoint1].append(endpoint2)

            if endpoint2 not in endpoints.keys():
                endpoints[endpoint2] = [endpoint1]
            else:
                if endpoint1 not in endpoints[endpoint2]:
                    endpoints[endpoint2].append(endpoint1)
    rate = 0
    for i in endpoints.keys():
        if i == '$':
            continue
        if len(endpoints[i]) == ploidy:
            rate += 1
    rate = round(rate / (len(endpoints.keys()) - 1), 4)
    return rate

blockdir = 'D:/InferAncestorGenome/' \
           'DCW_HN1_YMR_newpair/IAG/inputFiles/'
resultdir = 'D:/InferAncestorGenome/' \
            'DCW_HN1_YMR_newpair/IAG/outputFiles/'

filelist = [blockdir + 'PT.filter.block',blockdir + 'PS.filter.block']
prePTrate = calculateBPReuseRateForPrePT(filelist,3)

filelist = [blockdir + 'PS.filter.block',
            resultdir + '1_PT_pre/PT_preWGD_ancestor.block',
            resultdir + '2_PS_post/PR.filter.block.double']
postPSrate = calculateBPReuseRateForPostPs(filelist,3)

filelist = [resultdir + '2_PS_post/PS_postWGD_ancestor.block',blockdir + 'PR.filter.block']
prePSrate = calculateBPReuseRateForPrePS(filelist,3)
print('prePT:'+str(prePTrate))
print('postPS:'+str(postPSrate))
print('prePS:' + str(prePSrate))




