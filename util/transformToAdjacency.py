def transformToAdjacency(file_list,outfile):
    outfile = open(outfile, 'w')
    for i in file_list:
        adjacency_list = []
        with open(i) as df:
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
        for j in adjacency_list:
            outfile.write(j[0] + ' ' + j[1] + ' ')
        outfile.write('\n')
    outfile.close()

