import os
import re

from kn_tools.basic_tools import text_to_df
from kn_tools.rna_path_tools import get_all_paths

def go_to_bed(rna,raw):
    """Makes the BED Files for all possible paths and outputs to bed and bedinfo directories, will consider existing directories and enumerate, also labels with RNA Name

    Parameters:
        rna: (str) The name of the RNA in the name2 column
        raw: (str) File location of UCSC data
    """
    data_df = text_to_df(raw)
    paths = get_all_paths(rna,data_df,detail=True)
    ostrich = ''
    realinfo = ''
    ncount = 1
    
    fcount = 1
    while os.path.exists(rna + 'bed' + str(fcount) + '/'):
        fcount += 1
    os.mkdir(rna + 'bed' + str(fcount) + '/')
    os.mkdir(rna + 'bedinfo' + str(fcount) + '/')
    rnadir = rna + 'bed' + str(fcount) + '/'
    infodir = rna + 'bedinfo' + str(fcount) + '/'
    
    for path in paths:
        name = path[0]
        strand = path[1]
        chrom = path[2]
        
        if path[3] != ' ':
            pathdet = path[4:-2]
            realinfo +=  '\t'.join([path[3]] + map(str,path[-2:]))
            
            finfo = infodir + rna + '_info_' + str(ncount) + '.txt'
            f1 = open(finfo,'w')
            f1.write(realinfo)
            f1.close()
            realinfo = ''
            
        else:
            pathdet = path[4:]
            
        for exon in pathdet:
            start = str(exon[0])
            end = str(exon[1])
            info = [chrom,start,end,name,'score',strand]
            if strand == '+':
                ostrich += '\t'.join(info) + '\n'
            else:
                ostrich = '\t'.join(info) + '\n' + ostrich
            
        fname = rna + '_bed' + '_' + str(ncount) + '.txt'
        ncount += 1
        
        fname = rnadir + fname
        f = open(fname,'w')
        f.write(ostrich)
        f.close()
        ostrich = ''

def seq_index(ind,nodes,strand):
    """Convert ind from genomic to sequence coordinates based on nodes and strand
    
    Parameters:
        ind: (int) The genomic coordinate, usually a long number from UCSC browser
        nodes: (list of int) The genomic coordinates of the exons start/ends
        strand: (str) + or - depending on the transcript

    Returns:
        to_return: (int) The sequence coordinate
    """
    to_return = 0
    ind = int(ind)
    nodes = map(int,nodes)
    if strand == '+':
        for j in range(0,len(nodes)-1,2):
            if nodes[j+1] < ind:
                to_return += nodes[j+1] - nodes[j]
            elif nodes[j] < ind:
                to_return += ind - nodes[j]
                return to_return
            else:
                print 'Error: Invalid exon nodes'
                return None
            
    else:
        for j in range(len(nodes)-1,0,-2):
            if nodes[j-1] > ind:
                to_return += nodes[j] - nodes[j-1]
            elif nodes[j] >= ind:
                to_return += nodes[j] - ind
                return to_return
            else:
                print 'Error: Invalid exon nodes'
                return None


def fetch_coords(seq):
    """Gets the location of the start and stop codon

    Parameters:
        seq: (str) The string of A,T,G,C from the bed file

    Returns:
        se: (list of 2 int) first value is the start, second value is the stop
    """
    matches = list(re.finditer('(ATG)(...)*?(TAG|TAA|TGA).*?',seq))
    
    pairs = []
    
    starts = sorted(list(set([item.span(1)[0] for item in matches])))
    stops = sorted(list(set([item.span(3)[1] for item in matches])))
    
    temp_i = []
    temp_j = []
    temp_diffs = []
    
    for i in range(len(starts)):
        mod = starts[i] % 3
        
        j = 0
        while j >= 0 and j < len(stops):
            if stops[j] > starts[i]:
                if stops[j] % 3 == mod:
                    temp_i.append(i)
                    temp_j.append(j)
                    temp_diffs.append(stops[j]-starts[i])
                    j = len(stops) + 1
            j += 1
    
    ind = temp_diffs.index(max(temp_diffs))
    se = [starts[temp_i[ind]],stops[temp_j[ind]]]
    return se
        