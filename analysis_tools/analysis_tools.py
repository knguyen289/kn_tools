import os
import regex as re

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
                return 0
            
    else:
        for j in range(len(nodes)-1,0,2):
            if nodes[j-1] > ind:
                to_return += nodes[j] - nodes[j-1]
            elif nodes[j] > ind:
                to_return += nodes[j] - ind
                return to_return
            else:
                print 'Error: Invalid exon nodes'
                return 0

def regex_possible(bases):
    """Gets the location of all possible start and end codons

    Parameters:
        bases: (str) The string of A,T,G,C from the bed file

    Returns:
        atg_inds: (list of int) The possible cdsStarts in sequence coordinates
        stop_inds: (list of int) The possible cdsEnds in sequence coordinates
    """
    starts = list(re.finditer('(...)+?(ATG)(...)+?',bases,overlapped=True))
    stops = list(re.finditer('(...)+?(TAG|TAA|TGA)(...)+?',bases,overlapped=True))
    
    atg_inds = []
    for item in starts:
        atg_inds.append(item.span(2))
    atg_inds = list(set(atg_inds))

    stop_inds = []
    for item in stops:
        stop_inds.append(item.span(2))
    stop_inds = list(set(stop_inds))
    
    atg_inds = sorted([item[0] for item in atg_inds])
    stop_inds = sorted([item[1] for item in stop_inds])
    return atg_inds,stop_inds
        