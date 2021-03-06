import os
import regex as re
import copy
import math
import pandas as pd

from kn_tools.basic_tools import text_to_df
from kn_tools.rna_path_tools import get_all_paths

def go_to_bed(rna,raw):
    """Makes the BED Files for all possible paths and outputs to bed and bedinfo directories, will consider existing directories and enumerate, also labels with RNA Name

    Parameters:
        rna: (str) The name of the RNA in the name2 column
        raw: (str) File location of UCSC data

    Returns:
        fcount: (int) The ID of the directory
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
            start = str(int(exon[0]))
            end = str(int(exon[1]))
            info = [chrom,start,end,name,'score',strand]
            #if strand == '+':
            ostrich += '\t'.join(info) + '\n'
            #else:
                #ostrich = '\t'.join(info) + '\n' + ostrich
            
        fname = rna + '_bed' + '_' + str(ncount) + '.txt'
        ncount += 1
        
        fname = rnadir + fname
        f = open(fname,'w')
        f.write(ostrich)
        f.close()
        ostrich = ''

    return fcount

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
            elif nodes[j] <= ind:
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


def gene_index(ind,nodes,strand):
    """Convert ind from sequence to genomic coordinates based on nodes and strand
    
    Parameters:
        ind: (int) The sequence coordinate, usually a long number from UCSC browser
        nodes: (list of int) The genomic coordinates of the exons start/ends
        strand: (str) + or - depending on the transcript

    Returns:
        to_return: (int) The genomic coordinate
    """
    running = int(ind)
    nodes = map(int,nodes)
    
    if strand == '+':
        for i in range(0,len(nodes)-1,2):
            if nodes[i+1] - nodes[i] < running:
                running -= nodes[i+1] - nodes[i]
            else:
                to_return = running + nodes[i]
                return to_return
    else:
        for i in range(len(nodes)-1,0,-2):
            if nodes[i] - nodes[i-1] < running:
                running -= nodes[i] - nodes[i-1]
            else:
                to_return = nodes[i] - running
                return to_return



def fetch_coords(seq):
    """Gets the location of the start and stop codon

    Parameters:
        seq: (str) The string of A,T,G,C from the bed file

    Returns:
        se: (list of 2 int) first value is the start, second value is the stop
    """
    matches = list(re.finditer('(ATG)(...)*?(TAG|TAA|TGA).*?',seq,overlapped=True))
    
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

def set_flags(ss_df):
    """ Calculates number of bases to the last junction and exons from end in the startstop DataFrame
    
    Parameters:
        ss_df: (Pandas DataFrame) The outputted DataFrame from startstop.py with the added columns 'exonsFromEnd', 'basesFromJunct', 'exists', 'theorized_nmd', 'erroneous_nmd' if in last exon, 'basesFromJunct' is -1, 'exists', 'theorized_nmd', 'erroneous_nmd' are 0 and 1 boolean value
    
    Returns:
        new_df: (Pandas DataFrame) ss_df updated with the aforementioned data
    """
    new_df = copy.deepcopy(ss_df)
    for index,row in new_df.iterrows():
        seq_nodes = sorted(map(int,str(row.get_value('seqNodes')).split(',')))
        
        end = int(row.get_value('regexStop'))
        i1 = len(seq_nodes) - 1
        i2 = 0
        s = -2
        
        junct = i1+s/2
        dist = ''
        excount = 0

        for i in range(i1,i2,s):
            ex_test = sorted([seq_nodes[i],seq_nodes[i + s/2]])
            
            # If the sequence coordinates of the end lie in the sequence exon of the ends exon
            if ex_test[0] < end and end <= ex_test[1]:
                dist = abs(seq_nodes[junct] - end)
                if excount == 0:
                    dist = -1
                break
            excount += 1
        
        new_df.loc[index,'exonsFromEnd'] = excount
        new_df.loc[index,'basesFromJunct'] = dist

        if row.get_value('basesFromJunct') >= 50:
            new_df.loc[index,'theorized_nmd'] = 1
        new_df.loc[index,'erroneous_nmd'] = (row.get_value('exists') + row.get_value('theorized_nmd') + 1) % 2
    return new_df

def entropy(x):
    """Gets the entropy for a list of boolean

    Parameters:
        x: (list of int) List of boolean, 0 false and 1 true used for inclusion and exclusion of exons

    Returns:
        H_x: (float) Entropy calculated
    """
    p = x.count(1) / float(len(x))
    q = x.count(0) / float(len(x))
    if p == 0 or q == 0:
        return 0
    H_x = -p*math.log(p,2) -q*math.log(q,2)
    return H_x

def mutual_inf(pairs):
    """Gets the mutual information for a list of boolean pairs, uses I(X,Y) = H(X) + H(Y) - H(X,Y)

    Parameters:
        pairs: (list of int) List of boolean pairs, 0 false and 1 true used for inclusion and exclusion of exons

    Returns:
        MI: (float) Mutual information calculated
        detailed: (dictionary) Provides counts for each pairing (0,0), (1,1), (0,1), (1,0)
    """
    x = [item[0] for item in pairs]
    y = [item[1] for item in pairs]
    t = len(x)
    total = 0
    detailed = {}
    H_x = entropy(x)
    H_y = entropy(y)
    H_xy = 0
    for i in range(2):
        for j in range(2):
            p_xy = float(pairs.count((i,j)))/t
            c_xy = pairs.count((i,j))
            detailed[(i,j)] = c_xy
            if p_xy == 0:
                H_xy += 0
            else:
                H_xy += -p_xy*math.log(p_xy,2)
    MI = H_x + H_y - H_xy
    return MI,detailed

def get_splice_sites(mod_df):
    """Gets the splice site information about a mod dataframe
    
    Parameters:
        mod_df: (Pandas DataFrame) For one rna, the mod_df, must have strand and name2 columns
    
    Returns:
        ss_df: (Pandas DataFrame) Includes the unique splice site info for the rna
    """
    ss_list = []
    ss_set = []
    exon_cols = list(mod_df.columns)[5:-1]
    for index,row in mod_df.iterrows():
        mods = [row.get_value(col) for col in exon_cols]
        nmd_ind = row.get_value('NMD_ind')
        strand = row.get_value('strand')
        for i in range(len(mods)-1,-1,-1):
            if mods[i] != -1:
                if strand == '+':
                    coord = row.get_value('cdsEnd')
                else:
                    coord = row.get_value('cdsStart')
                if not (coord,exon_cols[i],mods[i],nmd_ind) in ss_set:
                    ss_list.append([coord,exon_cols[i],mods[i],nmd_ind])
                    ss_set.append((coord,exon_cols[i],mods[i],nmd_ind))
                break
    ss_df = pd.DataFrame(ss_list,columns=['splice_site','pexon','frame','NMD_ind'])
    ss_df.insert(0,'name2',row.get_value('name2'))
    ss_df.insert(1,'strand',row.get_value('strand'))
    return ss_df
        