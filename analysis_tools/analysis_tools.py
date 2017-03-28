import os

from kn_tools.basic_tools import *
from kn_tools.rna_path_tools import *

def go_to_bed(rna,data_df):
    """Makes the BED Files for all possible paths and outputs to bed and bedinfo directories, will consider existing directories and enumerate, also labels with RNA Name

    Parameters:
        rna: (str) The name of the RNA in the name2 column
        data_df: (Pandas DataFrame) The UCSC dataframe from text_to_df
    """
    bed = kn.get_all_paths(rna,data_df,detail=True)
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
    
    for path in bed:
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
            
        else:
            pathdet = path[4:]
            
        for exon in pathdet:
            start = str(exon[0])
            end = str(exon[1])
            info = [chrom,start,end,name,'score',strand]
            ostrich += '\t'.join(info) + '\n'
            
        fname = rna + '_bed' + '_' + str(ncount) + '.txt'
        ncount += 1
        
        fname = rnadir + fname
        f = open(fname,'w')
        f.write(ostrich)
        f.close()
        