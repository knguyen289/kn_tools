import pandas as pd
import copy
import numpy as np 
from scipy.spatial.distance import euclidean

def filter_df(df,indexes):
    """Returns a dataframe that includes info from indicated indexes
    
    Parameters:
        df: (Pandas DataFrame) The original DataFrame
        indexes: (list of str) The list of indexes
    """
    head = df.columns
    data = []
    for i in indexes:
        row = []
        for c in head:
            row.append(df.get_value(i,c))
        data.append(row)
    new_df = pd.DataFrame(data,columns=head)
    new_df.index = indexes
    return new_df

def redact_df(df,columns):
    """Returns a dataframe that includes info from indicated columns
    
    Parameters:
        df: (Pandas DataFrame) The original DataFrame
        indexes: (list of str) The list of column names
    """
    ind = df.index
    data = []
    for c in columns:
        vert = list(df[c])
        data.append(vert)
    new_df = pd.DataFrame(data)
    new_df = new_df.transpose()
    new_df.columns = columns
    new_df.index = ind
    return new_df

def get_certain_mef(df,mef):
    """Get certain mef or columns from a df that includes cytosolic, membrane, insoluble data

    Parameters:
        df: (Pandas DataFrame) The original DataFrame
        mef: (str) The name of the mef, prefix before _cyt, _mem, _ins"""
    types = ['cyt','mem','ins']
    cs = [mef + '_' + t for t in types]
    return redact_df(df,cs)

def exp_err_adj(c,m,i,p_c,p_m):
    """Experimental error adjustment
    
    Parameters:
        c = (float) Cytosolic experimental value
        m = (float) Membrane experimental value
        i = (float) Insoluble experimental value
        p_c = (float) Proportion of actual cytosolic obtained, between 0 and 1
        p_m = (float) Proportion of actual membrane obtained, between 0 and 1"""
    adj_c = c/p_c
    adj_m = m/p_m + (p_c-1)*adj_c
    adj_i = c + m + i - adj_c - adj_m

    return [adj_c,adj_m,adj_i]

def adj_df(df,mef,p_c,p_m):
    """Adjusts the df for experimental extraction error
    
    Parameters:
        df: (Pandas DataFrame) The experimentally derived DataFrame
        mef: (str) The MEF or prefix of columns [with suffix _ins and the like]
        p_c = (float) Proportion of actual cytosolic obtained, between 0 and 1
        p_m = (float) Proportion of actual membrane obtained, between 0 and 1"""
    return_df = copy.deepcopy(df)
    types = ['cyt','mem','ins',]
    cols = [mef + '_' + t for t in types]

    for index,row in df.iterrows():
        c = float(row.get_value(cols[0]))
        m = float(row.get_value(cols[1]))
        i = float(row.get_value(cols[2]))
        c,m,i = exp_err_adj(c,m,i,p_c,p_m)

        return_df.loc[index,cols[0]] = c
        return_df.loc[index,cols[1]] = m
        return_df.loc[index,cols[2]] = i
        
    return return_df
    
def text_to_df(filename,index=None):
    """Converts a text file that is tab delimited with first row being header to Pandas DataFrame

    Parameters:
        filename: (str) The location of the text file
        index: (str) The column that can be used as an index [Optional]
    """
    f = open(filename,'r')
    data = []
    head = []
    for line in f:
        line = line.rstrip()
        if len(head) == 0:
            head = line.split('\t')
            
        else:
            data.append(line.split('\t'))

    df = pd.DataFrame(data,columns=head)
    if index != None:
        df.index = df[index]
        df.drop(index,1,inplace=True)

    return df

def trieuclid(df1,df2,d3):
    """Gets a list of the perimeters of the triangle created by the gene locations of each gene in each of df1,2,3...Returns the list of distances and the ordered list of genes...MUST HAVE THE SAME MEF NAME

    Parameters:
        df1: (Pandas Dataframe) The first MEF DataFrame
        df2: (Pandas Dataframe) The second MEF DataFrame
        df3: (Pandas Dataframe) The third MEF DataFrame"""
    genes = list(df1.index)
    mef = df1.columns[0].split('_')[0] + '_' + df1.columns[0].split('_')[1]
    c = mef + '_cyt'
    m = mef + '_mem'
    i = mef + '_ins'
    
    incl_dist = []
    incl_genes = []
    for g in genes:
        vecs = []
        for df in [df1,df2,df3]:
            temp = [df.loc[g,c],df.loc[g,m],df.loc[g,i]]
            temp = [float(item) for item in temp]
            df.append(temp)

        if sum(vecs[0]) != 0 and sum(vecs[1]) != 0 and sum(vecs[2]) != 0:
            v1 = vecs[0]
            v2 = vecs[1]
            v3 = vecs[2]
            incl_dist.append(euclidean(v1,v2)+euclidean(v1,v3)+euclidean(v2,v3))
            incl_genes.append(g)
    
    return incl_dist,incl_genes
