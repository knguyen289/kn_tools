import pandas as pd
import copy

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
    types = ['cyt','mem','ins']
    cs = [mef + '_' + t for t in types]
    return redact_df(df,cs)

def text_to_df(filename,index=None):
    """Converts a text file that is tab delimited with first row being header to Pandas DataFrame

    Parameters:
        filename: (str) The location of the text file
        index: (str) The column that can be used as an index [Optional]

    Returns:
        df: (Pandas DataFrame) The converted DataFrame
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

