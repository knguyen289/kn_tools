import pandas as pd
import numpy as np
import networkx as nx
import copy as copy

def get_rna_dfs(rna,data_df):
	"""Filters the UCSC dataframe for a specified RNA in the name2 column, gives a list of one or two dataframes
	
	Parameters:
		rna: (str) The name of the RNA in the name2 column
		data_df: (Pandas DataFrame) The UCSC dataframe from text_to_df
	
	Returns:
		dfs: (list) A list of one or two filtered Pandas DataFrames, one for each strand
	"""
	temp_df = copy.deepcopy(data_df)
	temp_df = temp_df[temp_df['name2'] == rna]

	strands = list(set(temp_df['strand']))
	dfs = []
	for strand in strands:
		strand_df = temp_df[temp_df['strand'] == strand]
		dfs.append(strand_df)
	return dfs

def get_lookup1(df):
	"""Gets the first lookup table for the given RNA DataFrame produced by get_rna_dfs
	
	Parameters:
		df: (Pandas DataFrame) The RNA DataFrame produced by get_rna_dfs
	
	Returns:
		lookup_df: (Pandas DataFrame) The first lookup table for exon ID, has upper and lower bounds for each exon, all exons, and number of diff forms per exon
	"""
	#get all exons
	master_exons = []
	for index,row in df.iterrows():
		s = map(int,row.get_value('exonStarts')[:-1].split(','))
		e = map(int,row.get_value('exonEnds')[:-1].split(','))
		for i in range(len(s)):
			if not (s[i],e[i]) in master_exons:
				master_exons.append((s[i],e[i]))
	master_exons = np.array(master_exons,dtype = [('start',int),('end',int)])
	exons_arr = np.sort(master_exons,order=['start','end'])
	
	#check for overlap with exons
	ex_num = 1
	ex_data = []
	for i in range(len(exons_arr)):
		to_append = []
		if i != 0:
			if list(exons_arr[i-1]) + list(exons_arr[i]) == sorted(list(exons_arr[i-1]) + list(exons_arr[i])):
				ex_num += 1
		to_append.append(ex_num)
		to_append.append(exons_arr[i][0])
		to_append.append(exons_arr[i][1])
		ex_data.append(to_append)
		
	#if there is overlap, this will maintain that starts and ends are lower/upper bound
	to_df = []
	ex_inds = list(set([item[0] for item in ex_data]))
	for i in ex_inds:
		temp = [item for item in ex_data if item[0] == i]
		starts = [item[1] for item in temp]
		ends = [item[2] for item in temp]
		s = min(starts)
		e = max(ends)
		
		to_df.append([i,s,e,','.join(map(str,starts)),','.join(map(str,ends)),len(ends)])
	
	#create the dataframe
	lookup_df = pd.DataFrame(to_df,columns=['exon','lb','ub','starts','ends','#']).set_index('exon')
	return lookup_df

def get_lookup2(df):
	"""Gets the second lookup table for the given RNA DataFrame produced by get_rna_dfs
	
	Parameters:
		df: (Pandas DataFrame) The RNA DataFrame produced by get_rna_dfs
		
	Returns:
		lu2: (Pandas DataFrame) The second lookup table, masks the multiple splice sites as different exons
	"""
	lu1_df = get_lookup1(df)
	info = []
	
	for index,row in lu1_df.iterrows():
		if lu1_df.loc[index,'#'] == 1:
			info.append([str(index)] + map(int,[row.get_value('starts'),row.get_value('ends')]))
		else:
			starts = map(int,row.get_value('starts').split(','))
			ends = map(int,row.get_value('ends').split(','))
			uniq = sorted(list(set(starts+ends)))
			
			exons = [(starts[i],ends[i]) for i in range(len(starts))]
			
			exons = np.array(exons,dtype=[('s',int),('e',int)])
			exons = np.sort(exons,order=['s','e'])
			exons = list(exons)
			
			for i in range(len(exons)):
				temp = []
				s_ind = uniq.index(exons[i][0])+1
				e_ind = uniq.index(exons[i][1])+1
				
				temp.append(str(index)+'.'+str(s_ind)+'.'+str(e_ind))
				temp.append(exons[i][0])
				temp.append(exons[i][1])
				info.append(temp)
	lu2 = pd.DataFrame(info,columns=['exon','start','end'])
	lu2.index = np.arange(1, len(lu2) + 1)
	lu2.index.name = 'pseudoexon'
	return lu2

def get_pexons(df,lu_df=[]):
	"""Gets the psuedoexons for the given RNA DataFrame produced by get_rna_dfs
	
	Parameters:
		df: (Pandas DataFrame) The RNA DataFrame produced by get_rna_dfs
		lu_df: (Pandas DataFrame) If lookup table 2 has already been produced, you can input it [Optional]
	
	Returns:
		strand_nodes: (list) A list of the psuedoexons for conversion to a graph
	"""
	if len(lu_df) == 0:
		lu_df = get_lookup2(df)
	
	strand_nodes = []
	for index,row in df.iterrows():
		s = map(int,row.get_value('exonStarts')[:-1].split(','))
		e = map(int,row.get_value('exonEnds')[:-1].split(','))
		exons = [[s[i],e[i]] for i in range(len(s))]
		
		pnodes = []
		nodes = []
		
		for exon in exons:
			unset = True
			while unset:
				if len(pnodes) == 0:
					i1 = 1
				else:
					i1 = pnodes[-1]
				for i in range(i1,len(lu_df)+1):
					lu_ex = list(lu_df.loc[i])
					if lu_ex[1] == exon[0] and lu_ex[2] == exon[1]:
						pnodes.append(i)
						unset = False
						
		strand_nodes.append(pnodes)
		
	return strand_nodes

def get_graph(strand_nodes):
	"""Gets the graph for a given list of nodes
	
	Parameters:
		strand_nodes: (list) A list of the psuedoexons for conversion to a graph created by get_pexons
	
	Returns:
		G: (Networkx DiGraph) The digraph created by the list of nodes
		se_pairs: (list of lists) List of all possible start,end pairs in the graph
	"""
	G = nx.DiGraph()
	se_pairs = set([])
	for strand in strand_nodes:
		G.add_path(strand)
		se_pairs.add((min(strand),max(strand)))
	se_pairs = [list(item) for item in list(se_pairs)]
	return G,se_pairs

def get_all_paths(df):
    """Gets all paths with replacement of the pseudoexons, adds an '*' if it exists
    
    Parameters:
        df: (Pandas DataFrame) The RNA DataFrame produced by get_rna_dfs
    
    Returns:
        paths: (list) List of all paths produced
    """
    lu_df = get_lookup2(df)
    pnodes = get_pexons(df,lu_df)
    
    G,se_pairs = get_graph(pnodes)
    paths = []
    for se in se_pairs:
        for path in nx.all_simple_paths(G, source=se[0], target=se[1]):
            path = list(path)
            pathcount += 1
            if path in pnodes:
                paths.append(['*']+[lu_df.loc[item,'exon'] for item in path])
            else:
                paths.append([' ']+map(str,path))
                
    return paths

def print_all_paths(df):
	"""Prints all paths with replacement of the pseudoexons, adds an '*' if it exists
	
	Parameters:
		df: (Pandas DataFrame) The RNA DataFrame produced by get_rna_dfs
	"""
	lu_df = get_lookup2(df)
	pnodes = get_pexons(df,lu_df)
	
	G,se_pairs = get_graph(pnodes)
	pathcount = 0
	for se in se_pairs:
		for path in nx.all_simple_paths(G, source=se[0], target=se[1]):
			path = list(path)
			pathcount += 1
			if path in pnodes:
				path = ' | '.join([lu_df.loc[item,'exon'] for item in path])
				print '*',path
			else:
				path = '  '+' | '.join(map(str,path))
				print path
				
	print pathcount, 'Total Paths'

def print_rna_paths(rna,data_df):
	"""Prints all RNA paths for a given RNA in name2 of UCSC DataFrame
	
	Parameters:
		rna: (str) The name of the RNA in the name2 column
		data_df: (Pandas DataFrame) The UCSC dataframe from text_to_df
	"""
	dfs = get_rna_dfs(rna,data_df)
	if len(dfs) != 1:
		print 'Warning: Different strand directions may have exons with the same number that are not truly equivalent\n'
	else:
		print
	for df in dfs:
		strand = list(set(df['strand']))[0]
		pnodes = get_pexons(df)
		print 'Strand Direction:',strand
		if len(pnodes[0]) == 1:
			print 'No producable graph'
		else:
			print_all_paths(df)
		print
