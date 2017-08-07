import pandas as pd
import numpy as np
import networkx as nx
import copy as copy

def get_rna_dfs(rna,data_df):
	"""Filters the UCSC dataframe for a specified RNA in the name2 column, gives a list of one or two dataframes, drops non coding rna
	
	Parameters:
		rna: (str) The name of the RNA in the name2 column
		data_df: (Pandas DataFrame) The UCSC dataframe from text_to_df
	
	Returns:
		final_df: (Pandas DataFrame) Filtered DataFrame, returns None if there is more than one strand defined
	"""
	temp_df = data_df[data_df['name2'] == rna]

	strands = list(set(temp_df['strand']))
	if len(strands) >= 2:
		return None

	non_coding = []
	for index,row in temp_df:
		if int(temp_df.loc[index,'cdsStart']) == int(temp_df[index,'cdsEnd']):
			non_coding += index
	final_df = temp_df.drop(non_coding)

	return final_df

def get_lookup1(df):
	"""Gets the first lookup table for the given RNA DataFrame produced by get_rna_dfs. UPDATE 6/22/17 - Accomodated for minus strands

	Parameters:
		df: (Pandas DataFrame) The RNA DataFrame produced by get_rna_dfs
	
	Returns:
		lookup_df: (Pandas DataFrame) The first lookup table for exon ID, has upper and lower bounds for each exon, all exons, and number of diff forms per exon
	"""
	#get all exons
	master_exons = []
	strand = list(set(df['strand']))[0]
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
		
		if strand == '-':
			adj_ex_ind = len(ex_inds) - i + 1
		else:
			 adj_ex_ind = i

		to_df.append([adj_ex_ind,s,e,','.join(map(str,starts)),','.join(map(str,ends)),len(ends)])
	
	#create the dataframe
	lookup_df = pd.DataFrame(to_df,columns=['exon','lb','ub','starts','ends','#']).set_index('exon')
	lookup_df.insert(len(lookup_df.columns),'strand',strand)
	return lookup_df

def get_lookup2(df):
	"""Gets the second lookup table for the given RNA DataFrame produced by get_rna_dfs. UPDATE 6/22/17 - Accomodated for minus strands
	
	Parameters:
		df: (Pandas DataFrame) The RNA DataFrame produced by get_rna_dfs
		
	Returns:
		lu2: (Pandas DataFrame) The second lookup table, masks the multiple splice sites as different exons
	"""
	lu1_df = get_lookup1(df)
	info = []

	for index,row in lu1_df.iterrows():
		strand = row.get_value('strand')
		if lu1_df.loc[index,'#'] == 1:
			se = map(int,[row.get_value('starts'),row.get_value('ends')])
			info.append([str(index)] + se + [(se[1] - se[0]) % 3])
		else:
			starts = map(int,row.get_value('starts').split(','))
			ends = map(int,row.get_value('ends').split(','))
			
			if strand == '-':
				uniq = sorted(list(set(starts+ends)),reverse=True)
			else:
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
				temp.append((exons[i][1]-exons[i][0]) % 3)
				info.append(temp)
				
	lu2 = pd.DataFrame(info,columns=['exon','start','end','mod'])
	if strand == '-':
		lu2.index = range(len(lu2),0,-1)
	else:
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
		detailed: (list) includes the real cdsStart and cdsEnd for a strand
		names: (list) includes the strand name
		indexes: (list) list of the dataframe index
	"""
	if len(lu_df) == 0:
		lu_df = get_lookup2(df)
	
	detailed = []
	names = []
	strand_nodes = []
	indexes = []
	for index,row in df.iterrows():
		s = sorted(map(int,row.get_value('exonStarts')[:-1].split(',')))
		e = sorted(map(int,row.get_value('exonEnds')[:-1].split(',')))
		exons = [[s[i],e[i]] for i in range(len(s))]
		
		pnodes = []
		nodes = []
		lu_exons = [[lu_df.loc[j][1],lu_df.loc[j][2]] for j in range(1,len(lu_df)+1)]

		if row.get_value('strand') == '+':
			order = range(0,len(exons))
		else:
			order = range(len(exons)-1,-1,-1)
		for i in order:
			exon = exons[i]

			if exon in lu_exons:
				pex_ind = lu_exons.index(exon)+1
				pnodes.append(pex_ind)
			else:
				raise


		detailed.append(map(int,[df.loc[index,'cdsStart'],df.loc[index,'cdsEnd']]))
		names.append(df.loc[index,'name'])
		strand_nodes.append(pnodes)
		indexes.append(index)
		
	return strand_nodes,detailed,names,indexes

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

def get_all_paths(rna,data_df,detail=False):
	"""Gets all paths with replacement of the pseudoexons, adds an '*' if it exists
	
	Parameters:
		rna: (str) The name of the RNA in the name2 column
		data_df: (Pandas DataFrame) The UCSC dataframe from text_to_df
		detail: (boolean) False if output is nodes, True for start and ends [Optional]
	
	Returns:
		paths: (list) List of all paths produced, None if a strange rna
	"""
	df = get_rna_dfs(rna,data_df)
	if df is None:
		return None
	paths = []
	lu_df = get_lookup2(df)
	print lu_df
	pnodes,pdetailed,names,indexes = get_pexons(df,lu_df)
	index = list(df.index)[0]
	strand = df.loc[index,'strand']
	chrom = df.loc[index,'chrom']
	print pnodes
	G,se_pairs = get_graph(pnodes)
	
	for se in se_pairs:
		for path in nx.all_simple_paths(G, source=se[0], target=se[1]):
			path = list(path)
			print path
			if path in pnodes:
				if detail:
					ind = pnodes.index(path)
					irl = pdetailed[ind]
					realname = names[ind]
					paths.append([rna,strand,chrom,realname]+[[lu_df.loc[item,'start'],lu_df.loc[item,'end']] for item in path]+irl)

				else:
					paths.append([rna,strand,chrom,realname]+[lu_df.loc[item,'exon'] for item in path])
			else:
				if detail:
					paths.append([rna,strand,chrom,' ']+[[lu_df.loc[item,'start'],lu_df.loc[item,'end']] for item in path])
				else:
					paths.append([rna,strand,chrom,' ']+map(str,path))
				
	return paths

