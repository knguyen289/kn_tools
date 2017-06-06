# kn_tools
Kim Nguyen's tools for Wang Lab, includes tools for plotting and data management
## basic_tools/basic_tools.py

### filter_df()
#### Description:
* Returns a dataframe that includes info from indicated indexes

#### Parameters:
* **df:** *(Pandas DataFrame)* The original DataFrame

* **indexes:** *(list of str)* The list of indexes


### redact_df()
#### Description:
* Returns a dataframe that includes info from indicated columns

#### Parameters:
* **df:** *(Pandas DataFrame)* The original DataFrame

* **indexes:** *(list of str)* The list of column names


### get_certain_mef()
#### Description:
* Get certain mef or columns from a df that includes cytosolic, membrane, insoluble data

#### Parameters:
* **df:** *(Pandas DataFrame)* The original DataFrame


### exp_err_adj()
#### Description:
* Experimental error adjustment

#### Parameters:
* **c =** *(float)* Cytosolic experimental value

* **m =** *(float)* Membrane experimental value

* **i =** *(float)* Insoluble experimental value

* **p_c =** *(float)* Proportion of actual cytosolic obtained, between 0 and 1


### adj_df()
#### Description:
* Adjusts the df for experimental extraction error

#### Parameters:
* **df:** *(Pandas DataFrame)* The experimentally derived DataFrame

* **mef:** *(str)* The MEF or prefix of columns [with suffix _ins and the like]

* **p_c =** *(float)* Proportion of actual cytosolic obtained, between 0 and 1


### text_to_df()
#### Description:
* Converts a text file that is delimited with first row being header to Pandas DataFrame

#### Parameters:
* **filename:** *(str)* The location of the text file

* **index:** *(str)* The column that can be used as an index [Optional]

* **sep:** *(str)* The separator used in the file

#### Returns:
* **df:** *(Pandas Dataframe)* The desired dataframe

### trieuclid()
#### Description:
* Gets a list of the perimeters of the triangle created by the gene locations of each gene in each of df1,2,3...Returns the list of distances and the ordered list of genes...MUST HAVE THE SAME MEF NAME

#### Parameters:
* **df1:** *(Pandas Dataframe)* The first MEF DataFrame

* **df2:** *(Pandas Dataframe)* The second MEF DataFrame

## ternary_tools/ternary_tools.py

### get_mag()
#### Description:
* Gets the magnitude of a vector

#### Parameters:
* **x_vals:** *(list)* [x_0,x_1]

* **y_vals:** *(list)* [y_0,y_1]


### get_angle()
#### Description:
* Gets the degree angle of a vector

#### Parameters:
* **x_vals:** *(list)* [x_0,x_1]

* **y_vals:** *(list)* [y_0,y_1]


### get_points()
#### Description:
* Gets the euclidean coordinates of the ternary plot points along with the ordered gene names

#### Parameters:
* **info_df:** *(Pandas Dataframe)* has columns for the vertices of ternary

* **top_col:** *(str)* Name of top vertex

* **left_col:** *(str)* Name of left vertex

* **right_col:** *(str)* Name of right vertex


### plot_ternary()
#### Description:
* Plots a ternary plot

#### Parameters:
* **data_specs:** *(list of dictionaries)* Has details about each data to be plotted keys are 'Data','Color','Label'

* **title:** *(str)* Title of the plot

* **fig:** *(Matplotlib Figure)* Figure to plot on [Optional]

* **ax:** *(Matplotlib Axis)* Axis to plot on [Optional]

* **t_scale:** *(float)* Text scale up or down [Optional]

* **location:** *(str)* Where the plot should be saved, if it is to be saved [Optional]

* **t:** *(str)* Title of the top vertex [Optional]

* **l:** *(str)* Title of the left vertex [Optional]

* **r:** *(str)* Title of the right vertex [Optional]


### alter_length()
#### Description:
* Changes the length of the vector for a given proportion

#### Parameters:
* **x_vals =** *(list)* [x_0,x_1]

* **y_vals =** *(list)* [y_0,y_1[


### get_points_dyn()
#### Description:
* Gets the euclidean coordinates of the dynamic ternary plot vectors

#### Parameters:
* **info_df:** *(Pandas Dataframe)* has columns for the vertices of ternary

* **head:** *(str)* Name of the cell type that is the head of the vectors

* **tail:** *(str)* Name of the cell type that is the tail of the vectors

* **t_suff:** *(str)* Suffix of the column in info_df for top vertex data [Optional]

* **l_suff:** *(str)* Suffix of the column in info_df for left vertex data [Optional]

* **r_suff:** *(str)* Suffix of the column in info_df for right vertex data [Optional]


### plot_dyn_ternary()
#### Description:
* Plots a ternary plot

#### Parameters:
* **dyn_points:** *(list of lists)* Points from get_points_dyn

* **title:** *(str)* Title of the plot

* **fig:** *(Matplotlib Figure)* Figure to plot on [Optional]

* **ax:** *(Matplotlib Axis)* Axis to plot on [Optional]

* **t_scale:** *(float)* Text scale up or down [Optional]

* **cmap:** *(str)* Matplotlib colormap [Optional]

* **location:** *(str)* Where the plot should be saved, if it is to be saved [Optional]

* **t:** *(str)* Title of the top vertex [Optional]

* **l:** *(str)* Title of the left vertex [Optional]

* **r:** *(str)* Title of the right vertex [Optional]

* **spec_lab:** *(str)* Label for the angle spectrum [Optional]

## rna_plot_tools/rna_plot_tools.py

### am_gen()
#### Description:
* Generates the geometric mean of the proportion of the window exons take up for a given dataframe (which should represent one plot)

#### Parameters:
* **data_df:** *(Pandas DataFrame)* The DataFrame created using basic_tools text_to_df for the UCSC data, has a 'name2' column


### scale()
#### Description:
* Takes in a dataframe of ucsc info to be plotted, and creates the plottable segments without augmentation, returns the plottable dataframe

#### Parameters:
* **df:** *(Pandas DataFrame)* The DataFrame created by the filter on the UCSC data_df, has one bin, one strand direction, and one rna attributed to it


### augment()
#### Description:
* Takes in a dataframe of ucsc info to be plotted, and creates the plottable segments with augmentation, returns the plottable dataframe

#### Parameters:
* **df:** *(Pandas DataFrame)* The DataFrame created by the filter on the UCSC data_df, has one bin, one strand direction, and one rna attributed to it

* **am:** *(float)* The arithmetic mean generated by am_gen function, must be below 0.02 to be useful, used as the index to scale up exons


### plot_seg()
#### Description:
* Plots the exons on a single segment of the rna, each rna has 3 segments: 1 is before the cdsStart, 2 is between cdsStart and cdsEnd, 3 is past cdsEnd

#### Parameters:
* **ax:** *(Matplotlib Axis)* The axis to plot the seg onto

* **r:** *(int)* The y value to plot the seg onto

* **seg:** *(list)* The list of exon start and ends to plot

* **exons:** *(list of lists)* The list of exons for ex_i, the list has an element [s_i,e_i]

* **c:** *(str)* Valid matplotlib color to plot, for seg 2's it should be Crimson, for seg 1 and 3's it should be #8c0d26

* **w:** *(float)* The width of the line, for seg 2's it should be 20, for seg 1 and 3's it should be 10


### ucsc_plot()
#### Description:
* Produces all the plot for the UCSC data of a single rna, separates plots by bin and strand, augments if am is below 0.02

#### Parameters:
* **rna:** *(str)* The RNA name, must be one from the 'name2' column in the UCSC dataframe

* **data:** *(str)* The UCSC DataFrame created from the text_to_df function

* **fname:** *(str)* If file output to pdf, put a special fname for this rna, separate plots will be produced and labeled with corresponding bin and strand [Optional]

* **override:** *(boolean)* Set to True if you want the unscaled version even if it is under the 0.02 threshold [Optional]

* **excd:** *(str)* Matplotlib compatible color for exons in coding region [Optional]

* **extx:** *(str)* Matplotlib compatible color for exons out of coding region [Optional]

* **incd:** *(str)* Matplotlib compatible color for introns in coding region [Optional]

* **intx:** *(str)* Matplotlib compatible color for introns out of coding region [Optional]

## rna_path_tools/rna_path_tools.py

### get_rna_dfs()
#### Description:
* Filters the UCSC dataframe for a specified RNA in the name2 column, gives a list of one or two dataframes

#### Parameters:
* **rna:** *(str)* The name of the RNA in the name2 column

* **data_df:** *(Pandas DataFrame)* The UCSC dataframe from text_to_df

#### Returns:
* **dfs:** *(list)* A list of one or two filtered Pandas DataFrames, one for each strand

### get_lookup1()
#### Description:
* Gets the first lookup table for the given RNA DataFrame produced by get_rna_dfs

#### Parameters:
* **df:** *(Pandas DataFrame)* The RNA DataFrame produced by get_rna_dfs

#### Returns:
* **lookup_df:** *(Pandas DataFrame)* The first lookup table for exon ID, has upper and lower bounds for each exon, all exons, and number of diff forms per exon

### get_lookup2()
#### Description:
* Gets the second lookup table for the given RNA DataFrame produced by get_rna_dfs

#### Parameters:
* **df:** *(Pandas DataFrame)* The RNA DataFrame produced by get_rna_dfs

#### Returns:
* **lu2:** *(Pandas DataFrame)* The second lookup table, masks the multiple splice sites as different exons

### get_pexons()
#### Description:
* Gets the psuedoexons for the given RNA DataFrame produced by get_rna_dfs

#### Parameters:
* **df:** *(Pandas DataFrame)* The RNA DataFrame produced by get_rna_dfs

* **lu_df:** *(Pandas DataFrame)* If lookup table 2 has already been produced, you can input it [Optional]

#### Returns:
* **strand_nodes:** *(list)* A list of the psuedoexons for conversion to a graph
* **detailed:** *(list)* includes the real cdsStart and cdsEnd for a strand
* **names:** *(list)* includes the strand name

### get_graph()
#### Description:
* Gets the graph for a given list of nodes

#### Parameters:
* **strand_nodes:** *(list)* A list of the psuedoexons for conversion to a graph created by get_pexons

#### Returns:
* **G:** *(Networkx DiGraph)* The digraph created by the list of nodes
* **se_pairs:** *(list of lists)* List of all possible start,end pairs in the graph

### get_all_paths()
#### Description:
* Gets all paths with replacement of the pseudoexons, adds an '*' if it exists

#### Parameters:
* **rna:** *(str)* The name of the RNA in the name2 column

* **data_df:** *(Pandas DataFrame)* The UCSC dataframe from text_to_df

* **detail:** *(boolean)* False if output is nodes, True for start and ends [Optional]

#### Returns:
* **paths:** *(list)* List of all paths produced
## analysis_tools/analysis_tools.py

### go_to_bed()
#### Description:
* Makes the BED Files for all possible paths and outputs to bed and bedinfo directories, will consider existing directories and enumerate, also labels with RNA Name

#### Parameters:
* **rna:** *(str)* The name of the RNA in the name2 column

* **raw:** *(str)* File location of UCSC data


### seq_index()
#### Description:
* Convert ind from genomic to sequence coordinates based on nodes and strand

#### Parameters:
* **ind:** *(int)* The genomic coordinate, usually a long number from UCSC browser

* **nodes:** *(list of int)* The genomic coordinates of the exons start/ends

* **strand:** *(str)* + or - depending on the transcript

#### Returns:
* **to_return:** *(int)* The sequence coordinate

### gene_index()
#### Description:
* Convert ind from sequence to genomic coordinates based on nodes and strand

#### Parameters:
* **ind:** *(int)* The sequence coordinate, usually a long number from UCSC browser

* **nodes:** *(list of int)* The genomic coordinates of the exons start/ends

* **strand:** *(str)* + or - depending on the transcript

#### Returns:
* **to_return:** *(int)* The genomic coordinate

### fetch_coords()
#### Description:
* Gets the location of the start and stop codon

#### Parameters:
* **seq:** *(str)* The string of A,T,G,C from the bed file

#### Returns:
* **se:** *(list of 2 int)* first value is the start, second value is the stop
