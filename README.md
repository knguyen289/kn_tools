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
* Converts a text file that is tab delimited with first row being header to Pandas DataFrame

#### Parameters:
* **filename:** *(str)* The location of the text file
* **index:** *(str)* The column that can be used as an index [Optional]

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

## rna_tools/rna_tools.py

### rna_vis()
#### Description:
* Creates a bar visualization for an rna strand with UCSC data

#### Parameters:
* **rna:** *(str)* The RNA name in the 'name2' column of the UCSC data

* **data:** *(Pandas DataFrame)* The DataFrame created using basic_tools text_to_df for the UCSC data, has a 'name2' column

* **fname:** *(str)* The filename of the figure saved (no suffix) [Optional]

* **excd:** *(str)* HTML color name for exons within coding region (default: 'Crimson') [Optional]

* **extx:** *(str)* HTML color name for exons outside coding region (default: '#8c0d26') [Optional]

* **incd:** *(str)* HTML color name for introns within coding region (default: '#ffbab3') [Optional]

* **intx:** *(str)* HTML color name for introns outside coding region (default: 'MistyRose') [Optional]

