# kn_tools
Kim Nguyen's tools for Wang Lab, includes tools for plotting and data management
## ternary_tools/ternary_tools.py

### get_mag()
#### Description:
* Gets the magnitude of a vector

#### Parameters:
* x_vals: (list) [x_0,x_1]
* y_vals: (list) [y_0,y_1]

### get_angle()
#### Description:
* Gets the degree angle of a vector

#### Parameters:
* x_vals: (list) [x_0,x_1]
* y_vals: (list) [y_0,y_1]

### get_points()
#### Description:
* Gets the euclidean coordinates of the ternary plot points along with the ordered gene names

#### Parameters:
* info_df: (Pandas Dataframe) has columns for the vertices of ternary
* top_col: (str) Name of top vertex
* left_col: (str) Name of left vertex
* right_col: (str) Name of right vertex

### plot_ternary()
#### Description:
* Plots a ternary plot

#### Parameters:
* data_specs: (list of dictionaries) Has details about each data to be plotted
* keys are 'Data','Color','Label'
* title: (str) Title of the plot
* fig: (Matplotlib Figure) Figure to plot on [Optional]
* ax: (Matplotlib Axis) Axis to plot on [Optional]
* t_scale: (float) Text scale up or down [Optional]
* location: (str) Where the plot should be saved, if it is to be saved [Optional]
* t: (str) Title of the top vertex [Optional]
* l: (str) Title of the left vertex [Optional]
* r: (str) Title of the right vertex [Optional]

### alter_length()
#### Description:
* Changes the length of the vector for a given proportion

#### Parameters:
* x_vals = (list) [x_0,x_1]
* y_vals = (list) [y_0,y_1[

### get_points_dyn()
#### Description:
* Gets the euclidean coordinates of the dynamic ternary plot vectors

#### Parameters:
* info_df: (Pandas Dataframe) has columns for the vertices of ternary
* head: (str) Name of the cell type that is the head of the vectors
* tail: (str) Name of the cell type that is the tail of the vectors
* t_suff: (str) Suffix of the column in info_df for top vertex data [Optional]
* l_suff: (str) Suffix of the column in info_df for left vertex data [Optional]
* r_suff: (str) Suffix of the column in info_df for right vertex data [Optional]

## basic_tools/basic_tools.py

### filter_df()
#### Description:
* Returns a dataframe that includes info from indicated indexes

#### Parameters:
* df: (Pandas DataFrame) The original DataFrame
* indexes: (list of str) The list of indexes

### redact_df()
#### Description:
* Returns a dataframe that includes info from indicated columns

#### Parameters:
* df: (Pandas DataFrame) The original DataFrame
* indexes: (list of str) The list of column names

