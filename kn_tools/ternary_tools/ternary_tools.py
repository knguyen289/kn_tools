import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import math

def get_mag(x_vals,y_vals):
    """Gets the magnitude of a vector
    
    Parameters:
        x_vals: (list) [x_0,x_1]
        y_vals: (list) [y_0,y_1]
    """
    dx = x_vals[1] - x_vals[0]
    dy = y_vals[1] - y_vals[0]
    mag = (dx**2 + dy**2)**0.5
    return mag

def get_angle(x_vals,y_vals):
    """Gets the degree angle of a vector
    
    Parameters:
        x_vals: (list) [x_0,x_1]
        y_vals: (list) [y_0,y_1]
    """
    dx = x_vals[1] - x_vals[0]
    dy = y_vals[1] - y_vals[0]
    if dx == 0:
        if dy > 0:
            return 90
        else:
            return 270
    a = math.atan(dy/dx)*180/math.pi
    
    if dx > 0 and dy > 0:
        return a
    if dx < 0 and dy > 0:
        return 180 + a
    if dx < 0 and dy < 0:
        return 180 + a
    return 360 + a

def get_points(info_df,top_col,left_col,right_col):
    """Gets the euclidean coordinates of the ternary plot points along with the ordered gene names
    
    Parameters:
        info_df: (Pandas Dataframe) has columns for the vertices of ternary
        top_col: (str) Name of top vertex
        left_col: (str) Name of left vertex
        right_col: (str) Name of right vertex
    """
    #up
    top = [0,1]
    #right
    left = [((3**0.5)/2),-0.5]
    #left
    right = [(-(3**0.5)/2),-0.5]
    
    points = []
    genes = []
    for index,row in info_df.iterrows():
        t = float(row.get_value(top_col))
        l = float(row.get_value(left_col))
        r = float(row.get_value(right_col))
        total = t + l + r
        if min(t,l,r) != 0:
            t /= total
            l /= total
            r /= total
            row_sum = [top[0]*t + left[0]*l + right[0]*r,
                       top[1]*t + left[1]*l + right[1]*r]
            points.append(row_sum)
            genes.append(index)
            
    all_x = [item[0] for item in points]
    all_y = [item[1] for item in points]
    return [all_x,all_y],genes

def plot_ternary(data_specs,title,fig=None,ax=None,t_scale=1,location=None,t=r'Membrane',l=r'Insoluble',r=r'Cytosolic'):
    """Plots a ternary plot

    Parameters:
        data_specs: (list of dictionaries) Has details about each data to be plotted
                    keys are 'Data','Color','Label'
        title: (str) Title of the plot
        fig: (Matplotlib Figure) Figure to plot on [Optional]
        ax: (Matplotlib Axis) Axis to plot on [Optional]
        t_scale: (float) Text scale up or down [Optional]
        location: (str) Where the plot should be saved, if it is to be saved [Optional]
        t: (str) Title of the top vertex [Optional]
        l: (str) Title of the left vertex [Optional]
        r: (str) Title of the right vertex [Optional]
    """
    if fig == None and ax == None:
        fig,ax = plt.subplots(figsize=(14,7*(3**0.5)))
    numbers = ''
    
    for spec in data_specs:
        data = spec['Data']
        color = spec['Color']
        label = spec['Label']
        
        if label != None:
            ax.plot(data[0],
                    data[1],
                    marker='.',
                    c=color,
                    markeredgecolor=color,
                    linestyle='None',
                    label=label + ' (' + str(len(data[0])) + ' points)')
        
    #triangle
    ax.plot([((3**0.5)/2),0],[-0.5,1],c='black',linewidth=2*t_scale)
    ax.plot([(-(3**0.5)/2),0],[-0.5,1],c='black',linewidth=2*t_scale)
    ax.plot([((3**0.5)/2),(-(3**0.5)/2)],[-0.5,-0.5],c='black',linewidth=2*t_scale)
    ax.grid(b=False)
    ax.axis('off')
    ax.set_xlim([(-(3**0.5)/2)-0.25,((3**0.5)/2)+0.25])
    ax.set_ylim([-0.5-0.25,1+0.25])

    #labels
    ax.text(0,1.05,t,fontsize=28*t_scale,weight='demi',horizontalalignment='center')
    ax.text((3**0.5)/2,-0.6,l,fontsize=28*t_scale,weight='demi',horizontalalignment='center')
    ax.text(-(3**0.5)/2,-0.6,r,fontsize=28*t_scale,weight='demi',horizontalalignment='center')
    
    legend = ax.legend(fontsize=17*t_scale,loc=2,frameon=True)
    frame = legend.get_frame()
    frame.set_edgecolor('black')

    ax.set_title(title,fontsize=40*t_scale,weight='bold')
    if location != None:
        fig.savefig(location)

def alter_length(x_vals,y_vals,p):
    """Changes the length of the vector for a given proportion

    Parameters:
        x_vals = (list) [x_0,x_1]
        y_vals = (list) [y_0,y_1[
        p = (float) proportion ex: 10 becomes a tenth of length"""
    dx = float(x_vals[1] - x_vals[0])
    dy = float(y_vals[1] - y_vals[0])
    
    new_xs = [x_vals[0],x_vals[0] + dx/float(p)]
    new_ys = [y_vals[0],y_vals[0] + dy/float(p)]
    return [new_xs,new_ys]

def get_points_dyn(info_df,head,tail,t_suff='_mem',l_suff='_ins',r_suff='_cyt'):
    """Gets the euclidean coordinates of the dynamic ternary plot vectors
    
    Parameters:
        info_df: (Pandas Dataframe) has columns for the vertices of ternary
        head: (str) Name of the cell type that is the head of the vectors
        tail: (str) Name of the cell type that is the tail of the vectors
        t_suff: (str) Suffix of the column in info_df for top vertex data [Optional]
        l_suff: (str) Suffix of the column in info_df for left vertex data [Optional]
        r_suff: (str) Suffix of the column in info_df for right vertex data [Optional]
    """
    #up
    top = [0,1]
    #right
    left = [((3**0.5)/2),-0.5]
    #left
    right = [(-(3**0.5)/2),-0.5]
    
    points = []
    genes = []
    for index,row in info_df.iterrows():
        xyxy = []
        for p in [head,tail]:
            top_col = p + t_suff
            left_col = p + l_suff
            right_col = p + r_suff
            
            t = float(row.get_value(top_col))
            l = float(row.get_value(left_col))
            r = float(row.get_value(right_col))
            total = t + l + r
            if min(t,l,r) != 0:
                t /= total
                l /= total
                r /= total
                row_sum = [top[0]*t + left[0]*l + right[0]*r,
                           top[1]*t + left[1]*l + right[1]*r]
                xyxy.append(row_sum[0])
                xyxy.append(row_sum[1])
        if len(xyxy) == 4:
            points.append(xyxy)
            genes.append(index)
    all_x = [[item[0],item[2]] for item in points]
    all_y = [[item[1],item[3]] for item in points]
    return [all_x,all_y]

def plot_dyn_ternary(dyn_points,title,fig=None,ax=None,t_scale=1,cmap='hsv',location=None,t=r'Membrane',l=r'Insoluble',r=r'Cytosolic',spec_lab='Angle Spectrum'):
    """Plots a ternary plot

    Parameters:
        dyn_points: (list of lists) Points from get_points_dyn
        title: (str) Title of the plot
        fig: (Matplotlib Figure) Figure to plot on [Optional]
        ax: (Matplotlib Axis) Axis to plot on [Optional]
        t_scale: (float) Text scale up or down [Optional]
        cmap: (str) Matplotlib colormap [Optional]
        location: (str) Where the plot should be saved, if it is to be saved [Optional]
        t: (str) Title of the top vertex [Optional]
        l: (str) Title of the left vertex [Optional]
        r: (str) Title of the right vertex [Optional]
        spec_lab: (str) Label for the angle spectrum [Optional]
    """
    hsv = plt.get_cmap(cmap)
    cNorm = colors.Normalize(vmin=0,vmax=360)
    scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=hsv)
    
    if fig == None and ax == None:
        fig,ax = plt.subplots(figsize=(14,7*(3**0.5)))
    for i in range(len(dyn_points[0])):
        v = [dyn_points[0][i],dyn_points[1][i]]
        x1 = v[0][0]
        x2 = v[0][1]
        y1 = v[1][0]
        y2 = v[1][1]
        angle = get_angle(v[0],v[1])
        colorVal = scalarMap.to_rgba(angle)
        temp = list(colorVal)
        for i in range(3):
            temp[i] *= 0.9
        colorVal = list(temp)
        if get_mag(v[0],v[1]) > 0.05:
            v = alter_length([x1,x2],[y1,y2],15)
            x1 = v[0][0]
            x2 = v[0][1]
            y1 = v[1][0]
            y2 = v[1][1]
            dx = x2 - x1
            dy = y2 - y1
            ax.arrow(x1,y1,dx,dy,color=colorVal,width=0.0005,length_includes_head=True)
        else:
            ax.plot(x1,y1,marker='.',c=colorVal,markeredgecolor=colorVal,markersize=6)

    #triangle
    ax.plot([((3**0.5)/2),0],[-0.5,1],c='black',lw=2)
    ax.plot([(-(3**0.5)/2),0],[-0.5,1],c='black',lw=2)
    ax.plot([((3**0.5)/2),(-(3**0.5)/2)],[-0.5,-0.5],c='black',lw=2)
    ax.grid(b=False)
    ax.axis('off')
    ax.set_xlim([(-(3**0.5)/2)-0.2,((3**0.5)/2)+0.25])
    ax.set_ylim([-0.5-0.2,1+0.2])

    #labels
    ax.text(0,1.05,t,fontsize=28*t_scale,weight='demi',horizontalalignment='center')
    ax.text((3**0.5)/2,-0.6,l,fontsize=28*t_scale,weight='demi',horizontalalignment='center')
    ax.text(-(3**0.5)/2,-0.6,r,fontsize=28*t_scale,weight='demi',horizontalalignment='center')

    ax.text(0.9,0.48,spec_lab,fontsize=20*t_scale,weight='demi',horizontalalignment='center')
    for i in range(0,360,3):
        x = math.cos(math.radians(i))
        y = math.sin(math.radians(i))
        colorVal = scalarMap.to_rgba(i)
        temp = list(colorVal)
        for i in range(3):
            temp[i] *= 0.9
        colorVal = list(temp)
        ax.plot((x/12 + 0.9,x/6 + 0.9),(y/12 + 0.23,y/6+0.23),color=colorVal)

    ax.set_title(title,fontsize=40*t_scale,weight='bold')
    if location != None:
        fig.savefig(location)

