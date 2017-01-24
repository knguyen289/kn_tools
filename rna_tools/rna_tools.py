import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np 
import pandas as pd 
import copy

def rna_vis(rna,data,fname=None,excd='Crimson',extx='DARKRED',incd='PINK',intx='MistyRose'):
    """Creates a bar visualization for an rna strand with UCSC data

    Parameters:
        rna: (str) The RNA name in the 'name2' column of the UCSC data
        data: (Pandas DataFrame) The DataFrame created using basic_tools text_to_df for the UCSC data, has a 'name2' column
        fname: (str) The filename of the figure saved (must be .png) [Optional]
        excd: (str) HTML color name for exons within coding region (default: 'Crimson') [Optional]
        extx: (str) HTML color name for exons outside coding region (default: 'DARKRED') [Optional]
        incd: (str) HTML color name for introns within coding region (default: 'PINK') [Optional]
        intx: (str) HTML color name for introns outside coding region (default: 'MistyRose') [Optional]
    """
    sns.set_style('whitegrid')
    temp_df = copy.deepcopy(data)
    temp_df = temp_df[temp_df['name2'] == rna]
    try:
        fig,ax = plt.subplots()
        cols = [extx,intx,excd,incd]
        r = 0

        for index,row in temp_df.iterrows():
            tx = [int(row.get_value('txStart')),int(row.get_value('txEnd'))]
            cds = [int(row.get_value('cdsStart')),int(row.get_value('cdsEnd'))]
            if row.get_value('strand') == '+':
                sect = tx[1] - cds[1]
            else:
                sect = cds[0] - tx[0]

            s = row.get_value('exonEnds')[:-1].split(',')
            s = sorted([int(item) for item in s])
            e = row.get_value('exonStarts')[:-1].split(',')
            e = sorted([int(item) for item in e])
            nodes = sorted(s+e)
            
            exons = []
            for i in range(len(s)):
                exons.append(sorted([s[i],e[i]]))
                
            ch1 = sorted(list(set([tx[0]] + [item for item in nodes if item <= cds[0]] + [cds[0]])))
            ch2 = sorted(list(set([cds[0]] + [item for item in nodes if item > cds[0] and item <= cds[1]] + [cds[1]])))
            ch3 = sorted(list(set([cds[1]] + [item for item in nodes if item > cds[1]] + [tx[1]])))

            for i in range(len(ch1)-1):
                c = 1
                w = 5
                for exon in exons:
                    if ch1[i] >= exon[0] and ch1[i+1] <= exon[1]:  
                        c = 0
                        pass
                plt.hlines(r,ch1[i],ch1[i+1],colors=cols[c],linewidth=w)

            for i in range(len(ch2)-1):
                c = 3
                w = 5
                for exon in exons:
                    if ch2[i] >= exon[0] and ch2[i+1] <= exon[1]:  
                        c = 2
                        w = 15
                        pass
                plt.hlines(r,ch2[i],ch2[i+1],colors=cols[c],linewidth=w)

            for i in range(len(ch3)-1):
                c = 1
                w = 5
                for exon in exons:
                    if ch3[i] >= exon[0] and ch3[i+1] <= exon[1]:  
                        c = 0
                        pass
                plt.hlines(r,ch3[i],ch3[i+1],colors=cols[c],linewidth=w)
            r += 1

        if row.get_value('strand') == '+':
            ax.set_xlim([tx[0]-100,tx[1]+100])


        else:
            ax.set_xlim([tx[1]+100,tx[0]-100])

        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)
        ax.set_ylim([-1,len(temp_df.index)])
        ax.set_title(rna + ' (' + row.get_value('strand') + ') length: ' + str(tx[1] - tx[0]))
        fig.set_size_inches(12,len(temp_df.index)*0.5+1)
        plt.tight_layout()
        if fname == None:
            plt.show()
        else:
            plt.savefig(fname)
    except:
        plt.close()
        return 'Plot could not be produced'
