
import pandas as pd
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt

import os

def plot_clustered_stacked(dfall, labels=None, title="multiple stacked bar plot",  H="", alpha=1, **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
labels is a list of the names of the dataframe, used for the legend
title is a string for the title of the plot
H is the hatch used for identification of the different dataframe"""

    n_df = len(dfall)
    n_col = len(dfall[0].columns) 
    n_ind = len(dfall[0].index)
    axe = plt.subplot(111)

    for df in dfall : # for each data frame
        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ax=axe,
                      legend=False,
                      grid=False,
                      #width=[1,1,1,1,1,1,2],
                      **kwargs)  # make bar plots

    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x((rect.get_x() + 1 / float(n_df + 1) * (i) / float(n_col))+.03*int(i / n_col))
                rect.set_hatch(H * int(i / n_col)) #edited part  
                #rect.set_alpha(alpha/(1+int(i / n_col)))
                if j==len(h[i:i+n_col]):
                    rect.set_width(2 / float(n_df + 1)) 
                else:
                    rect.set_width(1 / float(n_df + 1))

    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    axe.set_xticklabels(df.index, rotation = 0)
    axe.set_title(title)

    # Add invisible data to add another legend
    n=[]        
    #for i in range(n_df):
    #    n.append(axe.bar(0, 0, color="gray", hatch=H * i))#, width=[1,1,1,1,1,1,2]))

    """
    handles, leg_labels = axe.get_legend_handles_labels()
    print(leg_labels)
    ordering = range(len(leg_labels))

    x, handles = zip(*sorted(zip(ordering, handles), key=lambda t: t[0], reverse=True))
    x, leg_labels = zip(*sorted(zip(ordering, leg_labels), key=lambda t: t[0], reverse=True))
    print(leg_labels)
    axe.legend(handles, leg_labels)
    """

    print(l)
    print(l[:n_col][::-1])
    l1 = axe.legend(h[:n_col][::-1], l[:n_col][::-1], loc=[.01, 0.01])

    if labels is not None:
        #print(labels)
        #raise(ValueError)
        l2 = plt.legend(n, labels, loc=[.3, 0.01]) 
    axe.add_artist(l1)
    return axe



#Import CSVs with counts


#then try it two ways - all asharpened together and all unsharpened together, and seperate probes together

#place counts at the top



def plot_bleeding(path):
    with open(path, 'r') as f:
        df = pd.read_csv(f)

    df = df.set_index('Description')
    for column in df.columns:
        if 'Unnamed' in column:
            df = df.drop(column, axis=1)
    #fig, ax = plt.subplots()


    df = df.replace(0,.001)

    sharp_df = pd.DataFrame(columns= df.columns)
    unsharp_df = pd.DataFrame(columns= df.columns)

    area_mapping = {
        'probe A': 'visAM',
        'probe B': 'visPM',
        'probe C': 'visP',
        'probe D': 'visLM',
        'probe E': 'visAL',
        'probe F': 'visRL',
        'All Probes': 'All 6 areas'
    }

    for index, row in df.iterrows():
        print(row)
        new_name = (' ').join(index.split(' ')[:2])
        try:
            new_name = area_mapping[new_name]
        except KeyError as E:
            pass
        row.name = new_name
        print(row.name)
        if 'unsharp' in index.lower():
            unsharp_df = unsharp_df.append(row)
        else:
            sharp_df = sharp_df.append(row)

    sharp_counts = sharp_df['N_insertions'].astype(int)
    sharp_df = sharp_df.drop('N_insertions', axis=1)
    #sharp_counts.append(sum(sharp_counts))
    unsharp_counts = unsharp_df['N_insertions'].astype(int)
   # unsharp_counts.append(sum(unsharp_counts))
    unsharp_df = unsharp_df.drop('N_insertions', axis=1)

    print(sharp_df.index)
    print(sharp_df)
    print(unsharp_df)

    ax = plot_clustered_stacked([sharp_df, unsharp_df],
                           ["sharpened", "unsharpened"],
                           cmap=plt.cm.cool)

    print(ax.get_xticklabels())
    
    ax.set_xticklabels(ax.get_xticklabels(), rotation = -45)
    ax.set_ylabel('Fraction of Insertions')

    """
    handles, labels = ax.get_legend_handles_labels()
    print(labels)
    ordering = rang(len(labels))

    x, handles = zip(*sorted(zip(ordering, handles), key=lambda t: t[0], reverse=True))
    x, labels = zip(*sorted(zip(ordering, labels), key=lambda t: t[0], reverse=True))
    ax.legend(handles, labels)
    """


    ##something like this to add the counts
    #ax.bar_label(p1, label_type='center')
    from matplotlib.container import BarContainer

    bars = [i for i in ax.containers if isinstance(i, BarContainer)]
    print(len(bars))
    for idx, bar in enumerate(bars):
        if idx==4:
            sharp_counts = [f'S\n({i})' for i in sharp_counts]
            ax.bar_label(bar, labels=sharp_counts, label_type='edge', size=7)
        if idx==9:
            unsharp_counts = [f'U\n({i})' for i in unsharp_counts]
            ax.bar_label(bar, labels=unsharp_counts, label_type='edge', size=7)
    print('##############################################################################')
    

    fig = ax.get_figure()
    plt.tight_layout()
    plt.show()
    #raise(ValueError)
    keyword = os.path.split(path)[1].split('_')[0]
    save_path = os.path.join(r'C:\Users\greggh\Documents\python_scripts\sharpened_probes\plots', keyword+'_bar_clustered.pdf')

    fig.savefig(save_path)

paths = [
    r"C:\Users\greggh\Documents\python_scripts\sharpened_probes\mean_fraction.csv",
    r"C:\Users\greggh\Documents\python_scripts\sharpened_probes\max_fraction.csv",
    ]

for path in paths:
    plot_bleeding(path)

rasie(ValueError)
#VV this sorts the columns so sharpened is next to eachother, it didn't look better
"""
ordering = [1,8,2,9,3,10,4,11,5,12,6,13,7,14]

df['ordering'] = ordering

df = df.sort_values(by=['ordering'])
df = df.drop('ordering', axis=1)
"""
print(df.index)
print(df.head(5))

ax = df.plot.bar(stacked=True)

handles, labels = ax.get_legend_handles_labels()
print(labels)
ordering = [1,2,3,4,5]

x, handles = zip(*sorted(zip(ordering, handles), key=lambda t: t[0], reverse=True))
x, labels = zip(*sorted(zip(ordering, labels), key=lambda t: t[0], reverse=True))
ax.legend(handles, labels)

#a = sorted(a, key=lambda x: x.modified, reverse=True)

#something like this should work if we want to put in the N somewhere in the middle? at the bottom?
"""
for idx, p in enumerate(ax.patches):
    print(idx)
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
"""



print(labels)

fig = ax.get_figure()
fig.savefig(r'C:\Users\greggh\Documents\python_scripts\sharpened_probes\plots\max_bar_plot.png')