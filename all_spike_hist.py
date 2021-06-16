
import os, glob
from matplotlib import pyplot as plt
import numpy as np
import re


save_fig_dir = r"\\allen\programs\braintv\workgroups\nc-ophys\corbettb\all_spike_hists\sahar"


dirpaths = [r"\\10.128.54.155\Data"]

failed = []
count = 0
for dirpath in dirpaths:
    for d,_,_ in os.walk(dirpath):
        stfiles = glob.glob(os.path.join(d, 'spike_times.npy'))
        count+=1
        if count>1000:
            print('walking')
            count=0
        if len(stfiles):
            print('found spike times', d)
            sts_file = stfiles[0]
            try:
                sts = np.load(sts_file).flatten()
                h, b = np.histogram(sts/30000., bins=np.arange(sts[-1]/30000.))
                fig, ax = plt.subplots()
                title = '_'.join(d.split('\\')).replace('.', '')
                fig.suptitle(title)
                ax.plot(b[:-1], h)
                save_path = os.path.join(save_fig_dir, title)
                print(save_path)
                fig.savefig(save_path)
                plt.close('all')
                #raise(ValueError)
            except Exception as E:
                #raise(E)
                print('Error plotting', str(sts_file))
                pass