
import os, glob
from matplotlib import pyplot as plt
import numpy as np
import re


save_fig_dir = r"\\allen\programs\braintv\workgroups\nc-ophys\corbettb\all_spike_hists\settle_data"


dirpaths = [r"\\10.128.54.20\sd8", r"\\10.128.54.20\sd8.2"]
sessions_to_run = []
"""
for d,_,_ in os.walk(dirpath):
    stfiles = glob.glob(os.path.join(d, 'spike_times.npy'))
    try:
        lower_or_settle = 'lower' in stfiles[0] or 'settle' in stfiles[0]
        if len(stfiles)>0 and lower_or_settle:
            sessions_to_run.append(glob.glob(os.path.join(d, 'spike_times.npy')))
    except Exception as E:
        pass
"""

for dirpath in dirpaths:
    for dirname in os.listdir(dirpath):
        lower_or_settle = 'lower' in dirname or 'settle' in dirname
        if lower_or_settle:
            sessions_to_run.append(os.path.join(dirpath, dirname))
#print(sessions_to_run)
#raise(ValueError)

def get_corresponding_probe_spike_times(sts_file):
    probe = sts_file[sts_file.find('probe'):sts_file.find('probe')+6]
    session_dir = sts_file
    while True:
        session_dir, sort_dirname = os.path.split(session_dir)
        if not(probe in session_dir):
            break
    session = get_corresponding_session(session_dir)
    session_name = os.path.split(session)[1]
    sorted_dir = os.path.join(session, session_name+'_'+probe+'_sorted')
    st_file = None
    for d,_,_ in os.walk(sorted_dir):

        stfiles = glob.glob(os.path.join(d, 'spike_times.npy'))
        if len(stfiles):
            print('found spike times', d)
            st_file = stfiles[0]
    return st_file, sort_dirname

 
def get_corresponding_session(dirpath):
    print('dirpath: ', dirpath)
    dirname = os.path.split(dirpath)[1]
    print('dirname: ', dirname)
    mouse_date = '_'.join(dirname.split('_')[1:])
    print('mouse_date: ', mouse_date)
    dir_list = []
    for dirpath in dirpaths:
        for dirname in os.listdir(dirpath):
            if mouse_date in dirname and not('lower' in dirname or 'settle' in dirname):
                dir_list.append(os.path.join(dirpath, dirname))
    print('dir_list: ', str(dir_list))
    assert(len(dir_list)==1)
    return dir_list[0]


include_samples = 8*30000*60

failed = []
get_mouse_id = lambda x: re.search('[0-9]{6}', x)
for ind, session_path in enumerate(sessions_to_run):
    #session_path = os.path.join(dirpath, s)
    print(session_path)
    if os.path.isdir(session_path):
        print('running')
        for dirname in os.listdir(session_path):
            if True:#'sorted' in dirname:
                print('sorted')
                sort_path = path = os.path.join(session_path, dirname)
                for d,_,_ in os.walk(sort_path):

                    stfiles = glob.glob(os.path.join(d, 'spike_times.npy'))
                    if len(stfiles):
                        print('found spike times', d)
                        sts_file = stfiles[0]
                        try:
                            
                            session_sts_path, sort_dirname = get_corresponding_probe_spike_times(sts_file)
                            print('session_sts_path: ', session_sts_path)
                            #mid = get_mouse_id(sts_file)
                            #mid = mid.group(0) if mid is not None else ind
                            #mid = '553963'#dirname
                            #print(mid)
                            if sts_file is not None:
                                sts = np.load(sts_file).flatten()
                                session_sts = np.load(session_sts_path).flatten()
                                max_st = sts[-1]
                                print('max_st: ', max_st)
                                print('max_st2: ', np.max(sts))
                                start_st = max_st - include_samples
                                print('start_st: ', start_st)

                                session_sts = session_sts+sts[-1]
                                sts = np.append(sts, session_sts)
                                sts = sts[sts>start_st]
      
                                print('here1')
                                h, b = np.histogram(sts/30000., bins=np.arange(sts[-1]/30000.))
                                fig, ax = plt.subplots()
                                title = sort_dirname+'_full'
                                fig.suptitle(title)
                                ax.plot(b[:-1], h)
                                fig.savefig(os.path.join(save_fig_dir, 'full', title))
                                plt.close('all')
                                print('here2')
                                end_st = max_st + include_samples
                                print(end_st)

                                #eight_sts = sts[sts>start_st]
                                print('max_st3: ', np.max(sts))
                                eight_sts = sts[sts<end_st]
                                print('max_st4: ', np.max(eight_sts))
                                #eight_sts = eight_sts - start_st
                                print('eight_stslen: ', len(eight_sts))
                                #eight_sts = []
                                #for st in sts:
                                #    if st> start_st and st< end_st:
                                #       eight_sts.append(st - start_st)
                                print('here3')
                                h, b = np.histogram(eight_sts/30000., bins=np.arange(np.max(eight_sts)/30000.))
                                fig, ax = plt.subplots()
                                title = sort_dirname+'_plus_minus_eight'
                                fig.suptitle(title)
                                ax.plot(b[:-1], h)
                                fig.savefig(os.path.join(save_fig_dir, 'plus_minus_eight', title))
                                plt.close('all')
                                print('here4')
                                eight_sts_settle_only = eight_sts[eight_sts<include_samples+start_st]
                                #for st in eight_sts:
                                #    if st< 8*30000:
                                #        eight_sts_settle_only.append(st)
                                print('here5')

                                h, b = np.histogram(eight_sts_settle_only/30000., bins=np.arange(np.max(eight_sts_settle_only)/30000.))
                                fig, ax = plt.subplots()
                                title = sort_dirname+'_settle_only'
                                fig.suptitle(title)
                                ax.plot(b[:-1], h)
                                save_path = os.path.join(save_fig_dir, 'settle_only', title)
                                print('saving_to: ', save_path)
                                fig.savefig(save_path)
                                plt.close('all')

                                raise(valueError)


                        except Exception as E:
                            #raise(E)
                            print('Error plotting', str(sort_path))
                            pass