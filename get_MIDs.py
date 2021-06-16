import pandas as pd
import os

base_path =r"C:\Users\svc_neuropix\Documents\python_scripts\siezure_analysis"
seizure_list_path = os.path.join(base_path, 'seizure_dataframe_20210527.csv')#r"\\10.128.50.20\sd7.2\siezure_analysis\siezure_mice.csv"

with open(seizure_list_path, 'r') as f:
	seizure_list = pd.read_csv(f)


MIDs = set(seizure_list['mid'])

for MID in MIDs:
		print(MID)