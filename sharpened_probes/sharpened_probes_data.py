
import pandas as pd
import numpy as np
from collections import defaultdict
import os

gregg_path = r"C:\Users\greggh\Documents\python_scripts\sharpened_probes\gregg_sharpened.csv"
tamina_path = r"C:\Users\greggh\Documents\python_scripts\sharpened_probes\tamina_sharpened.csv"

path_list = [gregg_path, tamina_path]
#sharpened_sessions_path = r"C:\Users\greggh\Documents\python_scripts\sharpened_probes\sharpened_sessions.csv"



#first check the data
pd.options.display.max_columns = None
pd.options.display.max_rows = None

def visual_check_path(csv_path):
	with open(csv_path, 'r') as f:
		df = pd.read_csv(f)
		visual_check(df)

def visual_check(df):
	print(df.columns)
	print(df.head())
	print(df.describe(include='all'))

def visual_checks(dfs):
	for df in dfs:
		visual_check(df)



#first combine them - try it with both mean and max - make sure it doens't treat absences (only one person did it ) as a 0
#open both CSVs
dfs = []
for path in path_list:
	with open(path, 'r') as f:
		df = pd.read_csv(f)
		dfs.append(df)

for idx, df in enumerate(dfs):
	for column_name in df.columns:
		#print(column_name)
		if 'ï»¿' in column_name:
			new_name = column_name[3:]
			print("renaming ", column_name, " to ", new_name )
			df = df.rename(columns={column_name: new_name})
			#print(df.columns)
		dfs[idx] = df

#visual_checks(dfs)

#grab all the session IDS - use set to make unique
session_IDs = set()
for df in dfs:
	these_ids = df['MID']
	#print(these_ids)
	session_IDs = session_IDs.union(set(these_ids))

print(session_IDs)

#Then fill in new dataframes
probes = 'ABCDEF'
def probe_name(probe):
	return 'probe'+probe

def blood_name(probe):
	return 'probe'+probe+'_blood'

def sharp_name(probe):
	return 'probe '+probe+' sharpened'

def unsharp_name(probe):
	return 'probe '+probe+' unsharpened'	

def mean_name(probe):
	return blood_name(probe)+'_mean'

def max_name(probe):
	return blood_name(probe)+'_max'

print()


columns = ['MID']
for probe in probes:
	columns.append(mean_name(probe))
	columns.append(max_name(probe))
	columns.append(sharp_name(probe))

combined_df = pd.DataFrame(columns= columns)

for MID in session_IDs:
	if not(np.nan is MID):
		#print(MID)
		row_dict = {}
		row_dict['MID'] = MID
		bleed_values = defaultdict(list)
		sharp_values = defaultdict(list)

		for df in dfs:
			idxs = None
			try:
				idxs = np.where(df["MID"] == MID)[0]
				#print(idxs)
			except Exception as E:
				pass
			else:
				for idx in idxs:
					row = df.iloc[idx, :]
					for probe in probes:
						#print(probe)
						#print(row.loc['Sharpened or Unsharpened'])
						#print(int(row.loc['Sharpened or Unsharpened'])
						#print(row.loc[probe_name(probe)])
						#print(row)
						try:
							bleed_values[probe].append(int(row.loc[probe_name(probe)]))
							#print('int value', int(row.loc[probe_name(probe)]))
							sharp_values[probe].append(int(row.loc['Sharpened or Unsharpened']))
						except Exception as E:
							pass
		for probe in probes:
			if sharp_values[probe]:
				#print(sharp_values)
				#print(bleed_values)
				#raise(ValueError)
				sharpened = np.median(sharp_values[probe])
				assert(sharpened == sharp_values[probe][0])
				row_dict[mean_name(probe)] = np.mean(bleed_values[probe])
				row_dict[max_name(probe)] = np.max(bleed_values[probe])
				row_dict[sharp_name(probe)] = sharpened
		combined_df = combined_df.append(row_dict, ignore_index=True)

visual_check(combined_df)

save_path = r"C:\Users\greggh\Documents\python_scripts\sharpened_probes\combined.csv"
combined_df.to_csv(save_path)



#then do the counts for all 6 probes, and then for all combined, sharpened and unsharpened
categories = ['None', 'Miniscule', 'Mild', 'Moderate', 'Severe']
bins = [0,.1,1.5,2.5,3.5,6]
columns = ['Description', 'N_insertions'] 
columns.extend(categories)

max_counts_df = pd.DataFrame(columns= columns)
mean_counts_df = pd.DataFrame(columns= columns)

df_dict = {
	'max': max_counts_df,
	'mean': mean_counts_df
}
func_dict = {
	'max': max_name,
	'mean': mean_name
}


for key, df in df_dict.items():
	total_sharp = {'Description': 'All Probes Sharpened'}
	total_unsharp = {'Description': 'All Probes Unsharpened'}
	for idx, category in enumerate(categories):
		total_sharp[category] = 0
		total_unsharp[category] = 0
	total_count_sharps = 0
	total_count_unsharps = 0

	for probe in probes:
		sharp_status = combined_df[sharp_name(probe)].astype(bool)
		value_name = func_dict[key](probe)
		sharp_values = combined_df[value_name][sharp_status]
		unsharp_values = combined_df[value_name][~sharp_status]
		#print(sharp_values)
		#print(unsharp_values)
		sharp_dict = {}
		unsharp_dict = {}
		sharp_dict['Description'] = sharp_name(probe)
		unsharp_dict['Description'] = unsharp_name(probe)
		sharp_value_counts = sharp_values.value_counts(bins=bins, sort=False) 
		total_probe_sharps = sum(sharp_value_counts)
		sharp_dict['N_insertions'] = total_probe_sharps
		total_count_sharps+= total_probe_sharps
		#print(sharp_value_counts)
		#print(total_sharps)
		unsharp_value_counts = unsharp_values.value_counts(bins=bins, sort=False) 
		total_probe_unsharps = sum(unsharp_value_counts)
		unsharp_dict['N_insertions'] = total_probe_unsharps
		total_count_unsharps+= total_probe_unsharps
		for idx, category in enumerate(categories):
			sharp_dict[category] = sharp_value_counts[bins[idx+1]]/total_probe_sharps
			total_sharp[category] = sharp_value_counts[bins[idx+1]]+ total_sharp[category]
			unsharp_dict[category] = unsharp_value_counts[bins[idx+1]]/total_probe_unsharps
			total_unsharp[category] = unsharp_value_counts[bins[idx+1]] + total_unsharp[category]
		print('########################################')
		print(sharp_dict)
		df = df.append(sharp_dict, ignore_index=True)
		df = df.append(unsharp_dict, ignore_index=True)
		print(df.head(5))

	total_sharp['N_insertions'] = total_count_sharps
	total_unsharp['N_insertions'] = total_count_unsharps
	for category in categories:
		total_sharp[category] = total_sharp[category]/total_count_sharps
		total_unsharp[category] = total_unsharp[category]/total_count_unsharps

	df = df.append(total_sharp, ignore_index=True)
	df = df.append(total_unsharp, ignore_index=True)
	df_dict[key] = df

#just save the counts in a sperate CSV (or2) so the plotting doesn't work with the data at all
for key, df in df_dict.items():
	save_path = os.path.join(r"C:\Users\greggh\Documents\python_scripts\sharpened_probes", key+"_fraction.csv")
	df.to_csv(save_path)



print('finished')
	






