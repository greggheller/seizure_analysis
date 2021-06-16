import csv
import pandas as pd
import os
import numpy as np

from matplotlib import pyplot as plt
from matplotlib import dates as mpldates

from datetime import datetime as dt
from statsmodels.stats.diagnostic import lilliefors
import seaborn as sns


base_path =r"C:\Users\svc_neuropix\Documents\python_scripts\siezure_analysis"
siezure_start_date = 20201202

analysis_start_date = 20190505
analysis_start_date_dt = dt.strptime(str(analysis_start_date), '%Y%m%d')

water_info_path = os.path.join(base_path, 'mouse_water_info_NP_rigs.csv') #r"\\10.128.50.20\sd7.2\siezure_analysis\mouse_water_info_NP_rigs.csv"

training_water_info_path = water_info_path = os.path.join(base_path, 'mouse_water_info_ALL_rigs.csv') 

mouse_lists = os.path.join(base_path, 'siezure_mice.csv')#r"\\10.128.50.20\sd7.2\siezure_analysis\siezure_mice.csv"

combined_stats_path = os.path.join(base_path, 'combined_mouse_stats.csv')#combined_stats_path = r"\\10.128.50.20\sd7.2\siezure_analysis\combined_mouse_stats.csv"

plots_path = os.path.join(base_path, 'plots')#r"\\10.128.50.20\sd7.2\siezure_analysis\plots"


date_columns = ['exp_date', 'Date of HP Surgery', 'Date of Birth', 'Date of WR']

date_diff_names = {
	'age_at_exp': ['exp_date', 'Date of Birth'],
	'age_at_surgery': ['Date of Birth', 'Date of HP Surgery'],
	'age_at_wr': ['Date of Birth', 'Date of WR'],
	'days_since_wr': ['exp_date', 'Date of WR'],
	'days_since_surgery': ['exp_date', 'Date of HP Surgery'],
	'days_of_recover_between_surgery_and_WR': ['Date of WR', 'Date of HP Surgery'],
}

def process_water_info():
	with open(water_info_path, 'r') as f:
		water_info = pd.read_csv(f)

	print(water_info.head())
	print(water_info.describe())
	twf = water_info.loc[:,'TWF']
	print(twf)
	mask = twf==0.85
	water_info = water_info.loc[mask, :]
	water_info['TotalWater'] = water_info.loc[:,'WS_ml'] + water_info.loc[:,'WE_ml']
	water_info['TotalWaterPercent'] = water_info.loc[:,'TotalWater']/water_info.loc[:,'BLW_g'] 
	print(water_info.head())
	print(water_info.describe())

	mouse_water_medians = water_info.groupby(["MID"]).median()
	print(mouse_water_medians.head())

	save_path = os.path.join(os.path.dirname(water_info_path), 'mouse_water_medians.csv')
	mouse_water_medians.to_csv(save_path)


def process_training_water_info():
	with open(training_water_info_path, 'r') as f:
		water_info = pd.read_csv(f)

	print(water_info.head())
	print(water_info.describe())
	#twf = water_info.loc[:,'TWF']
	#print(twf)
	#mask = twf==0.85
	#water_info = water_info.loc[mask, :]
	#water_info['TotalWater'] = water_info.loc[:,'WS_ml'] + water_info.loc[:,'WE_ml']
	#water_info['TotalWaterPercent'] = water_info.loc[:,'TotalWater']/water_info.loc[:,'BLW_g'] 
	#print(water_info.head())
	#print(water_info.describe())
	mouse_dates = water_info.loc[:, ('date_time','MID')]

	mouse_min_date = mouse_dates.groupby(["MID"]).min()
	mouse_min_date = mouse_min_date.rename(columns={'date_time':'Date of WR'})
	print(mouse_min_date.head())


	save_path = os.path.join(os.path.dirname(water_info_path), 'mouse_water_restrict_start_dates.csv')
	mouse_min_date.to_csv(save_path)

def get_exp_date(mouse_id, exp_dirs):
	mouse_dates = []
	for exp_dir in exp_dirs:
		for dirname in os.listdir(exp_dir):
			if str(mouse_id) in dirname and int(dirname.split('_')[1]) == int(mouse_id):
				mouse_dates.append(dirname)
	exp_session = min(mouse_dates)
	return get_date_from_session(exp_session)


def get_date_from_session(session):
	try:
		session_date = session.split('_')[2]
		#print(session_date)
		assert(len(session_date)==8)
	except Exception as E:
		session_date = None
	return session_date


def get_hand_off_date(mouse_id, hab_dirs):
	mouse_dates = []
	for exp_dir in hab_dirs:
		for dirname in os.listdir(exp_dir):
			if  str(mouse_id) in dirname and int(dirname.split('_')[1]) == int(mouse_id):
				mouse_dates.append(dirname)
	try:
		handoff_session = min(mouse_dates)
	except Exception as E:
		print('No habs found for '+str(mouse_id))
		handoff_session = None
	handoff_date = get_date_from_session(handoff_session)
	number_of_habs = len(mouse_dates)
	return handoff_date , number_of_habs


def gather_mouse_info_from_SDs():
	exp_paths = [r"\\10.128.50.20\sd7.2", r"\\10.128.50.20\sd7", r"\\10.128.50.43\sd6.3"]
	hab_paths = [os.path.join(path, 'habituation') for path in exp_paths]
	with open(mouse_lists, 'r') as f:
		mouse_info = pd.read_csv(f)

	columns = ['MID', 'exp_date', 'handoff_date', 'number_of_habs']
	sd_info = pd.DataFrame(columns=columns)	

	analyze_mice = mouse_info.loc[:, 'all mice']
	for mouse_id in analyze_mice:
		exp_date = get_exp_date(mouse_id, exp_paths)
		handoff_date, number_of_habs = get_hand_off_date(mouse_id, hab_paths)
		sd_info.loc[len(sd_info.index)] = [mouse_id, exp_date, handoff_date, number_of_habs] 

	save_path = os.path.join(os.path.dirname(water_info_path), 'mouse_sd_info.csv')
	sd_info.to_csv(save_path)


def combine_mouse_stats():

	file_names = ['mouse_sd_info.csv', 'mouse_water_medians.csv', 'Brain_Compression_Log.csv', 'surgeon_info.csv', 'mouse_water_restrict_start_dates.csv']

	#load the first csf as a DF
	filename1 = file_names.pop(0)
	first_path = os.path.join(os.path.dirname(water_info_path), filename1)
	with open(first_path, 'r') as f:
		combined_df = pd.read_csv(f)
	combined_df['MID']=combined_df['MID'].astype(int)

	#add all the other csvs to the end of the first
	print(combined_df.head())
	for filename in file_names:
		this_path = os.path.join(os.path.dirname(water_info_path), filename)
		with open(this_path, 'r') as f:
			this_df = pd.read_csv(f)
		this_df['MID']=this_df['MID'].astype(int)
		print(this_df.head())
		combined_df = combined_df.merge(this_df, how='left',  left_on='MID', right_on='MID')

	#drop all the "unnamed" 
	for column in combined_df.columns:
		#print(column)
		if 'unnamed' in column.lower():
			combined_df = combined_df.drop(column, 1)


	#convert date clumns to dt and plt
	for date_column in date_columns:
		date_str_format = '%m/%d/%Y'
		if date_column == 'exp_date':
			date_str_format = '%Y%m%d'
		if date_column == 'Date of WR':
			date_str_format = '%m/%d/%Y %H:%M'
		combined_df[date_column+'_dt'] =  pd.to_datetime(combined_df[date_column], format=date_str_format)
		combined_df[date_column+'_plt'] = mpldates.date2num(combined_df[date_column+'_dt'])


	#add columns for the time between important dates
	for name, pair in date_diff_names.items():
		col_1 = pair[0]+'_dt'
		col_2 = pair[1]+'_dt'
		combined_df[name] = abs(combined_df[col_1].sub(combined_df[col_2]).dt.days)

	#add the siezure column as bools
	with open(mouse_lists, 'r') as f:
		mouse_info = pd.read_csv(f)

	siezure_mice = mouse_info['siezure mice'].fillna(0).astype(int)
	#siezure_mice = siezure_mice_int.astype(str)
	#print(siezure_mice.dtype)
	#print(siezure_mice)
	#print(combined_df['MID'].dtype)
	#print('testing signle num 546512: ' , (546512 in list(siezure_mice)))
	combined_df['siezure'] = [value in list(siezure_mice) for value in combined_df['MID']]
	combined_df['siezure'] = combined_df['siezure'].astype(int)
	#print('SIEZURE MICE')
	#print('########################')
	#print(combined_df['siezure'])



	#convert thes strings to ints
	psuedo_numerical_cols = ['Edema (Swelling)', 'Hematoma', 'Sinus Bleed', 'Laceration']
	mapping = {
		'None': 0,
		'Mild': 1,
		'Moderate': 2,
		'Severe': 3
	}
	for column in psuedo_numerical_cols:
		num_values = []
		for value in combined_df[column]:
			num_value = np.nan
			try:
				num_value = mapping[value]
			except Exception as E:
				pass
			num_values.append(num_value)
		num_name = column+'_num'
		combined_df[num_name] = num_values


	#add summary statistics for surgical damage
	surg_num_cols = []
	for column in psuedo_numerical_cols:
		surg_num_cols.append(column+'_num')
	combined_df['surgery_damage_average'] = combined_df.loc[:, surg_num_cols].mean(axis=1)
	combined_df['surgery_damage_sum'] = combined_df.loc[:, surg_num_cols].sum(axis=1)
	combined_df['surgery_damage_max'] = combined_df.loc[:, surg_num_cols].max(axis=1)



	print(combined_df.head())
	save_path = os.path.join(os.path.dirname(water_info_path), 'combined_mouse_stats.csv')
	combined_df.to_csv(save_path)




def plots_over_time(combined_df, column_names, column_values, category_name=False):
	if category_name is False:
		category_name = ', '.join(column_names)
	for date_column in date_columns:
		fig, ax = plt.subplots()
		column_name = str(column_names)
		#print(combined_df.loc[:, 'exp_date'])
		try:
			date_column_plt = date_column+'_plt'
			date_column_dt = date_column+'_dt'
			combined_df[date_column_plt] = mpldates.date2num(combined_df[date_column_dt])
				#print(combined_df['exp_date_plt'])
			#raise(ValuError)
			#print(combined_df.loc[:, date_column_plt])
			#print(column_values)
			plt.plot_date(combined_df.loc[:, date_column_plt], column_values)
			plt.legend(column_names)
			ax.set_ylabel(category_name)
			ax.set_xlabel(date_column)
			ax.set_title(category_name+' over time ('+date_column+')')
			plt.xticks(rotation=30)

			save_path = os.path.join(plots_path, category_name+'_'+date_column+'.png')

			fig.savefig(save_path)
			

			days_to_smooth_list = [30]
			for days_to_smooth in days_to_smooth_list:
				fig, ax = plt.subplots()
				smoothed_column = []
				for date_in in combined_df.loc[:, date_column_dt]:
					in_range = []
					in_range_dates = []
					for idx, date_in2 in enumerate(combined_df.loc[:, date_column_dt]):
						diff_time = date_in2 - date_in
						diff_days = abs(diff_time.days)
						if diff_days < days_to_smooth:
							in_range.append(column_values[idx])
							in_range_dates.append(combined_df.loc[:, date_column_dt][idx])
					#print(in_range_dates)
					#print(in_range)
					smoothed_value = np.nanmean(in_range)
					#print(smoothed_value)
					smoothed_column.append(smoothed_value)
				#print('These are the column with values', str(len(smoothed_column)), smoothed_column)
				#(unique, counts) = numpy.unique(smoothed_column, return_counts=True)
				#print('num of smoothed_values', counts)
				#(unique, counts) = numpy.unique(combined_df.loc[:, 'exp_date_plt'], return_counts=True)
				#print('num of smoothed_dates', counts)
				#print('These are the dates', str(len(combined_df.loc[:, 'exp_date_plt'])), combined_df.loc[:, 'exp_date_plt'])
				plt.plot_date(combined_df.loc[:, date_column_plt], smoothed_column)
				plt.legend(column_names)
				ax.set_ylabel(str(days_to_smooth)+' Day Rolling average of '+category_name)
				ax.set_xlabel(date_column)
				ax.set_title(category_name+' (smoothed over '+date_column+' days)')
				plt.xticks(rotation=30)
				save_path = os.path.join(plots_path, category_name+'_time_smoothed_'+str(days_to_smooth)+'_'+date_column+'.png')
				fig.savefig(save_path)
		except Exception as E:
			print('failed to plot time for ', column_name)
			#raise(E)


def plot_siezure_plots_categorical(combined_df, column_name):

	column_names = combined_df[column_name].unique()
	column_values = pd.get_dummies(combined_df, columns=[column_name])
	#print(column_values.head())
	#print(column_values.columns)
	dummy_names = [column for column in column_values.columns if column_name+'_' in column]
	#print(dummy_names)
	column_values = column_values.loc[:, dummy_names]
	#print(column_values.head())
	#raise(ValueError)
	plots_over_time(combined_df, column_names, column_values, category_name=column_name)

	siezure_start_date_dt = dt.strptime(str(siezure_start_date), '%Y%m%d')







def plot_siezure_plots_numeric(combined_df, column_name):

	column_values = combined_df[column_name]
	plots_over_time(combined_df, [column_name], column_values)

	siezure_start_date_dt = dt.strptime(str(siezure_start_date), '%Y%m%d')
	before_siezure = []
	after_siezure = []
	#print(len(combined_df.loc[:, 'exp_date_dt']))
	#print(len(combined_df.loc[:, column_name]))
	#print(combined_df.loc[:, column_name])
	for idx, date_in in enumerate(combined_df.loc[:, 'exp_date_dt']):
		#print(idx)
		#print(date_in)
		#print(siezure_start_date_dt)
		diff_time = siezure_start_date_dt - date_in
		diff_days = diff_time.days
		value = combined_df.loc[:, column_name][idx]
		if diff_days > 0:
			before_siezure.append(value)
		else:
			after_siezure.append(value)


	def clean_list(in_list):
		#print('UNCLEAN: ', in_list)
		in_list = np.array(in_list)
		in_list = in_list[~pd.isnull(in_list)]
		#in_list = in_list[pd.isnull(in_list)]
		#for idx, element in enumerate(in_list):
		#	if not(element is None):
	#			remove_idxs.append(idx)
	#		else:
	#			try:
	#				if np.isnan(element):
	#					print('Found nan')
	#					in_list.remove(element)
	#			except Exception as E:
	#				pass
			#cleanedList = [x for x in in_list if str(x) != 'nan']
		#print('CLEAN: ', in_list)
		return in_list

	before_siezure = clean_list(before_siezure)
	after_siezure = clean_list(after_siezure)

	n, bins, patches = plt.hist(combined_df.loc[:, column_name])
	fig, ax = plt.subplots()
	plt.hist(after_siezure, bins, weights=np.ones(len(after_siezure)) / len(after_siezure), alpha=0.5, label='after DEC')
	plt.hist(before_siezure, bins, weights=np.ones(len(before_siezure)) / len(before_siezure), alpha=0.5, label='before DEC')
	plt.legend(loc='upper right')
	save_path = os.path.join(plots_path, column_name+'_date_hist.png')
	ax.set_ylabel('Percentage of date group (before DEC/after DEC)')
	ax.set_xlabel('Categories of: '+column_name)
	ax.set_title('Percentage of date group for '+column_name+' categories')
	plt.xticks(rotation=30)

	fig.savefig(save_path)


	with open(mouse_lists, 'r') as f:
		mouse_info = pd.read_csv(f)

	#print(mouse_info.head())

	#analyze_mice = mouse_info.loc[:, 'all mice']
	#print('Mice to Analyze:')
	#print(analyze_mice)

	siezure_mice = mouse_info['siezure mice']
	#print('Mice with Siezures:')
	#print(siezure_mice)
	siezure = []
	no_siezure = []
	for idx, MID in enumerate(combined_df.loc[:, 'MID']):
		#print(MID)
		#print(float(MID))
		value = combined_df.loc[:, column_name][idx]
		if float(MID) in list(siezure_mice):
			siezure.append(value)
		else:
			no_siezure.append(value)
	siezure = clean_list(siezure)
	no_siezure = clean_list(no_siezure)


	#print(no_siezure)
	n, bins, patches = plt.hist(list(siezure)+list(no_siezure))

	fig, ax = plt.subplots()
	siezure_n, bins2, patches = plt.hist(siezure, bins, weights=np.ones(len(siezure)) / len(siezure), alpha=0.5, label='siezure')
	no_siezure_n, bins2, patches = plt.hist(no_siezure, bins, weights=np.ones(len(no_siezure)) / len(no_siezure), alpha=0.5, label='no_siezure')
	plt.legend(loc='upper right')
	ax.set_ylabel('Percentage of siezure group (siezure/non-sizure)')
	ax.set_xlabel('Categories of: '+column_name)
	ax.set_title('Percentage of siezure group for '+column_name+' categories')
	plt.xticks(rotation=30)


	save_path = os.path.join(plots_path, column_name+'_siez_hist.png')
	fig.savefig(save_path)
	n, bins, patches = plt.hist(list(siezure)+list(no_siezure))

	fig, ax = plt.subplots()


	labels=[item for item in ax.get_xticklabels()]
	fig, ax = plt.subplots()

	siezure_n, bins2, patches = plt.hist(siezure, bins, alpha=0.5, label='siezure')
	no_siezure_n, bins2, patches = plt.hist(no_siezure, bins, alpha=0.5, label='no_siezure')

	#print('siez_n: ',siezure_n)
	#print('siez_n: ',no_siezure_n)
	#print('n: ',n)
	fraction_siez = siezure_n/n
	fraction_siez[fraction_siez==np.inf]=0
	fraction_siez[pd.isnull(fraction_siez)] = 0
	#print(fraction_siez)
	fraction_non = no_siezure_n/n
	fraction_non[pd.isnull(fraction_non)] = 0
	#print(fraction_non)
	#print(bins)

	
	fig, ax = plt.subplots()

	labels = []
	for idx in range(len(bins)-1):
		labels.append(np.round((bins[idx]+bins[idx+1])/2, 3))

	bar_placements = list(range(len(labels)))
	#print(labels)
	ax.bar(bar_placements, fraction_siez, width=.8, label='siezure')
	ax.bar(bar_placements, fraction_non, width=.8, bottom=fraction_siez,
	       label='no siezure')

	ax.set_xticks(bar_placements)
	ax.set_xticklabels(labels)

	ax.set_ylabel('Fraction of category with and without siezure')
	ax.set_xlabel('Categories of: '+column_name)
	ax.set_title('Fraction of siezures for '+column_name+' categories')
	ax.legend()
	plt.xticks(rotation=30)
	#fig = plt.figure()
	#plt.hist([siezure, no_siezure], bins, weights=np.ones(len(siezure)) / len(siezure), alpha=0.5, label='siezure')
	#plt.hist(, bins weights=np.ones(len(no_siezure)+len(no_siezure)) / len(no_siezure), alpha=0.5, label='no_siezure')
	#plt.legend(loc='upper right')
	save_path = os.path.join(plots_path, column_name+'_siez_stacked_bar.png')
	fig.savefig(save_path)

	#make a second one that splits it based on date - two subplots from both before and after date shift
	#gonna need to do something different if categorical anyway...


def load_df_with_dates(path):
	with open(path, 'r') as f:
		combined_stats = pd.read_csv(f)

	for column in combined_stats.columns:
		if '_dt' in column:
			col_vals = combined_stats[column]
			dates = [None]*len(col_vals)
			for idx, date_in in enumerate(col_vals):
				#print(date_in)
				try:
					dates[idx] =  dt.strptime(str(date_in), '%Y-%m-%d %H:%M:%S')
				except Exception as E:
					pass
			combined_stats[column] = dates

	return combined_stats


def produce_all_plots(combined_stats):

	for column in ['HP Weight before', 'age_at_exp', 'age_at_wr', 'age_at_surgery']:
		"""['siezure', 'TotalWaterPercent', 'BLW_g', 'number_of_habs', 'Compression ',  'Health ',
		'Edema (Swelling)_num', 'Hematoma_num', 'Sinus Bleed_num', 'Laceration_num', 
		'surgery_damage_average', 'surgery_damage_sum', 'surgery_damage_max', 'HP Work Station', 
		'HP Iso Level', 'HP Recovery', 'HP Iso Duration', 'Sex', 'LabTracks Group', 'HP Surgeon',
		'age_at_surgery', 'age_at_exp', 'days_since_surgery', #'np_rig', 'np_operator', 
		'days_since_wr', 'age_at_wr', 'days_of_recover_between_surgery_and_WR', ]:"""
		print('producing plots for '+ column)

		if pd.api.types.is_numeric_dtype(combined_stats[column]):
			print('numeric')
			plot_siezure_plots_numeric(combined_stats, column)
			
		else:
			print('categorical')
			plot_siezure_plots_categorical(combined_stats, column)
			


if __name__ == '__main__':
	#process_water_info()

	#process_training_water_info()

	#gather_mouse_info_from_SDs()

	combine_mouse_stats()

	combined_stats = load_df_with_dates(combined_stats_path)
	#print(combined_stats['exp_date_dt'])

	pre_covid_surgery = combined_stats[(combined_stats['Date of HP Surgery_dt'] < analysis_start_date_dt)]
	print('How many mice had siezures that recieved surgery before covid?')
	print(pre_covid_surgery['siezure'].sum())
	print('here is when we ran these expeirments:')
	print(pre_covid_surgery[['MID', 'exp_date', 'age_at_exp']])
	print('here are the ages:')
	print(pre_covid_surgery['age_at_exp'])

	combined_stats = combined_stats[(combined_stats['Date of HP Surgery_dt'] > analysis_start_date_dt)]
	save_path = os.path.join(os.path.dirname(water_info_path), 'filtered_combined_stats.csv')
	combined_stats.to_csv(save_path)

	#testing poinson - 
	#first order by date and grab the siezure column
	sorted_stats = combined_stats.loc[:,('MID', 'exp_date','siezure')]
	sorted_stats = sorted_stats.sort_values(by='exp_date')
	print(sorted_stats.head(50))
	#get indexes of siezures
	ordered_siezures = sorted_stats['siezure'].to_list()
	print(ordered_siezures)
	indexes = np.where(sorted_stats['siezure']==1)[0]
	print(indexes)
	#add 0 and the end - decided not to do this
	#subtract to get inter-event intervals
	intervals = np.abs(indexes[0:len(indexes)-1] - indexes[1:len(indexes)])
	print(intervals)

	#plot the distribution
	fig, ax = plt.subplots()
	plt.hist(intervals, bins= [1,2,3,4,5,6,7,8], alpha = .5)
	plt.legend(loc='upper right')
	dist_save_path = os.path.join(plots_path, 'dist_inter_ziezure_intervals.png')
	swapped = [1, 2, 4, 6, 2, 7, 1, 3, 2, 1, 1, 1, 1, 2]
	plt.hist(swapped, bins= [1,2,3,4,5,6,7,8], alpha = .5)
	ax.set_ylabel('Number of observances')
	ax.set_xlabel('inte siezure interval')
	ax.set_title('distrubution of inter-ziezure intervals')
	plt.xticks(rotation=30)
	fig.savefig(dist_save_path)

	ksstat, pvalue = lilliefors(intervals, dist='exp')
	print('test exp For intervals: ', ksstat, pvalue)
	ksstat, pvalue = lilliefors(swapped, dist='exp')
	print('test exp For swapped: ', ksstat, pvalue)

	ksstat, pvalue = lilliefors(intervals, dist='exp')
	print('test norm For intervals: ', ksstat, pvalue)
	ksstat, pvalue = lilliefors(swapped, dist='exp')
	print('test norm For swapped: ', ksstat, pvalue)

	for i in range(10):
		test = np.random.exponential(size=len(swapped))
		ksstat, pvalue = lilliefors(test, dist='exp')
		print('For test ', i, ': ', ksstat, pvalue)

	combined_stats = load_df_with_dates(save_path)
	combined_stats['Genotype'] = combined_stats['LabTracks Group']
	combined_stats['Genotype'] = combined_stats['Genotype'].replace('C57BL6J', 'C57BL6J')
	combined_stats['Genotype'] = combined_stats['Genotype'].replace('C57BL6J(NP)', 'C57BL6J')
	combined_stats['Genotype'] = combined_stats['Genotype'].replace('C57BL6J (NP)', 'C57BL6J')
	print('GENOTYPE#############3')
	print(combined_stats['Genotype'])
	combined_stats['Genotype_Sex'] = combined_stats['Genotype']+'_'+combined_stats['Sex']
	#print(combined_stats['exp_date_dt'])

	fig, ax = plt.subplots()

	dist_save_path = os.path.join(plots_path, 'baseline_hist_genotype_sex.png')

	#sns.scatterplot(x='age_at_wr', y='BLW_g', hue = 'LabTracks Group', data=combined_stats) #hue='label'Pedigree Name
	sns.histplot(x='BLW_g', hue = 'Genotype_Sex', data=combined_stats, multiple="stack")#LabTracks Group #stat='percent'?
	ax.set_title('Baseline weight vs age at water restriction')
	plt.xticks(rotation=30)
	fig.savefig(dist_save_path)
	produce_all_plots(combined_stats)


	