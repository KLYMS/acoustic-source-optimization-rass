import pandas as pd
import numpy as np
import warnings
import argparse
import re
warnings.filterwarnings('ignore') # ignore warnings in output 

output_filepath = '/home/murali/Documents/rass/automation/gen_data/'

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('file_path', type=str, help='Path to the CSV file')
args = parser.parse_args()
file_path = args.file_path

pattern_filename = r'/([^/]+)\.CSV$'
match_filename = re.search(pattern_filename, file_path)

if match_filename:
    file_name = match_filename.group(1)
    # print(file_name)
else:
    print("No match found")

df = pd.read_csv(file_path, skiprows = 6)
pd.set_option('display.max_columns', None)
df.columns = df.columns.str.strip()

wanted_cols = ['ObsTime', 'Height', 'WS', 'WD', 'Temp0']
data = df[wanted_cols]

# time 
data['Time'] = pd.to_datetime(data['ObsTime'])
data.drop(['ObsTime'], axis=1, inplace = True)

data['timedelta'] = data['Time'] - data.loc[0,'Time']
data['seconds'] = data['timedelta'].dt.total_seconds()
data.drop(['timedelta'],axis=1,inplace=True)

# cleaning after the max height and if less than 0
max_height_time = data.loc[data['Height'].idxmax(), 'Time']
data = data[data.sort_values(by = 'Time')['Time'] <=max_height_time]
idx = data.loc[data['Height']<0].index
data.drop(idx, inplace=True)

# cleaning the data
data.loc[data['Temp0'] == '-----','Temp0'] = np.nan
data['Temp0'] = data['Temp0'].astype(float)
data['Temp'] = data['Temp0'] + 273
data.drop(['Temp0'],axis=1, inplace=True)
data = data[~data['Temp'].isna()]

bin_size = 100

data['Height_bin'] = (data['Height'] // bin_size) * bin_size
grouped = data.groupby('Height_bin').mean().reset_index()
grouped['WD_'] = (grouped['WD'] + 180) % 360

saving_data = grouped[['Height','Temp','WS','WD']]
saving_data.to_csv( output_filepath + file_name + '_cleaned.csv', index=False)