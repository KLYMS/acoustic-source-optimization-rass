import numpy as np
import pandas as pd

df = pd.read_csv('/home/murali/Documents/rass/data/sim_data/cpp_sim/wTemp_wWind_cpp_alglib_CS_8.csv')

group_t = df.groupby('time')

for t in df['time'].unique():
    df_= group_t.get_group(t)
    df_.to_csv(f'/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_8/wTemp_wWind_cpp_alglib_CS_8_{t}.csv', index=False)

t_step = df['time'].unique()
np.savetxt('/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_8/time_step_cpp_alglib_CS_8.txt', t_step, fmt='%d')