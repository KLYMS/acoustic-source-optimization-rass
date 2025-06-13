import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
df_0 = pd.read_csv('/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_0/backscatter_coordinates_0_G250-r2.csv')

df_1 = pd.read_csv('/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_1/backscatter_coordinates_1_G250-r2.csv')
df_2 = pd.read_csv('/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_2/backscatter_coordinates_2_G250-r2.csv')
df_3 = pd.read_csv('/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_3/backscatter_coordinates_3_G250-r2.csv')
df_4 = pd.read_csv('/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_4/backscatter_coordinates_4_G250-r2.csv')

df_5 = pd.read_csv('/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_5/backscatter_coordinates_5_G250-r2.csv')
df_6 = pd.read_csv('/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_6/backscatter_coordinates_6_G250-r2.csv')
df_7 = pd.read_csv('/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_7/backscatter_coordinates_7_G250-r2.csv')
df_8 = pd.read_csv('/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_8/backscatter_coordinates_8_G250-r2.csv')

df =  pd.concat([df_0,df_1,df_2,df_3,df_4,df_5,df_6,df_7,df_8])

plt.figure(figsize=(10, 6))
plt.scatter(df['x_val'],df['z_val'], s=0.5)

dist = 35 * 1000
plt.xlim((-dist,dist))
plt.ylim((0,dist))
plt.title('N-S')
plt.xlabel("distance")
plt.ylabel("altitude")
plt.grid()
plt.show()

plt.figure(figsize=(10, 6))
plt.scatter(df['y_val'],df['z_val'], s=0.5)

dist = 35 * 1000
plt.xlim((-dist,dist))
plt.ylim((0,dist))
plt.title('E-W')
plt.xlabel("distance")
plt.ylabel("altitude")

plt.grid()
plt.show()