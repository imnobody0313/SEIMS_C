import os
import os.path
import numpy as np
import pandas as pd
from pymongo import MongoClient
# import db_import_field_arrays

# subbasin_param_csv = './subbasin.csv'

# db_import_field_arrays.import_array_to_mongodb(spatial_gfs, mask_array, mask_name)
# 
conn = MongoClient('127.0.0.1', 27017)
db = conn.hlg_hband_longterm_model
is_permafrost = [0,0,1,1,1,0,0,0,1,1,0,1]
CH_BNK_K = [0.00001,0.001,0.001,0.001,0.001,0.00001,0.00001,0.001,0.001,0.001,0.00001,0.001]
GW_SPYLD = [0.1468,0.029,0.2936,0.2936,0.2936,0.1468,0.1468,0.029,0.2936,0.2936,0.1468,0.2936]

for i in range(1,13):
	db["REACHES"].update({'SUBBASINID': i}, {'$set': {"is_permafrost": is_permafrost[i-1]}},False,True)
	db["REACHES"].update({'SUBBASINID': i}, {'$set': {"GW_SPYLD": GW_SPYLD[i-1]}},False,True)
	db["REACHES"].find_one_and_update({'SUBBASINID': i},{'$set': {'CH_BNK_K': CH_BNK_K[i-1]}})#0.95
