import pandas as pd
import os
from I_OpenAQ_Share import find_root_dir

# Define the file path
debug_mode = 0;
direc = find_root_dir(debug_mode);

# output file
out_file_max = f"{direc}/Analysis_Data/Ion_Chromatography/IC_MaxConc_By_Site.csv";
out_file_perc = f"{direc}/Analysis_Data/Ion_Chromatography/IC_99percentile_By_Site.csv";
ic_max = []
ic_perct = []

# read site sampling data
site_info = pd.read_excel(f"{direc}/Site_Sampling/Site_details.xlsx")
# loop througth site we have to get the master file and read the max conc of each IC species
for index, row in site_info.iterrows():
    site_code = row['Site_Code']

    master_file_name = f"{direc}/Analysis_Data/Master_files/{site_code}_master.csv"
    if os.path.exists(master_file_name):
        df = pd.read_csv(master_file_name,na_values=['NaN',' NaN','  NaN'])
        master = df.fillna(0) # replace nan with 0
        
        if index == 0:
            # Find columns that start with 'IC_'
            ic_columns = [col for col in master.columns if col.startswith('IC_')]
        
        # Calculate the maximum values 
        max_values = master[ic_columns].max().tolist()
        # Calculate the 90th percentile values 
        percentile_values = master[ic_columns].quantile(0.99).tolist()
        
        # Append the maximum values to the output list
        ic_max.append([site_code] + max_values)
        ic_perct.append([site_code] + percentile_values)

output_columns = ['Site Code'] + ic_columns

# Convert the output to a DataFrame and output as csv
output_df = pd.DataFrame(ic_max, columns=output_columns)
output_df.to_csv(out_file_max, index=False)

# Convert the output to a DataFrame and output as csv
output_df = pd.DataFrame(ic_perct, columns=output_columns)
output_df.to_csv(out_file_perc, index=False)

print('done')

