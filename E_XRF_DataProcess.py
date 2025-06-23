# PURPOSE: ------------------------------------------------------------
# Read and perform QA/QC on XRF data straight from the XRF
# system with no pre-processing; put XRF data and relevant flags in the
# Master data file for each SPARTAN site that has data contained in the XRF
# file; move XRF files from "Raw" to relevant site folder
#
# FILE REQUIREMENTS: ------------------------------------------------------
# Files in the raw file directory should be straight from the XRF system, unmodified
# 
# Input: -----------------------------------------------------------
# Analysis_Data/Filter_Masses/NEW_MTL/$SiteCode_CartID_{pre/post}_weight.xlsx
# 
# Output: ----------------------------------------------------------
# Analysis_Data/Filter_Masses/Masses_by_site/$SiteCode_masses.xlsx
#
# Suppliment: ------------------------------------------------------
# Acceptancrelease_testing_Summary_WashU.xlsx
# Analysis_Data/E-Logs/XRF E-Log.xlsx
# 
# To re-process a raw MTL file: ------------------------------------------- 
# Turn Redo_All_Archieve to 1
# 
# Created by: --------------------------------------------------------------
# Haihui Zhu, Nidhi Anchan, Crystal Weagle, Emme Le Roy, and Chris Oxford 
# May-2024

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import pandas as pd
import os
import math
from datetime import datetime
import logging
import spt_utils as su

# ============================
# Initial Setup and Configuration
# ============================
Redo_All_Archive = 0
debug_mode = 0
direc = su.find_root_dir(debug_mode)

# ============================
# Directory & File Paths
# ============================
direc_new = os.path.join(direc, 'Analysis_Data', 'XRF', 'NEW')
direc_archive = os.path.join(direc, 'Analysis_Data', 'Archived_Filter_Data', 'XRF_data')
direc_master = os.path.join(direc, 'Analysis_Data', 'Master_files')
log_file_path = os.path.join(direc, 'Public_Data/Data_Processing_Records/XRF', f"{datetime.now().strftime('%Y-%m-%d-%H%M%S')}_XRF_Record.txt")
elogfname = f'{direc}/Analysis_Data/E-Logs/XRF E-Log.xlsx' 

# TEST 
# direc_new = os.path.join('./release_test', 'XRF_NEW')
# direc_archive = os.path.join('./release_test','XRF_archive')
# direc_master = os.path.join('./release_test', 'Master_files')
# log_file_path = os.path.join('./release_test', f"{datetime.now().strftime('%Y-%m-%d-%H%M%S')}_XRF_Record.txt")
# elogfname = f'./release_test/XRF E-Log.xlsx' 

# Log file setup
su.setup_logging(log_file_path)  # Initialize logging
logging.info("E started")
logging.info(f'{datetime.now()} \n')

# Read the site info table
SiteInfoFile = f'{direc}/Site_Sampling/Site_details.xlsx'
site_details = pd.read_excel(SiteInfoFile)


# ============================
# Functions & Constants
# ============================
r = (21.2 / 10) / 2  # diameter is 21.2 mm, need to convert to cm and then to radius
filter_area = math.pi * r ** 2

# Column indices for metal data (element of interest is in brackets)
xrf_cols = ['C (Na)', 'C (Al)', 'C (Si)', 'C (S)', 'C (Cl)', 'C (K)', 'C (Ca)', 'C (Ti)', 
          'C (V)', 'C (Fe)', 'C (Zn)', 'C (Ce)', 'C (Pb)', 'C (As)', 'C (Co)', 
          'C (Cr)', 'C (Cu)', 'C (Mg)', 'C (Mn)', 'C (Ni)', 'C (Sb)', 'C (Rb)', 
          'C (Sr)', 'C (Cd)', 'C (Se)', 'C (Sn)']
metals = []
for col in xrf_cols:
    start_idx = col.find('(') + 1  # Start after '('
    end_idx = col.find(')')        # End at ')'
    metals.append(col[start_idx:end_idx].strip())  # Extract and clean the element symbol
    
# ============================
# Start processing 
# ============================
# move archived files to new (if enabled)
if Redo_All_Archive:
    logging.info(f'Copying archived XRF report for reprocessing')
    # step 1: read the list of files that do not need to be reprocessed
    skipfile_path = os.path.join(direc_archive, 'No_need_to_reprocess.txt')
    if os.path.exists(skipfile_path):
        with open(skipfile_path, 'r') as f:
            skip_list = [line.strip() for line in f if line.strip()]
    else:
        skip_list = []

    # step 2: moving other files to 'NEW'
    for site in site_details['Site_Code']:
        site_dir = os.path.join(direc_archive, site)
        if not os.path.exists(site_dir):
            continue

        for file in su.get_files(site_dir):
            if file in skip_list:
                continue

            source = os.path.join(site_dir, file)
            su.copy_file(source, direc_new)
            
    logging.info(f'Done moving.\n')

# %% --------- Find and read data files in "NEW" file directory ---------
files = su.get_files(direc_new)

for XRF_file in files:
    file_path = os.path.join(direc_new, XRF_file)
    logging.info(f"Processing file: {XRF_file}")
    
    # Read the XRF data using pandas
    df_xrf = pd.read_excel(file_path)
    XRFtitle = df_xrf.columns.tolist()
    
    # Find rows where 'Pos' is not NaN
    rows = df_xrf['Pos'].notna()
    df_xrf = df_xrf.loc[rows,:]
    xrf_IDs = df_xrf.loc[:, 'Ident']  # filter ID
    site_ID = xrf_IDs.str.slice(0, 4).unique()

    if len(site_ID) > 1:
        logging.error('File contain data from multiple sites. Please separate!')

    site_ID = site_ID[0]

    # ------- Find index for metals in XRF file -----
    index_metals = []
    for col in xrf_cols:
        index = XRFtitle.index(col)
        index_metals.append(index)

    # Isolate XRF data, units of ug/cm2
    xrf_dataB = df_xrf.iloc[:, index_metals].astype(float).values  # Assuming the data can be converted to float

    # Multiply by filter surface area to get in units of ug, multiply by 1000 to get in final units of ng
    xrf_data = (xrf_dataB * filter_area) * 1000

    # Read the master file and write XRF data to it
    master_file = os.path.join(direc_master, f'{site_ID}_master.csv')
    master_data = su.read_master(master_file)
  
    # Writing new XRF data to Master_XRF
    master_idx = [i for i, master_id in enumerate(master_data['FilterID']) if master_id in xrf_IDs.values]
    for mid, metal in enumerate(metals):
        master_data.loc[master_idx, f'{metal}_XRF_ng'] = xrf_data[:,mid]

    # Write master file
    su.write_master(master_file, master_data)
    logging.info(f'Done writing to {site_ID} master file')
    
    # Move raw data to archive
    destination = f'{direc_archive}/{site_ID}/'
    su.force_movefile(file_path, destination)

    # adding dates to elog for future reference 
    # (scirpt F finds MDL using analysis dates)
    city = site_details.loc[site_details['Site_Code']==site_ID, 'City'].values[0]
    sheet_name = f'{site_ID}_{city}'
    xls = pd.ExcelFile(elogfname)
    if sheet_name in xls.sheet_names:
        df_elog = pd.read_excel(elogfname, sheet_name=sheet_name)
    else:
        df_elog = pd.DataFrame(columns=[
            'Cartridge ID', 'Filter ID', 'Analysis ID', 
            'Analysis Date', 'Data Uploaded', 'Comments'
        ])
        
    # Extract Cartridge ID from file_path
    cartridge_id = os.path.basename(file_path).split('.')[0]
        
    # find if filters are in e log and update accordingly
    existing_filters = set(df_elog['Filter ID'])
    new_rows = []

    for _, row in df_xrf.iterrows():
        ident = row['Ident']
        time_val = row['Time']
        if ident in existing_filters:
            df_elog.loc[df_elog['Filter ID'] == ident, 'Analysis Date'] = time_val
        else:
            new_rows.append({
                'Cartridge ID': cartridge_id,
                'Filter ID': ident,
                'Analysis ID': None,
                'Analysis Date': time_val,
                'Data Uploaded': None,
                'Comments': None
            })
            
    # Add new rows to df_elog if needed
    if new_rows:
        df_elog = pd.concat([df_elog, pd.DataFrame(new_rows)], ignore_index=True)

    # save the updated elog data
    with pd.ExcelWriter(elogfname, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
        df_elog.to_excel(writer, sheet_name=sheet_name, index=False)
    
    
logging.info(f"END of XRF data processing for {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

