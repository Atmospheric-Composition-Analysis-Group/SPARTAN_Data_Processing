#PURPOSE: ------------------------------------------------------------
# Read and perform QA/QC on BC SSR
# 
# Input: -------------------------------------------------------------------
# Analysis_Data/Black_Carbon/SSR_by_site/[Site_codes]_SSR.xlsx
#
# Output: -------------------------------------------------------------------
# Analysis_Data/Master_files/$SiteCode_master.csv
# 
# Suppliment: ------------------------------------------------------
# Site_Sampling/Site_details.xlsx
#
# To re-process: -----------------------------------------------------------
# Set Clean_Up_Master to 1 and run this script. 
#
# Created by: --------------------------------------------------------------
# Crystal Weagle, 21 November 2018
# Revised by Nidhi Anchan and Haihui Zhu, May, 2025

import pandas as pd  
import os
import numpy as np 
from datetime import datetime
import logging
import spt_utils as su



# ============================
# Initial Setup and Configuration
# ============================
debug_mode = 0
direc = su.find_root_dir(debug_mode)

Clean_Up_Master = 1

# ============================
# Directory & File Paths
# ============================
direc_master = os.path.join(direc, 'Analysis_Data/Master_files')
direc_SSR = os.path.join(direc, 'Analysis_Data/Black_Carbon/SSR_by_site')
log_file_path = os.path.join(direc, 'Public_Data/Data_Processing_Records/BC_SSR',
                             datetime.now().strftime('%Y-%m-%d-%H%MSS') + '_BC_SSR_Record.txt')

su.setup_logging(log_file_path)  # Initialize logging
logging.info("C2 started")


# ============================
# Define Functions & Constants
# ============================
# Constants
SSR_threshold = [20, 90]
SSR_std = 1.5
SA = 3.142
ABS_coeff = 0.1 # Updated from 0.06 to 0.1
R_100 = 100
R_avg_stretch = 85
q = 1 / 1.5

# Site details
SiteInfoFile = f'{direc}/Site_Sampling/Site_details.xlsx'
site_details = pd.read_excel(SiteInfoFile)
site_codes = site_details.iloc[:, 0].values
site_cities = site_details.iloc[:, 2].values

# ============================
# Processing Loop
# ============================
for idx, site_code in enumerate(site_codes):

    # SSR file path
    ssr_file = os.path.join(direc_SSR, f"{site_code}_SSR.xlsx")
    if not os.path.exists(ssr_file):
        logging.info(f"[SKIP] No SSR file for {site_code}")
        continue
    else:
        logging.info(f"Processing site: {site_code} ({site_cities[idx]})")

    # Master file
    master_file = os.path.join(direc_master, f"{site_code}_master.csv")
    df_master = su.read_master(master_file)

    # Check for possible Filter ID column names
    if 'FilterID' not in df_master.columns:
        logging.warning(f"No recognized Filter ID column in {master_file}. Skipping.")
        continue

    if df_master is None:
        continue
    
    if Clean_Up_Master == 1:
        # Prepare for reprocessing
        df_master['BC_SSR_ug'] = np.nan
        
    idx_reprocess = df_master['BC_SSR_ug'].isna()

    if not idx_reprocess.any():
        logging.info(f"[SKIP] No SSR BC needed need to be added for {site_code}")
        continue

    # Read SSR data
    ssr_df = pd.read_excel(ssr_file, header=7)  #because in all excel sheet the Sample ID is in row number 8 of excel

    # Format filter IDs
    ssr_df['Formatted_ID'] = su.filter_id_format(ssr_df['Sample ID#'])

    # Match IDs
    matched = df_master[df_master['FilterID'].isin(ssr_df['Formatted_ID']) & idx_reprocess]

    R_avg = ssr_df['Mean Reading'].values
    R_std = ssr_df['Std. Dev.'].values
    analysis_ids = ssr_df['Formatted_ID'].values

    BC_vals = []
    R_toolow_count = 0

    for i, fid in enumerate(analysis_ids):
        if fid not in df_master['FilterID'].values:
            # This filter is not in the master file, skip it
            continue

        row = df_master[df_master['FilterID'] == fid].iloc[0]
        LotID = row['LotID']
        mass_type = row['Mass_type']
        R = R_avg[i]
        Rdev = R_std[i]

        if pd.isna(R):
            BC_vals.append(np.nan)
            continue

        # Apply blank correction
        if pd.notna(LotID) and LotID != "Unknown":
            # Use stretched teflon correction
            BC = (q * SA / ABS_coeff) * np.log(R_avg_stretch / R)
        else:
            # Standard mesh teflon
            BC = (q * SA / ABS_coeff) * np.log(R_100 / R)

        # QA/QC filtering
        if Rdev > SSR_std:
            BC = -899
        elif R > SSR_threshold[1]:
            BC = 0
        elif R < SSR_threshold[0]:
            R_toolow_count += 1

        if pd.isna(mass_type) or mass_type == 3 or mass_type == 4:
            BC = -899

        # assign the calculated BC value to master data
        df_master.loc[df_master['FilterID'] == fid, 'BC_SSR_ug'] = BC


    # write to master file
    su.write_master(master_file, df_master)


logging.info(f"SSR BC processing completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
