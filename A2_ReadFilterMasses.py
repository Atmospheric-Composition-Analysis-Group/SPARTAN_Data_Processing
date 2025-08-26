########################### CODE DESCRIPTION ################################
# PURPOSE: ------------------------------------------------------------
# Read filter masses (output from A1) and add to the master data file
# Perform QC on filter masses.
#
# Input: -------------------------------------------------------------------
# Analysis_Data/Filter_Masses/Masses_by_site/$SiteCode_masses.xlsx
#
# Output: -------------------------------------------------------------------
# Analysis_Data/Master_files/$SiteCode_master.csv
# 
# Suppliment: ------------------------------------------------------
# Site_Sampling/Site_details.xlsx
# 
# To re-process: -----------------------------------------------------------
# Delete wrong weights in Master Files and rerun this script. 
# If there are values in the mass column in the master file, it won't be replaced
#
# Created by: --------------------------------------------------------------
# Crystal Weagle, 14 May 2020
# Revised by Haihui Zhu, April 2025
# 

import pandas as pd
import numpy as np
import os
from datetime import datetime
import logging
import spt_utils as su  # common functions are defined here

# ============================
# Initial Setup and Configuration
# ============================

# Define the debug mode for finding the root directory
debug_mode = 0  # 0 for standard mode, 1 for debug mode

# Find the root directory of the project (e.g., parent folder for data files)
direc = su.find_root_dir(debug_mode)

# Constant
neg_threshold = -0.3

# ============================
# Directory & File Paths
# ============================
# Define paths for input, output, and archived data files
direc_new = os.path.join(direc, 'Analysis_Data/Filter_Masses/Masses_by_site')
direc_master = os.path.join(direc, 'Analysis_Data/Master_files')
log_file_path = os.path.join(direc, 'Public_Data/Data_Processing_Records/MTL_masses',
                             datetime.now().strftime('%Y-%m-%d-%H%MSS') + '_MTLmasses_Record.txt')

# Read the site info table
SiteInfoFile = f'{direc}/Site_Sampling/Site_details.xlsx'
site_details = pd.read_excel(SiteInfoFile)

# Initialize logging to track script operations
su.setup_logging(log_file_path)
logging.info("A2 started")

# ============================
# Define Functions
# ============================
def add_entry(master_data, filter,  mass, mass_type, sampling_mode, barcode, cartid, lotid):
    if filter in master_data['FilterID'].values:
        mask = master_data['FilterID'] == filter
        master_data.loc[mask, 'mass_ug'] = mass
        master_data.loc[mask, 'Filter_Barcode'] = barcode
        master_data.loc[mask, 'CartridgeID'] = cartid
        master_data.loc[mask, 'LotID'] = lotid
        master_data.loc[mask, 'projectID'] = sampling_mode
        master_data.loc[mask, 'Mass_type'] = mass_type
        
    else:
        # Create template for new row with type-appropriate defaults
        new_row = {col:  np.nan for col in master_data.columns}
        
        master_data = pd.concat(
            [master_data, pd.DataFrame([new_row])],
            ignore_index=True
        )
     
        master_data.iloc[-1, master_data.columns.get_loc('FilterID')] = filter
        master_data.iloc[-1, master_data.columns.get_loc('mass_ug')] = mass
        master_data.iloc[-1, master_data.columns.get_loc('Filter_Barcode')] = barcode
        master_data.iloc[-1, master_data.columns.get_loc('CartridgeID')] = cartid
        master_data.iloc[-1, master_data.columns.get_loc('LotID')] = lotid
        master_data.iloc[-1, master_data.columns.get_loc('projectID')] = sampling_mode
        master_data.iloc[-1, master_data.columns.get_loc('Mass_type')] = mass_type
    
    return master_data
    

# ============================
# Start processing 
# ============================

# Find & Read Data
for idx, site in enumerate(site_details['Site_Code']):
    
    sampling_mode = site_details['Sampling_Mode'][idx][0]

    # Check if the master file for the site already exists
    master_file = f'{direc_master}/{site}_master.csv'
    master_data = su.read_master(master_file)
   
    # looping through sites and add new data to master files
    filename = f'{direc_new}/MTL_weighing_WashU/{site}_MTL_masses.csv'
    if os.path.exists(filename):  # file found
        logging.info(f'Reading MTL weighing file for {site}')

        mtl_data = pd.read_csv(filename, dtype=str)
        mtl_data['Net_Weight_ug'] = pd.to_numeric(mtl_data['Net_Weight_ug'], errors='coerce')
        mtl_data['AnalysisID_formated'] = su.filter_id_format(mtl_data['AnalysisID'])
        
        for filter in mtl_data['AnalysisID_formated']:
            filter_mask_master = master_data['FilterID'] == filter
            filter_mask_mtl = mtl_data['AnalysisID_formated'] == filter
            
            # if filter in master_data but mass = nan
            if master_data.loc[filter_mask_master, 'mass_ug'].isna().any():
                # two possibilities:
                # 1. filter is a lab blank
                # 2. filter mass has not been added 
               
                mass = mtl_data.loc[filter_mask_mtl, 'Net_Weight_ug'].values[0]

                # mass_type: keep existing label unless lab blank logic applies
                mass_type = master_data.loc[filter_mask_master, 'Mass_type'].values[0]
                 # Copy IDs from MTL files
                barcode = mtl_data.loc[filter_mask_mtl, 'Filter_Barcode'].values[0]
                cartid  = mtl_data.loc[filter_mask_mtl, 'CartridgeID'].values[0]
                lotid   = mtl_data.loc[filter_mask_mtl, 'LotID'].values[0]
                
                # for negative data below threshold, label mass as nan, mass_type as 5.
                if mass < neg_threshold and pd.isna(mass_type):
                    # for such filters, mass_type is either nan or 5. Example: lab blank
                    # if nan, label as 5 (invalid measurement)
                    mass_type = 5
                    logging.info(f'{filter} net weight too negative, label as invalid (mass_type = 5)') 
                
                master_data = add_entry(master_data, filter,  mass, mass_type, sampling_mode, \
                    barcode, cartid, lotid)
                
                
            # if the filer not in master; 
            if not filter_mask_master.any():
                
                mass = mtl_data.loc[filter_mask_mtl, 'Net_Weight_ug'].values[0]
                mass_type = np.nan # initialize as nan. will be populated in script B. 
                
                # for negative data below threshold, label mass as nan, mass_type as 5.
                if mass < neg_threshold:
                    mass_type = 5
                    logging.info(f'{filter} net weight too negative, set to nan') 
                       
                master_data = add_entry(master_data, filter,  mass, mass_type, sampling_mode, \
                    mtl_data.loc[filter_mask_mtl, 'Filter_Barcode'].values[0], mtl_data.loc[filter_mask_mtl, 'CartridgeID'].values[0], mtl_data.loc[filter_mask_mtl, 'LotID'].values[0])
        # Sort the master data by FilterID and reset index
        master_data = master_data.sort_values(by="FilterID").reset_index(drop=True)
        # Write file 
        su.write_master(master_file, master_data)
        logging.info(f'Done adding filter mass for {site}')
        
    else:
        logging.info(f'No MTL weighing file exists for {site}')
        continue  # no post-weights available from MTL or Dal weighing

logging.info(f'END of filter mass processing for {datetime.now()}')
