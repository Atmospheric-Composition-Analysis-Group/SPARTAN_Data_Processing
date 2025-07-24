# 
# PURPOSE: ------------------------------------------------------------
# Read filter masses from MTL weighing system output ($SiteCode_CartID_pre/post_weight.xlsx)
# and write to site-specific file ($SiteCode_MTL_masses.csv). 
#
# FILE REQUIREMENTS: ------------------------------------------------------
# Filter mass data from MTL should be exported as a Test 
# Report ("pre-weights" or "pre- & post-weights") for each individual cart-
# ridge and the comma-delimited text file should be saved as a '.xlsx' file. 
# The full filterID (i.e.'CAHA-0105-3') should be used in the MTL system.
#
# Input: -----------------------------------------------------------
# Analysis_Data/Filter_Masses/NEW_MTL/$SiteCode_CartID_{pre/post}_weight.xlsx
# 
# Output: ----------------------------------------------------------
# Analysis_Data/Filter_Masses/Masses_by_site/$SiteCode_masses.xlsx
#
# Suppliment: ------------------------------------------------------
# Acceptance_Testing_Summary_WashU.xlsx
# 
# To re-process a raw MTL file: ------------------------------------------- 
# Simply move the raw file to Filter_Masses/NEW_MTL, and rerun.
# Warnings about replacing existing weights will be given, unless
# replace_warning is set to 0. 
# Replace_warning is designed to prevent the same (and wrong) 
# filter ID in multiple raw MTL files.
# 
# Created by: --------------------------------------------------------------
# Haihui Zhu, Crystal Weagle, Emme Le Roy, and Chris Oxford 
# Dec-24-2024

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
replace_warning = 1  # if you expect to have masses being replaced and don't want to see the warnings, turn this off (= 0). 

# ============================
# Directory & File Paths
# ============================
direc_new = os.path.join(direc, 'Analysis_Data', 'Filter_Masses', 'NEW_MTL')
direc_output = os.path.join(direc, 'Analysis_Data', 'Filter_Masses', 'Masses_by_site', 'MTL_weighing_WashU') 
direc_archive = os.path.join(direc, 'Analysis_Data', 'Archived_Filter_Data', 'MTL_masses')
log_file_path = os.path.join(direc, 'Public_Data/Data_Processing_Records/MTL_masses',
                                datetime.now().strftime('%Y-%m-%d-%H%MSS') + '_MTLmasses_Record.txt')

# Log file setup
su.setup_logging(log_file_path)  # Initialize logging
logging.info("A1 started")

# Read Acceptance Testing Summary File (to fill in CartridgeIDs and LotIDs)
Sum_lots = pd.read_excel(os.path.join(direc, 'Analysis_Data/Acceptance_Testing/Acceptance_Testing_Summary_WashU.xlsx'),header=None)

# ============================
# Define Functions
# ============================
def find_lot_cart(filter_IDs, Sum_lots):

    site_lotIDs = []
    site_cartridgeIDs = []
    
    for fid in filter_IDs:
        locations = find_string_in_df(Sum_lots, fid)
        if len(locations) == 1:
            site_lotIDs.append(Sum_lots.iloc[locations[0][0], locations[0][1] - 2])
            site_cartridgeIDs.append(Sum_lots.iloc[locations[0][0], locations[0][1] - 1])
        elif len(locations) == 0:
            logging.error(f'{fid} not found in the acceptence sheet!')
        elif len(locations) >1:
            logging.error(f'More than one {fid} found in the acceptence sheet!')
            
    return site_lotIDs, site_cartridgeIDs

def check_and_fill(existing_table, var, idx_massfile, new_value, filterID, replace_warning):
    old_value = existing_table.at[idx_massfile, var]

    # If new_value is numeric (int or float)
    if isinstance(new_value, (int, float)):
        if not np.isnan(new_value) and old_value != new_value:
            existing_table.at[idx_massfile, var] = new_value
            if replace_warning == 1:
                logging.info(f'{filterID} {var} old value {old_value} replaced with {new_value}')
    
    else:  # Assume it's a string or object
        # Don't update if old value is already informative (i.e., not 'U' or NaN)
        if pd.isna(old_value) or old_value in ['U', '1']:
            existing_table.at[idx_massfile, var] = new_value
            if replace_warning == 1:
                logging.info(f'{filterID} {var} old value {old_value} replaced with {new_value}')
    
    return existing_table

# Function to find the cell containing a specific string
def find_string_in_df(df, search_string):
    # Apply a mask where each cell is checked for the string
    mask = df.applymap(lambda x: search_string in str(x))
    # Using np.where to get the indices of True values in the mask
    row_indices, col_indices = np.where(mask)
    # Create a list of tuples (row, column) for each found location
    locations = [(row, col) for row, col in zip(row_indices, col_indices)]
    
    return locations


# ============================
# Start processing 
# ============================
# Read all MTL filter masses test reports in the 'MTL_NEW' directory
files = su.get_files(direc_new)

for data_file in files:
    
    filename = os.path.join(direc_new, data_file)
    logging.info(f'Reading {data_file}')
    
    all_data = pd.read_excel(filename, header=None) # read all info in the sheet
    
    # read weighting type (pre or post)
    locations = find_string_in_df(all_data, 'Pretest Date')
    pre_date = all_data.iloc[locations[0][0],locations[0][1]+1]
    locations = find_string_in_df(all_data, 'Posttest Date')
    post_date = all_data.iloc[locations[0][0],locations[0][1]+1]
    if pd.isna(post_date):  #post date is NaN
        postweights_avail = 0  # will be used later
        post_date = 'Unknown'
    else:  # There are post-weights
        postweights_avail = 1
    
    # find site code 
    locations = find_string_in_df(all_data, 'Test Code')
    site = all_data.iloc[locations[0][0],locations[0][1]+1][0:4]

    # Read filter mass
    locations = find_string_in_df(all_data, 'Position')
    filter_data_raw = pd.read_excel(filename,skiprows=locations[0][0] )
    # Filter rows where 'Position' starts with 'Phase #'
    filter_data_raw = filter_data_raw[filter_data_raw['Position'].fillna('').str.startswith('Phase #')]
    filter_data_raw = filter_data_raw.iloc[:8, :]
    filter_data_raw['Pretest Standard Deviation'] = pd.to_numeric(filter_data_raw['Pretest Standard Deviation'], errors='coerce').round(2)
    filter_data_raw['Posttest Standard Deviation'] = pd.to_numeric(filter_data_raw['Posttest Standard Deviation'], errors='coerce').round(2)
    raw_titles = filter_data_raw.columns.tolist()

    # check if standard deviation too high:
    if (filter_data_raw['Pretest Standard Deviation'] > 2.5).sum() > 0:
        logging.error('Standard deviation of preweight mass too high! Measure the mass again.')
    elif (filter_data_raw['Posttest Standard Deviation'] > 2.5).sum() > 0:
        logging.error('Standard deviation of postweight mass too high! Measure the mass again.')

    # Filter IDs and Analysis IDs 
    filter_IDs = filter_data_raw['Filter ID'].values
    analysis_IDs = [ids[:9] for ids in filter_IDs]

    # Find barcodes
    filter_barcodes = filter_data_raw['MTL Filter ID'].values  

    # Find lot IDs and Cart IDs
    site_lotIDs, site_cartridgeIDs = find_lot_cart(filter_IDs, Sum_lots)

    # Get mass information for this site from MTL report file
    preweigh_date = [pre_date] * 8
    preweight_avgs = filter_data_raw['Corrected Pretest'].values  # units = mg
    preweight_stds = filter_data_raw['Pretest Standard Deviation'].values  # units = ug
    preweight_reps = filter_data_raw[' Pretest Repetitions'].values  # unitless # Don't delete the space before the string. This is how the excel sheets are formated. 

    if postweights_avail == 1:
        postweight_date = [post_date] * 8
        postweight_avgs = filter_data_raw['Corrected Posttest'].values  # units = mg
        postweight_stds = filter_data_raw['Posttest Standard Deviation'].values  # unit = ug
        postweight_reps = filter_data_raw[' Posttest Repetitions'].values  # unitless
        netweight = np.round((np.array(postweight_avgs, dtype=np.float64) - 
                              np.array(preweight_avgs, dtype=np.float64)) * 1000, 1)  # multiply by 1000 to get in units of ug
    elif postweights_avail == 0:
        postweight_date = [post_date] * 8
        postweight_avgs = np.full(len(preweight_avgs), np.nan)
        postweight_stds = np.full(len(preweight_avgs), np.nan)
        postweight_reps = np.full(len(preweight_avgs), np.nan)
        netweight = np.full(len(preweight_avgs), np.nan)

    # combine all variable to a all_data in the same format as the mass file:
    newtable = pd.DataFrame({
        'FilterID': filter_IDs,
        'AnalysisID': analysis_IDs,
        'Filter_Barcode': filter_barcodes,
        'CartridgeID': site_cartridgeIDs,
        'LotID': site_lotIDs,
        'Preweigh_Date': preweigh_date,
        'Preweight_avg_mg': preweight_avgs,
        'Preweight_std_ug': preweight_stds,
        'Preweight_NumReps': preweight_reps,
        'Postweigh_Date': postweight_date,
        'Postweight_avg_mg': postweight_avgs,
        'Postweight_std_ug': postweight_stds,
        'Postweight_NumReps': postweight_reps,
        'Net_Weight_ug': netweight
    })
    
    # ---- Reading MTL MASS SHEET  ----
    mass_file = os.path.join(direc_output, f"{site}_MTL_masses.csv")
    
    if os.path.exists(mass_file): 
        
        existing_table = pd.read_csv(mass_file)
            
        # find filter IDs NOT in MTL mass file
        newIDs_idx = np.where(~np.isin(filter_IDs, existing_table['FilterID']))[0]  # len should be either 8 or 0
        
        if len(newIDs_idx) == 0: # every filter is in mtl mass sheet
            
            # udpate number 
            for ii in range(len(filter_IDs)):
                # find index of filters in existing_table
                idx_massfile = existing_table[existing_table['FilterID'] == filter_IDs[ii]].index[0]

                existing_table = check_and_fill(existing_table, 'AnalysisID', idx_massfile, analysis_IDs[ii], filter_IDs[ii], replace_warning)
                existing_table = check_and_fill(existing_table, 'Filter_Barcode', idx_massfile, filter_barcodes[ii], filter_IDs[ii], replace_warning)
                existing_table = check_and_fill(existing_table, 'CartridgeID', idx_massfile, site_cartridgeIDs[ii], filter_IDs[ii], replace_warning)
                existing_table = check_and_fill(existing_table, 'LotID', idx_massfile, site_lotIDs[ii], filter_IDs[ii], replace_warning)

                existing_table = check_and_fill(existing_table, 'Preweigh_Date', idx_massfile, preweigh_date[ii], filter_IDs[ii], replace_warning)
                existing_table = check_and_fill(existing_table, 'Preweight_avg_mg', idx_massfile, preweight_avgs[ii], filter_IDs[ii], replace_warning)
                existing_table = check_and_fill(existing_table, 'Preweight_std_ug', idx_massfile, preweight_stds[ii], filter_IDs[ii], replace_warning)
                existing_table = check_and_fill(existing_table, 'Preweight_NumReps', idx_massfile, preweight_reps[ii], filter_IDs[ii], replace_warning)
                
                if postweights_avail == 1:
                    existing_table = check_and_fill(existing_table, 'Postweigh_Date', idx_massfile, postweight_date[ii], filter_IDs[ii], replace_warning)
                    existing_table = check_and_fill(existing_table, 'Postweight_avg_mg', idx_massfile, postweight_avgs[ii], filter_IDs[ii], replace_warning)
                    existing_table = check_and_fill(existing_table, 'Postweight_std_ug', idx_massfile, postweight_stds[ii], filter_IDs[ii], replace_warning)
                    existing_table = check_and_fill(existing_table, 'Postweight_NumReps', idx_massfile, postweight_reps[ii], filter_IDs[ii], replace_warning)

                    existing_table = check_and_fill(existing_table, 'Net_Weight_ug', idx_massfile, netweight[ii], filter_IDs[ii], replace_warning)
        
        elif len(newIDs_idx) == 8:
            existing_table = pd.concat([existing_table, newtable], ignore_index=True)
            
        else:
            logging.info(f'Warning: {site_cartridgeIDs[0]} has less than 8 filters.\n')
            
    else:
        existing_table = newtable
        logging.info(f'No MTL mass file exists for {site}. Creating new MTL mass file.\n')
        
    # ---- Writing to MTL MASS SHEET  ----
    existing_table['Postweight_NumReps'] = existing_table['Postweight_NumReps'].fillna(0).astype(int)
    existing_table['Preweight_NumReps'] = existing_table['Preweight_NumReps'].fillna(0).astype(int)
    existing_table.to_csv(mass_file, index=False)
    logging.info(f'Finished writing to {site} MTL mass file \n')


    # Copy file to archived MTL masses folder for future reference
    file_destination = os.path.join(direc_archive, site)
    status = su.copy_file(filename, file_destination)
    if status is True:
        os.remove(filename)
        logging.info('File has been moved from NEW_MTL folder to Archives')
    else:
        logging.info('Failed to copy to Archives. Cannot delete file from NEW_MTL data folder. ')
