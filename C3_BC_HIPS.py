#
# PURPOSE: ------------------------------------------------------------
# This script processes Black Carbon (BC) data measured using the
# HIPS method. It performs quality control, calculates BC mass, and writes the
# processed data into corresponding master files for each sampling site.
# 
# Input: -----------------------------------------------------------
# Master files
# /Analysis_Data/HIPS/Reports/*.xlsx
# 
# Output: ----------------------------------------------------------
# Master files
#
# Suppliment: ------------------------------------------------------
# Site_Sampling/Site_details.xlsx
# 
# To re-process a HIPS file: ------------------------------------------- 
# Turn Redo_All_Archieve to 1
# 
# Created by: --------------------------------------------------------------
# Nidhi Anchan, Haihui Zhu, and Crystal Weagle
# May-2025

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
Redo_All_Archive = 0

# ============================
# Directory & File Paths
# ============================
direc = su.find_root_dir(debug_mode)
direc_master = os.path.join(direc, 'Analysis_Data', 'Master_files')
direc_HIPS = os.path.join(direc, 'Analysis_Data', 'HIPS', 'Reports')
direc_archive = os.path.join(direc, 'Analysis_Data', 'Archived_Filter_Data', 'HIPS')
log_file_path = os.path.join(direc, 'Public_Data/Data_Processing_Records/HIPS',
                                datetime.now().strftime('%Y-%m-%d-%H%MSS') + '_HIPS_Record.txt')

# # TEST Dir
# direc_master = os.path.join('./release_test',  'Master_files')
# direc_HIPS = os.path.join('./release_test', 'HIPS_Reports')
# direc_archive = os.path.join('./release_test', 'HIPS_archive')
# log_file_path = os.path.join('./release_test', 
#                                 datetime.now().strftime('%Y-%m-%d-%H%MSS') + '_HIPS_Record.txt')


# Log file setup
su.setup_logging(log_file_path)  # Initialize logging
logging.info("C3 started")

# Read the site info table
SiteInfoFile = f'{direc}/Site_Sampling/Site_details.xlsx'
site_details = pd.read_excel(SiteInfoFile)
Site_codes = site_details.iloc[:, 0].astype(str).tolist()


# ============================
# Define Functions
# ============================
def process_hips_file(fpath, sigma):
    fname = os.path.basename(fpath)
    try:
        if fname.endswith('.csv'):
            rawdata = pd.read_csv(fpath)
            logging.info(f'{fpath} read ')
        else:
            xl = pd.ExcelFile(fpath)
            if 'HIPS Results' not in xl.sheet_names:
                logging.warning(f"Sheet 'HIPS Results' not found in {fname}")
                return None
            rawdata = xl.parse('HIPS Results')
            logging.info(f'{fpath} read ')
    except Exception as e:
        logging.warning(f"Skipping {fname}, read error: {e}")
        return None
    colnames = rawdata.columns.str.lower()
    
    try:
        # Step 1: Make column name search case-insensitive
        lower_colnames = rawdata.columns.str.lower()

       # Step 2: First look for columns named exactly 'id' (case-insensitive)
        id_cols = [col for col in rawdata.columns if col.lower().strip() == 'id']

        # If more than one column named 'id' is found, fall back to 'filter id' or 'filterid'
        if len(id_cols) == 1:
            filter_id_cols = id_cols
        else:
            filter_id_cols = [col for col in rawdata.columns if col.lower().strip() in ['filter id', 'filterid']]

        # Step 2: Look for 'filter id' or 'filterid' specifically
        possible_names = ['filter id', 'filterid','id'] # revise to first try 'id' and for more than 2 columns try 'filter id' and 'filterid'
        filter_id_cols = [col for col in rawdata.columns if col.lower() in possible_names]

        # Step 3: Fallback or raise an error if not found
        if len(filter_id_cols) == 0:
            raise ValueError("No 'filter ID' column found in the input file.")
        elif len(filter_id_cols) > 1:
            print(f"Multiple candidate columns for filter ID found: {filter_id_cols}. Choosing the first.")
            
        # Step 4: Extract the actual filter ID values
        filter_id = rawdata[filter_id_cols[0]].astype(str)

        filter_id = su.filter_id_format(filter_id)
        tau = rawdata.loc[:, colnames.str.contains('tau')].iloc[:, 0]
        dep_area = rawdata.loc[:, colnames.str.contains('depositarea|deposit area')].iloc[:, 0]
    except Exception as e:
        logging.warning(f"Skipping {fname}, missing critical columns: {e}")
        return None

    filter_type = rawdata.loc[:, colnames.str.contains('filter type')].iloc[:, 0] if any(colnames.str.contains('filter type')) else pd.Series([''] * len(filter_id))
    comments = rawdata.loc[:, colnames.str.contains('comment')].iloc[:, 0] if any(colnames.str.contains('comment')) else pd.Series([''] * len(filter_id))

    remove_idx = filter_id[filter_id.str.contains('LB|FB', case=False, na=False)].index.tolist()
    if not filter_type.empty:
        remove_idx += filter_type[filter_type.str.contains('LB|FB', case=False, na=False)].index.tolist()
    
    drop_idx = list(set(remove_idx))

    filter_id = filter_id.drop(index=drop_idx).reset_index(drop=True)
    tau = tau.drop(index=drop_idx).reset_index(drop=True)
    dep_area = dep_area.drop(index=drop_idx).reset_index(drop=True)
    comments = comments.drop(index=drop_idx).reset_index(drop=True)
    logging.info(f"{fname}: {len(filter_id)} valid rows after removing blanks")


    if filter_id.empty:
        logging.warning(f"Skipping {fname}, all rows dropped after blank filtering")
        return None

    formatted_id, site_code_list = [], []
    fid_list = filter_id.tolist()
    formatted_id = su.filter_id_format(pd.Series(fid_list)).tolist()
    site_code_list = [fid[:4] for fid in formatted_id]


    bcmass = dep_area * tau / sigma

    return pd.DataFrame({
        'Site': site_code_list,
        'FormattedID': formatted_id,
        'BCMass': bcmass,
        'Comment': comments
    })


def update_master_file(master_data, group):
    for _, row in group.iterrows():
        fid = str(row['FormattedID']).strip()
        bcm = row['BCMass']
        comment = row['Comment']

        # match_idx = master_data[master_data['master_ids'].str.strip() == fid].index
        match_idx = master_data[master_data['FilterID'].astype(str).str.strip() == str(fid).strip()].index

        if match_idx.empty:
            logging.warning(f"{fid} not found in master file")
            continue

        i = match_idx[0]

        # Update BC_HIPS_ug column
        current_val = master_data.at[i, 'BC_HIPS_ug']
        if pd.isna(current_val) or abs(current_val - bcm) < 0.3:
            master_data.at[i, 'BC_HIPS_ug'] = bcm
        else:
            logging.info(f"Overwriting BC mass for {fid} from {current_val} to {bcm}")
            master_data.at[i, 'BC_HIPS_ug'] = bcm

        # Update flags
        for code, keywords in {
            'VIH': ['weight', 'inhomogenous', 'concentrated'],
            'FS': ['stretch', 'wrinkle'],
            'PH': ['pin', 'hole'],
            'VCT': ['contaminat', 'streak', 'particle', 'flecks'],
            'SGP': ['screen pattern', 'grid'],
            'MTL': ['metallic', 'reflect'],
            'SUS': ['suspect', 'compromise', 'questionable']
        }.items():
            if isinstance(comment, str) and any(kw in comment.lower() for kw in keywords):
                master_data.at[i, 'Flags'] = su.add_flag(master_data.at[i, 'Flags'], code)

    return master_data


# ============================
# Start processing 
# ============================
# Constants
SIGMA = 0.1 # absorption coefficient. We are recently changing it from 0.06 to 0.1. 

# Archive reprocessing (if enabled)
if Redo_All_Archive:
    logging.info(f'Copying archived HIPS report for reprocessing')
    skipfile_path = os.path.join(direc_archive, 'No_need_to_reprocess.txt')
    if os.path.exists(skipfile_path):
        with open(skipfile_path, 'r') as f:
            skip_list = [line.strip() for line in f if line.strip()]
    else:
        skip_list = []

    for site in Site_codes:
        site_dir = os.path.join(direc_archive, site)
        if not os.path.exists(site_dir):
            continue

        for file in su.get_files(site_dir):
            if file in skip_list:
                logging.info(f"Skipping file in skiplist: {file}")
                continue

            source = os.path.join(site_dir, file)
            # destination = os.path.join(direc_HIPS, file)
            try:
                su.copy_file(source, direc_HIPS)
            except Exception as e:
                logging.warning(f"Failed to copy {file}: {e}")

    # Clear HIPS values in master files
    logging.info(f'Resetting HIPS BC data in master files...')
    for site in Site_codes:
        master_path = os.path.join(direc_master, f"{site}_master.csv")
        master_data = su.read_master(master_path, create_new=True)
        if master_data is None:
            raise RuntimeError(f"Failed to load master data from {master_path}")
        master_data['BC_HIPS_ug'] = np.nan

        su.write_master(master_path, master_data)
    logging.info(f'Done.')

                
# Loop through each HIPS file in Reports
# Get all HIPS report files
all_files = su.get_files(direc_HIPS)

for fname in all_files:
    fpath = os.path.join(direc_HIPS, fname)
    result = process_hips_file(fpath, SIGMA)

    if result is None or result.empty:
        logging.info(f"{fname}: No valid data found.")
        continue

    # Drop invalid or blank 'Site' values
    result = result.dropna(subset=['Site'])
    result = result[result['Site'].astype(str).str.strip() != '']

    # Collect unique sites in this report
    site_codes_in_report = (
        result['Site'].astype(str).str.strip().dropna().unique().tolist()
    )
    logging.info(f"Detected sites in {fname}: {site_codes_in_report}")

    # Update master files for each site
    for site, group in result.groupby('Site'):
        site_code = site.strip()
        if site_code == '' or site_code.lower().startswith('nan'):
            logging.warning(f"Skipping invalid site name in {fname}.")
            continue

        master_path = os.path.join(direc_master, f"{site_code}_master.csv")
        master_data = su.read_master(master_path)

        if master_data is None or master_data.empty:
            logging.warning(f"{master_path} is empty. Skipping.")
            continue

        updated = update_master_file(master_data, group)
        su.write_master(master_path, updated)

        # Extract unique site codes from the current result DataFrame
        sites_in_file = result['Site'].dropna().astype(str).str.strip().unique().tolist()

        # Archive the file to site-specific folders (only for relevant sites)
        all_success = True
        archive_dest = os.path.join(direc_archive, site_code)
        os.makedirs(archive_dest, exist_ok=True)
        
        if su.copy_file(fpath, archive_dest):
            logging.info(f"Successfully archived {fpath} to {archive_dest}")
        else:
            logging.warning(f"Failed to archive {fpath} to {archive_dest}")
            all_success = False  # Track failure

    # Remove file only if all copies were successful
    if all_success and os.path.exists(fpath):
        try:
            os.remove(fpath)
            logging.info(f"Archived and removed {fpath}")
        except Exception as e:
            logging.error(f"Error removing {fpath}: {e}")
    elif not all_success:
        logging.warning(f"{fpath} was not removed due to archiving failures.")



                

logging.info("C3 HIPS BC processing complete.")


                

