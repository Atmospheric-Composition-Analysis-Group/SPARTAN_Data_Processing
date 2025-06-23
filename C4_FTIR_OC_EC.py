########################### CODE DESCRIPTION ################################
# PURPOSE: ------------------------------------------------------------
# Process FTIR OC/EC data and integrate into SPARTAN master files
#
# Input: -------------------------------------------------------------------
# Analysis_Data/FTIR/FTIR_data/spartan_FTIR_Batch_##.csv
# Analysis_Data/Master_files/$SiteCode_master.csv
#
# Output: -------------------------------------------------------------------
# Analysis_Data/Master_files/$SiteCode_master.csv
# 
# Suppliment: ---------------------------------------------------------------
# None
# 
# To re-process: -----------------------------------------------------------
# Set Redo_All_Archive to 1
#
# Created by: --------------------------------------------------------------
# Nidhi Anchan & Haihui Zhu, April 2025
# 

import os
import pandas as pd
import numpy as np
import shutil
from datetime import datetime
import logging
import spt_utils as su

# ============================
# Configuration
# ============================
debug_mode = 0
Redo_All_Archive = 0

# ============================
# Setup Directories (Original)
# ============================
direc = su.find_root_dir(debug_mode)
direc_master = os.path.join(direc, 'Analysis_Data', 'Master_files')
direc_FTIR = os.path.join(direc, 'Analysis_Data', 'FTIR', 'FTIR_data')
direc_archive = os.path.join(direc, 'Analysis_Data', 'Archived_Filter_Data', 'FTIR')
log_file_path = os.path.join(direc, 'Public_Data/Data_Processing_Records/FTIR',
                                datetime.now().strftime('%Y-%m-%d-%H%MSS') + '_FTIR_Record.log')

# ============================
# Setup Directories (TEST PATHS)
# ============================
# direc_FTIR = 'release_test/FTIR_data'
# direc_master = 'release_test/Master_files'
# direc_archive = 'release_test/FTIR_archive'
# log_file_path = os.path.join('release_test', datetime.now().strftime('%Y-%m-%d-%H%MSS') + '_FTIR_Record.log')

# Initialize logging
su.setup_logging(log_file_path)
logging.info("C4 started")


# ============================
# Optimized Functions
# ============================

def process_ftir_file(file_path):

    df = pd.read_csv(file_path)
    df = df[df['Parameter'].isin(['OC_ftir', 'EC_ftir', 'OM'])]

    if df.empty:
        logging.info(f"No OC or IC data found for {file_path}")
        return None

    logging.info(f'Processing {file_path}')  

    # Handle unusual MDLs (missing or Inf)
    df['MDL_ug_m3'] = df['MDL_ug_m3'].replace(np.inf, np.nan)

    # Create complete pivot table structure
    expected_params = ['OC_ftir', 'EC_ftir', 'OM']
    pivot_cols = pd.MultiIndex.from_product(
        [['MassLoading_ug', 'MDL_ug_m3'], expected_params],
        names=['Measurement', 'Parameter']
    )

    # Pivot with guaranteed column structure
    df_pivot = df.pivot_table(
        index=['Site', 'FilterId'],
        columns='Parameter',
        values=['MassLoading_ug', 'MDL_ug_m3'],
        aggfunc='first'
    ).reindex(columns=pivot_cols)  # Ensure all expected columns exist
    # Flatten column multi-index
    df_pivot.columns = [f"{meas}_{param}" for meas, param in df_pivot.columns]
    df = df_pivot.reset_index()
    # Fill MDL NaNs with column averages
    for param in expected_params:
        mdl_col = f'MDL_ug_m3_{param}'
        
        # Calculate mean ignoring NaNs
        col_mean = df[mdl_col].mean(skipna=True)
        
        # Fill NaNs if we have valid mean, else leave as NaN
        if not np.isnan(col_mean):
            df[mdl_col] = df[mdl_col].fillna(col_mean)
        else:
            logging.warning(f"No valid MDL values found for {param}, keeping NaNs")
    
    # lastly, set filter ID to standard format
    df['FilterId'] = su.filter_id_format(df['FilterId'])
    
    return df
    

def update_master_data(master_data, ftir_data):
    """Update master data with FTIR results"""
    updates = []
    for _, row in ftir_data.iterrows():
        mask = master_data['FilterID'] == row.FilterId
            
        # can go directly to master file for readability
        master_data.loc[mask, 'OC_FTIR_ug'] = row.MassLoading_ug_OC_ftir
        master_data.loc[mask, 'EC_FTIR_ug'] = row.MassLoading_ug_EC_ftir
        master_data.loc[mask, 'OM_FTIR_ug'] = row.MassLoading_ug_OM
        master_data.loc[mask, 'OC_FTIR_MDL'] = row.MDL_ug_m3_OC_ftir
        master_data.loc[mask, 'EC_FTIR_MDL'] = row.MDL_ug_m3_EC_ftir
        master_data.loc[mask, 'OM_FTIR_MDL'] = row.MDL_ug_m3_OM # currently don't have OM MDL
        
        # Check MDLs, add flag if data lower than MDL
        if row.MassLoading_ug_OC_ftir < row.MDL_ug_m3_OC_ftir:
            master_data.loc[mask, 'Flags'] = master_data.loc[mask, 'Flags'].apply(lambda x: su.add_flag(x, 'FTIR-OC<MDL'))
            logging.info(f"FTIR-OC for {row.FilterId} ({row.MassLoading_ug_OC_ftir}) is lower than MDL {row.MDL_ug_m3_OC_ftir}")
        if row.MassLoading_ug_EC_ftir < row.MDL_ug_m3_EC_ftir:
            master_data.loc[mask, 'Flags'] = master_data.loc[mask, 'Flags'].apply(lambda x: su.add_flag(x, 'FTIR-EC<MDL'))
            logging.info(f"FTIR-EC for {row.FilterId} ({row.MassLoading_ug_EC_ftir}) is lower than MDL {row.MDL_ug_m3_EC_ftir}")
        if row.MassLoading_ug_OM < row.MDL_ug_m3_OM: # currently don't have OM MDL
            master_data.loc[mask, 'Flags'] = master_data.loc[mask, 'Flags'].apply(lambda x: su.add_flag(x, 'FTIR-EC<MDL'))
            logging.info(f"FTIR-EC for {row.FilterId} ({row.MassLoading_ug_OM}) is lower than MDL {row.MDL_ug_m3_OM}")
            
    return master_data

# ============================ 
# Main Processing
# ============================

# Archive reprocessing (if enabled)
if Redo_All_Archive:
    try:
        with open(os.path.join(direc_archive, 'No_need_to_reprocess.txt'), 'r') as f:
            skip_files = set(line.strip() for line in f)
    except FileNotFoundError:
        skip_files = set()

    for file in su.get_files(direc_archive):
        if file not in skip_files:
            try:
                shutil.copy(os.path.join(direc_archive, file), direc_FTIR)
            except Exception as e:
                logging.warning(f"Archive copy failed: {file} - {str(e)}")

# Process FTIR files
for file in su.get_files(direc_FTIR):
    # not need to check if file starts with '.' or '~'
    # su.get_files handled this already
        
    file_path = os.path.join(direc_FTIR, file)
    ftir_data = process_ftir_file(file_path)
    
    if ftir_data is None:
        continue
        
    # Group by site for efficient processing
    for site, group in ftir_data.groupby('Site'):
        master_path = os.path.join(direc_master, f"{site}_master.csv")
        master_data = su.read_master(master_path) # su. read_master handles the situation when master_path does not exist
        updated_data = update_master_data(master_data, group)
        su.write_master(master_path, updated_data)
            
    # Archive processed file
    su.force_movefile(file_path, direc_archive)

logging.info(f"Processing complete.")
