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
#For older format files that are missing MassLoading_ug but have Concentration_ug_m3 and Volume_m3, we can backfill MassLoading_ug using the formula:
def backfill_massloading_from_conc_volume(df_full, df, file_path):
    """
    Backfill missing MassLoading_ug in filtered FTIR rows using:
        MassLoading_ug = Concentration_ug_m3 * Volume_m3

    Also updates the same rows in df_full and saves the corrected full file.
    """

    required_cols = ['Volume_m3', 'Parameter', 'MassLoading_ug', 'Concentration_ug_m3']
    missing_cols = [col for col in required_cols if col not in df.columns]

    if missing_cols:
        msg = (
            f"{os.path.basename(file_path)}: Missing required column(s): {missing_cols}. "
            f"This looks like an old or incompatible FTIR file format. Stopping script."
        )
        logging.error(msg)
        raise SystemExit(msg)

    # convert to numeric safely
    df['MassLoading_ug'] = pd.to_numeric(df['MassLoading_ug'], errors='coerce')
    df['Concentration_ug_m3'] = pd.to_numeric(df['Concentration_ug_m3'], errors='coerce')
    df['Volume_m3'] = pd.to_numeric(df['Volume_m3'], errors='coerce')

    missing_before = df['MassLoading_ug'].isna()

    if missing_before.any():
        missing_params = df.loc[missing_before, 'Parameter'].value_counts().to_dict()
        logging.warning(
            f"{os.path.basename(file_path)}: MassLoading_ug is empty for {missing_before.sum()} row(s), "
            f"Attempting backfill using Concentration_ug_m3 * Volume_m3. "
            f"Parameter breakdown = {missing_params}"
        )
    else:
        logging.info(
            f"{os.path.basename(file_path)}: MassLoading_ug already present for all OC/EC/OM rows."
        )

    fill_mask = (
        df['MassLoading_ug'].isna() &
        df['Concentration_ug_m3'].notna() &
        df['Volume_m3'].notna()
    )

    if fill_mask.any():
        df.loc[fill_mask, 'MassLoading_ug'] = (
            df.loc[fill_mask, 'Concentration_ug_m3'] *
            df.loc[fill_mask, 'Volume_m3']
        )

        # update same rows in full dataframe
        df_full.loc[df.index[fill_mask], 'MassLoading_ug'] = df.loc[fill_mask, 'MassLoading_ug']

        filled_params = df.loc[fill_mask, 'Parameter'].value_counts().to_dict()
        logging.info(
            f"{os.path.basename(file_path)}: Backfilled MassLoading_ug for {fill_mask.sum()} row(s); "
            f"parameter breakdown = {filled_params}"
        )

        try:
            df_full.to_csv(file_path, index=False)
            logging.info(f"{os.path.basename(file_path)}: Saved updated file after filling MassLoading_ug")
        except Exception as e:
            msg = (
                f"{os.path.basename(file_path)}: Backfilled MassLoading_ug in memory, "
                f"but could not save updated file: {e}. Stopping script."
            )
            logging.error(msg)
            raise SystemExit(msg)

    remaining_mask = df['MassLoading_ug'].isna()
    if remaining_mask.any():
        remaining_params = df.loc[remaining_mask, 'Parameter'].value_counts().to_dict()
        logging.warning(
            f"{os.path.basename(file_path)}: {remaining_mask.sum()} row(s) still have missing MassLoading_ug; "
            f"parameter breakdown = {remaining_params}"
        )

    return df

def process_ftir_file(file_path):

    df = pd.read_csv(file_path)
    # Standardize column names of MDL_ug_m3 and mass loading
    df.columns = (
        df.columns
          .str.strip()
          .str.replace(r'^MDL$', 'MDL_ug_m3', regex=True)
          .str.replace(r'^MassLoading$', 'MassLoading_ug', regex=True)
          .str.replace(r'^Concentration$', 'Concentration_ug_m3', regex=True)
          .str.replace(r'^Volume$', 'Volume_m3', regex=True)
    )

    required_cols = ['Site', 'FilterId', 'Volume_m3', 'Parameter', 'MassLoading_ug', 'Concentration_ug_m3', 'MDL_ug_m3']
    missing_cols = [col for col in required_cols if col not in df.columns]

    if missing_cols:
        msg = (
            f"{os.path.basename(file_path)}: Missing required column(s): {missing_cols}. "
            f"This looks like an old or incompatible FTIR file format. Stopping script."
        )
        logging.error(msg)
        raise SystemExit(msg)

    # Keep full copy so helper can save corrected full file without dropping other rows
    df_full = df.copy()

    df = df[df['Parameter'].isin(['OC_ftir', 'EC_ftir', 'OM'])].copy()

    if df.empty:
        logging.info(f"No OC or EC data found for {file_path}")
        return None

    logging.info(f'Processing {file_path}')

    # Handle unusual MDLs (missing or Inf)
    df['MDL_ug_m3'] = df['MDL_ug_m3'].replace(np.inf, np.nan)

    # Use helper to backfill MassLoading_ug where possible
    df = backfill_massloading_from_conc_volume(df_full, df, file_path)

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
