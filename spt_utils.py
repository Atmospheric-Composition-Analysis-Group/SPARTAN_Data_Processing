# supporting functions for SPARTAN scripts
# version 1.0.0
# created by Haihui Zhu,  Nidhi Anchan
# April 2025

import pandas as pd
import numpy as np
import os
from datetime import datetime, timedelta
import logging
import shutil
import re

def setup_logging(log_file_path):
    # setLevel(logging.DEBUG): Setting the level to DEBUG includes all messages 
    # at all severity levels (DEBUG, INFO, WARNING, ERROR, CRITICAL). 
    
    # Configure logging
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        handlers=[
                            logging.FileHandler(log_file_path),  # File handler
                            logging.StreamHandler()              # Console handler
                        ])

def find_root_dir(debug_mode):
    currentdir = os.getcwd()
     # For mac OS
    if 'Volumes' in currentdir and debug_mode == 0:
        direc='/Volumes/rvmartin/Active/SPARTAN-shared/'
        
    elif  'Volumes' in currentdir and debug_mode == 1:
        direc='/Volumes/rvmartin/Active/SPARTAN-shared/Public_Data/Scripts_dev/' 
    
    elif 'storage1/fs1' in currentdir and debug_mode == 0: 
        direc='/storage1/fs1/rvmartin/Active/SPARTAN-shared' 

    elif 'storage1/fs1' in currentdir and debug_mode == 1:
        direc='/storage1/fs1/rvmartin/Active/SPARTAN-shared/Public_Data/Scripts_dev/'   
    # for Windows OS
    elif '\\storage1.ris.wustl.edu'  in currentdir:
        direc= '\\\\storage1.ris.wustl.edu\\rvmartin\\Active\\SPARTAN-shared\\'
    elif 'Z:' in currentdir:
        direc= '\\\\storage1.ris.wustl.edu\\rvmartin\\Active\\SPARTAN-shared\\'
    # Generic Compute1 symlink (e.g., /spartan_shared)
    elif os.path.exists("/spartan_shared"):
        direc='/spartan_shared'
    else:
        raise ValueError('Cannot identify the root directory') 
        # direc='/Volumes/rvmartin/Active/SPARTAN-shared/'# test only
    return direc
  
def get_files(directory):
    """Get list of non-hidden .csv/.xlsx files after renaming .xls files to .xlsx."""
    # Step 1: Get all non-hidden .xls, .xlsx, and .csv files
    filelist = [f for f in os.listdir(directory) 
                if (f.endswith('.xls') or f.endswith('.xlsx') or f.endswith('.csv')) 
                and not f.startswith('.') and not f.startswith('~')]
    
    # # Step 2: Rename .xls files to .xlsx [script can't open xls]
    # for filename in filelist:
    #     if filename.endswith('.xls'):
    #         old_path = os.path.join(directory, filename)
    #         base_name = os.path.splitext(filename)[0]  # Remove extension
    #         new_filename = base_name + '.xlsx'
    #         new_path = os.path.join(directory, new_filename)
    #         os.rename(old_path, new_path)
    
    # # Step 3: Get updated list of non-hidden .csv and .xlsx files
    # filelist = [f for f in os.listdir(directory) 
    #                 if (f.endswith('.csv') or f.endswith('.xlsx')) 
    #                 and not f.startswith('.') and not f.startswith('~')]
    return filelist

def copy_file(filename, destination, screen_print=False):
    # create destination folder if not exist
    if not os.path.exists(destination):
        os.makedirs(destination)
    
    # copy specified filename to specified destination
    try:
        shutil.copy(filename, destination)  #copy2 tries to copy metadata as well.
        status = True
        if screen_print is True: # for debugging
            logging.info(f'{filename} copied to {destination}')
    
    except Exception as e:
        status = False
        logging.info(f'failed to copy {filename} to {destination} (spt_utils) - {str(e)}')
    
    return status


def filter_id_format(filter_ids):
    # Convert SITE_XXX to SITE_0XXX
    new_filter_ids = filter_ids.iloc[:].astype(str).apply(
        lambda x: x[:5] + '0' + x[5:8] if len(x) < 9 else x
    )
     
    # 'AEAZ-0010-1' to 'AEAZ-0010'
    new_filter_ids = new_filter_ids.iloc[:].astype(str).apply(
        lambda x: x[:9] if len(x) == 11 else x
    )
    
    # '10001-AEAZ-1T' or '10001-AEAZ-1N' (in {site}_dates_flows.xlsx)
    new_filter_ids = new_filter_ids.iloc[:].astype(str).apply(
        lambda x: x[6:10] + '-0' + x[2:5] if len(x) == 13 else x
    )
    
    
    return new_filter_ids

def add_flag(orig_flags, new_flag):
    
    # Convert NaN (float) to empty string
    if isinstance(orig_flags, float) and np.isnan(orig_flags):
        orig_flags = ''

    # remove leading and trailing blank
    orig_flags = orig_flags.strip()

    # ----- remove obsolete flags ----------------
    # Sx-{ion}: IC, correlation not exist. No longer labelled
    orig_flags = re.sub(r'\bSx-\w+\b\s*;?\s*', '', orig_flags).strip('; ')
    # FSx-{ion}: Not sure what it means
    orig_flags = re.sub(r'\bFSx-\w+\b\s*;?\s*', '', orig_flags).strip('; ')
    
    
    # ----- preventing repeating flags ----------------
    if new_flag:
        while new_flag in orig_flags:
            orig_flags = orig_flags.replace(new_flag, '')  # prevent repeated part

    
    # # ------ a special check for 'qc2-FSx-Li' [obsolete] -------------------------------
    # # when there is 'qc2-x-Li', it means the previous section erroneously 
    # # removed 'FS' in 'qc2-FSx-Li'.
    # if 'qc2-x-Li' in orig_flags:
    #     orig_flags = orig_flags.replace('qc2-x-Li', 'qc2-FSx-Li')

    # ------ adding new flag -------------------------------------------------
    if not new_flag:  # just in case new flag is empty
        flag_out = orig_flags
    else:
        flag_out = f"{orig_flags}; {new_flag}"

    # ------ Formatting flag -------------------------------------------------
    while "; ;" in flag_out:
        flag_out = flag_out.replace("; ;", ";")  # remove extra ';' in the middle of a flag
    while "  " in flag_out:
        flag_out = flag_out.replace("  ", " ")  # remove extra blank in the middle of a flag

    if flag_out:
        while flag_out and (flag_out[0] == ';' or flag_out[0] == ' '):  # prevent starting with ';' or a blank
            flag_out = flag_out[1:]

    # Remove the leading and trailing whitespace
    flag_out = flag_out.strip()

    return flag_out

def force_movefile(filename, file_destination):
    
    # Convert to absolute paths to handle relative paths and symlinks correctly
    filename = os.path.abspath(filename)
    file_destination = os.path.abspath(file_destination)
    
    dest_path = os.path.join(file_destination, os.path.basename(filename))
    
    # Attempt to move the file
    try:
        shutil.move(filename, dest_path)
         
    except Exception as e:
        error_msg = str(e)
        print(f'Moving not successful: {filename} \n msg:{error_msg} \n', end='')


# master file related functions:
def master_heading():
    
    dtype_mapping = {
        # First 7 columns
        'FilterID': 'object', 'Filter_Barcode': 'object', 'CartridgeID': 'object', 'LotID': 'object', 'projectID': 'object', 
        'hours_sampled':'float64', 'Mass_type':'Int64',
        
        # Date/time columns (8)
        'start_year':'Int64', 'start_month':'Int64', 'start_day':'Int64', 'start_hour':'float64',
        'stop_year':'Int64', 'stop_month':'Int64', 'stop_day':'Int64', 'stop_hour':'float64',
        
        # Mass/volume (2)
        'mass_ug':'float64', 'Volume_m3':'float64',
        
        # IC columns (13)
        'IC_F_ug_T':'float64', 'IC_Cl_ug_T':'float64', 'IC_NO2_ug_T':'float64', 'IC_Br_ug_T':'float64', 'IC_NO3_ug_T':'float64', 'IC_PO4_ug_T':'float64', 'IC_SO4_ug_T':'float64',
        'IC_Li_ug_T':'float64', 'IC_Na_ug_T':'float64', 'IC_NH4_ug_T':'float64', 'IC_K_ug_T':'float64', 'IC_Mg_ug_T':'float64', 'IC_Ca_ug_T':'float64',
        
        # ICP columns (21)
        'Li_ICP_ng':'float64', 'Mg_ICP_ng':'float64', 'Al_ICP_ng':'float64', 'P_ICP_ng':'float64', 'Ti_ICP_ng':'float64', 'V_ICP_ng':'float64', 'Cr_ICP_ng':'float64', 'Mn_ICP_ng':'float64',
        'Fe_ICP_ng':'float64', 'Co_ICP_ng':'float64', 'Ni_ICP_ng':'float64', 'Cu_ICP_ng':'float64', 'Zn_ICP_ng':'float64', 'As_ICP_ng':'float64', 'Se_ICP_ng':'float64', 'Ag_ICP_ng':'float64',
        'Cd_ICP_ng':'float64', 'Sb_ICP_ng':'float64', 'Ba_ICP_ng':'float64', 'Ce_ICP_ng':'float64', 'Pb_ICP_ng':'float64',
        
        # XRF columns (26)
        'Na_XRF_ng':'float64', 'Al_XRF_ng':'float64', 'Si_XRF_ng':'float64', 'S_XRF_ng':'float64', 'Cl_XRF_ng':'float64', 'K_XRF_ng':'float64', 'Ca_XRF_ng':'float64', 'Ti_XRF_ng':'float64',
        'V_XRF_ng':'float64', 'Fe_XRF_ng':'float64', 'Zn_XRF_ng':'float64', 'Ce_XRF_ng':'float64', 'Pb_XRF_ng':'float64', 'As_XRF_ng':'float64', 'Co_XRF_ng':'float64', 'Cr_XRF_ng':'float64',
        'Cu_XRF_ng':'float64', 'Mg_XRF_ng':'float64', 'Mn_XRF_ng':'float64', 'Ni_XRF_ng':'float64', 'Sb_XRF_ng':'float64', 'Rb_XRF_ng':'float64', 'Sr_XRF_ng':'float64', 'Cd_XRF_ng':'float64',
        'Se_XRF_ng':'float64', 'Sn_XRF_ng':'float64',
        
        # Carbon columns (4)
        'BC_SSR_ug':'float64', 'BC_HIPS_ug':'float64', 'EC_FTIR_ug':'float64', 'EC_FTIR_MDL':'float64', 'OC_FTIR_ug':'float64','OC_FTIR_MDL':'float64','OM_FTIR_ug':'float64','OM_FTIR_MDL':'float64',
        
        # IC Nylon columns (13)
        'IC_F_ug_N':'float64', 'IC_Cl_ug_N':'float64', 'IC_NO2_ug_N':'float64', 'IC_Br_ug_N':'float64', 'IC_NO3_ug_N':'float64', 'IC_PO4_ug_N':'float64', 'IC_SO4_ug_N':'float64',
        'IC_Li_ug_N':'float64', 'IC_Na_ug_N':'float64', 'IC_NH4_ug_N':'float64', 'IC_K_ug_N':'float64', 'IC_Mg_ug_N':'float64', 'IC_Ca_ug_N':'float64',
        
        # Final columns (2)
        'method_index':'Int64', 'Flags':'object'
    }
    
    columns = list(dtype_mapping.keys())
    
    return  columns, dtype_mapping


def write_master(master_file, master_data):
    # Standard version
    # master_data.to_csv(master_file, index=False, na_rep='NaN', float_format='%.3f')
    
    columns, dtype_mapping = master_heading()
    for col, dtype in dtype_mapping.items():
        if dtype == 'float64':
            master_data[col] = master_data[col].round(3)
    master_data.to_csv(master_file, index=False, na_rep='NaN')        
            
    # 'Light' versions:
    # Define column groups
    ic_columns = [col for col in master_data.columns if col.startswith('IC_')]
    icp_columns = [col for col in master_data.columns if '_ICP_ng' in col]
    xrf_columns = [col for col in master_data.columns if '_XRF_ng' in col]
    carbon_columns = [col for col in master_data.columns if col.startswith(('BC_', 'EC_', 'OC_', 'OM_'))]

    # Create versions by excluding specified groups
    ic_only = master_data.drop(columns=icp_columns + xrf_columns + carbon_columns)  # No ICP/XRF/Carbon
    xrf_only = master_data.drop(columns=ic_columns + carbon_columns + icp_columns)  # No ICP/IC/Carbon
    carbon_only = master_data.drop(columns=ic_columns + icp_columns + xrf_columns)  # No ICP/IC/XRF

    # split "master_file" to get the directory and filename
    original_dir, filename = os.path.split(master_file)
    
    # Ensure sub-master folders exist  
    for subfolder in ['only_carbon', 'only_ic', 'only_xrf']:
        os.makedirs(os.path.join(original_dir, subfolder), exist_ok=True)

    # Save versions (CSV format shown - modify as needed)
    ic_only.to_csv(f'{original_dir}/only_ic/{filename}', index=False)
    xrf_only.to_csv(f'{original_dir}/only_xrf/{filename}', index=False)
    carbon_only.to_csv(f'{original_dir}/only_carbon/{filename}', index=False)
   

def create_master(master_file):
    # Define column names in exact order with validation
    columns, dtype = master_heading()

    # Create empty DataFrame with headers
    df = pd.DataFrame(columns=columns)

    # Write to CSV (will overwrite existing file)
    write_master(master_file,df)
    
    logging.info(f"Successfully created {master_file}.")
    return df


def backup_master(master_file):
    # This function creates a backup of the master file in a dated folder.
    # It also moves old backup folders to an Archives folder if they are older than 1 year.
    
    # Split the master file path into directory and filename
    master_path, fname = os.path.split(master_file)
    
    # Define the Backups and Archives directory paths
    backups_dir = os.path.join(master_path, 'Backups')
    archives_dir = os.path.join(backups_dir, 'Archives')
    
    current_date = datetime.now()
    
    # Process existing backup folders to move old ones to Archives
    for folder_name in os.listdir(backups_dir):
        folder_path = os.path.join(backups_dir, folder_name)
        # Check if it's a backup folder and not the Archives directory
        if os.path.isdir(folder_path) and folder_name.startswith('Backup_'):
            try:
                # Extract date from folder name
                date_str = folder_name.split('_')[1]
                folder_date = datetime.strptime(date_str, '%Y-%m-%d')
                # Check if the folder is older than 1 year
                if (current_date - folder_date) > timedelta(days=365):
                    dest_path = os.path.join(archives_dir, folder_name)
                    shutil.move(folder_path, dest_path)
            except (IndexError, ValueError):
                # Skip folders with invalid date formats
                continue
    
    # Create today's backup folder
    today_str = current_date.strftime('%Y-%m-%d')
    today_backup_dir = os.path.join(backups_dir, f'Backup_{today_str}')
    os.makedirs(today_backup_dir, exist_ok=True)
    
    # Check if the file is already in today's backup directory
    dest_file = os.path.join(today_backup_dir, fname)
    if not os.path.exists(dest_file):
        shutil.copy2(master_file, dest_file)
 
        
def read_master (master_file, create_new=False ):
    
    logging.info(f'reading {master_file}')
    
    if os.path.exists(master_file):
        # reading headings
        titles, dtype = master_heading()
        
        # make a copy just in case we need to go back 
        backup_master(master_file)
    
        df = pd.read_csv(master_file)
        
        # If there are new columns added to titles, master files will be updated:
        for col in titles:
            if col not in df.columns:
                df[col] = np.nan
        
        # sort columns order
        df = df[titles]
        
        # Get all columns that should be Int64 # NOTE: could be deleted in the future (when all files are in the right format)
        for col, dty in dtype.items():
            if dty == 'Int64':
                df[col] = pd.to_numeric(df[col], errors='coerce').round().astype('Int64')
            if dty == 'float64':
                df[col] = pd.to_numeric(df[col], errors='coerce').astype('float64')
                
        write_master(master_file,df)
        
        # Read raw data with type preservation
        df = pd.read_csv(master_file,  dtype=dtype)
        df = df.replace('NaN', np.nan)
        df['Flags'] = df['Flags'].fillna('') 
            
                
        if df.empty:
            return df

        # Formate FilterIDs to SITE_XXXX
        df['FilterID'] = df.iloc[:, 0].astype(str).apply(
            lambda x: x[:5] + '0' + x[5:8] if len(x) < 9 else x
        )

        # # Convert numeric columns to NaN-filled floats
        # numeric_cols = [col for col in titles if col not in [
        #     'FilterID', 'Filter_Barcode', 'CartridgeID', 'LotID', 'Flags', 'projectID'
        # ]]
        # df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce')

        # # LotID should be str, not int or float
        # df['LotID'] = df['LotID'].apply(
        #     lambda x: 'NaN' if pd.isna(x) else str(int(x)) if isinstance(x, float) else str(x)
        # )
        
        # Reorder columns and filter to expected structure
        master_data = df[titles].copy()
                   
    else:
        if create_new is True:
            logging.info(f'{master_file} not found. Creating one.')
            master_data = create_master(master_file)
        else:
            logging.info(f'{master_file} not found.Skipping.')
            master_data = None
            
        
    return master_data

