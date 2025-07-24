# PURPOSE: ------------------------------------------------------------
# 1. read and perform QA/QC on IC data straight from the IC
# system with no pre-processing; 
# 2. put IC data and relevant flags in the Master data file for each 
# SPARTAN site that has data contained in the IC file.
# 3. move IC files from "Raw" to a relevant archive folder
#
# FILE REQUIREMENTS: ------------------------------------------------------
# directly from IC without any manual processing
# 
# Input: -----------------------------------------------------------
# Analysis_Data/Ion_Chromatography/NEW/
# Analysis_Data/Master_files 
# Analysis_Data/E-Logs/IC E-Log.xlsx
# Analysis_Data/Ion_Chromatography/Syringe_filter_contamination/cation/IC_low_calibration_xxx
# Analysis_Data/Ion_Chromatography/MDL
# 
# Output: ----------------------------------------------------------
# Analysis_Data/Master_files
# Analysis_Data/Archived_Filter_Data/IC_data/
#
# Suppliment: ------------------------------------------------------
# Acceptance_Testing_Summary_WashU.xlsx
# 
# To re-process IC files: ------------------------------------------- 
# Turn Redo_All_Archieve to 1
# 
# Created by: --------------------------------------------------------------
# Haihui Zhu, Crystal Weagle, Emme Le Roy, and Chris Oxford 
# June-2025

import pandas as pd
import os
import numpy as np
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
direc = su.find_root_dir(debug_mode)
direc_new = os.path.join(direc, "Analysis_Data/Ion_Chromatography/NEW")
direc_master = os.path.join(direc, "Analysis_Data/Master_files")
direc_archive = os.path.join(direc, "Analysis_Data/Archived_Filter_Data/IC_data/")
log_file_path = os.path.join(direc,'Public_Data/Data_Processing_Records/IC', datetime.now().strftime('%Y-%m-%d-%H%MSS') + '_IC_Record.txt')

# # TEST 
# direc_new = os.path.join('./release_test',"IC_NEW")
# direc_archive = os.path.join('./release_test',"IC_archive/")
# direc_master = os.path.join('./release_test',"Master_files")
# log_file_path = os.path.join('./release_test', datetime.now().strftime('%Y-%m-%d-%H%MSS') + '_IC_Record.txt')

# Log file setup
su.setup_logging(log_file_path)  # Initialize logging
logging.info("D started")

# Read the site info table
SiteInfoFile = f'{direc}/Site_Sampling/Site_details.xlsx'
site_details = pd.read_excel(SiteInfoFile)

# reading elog (contains extraction volumes for individual filters)
elogfname = f'{direc}/Analysis_Data/E-Logs/IC E-Log.xlsx' 

# read the calibration curves for low concentrations
low_conc_curve_file = f'{direc}/Analysis_Data/Ion_Chromatography/Syringe_filter_contamination/low_conc_curve_all.xlsx'

# read MDLs 
# only use teflon since MDL are only used for quality checks - if H2O samples has high ion concentration
# sometimes teflon and nylong are mixed in one batch, can't separate H2O samples for two MDLs
mdl_path = f'{direc}/Analysis_Data/Ion_Chromatography/IC_Teflon detection limit.xlsx'
mdl = pd.read_excel(mdl_path, nrows=2,index_col=0) 

# ============================
# Define Functions & constants
# ============================
# Correlation coefficient of a calibration curve should be higher than this value for the IC data to be consider valid.
CorrCoeff_threshold = 99.5;

# all ions:
ion_abb = {
    'Fluoride': 'F',
    'Chloride': 'Cl',
    'Nitrite': 'NO2',
    'Bromide': 'Br',
    'Nitrate': 'NO3',
    'Phosphate': 'PO4',
    'Sulfate': 'SO4',

    'Lithium': 'Li',
    'Sodium': 'Na',
    'Ammonium': 'NH4',
    'Potassium': 'K',
    'Magnesium': 'Mg',
    'Calcium': 'Ca',
    
    'HMS': 'HMS'
}
# concentration of ions in a standard sample should be within +/-10 % of the following values
ion_QC_conc = {'F':0.25,'Cl':0.375 ,'NO2':1.25 ,'Br':1.25 ,'NO3':1.25 ,'PO4':1.875 ,'SO4':1.875, \
    'Li':0.25 ,'Na':1.00 ,'NH4':1.25,'K':2.5,'Mg':1.25,'Ca':2.5};

def find_row_col(df, pattern):
    mask = df.applymap(lambda x: pattern in str(x).strip().lower())
    
    stacked = mask.stack()
    true_entries = stacked[stacked]  # Filter to True values
    if true_entries.empty:
        raise ValueError(f"{pattern} not found in IC input data.")
    first_index_tuple = true_entries.index[0]  # Get first (row, column) MultiIndex tuple
    row_idx = first_index_tuple[0]  # First level = row index
    col_idx = first_index_tuple[1]  # Second level = column index
    return row_idx, col_idx
    
def read_calib_curve(file_path, sheet_name):
    raw_calibration = pd.read_excel(file_path, sheet_name=sheet_name)
    
    # check if the data is from WashU. Stop when it isn't
    # Locate the cell containing 'Operator:'
    row_idx, col_idx = find_row_col(raw_calibration, 'operator')
    next_col = raw_calibration.columns[raw_calibration.columns.get_loc(col_idx) + 1]
    operator_value = raw_calibration.at[row_idx, next_col]
    if operator_value != 'SPARTAN':
        print(f'file operator: {operator_value}')
        exit()
            
    # find the date
    row_idx, col_idx = find_row_col(raw_calibration, 'date')
    next_col = raw_calibration.columns[raw_calibration.columns.get_loc(col_idx) + 1]
    extraction_date = raw_calibration.at[row_idx, next_col]
    
    # find if cubic curve exists in this sheet:
    cubic_idx, _ = find_row_col(raw_calibration, 'calibration summary')
    # cubic_idx = raw_calibration.index[raw_calibration.apply(lambda row: 'Calibration Summary' in str(row.values), 1).min()][0] 
    cubic = raw_calibration.iloc[cubic_idx, 1] 
            
    # Find start and end indices
    mask = raw_calibration.apply(lambda row: any('Peak Name' in str(cell) for cell in row), axis=1)
    if not mask.any():
        mask = raw_calibration.apply(lambda row: any('Name' in str(cell) for cell in row), axis=1)
        idx = mask.idxmax()  # First row where 'Name' is found
        # Replace 'Name' cells with 'Peak Name'
        raw_calibration.loc[idx] = raw_calibration.loc[idx].apply(
            lambda cell: 'Peak Name' if 'Name' in str(cell) else cell
        )
    start_idx = mask.idxmax()
    
    mask = raw_calibration.apply(lambda row: any('AVERAGE' in str(cell) for cell in row), axis=1)
    if not mask.any():
        end_idx = start_idx+11
    else: 
        end_idx = mask.idxmax()-1
        
    
    # Extract and clean data
    calib_curves = raw_calibration.loc[start_idx:end_idx].copy()
    calib_curves.columns = calib_curves.iloc[0].values # Set row6 as column names
    calib_curves = calib_curves.iloc[1:]         # Remove header row
    calib_curves.columns = calib_curves.columns.str.rstrip() # remove tailing space in column name 
    
    # Filter out non-data rows (NaN or 'CD' in first column)
    calib_curves = calib_curves[calib_curves.iloc[:, 0].notna() & (calib_curves.iloc[:, 0] != 'CD')].reset_index(drop=True)

    return calib_curves, extraction_date, cubic

def read_panel(df, amount_idx, sum_idx):

    # Calculate start and end indices
    start_idx = amount_idx + 4
    end_idx = sum_idx

    # Extract the data block
    data_block = df.iloc[start_idx:end_idx,:].copy()

    # Get the cation headers from the row above the data block
    cation_headers = df.iloc[amount_idx + 3, 2:10].tolist()

    # Create new column names
    new_columns = ['index', 'samplename'] + cation_headers

    # Select only the relevant columns (first 8 columns)
    data_block.columns = new_columns

    # Convert numeric columns to float, handle 'n.a.' as NaN
    for col in cation_headers:
        data_block[col] = pd.to_numeric(data_block[col], errors='coerce')

    # Reset index for cleaner output
    return data_block.reset_index(drop=True)

def read_conc_and_signal(file_path, sheet_name):
    df = pd.read_excel(file_path, sheet_name=sheet_name)
    
    # Find the row index with 'Amount'
    amount_idx, _ = find_row_col(df, 'amount')

    # Find the index of first row with 'Sum:'
    sum_idx, _ = find_row_col(df, 'sum')
    
    # Read the amount (concenctration) panel:
    conc = read_panel(df, amount_idx, sum_idx)
    
    # Find the row index with 'Area'
    area_idx, _ = find_row_col(df, 'area')

    # The next 'sum' should be somewhere near area_idx+sum_idx:
    sum_idx2, _ = find_row_col(df.iloc[area_idx:area_idx+sum_idx,:], 'sum')
    
    # Read the area panel:
    area = read_panel(df, area_idx, sum_idx2)
    
    return conc, area

def calc_conc(calib, calib_low_conc, selec_area, apply_low_conc_curve, cubic):
    
    """Solve for concentration using calibration curve parameters"""
    if 'Cubic' in calib['Cal.Type']: # this is not a robust check since some old files has 'cubic' in 'cal type' but don't provide cubic coeff

        if  cubic == 1: 
            if abs(selec_area)>=0:
                # Cubic equation: area = Offset + Slope*conc + Curve*conc² + Cubic*conc³
                # Rearrange to: Cubic*conc³ + Curve*conc² + Slope*conc + (Offset - area) = 0
                coeffs = [calib['Cubic'], calib['Curve'], calib['Slope'], (calib['Offset'] - selec_area)]
                # Find roots of cubic polynomial
                roots = np.roots(coeffs)
                
                # Filter real positive roots
                real_roots = [r.real for r in roots if np.isreal(r)] # do not need to find the positive root. For H2O samples, conc can be nagative
                if  len(real_roots) == 1:
                    conc = real_roots[0] 
                else:
                    real_positve_roots = [r for r in real_roots if r>0]
                    
                    if len(real_positve_roots) > 0:
                        conc = min(real_positve_roots)  # pick the minimum positive value
                    else: 
                        conc = max(real_roots) # pick the maximum nagative value
            else:
                conc = np.nan
        else:
            conc = np.nan
                     
    elif 'Lin' in calib['Cal.Type']:
        # Linear equation: area = Offset + Slope * conc
        conc = (selec_area - calib['Offset']) / calib['Slope']
        
        if apply_low_conc_curve == 1: # can apply low_conc curve for filtees with tiny areas
            if conc < calib_low_conc['conc'].max()*0.85: # notably lower than the high end of the curve
                conc = np.interp(selec_area, calib_low_conc['signal'], calib_low_conc['conc'])
                       
    return conc
 
 
def get_volume_mapping(elogfname, site):
    excel = pd.ExcelFile(elogfname, engine='openpyxl')
    sheet_name = next((s for s in excel.sheet_names if site in s), None)
    if sheet_name is None:
        raise ValueError(f"No sheet found for site '{site}' in elog")
    
    elog = pd.read_excel(elogfname, sheet_name=sheet_name)
    elog.columns = elog.columns.str.rstrip()
    volume_map = elog.set_index('Analysis ID')['Volume (mL)'].to_dict()
    
    # Add reformatted entries (remove leading zeros after '-')
    for fid in list(volume_map.keys()):
        if isinstance(fid, str):
            parts = fid.split('-', 1)
            if len(parts) == 2 and parts[1].startswith('0'):
                reformatted = f"{parts[0]}-{parts[1][1:]}"
                if reformatted not in volume_map:  # Avoid overwriting
                    volume_map[reformatted] = volume_map[fid]
    return volume_map
  
      
def update_flag(tflag, newstr):
    if tflag == '':
        udpatedflag = tflag
    else:
        udpatedflag = [segment.strip() + newstr for segment in tflag.split(';')]
        udpatedflag = '; '.join(udpatedflag)
    return udpatedflag
    



    
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

# move archived files to new (if enabled)
if Redo_All_Archive:
    logging.info(f'Copying archived IC report for reprocessing')
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
                # logging.info(f"Skipping file in skiplist: {file}") # not really need this printout
                continue

            source = os.path.join(site_dir, file)
            su.copy_file(source, direc_new)
            
    logging.info(f'Done restoring archived data.\n')

files = su.get_files(direc_new)
if not files:
    logging.info("No new files to process.")
    exit()

for file in files:
    file_path = os.path.join(direc_new, file)
    logging.info(f"Processing file: {file}")

    # Read sheet names
    excel = pd.ExcelFile(file_path)
    sheets = excel.sheet_names
    sheet_summary = next((s for s in sheets if 'Summary' in s), None)
    sheet_calibration = next((s for s in sheets if 'Calibration' in s), None)
    
    if not sheet_summary or not sheet_calibration:
        logging.error(f"Missing required sheets in {file}")
        continue

    # get calibration curves 
    calib_curve, extraction_date, cubic = read_calib_curve(file_path, sheet_calibration)
    
    # Read IC raw data
    all_conc, all_area = read_conc_and_signal(file_path, sheet_name=sheet_summary)
    
    # find if there is a peak named HMS, if yes, do not apply the low concentration curve
    if 'HMS' in calib_curve['Peak Name'].str.strip().values:
        apply_low_conc_curve = 0
        calib_curve = calib_curve[calib_curve['Peak Name'] != 'HMS']  # remove HMS from the calibration curve
    else:
        apply_low_conc_curve = 1 # use to indicate if low concentration curve applies

    tflag = '' # initiate flag
    # Loop through ions in the calib curve:
    for calib_idx, calib in  calib_curve.iterrows():
        ion_full = calib['Peak Name'] # find the full name of this ion
        if ion_full == 'Fluoride':
            logging.info(f'skipping {ion_full} data')
            # remove flouride in calib_curve
            calib_curve = calib_curve[calib_curve['Peak Name'] != 'Fluoride'] 
            continue # no need to calc concentration for flouride
        
        logging.info(f'processing {ion_full} data')
        
        ion = ion_abb[ion_full] # find the abbreviation of this ion
        
        # skip if no calibration curve:
        if pd.isna(calib['Coeff.Det.']) or calib['Coeff.Det.']=='n.a.':
            logging.warning(f'No calibration for {ion}')
            continue
        
        # read low conc calibration curve:
        if ion == 'NH4': # NH4 do not have a low conc curve
            calib_low_conc = None
        else:
            calib_low_conc = pd.read_excel(low_conc_curve_file, sheet_name=ion)
    
        #### calculate ion concenration ###
        select_area = all_area[ion_full] # select the column for this ion
        for idx, area_val in select_area.items():
            # replace concentration with our calculated values
            all_conc.loc[idx, ion_full]= calc_conc(calib, calib_low_conc, area_val, apply_low_conc_curve, cubic)
    
        ### QC and Flags ###
        # check correlation of calibration curve, flag if lower than CorrCoeff_threshold
        if calib['Coeff.Det.'] < CorrCoeff_threshold:
            logging.warning(f'{ion} curve correlation coefficient lower than {CorrCoeff_threshold}')
            tflag = su.add_flag(tflag, f'S-{ion}')
        
        # check if the QC STD 1.25 concentration is within +/-10 % of the specified concentration for each anion 
        # (data before May 2018 will not have these)
        qc_samples = all_conc.loc[all_conc['samplename'].str.strip().str.upper().str.contains('QC 1.25'), :]
        qc_samples = qc_samples[ion_full]
        qc_avg = qc_samples.mean() 
        if qc_avg > 1.1*ion_QC_conc[ion] or qc_avg < 0.9*ion_QC_conc[ion] : 
            logging.warning(f'Avg QC STD 1.25 for {ion} exceeds +/- 10 %% of specified concentration')
            tflag = su.add_flag(tflag, f'qc1-{ion}') 
    
        # H2O should be within 10xMDL (use areas here)
        h2o_samples = all_area[all_area['samplename'].str.strip().str.upper() == 'H2O']
        h2o_samples = h2o_samples[ion_full]
        h2o_avg = h2o_samples.mean() 
        if extraction_date < datetime(2021,11,25):
            mdl_values = mdl.loc['Before_2021_Nov25',f'IC_area_{ion}']
        else:
            mdl_values = mdl.loc['After_2021_Nov25',f'IC_area_{ion}']
        
        if h2o_avg > (10 * mdl_values): 
            logging.warning(f'Average H2O for {ion} is > 10*MDL')
            tflag = su.add_flag(tflag, f'qc2-{ion}') 
        
        
    # find sites in this batch:
    pattern = r'^([A-Za-z]{4})[-]\d{4}(?:\w+)?$'
    site_samples = all_conc[all_conc['samplename'].str.match(pattern, na=False)]
    site_unique = site_samples['samplename'].str.extract(pattern, expand=False).str.upper().unique()
    # prepare to move IC file to archive
    file_move_status = {}
            
    # loop through all sites:
    for site in site_unique:
        if site in site_details['Site_Code'].tolist():
            # Precompute volume mapping
            volume_map = get_volume_mapping(elogfname, site)
            
            master_file = os.path.join(direc_master, f'{site}_master.csv')
            master_data = su.read_master(master_file)
            
            # Filter out lab blanks and non-site samples
            valid_samples = site_samples[
                site_samples['samplename'].str.contains(site) &
                ~site_samples['samplename'].str.contains('LB') &
                ~site_samples['samplename'].str.lower().str.contains('lab')
            ].copy()

            # Classify filters and match volumes
            valid_samples['filter_type'] = np.where(
                valid_samples['samplename'].str.endswith('N'), 'nylon', 'teflon'
            )
            valid_samples['volume'] = valid_samples['samplename'].map(volume_map)
            if valid_samples['volume'].isna().any():
                missing = valid_samples.loc[valid_samples['volume'].isna(), 'samplename'].tolist()
                raise ValueError(f"Filters not found in elog: {missing}")

            # compute masses for all ions
            ion_cols = calib_curve['Peak Name'].tolist()
            for ion_full in ion_cols:
                ion = ion_abb[ion_full]
                valid_samples[f'mass_{ion}'] = valid_samples['volume'] * valid_samples[ion_full]
            
            
            # Process nylon filters (grouped by truncated ID)
            nylon_samples = valid_samples[valid_samples['filter_type'] == 'nylon'].copy()
            nylon_samples['FilterID'] = nylon_samples['samplename'].str[:9]

            if not nylon_samples.empty:
                # Aggregate ion masses for each truncated FilterID
                mass_cols = [f'mass_{ion_abb[ion]}' for ion in ion_cols]
                nylon_agg = nylon_samples.groupby('FilterID')[mass_cols].first().reset_index()
                
                # Merge to update master_data (nylon masses)
                for ion_full in ion_cols:
                    ion = ion_abb[ion_full]
                    
                    # Prepare the update data: rename column and index by FilterID
                    update_col = f'IC_{ion}_ug_N'
                    nylon_update = nylon_agg[['FilterID', f'mass_{ion}']].rename(
                        columns={f'mass_{ion}': update_col}
                    ).set_index('FilterID')

                    # Set FilterID as index in master_data (temporarily)
                    master_data.set_index('FilterID', inplace=True)

                    # Update only matching rows
                    master_data.update(nylon_update)

                    # Restore FilterID as a column if needed
                    master_data.reset_index(inplace=True)
            
                
                # Update flags for nylon
                nylon_filters = nylon_samples['FilterID'].unique()
                updatedflag = update_flag(tflag, '-T')
                for nidx in nylon_filters:
                    # Ensure FilterID exists in master_data
                    # If not, add a new entry with NaN mass and empty fields
                    # This is to ensure that we do not miss any FilterID in the master file
                    if nidx not in master_data['FilterID'].values:
                        print(f"[INFO] Adding missing FilterID: {nidx}")
                        master_data = add_entry(master_data, nidx, mass=np.nan, mass_type=np.nan, sampling_mode=np.nan,
                                                barcode=np.nan, cartid=np.nan, lotid=np.nan) #Script D just adds filterID, other data is added by A2 and B

                
                    master_data.loc[master_data['FilterID']==nidx, 'Flags'] = su.add_flag(master_data.loc[master_data['FilterID']==nidx, 'Flags'].values[0], updatedflag)
                

            # Process teflon filters
            teflon_samples = valid_samples[valid_samples['filter_type'] == 'teflon'].copy()
            teflon_samples['FilterID'] = teflon_samples['samplename']

            if not teflon_samples.empty:
                # Update master_data (teflon masses)
                for ion_full in ion_cols:
                    ion = ion_abb[ion_full]
                    
                    # Prepare the update data: rename column and index by FilterID
                    update_col = f'IC_{ion}_ug_T'
                    teflon_update = teflon_samples[['FilterID', f'mass_{ion}']].rename(
                        columns={f'mass_{ion}': update_col}
                    ).set_index('FilterID')

                    # Set FilterID as index in master_data (temporarily)
                    master_data.set_index('FilterID', inplace=True)

                    # Update only matching rows
                    master_data.update(teflon_update)

                    # Restore FilterID as a column if needed
                    master_data.reset_index(inplace=True)
            
                    # mass_map = teflon_samples.set_index('FilterID')[f'mass_{ion}'].to_dict()
                    # master_data[f'IC_{ion}_ug_T'] = master_data['FilterID'].map(mass_map).fillna(master_data[f'IC_{ion}_ug_T'])
                    
                        
                # Update flags for teflon
                teflon_filters = teflon_samples['FilterID'].unique()
                updatedflag = update_flag(tflag, '-T')
                for tidx in teflon_filters:
                    if tidx not in master_data['FilterID'].values:
                        print(f"[INFO] Adding missing FilterID: {tidx}")
                        master_data = add_entry(master_data, tidx, mass=np.nan, mass_type=np.nan, sampling_mode=np.nan,
                                                barcode=np.nan, cartid=np.nan, lotid=np.nan) #Script D just adds filterID, other data is added by A2 and B

                    master_data.loc[master_data['FilterID']==tidx, 'Flags'] = su.add_flag(master_data.loc[master_data['FilterID']==tidx, 'Flags'].values[0], updatedflag)

            # Sort all rows by FilterID before saving
            master_data.sort_values(by='FilterID', inplace=True) 

            # write to master file
            su.write_master(master_file, master_data)
            logging.info(f'Done writing IC data to {site} master file.')
            
            # Archive IC data files
            if  site not in file_move_status:
                destination = f'{direc_archive}/{site}/'
                file_move_status[site] = su.copy_file(file_path, destination)

    if any(file_move_status.values()):
        os.remove(file_path)
        logging.info('Moved IC file to archive.\n')
 
logging.info("Complete!")
