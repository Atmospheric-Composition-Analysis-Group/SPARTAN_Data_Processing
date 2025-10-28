########################### CODE DESCRIPTION ################################
# PURPOSE: ------------------------------------------------------------
# Read and perform QA/QC on flow data for each cartridge
# Assign mass type, method code, sampling hours to each filter
# 
# Input: -------------------------------------------------------------------
# Site_Sampling/$SiteCodes_$SiteCities/Cartridge_Data/$SiteCodes_dates_flows.xlsx"
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
# Crystal Weagle, 4 October 2018
# Revised by Haihui Zhu, April 2025

from datetime import datetime
import pandas as pd
import numpy as np
import os
import logging
import spt_utils as su

# ============================
# Initial Setup and Configuration
# ============================
Clean_Up_Master = 0  # set to 1 if want to reprocess all dates and volumes

# Setup directories 
debug_mode = 0
direc = su.find_root_dir(debug_mode)

# ============================
# Directory & File Paths
# ============================
direc_sampling = os.path.join(direc, 'Site_Sampling')
direc_master = os.path.join(direc, 'Analysis_Data/Master_files')
log_file_path = os.path.join(direc, 'Public_Data/Data_Processing_Records/Flow_Processing', 
                          f"{datetime.now().strftime('%Y-%m-%d-%H%M%S')}_Flow_Processing_Record.txt")

# TEST dirs
# direc_sampling = os.path.join('./B_test', 'Site_Sampling')
# direc_master = os.path.join('./release_test', 'Master_files')
# log_file_path = os.path.join('./release_test',  
#                           f"{datetime.now().strftime('%Y-%m-%d-%H%M%S')}_Flow_Processing_Record.txt")

# Log file setup
su.setup_logging(log_file_path)  # Initialize logging
logging.info("B started")

# Read the site info table
SiteInfoFile = f'{direc}/Site_Sampling/Site_details.xlsx'
site_details = pd.read_excel(SiteInfoFile)

# ============================
# Define Functions
# ============================
# None

# ============================
# Start processing 
# ============================
# Note before processing:
# Mass type code
# 0 = cartridge (traveling) blank
# 1 = PM2.5
# 2 = PM10
# 3 = PMcoarse
# 4 = saturated nuclepore, filter measurements invalid [Not found at any site. The data might be too old to be used.]
# 5 = negative mass, filter measurements invalid
# 6 = invalid flow, filter measurements invalid


# Read dates, hours, and volumes. Determine masstype and methods

# # Uncomment this if start from a specific site: 
# site_list = site_details['Site_Code'].tolist()
# start_site = 'IDBD' 
# start_idx = site_list.index(start_site)
# for idx in range(start_idx, len(site_list)):
#     site = site_details['Site_Code'][idx]
    
# Uncomment this if start from begining 
for idx, site in enumerate(site_details['Site_Code']):

    filename = os.path.join(direc_sampling, 
        f"{site}_{site_details['City'][idx]}/Cartridge_Data/{site}_dates_flows.xlsx")
    
    # Check if a site flow sheet exists
    if os.path.exists(filename):  # file found
        
        logging.info(f'Reading file for {site}')
        flow_data = pd.read_excel(filename)
        flow_data = flow_data.replace('NaN', np.nan)
        flow_data = flow_data.replace('NAN', np.nan)
        flow_filterid_origial = flow_data['Filter_ID'] # Will need it later
         
        # ---- make sure filter ID is in standard format (XXXX_NNNN) ----
        flow_data['Analysis_ID'] = su.filter_id_format(flow_data['Analysis_ID'])
        
        # correct heading name: Flow_enternal_end shoule be Flow_external_end
        if 'Flow_enternal_end' in flow_data.columns:
            flow_data.rename(columns={'Flow_enternal_end': 'Flow_external_end'}, inplace=True)
        if 'Hours_sampled' in flow_data.columns:
            flow_data.rename(columns={'Hours_sampled': 'hours_sampled'}, inplace=True)

        # ---- read master files ----
        master_file = os.path.join(direc_master, f"{site}_master.csv")
        master_data = su.read_master(master_file)
        
        if master_data is None:  # skip if no filter added (A2 should add filter info to master before running B)
            continue

        # when Clean Up Master = 1, NaN all relevant data so that they can be added later (i.e. updated)
        if Clean_Up_Master == 1:
            master_data['Volume_m3'] = np.nan
            master_data['hours_sampled'] = np.nan
            master_data.loc[master_data['Mass_type'] != 5, 'Mass_type'] = np.nan
            master_data['method_index'] = np.nan
            
            for col in master_data.columns:
                if 'start' in col or 'stop' in col:
                    master_data[col] = np.nan


        for idx, tfilter in enumerate(master_data['FilterID']):
            
            if pd.isna(master_data.loc[idx, 'hours_sampled']): # sampling hour not added, search in flow_data
                
                mask_flow = flow_data['Analysis_ID']== tfilter  # find the corresponding filter position in flow_data

                if mask_flow.any(): # if filter found in flow_data 
                    
                    tflow_SSe_ID = flow_data.loc[mask_flow, 'SSe_ID'].values[0]

                    # Minimal strict check: stop if any alphabetic character is present in SSe_ID
                    if isinstance(tflow_SSe_ID, str) and any(c.isalpha() for c in tflow_SSe_ID):
                        raise ValueError(f"{site}: SSe_ID must be numeric; got {tflow_SSe_ID!r} for FilterID={tfilter}")

                    
                    tMaster_LotID     = master_data['LotID'][idx]
                    tMaster_ProjectID = master_data['projectID'][idx]
                    
                    # Mass Type Criterion 1: 
                    # mass_type = 3 (PM coarse) for nucleopore filters 
                    if flow_filterid_origial[mask_flow].values[0][-1] == 'N':
                        master_data.loc[idx, 'Mass_type'] = 3 

                    # Mass Type Criterion 2: 
                    # mass_type = 1 (PM25) for external start flow rate = 5
                    # mass_type = 2 (PM10) for external start flow rate = 1.5
                    col_ext_start = 1
                    targets = [5, 1.5, 0]  # possible target flow rates
                    types = [1, 2, 0]
                    if pd.isna(master_data.loc[idx, 'Mass_type']):
                        exFlow = flow_data.loc[mask_flow, 'Flow_external_start'].values[0]
                        fidx = np.argmin(np.abs(np.array(targets) - exFlow))  # find which flow rate is closer  
                        master_data['Mass_type'][idx] = types[fidx]  # sets mass type

                    # Set Method Index
                    method_index = np.nan
                    if tflow_SSe_ID / 5000 < 1:  # SS4
                        if site == 'SGSU':  # OLD SGSU sampler is SS4 but has target flow rate of 5 LPM
                            target = 5  
                            method_index = 2
                        else:
                            target = 4
                            method_index = 1
                        
                        r_SF = flow_data.loc[mask_flow, 'Flow_external_start'].values[0] * np.array([0.9, 1.1])  # 10% range of external start flow rate
                        r_target = target * np.array([0.9, 1.1])  # range of acceptable start flow rates
                    
                    elif master_data.loc[idx, 'Mass_type'] == 2: # PM10
                        
                        r_SF = flow_data.loc[mask_flow, 'Flow_external_start'].values[0] * np.array([0.8, 1.2])  # range of acceptable flow rates
                        r_target = np.array([1.2, 1.8]) 
                        if tMaster_LotID == 'Unknown':
                            method_index = 1
                        elif tMaster_ProjectID == 'S':
                            method_index = 2
                        elif tMaster_ProjectID == 'M':
                            method_index = 3
                            
                    elif master_data.loc[idx, 'Mass_type'] == 1: # PM25
                        
                        target = 5
                        r_SF = flow_data.loc[mask_flow, 'Flow_external_start'].values[0] * np.array([0.9, 1.1])  # 10% range of external start flow rate
                        r_target = target * np.array([0.9, 1.1])  
                        if tMaster_LotID == 'Unknown':
                            method_index = 3
                        elif tMaster_ProjectID == 'S':
                            method_index = 4
                        elif tMaster_ProjectID == 'M':
                            method_index = 5
                            
                    elif np.isnan(tflow_SSe_ID):
                        r_target = np.array([0, 0])
                        r_SF = np.array([0, 0])
                        method_index = np.nan
                        logging.info(f'No SSe ID found for flows of {tfilter}')

                    # Mass Type Criterion 3: 
                    # mass_type = 0 (field blank) if year/volume/any rate == 0
                    if flow_data.loc[mask_flow,'start_year'].values[0] == 0 and flow_data.loc[mask_flow,'volume_m3'].values[0] == 0:
                        master_data.loc[idx, 'Mass_type'] = 0
                        
                    elif master_data.loc[idx, 'Mass_type'] == 5 or master_data.loc[idx, 'Mass_type'] == 3:
                        pass
                    
                    # Mass Type Criterion 4: 
                    # mass_type = 6 (invalid) if volume = nan, or flow rate out of target range
                    elif np.isnan(flow_data.loc[mask_flow,'volume_m3'].values[0]):  # adding a check: if volume is nan, ignore flow rates and assign mass type as invalid.
                        master_data.loc[idx, 'Mass_type'] = 6
                        
                    else:
                        # CHECK: External Start Flow within 10% of Target Flow Rate
                        if flow_data.loc[mask_flow, 'Flow_external_start'].values[0]  > r_target[1] or \
                            flow_data.loc[mask_flow, 'Flow_external_start'].values[0]  < r_target[0]:
                                
                            logging.info(f'WARNING: external start flow rate not within 10% of target flow rate for {tfilter}')
                            master_data.loc[idx, 'Mass_type'] = 6 

                        # CHECK: External End Flow within 10% of External Start Flow Rate
                        if flow_data.loc[mask_flow, 'Flow_external_end'].values[0]  > r_SF[1] or \
                            flow_data.loc[mask_flow, 'Flow_external_end'].values[0]  < r_SF[0]:
                                
                            logging.info(f'WARNING: external end flow rate not within 10% of external start flow rate, indicates clogging for {tfilter}')
                            master_data.loc[idx, 'Mass_type'] = 6 

                        # CHECK: internal start flow rate must be within 10 % of external start flow rate
                        if flow_data.loc[mask_flow, 'Flow_internal_start'].values[0]  != -899:  # internal flow rate for CASH never work, skip CASH to do the check
                            if flow_data.loc[mask_flow, 'Flow_internal_start'].values[0]  > r_SF[1] or \
                                flow_data.loc[mask_flow, 'Flow_internal_start'].values[0]  < r_SF[0]:
                                master_data.loc[idx, 'Mass_type'] = 6 
                                logging.info(f'WARNING: internal start flow rate invalided for {tfilter} reported by instrument')
                                
                         
                    # Assign values using .loc for unambiguous indexing
                    master_data.loc[idx, 'hours_sampled'] = flow_data.loc[mask_flow, 'hours_sampled'].values[0]
                    master_data.loc[idx, 'start_year'] = flow_data.loc[mask_flow, 'start_year'].values[0]
                    master_data.loc[idx, 'start_month'] = flow_data.loc[mask_flow, 'start_month'].values[0]
                    master_data.loc[idx, 'start_day'] = flow_data.loc[mask_flow, 'start_day'].values[0]
                    master_data.loc[idx, 'start_hour'] = flow_data.loc[mask_flow, 'start_hour'].values[0]
                    master_data.loc[idx, 'stop_year'] = flow_data.loc[mask_flow, 'stop_year'].values[0]
                    master_data.loc[idx, 'stop_month'] = flow_data.loc[mask_flow, 'stop_month'].values[0]
                    master_data.loc[idx, 'stop_day'] = flow_data.loc[mask_flow, 'stop_day'].values[0]
                    master_data.loc[idx, 'stop_hour'] = flow_data.loc[mask_flow, 'stop_hour'].values[0]
                    master_data.loc[idx, 'Volume_m3'] = flow_data.loc[mask_flow, 'volume_m3'].values[0]  # Writes volumes to master file
                    master_data.loc[idx, 'method_index'] = method_index
                    
                    # Adding flags from flow file
                    if not isinstance(flow_data.loc[mask_flow,'flags'].values[0], (float, np.floating)):
                        master_data.loc[idx, 'Flags'] = su.add_flag(master_data.loc[idx, 'Flags'], flow_data.loc[mask_flow,'flags'].values[0])
                    
                    # Adding additional flags if mass_type is 5 or 6: 
                    if master_data.loc[idx, 'Mass_type']== 5:
                        master_data.loc[idx, 'Flags'] =  su.add_flag(master_data.loc[idx, 'Flags'], 'Negative Mass')
                    elif master_data.loc[idx, 'Mass_type'] == 6:
                        master_data.loc[idx, 'Flags'] =  su.add_flag(master_data.loc[idx, 'Flags'], 'Invalid Flow')
                    elif master_data.loc[idx, 'Mass_type'] == 1 | master_data.loc[idx, 'Mass_type'] == 2:
                        logging.info(f'updating {tfilter}')
                        

        # Write to Master Data File  # NOTE_HAIHUI: Relace with df.to_csv
        su.write_master(master_file, master_data)
        
    else:
        logging.info(f'No flow file exists for {site}')
        print(f'{filename}')
        
logging.info('Finished reading flow dates')
