########################### CODE DESCRIPTION ################################
# PURPOSE: ------------------------------------------------------------
# Produce public files after applying apply data validations on master data files 
# STEPs include: 
# 1. Apply a default blank correction for metals determined by ICP-MS
# 2. Search for any flag or signs of invalid sampling (pinholes, bad flow
# rate, etc) and remove the filter for further analysis 
# 3. Calculate BC using SSR when HIPS measurement is not available. 
# 4. Uses flags from IC and the 'Bad_sodium' list to remove invalid IC measurements
# 5. Check sum of measured species vs. mass collected and mark filter as invalid if sum of species is greater than collected mass
# 6. Reconstruct fine mass using measured constituents
# 7. Determine filter-specific kappa values based on measured constituents
# 8. Generate data products for research and posting on SPARTAN website.
#    Data products generated include:
#     - all measured components (dry, 0% RH), with any necessary corrections applied for PM2.5, PM10, and when applicable, PMcoarse
#     - reconstructed fine mass with residual matter and PBW at 35 %RH
# 
# Input: -------------------------------------------------------------------
# Analysis_Data/Master_files/$SiteCode_master.csv
#
# Output: -------------------------------------------------------------------
# everything under:
# Public_Data/Chemical_Filter_Data/ 
# Public_Data/RCFM/
# 
# Suppliment: ------------------------------------------------------
# Site_Sampling/Site_details.xlsx
# 
# Version Update Notes (Older data versions can be found in Archive folder): -------------------------------------------------------
# Version 2.0 (March 2021): Initial public data release
# Update from version 3.0 to 3.1 (Oct 2025): Dr. Xuan Liu's updates dust attenuation correction for Al and Si
# Version 3.1.1 (Dec 2025): Bug Fix to avoid adding OC abd EC_FTR_MDL values again as a separate row
#
# Created by: --------------------------------------------------------------
# Crystal Weagle, 22 January 2019
# Revised by Haihui Zhu, April 2025, Nidhi Anchan, Oct 2025
##############################################################################
import os, sys
import logging, warnings
from datetime import datetime, date
import matplotlib as mpl
import pandas as pd

# Silence Matplotlib/fontTools warnings
if hasattr(mpl, "set_loglevel"):
    mpl.set_loglevel("warning")
logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("fontTools").setLevel(logging.WARNING)
logging.getLogger("fontTools.subset").setLevel(logging.WARNING)

# ---------------------------------------------------------------------------
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import spt_utils as su
import f_utils as fu

# ===============================
# Initial Setup and Configuration
# ===============================
debug_mode = 0
direc = su.find_root_dir(debug_mode)
DataVersion = 'Data version 3.1.1'  # will be in the first line of the public files

# ============================
# Directory & File Paths
# ============================
direc_master = os.path.join(direc, 'Analysis_Data', 'Master_files')
direc_output = os.path.join(direc, 'Public_Data')
log_file_path = os.path.join(direc, 'Public_Data/Data_Processing_Records/Master_File',
                                datetime.now().strftime('%Y-%m-%d-%H%MSS') + '_MasterFile_Processing_Record.log')
Sampling_Parameters_Methods = os.path.join(direc,'SOPs/Public SOPs/Sampling_Parameters_Methods_2.4.xlsx')

# # TEST 
# direc_master = os.path.join('release_test', 'Master_files')
# direc_output = os.path.join('release_test', 'Public_Data')
# log_file_path = os.path.join('release_test',  datetime.now().strftime('%Y-%m-%d-%H%MSS') + '_MasterFile_Processing_Record.txt')
# Sampling_Parameters_Methods = os.path.join('release_test/Sampling_Parameters_Methods_2.4.xlsx')

# Log file setup
su.setup_logging(log_file_path)  # Initialize logging
logging.info("F started")
logging.info(f'{datetime.now()} \n')

# Read the site info table
site_details = pd.read_excel(f'{direc}/Site_Sampling/Site_details.xlsx')

# define IC elog and MDL
icelog = f'{direc}/Analysis_Data/E-Logs/IC E-Log.xlsx' 

ic_mdl_tef_path = f'{direc}/Analysis_Data/Ion_Chromatography/IC_Teflon detection limit.xlsx'
ic_mdl_tef = pd.read_excel(ic_mdl_tef_path, header=4, index_col=0) # unit = ug/filter

ic_mdl_nyl_path = f'{direc}/Analysis_Data/Ion_Chromatography/IC_Nylon detection limit.xlsx'
ic_mdl_nyl = pd.read_excel(ic_mdl_nyl_path, header=0, skiprows=range(1, 6), nrows=1, index_col=0) 

# define XRF elof and MDL
xrfelog = f'{direc}/Analysis_Data/E-Logs/XRF E-Log.xlsx' 
xrf_mdl_path = f'{direc}/Analysis_Data/XRF/SPARTAN_XRF_MDL_Unc_wide.xlsx'
xrf_mdl = pd.read_excel(xrf_mdl_path)

# Load Sampling Methods 
# Sampling_Parameters_Methods = os.path.join(direc, 'SOPs', 'Public SOPs', 'Sampling_Parameters_Methods_2.4.xlsx')
water_content = pd.read_excel(Sampling_Parameters_Methods, sheet_name='Water_content', header=0, usecols="A:D",nrows=1)
RCFM_parameters = pd.read_excel(Sampling_Parameters_Methods, sheet_name='PM2.5 RCFM', header=0)
PM25_mass_para = pd.read_excel(Sampling_Parameters_Methods, sheet_name='PM2.5 mass', header=0)
PM25_IC_para = pd.read_excel(Sampling_Parameters_Methods, sheet_name='PM2.5 water-soluble ions', header=0)
PM25_carbon_para = pd.read_excel(Sampling_Parameters_Methods, sheet_name='PM2.5 BC_OC_EC', header=0)
Flag_parameters = pd.read_excel(Sampling_Parameters_Methods, sheet_name='Flags', header=0)

# Load the list of bad sodium cartridges 
BadNa = os.path.join(direc, 'Analysis_Data', 'Ion_Chromatography', 'Bad_sodium_data_in_2017.xlsx')
BadNa = pd.read_excel(BadNa, header=0)
BadNa = BadNa['CartID'].values


# Additional output initialization 
# BC source summary
bc_table_path = os.path.join(direc, 'Analysis_Data/Black_Carbon/BC_availability.xlsx')
BC_table = pd.DataFrame(
    index=site_details['Site_Code'],
    columns=[
        'HIPS_Measured_BC', 
        'HIPS_Est_BC', 
        'HIPS_total_BC', 
        'SSR_total_BC'
    ],
    dtype=float
)

# Ref: Mass type code
# 0 = cartridge (traveling) blank
# 1 = PM2.5
# 2 = PM10
# 3 = PMcoarse
# 4 = saturated nuclepore, filter measurements invalid [Not found at any site. The data might be too old to be used.]
# 5 = negative mass, filter measurements invalid
# 6 = invalid flow, filter measurements invalid

# ============================
# Start processing 
# ============================
# %% Read master files and process public files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# #Uncomment this if start from a specific site: 
# site_list = site_details['Site_Code'].tolist()
# start_site = 'ETAD' 
# start_idx = site_list.index(start_site)
# for idx in range(start_idx, len(site_list)):
#     site = site_details['Site_Code'][idx]
    
#Uncomment this if start from begining 
for idx, site in enumerate(site_details['Site_Code']): 
    
    # Read Master Files & IC_MDL files
    master_file = os.path.join(direc_master, f'{site}_master.csv')
    master_data = su.read_master(master_file)
    Titles, dtype = su.master_heading()
    
    if master_data.empty:
        continue
    
    # extrate site info 
    tsite_info = site_details.iloc[idx]
    
    # Read IC elog data
    ic_elog = fu.get_ic_elog_dates(icelog, site)
    
    # Read XRF elog data
    city = tsite_info['City']
    xrf_elog = fu.get_xrf_elog_dates(xrfelog, site, city)
    

    # ============================
    # Start data validation
    # ============================
    # STEP 1: For ICP-MS measurements, subtract field blank values in PM2.5 and PM10 filters
    # for more information, go to the function 'subtract_field_blank' in f_utils.py
    master_data = fu.icp_subtract_field_blank(master_data)
    
    # STEP 2: Exclude bad filters or filter with invalid sampling hour, method, or volume
    # for more information, go to the function 'remove_invalid_filters' in f_utils.py
    master_data = fu.remove_invalid_filters(master_data)
    
    # STEP 3: Estimate HIPS-BC when HIPS-BC not available
    # for more information, go to the function 'process_BC_data' in f_utils.py
    master_data, BC_table = fu.process_BC_data(master_data, site, BC_table)
    
    # STEP 4: Mask out invalid IC data: bad correlation coefficient 
    # for more information, go to the function 'mask_invalid_ic_data' in f_utils.py
    master_data = fu.mask_invalid_ic_data(master_data)
    
    # STEP 5: Mask out invalid IC data: bad Na conc in blanks
    # find cartridges with bad Na (Na too high in blanks)
    Ind = np.where(np.isin(BadNa, site))[0]
    if len(Ind) > 0:
        for ii in Ind:
            tCartridge = BadNa[ii]
            MasterInd = np.where(master_data['CartridgeID'] == tCartridge)[0]
            master_data.iloc[MasterInd, 'IC_Na_ug_T'] = np.nan  

    # STEP 6: Convert mass to concentration (ug/m3)
    master_data = fu.convert_mass_to_concentration(master_data)

    # STEP 7: Compare filter total mass to sum of measured species
    # remove filter if sum of species greater than collected mass by more than 10%
    # Mass = Metals + IC (without F, Cl, K, Mg, Ca) + BC
    master_data = fu.check_spec_sum(master_data)


    # ===================================
    # Prepare for generating public files
    # ===================================
    # STEP 1: Sort data by date
    master_data = master_data.sort_values(
        by=['start_year', 'start_month', 'start_day', 'start_hour'],
        kind='mergesort',  # Preserves original order for tied rows 
        ascending=[True, True, True, True],  # Chronological order
        na_position='last'  # Place rows with missing dates at the end
    ).reset_index(drop=True)  # Reset index for clean ordering
        
    # STEP 2: Sort into PM2.5, PM10, and PMc
    pm25_data = master_data[master_data['Mass_type'] == 1].copy()
    pm10_data = master_data[master_data['Mass_type'] == 2].copy()
    # pmc_data = master_data[master_data['Mass_type'] == 3].copy() # no longer update pm coarse
    
    if len(pm25_data) == 0:  # No PM2.5 data then no need to go further
        logging.info(f'No valid PM25 data for {site}')
        continue

    # --- Handle pending-weight filters gracefully (no deletions) ---
    ALLOW_PENDING_MASS = True   # Set False to enforce hard-stop/exit # If check is needed # Sanity check: PM2.5 mass should not be nan
                                                                                                  # This should not happen if step 2 in quality checks is done correctly
    VOLUME_COL = "Volume_m3"     # <-- change if your volume column differs

    # Step 1: rows where PM2.5 mass is NaN
    is_mass_nan = pm25_data['mass_ug'].isna()

    if is_mass_nan.any():
        # Step 2: among those, check if volume exists (>0)
        if VOLUME_COL in pm25_data.columns:
            has_volume = pm25_data[VOLUME_COL].notna() & (pm25_data[VOLUME_COL] > 0)
        else:
            logging.warning(
                f"Volume column '{VOLUME_COL}' not found; cannot verify volume for pending-weight detection."
            )
            has_volume = pd.Series(False, index=pm25_data.index)

        # Pending := PM2.5 NaN AND volume present (>0)
        pending = is_mass_nan & has_volume

        if pending.any():
            flagged_ids = pm25_data.loc[pending, 'FilterID'].astype(str).tolist()

            if ALLOW_PENDING_MASS:
                # Just warn; DO NOT delete or filter rows
                logging.warning(
                    f"{len(flagged_ids)} PM2.5 filters flagged: Volume exists but PM2.5 doesn't "
                    f"(likely in shipment). FilterIDs: {', '.join(flagged_ids)}"
                )
            else:
                # Strict mode: fail fast and say why
                logging.error(
                    "Detected filters where volume exists but PM2.5 is NaN (likely in shipment). "
                    f"ALLOW_PENDING_MASS=False, so exiting. FilterIDs: {', '.join(flagged_ids)}"
                )
                exit(1)
 

    # STEP 3: Apply dust attenuation
    # Sub-step 1: calculate dust loading for PM2.5 and PM10
    # Sub-step 2: correct Al and Si concentrations based on dust loading
    
    # ---- Sub-step 1: dust loading -------------
    # if XRF data available, define soil (aka mineral dust) by combination of {Al, Si, Ca, Ti, Fe} following Xuan Liu's eqn:
    # [dust] = [1.89Al×(1+MAL)+2.14Si+1.40Ca+1.36Fe+1.67Ti]×CF
    # if no XRF, only ICP-MS data available, soil is defined as:
    # [dust] = 10*([Al] + [Fe] + [Mg]), if no Fe use: [dust] = 10*([Mg] + 30*[Ti] + [Al])
    dust25_ugm3, dust_pm25filters = fu.get_dust(pm25_data, tsite_info) #soil,soilxvol
    dust10_ugm3, dust_pm10filters = fu.get_dust(pm10_data, tsite_info)
    
    # ---- Sub-step 2: Al and Si attenuation ----
    # Xuan Liu, Jan, 2025: 
    # ------------------ PM2.5 -------------------
    # A = 0.78 - 8.6e-4 * dust_loading + 4.0e-7 * dust_loading^2
    # [Si_new] = [Si] / A
    # [Al_new] = [Al] * 0.77 / A
    # The additional coefficient of 0.77 is included here to exclude the impact of attenuation included in the Al calibration curve.
    # ------------------ PM10 --------------------
    # A = (0.78 - 8.6e-4 * dust_loading + 4.0e-7 * dust_loading^2) * 0.91
    # The ratio of 0.91 is used to calculate the PM10 attenuation based on the PM2.5 attenuation.
    # [Si_new] = [Si] / A
    # [Al_new] = [Al] * 0.77 / A
    # 
    # The dust_loading used in the equations needs to be in the form of μg/cm2. The deposition area used to calculate dust mass loading is 3.53 cm2.
    # v3.0 data Code for the above calculations:
    # # PM2.5
    # dust_loading = dust_pm25filters / 3.53  # unit ug/cm2
    # A = 0.78 - 8.6e-4 * dust_loading + 4.0e-7 * dust_loading**2
    # pm25_data.loc[:, 'Si_XRF_ng'] = pm25_data.loc[:, 'Si_XRF_ng'].values.reshape(-1, 1) / A   
    # pm25_data.loc[:, 'Al_XRF_ng'] = pm25_data.loc[:, 'Al_XRF_ng'].values.reshape(-1, 1) * 0.77 / A    
    
    # # PM10  # The ratio of 0.91 is used to calculate the PM10 attenuation based on the PM2.5 attenuation.
    # dust_loading = dust_pm10filters / 3.53  # unit ug/cm2
    # A = 0.78 - 8.6e-4 * dust_loading + 4.0e-7 * dust_loading**2 * 0.91
    # pm10_data.loc[:, 'Si_XRF_ng'] = pm10_data.loc[:, 'Si_XRF_ng'].values.reshape(-1, 1) / A   
    # pm10_data.loc[:, 'Al_XRF_ng'] = pm10_data.loc[:, 'Al_XRF_ng'].values.reshape(-1, 1) * 0.77 / A    

    # Update Xuan Liu: Oct 14, 2025: v3.1 data
    # ------------------ PM2.5 -------------------
    # A = (0.78 - 8.6e-4 * dust_loading + 4.0e-7 * dust_loading^2) * f
    # [Si_new] = [Si] / A
    # [Al_new] = [Al] * 0.77 / A
    # The additional coefficient of 0.77 is used to exclude the impact of attenuation introduced by non-representative standards in Al calibration.
    
    # ------------------ PM10 --------------------
    # A = (0.78 - 8.6e-4 * dust_loading + 4.0e-7 * dust_loading^2) * 0.92 * f
    # [Si_new] = [Si] / A
    # The ratio of 0.92 is used to calculate the PM10 attenuation based on the PM2.5 attenuation.
    # [Al_new] = [Al] * 0.77 / A
    
    # The scaling factor is calculated as:
    # f = w_dust + (1 - w_dust) * 1.2
    # w_dust is the dust mass fraction. The dust_loading used in the equations needs to be in the form of μg/cm2. The deposition area used to calculate dust mass loading is 3.53 cm2.
    
    # Calculate dust mass fraction
    pm25_total = np.asarray(pm25_data['mass_ug'].to_numpy(), dtype=float).ravel() # ug/m3
    pm25_dust  = np.asarray(dust25_ugm3, dtype=float).ravel() #ug/m3
    w_dust25 = np.zeros_like(pm25_total, dtype=float)
    np.divide(pm25_dust, pm25_total, out=w_dust25, where=pm25_total > 0)
    #w_dust25 = np.clip(np.nan_to_num(w_dust25), 0.0, 1.0)

    pm10_total = np.asarray(pm10_data['mass_ug'].to_numpy(), dtype=float).ravel()
    pm10_dust  = np.asarray(dust10_ugm3, dtype=float).ravel()
    w_dust10 = np.zeros_like(pm10_total, dtype=float)
    np.divide(pm10_dust, pm10_total, out=w_dust10, where=pm10_total > 0)
    #w_dust10 = np.clip(np.nan_to_num(w_dust10), 0.0, 1.0)
    
    # Scaling Factor f
    f25 = w_dust25 + (1 - w_dust25) * 1.2
    f10 = w_dust10 + (1 - w_dust10) * 1.2
 
    # PM2.5
    dust_loading = np.asarray(dust_pm25filters, dtype=float).ravel() / 3.53  # ug/cm2 (force 1-D)
    A = (0.78 - 8.6e-4 * dust_loading + 4.0e-7 * dust_loading**2) * f25        
    A = np.asarray(A, dtype=float).ravel()                                   
    pm25_data.loc[:, 'Si_XRF_ng'] = pm25_data['Si_XRF_ng'].to_numpy(dtype=float) / A
    pm25_data.loc[:, 'Al_XRF_ng'] = pm25_data['Al_XRF_ng'].to_numpy(dtype=float) * 0.77 / A


    # PM10  
    dust_loading = np.asarray(dust_pm10filters, dtype=float).ravel() / 3.53   # ug/cm2 (force 1-D)
    A = (0.78 - 8.6e-4 * dust_loading + 4.0e-7 * dust_loading**2) * 0.92 * f10  
    A = np.asarray(A, dtype=float).ravel()                                     

    pm10_data.loc[:, 'Si_XRF_ng'] = pm10_data['Si_XRF_ng'].to_numpy(dtype=float) / A
    pm10_data.loc[:, 'Al_XRF_ng'] = pm10_data['Al_XRF_ng'].to_numpy(dtype=float) * 0.77 / A

    # ===================================
    # OUTPUT FILE 1: {Site}_speciation.csv
    # =================================== 
    fu.populate_SpecChem(pm25_data,Sampling_Parameters_Methods,tsite_info, direc_output, ic_mdl_tef, ic_mdl_nyl, ic_elog, xrf_mdl, xrf_elog, DataVersion, masstype='pm25')
    
    fu.populate_SpecChem(pm10_data,Sampling_Parameters_Methods,tsite_info, direc_output, ic_mdl_tef, ic_mdl_nyl, ic_elog, xrf_mdl, xrf_elog, DataVersion, masstype='pm10')
    
    # no longer update pm coarse, this line won't work:
    # pmc_chemspec = fu.populate_SpecChem(pmc_data,Sampling_Parameters_Methods,tsite_info, masstype='pmc') 
      

    # ===================================
    # Output FILE 2: Reconstruct PM2.5
    # =================================== 
    # List of components that will be collected in the RCFM_condense (not all are reported to the public)
    comp_list = ['FilterID','Year','Start_Date', 'End_Date', 'Filter PM2.5 Mass','Sulfate','Nitrate','Ammonium','Sea Salt','Fine Soil',
                 'Equivalent BC PM2.5','Trace Element Oxides','Organic Carbon','Organic Matter', 'Residual Matter',
                 'vol_growth','Particle Bound Water','Potassium', 'Chlorine','Sodium']
    # List of components in RCFM_condense that we don't report
    DoNotReport_List = ['Chlorine', 'Potassium', 'Sodium', 'vol_growth', 'Particle Bound Water', 'FilterID', 'Year', 'Start_Date', 'End_Date']
    
    # Note: all species have been converted to ug/m3 or ng/m3 in Quality Checks STEP 6.
    
    # Initialize a rcfm DataFrame for figures
    rcfm_df = pd.DataFrame(columns=comp_list)

    rcfm_titles = ['Site_Code', 'Latitude', 'Longitude', 'Elevation_meters', 'Filter_ID',
            'Start_Year_local', 'Start_Month_local', 'Start_Day_local', 'Start_hour_local',
            'End_Year_local', 'End_Month_local', 'End_Day_local', 'End_hour_local',
            'Hours_sampled', 'Parameter_Code', 'Parameter_Name', 'Value', 'Units',
            'Method_Code', 'Collection_Description', 'Analysis_Description', 'Conditions', 'Flag']
    # Initialize another rcfm DataFrame for output (with sampling parameters)
    rcfm_out = pd.DataFrame(columns=rcfm_titles)

    for idx, row in pm25_data.iterrows():
        
        # read info that might needed for method parameter identification
        method_idx = int(row['method_index'])-1 # iloc start from 0
        if row['projectID'] == 'S':
            sam_mode = 'SPARTAN'
        elif row['projectID'] == 'M':
            sam_mode = 'MAIA'
        
        nylon_exist = pd.notna(row.get('IC_NO3_ug_N'))
        xrf_exist = pd.notna(row.get('Al_XRF_ng'))
        icp_exist = pd.notna(row.get('Al_ICP_ng'))
        
        rcfm = {} # initiate a dictionary for *dry* components
        rcfm['FilterID'] = row['FilterID']
        rcfm['Year'] = row['start_year']
        hour = min(23, int(round(row['start_hour'])))
        rcfm['Start_Date'] = datetime(int(row['start_year']), int(row['start_month']), int(row['start_day']), hour)
        hour = min(23, int(round(row['stop_hour'])))
        rcfm['End_Date'] = datetime(int(row['stop_year']), int(row['stop_month']), int(row['stop_day']), hour)
        
        for component in comp_list:
            
            if component == 'Filter PM2.5 Mass':
                if nylon_exist:
                    rcfm[component] = row['mass_ug'] + row['IC_NO3_ug_N']* water_content['Inorganics'].values[0]
                    
                    method_params = PM25_mass_para[
                        (PM25_mass_para['Collection Description'].str.contains('Nylon')) &
                        (PM25_mass_para['Sampling Mode'] == sam_mode )
                    ].iloc[0] # Turning dataframe to series 
                    
                else:
                    rcfm[component] = row['mass_ug']
                    method_params = PM25_mass_para.iloc[method_idx]
                    
            elif component == 'Sulfate':
                rcfm[component] = row['IC_SO4_ug_T']
                 
                method_params = PM25_IC_para[
                    (PM25_IC_para['Parameter'].str.contains('Sulfate', case=False)) & 
                    (PM25_IC_para['Collection Description'].str.contains('Teflon')) &
                    (PM25_IC_para['Sampling Mode'] == sam_mode) 
                ]
                if not method_params.empty and sam_mode == 'SPARTAN':
                    method_params = method_params.iloc[method_idx] 
                    
                elif not method_params.empty and sam_mode == 'MAIA':
                    method_params = method_params.iloc[0]
                else:
                    id = row['FilterID']
                    logging.error('cannot find method for {id}, exiting')
                    exit()
            
            elif component == 'Nitrate':
                if nylon_exist:
                    rcfm[component] = row['IC_NO3_ug_T'] + row['IC_NO3_ug_N']
                    
                    method_params = PM25_IC_para[
                        (PM25_IC_para['Parameter'].str.contains('Nitrate', case=False)) & 
                        (PM25_IC_para['Collection Description'].str.contains('Nylon')) &
                        (PM25_IC_para['Sampling Mode'] == sam_mode) 
                    ].iloc[0]
                    
                    method_params['Collection Description'] += '; Added to a Teflon filter mass'
                    
                else:
                    rcfm[component] = row['IC_NO3_ug_T']   
                 
                    method_params = PM25_IC_para[
                        (PM25_IC_para['Parameter'].str.contains('Nitrate', case=False)) & 
                        (PM25_IC_para['Collection Description'].str.contains('Teflon')) &
                        (PM25_IC_para['Sampling Mode'] == sam_mode) 
                    ]
                    if not method_params.empty and sam_mode == 'SPARTAN':
                        method_params = method_params.iloc[method_idx] 
                        
                    elif not method_params.empty and sam_mode == 'MAIA':
                        method_params = method_params.iloc[0]
                    else:
                        id = row['FilterID']
                        logging.error('cannot find method for {id}, exiting')
                        exit()
        
            elif component == 'Ammonium':
                
                if nylon_exist:
                    rcfm[component] = row['IC_NH4_ug_T'] + row['IC_NH4_ug_N']*0.29
                    
                    method_params = RCFM_parameters[
                        (RCFM_parameters['Parameter']=='Ammonium') & 
                        (RCFM_parameters['Collection Description'].str.contains('Telfon filter plus')) 
                    ].iloc[0]
                    
                else:   
                    rcfm[component] = row['IC_NH4_ug_T'] 
                    
                    method_params = PM25_IC_para[
                        (PM25_IC_para['Parameter'].str.contains('Ammonium', case=False)) & 
                        (PM25_IC_para['Collection Description'].str.contains('Teflon')) &
                        (PM25_IC_para['Sampling Mode'] == sam_mode) 
                    ]
                    if not method_params.empty and sam_mode == 'SPARTAN':
                        method_params = method_params.iloc[method_idx] 
                        
                    elif not method_params.empty and sam_mode == 'MAIA':
                        method_params = method_params.iloc[0]
                    else:
                        id = row['FilterID']
                        logging.error('cannot find method for {id}, exiting')
                        exit()
                
            elif component == 'Sea Salt':
                # Sea salt(dry) = (2.54*[Na+] - 0.1[Al]) or (1.65[Cl] - 0.1[Al])
                
                # default method:
                method_params = RCFM_parameters[
                        (RCFM_parameters['Parameter']=='Sea Salt') & 
                        (RCFM_parameters['Analysis Description'].str.contains('XRF')) 
                    ].iloc[0]
                
                if xrf_exist:
                    rcfm[component] = 2.54 * row['IC_Na_ug_T'] - 0.1 * row['Al_XRF_ng']/1000
                    
                elif icp_exist:
                    rcfm[component] = 2.54 * row['IC_Na_ug_T'] - 0.1 * row['Al_ICP_ng']/1000
                    # if use ICP, update method_params
                    method_params = RCFM_parameters[
                            (RCFM_parameters['Parameter']=='Sea Salt') & 
                            (RCFM_parameters['Analysis Description'].str.contains('ICP-MS')) 
                        ].iloc[0]
                else:
                    rcfm[component] = 2.54 * row['IC_Na_ug_T']
                
                if rcfm[component] < -0.3: # nan out if too negative
                    rcfm[component] = np.nan
                
            elif component == 'Fine Soil':
                soil, _ = fu.get_dust(row.to_frame().T, tsite_info)
                rcfm[component] = soil[0][0]
                
                if xrf_exist:
                    method_params = RCFM_parameters[
                            (RCFM_parameters['Parameter']=='Fine Soil') & 
                            (RCFM_parameters['Analysis Description'].str.contains('XRF')) 
                        ].iloc[0]
                else:
                    method_params = RCFM_parameters[
                            (RCFM_parameters['Parameter']=='Fine Soil') & 
                            (RCFM_parameters['Analysis Description'].str.contains('ICP-MS')) 
                        ].iloc[0]
                
            elif component == 'Trace Element Oxides':
                rcfm[component] = fu.get_teo(row.to_frame().T, xrf_exist, icp_exist)
                
                if xrf_exist:
                    method_params = RCFM_parameters[
                            (RCFM_parameters['Parameter']=='Trace Element Oxides') & 
                            (RCFM_parameters['Analysis Description'].str.contains('XRF')) 
                        ].iloc[0]
                else:
                    method_params = RCFM_parameters[
                            (RCFM_parameters['Parameter']=='Trace Element Oxides') & 
                            (RCFM_parameters['Analysis Description'].str.contains('ICP-MS')) 
                        ].iloc[0]
                
            elif component == 'Equivalent BC PM2.5':
                rcfm[component] = row['BC_HIPS_ug']

                carbon_para = PM25_carbon_para
                    
                if 'HIPS-BC estimated using SSR-BC' in row['Flags']:
                    method_params = carbon_para[
                            (carbon_para['Analysis Description'].str.contains('Estimated')) &
                            (carbon_para['Sampling Mode'] == sam_mode) 
                    ].iloc[0]
                    
                else:
                    method_params = carbon_para[
                            (carbon_para['Analysis Description'] == 'HIPS') &
                            (carbon_para['Sampling Mode'] == sam_mode) 
                        ].iloc[0]
                    
            elif component == 'Organic Carbon':
                rcfm[component] = row['OC_FTIR_ug']
                
                carbon_para = PM25_carbon_para
                
                method_params = carbon_para[
                        (carbon_para['Parameter'].str.contains('OC', case=False)) & 
                        (carbon_para['Analysis Description'] == 'FTIR') &
                        (carbon_para['Sampling Mode'] == sam_mode) 
                    ].iloc[0]
                
            elif component == 'Organic Matter':
                rcfm[component] = row['OM_FTIR_ug']
                
                carbon_para = PM25_carbon_para
                
                method_params = carbon_para[
                        (carbon_para['Parameter'].str.contains('OM', case=False)) & 
                        (carbon_para['Analysis Description'] == 'FTIR') &
                        (carbon_para['Sampling Mode'] == sam_mode) 
                    ].iloc[0]
                    
            elif component == 'Residual Matter':
                # 1. only proceed if major components available
                if pd.notna(rcfm['Sulfate']) & pd.notna(rcfm['Nitrate']) & pd.notna(rcfm['Ammonium']) &\
                    pd.notna(rcfm['Fine Soil']): 
                    
                    # 2. calc wet mass    
                    inorgic_wet = (rcfm['Sulfate'] + rcfm['Nitrate'] + rcfm['Ammonium']+   np.nan_to_num(row['IC_K_ug_T'])) * (1 + water_content['Inorganics'].values)
                    
                    organic_wet =  np.nan_to_num(rcfm['Organic Carbon']) * (1 + water_content['RM'].values)
                    
                    ss_wet =  np.nan_to_num(rcfm['Sea Salt']) * (1 + water_content['SS'].values)
                    
                    # 3. calc residual
                    rm = rcfm['Filter PM2.5 Mass'] - inorgic_wet - organic_wet - ss_wet - np.nan_to_num(rcfm['Equivalent BC PM2.5']) -  np.nan_to_num(rcfm['Trace Element Oxides']) - rcfm['Fine Soil']
                    
                    rcfm[component] = rm[0]
                    
                    method_params = RCFM_parameters[ (RCFM_parameters['Parameter']=='Residual Matter')  ].iloc[0]
                    
                    # prepare water for vol_growth
                    inorgic_water = (rcfm['Sulfate'] + rcfm['Nitrate'] + rcfm['Ammonium']+   np.nan_to_num(row['IC_K_ug_T'])) * water_content['Inorganics'].values
                    
                    organic_water =  np.nan_to_num(rcfm['Organic Carbon']) * water_content['RM'].values
                    
                    ss_wet_water =  np.nan_to_num(rcfm['Sea Salt']) * water_content['SS'].values
                    
                    water_mass = inorgic_water + organic_water + ss_wet_water
                    water_mass = water_mass[0]
                    
                else:
                    rcfm[component] = np.nan
                    water_mass = np.nan
                    
            elif component == 'Particle Bound Water':
                # similar to wet mass but without adding '1' and without worrying about nan
                # not reported to the public files
                inorgic_water = water_content['Inorganics'].values * ( np.nan_to_num(rcfm['Sulfate']) + np.nan_to_num(rcfm['Nitrate']) +\
                    np.nan_to_num(rcfm['Ammonium']) +  np.nan_to_num(row['IC_K_ug_T'])) 
                
                organic_water = water_content['RM'].values * np.nan_to_num(rcfm['Organic Carbon']) 
                
                ss_water = water_content['SS'].values * np.nan_to_num(rcfm['Sea Salt']) 
                
                # Particle bounded water is the sum of water from above
                pbw = inorgic_water  + organic_water + ss_water
                rcfm[component] = pbw[0]
                
            elif component == 'vol_growth':
                # Species     = [SNA   RM    OC  NaCl Soil BC   TEO  ]
                # Density_Dry = [1.70, 1.3, 1.3, 2.2, 2.6, 1.8, 2.5  ]  # same as in GEOS-Chem (except for TEO)
                # vol_wet = vol_dry + vol_water
                # vol_growth = vol_wet/vol_dry = 1 + vol_water/vol_dry (volumn growth factor)
                vol_dry = (rcfm['Sulfate']/1.7 + rcfm['Nitrate']/1.7 + rcfm['Ammonium']/1.7 + rcfm['Fine Soil']/2.6 +  \
                            np.nan_to_num(rcfm['Organic Carbon'])/1.3+ np.nan_to_num(rcfm['Sea Salt'])/2.2 +\
                            np.nan_to_num(rcfm['Equivalent BC PM2.5'])/1.8 + np.nan_to_num(rcfm['Trace Element Oxides'])/2.5)
                 
                # if any of SNA and soil is nan, vol_water will be nan. So it is fine to have vol_dry being nan too. 
                vol_water = water_mass/1.0
                vol_growth = 1+ vol_water/vol_dry 

                rcfm[component] = vol_growth

            # The following two ions are not reported in RCFM, just for internal figure purpose
            elif component == 'Potassium':
                rcfm[component] = row['IC_K_ug_T']
                
            elif component == 'Chlorine':
                rcfm[component] = row['IC_Cl_ug_T']
                
            elif component == 'Sodium':
                rcfm[component] = row['IC_Na_ug_T']
                
                
            # Now adding method information and copy to rcfm_df
            if pd.notna(rcfm[component]) and component not in DoNotReport_List:
                new_row = {}
                
                if component == 'Trace Element Oxides':
                    new_row['Value'] = round(rcfm[component] * 1000, 2) # reported as ng/m3
                else: 
                    # adding value of the measurement
                    new_row['Value'] = round(rcfm[component], 2)
                
                # basic site information
                new_row['Site_Code'] = tsite_info['Site_Code']
                new_row['Latitude'] = tsite_info['Latitude']
                new_row['Longitude'] = tsite_info['Longitude']
                new_row['Elevation_meters'] = tsite_info['Elevation_meters']
                
                # filter and time information from pm25_data
                new_row['Filter_ID'] = row['FilterID']  # Assuming column is named FilterID
                new_row['Start_Year_local'] = row['start_year']
                new_row['Start_Month_local'] = row['start_month']
                new_row['Start_Day_local'] = row['start_day']
                new_row['Start_hour_local'] = row['start_hour']
                new_row['End_Year_local'] = row['stop_year']
                new_row['End_Month_local'] = row['stop_month']
                new_row['End_Day_local'] = row['stop_day']
                new_row['End_hour_local'] = row['stop_hour']
                new_row['Hours_sampled'] = row['hours_sampled']
                
                
                # Fill method-related information
                new_row['Parameter_Code'] = method_params['Parameter Code']
                new_row['Parameter_Name'] = method_params['Parameter']
                new_row['Units'] = method_params['Units']
                new_row['Collection_Description'] = method_params['Collection Description']
                new_row['Analysis_Description'] = method_params['Analysis Description']
                new_row['Method_Code'] = method_params['Method Code']
                new_row['Conditions'] = 'Ambient local'
                
                # adding flag
                new_row['Flag'] = fu.copyflag( row['Flags'], Flag_parameters)
                
                # Append the new row to RCFM output dataframe (will save to csv later)
                rcfm_out = pd.concat([rcfm_out, pd.DataFrame([new_row])], ignore_index=True)
        
        # append to the rcfm_df as the condensed version of RCFM (internal use only)
        rcfm_df = pd.concat([rcfm_df, pd.DataFrame([rcfm])], ignore_index=True)
    
    # save RCFM output CSV file
    fname = f'{direc_output}/RCFM/{tsite_info.Site_Code}_PM25_RCFM.csv'
    # generate metadata lines
    today = date.today().isoformat()
    metadata = [
        f'# File Updated: {today}',
        f'# {DataVersion}',
        f'# Site: {tsite_info.City}, {tsite_info.Country}'
    ]

    # Write metadata + DataFrame
    with open(fname, 'w') as f:
        for line in metadata:
            f.write(line + '\n')
        rcfm_out.to_csv(f, index=False)
        
    logging.info(f'{fname} saved') 
    
    # save a internal used RCFM summary file
    fname = f'{direc_output}/RCFM/RCFM_condense/{tsite_info.Site_Code}_PM25_RCFM_PieChartData.csv'
    rcfm_df.to_csv(fname, index=False, float_format='%.2f')
    logging.info(f'{fname} saved')      

    # ===============================================================
    # FIGURE 1: Time series and site average PM2.5 for website
    # ===============================================================
    savedir = f'{direc_output}/Chemical_Filter_Data/Plots/Filter_PM25_plots'
    fname = f"{tsite_info.Site_Code}_PM25_SiteAvg.png"
    fu.timeseries_pm25(rcfm_df,fname, savedir, tsite_info.City)
    
    # =================================================================
    # FIGURE 2: Bar chart - all data with dates or filter IDs as labels
    # =================================================================
    savedir = f'{direc_output}/Chemical_Filter_Data/Plots/'
    fname = f"{tsite_info.Site_Code}_PM25_Spec_all.png"
    fu.bar_all(rcfm_df,fname, savedir, tsite_info.City)
    
    # =================================================================
    # FIGURE 3: Bar chart -  most recent and complete data for website
    # =================================================================
    savedir = f'{direc_output}/Chemical_Filter_Data/Plots/'
    fname = f"{tsite_info.Site_Code}_PM25_Spec_website.png"
    fu.bar_website(rcfm_df,fname, savedir, tsite_info.City)
    
    # =================================================================
    # FIGURE 4: RCFM Pie chart - all data
    # =================================================================
    # There are several versions of pie charts:
    # no RM version
    savedir = f'{direc_output}/RCFM/Pie_spec_plots'
    fname = f"{tsite_info.Site_Code}_PM25_RCFMavg"
    fu.rcfm_pie(rcfm_df,fname, savedir, tsite_info.City)
    
    # with RM version
    savedir = f'{direc_output}/RCFM/Pie_spec_plots_withRM'
    fname = f"{tsite_info.Site_Code}_PM25_RCFMavg"
    fu.rcfm_wRM(rcfm_df,fname,savedir, tsite_info.City)
    
    # OM + RM, no OC
    # savedir = f'{direc_output}/RCFM/Pie_spec_plots_OM'
    # fname = f"{tsite_info.Site_Code}_PM25_RCFMavg"
    # fu.rcfm_OM_wRM(rcfm_df,fname,savedir,tsite_info.City)
    
    # =================================================================
    # FIGURE 5: RCFM Pie chart - 2020 and forward, with Cl and K
    # =================================================================
    # select filters by year
    rcfm_washu = rcfm_df[(rcfm_df['Year']>= 2020)]
     
    # no RM version
    savedir = f'{direc_output}/RCFM/Pie_spec_plots/MTL_filters'
    fname = f"{tsite_info.Site_Code}_PM25_RCFMavg"
    fu.rcfm_with_K_Cl(rcfm_washu,fname,savedir,tsite_info.City)
    
    # with RM version
    savedir = f'{direc_output}/RCFM/Pie_spec_plots_withRM/MTL_filters'
    fname = f"{tsite_info.Site_Code}_PM25_RCFMavg"
    fu.rcfm_with_K_Cl_wRM(rcfm_washu,fname, savedir,tsite_info.City)
    
    # OM + RM, no OC
    # savedir = f'{direc_output}/RCFM/Pie_spec_plots_OM/MTL_filters'
    # fname = f"{tsite_info.Site_Code}_PM25_RCFMavg"
    # fu.rcfm_OM_with_K_Cl_wRM(rcfm_df,fname,savedir, tsite_info.City)
    
    plt.close('all')
    logging.info(f'Finished processing data for {site} \n\n')




# =================================================================
#  FINALLY: Print out data summary tables  
# =================================================================
# BC summary table
BC_table.reset_index().rename(columns={'index': 'Site_Code'}).to_excel( bc_table_path, index=False )

logging.info(f'END of Master File processing for {datetime.now()}')
