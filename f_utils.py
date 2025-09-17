# supporting functions for SPARTAN scripts
# version 1.0.0
# created by Haihui Zhu,  Nidhi Anchan
# April 2025
import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from datetime import date
import logging
import spt_utils as su
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")
logging.getLogger('matplotlib').setLevel(logging.WARNING)  # Suppress Matplotlib Logs
from matplotlib import font_manager as fm
from pathlib import Path

# ---- Register user-level Arial fonts (no admin required) ----
AR_DIR = Path("/storage1/fs1/rvmartin/Active/SPARTAN-shared/Public_Data/Scripts/.fonts")

# Be flexible about file names (Arial.ttf vs arial.ttf etc.)
candidates = []
candidates += list(AR_DIR.glob("Arial*.ttf"))
candidates += list(AR_DIR.glob("arial*.ttf"))
candidates += list(AR_DIR.glob("Arial*.ttc"))
candidates += list(AR_DIR.glob("arial*.ttc"))

found = 0
for p in candidates:
    try:
        fm.fontManager.addfont(str(p))
        found += 1
    except Exception as e:
        logging.warning(f"Could not add font {p}: {e}")

if found == 0:
    logging.warning(f"No Arial fonts found under {AR_DIR}; will fall back to DejaVu Sans.")

# ---- Use Arial for plots; keep text editable in PDF/SVG ----
mpl.rcParams.update({
    "text.usetex"  : False,     # TeX would outline text
    "pdf.fonttype" : 42,        # embed TrueType (editable in AI)
    "svg.fonttype" : "none",    # keep <text> nodes in SVG
    "ps.useafm"           : False,
    "pdf.use14corefonts"  : False,
    "font.family"  : "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],  # Arial first, fallback ok
})

def get_ic_elog_dates(icelog, site):
    excel = pd.ExcelFile(icelog, engine='openpyxl')
    sheet_name = next((s for s in excel.sheet_names if site in s), None)
    if sheet_name is None:
        logging.error(f"No sheet found for site '{site}' in elog")
        return None
    
    elog = pd.read_excel(icelog, sheet_name=sheet_name)
    elog.columns = elog.columns.str.rstrip()
    elog['Filter ID'] = elog['Filter ID'].str.strip() # not sure why the above line didn't work for Filter ID
    # reformat filter ID:
    elog['Filter ID'] = su.filter_id_format(elog['Filter ID'])
    # Convert column to datetime, forcing errors to NaT
    elog['Extraction Date'] = pd.to_datetime(elog['Extraction Date'], errors='coerce')
    # Replace NaT (invalid entries) with the default date
    elog['Extraction Date'] = elog['Extraction Date'].fillna(pd.Timestamp('2014-01-01'))
    # Then convert to dict with Analysis ID as key
    elog_selected = elog.set_index('Filter ID')['Extraction Date'].to_dict()

    return elog_selected
    
    
def get_xrf_elog_dates(xrfelog, site_ID, city):
    sheet_name = f'{site_ID}_{city}'
    xls = pd.ExcelFile(xrfelog)
    if sheet_name in xls.sheet_names:
        elog = pd.read_excel(xrfelog, sheet_name=sheet_name)
        elog.columns = elog.columns.str.rstrip()
        elog_selected = elog.set_index('Cartridge ID')['Analysis Date'].to_dict()
    else:
        elog_selected = None # for sites without any XRF data. This won't be used. 
        
    return elog_selected
    
    
def icp_subtract_field_blank(master_data):
    
    # Create subgroups for cartridge blocks (split every 8 rows)
    master_data["block"] = master_data.groupby("CartridgeID").cumcount()
    master_data["subgroup"] = (
        (master_data["CartridgeID"] != master_data["CartridgeID"].shift()) | 
        (master_data["block"] % 8 == 0)
    ).cumsum()
    master_data = master_data.drop(columns=["block"])
    
    # Identify ICP columns (assuming they start with 'ICP_')
    icp_cols = [col for col in master_data.columns if col.startswith('ICP_')]
    
    # Process each subgroup separately
    for subgroup_id, subgroup in master_data.groupby('subgroup'):
        # if full of nan, skip
        if subgroup['Volume_m3'].isna().all():
            continue
        
        # Find field blank (Mass_type=0) in this subgroup
        blanks = subgroup[subgroup['Mass_type'] == 0]
        
        if blanks.empty:
            logging.info(f"No field blank found in for cartridge {subgroup['CartridgeID'].iloc[0]}")
            continue
            
        # Get first blank (assuming one per subgroup)
        blank_values = blanks.iloc[0][icp_cols]
        
        # Get sample indices to adjust (Mass_type 1/2)
        sample_indices = subgroup.index[subgroup['Mass_type'].isin([1, 2])]  
        
        # Subtract blank values from samples in MAIN DataFrame
        master_data.loc[sample_indices, icp_cols] -= blank_values.values
        
    return master_data.drop(columns=['subgroup'])

def remove_invalid_filters(master_data):
    
    # Create a copy to avoid modifying the original DataFrame
    filtered_master_data = master_data.copy()
    
    # Condition 1: Exclude rows where 'Flags' contains exclusion phrase
    flag_mask = filtered_master_data['Flags'].str.contains(
        'Need to be excluded', 
        na=False  # Treat NaN/None as False
    )
    
    # Condition 2: Exclude rows with bad volume or sampling hours
    volume_hours_mask = (
        filtered_master_data[['Volume_m3', 'hours_sampled']].isna().any(axis=1) |  # NaN check
        (filtered_master_data['Volume_m3'] == 0) |  # Zero volume check
        (filtered_master_data['hours_sampled'] == 0)  # Zero hours check
    )
    
    # Condition 3: Exclude rows with missing method index
    method_mask = filtered_master_data['method_index'].isna()
    
    # Combine all exclusion masks
    combined_mask = flag_mask | volume_hours_mask | method_mask
    
    # Apply inverse mask to keep valid rows
    filtered_master_data = filtered_master_data[~combined_mask]
    
    return filtered_master_data.reset_index(drop=True)

def process_BC_data(master_data, site, BC_table):

    # 1. Replace -899 with NaN in BC_SSR_ug
    master_data['BC_SSR_ug'] = master_data['BC_SSR_ug'].replace(-899, np.nan)
    
    # 2. count initial BC data
    valid_mass = master_data['Mass_type'].isin([1, 2]) 
    BC_table.loc[site]['HIPS_Measured_BC'] = master_data.loc[valid_mass & master_data['BC_HIPS_ug'].notna()].shape[0]
    
    # 3. Estimate missing HIPS values and add flags
    mask = (
        master_data['Mass_type'].isin([1, 2]) &
        master_data['BC_SSR_ug'].notna() & 
        (master_data['BC_SSR_ug'] > 0) & 
        master_data['BC_HIPS_ug'].isna()
    )
    # 3.1 count Est HIPS data
    BC_table.loc[site]['HIPS_Est_BC'] = master_data.loc[mask].shape[0]
    
    if mask.any():
        # Coefficients for estimation
        a1 = 2163.54158502065
        b1 = 0.00198619072508564
        c1 = -0.00467334844964624
        a2 = 262.938831079902
        b2 = 0.0123767218647736
        c2 = 3.10315712581732
        
        # Calculate HIPS estimates
        SSR = master_data.loc[mask, 'BC_SSR_ug']
        HIPS = a1 * np.sin(b1 * SSR + c1) + a2 * np.sin(b2 * SSR + c2)
        
        # Update values and flags
        master_data.loc[mask, 'BC_HIPS_ug'] = HIPS
        master_data['Flags'] = master_data['Flags'].fillna('')
        master_data.loc[mask, 'Flags'] += 'HIPS-BC estimated using SSR-BC; '
    
    # 4. Count valid filters
    valid_mass = master_data['Mass_type'].isin([1, 2])
    BC_table.loc[site]['HIPS_total_BC'] = master_data.loc[valid_mass & master_data['BC_HIPS_ug'].notna()].shape[0]
    BC_table.loc[site]['SSR_total_BC'] = master_data.loc[valid_mass & master_data['BC_SSR_ug'].notna()].shape[0]
    
    # 5. Check for negative SSR values
    mask = master_data['BC_SSR_ug'] < 0
    # if mask.any():
    #     logging.info(f"Negative BC_SSR_ug values found for {master_data.loc[mask, 'Filter_ID'].values}")
    
    return master_data, BC_table

def mask_invalid_ic_data(master_data):
    # Set up the ions list
    ions = {'F','Cl','NO2','Br','NO3','PO4','SO4','Li','Na','NH4','K','Mg','Ca'}

    for ion in ions:
        # Create column names
        t_col = f'IC_{ion}_ug_T'
        n_col = f'IC_{ion}_ug_N'
        
        # 1. Handle S-{ion} and Sx-{ion} flags for T column
        t_mask = master_data['Flags'].str.contains(
            fr'\b(?:S-{ion}|Sx-{ion})\b', 
            na=False, 
            regex=True
        )
        
        if t_col in master_data.columns and t_mask.any():
            master_data.loc[t_mask, t_col] = np.nan
        
        # 2. Handle SN-{ion} and SxN-{ion} flags for N column
        n_mask = master_data['Flags'].str.contains(
            fr'\b(?:SN-{ion}|SxN-{ion})\b', 
            na=False, 
            regex=True
        )
        if n_col in master_data.columns and n_mask.any():
            master_data.loc[n_mask, n_col] = np.nan
            
    return master_data
  
def convert_mass_to_concentration(df):
    # Get positions of the boundary columns
    vol_col = 'Volume_m3'
    method_col = 'method_index'
    vol_idx = df.columns.get_loc(vol_col) + 1  # Start after Volume_m3
    method_idx = df.columns.get_loc(method_col)  # End before method_index
  
    # Get columns between Volume_m3 and method_index
    cols_to_convert = df.columns[vol_idx:method_idx]
    
    # Convert mass loading to concentration
    df[cols_to_convert] = df[cols_to_convert].div(df[vol_col], axis=0)
    
    # don't forget the mass column:
    df['mass_ug'] = df['mass_ug'].div(df[vol_col], axis=0)
    
    return df
  
def check_spec_sum(master_data):
    """Compare sum of chemical components with total filter mass. """
    """Remove filters when sum of chemical components is higher than total filter mass. """
    # Major IC columns 
    initial_targets = [
        'IC_NO2_ug_T', 'IC_Br_ug_T', 'IC_NO3_ug_T', 'IC_PO4_ug_T', 
        'IC_SO4_ug_T', 'IC_Li_ug_T', 'IC_Na_ug_T', 'IC_NH4_ug_T'
    ]
    

    # Include metal but not Na and S (they can come from IC)
    xrf_cols = [col for col in master_data.columns if '_XRF_' in col and not any(excl in col for excl in ['Na_', 'S_'])]
    icp_cols = [col for col in master_data.columns if '_ICP_' in col and not any(excl in col for excl in ['Na_', 'S_'])]
    # choices of carbon
    ftir_cols = ['OC_FTIR_ug', 'EC_FTIR_ug']
    bc_col = ['BC_HIPS_ug']
    
    # convert units of metals from ng to ug
    master_data[xrf_cols] = master_data[xrf_cols] / 1000
    master_data[icp_cols] = master_data[icp_cols] / 1000

    for idx, row in master_data.iterrows():
    
        targets = set(initial_targets)
        
        # Add XRF/ICP columns
        if pd.notna(row.get('Al_XRF_ng')): 
            targets.update(xrf_cols) # adding xrf if data available
        else:
            targets.update(icp_cols) # only use icp when xrf not available
        
        # Add FTIR or HIPS columns
        if pd.notna(row.get('OC_FTIR_ug')):
            targets.update(ftir_cols) # adding ftir when data available
        else:
            targets.update(bc_col) # use HIPS (or SSR inferred HIPS) when there isn't FTIR
        
        # Filter valid columns
        valid_cols = [col for col in targets if col in master_data.columns]
        
        # Calculate sum with NaN handling
        sum_val = row[valid_cols].fillna(0).sum()
        
        # Mass comparison logic
        mass_ug = row.get('mass_ug')
        if pd.notna(mass_ug):
            if sum_val > mass_ug:
                logging.info(f"Filter {row.get('FilterID', '')}: sum ({sum_val:.2f}) > mass ({mass_ug:.2f})")
                
                if sum_val > 1.1 * mass_ug:
                    master_data = master_data.drop(index=idx)
                    logging.info(f"Filter mass for {row.get('FilterID', '')} is invalid: exceeding 10% allowance.")
                
    # convert units of metals back to ng
    master_data[xrf_cols] = master_data[xrf_cols] * 1000
    master_data[icp_cols] = master_data[icp_cols] * 1000
    
    return master_data


# ========= Output Files =========
def get_dust(pm_data, site_details):
    
    # if XRF data available, define soil (aka mineral dust) by combination of {Al, Si, Ca, Ti, Fe} following Xuan Liu's eqn:
    # [SOIL] = [1.89Al×(1+MAL)+2.14Si+1.40Ca+1.36Fe+1.67Ti]×CF
    # if no XRF, only ICP-MS data available, soil is defined as:
    # [SOIL] = 10*([Al] + [Fe] + [Mg]), if no Fe use: [SOIL] = 10*([Mg] + 30*[Ti] + [Al])
    
    soil_elements = ['Al_XRF_ng', 'Si_XRF_ng', 'Ca_XRF_ng', 'Fe_XRF_ng', 'Ti_XRF_ng']
    xrf_soil_factors = [1.89, 2.14, 1.40, 1.36, 1.67]
    
    xrf_soil_ele_idx = []
    filtered_factors = []
    
    MAL = site_details['MAL'] 
    CF = site_details['CF'] 
    if np.isnan(MAL):
        raise ValueError(f'MAL and CF not found for {site_details.Site_Code} in Site_Details.xlsx')

    Soil = np.full((pm_data.shape[0], 1), np.nan)  # unit ug/m3
    
    for idx, (df_idx, row) in enumerate(pm_data.iterrows()):
        if pd.notna(row.get('Al_XRF_ng')):  # there is XRF data
            # ensure order correct:
            for i, elem in enumerate(soil_elements):
                if elem == 'Al_XRF_ng':
                    tSoil = row.get(elem) * xrf_soil_factors[i] * (1 + MAL)
                else:
                    tSoil += row.get(elem) * xrf_soil_factors[i]
                    
            Soil[idx] = tSoil * CF / 1000  # convert to ug/m3
            
        elif pd.notna(row.get('Al_ICP_ng')):  # no XRF, but there is ICP-MS data
            if pd.notna(row.get('Fe_ICP_ng')):  # Fe is available in ICP-MS data
                Soil[idx] = (10 * (row.get('Al_ICP_ng') + row.get('Fe_ICP_ng')+ row.get('Mg_ICP_ng'))) / 1000  # convert to ug/m3  
            else:  # no Fe in ICP-MS data
                Soil[idx] = (10 * (row.get('Al_ICP_ng') + 30*row.get('Ti_ICP_ng')+ row.get('Mg_ICP_ng'))) / 1000  # convert to ug/m3  
        else:  # no ICP-MS or XRF data
            Soil[idx] = np.nan

    SoilxVol = Soil * pm_data.loc[:, 'Volume_m3'].values.reshape(-1, 1)  # unit ug
    
    return Soil, SoilxVol

def copyflag(fullflag, Flag_parameters):
    publicflags = Flag_parameters['Flag'].values
    flagout = ""
    # Check each public flag
    for flag in publicflags:
        if flag in fullflag:  # Check if this public flag exists in fullflag
            if flagout:  # If flagout already has content
                flagout += "; " + flag
            else:
                flagout = flag
    return flagout
    
def populate_SpecChem(pm_data, Sampling_Parameters_Methods, site_info, direc_output, ic_mdl_tef, ic_mdl_nyl, ic_elog, xrf_mdl, xrf_elog,DataVersion, masstype):
    logging.info(f'Writing to {masstype} SpecChem CSV file')
  
    ChemSpec_Titles = [
        'Site_Code', 'Latitude', 'Longitude', 'Elevation_meters', 'Filter_ID',
        'Start_Year_local', 'Start_Month_local', 'Start_Day_local', 'Start_hour_local',
        'End_Year_local', 'End_Month_local', 'End_Day_local', 'End_hour_local',
        'Hours_sampled', 'Parameter_Code', 'Parameter_Name', 'Value', 'Units',
        'Analytical_MDL', 'MDL', 'UNC', 'Method_Code', 'Collection_Description',
        'Analysis_Description', 'Conditions', 'Flag'
    ]
    # Initialize empty DataFrame with specified columns
    pm25_chemspec = pd.DataFrame(columns=ChemSpec_Titles)

    # components contained in PM25_data_public matrix for searching the parameters sheet
    element_codes = {'SO4': 'Sulfate', 'NO3': 'Nitrate', 'NH4': 'Ammonium', 'Na': 'Sodium', 'PO4': 'Phosphate',
        'NO2': 'Nitrite', 'Br': 'Bromide', 'K': 'Potassium', 'Mg': 'Magnesium', 'Ca': 'Calcium', 'Li': 'Lithium', 
        'Al': 'Aluminum', 'P': 'Phosphorus', 'Ti': 'Titanium', 'V': 'Vanadium', 'Cr': 'Chromium', 'Mn': 'Manganese', 
        'Fe': 'Iron', 'Co': 'Cobalt', 'Ni': 'Nickel', 'Cu': 'Copper', 'Zn': 'Zinc', 'As': 'Arsenic', 'Se': 'Selenium', 
        'Au': 'Gold', 'Cd': 'Cadmium', 'Sb': 'Antimony', 'Ba': 'Barium', 'Ce': 'Cerium', 'Pb': 'Lead', 'Si': 'Silicon', 
        'S': 'Sulfur', 'Cl': 'Chlorine', 'Rb': 'Rubidium', 'Sr': 'Strontium', 'Sn': 'Tin'}
    
    PM25_mass_para = pd.read_excel(Sampling_Parameters_Methods, sheet_name='PM2.5 mass', header=0)
    PM25_IC_para = pd.read_excel(Sampling_Parameters_Methods, sheet_name='PM2.5 water-soluble ions', header=0)
    PM25_metals_para = pd.read_excel(Sampling_Parameters_Methods, sheet_name='PM2.5 trace elements', header=0)
    PM25_carbon_para = pd.read_excel(Sampling_Parameters_Methods, sheet_name='PM2.5 BC_OC_EC', header=0)

    PM10_mass_para = pd.read_excel(Sampling_Parameters_Methods, sheet_name='PM10 mass', header=0)
    PM10_IC_para = pd.read_excel(Sampling_Parameters_Methods, sheet_name='PM10 water-soluble ions', header=0)
    PM10_metals_para = pd.read_excel(Sampling_Parameters_Methods, sheet_name='PM10 trace elements', header=0)
    PM10_carbon_para = pd.read_excel(Sampling_Parameters_Methods, sheet_name='PM10 BC_OC_EC', header=0)

    Flag_parameters = pd.read_excel(Sampling_Parameters_Methods, sheet_name='Flags', header=0)
    
    # Loop through each row in pm data
    for idx, row in pm_data.iterrows():
        
        for col in pm_data.columns:
            
            if pd.notna(row.get(col)): # only report data that are not a nan
                
                # Create a new row dictionary
                new_row = {}
                
                # read info that might needed for method parameter identification
                method_idx = int(row['method_index'])-1 # iloc start from 0
                
                if row['projectID'] == 'S':
                    sam_mode = 'SPARTAN'
                elif row['projectID'] == 'M':
                    sam_mode = 'MAIA'
                
                # Get method parameters 
                if 'mass_ug' in col: 
                    if masstype == 'pm25':
                        method_params = PM25_mass_para.iloc[method_idx]
                    elif masstype == 'pm10':
                        method_params = PM10_mass_para.iloc[method_idx]
                    else:
                        logging.ERROR(f'mass type error! current mass type: {masstype}, please use one of these: pm25, pm10')
                    
                    # need to add flag if there is collocated Nylon filter 
                    # flag = 'Collocated Nylon filter exists, mass not added' 
                    
                     
                elif '_XRF_' in col:
                    # identify element
                    ele = col.split('_')[0]
                    ele_ful = element_codes[ele]
                    
                    if masstype == 'pm25':
                        metals_para = PM25_metals_para
                    elif masstype == 'pm10':
                        metals_para = PM10_metals_para
                        
                        
                    method_params = metals_para[
                        (metals_para['Parameter'].str.contains(ele_ful, case=False)) & 
                        (metals_para['Analysis Description'].str.contains('ED-XRF')) &
                        (metals_para['Sampling Mode'] == sam_mode )
                    ].iloc[0]
                    
                    # MDL & UNC 
                    if row['CartridgeID'] in xrf_elog.keys():
                        elog_date = xrf_elog[row['CartridgeID']]
                        future_mdl = xrf_mdl[xrf_mdl['Date'] > elog_date]
                        
                        # Find the minimum date
                        if not future_mdl.empty:
                            next_date = future_mdl['Date'].min() # nearest future mdl data
                            this_mdl = xrf_mdl[xrf_mdl['Date'] == next_date]
                        else:
                            this_mdl = xrf_mdl[xrf_mdl['Date'] == xrf_mdl['Date'].max()] # nearest past mdl data
                        
                        new_row['Analytical_MDL'] = this_mdl.loc[this_mdl['Statistics'] =='AnalyticalMDL (ug/cm2)',ele].values[0]
                        new_row['MDL'] = this_mdl.loc[this_mdl['Statistics'] =='MDL (ug/cm2)',ele].values[0]
                        new_row['UNC'] = this_mdl.loc[this_mdl['Statistics'] =='Uncertainty (ug/cm2)',ele].values[0] 
                    
                    else:
                        # some xrf data are archived due to quality issue. The values are still in master files but not considered for MDL
                        # cartid_temp = row['CartridgeID']
                        # value_temp = row[col]
                        # logging.warning(f'{cartid_temp} not found in XRF elog file, not MDL can be reported. column: {col} xrf value: {value_temp}')
                        new_row['Analytical_MDL'] =''
                        new_row['MDL']=''
                        new_row['UNC'] =''
                        
                        
                elif '_ICP_' in col:
                    # identify element
                    ele = col.split('_')[0]
                    if ele in element_codes: # there are icp element that are not in method list (e.g. Silver)
                        ele_ful = element_codes[ele]
                        
                        if masstype == 'pm25':
                            metals_para = PM25_metals_para
                        elif masstype == 'pm10':
                            metals_para = PM10_metals_para
                            
                        method_params = metals_para[
                                (metals_para['Parameter'].str.contains(ele_ful, case=False)) & 
                                (metals_para['Analysis Description'].str.contains('ICP-MS')) 
                            ]
                        
                        if not method_params.empty:
                            method_params = method_params.iloc[method_idx] 
                        else:
                            # In case there are ICP elements that not reported (not avail in metals_para)
                            continue
                        
                elif 'IC_' in col:
                    # identify element
                    ele = col.split('_')[1]
                    if ele in element_codes: # there are icp element that are not in method list (e.g. Silver)
                        ele_ful = element_codes[ele]
                        
                        if masstype == 'pm25':
                            ic_para = PM25_IC_para
                        elif masstype == 'pm10':
                            ic_para = PM10_IC_para
                        
                        if col.endswith('_T'):
                            method_params = ic_para[
                                (ic_para['Parameter'].str.contains(ele_ful, case=False)) & 
                                (ic_para['Collection Description'].str.contains('Teflon')) &
                                (ic_para['Sampling Mode'] == sam_mode) 
                            ]
                            
                            if not method_params.empty and sam_mode == 'SPARTAN':
                                method_params = method_params.iloc[method_idx] 
                                
                            elif not method_params.empty and sam_mode == 'MAIA':
                                method_params = method_params.iloc[0]
                            else:
                                # There are IC elements that not reported (not avail in IC_para)
                                continue
                            
                            # MDL & UNC 
                            if ic_elog[row['FilterID']] < datetime(2021,11,25):
                                mdl_values = ic_mdl_tef.loc['Before_2021_Nov25_ug',f'IC_conc_{ele}']
                            else:
                                mdl_values = ic_mdl_tef.loc['After_2021_Nov25_ug',f'IC_conc_{ele}']
                            
                            new_row['Analytical_MDL'] = ''
                            new_row['MDL'] = mdl_values
                            new_row['UNC'] = ''
                            
                        elif col.endswith('_N'):
                            if ele == 'NH4': # report estimated NH4 using nitrate
                                method_params = ic_para[
                                    (ic_para['Parameter'].str.contains(ele_ful, case=False)) & 
                                    (ic_para['Collection Description'].str.contains('Estimated Ammonium')) &
                                    (ic_para['Sampling Mode'] == sam_mode) 
                                ]
                                method_params = method_params.iloc[0]
                            else:
                                method_params = ic_para[
                                    (ic_para['Parameter'].str.contains(ele_ful, case=False)) & 
                                    (ic_para['Collection Description'].str.contains('Nylon')) &
                                    (ic_para['Sampling Mode'] == sam_mode) 
                                ]
                                # Check if any rows found match the criteria
                                if not method_params.empty:
                                    method_params = method_params.iloc[0]
                                else:
                                    # There are IC elements that not reported (not avail in IC_para)
                                    continue
                                
                            # MDL & UNC 
                            mdl_values = ic_mdl_nyl[f'IC_{ele}_ug_N'].values[0]
                            
                            new_row['Analytical_MDL'] = ''
                            new_row['MDL'] = mdl_values
                            new_row['UNC'] = ''
                    
                    
                elif '_FTIR_' in col:
                    ele = col.split('_')[0]
                    
                    if masstype == 'pm25':
                        carbon_para = PM25_carbon_para
                    elif masstype == 'pm10':
                        carbon_para = PM10_carbon_para
                        
                    method_params = carbon_para[
                            (carbon_para['Parameter'].str.contains(ele, case=False)) & 
                            (carbon_para['Analysis Description'] == 'FTIR') &
                            (carbon_para['Sampling Mode'] == sam_mode) 
                        ].iloc[0]
                    
                    new_row['Analytical_MDL'] = ''
                    new_row['MDL'] = row[f'{ele}_FTIR_MDL']
                    new_row['UNC'] = ''
            
                elif '_SSR_' in col:
                    ele = col.split('_')[0]
                    
                    if masstype == 'pm25':
                        carbon_para = PM25_carbon_para
                    elif masstype == 'pm10':
                        carbon_para = PM10_carbon_para
                        
                    method_params = carbon_para[
                            (carbon_para['Parameter'].str.contains(ele, case=False)) & 
                            (carbon_para['Analysis Description'].str.contains('Smoke Stain Reflectometer')) 
                        ]
                    
                    method_params = method_params.iloc[method_idx]
                        
                elif 'BC_HIPS_ug' in col:
                    
                    if masstype == 'pm25':
                        carbon_para = PM25_carbon_para
                    elif masstype == 'pm10':
                        carbon_para = PM10_carbon_para
                        
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
                else: # skip columns about site/filter info 
                    continue
                
                # adding value of the measurement
                new_row['Value'] = round(row[col], 2)
                
                # basic site information
                new_row['Site_Code'] = site_info['Site_Code']
                new_row['Latitude'] = site_info['Latitude']
                new_row['Longitude'] = site_info['Longitude']
                new_row['Elevation_meters'] = site_info['Elevation_meters']
                
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
                new_row['Flag'] = copyflag( row['Flags'], Flag_parameters)
                
                # Append the new row to pm25_chemspec
                pm25_chemspec = pd.concat([pm25_chemspec, pd.DataFrame([new_row])], ignore_index=True)
    
    fname = f'{direc_output}/Chemical_Filter_Data/{masstype.upper()}/{site_info.Site_Code}_{masstype.upper()}_speciation.csv'
    # generate metadata lines
    today = date.today().isoformat()
    metadata = [
        f'# File Updated: {today}',
        f'# {DataVersion}',
        f'# Site: {site_info.City}, {site_info.Country}'
    ]
    # Write metadata + DataFrame
    with open(fname, 'w') as f:
        for line in metadata:
            f.write(line + '\n')
        pm25_chemspec.to_csv(f, index=False)

    logging.info(f'{fname} saved')

def get_teo(pm25_data, xrf_exist,icp_exist):
    # TEO ratios if ICP-MS = (1.47[V] + 1.27[Ni] + 1.25[Cu] + 1.24[Zn] + 1.32[As] + 1.2[Se] + 1.07[Ag] + 1.14[Cd] + 1.2[Sb] + 1.12[Ba] + 1.23[Ce] + 1.08[Pb])
    # TEO_ratio_ICP = [1.47, 1.27, 1.25, 1.24, 1.32, 1.2, 1.07, 1.14, 1.2, 1.12, 1.23, 1.08]
    
    # TEO ratios if XRF = (1.79[V] + 1.69[Cr] + 1.63[Mn] + 1.34[Co] + 1.27[Ni] + 1.25[Cu] + 1.24[Zn] + 1.43[As] + 1.41[Se] + 1.09[Rb] + 1.18[Sr] + 1.14[Cd] + 1.20[Sn] + 1.26[Sb] + 1.20[Ce] + 1.12[Pb])
    # TEO_ratio_XRF = [1.79, 1.69, 1.63, 1.34, 1.27, 1.25, 1.24, 1.43, 1.41, 1.09, 1.18, 1.14, 1.20, 1.26, 1.20, 1.12]
    
    # Define coefficients and element mappings
    ICP_COEFFICIENTS = {
        'V': 1.47,
        'Ni': 1.27,
        'Cu': 1.25,
        'Zn': 1.24,
        'As': 1.32,
        'Se': 1.2,
        'Ag': 1.07,
        'Cd': 1.14,
        'Sb': 1.2,
        'Ba': 1.12,
        'Ce': 1.23,
        'Pb': 1.08
    }

    XRF_COEFFICIENTS = {
        'V': 1.79,
        'Cr': 1.69,
        'Mn': 1.63,
        'Co': 1.34,
        'Ni': 1.27,
        'Cu': 1.25,
        'Zn': 1.24,
        'As': 1.43,
        'Se': 1.41,
        'Rb': 1.09,
        'Sr': 1.18,
        'Cd': 1.14,
        'Sn': 1.20,
        'Sb': 1.26,
        'Ce': 1.20,
        'Pb': 1.12
    }
    teo = 0
    
    if xrf_exist:
        for element, coeff in XRF_COEFFICIENTS.items():
            col_name = f"{element}_XRF_ng"
            if col_name in pm25_data.columns:
                teo += pm25_data[col_name].fillna(0).values * coeff
                teo = teo[0]/1000

    elif icp_exist:
        for element, coeff in ICP_COEFFICIENTS.items():
            col_name = f"{element}_ICP_ng"
            if col_name in pm25_data.columns:
                teo += pm25_data[col_name].fillna(0).values * coeff
                teo = teo[0]/1000
    else:
        teo = np.nan
    
    
    return teo


# ========= FIGURES - PM25 time series =========
def save_fig(output_path, figname, close=True):
    os.makedirs(output_path, exist_ok=True)
    figname = os.path.join(output_path, figname)
    plt.savefig(figname, dpi=300, bbox_inches='tight')
    if close:
        plt.close()
    logging.info(f'{figname} saved')


def timeseries_pm25(rcfm_df,fname, savedir, Site_cities):
    PM25_data = rcfm_df['Filter PM2.5 Mass'].values
    # Convert dates to datetime objects
    start_date = rcfm_df['Start_Date']
    end_date = rcfm_df['End_Date']
    
    # Create full date range
    plotting_dates = pd.date_range(start_date.iloc[0] , end_date.iloc[-1])
    
    # Initialize plotting array
    PM_Plotting = pd.DataFrame({'Date': plotting_dates, 'PM25': np.nan}).set_index('Date')

    # Fill PM25 values for each sampling period
    for idx, row in rcfm_df.iterrows():
        start_dt = row.Start_Date
        end_dt = row.End_Date
        PM_Plotting.loc[start_dt:end_dt, 'PM25'] = row['Filter PM2.5 Mass']

    # Calculate average
    Site_PM_Avg = PM_Plotting['PM25'].mean(skipna=True)
    Avg_Line = np.full(len(PM_Plotting), Site_PM_Avg)
    Day_num = len(PM_Plotting)

    # Create plot
    plt.figure(figsize=(6.5, 3.6))  # Approx 500x288 points (1 inch = 72 points)
    ax = plt.gca()

    # Plot data
    plt.scatter(PM_Plotting.index, PM_Plotting['PM25'], marker='s', s=20, linewidth=1)

    # Set x-ticks based on duration
    if Day_num > 365*3:
        months = round( Day_num/(365*3) )
    else:
        months = 1 
        
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=months))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%Y'))

    # Formatting
    plt.ylim(bottom=0)
    plt.ylabel('PM$_{2.5}$ Concentration (μg/m$^3$)', fontweight='bold', fontsize=9)
    plt.xlabel('Date (MMM-YYYY)', fontweight='bold', fontsize=9)
    plt.title(f'Filter-based PM$_{{2.5}}$, {Site_cities}', fontsize=11)
    
    # Add average line
    plt.plot(PM_Plotting.index, Avg_Line, '--', linewidth=1.5, color='chocolate')
    
    # Add legend and grid
    plt.legend(['Filter Data', 'Site Average'], fontsize=8, loc='best')
    plt.grid(True)
    
    # Rotate date labels
    plt.xticks(rotation=45, ha='right', fontsize=7)
    plt.yticks(fontsize=7)
    
    # Adjust layout and save
    plt.tight_layout()
    
    # Save figure
    save_fig(savedir, fname)
    

# ========= FIGURES - Bar charts =========
def make_bar(rcfm_df, spec_mapping, colors, figname, savedir, city, website=False):
    
    spec_order = list(spec_mapping.keys())
    species_cols = list(spec_mapping.values())
    n_bars = len(rcfm_df)
    
    # Calculate means
    means = {spec: rcfm_df[spec].mean() for spec in species_cols}
    
    plt.figure(figsize=(16, 6))
    # Create stacked bars using dataframe columns
    bottom = np.zeros(n_bars)
    bars = []
    
    for idx, spec in enumerate(species_cols):
        bar = plt.bar(
            range(n_bars), 
            rcfm_df[spec].fillna(0),  # Handle any remaining NaNs
            bottom=bottom,
            color=colors[idx,:],
            edgecolor=colors[idx,:],
        )
        bars.append(bar)
        bottom += rcfm_df[spec].fillna(0).values

    # Format plot using dataframe index for x-axis
    ax = plt.gca()
    ax.set_xlim(0.5, n_bars + 0.5)
    ax.set_xticks(range(n_bars))
    
    # First version: filter IDs as labels
    ax.set_xticklabels(rcfm_df['FilterID'], rotation=45, ha='right', fontsize=6)
    ax.tick_params(axis='both', labelsize=6)
    ax.set_xlabel('Filter ID', fontsize=8, fontweight='bold')
    ax.set_ylabel('Attributed Concentration (μg/m$^{\mathbf{3}}$)', fontsize=8, fontweight='bold')
    plt.title(f'{city} Chemical Speciation', fontsize=10)
    
    # Create legend using actual column names
    legend = plt.legend(bars, species_cols, fontsize=6, loc='best')
    plt.grid(False)

    # Save filter ID version
    if website is False:
        output_path = f'{savedir}Bar_spec_plots_filterlabels'
        save_fig(output_path, figname, close=False)

    # Second version: dates as labels
    date_labels = pd.to_datetime(rcfm_df['End_Date']).dt.strftime('%m-%d-%Y')
    tick_positions = np.arange(len(rcfm_df))
    # Set positions and labels
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(date_labels, rotation=45, ha='right', fontsize=6)  # Use actual datetime objects
    ax.set_xlabel('Filter sampling end date (MM-DD-YYYY)', fontsize=8, fontweight='bold')

    # Save date version
    if website is False:
        output_path = f'{savedir}Bar_spec_plots_dates'
    else:
        output_path = f'{savedir}Bar_spec_website'
        
    save_fig(output_path, figname)
        
def get_bar_spec():
        # Define species mapping, colors and order
    spec_mapping = {
        'Sulfate': 'Sulfate',
        'Ammonium': 'Ammonium',
        'Nitrate': 'Nitrate',
        'Sea Salt': 'Sea Salt',
        'Dust': 'Fine Soil',
        'Trace Element Oxides': 'Trace Element Oxides',
        'Black Carbon': 'Equivalent BC PM2.5',
        'Water': 'Particle Bound Water',
        'Organic Carbon': 'Organic Carbon',
        'Residual Matter': 'Residual Matter',
    }

    colors = np.array([
        [237, 48, 41],    # Sulfate (red)
        [240, 103, 166],  # Ammonium (pink)
        [245, 126, 32],   # Nitrate (orange)
        [57, 84, 165],    # Sea Salt (blue)
        [252, 238, 30],   # Fine Soil (yellow)
        [128, 130, 133],  # TEO (grey)
        [35, 31, 32],     # BC (black)
        [109, 207, 246],  # PBW (water/blue)
        [55, 98, 60],     # OC (dark green)
        [80, 184, 72],    # OM (green)
    ]) / 255
    return spec_mapping, colors

def bar_all(rcfm_df,fname, savedir, city):
    spec_mapping, colors = get_bar_spec()
    
    make_bar(rcfm_df, spec_mapping, colors, fname, savedir, city)

def bar_website(rcfm_df,fname, savedir, city):
    if len(rcfm_df) < 6:
        return
    elif len(rcfm_df) < 51:
        rcfm_plot = rcfm_df
    else:
        rcfm_plot = rcfm_df[-50:]
        
    spec_mapping, colors = get_bar_spec()
    
    make_bar(rcfm_plot, spec_mapping, colors, fname, savedir, city, website=True)


# ========= FIGURES - Pie charts =========
def make_pie(rcfm_df, spec_mapping, colors, figname, savedir, city):
    
    spec_order = list(spec_mapping.keys())
    
    # Calculate means
    means = {spec: rcfm_df[spec_mapping[spec]].mean() for spec in spec_order}

    # Handle negative values
    negative_notes = []
    for spec in spec_order:
        if means[spec] < 0:
            negative_notes.append(f"{spec} = {means[spec]:.2f} and is set to 0 in the chart.")
            means[spec] = 0

    # Prepare data for pie chart
    values = [means[spec] for spec in spec_order]
    labels = [f"{spec}: {v:.1f}" for spec, v in zip(spec_order, values)]
    
    # Convert NaNs to 0 and ensure positive values
    values = [max(0, np.nan_to_num(v)) for v in values]

    # Check if all values are zero
    if sum(values) == 0:
        logging.info(f"Skipping pie chart for {figname} - all values zero")
        plt.close()
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(5, 4))

    # Create pie chart
        # Create pie chart without percentage labels
    wedges, texts = ax.pie(
        values, colors=colors, startangle=90,
        wedgeprops={'width': 0.5}
    )

    # Add center circle and mass text
    ax.add_artist(plt.Circle((0, 0), 0.25, fc='white'))
    ax.text(0, 0, f"{rcfm_df['Filter PM2.5 Mass'].mean():.0f}",
            ha='center', va='center', fontsize=42, fontweight='bold')

    # Add legend on the right
    ax.legend(wedges, labels,
              loc="center left",
              bbox_to_anchor=(1, 0, 0.5, 1))

    # Add negative value notes
    if negative_notes:
        plt.figtext(0.7, 0.01, "\n".join(negative_notes),
                    ha='center', fontsize=10, color='red')
    
    # Add a title of City name
    plt.title(city, fontsize=36, fontweight='bold')
    
   # PDF (best for Illustrator; keeps live text if the font is installed)
    os.makedirs(savedir, exist_ok=True)
    plt.savefig(f"{savedir}/{figname}.pdf",
                format="pdf", dpi=300, bbox_inches="tight",
                transparent=False)    # png version
    
     # png version
    figname = f'{figname}.png'
    save_fig(savedir, figname)
 
    
def rcfm_pie(rcfm_df, figname,savedir, city):
    # Define species mapping, colors and order
    spec_mapping = {
        'BC': 'Equivalent BC PM2.5',
        'TEO': 'Trace Element Oxides',
        'Fine Soil': 'Fine Soil',
        'Sea Salt': 'Sea Salt',
        'Nitrate': 'Nitrate',
        'Ammonium': 'Ammonium',
        'Sulfate': 'Sulfate',
        'PBW': 'Particle Bound Water',
        'OC': 'Organic Carbon',
    }

    colors = np.array([
        [35, 31, 32],     # BC (black)
        [128, 130, 133],  # TEO (grey)
        [252, 238, 30],   # Fine Soil (yellow)
        [57, 84, 165],    # Sea Salt (blue)
        [245, 126, 32],   # Nitrate (orange)
        [240, 103, 166],  # Ammonium (pink)
        [237, 48, 41],    # Sulfate (red)
        [109, 207, 246],  # PBW (water/blue)
        [80, 184, 72],     # OC (green)
    ]) / 255
     
    make_pie(rcfm_df, spec_mapping, colors, figname, savedir, city)


def rcfm_wRM(rcfm_df, figname, savedir, city):
    # Define species mapping, colors and order
    spec_mapping = {
        'BC': 'Equivalent BC PM2.5',
        'TEO': 'Trace Element Oxides',
        'Fine Soil': 'Fine Soil',
        'Sea Salt': 'Sea Salt',
        'Nitrate': 'Nitrate',
        'Ammonium': 'Ammonium',
        'Sulfate': 'Sulfate',
        'PBW': 'Particle Bound Water',
        'OC': 'Organic Carbon',
        'RM': 'Residual Matter',
    }

    colors = np.array([
        [35, 31, 32],     # BC (black)
        [128, 130, 133],  # TEO (grey)
        [252, 238, 30],   # Fine Soil (yellow)
        [57, 84, 165],    # Sea Salt (blue)
        [245, 126, 32],   # Nitrate (orange)
        [240, 103, 166],  # Ammonium (pink)
        [237, 48, 41],    # Sulfate (red)
        [109, 207, 246],  # PBW (water/blue)
        [55, 98, 60],     # OC (dark green)
        [80, 184, 72],    # OM (green)
    ]) / 255
     
    make_pie(rcfm_df, spec_mapping, colors, figname, savedir, city)


def rcfm_OM_wRM(rcfm_df, figname, savedir, city):
    # Define species mapping, colors and order
    spec_mapping = {
        'BC': 'Equivalent BC PM2.5',
        'TEO': 'Trace Element Oxides',
        'Fine Soil': 'Fine Soil',
        'Sea Salt': 'Sea Salt',
        'Nitrate': 'Nitrate',
        'Ammonium': 'Ammonium',
        'Sulfate': 'Sulfate',
        'PBW': 'Particle Bound Water',
        'OM': 'Organic Matter',
        'RM': 'Residual Matter',
    }

    colors = np.array([
        [35, 31, 32],     # BC (black)
        [128, 130, 133],  # TEO (grey)
        [252, 238, 30],   # Fine Soil (yellow)
        [57, 84, 165],    # Sea Salt (blue)
        [245, 126, 32],   # Nitrate (orange)
        [240, 103, 166],  # Ammonium (pink)
        [237, 48, 41],    # Sulfate (red)
        [109, 207, 246],  # PBW (water/blue)
        [55, 98, 60],     # OM (dark green)
        [80, 184, 72],    # RM (green)
    ]) / 255
     
    make_pie(rcfm_df, spec_mapping, colors, figname, savedir, city)

       
def rcfm_with_K_Cl(rcfm_df, figname, savedir, city):
    # Define species mapping, colors and order
    spec_mapping = {
        'BC': 'Equivalent BC PM2.5',
        'TEO': 'Trace Element Oxides',
        'Fine Soil': 'Fine Soil',
        'Na': 'Sodium',
        'Cl': 'Chlorine',
        'K': 'Potassium',
        'Nitrate': 'Nitrate',
        'Ammonium': 'Ammonium',
        'Sulfate': 'Sulfate',
        'PBW': 'Particle Bound Water',
        'OC': 'Organic Carbon',
    }

    colors = np.array([
        [35, 31, 32],     # BC (black)
        [128, 130, 133],  # TEO (grey)
        [252, 238, 30],   # Fine Soil (yellow)
        [57, 84, 165],    # Sea Salt (blue)
        [56, 170, 165],   # Cl
        [180, 132, 212],   # purple K
        [245, 126, 32],   # Nitrate (orange)
        [240, 103, 166],  # Ammonium (pink)
        [237, 48, 41],    # Sulfate (red)
        [109, 207, 246],  # PBW (water/blue)
        [80, 184, 72],     # OC (green)
    ]) / 255
    
    make_pie(rcfm_df, spec_mapping, colors, figname, savedir, city )
    
           
def rcfm_with_K_Cl_wRM(rcfm_df, figname, savedir, city):
    # Define species mapping, colors and order
    spec_mapping = {
        'BC': 'Equivalent BC PM2.5',
        'TEO': 'Trace Element Oxides',
        'Fine Soil': 'Fine Soil',
        'Na': 'Sodium',
        'Cl': 'Chlorine',
        'K': 'Potassium',
        'Nitrate': 'Nitrate',
        'Ammonium': 'Ammonium',
        'Sulfate': 'Sulfate',
        'PBW': 'Particle Bound Water',
        'OC': 'Organic Carbon',
        'RM': 'Residual Matter',
    }

    colors = np.array([
        [35, 31, 32],     # BC (black)
        [128, 130, 133],  # TEO (grey)
        [252, 238, 30],   # Fine Soil (yellow)
        [57, 84, 165],    # Sea Salt (blue)
        [56, 170, 165],   # Cl
        [180, 132, 212],   # purple K
        [245, 126, 32],   # Nitrate (orange)
        [240, 103, 166],  # Ammonium (pink)
        [237, 48, 41],    # Sulfate (red)
        [109, 207, 246],  # PBW (water/blue)
        [55, 98, 60],     # OC (dark green)
        [80, 184, 72],     # RM (green)
    ]) / 255
    

    make_pie(rcfm_df, spec_mapping, colors, figname, savedir, city )
    
           
def rcfm_OM_with_K_Cl_wRM(rcfm_df, figname, savedir, city):
    # Define species mapping, colors and order
    spec_mapping = {
        'BC': 'Equivalent BC PM2.5',
        'TEO': 'Trace Element Oxides',
        'Fine Soil': 'Fine Soil',
        'Na': 'Sodium',
        'Cl': 'Chlorine',
        'K': 'Potassium',
        'Nitrate': 'Nitrate',
        'Ammonium': 'Ammonium',
        'Sulfate': 'Sulfate',
        'PBW': 'Particle Bound Water',
        'OM': 'Organic Matter',
        'RM': 'Residual Matter',
    }

    colors = np.array([
        [35, 31, 32],     # BC (black)
        [128, 130, 133],  # TEO (grey)
        [252, 238, 30],   # Fine Soil (yellow)
        [57, 84, 165],    # Sea Salt (blue)
        [56, 170, 165],   # Cl
        [180, 132, 212],   # purple K
        [245, 126, 32],   # Nitrate (orange)
        [240, 103, 166],  # Ammonium (pink)
        [237, 48, 41],    # Sulfate (red)
        [109, 207, 246],  # PBW (water/blue)
        [55, 98, 60],     # OM (dark green)
        [80, 184, 72],     # RM (green)
    ]) / 255

    make_pie(rcfm_df, spec_mapping, colors, figname, savedir, city )