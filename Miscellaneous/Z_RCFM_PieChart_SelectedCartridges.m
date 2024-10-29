%% %%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: to read master data files, examine mass and ions and make pie charts for RCFM
% Data validation steps are similar to script F but no filter will be excluded due to quality issues 
% No public files will be created. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc
addpath('../UtilityFunctions')

%% %%%%%%% USER SWITCHES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DataVersion = 'Data version 2.4'; % will be in the first line of the public files

% Setup directories 
debug_mode = 0;
direc = find_root_dir(debug_mode);

direc_master = strcat(direc,'/Analysis_Data/Master_files');
direc_output =  strcat(direc,'/Public_Data');

% ---- the diary function saves the processing history into monthly records ----
diary(sprintf('%s/Public_Data/Data_Processing_Records/Master_File/%s_MasterFile_Processing_Record.txt',direc, datestr(now,'yyyy-mm-dd-HHMMSS')))
fprintf('%s \n', datestr(now))

% ---- Load Site Details and Sampling Methods ----
Sampling_Parameters_Methods = strcat(direc,'/SOPs/Public SOPs/Sampling_Parameters_Methods_2.3.xlsx');

PM25_mass_para = readtable(Sampling_Parameters_Methods,'Sheet','PM2.5 mass','PreserveVariableNames',true);
PM25_IC_para   = readtable(Sampling_Parameters_Methods,'Sheet','PM2.5 water-soluble ions','PreserveVariableNames',true);
PM25_metals_para = readtable(Sampling_Parameters_Methods,'Sheet','PM2.5 trace elements','PreserveVariableNames',true);
PM25_carbon_para = readtable(Sampling_Parameters_Methods,'Sheet','PM2.5 BC_OC_EC','PreserveVariableNames',true);

PM10_mass_para = readtable(Sampling_Parameters_Methods,'Sheet','PM10 mass','PreserveVariableNames',true);
PM10_IC_para   = readtable(Sampling_Parameters_Methods,'Sheet','PM10 water-soluble ions','PreserveVariableNames',true);
PM10_metals_para = readtable(Sampling_Parameters_Methods,'Sheet','PM10 trace elements','PreserveVariableNames',true);
PM10_carbon_para = readtable(Sampling_Parameters_Methods,'Sheet','PM10 BC_OC_EC','PreserveVariableNames',true);

PMc_parameters  = readtable(Sampling_Parameters_Methods,'Sheet','PMc','PreserveVariableNames',true);
RCFM_parameters = readtable(Sampling_Parameters_Methods,'Sheet','PM2.5 RCFM','PreserveVariableNames',true);

Flag_parameters = readtable(Sampling_Parameters_Methods,'Sheet','Flags','PreserveVariableNames',true);

site_details = readtable(strcat(direc,'/Site_Sampling/Site_details.xlsx'),'PreserveVariableNames',true);
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));
Site_countries = table2array(site_details(:,2));
latitudes = table2array(site_details(:,5));
longitudes = table2array(site_details(:,6));
elevations = table2array(site_details(:,7));

% ---- addtional output initialization ----
HIPS_Mesured_BC = nan(size(Site_codes));
HIPS_Est_BC = nan(size(Site_codes));
HIPS_total_BC = nan(size(Site_codes));
SSR_total_BC = nan(size(Site_codes));

% calculate reportable data #:
DataNum_tab_rows = {'City','PM2.5','Water','SO4','NH4','NO3','SS','Dust','TEO','BC','OC','RM'};
tSiteDataNum = nan(length(Site_codes),length(DataNum_tab_rows)-1);

% ---- Load the MDL & Unc matrix for XRF ----  
MDLfname = sprintf('%s/Analysis_Data/XRF/MDL_Unc_Ref.mat',direc);
load(MDLfname,'Ref_Dates','Ref_Values','FilterGroup')

% ---- Load MAL and CF for dust estimation
MAL_all = site_details.MAL;
CF_all = site_details.CF;

% ---- Load the list of bad sodium cartridges ----
BadNa = sprintf('%s/Analysis_Data/Ion_Chromatography/Bad_sodium_data_in_2017.xlsx',direc);
opts = detectImportOptions(BadNa);
opts.VariableNames = {'CartID'};
BadNa = table2array(readtable(BadNa,opts));

% Mass type variable meaning:
% 0 = cartridge (traveling) blank
% 1 = PM2.5
% 2 = PM10
% 3 = PMcoarse
% 4 = saturated nuclepore, filter measurements invalid
% 5 = negative mass, filter measurements invalid
% 6 = invalid flow, filter measurements invalid

% fileID = fopen('Zero_values.txt', 'w');
% if fileID == -1
%     error('File could not be opened for appending');
% end
% fprintf(fileID,'Label    S_year S_month S_day S_hour E_year E_month E_day E_hour Parameter Value  Description\n');
% fclose(fileID);

%% Read master files and process public files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  for loc =  1:numel(Site_codes)
 for loc =  find(ismember(Site_codes,'CHTS'))
    %% Read Master Files
    master_file = sprintf('%s/%s_master.csv',direc_master,Site_codes{loc});
    [Titles,Master_IDs,   Master_Barcodes, Master_CartridgeIDs, Master_LotIDs, Master_projectIDs,  Master_hours,  Master_masstype, ...
            Master_dates, Master_mass,     Master_IC,           Master_ICP,    Master_XRF,...
            Master_carbon, Master_Method,  Master_flags] = ReadMaster(master_file,Site_codes{loc} );
    
    Vol_Col = 2; % 2nd col in mass matrix is volume
    IC_title = Titles(contains(Titles,'IC_')); 
    ICP_title = Titles(contains(Titles,'ICP_')); 
    XRF_title = Titles(contains(Titles,'XRF_')); 

    if isempty(Master_mass)
        continue
    end
    
    % Read IC MDL Files
    ic_mdl_file = sprintf('%s/Analysis_Data/Ion_Chromatography/MDL/%s_IC_MDL.xlsx',direc,Site_codes{loc});
    if exist(ic_mdl_file,'file')
        ic_mdl = readtable(ic_mdl_file);
        ic_mdl_filterid = ic_mdl.Filter_ID;
        ic_mdl_titles = ic_mdl.Properties.VariableNames;
        ic_mdl_titles(6:end) = {'Fluoride' 'Chloride' 'Nitrite'  'Bromide' 'Nitrate' 'Phosphate' 'Sulfate'...
            'Lithium' 'Sodium' 'Ammonium' 'Potassium' 'Magnesium' 'Calcium'};
    end
    %% Data validation step 1: Corrections for ICP-MS measurement
    %  Prcedure:
    %  1. sampled filters - blank = corrected trace elements;
    
    Cart_ALL = unique(Master_CartridgeIDs); % IDs of cartridges in master file
    
    for i = 1:length(Cart_ALL)
        Cart_IND = find(strcmp(Master_CartridgeIDs, Cart_ALL(i)));
        Cart_BL = find(Master_masstype(Cart_IND) == 0); % finds the blanks for the cartridge
        if ~isempty(Cart_BL)
            if length(Cart_IND)  == 16 && isempty(Cart_BL) == 0 % for cartridges with 2 blanks, one teflon and one nuclepore, in positions 8 and 16
                Cart_BL = [8 16];
                Master_ICP(Cart_IND(Cart_BL(1)-7):Cart_IND(Cart_BL(1) - 1),:) = Master_ICP(Cart_IND(Cart_BL(1)-7):Cart_IND(Cart_BL(1) - 1),:) - Master_ICP(Cart_IND(Cart_BL(1)),:); % corrects PM2.5
                Master_ICP(Cart_IND(Cart_BL(2)-7):Cart_IND(Cart_BL(2) - 1),:) = Master_ICP(Cart_IND(Cart_BL(2)-7):Cart_IND(Cart_BL(2) - 1),:) - Master_ICP(Cart_IND(Cart_BL(2)),:); % corrects PMcoarse
                
            elseif ismember(Site_codes{loc} ,{'SGSU'}) && Cart_BL == 8 % will only be one blank (teflon) and in position 8
                Master_ICP(Cart_IND(Cart_BL-7):Cart_IND(Cart_BL - 1),:) = Master_ICP(Cart_IND(Cart_BL-7):Cart_IND(Cart_BL - 1),:) - Master_ICP(Cart_IND(Cart_BL),:); % only PM2.5 at this site
                
            elseif length(Cart_IND) == 8 && isempty(Cart_BL) == 0 % SS5 station with blanks in filter slot 7
                Master_ICP(Cart_IND(Cart_BL(end)-6):Cart_IND(Cart_BL(end) - 1),:) = Master_ICP(Cart_IND(Cart_BL(end)-6):Cart_IND(Cart_BL(end) - 1),:) - Master_ICP(Cart_IND(Cart_BL(end)),:); % corrects PM2.5 filters
                Master_ICP(Cart_IND(Cart_BL(end)+1),:) = Master_ICP(Cart_IND(Cart_BL(end)+1),:) - Master_ICP(Cart_IND(Cart_BL(end)),:); % corrects PM10 filter
                
            end
        end
        clear Cart_BL Cart_IND
    end
    
    clear Cart_ALL

%     for i = 1:size(Master_ICP,2) % ICP-MS trace metal correction
%         Master_ICP(Master_ICP(:,i) <0,i) = 0; % set any negative metals to zero (due to correction above)
%     end

    %% Data validation step 2 [SKIP]: Exclude bad filters or filter with invalid sampling hour, method, or volume
    %{
    Ind = find(contains(Master_flags,{'Need to be excluded'})==1);
    if ~isempty(Ind)
        fprintf('Following filters excluded for public files as indicated by filter flags:\n')
        disp(Master_IDs(Ind))

          [ Master_IDs, Master_Barcodes, Master_CartridgeIDs, Master_LotIDs, Master_projectIDs,...
            Master_hours, Master_masstype, Master_dates, Master_mass, Master_IC, ...
            Master_ICP, Master_XRF, Master_carbon, Master_Method, Master_flags] = ...
            RemoveRows (Ind, ...
            Master_IDs, Master_Barcodes, Master_CartridgeIDs, Master_LotIDs, Master_projectIDs,...
            Master_hours, Master_masstype, Master_dates, Master_mass, Master_IC, ...
            Master_ICP, Master_XRF, Master_carbon, Master_Method, Master_flags);
    end

    % find valid cols: vol exist and sampled hr ~=0
    vol_idx = find(Master_mass(:,Vol_Col) ~= 0 & ~isnan(Master_mass(:,Vol_Col)) & Master_hours ~=0 & ~isnan(Master_hours));
     
    [Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
    Master_hours, Master_masstype,   Master_dates,        Master_mass,   Master_IC,    ...
    Master_ICP,   Master_XRF,        Master_carbon,       Master_Method, Master_flags]...
     = SortData (vol_idx,...
    Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
    Master_hours, Master_masstype,   Master_dates,        Master_mass,   Master_IC,    ...
    Master_ICP,   Master_XRF,        Master_carbon,       Master_Method, Master_flags);

    clear vol_idx
 
    % remove method is NaN (SSe_ID is NaN, refer to script B)
    non_nan_idx = find(~isnan(Master_Method(:,1)));
 
    [Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
    Master_hours, Master_masstype,   Master_dates,        Master_mass,   Master_IC,    ...
    Master_ICP,   Master_XRF,        Master_carbon,       Master_Method, Master_flags]...
     = SortData (non_nan_idx,...
    Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
    Master_hours, Master_masstype,   Master_dates,        Master_mass,   Master_IC,    ...
    Master_ICP,   Master_XRF,        Master_carbon,       Master_Method, Master_flags);

    clear non_nan_idx
 %}
    %% Data validation step 3: Calculate BC when HIPS-BC not available
    %  Set SSR Bc = -899 to NaN, set negative BC to 0.
    Master_carbon(Master_carbon(:,1) == -899,1) = NaN;
    Master_carbon(Master_carbon(:,1) < 0, 1) = 0; % set any negative BC to zero due to using blank R in calculation [rare]
    
    Ind = find( isnan(Master_carbon(:,2)) &  Master_masstype<3 & Master_carbon(:,1)>0); % HIPS not available & Teflon filters & SSR-BC available
    HIPS_Mesured_BC(loc,1) = sum(~isnan(Master_carbon(:,2)));
    HIPS_Est_BC (loc,1)    = length(Ind);
    % apply the curve:
    a1=2163.54158502065;
    b1=0.00198619072508564;
    c1=-0.00467334844964624;
    a2=262.938831079902;
    b2=0.0123767218647736;
    c2=3.10315712581732;
    Master_carbon(Ind,2) =  a1*sin(b1* Master_carbon(Ind,1) + c1) + a2*sin(b2* Master_carbon(Ind,1) + c2); % SSR-BC unit = ug; HIPS-BC unit = ug;
    % add flag to those filters
    for ii = 1:length(Ind)
        Master_flags{Ind(ii)} = AddFlag (Master_flags{Ind(ii)}, 'HIPS-BC estimated using SSR-BC');
    end

    HIPS_total_BC(loc,1) = sum( ~isnan(Master_carbon(:,2)) );
    SSR_total_BC(loc,1)  = sum( Master_masstype<3 & Master_carbon(:,1)>0 ); 

    %% Data validation step 4: Mask out invalid IC data: bad correlation coefficient or bad Na conc in blanks
    
    % ----- Mark data as NaN when flags indicate invalid data -----
    ions = {'F','Cl','NO2','Br','NO3','PO4','SO4','Li','Na','NH4','K','Mg','Ca'};
    
    for k = 1:length(ions)
        for j = 1:size(Master_IC,1) 
            if contains(Master_flags(j,:),sprintf('S-%s',ions{k})) == 1 % the correlation coefficient for species i was below the required limit.
                Master_IC(j,k) = NaN; %set value to NaN because violated QA protocol
            end
            if contains(Master_flags(j,:),sprintf('Sx-%s',ions{k})) == 1 % The correlation coefficient for species i was missing from IC data file
                Master_IC(j,k) = NaN; %set value to NaN because violated QA protocol
            end
%            if contains(Master_flags(j,:),sprintf('qc1-%s',ions{i})) == 1% The average quality control standard check for species i exceeded +/- 10 % of specified concentration.
%                Master_IC(j,k) = NaN; %set value to NaN because violated QA protocol
%            end
        end
    end
    
    % find cartridges with bad Na (Na too high in blanks)
    Ind = find(contains(BadNa,Site_codes{loc}));
    if ~isempty(Ind)
        NaInd = find(contains(IC_title,'IC_Na'));
        for ii = 1:length(Ind)
            tCartridge = BadNa{Ind(ii)};
            MasterInd = find(ismember(Master_CartridgeIDs,tCartridge));
            Master_IC(MasterInd,NaInd) = NaN; % instead of NaN, using 0 so that no SS will be reported but ASO4 won't be affected.
        end
    end
    clear j k

    
    %% Data validation step 5: Sort data according to dates  
%     
%     % sort dates appropriately 
%     temp_dates = datenum(Master_dates(:,1),Master_dates(:,2),Master_dates(:,3),0,0,0); % create datenum for sorting data based on start dates
%     [dates_sort, idx] = sort(temp_dates);
% 
%     [Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
%     Master_hours, Master_masstype,   Master_dates,        Master_mass,  Master_IC,    ...
%     Master_ICP,   Master_XRF,        Master_carbon,       Master_Method,Master_flags]...
%      = SortData (idx,...
%     Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
%     Master_hours, Master_masstype,   Master_dates,        Master_mass,  Master_IC,    ...
%     Master_ICP,   Master_XRF,        Master_carbon,       Master_Method,Master_flags);
% 
%     clear temp_dates idx dates_sort

    %% Data validation step 6: Convert mass to concentration (ug/m3); change unit of metals from ng to ug 
    if isempty(Master_mass) == 1
        disp('No Mass data detected')
        continue
    else
        % converting mass to concentration
        Master_dates = round(Master_dates); % need to round to nearest hour
        Master_mass(:,1) = Master_mass(:,1) ./Master_mass(:,Vol_Col);
        Master_IC = Master_IC./Master_mass(:,Vol_Col);
        Master_ICP = Master_ICP./Master_mass(:,Vol_Col);
        Master_XRF = Master_XRF./Master_mass(:,Vol_Col);
        Master_carbon = Master_carbon./Master_mass(:,Vol_Col);
    end
    % convert units of metals from ng to ug
    Master_XRF = Master_XRF./1000; Master_ICP = Master_ICP./1000;


    %% Data validation step 7: Compare mass to sum of measured species 
    % remove filter if sum of species greater than collected mass by more than 10%
    % Mass =  Metals + IC (without F, Cl, K, Mg, Ca) + BC 
    IC_col = find(contains(IC_title,{'IC_NO2_ug','IC_Br_ug','IC_NO3_ug',...
             'IC_PO4_ug','IC_SO4_ug','IC_Li_ug','IC_Na_ug','IC_NH4_ug'}));
     
    k = 0;
    
    for i = 1:length(Master_IDs)
        massi = Master_mass(i,1); 

        Master_sumspec = SumSpec(i,Master_XRF,Master_ICP,Master_IC,Master_carbon,IC_col);
        
        if Master_sumspec > massi
            fprintf('Sum of species (%.1f %s) for filter %s greater than collected mass (%.1f %s) \n', Master_sumspec,'ug/m3',Master_IDs{i},massi,'ug/m3')
        end
        
        if (Master_sumspec - massi) > 0.1*massi % test if difference between sum of species and collected mass is > 10% of collected mass
            k = k +1;
            fprintf('Filter mass for %s is invalid: exceeding 10%% allowance.\n', Master_IDs{i})
            filter_index(k) = i;
        end
    end
    
    if exist('filter_index','var') ~=0 % remove rows containing invalid mass data 
        [Master_IDs, Master_Barcodes,Master_CartridgeIDs, Master_LotIDs, Master_projectIDs,...
         Master_hours,  Master_masstype,Master_dates,        Master_mass,   Master_IC,  ...
         Master_ICP,    Master_XRF,     Master_carbon,       Master_Method, Master_flags] = ...
         RemoveRows (filter_index, ...
         Master_IDs, Master_Barcodes,Master_CartridgeIDs, Master_LotIDs, Master_projectIDs,...
         Master_hours,  Master_masstype,Master_dates,        Master_mass,   Master_IC,  ...
         Master_ICP,    Master_XRF,     Master_carbon,       Master_Method, Master_flags);
    end

    clear Master_sumspec massi i k filter_index 
    
    %% Formating: Public data, change the unit of metals back to ng 

    SNAcol  =  [find(contains(IC_title,'IC_SO4_ug')==1);
                find(contains(IC_title,'IC_NO3_ug')==1);
                find(contains(IC_title,'IC_NH4_ug')==1)];

    otherICcol=[find(contains(IC_title,'IC_Na_ug')==1);
                find(contains(IC_title,'IC_PO4_ug')==1);
                find(contains(IC_title,'IC_NO2')==1);
                find(contains(IC_title,'IC_Br')==1);
                find(contains(IC_title,'IC_K_ug')==1);
                find(contains(IC_title,'IC_Mg')==1);
                find(contains(IC_title,'IC_Ca')==1)];

    % Rearrange data to write to public (exclude IC_F, IC_Cl, IC_Li)
    Public_data = [ Master_mass(:,1), Master_IC(:,SNAcol), Master_carbon(:,2), Master_IC(:,otherICcol),...
                    Master_ICP, Master_XRF, Master_carbon(:,4), Master_Method];
    Public_data(:,[1:4 6:10 12]) = round(Public_data(:,[1:4 6:10 12]),2); % mass and water-soluble ions minus Mg round to 2 digits after decimal
    Public_data(:,5) = round(Public_data(:,5),2); % BC round to 1 digit after decimal
    Public_data(:,11) = round(Public_data(:,11)*1000,1); % convert Mg to nanograms/m3 and round to one digit after the decimal
    Public_data(:,13:end-2) = round(Public_data(:,13:end-2)*1000,2); % convert trace elements back to nanograms/m3 and round to one digit after the decimal (at ~ line 219 it was converted from ng to ug)
    
    % Titles for the public data:
    Title_public = {'PM25','IC_SO4_ug','IC_NO3_ug','IC_NH4_ug','BC_ug','IC_Na_ug','IC_PO4_ug','IC_NO2_ug','IC_Br_ug','IC_K_ug','IC_Mg_ng','IC_Ca_ug'};
    % column 13 to later will change if there is change in columns in master file
    Title_public(end+1:end+length(ICP_title)) = ICP_title;
    Title_public(end+1:end+length(XRF_title)) = XRF_title;
    Title_public(end+1:end+2) = {'OC_FTIR_ug','method_index'}; % 'BC_SSR_ug','EC_FTIR_ug' are removed from public data
    
    
   
    %% assign all filters as PM2.5. Will pick out PM10 and FB later
    PM25_index = 1:length(Master_masstype);
    
    if isempty(PM25_index) % No PM2.5 data then no need to go further 
        fprintf('No valid PM25 data for %s\n',Site_codes{loc})
        continue
    end

    % PM2.5
    [PM25_labels, PM25_Barcodes, PM25_cartridgeIDs, PM25_LotIDs, PM25_projectIDs,...
     PM25_hours,  PM25_masstype, PM25_dates,        PM25_Mass,   PM25_IC,   ...
     PM25_ICP,    PM25_XRF,      PM25_Carbon,       PM25_Method, PM25_flags]...
     = SortData (PM25_index,...
    Master_IDs,   Master_Barcodes,   Master_CartridgeIDs, Master_LotIDs,       Master_projectIDs, ...
    Master_hours, Master_masstype,   Master_dates,        Master_mass,  Master_IC,    ...
    Master_ICP,   Master_XRF,        Master_carbon,       Master_Method,Master_flags);
    
    PM25_data_public = Public_data(PM25_index,:);

    %% Write to public file 1: XXXX_speciation.csv
    % components contained in PM25_data_public matrix for searching the parameters sheet
    spec_components = {'Filter PM2.5 mass'...
        'Sulfate' 'Nitrate' 'Ammonium' 'Equivalent BC' 'Sodium' ... % SNA + HIPS_BC + Na 2-6
        'Phosphate' 'Nitrite' 'Bromide' 'Potassium','Magnesium' 'Calcium'... % IC 7-12
        'Lithium' 'Magnesium' 'Aluminum' 'Phosphorus' 'Titanium' ... % ICP 13-17
        'Vanadium' 'Chromium' 'Manganese' 'Iron' 'Cobalt' ... % ICP 18-22
        'Nickel' 'Copper' 'Zinc' 'Arsenic' 'Selenium' 'Gold' ... % ICP 23-28
        'Cadmium' 'Antimony' 'Barium' 'Cerium' 'Lead' ... % ICP 29-33
        'Sodium' 'Aluminum' 'Silicon' 'Sulfur' 'Chlorine' 'Potassium' ... % XRF 34-39
        'Calcium' 'Titanium' 'Vanadium' 'Iron' 'Zinc' 'Cerium' 'Lead' ... % XRF 40-46
        'Arsenic' 'Cobalt' 'Chromium' 'Copper' 'Magnesium' 'Manganese'... % XRF 47-52
        'Nickel' 'Antimony' 'Rubidium' 'Strontium' 'Cadmium' 'Selenium' 'Tin' ... % XRF 53-59
        'OC'}; % FTIR 60
    % PMc_data_public contains 4 more cols [2 carbons 1 method]


    %% Formating: Reconstruct PM2.5  
    
     missing_PM25 = find( isnan(PM25_data_public(:,1)) | PM25_data_public(:,1)<=0 );
    
     temp = NaN.*PM25_hours;
    [PM25_labels, PM25_Barcodes,PM25_cartridgeIDs, PM25_LotIDs, PM25_projectIDs,...
     PM25_hours,  PM25_dates,   PM25_flags,  PM25_data_public,  ...
     ~,  ~,    ~,     ~,   ~,    ~ ] = ...
     RemoveRows (missing_PM25, ...
     PM25_labels, PM25_Barcodes,PM25_cartridgeIDs, PM25_LotIDs, PM25_projectIDs,...
     PM25_hours,  PM25_dates,   PM25_flags,  PM25_data_public, ...
     temp, temp, temp, temp, temp, temp );

           % Species  = [ASO4 NH4NO3 RM   OC NaCl Soil BC  TEO  Na2SO4]
              kappa   = [0.61 0.61 0.10 0.10 1.5  0.01 0.01 0.01 0.68];   % GC: need to double check 
              Density = [1.70 1.70 1.30 0.10 2.2  2.2  1.8  2.5  2.67];   % GC 
              RH = 35;
              MGrowth = 1+kappa./Density.*RH/(100-RH); % mass growth factor
              MGrowth (1:2) = 1.1 ; % halve growth factor of SNA to be consistent with GC simulation. Ref:
              % http://wiki.seas.harvard.edu/geos-chem/index.php/Particulate_matter_in_GEOS-Chem

    % ==== Na2SO4 (sea-salt associated sulfate) ====
    NaSO4_dry=0.18.*PM25_data_public(:,6); % 0.18*[Na+]
    NaSO4_wet=NaSO4_dry*MGrowth(8);
    
     % === ANO3 ====
    ANO3_dry = nan(size(PM25_data_public,1),1);
    Ind = PM25_data_public(:,4)>0.29*PM25_data_public(:,3); % if NH4 is sufficient
    ANO3_dry(Ind) = 1.29*PM25_data_public(Ind,3); % [ANO3(dry)]=1.29*[NO3]
    % NH4 insufficient or one of them is NaN
    ANO3_dry(~Ind) = sum(PM25_data_public(~Ind,[3 4]),2,'omitnan');  % [ANO3(dry)]= [NH4] + [NO3]
    % mask out when both NH4 and NO3 are NaN
    ANO3_dry(sum(isnan(PM25_data_public(:,[3 4])),2)==2) = NaN;
    ANO3_wet=ANO3_dry*MGrowth(2);
    
    % ==== ASO4 ====
    SO4inNaSO4 = 0.12*PM25_data_public(:,6); % "0.12" factor from Henning et al 2003 JGR 108(D1) pp4030, doi:10.1029/2002JD002439
    SO4inNaSO4(isnan(SO4inNaSO4)) = 0;
    % SO4_NSS = [SO4] - 0.12*[Na+]
    SO4_NSS  = sum( [PM25_data_public(:,2) -SO4inNaSO4],2); clear SO4inNaSO4
    SO4_NSS(SO4_NSS<0)=0; % SO4_NSS should never be negative
    %[ASO4(dry)] = [SO42-] + [NH4+]-(18/62*[NO3-])
    NH4_remain = sum([PM25_data_public(:,4) -0.29*PM25_data_public(:,3)],2,'omitnan'); % remaining NH4 after bonded with NIT; if no NO3 (NaN value means did not pass QA/QC), will give NaN
    NH4_remain(NH4_remain<0) = 0; % remaining NH4 should never be negative. In a NIT exceeding event, HNO3 exist.
    ASO4_dry = SO4_NSS+NH4_remain; % low NIT condition, SO4 exists as (NH4)2SO4; high NIT condition, it can be NH4HSO4 or even H2SO4.
    ASO4_wet = ASO4_dry*MGrowth(1); 
    clear NH4_remain
    
    % ==== Sea Salt ====
    % Sea salt(dry) = (2.54*[Na+] - 0.1[Al]) or (1.65[Cl] - 0.1[Al])
    for i = 1:size(PM25_data_public,1)
        if ~isnan(PM25_data_public(i,contains(Title_public,'Al_XRF')))  % check for presence of XRF data first, column 35 is aluminum
            SSalt_dry(i,1) = 2.54*PM25_data_public(i,6)-0.1*(PM25_data_public(i,contains(Title_public,'Al_XRF'))/1000);
        elseif ~isnan(PM25_data_public(i,contains(Title_public,'Al_ICP'))) % check for presence of Al from ICP-MS
            SSalt_dry(i,1) = 2.54*PM25_data_public(i,6)-0.1*(PM25_data_public(i,contains(Title_public,'Al_ICP'))/1000);
        else % no Al from ICP or XRF, therefore do not correct with Al
            SSalt_dry(i,1) = 2.54*PM25_data_public(i,6);
        end
    end
    clear i
    SSalt_dry(SSalt_dry<-0.3)=NaN;  % nan out if too negative

    SSalt_wet = SSalt_dry*MGrowth(4);
     
    % ==== Dust & TEO ====
    % if XRF data available, define soil (aka mineral dust) by combination of {Al, Si, Ca, Ti, Fe} following Xuan Liu's eqn:
    % [SOIL] = [1.89Al×(1+MAL)+2.14Si+1.40Ca+1.36Fe+1.67Ti]×CF

    % if no XRF, only ICP-MS data available, soil is defined as:
    % [SOIL] = 10*([Al] + [Fe] + [Mg]), if no Fe use: [SOIL] = 10*([Mg] + 30*[Ti] + [Al])
    
    % soil is assumed to be hydrophobic, therefore no water uptake
    xrf_soil_factors = [1.89 2.14 1.40 1.36 1.67];
    Soil = nan(size(PM25_data_public,1),1);

    MAL = MAL_all(loc,1);
    CF =  CF_all(loc,1);
    if isempty(MAL)
        error('MAL and CF not found for %s in Site_Details.xlsx',Site_codes{loc})
    end

    % TEO ratios if ICP-MS = (1.47[V] + 1.27[Ni] + 1.25[Cu] + 1.24[Zn] + 1.32[As] + 1.2[Se] + 1.07[Ag] + 1.14[Cd] + 1.2[Sb] + 1.12[Ba] + 1.23[Ce] + 1.08[Pb]);
    TEO_ratio_ICP = [1.47 1.27 1.25 1.24 1.32 1.2 1.07 1.14 1.2 1.12 1.23 1.08];
    TEO_title_ICP =[find(contains(Title_public,'V_ICP')==1),find(contains(Title_public,'Ni_ICP')==1),find(contains(Title_public,'Cu_ICP')==1),...
                    find(contains(Title_public,'Zn_ICP')==1),find(contains(Title_public,'As_ICP')==1),find(contains(Title_public,'Se_ICP')==1),...
                    find(contains(Title_public,'Ag_ICP')==1),find(contains(Title_public,'Cd_ICP')==1),find(contains(Title_public,'Sb_ICP')==1),...
                    find(contains(Title_public,'Ba_ICP')==1),find(contains(Title_public,'Ce_ICP')==1),find(contains(Title_public,'Pb_ICP')==1)];
    % TEO ratios if XRF = (1.79[V] + 1.69[Cr] + 1.63[Mn] + 1.34[Co] + 1.27[Ni] + 1.25[Cu] + 1.24[Zn] + 1.43[As] + 1.41[Se] + 1.09[Rb] + 1.18[Sr] + 1.14[Cd] + 1.20[Sn] + 1.26[Sb] + 1.20[Ce] + 1.12[Pb]);
    TEO_ratio_XRF = [1.79 1.69 1.63 1.34 1.27 1.25 1.24 1.43 1.41 1.09 1.18 1.14 1.20 1.26 1.20 1.12];
    TEO_title_XRF =[find(contains(Title_public,'V_XRF')==1),find(contains(Title_public,'Cr_XRF')==1),find(contains(Title_public,'Mn_XRF')==1),...
                    find(contains(Title_public,'Co_XRF')==1),find(contains(Title_public,'Ni_XRF')==1),find(contains(Title_public,'Cu_XRF')==1),...
                    find(contains(Title_public,'Zn_XRF')==1),find(contains(Title_public,'As_XRF')==1),find(contains(Title_public,'Se_XRF')==1),...
                    find(contains(Title_public,'Rb_XRF')==1),find(contains(Title_public,'Sr_XRF')==1),find(contains(Title_public,'Cd_XRF')==1),...
                    find(contains(Title_public,'Sn_XRF')==1),find(contains(Title_public,'Sb_XRF')==1),find(contains(Title_public,'Ce_XRF')==1),...
                    find(contains(Title_public,'Pb_XRF')==1)];
    TEO = nan(size(PM25_data_public,1),1);
    
    
    % Start to calculate soil and TEO
    no_metals = 0;
    method_metals = ones(size(PM25_data_public,1),1); % default is ICP-MS (ICP-MS = 1, XRF = 2)
    
    for i = 1:size(PM25_data_public,1)

        if ~isnan(PM25_data_public(i,34)) % there is XRF data
            tSoil = xrf_soil_factors.*PM25_data_public(i,[35 36 40 41 43]);
            tSoil(1) = tSoil(1).*(1+MAL); 
            tSoil = tSoil.*CF;
            Soil(i,1) = (sum(tSoil,2,'omitnan'))/1000; clear tSoil % convert to ug/m3 from ng/m3
            
            TEO(i,1) = (sum(TEO_ratio_XRF.*PM25_data_public(i,TEO_title_XRF),2,'omitnan'))/1000;
            
            method_metals(i) = 2;

        elseif ~isnan(PM25_data_public(i,15)) % no XRF, but there is ICP-MS data
            if PM25_data_public(i,21) > 0 % Fe is available in ICP-MS data
                Soil(i,1) =(10*(PM25_data_public(i,14) + PM25_data_public(i,21) + PM25_data_public(i,15)))/1000; % convert to ug/m3 from ng/m3
            elseif PM25_data_public(i,21) == 0 || isnan(PM25_data_public(i,21)) % no Fe in ICP-MS data
                Soil(i,1) = (10*(PM25_data_public(i,14) +  30*PM25_data_public(i,17) + PM25_data_public(i,15)))/1000 ;% convert to ug/m3 from ng/m3
            end
            
            TEO(i,1)=(sum(TEO_ratio_ICP.*PM25_data_public(i,TEO_title_ICP),2,'omitnan'))/1000; % divide by 1000 to convert to ug/m3 from ng/m3
        
        else % no ICP-MS or XRF data
            no_metals = no_metals+1; % counts number of filters that do not have ICP-MS of XRF data available
            Soil(i,1) = nan;
            TEO(i,1)= nan;
        end
    end
    
    if no_metals > 0
        fprintf('There are %d filters out of %d that do not have metal data \n',no_metals, size(PM25_data_public,1))
    end
    clear no_metals
    % ==== BC ====
%     BC = PM25_data_public(:,5); % non-hygroscopic
    BC = PM25_data_public(:,5)*0.06./0.1; % assume sigma 0.06 is too low. scale to sigma = 0.10

    % ==== OC ====
    % Need to convert to OM?
    OC_dry = PM25_data_public(:,end-1); % FTIR
    OC_wet = OC_dry.*MGrowth(3);
   
    % ==== RESIDUAL MASS (Organic Matter) ====
    merged_data_wet = [PM25_data_public(:,1) ASO4_wet ANO3_wet SSalt_wet Soil BC OC_wet TEO NaSO4_wet];
    colmask = any(isnan(merged_data_wet(:,[1 2 3 5])),2); % if any of [ASO4_wet ANO3_wet Soil] is nan, do not calculate RM
    
    RM_wet=merged_data_wet(:,1)-nansum(merged_data_wet(:,2:end),2);
    RM_wet(colmask == 1) = NaN; % if missing input inorganic, do not calculate RM
    RM_dry=RM_wet./MGrowth(3); % Dry residual matter
    
    Too_neg_OM=find(RM_wet./merged_data_wet(:,1)<-0.1); % set to NaN if negative by more than 10% of total wet (35% RH) PM2.5 mass
    RM_dry(Too_neg_OM,1)=NaN;
    RM_wet(Too_neg_OM,1)=NaN;
    clear too_neg_OM mask colmask index_OM merged_data_wet
    
    % ==== PBW ====
    % PBW = particle-bound water (mass conc., ug/m3) [currently not reported in RCFM files]
    PBW = sum( [ ASO4_wet -ASO4_dry ANO3_wet -ANO3_dry RM_wet -RM_dry OC_wet -OC_dry SSalt_wet -SSalt_dry NaSO4_wet -NaSO4_dry ] ,2,'omitnan');
    
    % ==== Kappa mix ====
    % dry volume
    Vol_mix = [ASO4_dry  ANO3_dry  RM_dry  OC_dry  SSalt_dry  Soil  BC  TEO  NaSO4_dry]./Density;
    
    kappa(4) = kappa(4)/2; % reducing kappa for NaCl by half to retain PBW at 0% RH
    
    for i=1:size(Vol_mix,1)
        kappa_mix_spec(i,:)=Vol_mix(i,:).*kappa./sum(Vol_mix(i,:),2,'omitnan');
    end
    kappa_mix=sum(kappa_mix_spec,2,'omitnan');
    kappa_mix(kappa_mix==0)=NaN;
    kappa_mix(kappa_mix>0.9)=NaN;
    kappa_mix(isnan(kappa_mix_spec(:,2))) = NaN; % if missing SO4 we are missing that and something else, therefore set to NaN
   
    species_plot = [PBW PM25_data_public(:,[2 4 3 6])  Soil TEO BC OC_dry RM_dry]; % order: {'water' 'SO4','NH4','NO3','SS','Dust','TEO','BC','OC','RM'};
       
    %% Additional Pie Charts: Cartridge specific pie charts
  
    infname =  sprintf('%s/Public_Data/Chemical_Filter_Data/Plots/Pie_spec_plots/ByCartridge/Cartridges_need_pie_charts.xlsx',direc);
    sheets = sheetnames(infname);
    if sum(contains(sheets,Site_codes{loc}))>0 % there is a sheet for this site
        cartlist = readtable(infname,'Sheet',Site_codes{loc});
        cartlist = table2array(cartlist);
        ind = find(ismember(PM25_cartridgeIDs,cartlist)==1);

        if ~isempty(ind)
            % sort mass type
            if mod(length(ind),8)==0
                fb = 7*(1:length(ind)/8); % the 7th filter is FP
                pm10 = 8*(1:length(ind)/8); % the 8th filter is FP
                ind([fb pm10]) = [];
            else
                error('Cartridge filter number is not 8.')
            end



            PM25_total = mean(PM25_data_public(ind,1),'omitnan');
            Cl =  mean(PM25_IC(ind,contains(IC_title,'IC_Cl')),'omitnan');
            K =  mean(PM25_IC(ind,contains(IC_title,'IC_K')),'omitnan');
            N = length(ind); % number of filters included

            Pie_made = PM25_RCFMavg_pie(PM25_total, species_plot(ind,:) , Site_cities,loc); % function for making pie chart
            if Pie_made == 1
                % add note of filter number
                notes = sprintf('Num of Filters = %d',N);
                text(-0.3, -0.1, notes,'units','normalized')

                fname = sprintf('%s/Public_Data/Chemical_Filter_Data/Plots/Pie_spec_plots/ByCartridge/%s_PM25_RCFMavg',direc,Site_codes{loc});
                saveas(gcf,sprintf('%s.png',fname))
                print(sprintf('%s.eps',fname),'-depsc')
                close all
            end

            % making another pie chart that includes Cl
            Pie_made = PM25_RCFMavg_pie_with_Cl(PM25_total, species_plot(ind,:) , Cl,K, Site_cities,loc); % function for making pie chart
            if Pie_made == 1
                % add note of filter number
                notes = sprintf('Num of Filters = %d',N);
                text(-0.1, -0.1, notes,'units','normalized')
                saveas(gcf,sprintf('%s_with_Cl.png',fname))
                print(sprintf('%s_with_Cl.eps',fname),'-depsc')
                close all


                % save data to a xlsx
                PM25 = PM25_data_public(ind,1);
                SO4 = PM25_data_public(ind,2);
                NH4 = PM25_data_public(ind,4);
                NO3 = PM25_data_public(ind,3);
                Na = PM25_data_public(ind,6);
                Cl =  PM25_IC(ind,contains(IC_title,'IC_Cl'));
                K  =  PM25_IC(ind,contains(IC_title,'IC_K'));
                savedata_dry(PM25_labels(ind,:), PM25, BC(ind,:), TEO(ind,:), Soil(ind,:), Na, NO3, NH4, SO4, PBW(ind,:), OC_dry(ind,:), RM_dry(ind,:),Cl,K,fname)
                savedata_wet(PM25_labels(ind,:), PM25, BC(ind,:), TEO(ind,:), Soil(ind,:), ANO3_wet(ind,:), ASO4_wet(ind,:), NaSO4_wet(ind,:), OC_wet(ind,:), RM_wet(ind,:),SSalt_wet(ind,:),fname)
                clear  PM25 SO4 NO3 NH4 Na Cl

            end
        end
    end

    %%
    close all
    clear *_dry *_wet BC *_index  colmask nan_idx  Titles  RFM_* Master_* SO4_NSS Soil SSalt_water start_date end_date TEO  Vol_mix 
    clear fileID i IN_tot_wet index_OM kappa_mix kappa_mix_spec PBW Too_neg_OM missing_PM25 PM*_index  PM25_dates Cart_ALL
    clear  species  PM*_data_public Master_CartridgeIDs h method_metals ic_mdl ic_mdl_*
    clear PM*_labels PM*_Barcodes PM*_cartridgeIDs PM*_LotIDs PM*_projectIDs PM*_masstype PM*_dates PM*_Mass PM*_Carbon PM*_Method PM*_flags PM*_ICP PM*_XRF PM*_IC
    
    fprintf('Finished processing data for %s \n\n', Site_cities{loc})
    
 end

%% Print out data summary tables
% print out the BC availability table
BC_table = table(Site_codes, HIPS_Mesured_BC, HIPS_Est_BC, HIPS_total_BC, SSR_total_BC);
writetable(BC_table,sprintf('%s/Analysis_Data/Black_Carbon/BC_availability.xlsx',direc));

%  print out reportable data num 
T = table(Site_cities, tSiteDataNum(:,1),tSiteDataNum(:,2),tSiteDataNum(:,3),tSiteDataNum(:,4),tSiteDataNum(:,5),...
                       tSiteDataNum(:,6),tSiteDataNum(:,7),tSiteDataNum(:,8),tSiteDataNum(:,9),tSiteDataNum(:,10),tSiteDataNum(:,11),'VariableNames',DataNum_tab_rows);
writetable(T,sprintf('%s/Analysis_Data/Reportable_Data_Count_New.xlsx',direc));


fprintf('END of Master File processing for %s \n', datestr(now))

diary off

%% FUNCTIONS
function savedata_wet(PM25_labels, PM25, BC, TEO, Soil, ANO3_wet, ASO4_wet, NaSO4_wet, OC_wet, RM_wet,SSalt_wet,fname)
        T = table(PM25_labels, PM25, BC, TEO, Soil, ANO3_wet, ASO4_wet, NaSO4_wet, OC_wet, RM_wet, SSalt_wet);
        tablename = sprintf('%s.xlsx',fname);
        writetable(T,tablename,'Sheet','wet')
end

function savedata_dry(PM25_labels, PM25, BC, TEO, Soil, Na, NO3, NH4, SO4, PBW, OC_dry, RM_dry,Cl,K,fname) 
        T = table(PM25_labels, PM25, BC, TEO, Soil, Na, NO3, NH4, SO4, PBW,Cl, K, OC_dry, RM_dry);
        tablename = sprintf('%s.xlsx',fname);
        delete(tablename) % delete existing file 
        writetable(T,tablename,'Sheet','dry')
end

function [tAnaMDL,tMDL,tUNC ]= findMDLUnc (tPM25_labels,ThisElement,Mass_conc, Vol, Ref_Dates,Ref_Values,FilterGroup)

    Titles = Ref_Values.Properties.VariableNames;
    Ref_Elements = table2array(Ref_Values(:,contains(Titles,'Element')));
%     Vol          = table2array(Ref_Values(1,contains(Titles,'Volume')));
    Unc_Pro_Vol  = 0.01*table2array(Ref_Values(1,contains(Titles,'Unc_proportion_volume')));

    for ii = 1:size(tPM25_labels,1)
        ThisFilter = tPM25_labels(ii,:);
        DateInd = zeros(size(FilterGroup,2)-1,1);

        for dd = 1:length(DateInd)
            DateInd(dd) = sum(contains(FilterGroup{2,dd},ThisFilter));
        end
        if sum(DateInd)>0
            DateInd(DateInd>1)=1; % raw data are duplicated in Archive
            MDL_Ind = find( contains(Ref_Elements,ThisElement) & Ref_Dates == FilterGroup{1,DateInd==1} );
            tAnaMDL(ii,1) = table2array(Ref_Values(MDL_Ind,contains(Titles,'Analytical_MDL'))) ./Vol(ii);
            tMDL(ii,1)    = table2array(Ref_Values(MDL_Ind,ismember(Titles,'MDL (ng)')))       ./Vol(ii);
            UncA          = table2array(Ref_Values(MDL_Ind,contains(Titles,'Unc_add')));
            UncP          = 0.01*table2array(Ref_Values(MDL_Ind,contains(Titles,'Unc_proportion_mass')));
            tUNC(ii,1) = sqrt( (UncA./Vol(ii))^2 + (UncP^2 + Unc_Pro_Vol^2)*Mass_conc(ii)^2 );

        else
            tAnaMDL(ii,1) = NaN;
            tMDL(ii,1) = NaN;
            tUNC(ii,1) = NaN;
        end
    end
    tAnaMDL = round(tAnaMDL,2);
    tMDL = round(tMDL,2);
    tUNC = round(tUNC,2);
end

function [PM25_labels,PM25_Barcodes,PM25_cartridgeIDs,PM25_LotIDs,PM25_projectIDs,...
        PM25_hours,PM25_masstype,PM25_dates,PM25_mass,PM25_IC,PM25_ICP,PM25_XRF,...
        PM25_carbon,PM25_Method,PM25_flags] = ...
        SortData (PM25_index, Master_IDs, Master_Barcodes,CartridgeIDs_master, LotIDs, ...
        projectIDs_master,Master_hours, Master_masstype, Master_dates,  ...
        Master_mass,  Master_IC,  Master_ICP,  Master_XRF, Master_carbon, ...
        Master_Method,Master_flags)

        PM25_labels   = Master_IDs(PM25_index,:);
        PM25_Barcodes = Master_Barcodes(PM25_index,:);
        PM25_cartridgeIDs = CartridgeIDs_master(PM25_index);
        PM25_LotIDs     = LotIDs(PM25_index,:);
        PM25_projectIDs = projectIDs_master(PM25_index);
        PM25_hours = Master_hours(PM25_index,:);
        
        PM25_masstype  = Master_masstype(PM25_index,:);
        PM25_dates  = Master_dates(PM25_index,:);
        PM25_mass  = Master_mass(PM25_index,:);
        PM25_IC  =   Master_IC(PM25_index,:);
        PM25_ICP  =  Master_ICP(PM25_index,:); 
        PM25_XRF  =  Master_XRF(PM25_index,:); 
        PM25_carbon  = Master_carbon(PM25_index,:);
        PM25_Method  = Master_Method(PM25_index,:);
        
        PM25_flags = Master_flags(PM25_index,:);
end

function [PM25_labels,PM25_Barcodes,PM25_cartridgeIDs,PM25_LotIDs,PM25_projectIDs,...
        PM25_hours,PM25_masstype,PM25_dates,PM25_mass,PM25_IC,PM25_ICP,PM25_XRF,...
        PM25_carbon,PM25_Method,PM25_flags] = ...
        RemoveRows (filter_index, ...
        PM25_labels,PM25_Barcodes,PM25_cartridgeIDs,PM25_LotIDs,PM25_projectIDs,...
        PM25_hours,PM25_masstype,PM25_dates,PM25_mass,PM25_IC,PM25_ICP,PM25_XRF,...
        PM25_carbon,PM25_Method,PM25_flags)
        
        PM25_hours = removerows(PM25_hours, filter_index);
        PM25_labels = removerows(PM25_labels, filter_index);
        PM25_cartridgeIDs = removerows(PM25_cartridgeIDs, filter_index);
        PM25_projectIDs = removerows(PM25_projectIDs, filter_index);
        PM25_Barcodes = removerows(PM25_Barcodes, filter_index);
        PM25_LotIDs   = removerows(PM25_LotIDs, filter_index);
        
        PM25_masstype  = removerows(PM25_masstype, filter_index);
        PM25_dates  = removerows(PM25_dates, filter_index);
        PM25_mass  = removerows(PM25_mass, filter_index);
        PM25_IC  =   removerows(PM25_IC, filter_index);
        PM25_ICP  =  removerows(PM25_ICP, filter_index);
        PM25_XRF  =  removerows(PM25_XRF, filter_index);
        PM25_carbon  = removerows(PM25_carbon, filter_index);
        PM25_Method  = removerows(PM25_Method, filter_index);
        
        PM25_flags = removerows(PM25_flags, filter_index);
        
end

function sumspec = SumSpec(i,XRF,ICP,IC,carbon,IC_col)

    carbon_sum = sum(carbon(i,[1 4]),2,'omitnan'); % Sum of HIPS BC and FTIR OC (if available)

    if ~isnan(XRF(i,2))  % check for presence of XRF data (Al here) first and sum only XRF metals
        % Metals + IC (without F, Cl, K, Mg, Ca) + BC
        spec = [sum(XRF(i,[2:3,5:end]),2,'omitnan')   sum(IC(i,IC_col),2,'omitnan')   carbon_sum];
        % currently don't use Na, S from XRF; ==> use IC_Na IC_SO4 already

    elseif ~isnan(ICP(i,3))  % check for presence of Al from ICP-MS and sum only ICP-MS metals
        spec = [sum(ICP(i,:),2,'omitnan')   sum(IC(i,IC_col),2,'omitnan')   carbon_sum];
     
    else % no Al from ICP or XRF, therefore sum xrf as default
        spec = [sum(XRF(i,[2:3,5:end]),2,'omitnan')   sum(IC(i,IC_col),2,'omitnan')   carbon_sum];
    
    end
    sumspec = sum(spec, 2,'omitnan');
end

function out = copyflag(cart_flags,Flag_parameters)
    
    publicflags = Flag_parameters.Flag;
    out = cart_flags;
    % remove flags
    for i = 1:length(cart_flags)
        out{i}='';
    end
    % add filter flags
    for i = 1:length(publicflags)
        flag_ind = find(contains(cart_flags,publicflags{i})==1);

        for j = 1:length(flag_ind)
            out{flag_ind(j)} = AddFlag( out{flag_ind(j)} , publicflags{i} );
        end

    end

end
