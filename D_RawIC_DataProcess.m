 %% %%%%%%%%%%%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: read and perform QA/QC on IC data straight from the IC
% system with no pre-processing; put IC data and relevant flags in the
% Master data file for each SPARTAN site that has data contained in the IC
% file; move IC files from "Raw" to relevant site folder

% Written by: Crystal Weagle
% Created: 25 May 2018

% EDITS: 
% Jan. 13 2020 CW - Revision 2 created on January 13 2020 to incorporate changes
% relevant to the new IC systems at Washington University in St. Louis.
% Updates from Revision 1 include new MDL's and addition of checking that
% sample concentrations are above MDL, if not they are set to zero. 

% Dec. 8 2021 Haihui Zhu - Revision 3 
% 1. Applied ReadMaster and WriteToMaster functions; 
% 2. Instead of reading data in unit of ug/L, this version read data in
%    uS min and calculate conc. in ug/L using the calibration data; 
%    ***Equation: signal = slope*conc+offset ==> conc = (signal-offset)/slope
% 3. ***This version will set all the F conc. into NaN due to an issue from
% the instrument. It should be corrected in the future. 
% 4. Applied the function 'AddFlag' for easier and more robust Flag adding to the master files.
% 5. Repalced 'xlsread' with 'readtable';
% 6. Added warning and error information when IC new data is uncompleted.
% 7. Skip deleting new file when (unlikely) new samples are not all assigned to master files. 

% March 18 2023 Haihui Zhu - Revision 4
% 1. Remove the 0 cap when the measurement is under MDL.
% 2. Added the AddFlagIC function for easy adding of flags.  

% Oct. 05 2023 Haihui Zhu - Revision 5
% 1.	Fixed error in NH3 calculation
% 2.	Applied calibration curves at low concentration for cation. 
% 3.	Added the feature to reprocessing all the WashU IC script (IC from Dal are rejected)
% 4.	Added the feature to collect all the IC signal area and save them to Ion_Chromatography/Area_Data

% Apr. 2 2024 Haihui Zhu - Revision 6
% Added calcualtion of the MDL for each filter and save them to
% 'Ion_Chromatography/MDL'
% NOTE:
%     The MDL are the limit for IC area (Us*min), it is non-negative
%     Beucase conc = (signal-offset)/slope, offset are mostly negative, the
%     MDL reported would be most likely positive, even MDL area is 0. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear ; clc
warning('off','backtrace')
addpath('UtilityFunctions')

%% Directories and User Switches

Redo_All_Archieve = 0; % set to 1 if want to re-process all archived data
Re_Collect_IC_areas = 0; % set to 1 if want to re-collect all area data. Note: only work when Redo_All_Archieve = 1

% Setup directories 
debug_mode = 0;
direc = find_root_dir(debug_mode);

direc_new = strcat(direc,'/Analysis_Data/Ion_Chromatography/NEW');
direc_archive = strcat(direc,'/Analysis_Data/Archived_Filter_Data/IC_data/');
direc_master = strcat(direc,'/Analysis_Data/Master_files');
direc_sampling = strcat(direc,'/Site_Sampling'); 
direc_area = strcat(direc,'/Analysis_Data/Ion_Chromatography/Area_Data'); % directory to save the IC signal area 
direc_mdl = strcat(direc,'/Analysis_Data/Ion_Chromatography/MDL'); % directory to save MDL data


% The diary function saves the processing history into an annual record
diary(sprintf('%s/Public_Data/Data_Processing_Records/IC/%s_IC_Record.txt',direc,datestr(now,'yyyy-mm-dd-HHMMSS')))
fprintf('%s \n', datestr(now))

Elog_filename = strcat(direc,'/Analysis_Data/E-Logs/IC E-Log.xlsx'); % contains extraction volumes for individual filters
[elog_status, Sheets_elog] = xlsfinfo(Elog_filename);

% low concentration calibration curve
IC_low_cat = load(strcat(direc,'/Analysis_Data/Ion_Chromatography/Syringe_filter_contamination/cation/IC_low_calibration.mat'),'Ca','K','Mg','Li','Na'); 
IC_low_an = load(strcat(direc,'/Analysis_Data/Ion_Chromatography/Syringe_filter_contamination/anion/IC_low_calibration_A.mat'),'Br','Cl','F','NO2','NO3','PO4','SO4'); 

% MDL reference file
mdl_in = strcat(direc_mdl,'/IC_LQL_Area_ToScript.xlsx');
mdl_ref = readtable(mdl_in);
mdl_ref_title = mdl_ref.Properties.VariableNames;
mdl_ref_dates = mdl_ref.Start_date;
mdl_ref = table2array(mdl_ref(:,2:end));

% file name to print filter info when IC value is 0 or nan
ic_zero_fname = 'IC_zero.txt';

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
CorrCoeff_threshold = 99.5;

%-------------   SITE DETAILS   --------------
site_details = readtable(strcat(direc_sampling,'/Site_details.xlsx'),'PreserveVariableNames',true);
Site_codes = table2array(site_details(:,1));
Site_cities = table2array(site_details(:,3));


% QC_conc = STD 1.25 ONLY TOUCH IF STD CHANGES
anions = {'F','Cl','NO2','Br','NO3','PO4','SO4'};
anion_QC_conc = [0.25 0.375 1.25 1.25 1.25 1.875 1.875];
cations = {'Li','Na','NH4','K','Mg','Ca'};
cation_QC_conc = [0.25 1.00 1.25 2.5 1.25 2.5];

% ---------- Minimum Detection Limits ------------ %% NEEDS TO BE UPDATED
% MDL recalulated based on IC runs in January 2020
anions_MDL = [0.005 0.005 0.004 0.01 0.003 0.035 0.006]; % ug/mL
% anions order = {'F','Cl','NO2','Br','NO3','PO4','SO4'};
cations_MDL = [0.0007 0.002 0.003 0.005 0.004 0.006]; % ug/mL
% cations order = {'Li','Na','NH4','K','Mg','Ca'};



%% %%%%%%%%%%%% Find files in either "NEW" dir or 'Archived' %%%%%%%%%%%%%%
if Redo_All_Archieve == 1

    % read the list of files that do not need to reprocess:
    fname = strcat(direc_archive,'No_need_to_reprocess.txt');
    fileID = fopen(fname, 'r');
    if fileID == -1
        error('File could not be opened for appending');
    end
    SkipList = {};
    while ~feof(fileID)
        line = fgetl(fileID);
        SkipList{end+1} = line;
    end
    fclose(fileID);

    % moving archive data to the new dir
    for loc = 1:length(Site_codes)

        ardir = strcat(direc_archive,Site_codes{loc});
        files = getFiles(ardir);

        % new file 
        new_files = getFiles(direc_new);

        for fid = 1:length(files)
            tfile = strcat(ardir,'/',files{fid,1});
            if any(contains(SkipList,files{fid,1})) % skip those in SkipList
                continue
            elseif any(contains(new_files,files{fid,1})) % skip copy one file multiple times - might make the script faster
                continue
            else
                [SUCCESS,MESSAGE,MESSAGEID] = copyfile(tfile, direc_new);
                if SUCCESS ~= 1
                    disp(MESSAGE)
                end
            end

        end

        % create new signal file
        if Re_Collect_IC_areas == 1
            % read filter infomation from master file
            master_file = sprintf('%s/%s_master.csv',direc_master,Site_codes{loc});

            [Titles,Filter_ID, Barcode, CartridgeID, LotID, projectID,Sample_hour, masstype, ...
                Master_dates, Master_mass, Master_IC, Master_ICP, Master_XRF,...
                Master_carbon, Mater_Nylon, Master_Method, Master_flags] = ReadMaster(master_file, Site_codes{loc});

            % delete if there is an exiting file & create a new one using 'read_and_update'
            area_file = sprintf('%s/%s_IC_Area.xlsx',direc_area,Site_codes{loc});
            if exist(area_file,'file')
                delete(area_file) % delte the existing one
            end
            area_table = read_and_update(area_file,Filter_ID,CartridgeID,Master_dates,masstype);
            
            % save area file with no area data
            writetable(area_table,area_file)

            clear Titles Filter_ID Barcode CartridgeID  LotID projectID Sample_hour  masstype ...
                Master_dates  Master_mass  Master_IC  Master_ICP  Master_XRF ...
                Master_carbon  Master_Method  Master_flags IC_area_* start_* stop_* newtable
            
        end
    end
    fprintf('\nDone copying archived IC raw file to NEW.\n\n')

    % creat a new IC_zero file
    delete(ic_zero_fname)
    fileID = fopen(ic_zero_fname, 'w');
    if fileID == -1
        error('File could not be opened for appending');
    end
    fprintf(fileID, 'Filter ID    date     Species   Value  \n');
    fclose(fileID);


end 

% now read files in the raw file directory
files = getFiles(direc_new);
 
%% %%%%%%%%%%%% Processing each file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for IC_file = 1:length(files)
    
    % Do not read hidden files
    char_file = char(files(IC_file));
    if char_file(1)=="."
        continue
    end
    clear char_file
    
    filename = sprintf('%s/%s',direc_new,files{IC_file,1});
    fprintf('\nReading %s \n', files{IC_file,1});
    
    Sheets = sheetnames(filename);
    
    for i = 1:length(Sheets)
        if contains(Sheets(i),'Summary') 
            SheetName = char(Sheets(i));
        elseif contains(Sheets(i),'Calibration')
            cal_SheetName = char(Sheets(i));
        end
    end
    % read Calibration data
    [raw_cal,CubicIndicator] = ReadRawCal(filename,cal_SheetName);
    % Do not reprocess IC from Dal
    [I , J] = findIndexContaining('Operator:', raw_cal); 
    Operator = raw_cal(I,J+1);
    if ~matches(Operator,'SPARTAN') % a robust check
        fprintf('Skip this file from Dalhousie, Operator: %s\n',Operator )
        
        % add this file to the 'No_need_to_reprocess.txt'
        fname = strcat(direc_archive,'No_need_to_reprocess.txt');
        fileID = fopen(fname, 'a');
        if fileID == -1
            error('File could not be opened for appending');
        end
        fprintf(fileID, '%s\n',files{IC_file,1});
        fclose(fileID);
        
        % remove that from the 'NEW' folder
        delete(filename)
        continue
    end
    % find the date of operation
    [I , J] = findIndexContaining('Date', raw_cal); 
    oprt_date = raw_cal(I,J+1); % should be class of datetime
    dateind = find(oprt_date-mdl_ref_dates>0); % find the mdl_ref that measured before this IC data
    dateind = dateind(end); % the last one is the most recent mdl_ref

    % Read IC raw data
    [data_pre,labels_all,data_all,data_type]= ReadRawIC(filename,SheetName,files{IC_file,1});
    % data_type:
    % 1 = anion data
    % 2 = cation data
    
    clear cal_SheetName SheetName Sheets i 
    test = 0; % zero means the file does not have any test data contained, there is a loop below that changes this to 1 if test data is detected and will save the file in a designated folder
   
    % create an empty flag matrix for Calibration conditions 
    calib_flags = '';
 
    
    %% %%%%%%%%%%%%%%%%%%% BEGIN ANION SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if data_type == 1 % anion data file
        
        % Columns are:
        % 1 = Fluoride
        % 2 = Chloride
        % 3 = Nitrite
        % 4 = Bromide
        % 5 = Nitrate
        % 6 = Phosphate
        % 7 = Sulfate
        
        %% ------ get correlation coefficients from calibration curve and flag if below 0.98-------
        [~ , J] = findIndexContaining('Coeff.', raw_cal); 
        % this section finds the row where each correlation coefficient is found. 
        k = 19; 
        [F,Fj]    =  FindIdx(raw_cal,k,'Fluoride'); clear Fj % column vector not needed)
        [Cl,Clj]  =  FindIdx(raw_cal,k, 'Chloride'); clear Clj
        [NO2,NO2j]=  FindIdx(raw_cal,k,'Nitrite'); clear NO2j
        [Br,Brj]  =  FindIdx(raw_cal,k,'Bromide'); clear Brj
        [NO3,NO3j] =  FindIdx(raw_cal,k, 'Nitrate'); clear NO3j
        [PO4,PO4j] =  FindIdx(raw_cal,k,'Phosphate'); clear PO4j
        [SO4,SO4j] =  FindIdx(raw_cal,k,'Sulfate'); clear SO4j
        clear k
        
        disp('Starting quality assurance section')
        disp('Checking calibration data')
        
        % calculate ion conc. in ug/L 
        if sum(contains(raw_cal,'HMS'),'all') == 0 
            fprintf('No HMS calibration for this file.\n')
            data_post(:,1) = getpostdata(data_pre(:,1),raw_cal(F,:),  CubicIndicator, data_all(:,1),IC_low_an,'F');
            data_post(:,2) = getpostdata(data_pre(:,2),raw_cal(Cl,:), CubicIndicator, data_all(:,2),IC_low_an,'Cl');
            data_post(:,3) = getpostdata(data_pre(:,3),raw_cal(NO2,:),CubicIndicator, data_all(:,3),IC_low_an,'NO2');
            data_post(:,4) = getpostdata(data_pre(:,4),raw_cal(Br,:), CubicIndicator, data_all(:,4),IC_low_an,'Br');
            data_post(:,5) = getpostdata(data_pre(:,5),raw_cal(NO3,:),CubicIndicator, data_all(:,5),IC_low_an,'NO3');
            data_post(:,6) = getpostdata(data_pre(:,6),raw_cal(PO4,:),CubicIndicator, data_all(:,6),IC_low_an,'PO4');
            data_post(:,7) = getpostdata(data_pre(:,7),raw_cal(SO4,:),CubicIndicator, data_all(:,7),IC_low_an,'SO4');
        else 
            ind = find(contains(raw_cal,'HMS')); 
            if str2num(raw_cal(ind(1),4)) == 0 % contain the string 'HMS' but the calibration is nan
                fprintf('No HMS calibration for this file.\n')
                data_post(:,1) = getpostdata(data_pre(:,1),raw_cal(F,:),  CubicIndicator, data_all(:,1),IC_low_an,'F');
                data_post(:,2) = getpostdata(data_pre(:,2),raw_cal(Cl,:), CubicIndicator, data_all(:,2),IC_low_an,'Cl');
                data_post(:,3) = getpostdata(data_pre(:,3),raw_cal(NO2,:),CubicIndicator, data_all(:,3),IC_low_an,'NO2');
                data_post(:,4) = getpostdata(data_pre(:,4),raw_cal(Br,:), CubicIndicator, data_all(:,4),IC_low_an,'Br');
                data_post(:,5) = getpostdata(data_pre(:,5),raw_cal(NO3,:),CubicIndicator, data_all(:,5),IC_low_an,'NO3');
                data_post(:,6) = getpostdata(data_pre(:,6),raw_cal(PO4,:),CubicIndicator, data_all(:,6),IC_low_an,'PO4');
                data_post(:,7) = getpostdata(data_pre(:,7),raw_cal(SO4,:),CubicIndicator, data_all(:,7),IC_low_an,'SO4');
            else
                fprintf('HMS included in calibration, no low conc curve applied.\n')
                data_post(:,1) = getpostdata(data_pre(:,1),raw_cal(F,:),  CubicIndicator, data_all(:,1),IC_low_an,'N/A');
                data_post(:,2) = getpostdata(data_pre(:,2),raw_cal(Cl,:), CubicIndicator, data_all(:,2),IC_low_an,'N/A');
                data_post(:,3) = getpostdata(data_pre(:,3),raw_cal(NO2,:),CubicIndicator, data_all(:,3),IC_low_an,'N/A');
                data_post(:,4) = getpostdata(data_pre(:,4),raw_cal(Br,:), CubicIndicator, data_all(:,4),IC_low_an,'N/A');
                data_post(:,5) = getpostdata(data_pre(:,5),raw_cal(NO3,:),CubicIndicator, data_all(:,5),IC_low_an,'N/A');
                data_post(:,6) = getpostdata(data_pre(:,6),raw_cal(PO4,:),CubicIndicator, data_all(:,6),IC_low_an,'N/A');
                data_post(:,7) = getpostdata(data_pre(:,7),raw_cal(SO4,:),CubicIndicator, data_all(:,7),IC_low_an,'N/A');
            end
            clear ind
        end
    
        % calculate MDL
           mdl_calc = nan(1,size(data_pre,2));
        if sum(contains(raw_cal,'HMS'),'all') == 0
            fprintf('No HMS calibration for this file.\n')
            mdl_calc(:,1) = getpostdata(mdl_ref(dateind,1),raw_cal(F,:),  CubicIndicator, NaN, IC_low_an,'F');
            mdl_calc(:,2) = getpostdata(mdl_ref(dateind,2),raw_cal(Cl,:), CubicIndicator, NaN, IC_low_an,'Cl');
            mdl_calc(:,3) = getpostdata(mdl_ref(dateind,3),raw_cal(NO2,:),CubicIndicator, NaN, IC_low_an,'NO2');
            mdl_calc(:,4) = getpostdata(mdl_ref(dateind,4),raw_cal(Br,:), CubicIndicator, NaN, IC_low_an,'Br');
            mdl_calc(:,5) = getpostdata(mdl_ref(dateind,5),raw_cal(NO3,:),CubicIndicator, NaN, IC_low_an,'NO3');
            mdl_calc(:,6) = getpostdata(mdl_ref(dateind,6),raw_cal(PO4,:),CubicIndicator, NaN, IC_low_an,'PO4');
            mdl_calc(:,7) = getpostdata(mdl_ref(dateind,7),raw_cal(SO4,:),CubicIndicator, NaN, IC_low_an,'SO4');
        else
            ind = find(contains(raw_cal,'HMS'));
            if str2num(raw_cal(ind(1),4)) == 0 % contain the string 'HMS' but the calibration is nan
                mdl_calc(:,1) = getpostdata(mdl_ref(dateind,1),raw_cal(F,:),  CubicIndicator, NaN, IC_low_an,'F');
                mdl_calc(:,2) = getpostdata(mdl_ref(dateind,2),raw_cal(Cl,:), CubicIndicator, NaN, IC_low_an,'Cl');
                mdl_calc(:,3) = getpostdata(mdl_ref(dateind,3),raw_cal(NO2,:),CubicIndicator, NaN, IC_low_an,'NO2');
                mdl_calc(:,4) = getpostdata(mdl_ref(dateind,4),raw_cal(Br,:), CubicIndicator, NaN, IC_low_an,'Br');
                mdl_calc(:,5) = getpostdata(mdl_ref(dateind,5),raw_cal(NO3,:),CubicIndicator, NaN, IC_low_an,'NO3');
                mdl_calc(:,6) = getpostdata(mdl_ref(dateind,6),raw_cal(PO4,:),CubicIndicator, NaN, IC_low_an,'PO4');
                mdl_calc(:,7) = getpostdata(mdl_ref(dateind,7),raw_cal(SO4,:),CubicIndicator, NaN, IC_low_an,'SO4');
            else
                mdl_calc(:,1) = getpostdata(mdl_ref(dateind,1),raw_cal(F,:),  CubicIndicator, NaN, IC_low_an,'N/A');
                mdl_calc(:,2) = getpostdata(mdl_ref(dateind,2),raw_cal(Cl,:), CubicIndicator, NaN, IC_low_an,'N/A');
                mdl_calc(:,3) = getpostdata(mdl_ref(dateind,3),raw_cal(NO2,:),CubicIndicator, NaN, IC_low_an,'N/A');
                mdl_calc(:,4) = getpostdata(mdl_ref(dateind,4),raw_cal(Br,:), CubicIndicator, NaN, IC_low_an,'N/A');
                mdl_calc(:,5) = getpostdata(mdl_ref(dateind,5),raw_cal(NO3,:),CubicIndicator, NaN, IC_low_an,'N/A');
                mdl_calc(:,6) = getpostdata(mdl_ref(dateind,6),raw_cal(PO4,:),CubicIndicator, NaN, IC_low_an,'N/A');
                mdl_calc(:,7) = getpostdata(mdl_ref(dateind,7),raw_cal(SO4,:),CubicIndicator, NaN, IC_low_an,'N/A');
            end
            clear ind
        end

        % Set F concentration into NaN before the issue in IC is fixed
        data_post(:,1) = NaN; 
        % Set F MLD to NaN too. 
        mdl_calc(:,1) = NaN;
        
        % Replace data_all by data_post
        data_all = data_post;
        
        % Check Fluoride Calibration
        if ~isempty(F) 
            if str2double(raw_cal(F,J)) < CorrCoeff_threshold
                disp('WARNING: STD correlation for Fluoride below 0.995')
                calib_flags = AddFlagIC(calib_flags,'S-F');
            end
        else 
            calib_flags = AddFlagIC(calib_flags,'Sx-F');
            disp('ERROR: no correlation coefficient for Fluoride - please re-export data from IC computer')
        end
        
        % Check Chloride Calibration 
        if ~isempty(Cl) 
            if str2double(raw_cal(Cl,J)) < CorrCoeff_threshold
                disp('WARNING: STD correlation for Fluoride below 0.995')
                calib_flags = AddFlagIC(calib_flags,'S-Cl');
            end
        else
            calib_flags = AddFlagIC(calib_flags,'Sx-Cl');
            disp('ERROR: no correlation coefficient for Chloride - please re-export data from IC computer')
        end
        
        % Check Nitrite Calibration
        if ~isempty(NO2) 
            if str2double(raw_cal(NO2,J)) < CorrCoeff_threshold
                disp('WARNING: STD correlation for Nitrite below 0.995')
                calib_flags = AddFlagIC(calib_flags,'S-NO2');
            end
        else
            calib_flags = AddFlagIC(calib_flags,'Sx-NO2');
            disp('ERROR: no correlation coefficient for Nitrite - please re-export data from IC computer')
        end
        
        % Check Bromide Calibration 
        if ~isempty(Br) 
            if str2double(raw_cal(Br,J)) < CorrCoeff_threshold
                disp('WARNING: STD correlation for Bromide below 0.995')
                calib_flags = AddFlagIC(calib_flags,'S-Br');
            end
        else
            disp('ERROR: no correlation coefficient for Bromide - please re-export data from IC computer')
            calib_flags = AddFlagIC(calib_flags,'Sx-Br');
        end

        % Check Nitrate Calibration 
        if ~isempty(NO3)
            if str2double(raw_cal(NO3,J)) < CorrCoeff_threshold
                disp('WARNING: STD correlation for Nitrate below 0.995')
                calib_flags = AddFlagIC(calib_flags,'S-NO3');
            end
        else
            disp('ERROR: no correlation coefficient for Nitrate - please re-export data from IC computer')
            calib_flags = AddFlagIC(calib_flags,'Sx-NO3');
        end
        
        % Check Phosphate Calibration 
        if ~isempty(PO4) 
            if str2double(raw_cal(PO4,J)) < CorrCoeff_threshold
                disp('WARNING: STD correlation for Phosphate below 0.995')
                calib_flags = AddFlagIC(calib_flags,'S-PO4');
            end
        else
            disp('ERROR: no correlation coefficient for Phosphate - please re-export data from IC computer')
            calib_flags = AddFlagIC(calib_flags,'Sx-PO4');
        end

        % Check Sulfate Calibration 
        if ~isempty(SO4) 
            if str2double(raw_cal(SO4,J)) < CorrCoeff_threshold
                disp('WARNING: STD correlation for Sulfate below 0.995')
                calib_flags = AddFlagIC(calib_flags,'S-SO4');
            end
        else 
            disp('ERROR: no correlation coefficient for Sulfate - please re-export data from IC computer')
            calib_flags = AddFlagIC(calib_flags,'Sx-SO4');
        end
       
        %--------- index different data types -------------
        H2O_ind = find(contains(labels_all,'H2O'));
        Samples_ind = find(contains(labels_all,'-')); % containing blanks, will be removed later
        STD_ind = find(contains(labels_all,'STD'));
        LB_ind = find(contains(labels_all,'-LB')); % index for lab blanks (started this Fall 2017)
        if isempty(LB_ind) == 1
            LB_ind = find(contains(labels_all,{'LAB','BLANK'}));
        end
        
        if ~isempty(LB_ind) % checks if lab blanks in run - no lab blanks prior to Fall 2017
            for i = 1:length(LB_ind)
                ind(i) = find(Samples_ind == LB_ind(i));
            end
            Samples_ind = removerows(Samples_ind, ind); % removes lab blanks from index to samples
        end

       % separate filters end with 'N' as they are not for master files
        if_nylon = zeros(size(Samples_ind));
        for i = 1:length(Samples_ind)
            tfilter = labels_all{Samples_ind(i)};
            if tfilter(end) == 'N'
                if_nylon(i) = 1;
            end
        end
        ind = find(if_nylon==1);
        nylon_ind = Samples_ind(if_nylon==1);
        Samples_ind = removerows(Samples_ind, ind);  clear if_nylon i ind


        % removes double dashed samples from index to samples
        num_dash = zeros(size(Samples_ind));
        for i = 1:length(Samples_ind)
            num_dash(i) = count(labels_all(Samples_ind(i)),"-"); 
        end
        extra_dash = find(num_dash > 1);
        Samples_ind = removerows(Samples_ind, extra_dash); clear extra_dash num_dash ind
        
        % ----Get data for different data type (standards, samples, lab blanks, and waters) ------
        data_STD = data_all(STD_ind,:);
        data_H2O = data_all(H2O_ind(2:end),:); % first water is often not representative of run, therefore start at second water
        data_teflon = data_all(Samples_ind,:);
        data_nylon = data_all(nylon_ind,:);
%         data_samples(isnan(data_samples)) = 0; % if there is no ion concentration found, set to zero rather than NaN in file. NaN indicates a missing number
        data_LB = data_all(LB_ind,:);
        area_samples = data_pre(Samples_ind,:);

        labels_stds = labels_all(STD_ind);
        labels_teflon = labels_all(Samples_ind);
        labels_nylon = labels_all(nylon_ind);
        labels_nylon2 = labels_nylon;
        for i = 1:length(labels_nylon)
            tlabel = char(labels_nylon(i));
            tlabel = tlabel(1:end-1); % remove 'N'
            labels_nylon(i) = tlabel;
        end

        % identify sites 
        labels_char = char(labels_all([Samples_ind; nylon_ind]));
        site_IDs = unique(labels_char(:,1:4),'rows');
        labels_LB = labels_all(LB_ind);
        
        clear *_ind num_cal raw_cal txt_cal labels_teflon_char
        MDL_flags = cell(size(data_teflon,1),1);

        
        %-------- This section check if the QC STD 1.25 concentration is within +/-10 % of the specified concentration for each anion -------
        % (data before May 2018 will not have these)
        QC_STD_ind = find(contains(labels_stds,'1.25'));
        
        if ~isempty(QC_STD_ind)
            data_STD_QC = data_STD(QC_STD_ind, :);
            QC_avgpercentdiff = 100*((abs(mean(data_STD_QC,1) - anion_QC_conc))./anion_QC_conc);
            for i = 1:length(QC_avgpercentdiff)
                if QC_avgpercentdiff(i) > 10
                    fprintf('WARNING: Avg QC STD 1.25 for %s exceeds +/- 10 %% of specified concentration \n', anions{i})
                    calib_flags = AddFlagIC(calib_flags,sprintf('qc1-%s',anions{i}));
                end
            end
            clear i 
        end
        
        %-------- This section check if average anion concentrations in waters is < 10 * MDL  -------
        H2O_avg = mean(data_H2O,1,'omitnan');
        for i = 1:size(H2O_avg,2)
            if H2O_avg(i) > 10*anions_MDL(i)
                fprintf('WARNING: Average H2O for %s is > 10*MDL \n', anions{i})
                calib_flags = AddFlagIC(calib_flags,sprintf('qc2-%s',anions{i}));
            end
        end
        clear i 
        
        % combine flags for writing to data file
        if size(calib_flags,1) > 0
            for t = 1:length(data_teflon)
                MDL_flags{t} = AddFlagIC(MDL_flags{t},calib_flags);
            end
            fprintf('Anion flags for this file are: %s \n',calib_flags)
        else
            disp('No flags for this file')
        end
        clear calib_flags
        
        disp('All quality assurance for this file has been completed')
        
        %% -------- This section checks for existance of a master data file for each site in the IC data file ----------
        disp('Looking for master files')
        
        NN = 0; % mark any unassigned IC sample (i.e. not found in master file) 
        
        for i = 1:size(site_IDs,1)
            
            master_file = sprintf('%s/%s_master.csv',direc_master,site_IDs(i,:));
            
           [Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs, projectIDs_master,Master_hours, Master_masstype, ...
            Master_dates, Master_mass, Master_IC, Master_ICP, Master_XRF,...
            Master_carbon, Master_Nylon, Master_Method, Master_flags] = ReadMaster(master_file,site_IDs(i,:));

           if ~isempty(Master_IDs)

               % read area file 
               area_file = sprintf('%s/%s_IC_Area.xlsx',direc_area,site_IDs(i,:));
               area_table = read_and_update(area_file,Master_IDs,CartridgeIDs_master,Master_dates,Master_masstype);
               % read mdl file
               mdl_file = sprintf('%s/%s_IC_MDL.xlsx',direc_mdl,site_IDs(i,:));
               mdl_table = read_and_update(mdl_file,Master_IDs,CartridgeIDs_master,Master_dates,Master_masstype);

               samples_ICfile_ind = find(contains(labels_teflon,site_IDs(i,:))); % index of samples for site_ID(i) in IC file
                
               % find the index of teflon filters in the master file
                masterID_teflon_ind = [];
                for k = 1:length(samples_ICfile_ind)
                    if ~isempty( find(matches(Master_IDs, labels_teflon(samples_ICfile_ind(k))),1) ) 
                        % filter ID found in Master_IDs
                        masterID_teflon_ind(k) = find(matches(Master_IDs, labels_teflon(samples_ICfile_ind(k)))); % finds the row index in the master file for samples in IC file
                        teflon_IC_Master(k) = samples_ICfile_ind(k);
                    else
                        % filter ID not found: could be filter ID format
                        % mismatch, try fix it
                        tfilterid = char(labels_teflon(samples_ICfile_ind(k))) ;
                        if length(tfilterid) == 8
                            tfilterid2 = strcat(tfilterid(1:5),'0',tfilterid(6:8));
                            masterID_teflon_ind(k) = find(matches(Master_IDs, tfilterid2)); % finds the row index in the master file for samples in IC file
                            teflon_IC_Master(k) = samples_ICfile_ind(k);
                        else
                            tfilter = labels_teflon(samples_ICfile_ind(k));

                            warning( '%s not found in Master File\n',tfilter)
                            NN = NN +1;

                        end
                    end
                end
         
                % find the index of nylon filters in the master file
                samples_ICfile_ind = find(contains(labels_nylon,site_IDs(i,:))); % index of samples for site_ID(i) in IC file
                
                 masterID_nylon_ind = [];
                for k = 1:length(samples_ICfile_ind)
                    if ~isempty( find(matches(Master_IDs, labels_nylon(samples_ICfile_ind(k))),1) ) 
                        % filter ID found in Master_IDs
                        masterID_nylon_ind(k) = find(matches(Master_IDs, labels_nylon(samples_ICfile_ind(k)))); % finds the row index in the master file for samples in IC file
                        nylon_IC_Master(k) = samples_ICfile_ind(k);
                    else
                        % filter ID not found: could be filter ID format
                        % mismatch, try fix it
                        tfilterid = char(labels_nylon(samples_ICfile_ind(k))) ;
                        if length(tfilterid) == 8
                            tfilterid2 = strcat(tfilterid(1:5),'0',tfilterid(6:8));
                            masterID_nylon_ind(k) = find(matches(Master_IDs, tfilterid2)); % finds the row index in the master file for samples in IC file
                            nylon_IC_Master(k) = samples_ICfile_ind(k);
                        else
                            tfilter = labels_nylon(samples_ICfile_ind(k));

                            warning( '%s not found in Master File\n',tfilter)
                            NN = NN +1;

                        end
                    end
                end

                % get extracton volume
                elog_index = find(contains(Sheets_elog,site_IDs(i,:)) == 1);
                Elog_SheetName = char(Sheets_elog(elog_index));
                cell_elog=readcell(Elog_filename,'Sheet',Elog_SheetName);

                % nylon data
                if ~isempty(masterID_nylon_ind)
                    % nan out -xxxxx before writing to master files
                    temp = data_nylon(nylon_IC_Master,:);
                    temp(temp<-100) = nan;
                    data_nylon(nylon_IC_Master,:) = temp;

                    [extraction_volumes_N,extraction_dates_N] = GetVolume(cell_elog,labels_nylon2(samples_ICfile_ind));
                    Master_Nylon(masterID_nylon_ind,1:7) = data_nylon(nylon_IC_Master,:).*extraction_volumes_N;

                    % write to master file
                    WriteToMaster ( Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs, projectIDs_master,...
                                    Master_hours, Master_masstype, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
                                    Master_carbon,Master_Nylon, Master_Method, Master_flags,...
                                    direc_master,site_IDs(i,:))
                end

                % proceed to write data to master_IC and master_nylon
                if ~isempty(masterID_teflon_ind) 

                    [extraction_volumes,extraction_dates] = GetVolume(cell_elog,Master_IDs(masterID_teflon_ind));
                    
                    % find any zeros and print to 'IC_zero.txt'
                    find_ic_zeros(data_teflon(teflon_IC_Master,:), Master_IDs(masterID_teflon_ind), ...
                        extraction_dates, Master_masstype(masterID_teflon_ind), ic_zero_fname, data_type)
                    % nan out -xxxxx before writing to master files
                    temp = data_teflon(teflon_IC_Master,:);
                    temp(temp<-100) = nan;
                    data_teflon(teflon_IC_Master,:) = temp;


                    % adding IC data to Master_IC
                    Master_IC(masterID_teflon_ind,1:7) = data_teflon(teflon_IC_Master,:).*extraction_volumes;


                    % adding area to area file
                    for ii = 1:7
                        if size(area_table(masterID_teflon_ind,5+ii)) == size(area_samples(teflon_IC_Master,ii))
                            area_table(masterID_teflon_ind,5+ii) = table(area_samples(teflon_IC_Master,ii));
                        else
                            error('size of area data do not match')
                        end
                    end

                    % adding mdl to mdl file
                    for ii = 1:7
                        mdl_table(masterID_teflon_ind,5+ii) = table(mdl_calc(1,ii).*extraction_volumes);
                    end

                    % adding flag to master_initial_flag
                    if exist('MDL_flags','var')
                    for k = 1:numel(masterID_teflon_ind)
                        Master_flags{masterID_teflon_ind(k)} = AddFlag( Master_flags{masterID_teflon_ind(k)} , MDL_flags{teflon_IC_Master(k)} );
                    end
                    end


                    % --- finish this site ---
                    % write to master file
                    WriteToMaster ( Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs, projectIDs_master,...
                                    Master_hours, Master_masstype, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
                                    Master_carbon,Master_Nylon, Master_Method, Master_flags,...
                                    direc_master,site_IDs(i,:))
                    % write to area file
                    writetable(area_table,area_file)
                    % write to mdl file
                    writetable(mdl_table,mdl_file)
                else
                    fprintf('Samples in IC file not found in %s Master file \n', site_IDs(i,:));
                    disp('Moving on to next site')
                end
            end  
            
            clear Master_data Master_data_initial txt_elog raw_elog num_elog extraction_volumes LotIDs_master Master_IC
            clear Master_IDs masterID_teflon_ind masterID_nylon_ind samples_ICfile_ind teflon_IC_Master nylon_IC_Master CartridgeIDs_master LotIDs
            clear Master_flags Master_flags_initial master_file hours_sampled Master_hours projectIDs_master Master_Barcodes 
        end
        
        
        %% copy and archive IC raw files
        for i = 1:size(site_IDs,1)
            file_destination = strcat(direc_archive,sprintf('%s/',site_IDs(i,:)));
            status(i) = CopyFile(filename,file_destination);
            clear file_destination
        end
        clear i
        
        if test == 1
            file_destination = strcat(direc_archive,'Test_data/');
            status = CopyFile(filename,file_destination);
            clear file_destination
            disp('Test data detected, file has been saved to test data folder')
        end
        
        if isempty(find(status == 0) == 1) && NN == 0 
            delete(filename)
            disp('File has been deleted from raw IC data folder')
        elseif test == 1
            delete(filename)
        else
            if NN == 0
                disp('Cannot delete file from raw IC data folder as it has not been moved to all site folders')
            else
                disp('Cannot delete file from raw IC data folder as some samples have not been assigned to master files')
            end
        end
        clear status 
        
    end
    
    %%%%%%%%%%%%%%%%%%%% END OF ANION SECTION OF SORTING AND QC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% %%%%%%%%%%%%%%%%%%% BEGIN CATION SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if data_type == 2 % cation data file
        
        % columns are:
        % 1 = Lithium
        % 2 = Sodium
        % 3 = Ammonium
        % 4 = Potassium
        % 5 = Magnesium
        % 6 = Calcium
        
        % ------ get correlation coefficients from calibration curve and flag if below 0.98-------
        [~ , J] = findIndexContaining('Coeff.', raw_cal); 
        % this section finds the row where each correlation coefficient is found. 
        k = 19; 
        [Li,~] =  FindIdx(raw_cal,k,'Lithium');  % column vector not needed)
        [Na,~] =  FindIdx(raw_cal,k,'Sodium');
        [NH4,~]=  FindIdx(raw_cal,k,'Ammonium'); 
        [K,~]  =  FindIdx(raw_cal,k,'Potassium'); 
        [Mg,~] =  FindIdx(raw_cal,k,'Magnesium'); 
        [Ca,~] =  FindIdx(raw_cal,k,'Calcium');   clear k
        disp('Starting quality assurance section')
        disp('Checking calibration data')
        
        % calculate ion conc. in ug/L 
        data_post(:,1) = getpostdata(data_pre(:,1),raw_cal(Li,:), CubicIndicator, data_all(:,1),IC_low_cat,'Li'); 
        data_post(:,2) = getpostdata(data_pre(:,2),raw_cal(Na,:), CubicIndicator, data_all(:,2),IC_low_cat,'Na'); 
        data_post(:,3) = getpostdata(data_pre(:,3),raw_cal(NH4,:),CubicIndicator, data_all(:,3),IC_low_cat,'NH4'); 
        data_post(:,4) = getpostdata(data_pre(:,4),raw_cal(K,:), CubicIndicator, data_all(:,4),IC_low_cat,'K'); 
        data_post(:,5) = getpostdata(data_pre(:,5),raw_cal(Mg,:),CubicIndicator, data_all(:,5),IC_low_cat,'Mg'); 
        data_post(:,6) = getpostdata(data_pre(:,6),raw_cal(Ca,:),CubicIndicator, data_all(:,6),IC_low_cat,'Ca'); 
       
        % calculate mdl
        mdl_calc = nan(1,size(data_post,2));
        mdl_calc(:,1) = getpostdata(mdl_ref(dateind,1),raw_cal(Li,:), CubicIndicator, NaN,IC_low_cat,'Li'); 
        mdl_calc(:,2) = getpostdata(mdl_ref(dateind,2),raw_cal(Na,:), CubicIndicator, NaN,IC_low_cat,'Na'); 
        mdl_calc(:,3) = getpostdata(mdl_ref(dateind,3),raw_cal(NH4,:),CubicIndicator, NaN,IC_low_cat,'NH4'); 
        mdl_calc(:,4) = getpostdata(mdl_ref(dateind,4),raw_cal(K,:), CubicIndicator, NaN,IC_low_cat,'K'); 
        mdl_calc(:,5) = getpostdata(mdl_ref(dateind,5),raw_cal(Mg,:),CubicIndicator, NaN,IC_low_cat,'Mg'); 
        mdl_calc(:,6) = getpostdata(mdl_ref(dateind,6),raw_cal(Ca,:),CubicIndicator, NaN,IC_low_cat,'Ca'); 
       
        % Replace data_all by data_post
        data_all = data_post;
        
        % Check Lithium Calibration
        if isempty(Li) == 0
            if str2double(raw_cal(Li,J)) < CorrCoeff_threshold
                disp('WARNING: STD correlation for Lithium below 0.995')
                calib_flags = AddFlagIC(calib_flags,'S-Li');
            end
        elseif isempty(Li) == 1 
            disp('ERROR: no correlation coefficient for Lithium - please re-export IC data')
            calib_flags = AddFlagIC(calib_flags,'Sx-Li');
        end
        
        % Check Sodium Calibration
        if isempty(Na) == 0 
            if str2double(raw_cal(Na,J)) < CorrCoeff_threshold
                disp('WARNING: STD correlation for Sodium below 0.995')
                calib_flags = AddFlagIC(calib_flags,'S-Na');
            end
        elseif isempty(Na) == 1 
            disp('ERROR: no correlation coefficient for sodium - please re-export IC data')
            calib_flags = AddFlagIC(calib_flags,'Sx-Na');
        end
        
        % Check Ammonium Calibration
        if isempty(NH4) == 0 
            if str2double(raw_cal(NH4,J)) < CorrCoeff_threshold
                disp('WARNING: STD correlation for ammonium below 0.995')
                calib_flags = AddFlagIC(calib_flags,'S-NH4');
            end
        elseif isempty(NH4) == 1 
            disp('ERROR: no correlation coefficient for ammonium - please re-export data from IC computer')
            calib_flags = AddFlagIC(calib_flags,'Sx-NH4');
        end
        
        % Check Potassium Calibration
        if isempty(K) == 0 
            if str2double(raw_cal(K,J)) < CorrCoeff_threshold
                disp('WARNING: STD correlation for Potassium below 0.995')
                calib_flags = AddFlagIC(calib_flags,'S-K');
            end
        elseif isempty(K) == 1 
            disp('ERROR: no correlation coefficient for Potassium - please re-export data from IC computer')
            calib_flags = AddFlagIC(calib_flags,'Sx-K');
        end
        
        % Check Magnesium Calibration
        if isempty(Mg) == 0 
            if str2double(raw_cal(Mg,J)) < CorrCoeff_threshold
                disp('WARNING: STD correlation for Magnesium below 0.995')
                calib_flags = AddFlagIC(calib_flags,'S-Mg');
            end
        elseif isempty(Mg) == 1 
            disp('ERROR: no correlation coefficient for Magnesium - please re-export data from IC computer')
            calib_flags = AddFlagIC(calib_flags,'Sx-Mg');
        end
        
        % Check Calcium Calibration
        if isempty(Ca) == 0 
            if str2double(raw_cal(Ca,J)) < CorrCoeff_threshold
                disp('WARNING: STD correlation for Calcium below 0.995')
                calib_flags = AddFlagIC(calib_flags,'S-Ca');
            end
        elseif isempty(Ca) == 1 
            disp('ERROR: no correlation coefficient for Calcium - please re-export data from IC computer')
            calib_flags = AddFlagIC(calib_flags,'Sx-Ca');
        end
      
        
        %--------- index different data types -------------
        H2O_ind = find(contains(labels_all,'H2O'));
        Samples_ind = find(contains(labels_all,'-'));
        STD_ind = find(contains(labels_all,'STD'));
        LB_ind = find(contains(labels_all,'-LB')); % index for lab blanks (started this Fall 2017)

        if isempty(LB_ind) == 1
            LB_ind = find(contains(labels_all,'LAB'));
        end
        if isempty(LB_ind) == 1
            LB_ind = find(contains(labels_all,'BLANK'));
        end
        
        if isempty(LB_ind) == 0 % checks if lab blanks in run, if not in loop script will stop - no lab blanks prior to Fall 2017
            for i = 1:length(LB_ind)
                ind(i) = find(Samples_ind == LB_ind(i));
            end
            Samples_ind = removerows(Samples_ind, ind); % removes lab blanks from index to samples
        end
        
        % separate filters end with 'N' as they are not for master files
        if_nylon = zeros(size(Samples_ind));
        for i = 1:length(Samples_ind)
            tfilter = labels_all{Samples_ind(i)};
            if tfilter(end) == 'N'
                if_nylon(i) = 1;
            end
        end
        ind = find(if_nylon==1);
        nylon_ind = Samples_ind(if_nylon==1);
        Samples_ind = removerows(Samples_ind, ind);  clear if_nylon i ind

        % remove rows with extra dashes
        num_dash = zeros(size(Samples_ind));
        for i = 1:length(Samples_ind)
            num_dash(i) = count(labels_all(Samples_ind(i)),"-");
        end
        extra_dash = find(num_dash > 1);
        Samples_ind = removerows(Samples_ind, extra_dash); clear extra_dash ind
        
        % ----Get data for different data type (standards, samples, lab blanks, and waters) ------
        data_STD = data_all(STD_ind,:);
        data_H2O = data_all(H2O_ind(2:end),:); % first water is often not representative of run, therefore start at second water
        data_teflon = data_all(Samples_ind,:);
        data_nylon = data_all(nylon_ind,:);
%         data_samples(isnan(data_samples)) = 0; % if there is no ion concentration found, set to zero rather than NaN in file. NaN indicates a missing number
        data_LB = data_all(LB_ind,:);
        area_samples = data_pre(Samples_ind,:);

        labels_stds = labels_all(STD_ind);
        labels_teflon = labels_all(Samples_ind);

        labels_nylon = labels_all(nylon_ind);
        labels_nylon2 = labels_nylon;
        for i = 1:length(labels_nylon)
            tlabel = char(labels_nylon(i));
            tlabel = tlabel(1:end-1); % remove 'N'
            labels_nylon(i) = tlabel;
        end

        % identify sites 
        labels_char = char(labels_all([Samples_ind; nylon_ind]));
        site_IDs = unique(labels_char(:,1:4),'rows');
        labels_LB = labels_all(LB_ind);

        clear *_ind num_cal raw_cal txt_cal 
        MDL_flags = cell(size(data_teflon,1),1);
        
        % This check is remove as we will report MDL in the public files.
%         % --------- Check that data is above MDL, if not, set to 0 ----------
%         MDL_flags = cell(size(data_samples,1),1);
%         for i = 1:size(cations_MDL,2)
%             idx = find(data_samples(:,i) < cations_MDL(i));
%             if ~isempty(idx)
%                 fprintf('FYI: there are %d of %d samples below MDL for %s  \n', length(idx), size(data_samples,1), cations{i})
%             end
%             clear idx
%         end

        %-------- This section check if the QC STD 1.25 concentration is within +/-10 % of the specified concentration for each anion -------
        % (data before May 2018 will not have these)
        QC_STD_ind = find(contains(labels_stds,'1.25'));
        
        if isempty(QC_STD_ind) == 0
            data_STD_QC = data_STD(QC_STD_ind, :);
            QC_avgpercentdiff = 100*((abs(mean(data_STD_QC,1) - cation_QC_conc))./cation_QC_conc);
            for i = 1:length(QC_avgpercentdiff)
                if QC_avgpercentdiff(i) > 10
                    fprintf('WARNING: Avg QC STD 1.25 for %s exceeds +/- 10 %% of specified concentration \n', cations{i})
                    
                    calib_flags = AddFlagIC(calib_flags,sprintf('qc1-%s',cations{i}) );
                end
            end
            clear i 
        end
        
        %-------- This section check if average anion concentrations in waters is < 10 * MDL  -------
        H2O_avg = mean(data_H2O,1,'omitnan');
        for i = 1:size(H2O_avg,2)
            if H2O_avg(i) > 10*cations_MDL(i)
                fprintf('WARNING: Average H2O for %s is > 10*MDL \n', cations{i})
                calib_flags = AddFlagIC(calib_flags,sprintf('qc2-%s',cations{i}) );
            end
        end
        clear i 
        
        % Create proper flag matrix for writing to data file
        if size(calib_flags,1) > 0
            for t = 1:length(data_teflon)
                MDL_flags{t} = AddFlagIC(MDL_flags{t},calib_flags);
            end
            fprintf('Cation flags for this file are: %s \n',calib_flags)
        else
            disp('No flags for this file')
        end
        
        disp('All quality assurance for this file has been completed')
        
        % -------- This section checks for existance of a master data file for each site in the IC data file ----------
        disp('Looking for master files')
        NN = 0;
        for i = 1:size(site_IDs,1)
            master_file = sprintf('%s/%s_master.csv',direc_master,site_IDs(i,:));
             
            [Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs, projectIDs_master,Master_hours, Master_masstype, ...
                Master_dates, Master_mass, Master_IC, Master_ICP, Master_XRF,...
                Master_carbon, Master_Nylon, Master_Method, Master_flags] = ReadMaster(master_file,site_IDs(i,:));

            if ~isempty(Master_IDs)

               % read area file 
               area_file = sprintf('%s/%s_IC_Area.xlsx',direc_area,site_IDs(i,:));
               area_table = read_and_update(area_file,Master_IDs,CartridgeIDs_master,Master_dates,Master_masstype);
               % read mdl file
               mdl_file = sprintf('%s/%s_IC_MDL.xlsx',direc_mdl,site_IDs(i,:));
               mdl_table = read_and_update(mdl_file,Master_IDs,CartridgeIDs_master,Master_dates,Master_masstype);

                % find index of teflon filter in master file
                samples_ICfile_ind = find(contains(labels_teflon,site_IDs(i,:))); % index of samples for site_ID(i) in IC file

                masterID_teflon_ind = [];
                for k = 1:length(samples_ICfile_ind)
                    if isempty( find(matches(Master_IDs, labels_teflon(samples_ICfile_ind(k))),1) ) == 0
                        masterID_teflon_ind(k) = find(matches(Master_IDs, labels_teflon(samples_ICfile_ind(k)))); % finds the row index in the master file for samples in IC file
                        teflon_IC_Master(k) = samples_ICfile_ind(k);
                    else
                        tfilterid = char(labels_teflon(samples_ICfile_ind(k))) ;
                        if length(tfilterid) == 8
                            tfilterid2 = strcat(tfilterid(1:5),'0',tfilterid(6:8));
                            masterID_teflon_ind(k) = find(matches(Master_IDs, tfilterid2)); % finds the row index in the master file for samples in IC file
                            teflon_IC_Master(k) = samples_ICfile_ind(k);
                        else
                            warning( '%s not found in Master File\n',labels_teflon(samples_ICfile_ind(k)) )
                            NN = NN + 1;
                        end
                    end
                end

                % find index of nylon filter in master file
                samples_ICfile_ind = find(contains(labels_nylon,site_IDs(i,:))); % index of samples for site_ID(i) in IC file

                masterID_nylon_ind = [];
                for k = 1:length(samples_ICfile_ind)
                    if isempty( find(matches(Master_IDs, labels_nylon(samples_ICfile_ind(k))),1) ) == 0
                        masterID_nylon_ind(k) = find(matches(Master_IDs, labels_nylon(samples_ICfile_ind(k)))); % finds the row index in the master file for samples in IC file
                        nylon_IC_Master(k) = samples_ICfile_ind(k);
                    else
                        tfilterid = char(labels_nylon(samples_ICfile_ind(k))) ;
                        if length(tfilterid) == 8
                            tfilterid2 = strcat(tfilterid(1:5),'0',tfilterid(6:8));
                            masterID_nylon_ind(k) = find(matches(Master_IDs, tfilterid2)); % finds the row index in the master file for samples in IC file
                            nylon_IC_Master(k) = samples_ICfile_ind(k);
                        else
                            warning( '%s not found in Master File\n',labels_nylon(samples_ICfile_ind(k)) )
                            NN = NN + 1;
                        end
                    end
                end
                       
                % get extraction volume
                elog_index = find(contains(Sheets_elog,site_IDs(i,:)) == 1);
                Elog_SheetName = char(Sheets_elog(elog_index));
                cell_elog=readcell(Elog_filename,'Sheet',Elog_SheetName);
                
                if isempty(masterID_nylon_ind) == 0
                    % nan out -xxxxx before writing to master files
                    temp = data_nylon(nylon_IC_Master,:);
                    temp(temp<-100) = nan;
                    data_nylon(nylon_IC_Master,:) = temp;

                    [extraction_volumes_N, extraction_dates_N] = GetVolume(cell_elog,labels_nylon2(samples_ICfile_ind));
                    % adding IC data to Master_IC
                    Master_Nylon(masterID_nylon_ind, 8:13) = data_nylon(nylon_IC_Master,:).*extraction_volumes_N;
                    
                    % save to master file
                    WriteToMaster( Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs, projectIDs_master,...
                        Master_hours, Master_masstype, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
                        Master_carbon,Master_Nylon, Master_Method, Master_flags,...
                        direc_master,site_IDs(i,:))

                end
                       
                if isempty(masterID_teflon_ind) == 0
                    [extraction_volumes, extraction_dates] = GetVolume(cell_elog,Master_IDs(masterID_teflon_ind));
                    
                    % find any zeros and print to 'IC_zero.txt'
                    find_ic_zeros(data_teflon(teflon_IC_Master,:), Master_IDs(masterID_teflon_ind), ...
                        extraction_dates, Master_masstype(masterID_teflon_ind), ic_zero_fname, data_type)
                    % nan out -xxxxx before writing to master files
                    temp = data_teflon(teflon_IC_Master,:);
                    temp(temp<-100) = nan;
                    data_teflon(teflon_IC_Master,:) = temp;

                    % adding IC data to Master_IC
                    Master_IC(masterID_teflon_ind, 8:13) = data_teflon(teflon_IC_Master,:).*extraction_volumes;
        
                    % adding area to area file
                    for ii = 1:6
                        if size(area_table(masterID_teflon_ind,5+7+ii)) == size(area_samples(teflon_IC_Master,ii))
                            area_table(masterID_teflon_ind,5+7+ii)= table(area_samples(teflon_IC_Master,ii));
                        else
                            error('size of area data do not match')
                        end
                    end

                    % adding mdl to mdl file
                    for ii = 1:6
                        mdl_table(masterID_teflon_ind,5+7+ii) = table(mdl_calc(1,ii).*extraction_volumes);
                        
                    end

                    % adding flag to master_initial_flag
                    if exist('Flags','var')
                    for k = 1:numel(masterID_teflon_ind)
                        Master_flags{masterID_teflon_ind(k)}  = AddFlag( Master_flags{masterID_teflon_ind(k)} , MDL_flags{teflon_IC_Master(k)} );
                    end
                    end
                    clear ind index filler
                    
                    % --- finish this site ---
                    % save to master file
                    WriteToMaster( Titles,Master_IDs, Master_Barcodes, CartridgeIDs_master, LotIDs, projectIDs_master,...
                             Master_hours, Master_masstype, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
                             Master_carbon,Master_Nylon, Master_Method, Master_flags,...
                             direc_master,site_IDs(i,:))
                        
                    % save to area file
                    writetable(area_table,area_file)

                    % save to mdl file
                    writetable(mdl_table,mdl_file)

                    fprintf('Finished writing cation data to %s master data file \n', site_IDs(i,:))
                    
                else
                    fprintf('Samples in IC file not found in %s Master file \n', site_IDs(i,:));
                    disp('Moving on to next site')
                end

                clear   Master_flags Master_IDs masterID_teflon_ind samples_ICfile_ind area_table Master_IC
                 
            end
            
            clear   txt_elog raw_elog num_elog extraction_volumes LotIDs_master 
            clear Master_IDs masterID_teflon_ind masterID_nylon_ind samples_ICfile_ind teflon_IC_Master nylon_IC_Master CartridgeIDs_master LotIDs
            clear Master_flags Master_flags_initial master_file hours_sampled Master_hours projectIDs_master Master_Barcodes 
            
        end
        
        for i = 1:size(site_IDs,1)
            file_destination = strcat(direc_archive,sprintf('%s/',site_IDs(i,:)));
            status = CopyFile(filename,file_destination);
            clear file_destination
        end
        clear i
        
        if test == 1
            file_destination = strcat(direc_archive,'Test_data/');
            status = CopyFile(filename,file_destination);
            clear file_destination
            disp('Test data detected, file has been saved to test data folder')
        end
        
        if isempty(find(status == 0) == 1) && NN == 0
            delete(filename)
            disp('File has been deleted from raw IC data folder')
        else
            if NN == 0
                disp('Cannot delete file from raw IC data folder as it has not been moved to all site folders')
            else
                disp('Cannot delete file from raw IC data folder as some samples have not been assigned to master files')
            end
        end
        clear status
        
    end
    
    %%%%%%%%%%%%%%%%%%%%% END OF CATION SECTION OF SORTING AND QC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear cal_corrcoeff data_all data_H2O data_LB data_teflon data_STD data_STD_QC  fileID calib_flags h H2O_avg 
    clear k labels_all labels_LB  labels_teflon labels_stds Master_data Master_data_initial
    clear Master_IDs masterID_teflon_ind QC_avgpercentdiff QC_STD_ind Flags samples_ICfile_ind site_IDs data_post data_pre
    clear t time_now Titles J F Cl Br NO2 NO3 PO4 SO4 Master_flags Master_flags_initial master_file num_dash status
    
    
end

fprintf('END of IC data processing for %s \n', datestr(now))
diary off

disp('Program finished')
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_post= getpostdata(data_pre,raw_cal, CubicIndicator, data_all,IC_low_cat,IonID)
    
if isempty(raw_cal)
   data_post = NaN.*data_pre;
   data_post(isnan(data_post)) = -77777; % no calibration curve

else
    if contains(raw_cal{3},'Lin') % linear regression
        offset = str2double(raw_cal{5});
        slope = str2double(raw_cal{6});
        data_post =  (data_pre-offset)./slope;

        data_post(isnan(data_post)) = -88888; % no signal

    elseif contains(raw_cal{3},'Cubic')
    
        if CubicIndicator == 1
            b0 = str2double(raw_cal{5});
            b1 = str2double(raw_cal{6});
            b2 = str2double(raw_cal{7});
            b3 = str2double(raw_cal{9});
                    
            for i = 1:length(data_pre)
                if ~isnan(data_pre(i,1))
                    syms x
                    eqn = data_pre(i,1) == b0+b1*x+b2*x^2+b3*x^3;
                    a = double( vpasolve(eqn,x,[0 4.51]) );
                    
                    if ~isempty(a)
                        data_post(i,1) =  double( vpasolve(eqn,x,[0 4.51]) );  % the upper limite of NH4 conc is 3 ug/mL. but for other species it could be 4.5 ug/mL 

                    else
                        % see if there is a solution outside of the up limit of current calibration curve. 
                        % Very high NH4 concentration. Using a specific coefficient set (from Chris/Haihui, ammonium.xlsx):
                        b00 = 0.0645945;
                        b11 = 0.241887;
                        b22 = -0.006855;
                        b33 = 0.00011;
                        eqn = data_pre(i,1) == b00+b11*x+b22*x^2+b33*x^3;
                        a = double( vpasolve(eqn,x,[0 31]) );
                        if ~isempty(a) 
                            data_post(i,1) =  double( vpasolve(eqn,x,[0 31]) ); % these coefficients valid up to 31 mg/L
                        else % there is not solution (unlikely to reach this point)
                            data_post(i,1) =  -99999; % no solution to the signal
                        end

                    end
                else
                    data_post(i,1) = -88888; % no signal
                end
            end
         
        else 
            data_post =  data_all;
            data_post(isnan(data_post)) = -66666; % cubic curve, no enough coefficient
        end
        
    end

% Commented this out since we want to leave the decision to the users:
%     data_post(data_post<0) = 0; % there should not be any negative conc. 


    % *********** Haihui Zhu, Oct 2023 ************************************
    % For very low concentration (e.g. field blank), apply the low concentration calibration curves

    % NOTE: using interp1, instead of curve fitting, to estimate ion
    % concentration
    switch IonID
        case 'Ca'
            data_post = check_for_low_conc(IC_low_cat.Ca, data_post, data_pre);
        case 'Na'
            data_post = check_for_low_conc(IC_low_cat.Na, data_post, data_pre);
        case 'Mg'
            data_post = check_for_low_conc(IC_low_cat.Mg, data_post, data_pre);
        case 'K'
            data_post = check_for_low_conc(IC_low_cat.K, data_post, data_pre);
        case 'Li'
            data_post = check_for_low_conc(IC_low_cat.Li, data_post, data_pre);
        case 'Br'
            data_post = check_for_low_conc(IC_low_cat.Br, data_post, data_pre);
        case 'Cl'
            data_post = check_for_low_conc(IC_low_cat.Cl, data_post, data_pre);
        case 'F'
            data_post = check_for_low_conc(IC_low_cat.F, data_post, data_pre);
        case 'NO2'
            data_post = check_for_low_conc(IC_low_cat.NO2, data_post, data_pre);
        case 'NO3'
            data_post = check_for_low_conc(IC_low_cat.NO3, data_post, data_pre);
        case 'PO4'
            data_post = check_for_low_conc(IC_low_cat.PO4, data_post, data_pre);
        case 'SO4'
            data_post = check_for_low_conc(IC_low_cat.SO4, data_post, data_pre);

    end

    % *********************************************************************
end

end

function data_post = check_for_low_conc(curv, data_post, data_pre)
% if the concentration is lower than 85% of the upper limit of the low_conc_curve, apply the low_conc_curve to get a more accurate conc.  
ind = find(data_post< 0.85*curv.conc(end) & data_post > -100);

if ~isempty(ind)
    data_post(ind) = interp1(curv.signal, curv.conc, data_pre(ind));
end

end

function [data_pre,labels_all,data_all,data_type]= ReadRawIC(filename,SheetName,file)

    if contains(SheetName,'ANION')
        data_type = 1;
        % Setup the Import Options and import the data
        opts = spreadsheetImportOptions("NumVariables", 9);
        opts.DataRange = "A1:I999"; % number is large enough to include all the data

    elseif contains(SheetName,'CATION')
        data_type = 2;
        % Setup the Import Options and import the data
        opts = spreadsheetImportOptions("NumVariables", 8);
        opts.DataRange = "A1:H999"; % number is large enough to include all the data
    end

    % Specify sheet and range
    opts.Sheet = SheetName;

    % Specify column names and types
    opts.VariableNames = ["No", "Name", "Amount", "Amount1", "Amount2", "Amount3", "Amount4", "Amount5", "Amount6"];
    opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double"];

    % Specify variable properties
    opts = setvaropts(opts, "Name", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "Name", "EmptyFieldRule", "auto");

    % Import the data
    initialdata = readtable(filename, opts, "UseExcel", false);

    labels_i = table2array(initialdata(:,2));
    data_i = table2array(initialdata(:,3:end));

    %------ find location of data and labels within the txt and num matrices----
    ind1 = find(contains(labels_i,'Name')) + 4; % ind1 is row location of first sample in units ug/mL and in uS*min
    ind2 = find(contains(labels_i,'Sum:')) - 1; % ind2 is row location of last sample in units ug/mL and in uS*min
    labels_all = labels_i(ind1(1):ind2(1));  % label_all is all data that will need to be reported
    data_all   = data_i  (ind1(1):ind2(1), :);
    labels_pre = labels_i(ind1(2):ind2(2)); % sometimes labels_pre (sample in the signal panel) can be more than label_all.
    data_pre   = data_i  (ind1(2):ind2(2),:);
    % Need to filter out the extra data in labels_pre and data_pre.
    if length(labels_all) ~= length(labels_pre)
    warning('There are more signal data than exported data in %s\n         Ignore the last %d samples in signal penal',...
             file,length(labels_pre)-length(labels_all))

    % remove extra data in the signal panel
    data_pre = data_pre(1:length(labels_all),:);

    end

end

function [raw_cal,CubicIndicator]= ReadRawCal(filename,SheetName)

% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 9);

% Specify sheet and range
opts.Sheet = SheetName;
opts.DataRange = "A1:I150";

% Specify column names and types
opts.VariableNames = ["VarName1", "CalibrationBatchReport", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify variable properties
opts = setvaropts(opts, ["VarName1", "CalibrationBatchReport", "VarName3"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName1", "CalibrationBatchReport", "VarName3"], "EmptyFieldRule", "auto");

% Import the data
raw_cal = readtable(filename, opts, "UseExcel", false);

% if CubicIndicator = 1, we could have a cubic relationship and enough
% coefficient; if it = 0 or NaN, we won't have enough coefficient
    % CubicIndicator = str2double(raw_cal(6,2));
% try a more robust way to find CubicIndicator
raw_cal_name = table2cell(raw_cal(:,1));
for i = 1:length(raw_cal_name)
    ind = contains(raw_cal_name{i},'Calibration Summary');
    if ind == 1
        CubicIndicator = str2double(table2array(raw_cal(i,2)));
    end
end
raw_cal = table2array(raw_cal);

end

function files = getFiles(directory)
    % retrieve all files in specified directory
    names = dir(sprintf('%s/*.xls',directory));
    names2 = dir(sprintf('%s/*.xlsx',directory));
    
    files1 = {names.name}; 
    files2 = {names2.name};
    
    files = [files1';files2'];
end

function [row, col] = findIndexContaining(string, file)
     % find row and column index of string in file
     index = find(contains(file,string));
     [row, col] = ind2sub(size(file),index); clear index
end

function [ I, J ] = FindIdx(raw_cal,k,SpecName)

% [I,J] = ind2sub( size(raw_cal),find(contains(raw_cal(1:k,1:size(raw_cal,2)),SpecName)) );
[I,J] = ind2sub( size(raw_cal),find(contains(raw_cal(1:k,1),SpecName)) );

if isempty(I)    
    warning('Cannot find %s correlation coefficient in the calibration sheet. Check species name or the sheet\n',SpecName)
end

end

function flags = AddFlagIC(oldflag,newflag)

if isempty(oldflag) 
    flags = newflag;
else 
    flags = sprintf('%s %s',oldflag,newflag);
end

end

function  out_table = read_and_update(in_file,Filter_ID,CartridgeID,Master_date_hour,masstype)
% start date & time
master_hour = floor(Master_date_hour(:,4));
master_minute = floor(60*(Master_date_hour(:,4)-master_hour));
master_second = floor(  60*(60*(Master_date_hour(:,4)-master_hour)-master_minute)   );
start_datetime = datetime([Master_date_hour(:,1:3), master_hour,master_minute,master_second]);
% stop date & time
master_hour = floor(Master_date_hour(:,8));
master_minute = floor(60*(Master_date_hour(:,8)-master_hour));
master_second = floor(  60*(60*(Master_date_hour(:,8)-master_hour)-master_minute)   );
stop_datetime = datetime([Master_date_hour(:,5:7), master_hour,master_minute,master_second]);
           
if exist(in_file,'file') 
    out_table = readtable(in_file);

    if length(out_table.Filter_ID) < length(Filter_ID) % area file has less filters than the master file
        F = nan(size(start_datetime));
        Cl = nan(size(start_datetime));
        NO2 = nan(size(start_datetime));
        Br = nan(size(start_datetime));
        NO3 = nan(size(start_datetime));
        PO4 = nan(size(start_datetime));
        SO4 = nan(size(start_datetime));
        Li = nan(size(start_datetime));
        Na = nan(size(start_datetime));
        NH4 = nan(size(start_datetime));
        K = nan(size(start_datetime));
        Mg = nan(size(start_datetime));
        Ca = nan(size(start_datetime));

        % fill out the existing values
        for ii = 1:length(Filter_ID)
            ind = find(ismember(out_table.Filter_ID,Filter_ID{ii} ));
            if ~isempty(ind)
                F(ii)   = out_table.F(ind);
                Cl(ii)  = out_table.Cl(ind);
                NO2(ii) = out_table.NO2(ind);
                Br(ii)  = out_table.Br(ind);
                NO3(ii) = out_table.NO3(ind);
                PO4(ii) = out_table.PO4(ind);
                SO4(ii) = out_table.SO4(ind);
                Li(ii)  = out_table.Li(ind);
                Na(ii)  = out_table.Na(ind);
                NH4(ii) = out_table.NH4(ind);
                K(ii)   = out_table.K(ind);
                Mg(ii)  = out_table.Mg(ind);
                Ca(ii)  = out_table.Ca(ind);
            end
        end


        % write to new IC area file
        out_table = table(Filter_ID,CartridgeID,masstype,... 1-3
            start_datetime,stop_datetime,... 4-5
            F,Cl,NO2,Br,NO3,PO4,SO4,... 6-12
            Li,Na,NH4,K,Mg,Ca); % 13-18
    end

else % not existing area/mdl file, create one

    F = nan(size(start_datetime));
    Cl = nan(size(start_datetime));
    NO2 = nan(size(start_datetime));
    Br = nan(size(start_datetime));
    NO3 = nan(size(start_datetime));
    PO4 = nan(size(start_datetime));
    SO4 = nan(size(start_datetime));
    Li = nan(size(start_datetime));
    Na = nan(size(start_datetime));
    NH4 = nan(size(start_datetime));
    K = nan(size(start_datetime));
    Mg = nan(size(start_datetime));
    Ca = nan(size(start_datetime));

    % write to new IC area file
    out_table = table(Filter_ID,CartridgeID,masstype,... 1-3
        start_datetime,stop_datetime,... 4-5
        F,Cl,NO2,Br,NO3,PO4,SO4,... 6-12
        Li,Na,NH4,K,Mg,Ca); % 13-18
end


end

function [extraction_volumes, extraction_dates] = GetVolume(cell_elog,Master_IDs)
elog_ID = cell_elog(:,3);
for i = 1:length(elog_ID)
   if ismissing( elog_ID{i} )
       elog_ID{i} = 'Missing';
   end
end

% find the filter ind in elog file
Ind = nan(length(Master_IDs),1);

for ii = 1:length(Master_IDs)
    a = find( matches( elog_ID, Master_IDs(ii) ) );
    if length(a)==1
        Ind(ii) = a;
    elseif length(a) > 1
        error('%d record found in elog for %s',length(a),Master_IDs(ii))
    else
        tfilter = char(Master_IDs(ii));
        tfilter2(1:5) = tfilter(1:5);
        tfilter2(6:8) = tfilter(7:9);
        a = find( matches( elog_ID, tfilter2 ) );
        if ~isempty(a)
            Ind(ii) = a;
        else
            error('%s extraction vol not found in elog',tfilter2)
        end
    end
end


% selected volumn, IDs, and dates
volumes = cell_elog(Ind,6) ;
IDs = cell_elog(Ind,3);
dates = cell_elog(Ind,5);
N = 0;
for i = 1:length(volumes)
    if ~isa(volumes{i},'double')
        warning('Extraction Volumes for %s is %s.\n ',IDs{i},volumes{i})
        N = N+1;
    end
end

if N>0
    error('Check E-log for invalid volume.\n')
else
    extraction_volumes = cell2mat( volumes );
    extraction_dates   = dates ;
end

end


