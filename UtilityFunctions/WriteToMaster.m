function [] = WriteToMaster( Titles,Master_IDs, Master_Barcodes, Master_CartridgeIDs, Master_LotIDs, Master_ProjectIDs,...
                             Master_hours, Master_masstype, Master_dates, Master_mass,Master_IC, Master_ICP, Master_XRF,...
                             Master_carbon, Master_Nylon, Master_Method, Master_flags,...
                             direc_master,site_ID)
% This function will write all the input data into the master file 
% Inputs:
% direc_master should be the dir for the master file
% Site ID shoule be the 4 letter Site Code.
% Other inputs (master data) should be same as listed in the ReadMaster
% function

% Written by: Haihui Zhu
% Created: 2021-08-25


% Outputs:
% Titles = 96 x 1 cell 
% Master_IDs = N x 1 cell 
if size(Master_IDs,2) ~= 1
    error('Master ID column format wrong: should contain 1 column, but it has %d colums.', size(Master_IDs,2))
end
% Master_Barcodes = N x 1 cell
if size(Master_Barcodes,2) ~= 1
    error('Master barcode column format wrong: should contain 1 column, but it has %d colums.', size(Master_Barcodes,2))
end
% Master_CartridgeIDs = N x 1 cell
if size(Master_CartridgeIDs,2) ~= 1
    error('Master cartridge column format wrong: should contain 1 column, but it has %d colums.', size(Master_CartridgeIDs,2))
end
% Master_LotIDs = N x 1 cell
if size(Master_LotIDs,2) ~= 1
    error('Master lot id column format wrong: should contain 1 column, but it has %d colums.', size(Master_LotIDs,2))
end
% Master_ProjectIDs = N x 1 cell
if size(Master_ProjectIDs,2) ~= 1
    error('Master project id column format wrong: should contain 1 column, but it has %d colums.', size(Master_ProjectIDs,2))
end
% Master_hours = N x 1 mat; hour sampled
if size(Master_hours,2) ~= 1
    error('Master hours format wrong: should contain 1 column, but it has %d colums.', size(Master_hours,2))
end
% Master_masstype = N x 1 mat 
if size(Master_masstype,2) ~= 1
    error('Master masstype format wrong: should contain 1 column, but it has %d colums.', size(Master_masstype,2))
end
% Master_dates = N x 8 mat; Start year month day and End year month day
if size(Master_dates,2) ~=8
    error('Master lot id column format wrong: should contain 8 column, but it has %d colums.', size(Master_dates,2))
end
% Master_mass = N x 2 mat; 'mass_ug','Volume_m3'
if size(Master_mass,2) ~= 2
    error('Master mass format wrong: should contain 2 column, but it has %d colums.', size(Master_mass,2))
end
% Master_IC = N x 13 mat
if size(Master_IC,2) ~= 13
    error('Master IC format wrong: should contain 13 column, but it has %d colums.', size(Master_IC,2))
end
% Master_ICP = N x 21 mat
if size(Master_ICP,2) ~= 21
    error('Master ICP format wrong: should contain 21 column, but it has %d colums.', size(Master_ICP,2))
end
% Master_XRF = N x 26 mat
if size(Master_XRF,2) ~= 26
    error('Master XRF format wrong: should contain 26 column, but it has %d colums.', size(Master_XRF,2))
end
% Master_carbon = N x 4 mat; 'BC_SSR_ug', 'BC_HIPS_ug','EC_FTIR_ug','OC_FTIR_ug'
if size(Master_carbon,2) ~= 4
    error('Master carbon format wrong: should contain 4 column, but it has %d colums.', size(Master_carbon,2))
end
% Master_Nylon = N x 13 mat
if size(Master_Nylon,2) ~= 13
    error('Master Nylon format wrong: should contain 13 column, but it has %d colums.', size(Master_Nylon,2))
end
% Master_Method = N x 1 mat
if size(Master_Method,2) ~= 1
    error('Master method format wrong: should contain 1 column, but it has %d colums.', size(Master_Method,2))
end
% Master_flags = N x 1 cell
if size(Master_flags,2) ~= 1
    error('Master flag format wrong: should contain 1 column, but it has %d colums.', size(Master_flags,2))
end


% If the script reaches here, the test write is successful. Can start to write to master files:

    fileID = fopen(sprintf('%s/%s_master.csv',direc_master,site_ID),'w');
    if fileID == -1
        error('Cannot open master file: %s',site_ID)
    end
    % Print Titles
    printTitles(fileID, Titles)

    % Print Data
    for h=1:size(Master_IDs,1) % 7 cols                 % 8 cols               % 2 cols                     % 13 cols                                                   % 21 cols                                                                              % 26 cols                        % 4 cols                       % 13 cols                  % 2 cols
        fprintf(fileID,'%s,%s,%s,%s,%s,%f,%6.1f,   %f,%f,%f,%f,%f,%f,%f,%f,   %6.2f,%f,  %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,%f,%f,   %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,    %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,%f,%f,%f,%f,%f, %f,%f,%f,%f,    %f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,%f,%f,   %f,%s\n',...
            char(Master_IDs(h,:)), char(Master_Barcodes(h,:)), char(Master_CartridgeIDs(h,:)), char(Master_LotIDs(h,:)), char(Master_ProjectIDs(h,:)),...
            Master_hours(h),       Master_masstype(h,:),       Master_dates(h,:),        Master_mass(h,:),         Master_IC(h,:), ...
            Master_ICP(h,:),       Master_XRF(h,:),            Master_carbon(h,:),       Master_Nylon(h,:),           Master_Method(h,:), char(Master_flags(h,:)));
    end

    fclose(fileID);
    fprintf('Finished writing to %s master file \n\n', site_ID)

end
