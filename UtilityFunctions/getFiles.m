function files = getFiles(directory)
    % retrieve all files in specified directory
    names = dir(sprintf('%s/*.xlsx',directory));
    files = {names.name}; 
    files = files';

    names = dir(sprintf('%s/*.csv',directory));
    files2 = {names.name}; 
    files2 = files2';

    if ~isempty(files2)
        for ff = 1:length(files2)
            files{end+1} = files2{ff};
        end
    end

end
