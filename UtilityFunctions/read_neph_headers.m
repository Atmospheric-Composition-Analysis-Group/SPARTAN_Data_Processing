function [headers, opts] = read_neph_headers(filename)


opts = detectImportOptions(filename);
if length(opts.VariableNames) == 30 % this means the script might have ingored empty columns at the end (CR columns)
    opts = delimitedTextImportOptions("NumVariables", 36);
end
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];
opts.VariableNamingRule = 'preserve';

rawtable = readtable(filename, opts);

headers = rawtable.Properties.VariableNames;