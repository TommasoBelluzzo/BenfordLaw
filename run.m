warning('off','all');

close('all');
clearvars();
clc();

[path,~,~] = fileparts(mfilename('fullpath'));

if (~strcmpi(path(end),filesep()))
    path_base = [path filesep()];
end

paths_base = genpath(path);
addpath(paths_base);

file = fullfile(path,strrep('Datasets\Example_Small.xlsx','\',filesep()));
data = importdata(file);
benford_analyse(data);

rmpath(paths_base);
