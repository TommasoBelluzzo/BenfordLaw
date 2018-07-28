warning('off','all');

close('all');
clearvars();
clc();

[path_base,~,~] = fileparts(mfilename('fullpath'));

if (~endsWith(path_base,filesep()))
    path_base = [path_base filesep()];
end

paths_base = genpath(path_base);
addpath(paths_base);

file = fullfile(path_base,strrep('Datasets\Example_Small.xlsx','\',filesep()));
data = importdata(file);
benford_analyse(data);

rmpath(paths_base);