warning('off','all');

close('all');
clearvars();
clc();

[path,~,~] = fileparts(mfilename('fullpath'));
paths = genpath(path);
addpath(paths);

imp = importdata('\Datasets\Example.xlsx');
data = imp.data;

benford_analyse(data);

rmpath(paths);