warning('off','all');

close('all');
clearvars();
clc();

[path,~,~] = fileparts(mfilename('fullpath'));
paths = genpath(path);
addpath(paths);

%imp = importdata('\Datasets\Example.xlsx');
%data = imp.data;

data = importdata('\Datasets\CP.xlsx');
%data = benford_random([1000 1],50000);

benford_analyse(data);

rmpath(paths);