% Author: Mostaquimur Rahman, University of Bristol, UK
% (ar15645@bristol.ac.uk)
% This script executes the groundwater model code groundwater.m

clear
close all
clc

topography = load('test_topography.txt'); % topographic height (m)
topo_resolution = 1000; % topographic resolution (m)
recharge = load('test_recharge.txt'); % groundwater recharge (m/d)
Sy = load('test_Sy.txt'); % specific yield (-)
T = load('test_T.txt'); % transmissivity (m2/d)
ini_WTD = 100; % initial groundwater table depth (m below surface)
dt = 1; % time step (d)
total_sim_steps = 1000; % total simulation steps
r_inac = []; % row of inactive (sea) cells
c_inac = []; % column of inactive (sea) cell
sea_level = 0; % sea level hight (m)

Hold = topography - ini_WTD; % initial groundwater level (m)
for t = 1:total_sim_steps   
    % call groundwater model   
    [Hnew,runoff_md] = groundwater(topography,topo_resolution, ...
                                   recharge,Sy,T,Hold,r_inac,c_inac,dt, sea_level);
    Hold = Hnew;
    head_mOB = Hnew;   
    disp(['Timestep ' num2str(t) ' complete'])
end
% outputs: head_mOB = groundwater head elevation (m)
%          runoff_md = runoff in m/d
save('results.mat','head_mOB','runoff_md')

% plot groundwater table depth
figure
imagesc(topography-head_mOB)
colorbar