clear all
close all
clc

%%

load('/Home/siv29/dsa030/Desktop/GEANT4/TGF-TEB-Propagation-Geant4/tests/data_t50_rd_15km.mat');

histogram('BinEdges',bins_rad_dist,'BinCounts',t50,'DisplayStyle','stairs','LineWidth',2);
hold on

load('/Home/siv29/dsa030/Desktop/GEANT4/TGF-TEB-Propagation-Geant4/tests/data_t50_rd_12km.mat');

histogram('BinEdges',bins_rad_dist,'BinCounts',t50,'DisplayStyle','stairs','LineWidth',2);
hold on

legend('source at 15 km','source at 12 km')

xlabel('TGF ISS radial distance (km)')
ylabel('T_{50} duration (micro-second)')