%To get data from expt readings and convert to lattice data
%function-> ExptData
%Input-> ExptNum, Experiment number 2 to 10,
%fileNum, file number for each expt.; zthresh: threshold on voltage signal
%Output-> lattice and (nev) 1D time series of point-process neural data
function [lattice,nev]=ExptData(ExptNum,fileNum,zthresh)
if exist('brainWindow')==0
    load('brain_window.mat');
end
fname = sprintf('%s%03d%s%d%s','Exp',ExptNum,'_',fileNum,'Data.mat');
load(fname);
% point process binary data and thresholded binary data
[ pp, thresh ] = ppmaker( ratioSequenceFiltered, zthresh);
pp=single(pp(:,:,201:end));
l=192;b=128;
Time=size(pp,3);
lattice=false(l,b,Time);
nev=false(l*b,Time);
% filtering the data by brain window
for t=1:Time
    pp_step1=logical(pp(:,:,t).*brainWindow);
    pp_step2=pp_step1(21:212,100:227);%;
    lattice(:,:,t)=pp_step2;
    % converting lattice to 1D data
    nev(:,t)=pp_step2(:);
end