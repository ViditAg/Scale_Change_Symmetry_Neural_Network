function [ pp, thresh ] = ppmaker( data, zthresh )
%PPMAKE Transform data into a point process
%   data is an (x,y,t) timeseries
%   zthresh is the number of SD to threshold above or below
%
%  Returns a thresholded version and a point process version.
%  
%  NOTE: if zthresh is negative ALL DATA WILL BE INVERTED (data=-data)
%  to facilitate thresholding and also subsequent calculations of 
%  magnitude
%
%  NOTE: amended by Greg to be more flexible
%
%  NOTE: "point process" has subtly different meaning in Enzo's paper,
%  i.e. whether a voxel is 'on' if it CROSSES from < to > 1SD, or whether
%  we just mean >1SD
%
%  NOTE: other versions of this stuff exist (see Eriks email with nLFP
%  stuff)
if nargin <2
    zthresh = 1;
end
   
    thresh = zeros(size(data));
    zimg = zscore(data, 0, 3); % transform to zstats voxelwise
    
    for t=1:size(data,3)
        if zthresh < 0
            thresh(:,:,t) = (zimg(:,:,t) < zthresh); 
        else
            thresh(:,:,t) = (zimg(:,:,t) > zthresh); 
        end
    end
   
    % point process forces only the 'moment' threshold crossing occurs
    % to be on
     pp=cat(3, thresh(:,:,1), ...
            and( thresh(:,:,2:end), not(thresh(:,:,1:(end-1))))); 
    
    % at this point both thresh and pp are binary values, with 1=on, 0=off
    % So that we can recover magnitude information we now multiply these
    % voxelwise by the original data, so that the magnitude information is 
    % preserved WHERE the threshold has been crossed
   % pp = data .* pp;
   % thresh = data .* thresh;
end

