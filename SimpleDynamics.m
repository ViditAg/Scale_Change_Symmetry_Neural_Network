%To calculate the dynamics of simple model
%Input-> l,b: lattice dimensions, Time: simulation time
%C: Coupling strength, p: probability of external activation
%Output-> lat: Lattice activity data as time series
%spkcnts: Spike counts time series
function [lat,spkcnts]=SimpleDynamics(l,b,Time,C,p)
nb_Con=[0 1 0; 1 1 1;0 1 0]; % Nearest-neighbor Connectivity
lat=false(l,b,Time); % neural network on 2D lattice (lXb) evolve over time
lat(:,:,1) = rand(l,b)<1; % initialize neural activity data matrix
spkcnts = single(zeros(l,b,Time)); % variable to store total spike input to a neuron
for t=1:Time-1
    % active neighbors based on the network activity
    spkcnts(:,:,t)=conv2(squeeze(single(lat(:,:,t))),nb_Con,'same');
    % activation step
    lat(:,:,t+1)=(1-(1-C).^spkcnts(:,:,t)*(1-p))>rand(l,b);
end