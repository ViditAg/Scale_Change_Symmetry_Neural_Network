%To calculate the dynamics of Realistic Model
%Input-> N: network size, Time: simulation time
%B: Connectivity Matrix
%Output-> nev: neural activity data as time series
%spkcnts: synaptic input time series
function [nev,syn_input]=RealisticDynamics(N,Time,B)
p_ext=5e-6;% external noise
Tau=80;% adaptation time constant
ref_parameter=1e-3; % refractory parameter
% Initialize variables
nev = false(N,Time); % variable to store neuronal activity
h=zeros(N,Time);% history time-series
syn_input = zeros(N,Time-1); % variable to store total spike input to a neuron
nev(:,1) = rand(N,1)<1; % 1% activity
% initial simulation without adaptation
for t=1:490
    syn_input(:,t)=(B*nev(:,t))+p_ext; % estimate synaptic input
    nev(:,t+1)=syn_input(:,t)>rand(N,1); % activation based on synaptic input
end
% adaptation introduced after sufficient history is present
for t=491:Time-1
    s_tau=sum(nev(:,t-Tau+1:t),2); 
    h(:,t)=s_tau+(s_tau==0); % estimate history of a neuron
    syn_input(:,t)=((B*nev(:,t))./h(:,t))+p_ext; % synaptic input with adaptation
    syn_input(nev(:,t),t)=ref_parameter*syn_input(nev(:,t),t); % add refractory-ness
    nev(:,t+1)=syn_input(:,t)>rand(N,1); % activation based on final synaptic input
end