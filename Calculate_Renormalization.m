task=1;
subtask=1;
rng('shuffle');
nb_Con=[0 1 0; 1 1 1;0 1 0];
if task==1
    % Simple Model Avg S
    % for spatiotemp. activity data call SimpleDynamics at C=0.2;0.23;0.3
    C=0.15:0.01:0.35;
    % Sandpile model simulation
    dim=400;%200
    Time =1.6102e4; % simulation time
    Trans=100;
    %dimensions of initial lattice
    l = dim; b = dim; N= l*b;
    % renormalization block length
    %Dynamical parameters
    p = 0.001;
    for Cstep=1:length(C)
        [lattice,Spike_Counts]=SimpleDynamics(l,b,Time,C(Cstep),p);
        AvgS(Cstep)=sum(sum(lattice(:,:,Trans+1:end)))/N;
    end
elseif task==2
    % Realistic Model mean correlation
    % for spatiotemp. activity data call RealisticDynamics at I=0.01;0.65;2.0
    Time = 1.65e4; % simulation time
    dim=160;
    l=dim;b=dim;N=l*b;
    % Connectivity Matrix for realistic Model
    B=RealConnect(l,b);
    I=0.01:0.1:2.01;
    for Istep=1:length(I)
        % tuning inhibition
        B_prime=B;
        B_prime(:,Ilist)=I(Istep)*B_prime(:,Ilist);
        [nev,synaptic_input]=RealisticDynamics(N,Time,B_prime);
        subNev_spk=nev(:,501:end)';
        CorrNev_spk=corr(subNev_spk);
        MeanCor(Istep)= nanmean(abs(CorrNev_spk(:)));
        S_spk(Istep)=mean(nev,1);
    end
elseif task==3
    % Expt Data lattice and mean correlation
    % ExptNum: experiment number runs from 2 to 10; increasing time till
    % drug
    %FileNum: Multiple readings recorded at each ExptNum
    ExptNum=2;FileNum=1;
    load('brain_window.mat');
    [lattice,nev]=ExptData(ExptNum,fileNum,zthresh,brainWindow);
    S=squeeze(sum(sum(lattice,1),2));
    % point process data correlation
    subNev=nev';
    CorrNev=corr(subNev);
    CorrNev(isnan(CorrNev))=0;
    MeanCor= nanmean(abs(CorrNev(:)));
elseif task==4
    if subtask==1
        % renormalization step
        % done at each point in parameter space over 100 random realizations
        % use the following piece to create data for 100 runs of same code
        % combine all the data by running Zeta_calculation.m
        % final answer calculated as zeta See equation.
        % Simple Model simulation
        dim=400;%200
        Time =1.6102e4; % simulation time
        Trans=100;
        %dimensions of initial lattice
        l = dim; b = dim; N= l*b;
        % renormalization block length
        %Dynamical parameters
        p = 0.001;%0.01;0.0001;
        C=0.23;%0.15:0.01:0.35;
        % transformational scheme parameters
        r=8;%4;16; % spatial dimension of transformational block
        time_r=1;%2;4;8;16; % temporal dimension of transformational block
        k=1:5:101; % steepness of transformation function f(S_b)
        x0=0.01:0.02:1; % mid point of transfromation function f(S_b)
        h_nume_k1=zeros(length(bins),length(k),length(x0));
        h_deno_k1=zeros(length(bins),length(k),length(x0));
        [Lattice,Spike_Counts]=SimpleDynamics(l,b,Time,C,p);
    elseif subtask==2
        Time = 1.65e4; % simulation time
        dim=160;
        l=dim;b=dim;N=l*b;
        nb_Con=[0 1 0; 1 1 1;0 1 0]; % connectivity nearest-neighbor for activation probability
        %Dynamical parameters
        I=0.6;%0.01:0.1:2.01
        % transformational scheme parameters
        r=8;%4;16; % spatial dimension of transformational block
        time_r=1;%2;4;8;16; % temporal dimension of transformational block
        k=1:5:101; % steepness of transformation function f(S_b)
        x0=0.01:0.02:1; % mid point of transfromation function f(S_b)
        % Connectivity Matrix for realistic Model
        B=RealConnect(l,b);
        % tuning inhibition
        B(:,Ilist)=I*B(:,Ilist);
        % initial simulation without adaptation
        [nev,synaptic_input]=RealisticDynamics(N,Time,B);
        % to use the spiking data for zeta calculation
        Lattice=nev_to_lattice(nev(:,501:end),l,b); % convert data dimensions
      % to use the continous synaptic data for zeta calculation
% %     Syn_Lattice=nev_to_lattice(synaptic_input(:,501:end),l,b); % convert data dimensions
% %     zthresh=0.5;
% %     Lattice=ppmaker(Syn_Lattice , zthresh);
      % active neighbors based on the network activity
        Spike_Counts=convn(squeeze(single(Lattice(:,:,1:end-1))),nb_Con,'same');
    elseif subtask==3
        ExptNum=2;FileNum=1;
        load('brain_window.mat');
        [Lattice,nev]=ExptData(ExptNum,fileNum,zthresh,brainWindow);
        Spike_Counts=convn(squeeze(single(Lattice(:,:,1:end-1))),nb_Con,'same'); 
    end
    [h_nume,h_deno]=nt_counts(Spike_Counts(:,:,Trans+1:end-1),Lattice(:,:,Trans+2:end));
    for ii=1:length(k)
        for jj=1:length(x0)
            [h_nume_k1(:,ii,jj),h_deno_k1(:,ii,jj)]=Renormalization(lattice(:,:,Trans+1:end),r,time_r,k(ii),x0(jj),nb_Con);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%zeta_min calculation%%%%%%%%%%%%%%%%
    mean_h_nume=squeeze(mean(h_nume,2));
    mean_h_deno=squeeze(mean(h_deno,2));
    mean_h_nume_k1=squeeze(mean(h_nume_k1,2));
    mean_h_deno_k1=squeeze(mean(h_deno_k1,2));
    [Zeta_min,k_zetamin,x0_zetamin]=Zeta_calculation(mean_h_nume,mean_h_deno,mean_h_nume_k1,mean_h_deno_k1,k,x0,task);
elseif task==5
    % avalanche distribution and kappa calculation
    % for realistic model model resuse the second half
    % with Lattice generated from codes already explained
    C=0.15:0.01:0.25;
    % Sandpile model simulation
    dim=400;%200
    Time =1.6102e4; % simulation time
    Trans=100;
    %dimensions of initial lattice
    l = dim; b = dim; N= l*b;
    % renormalization block length
    %Dynamical parameters
    p = 0.001;
    for Cstep=1:length(C)
        [lattice,~]=SimpleDynamics(l,b,Time,C(Cstep),p);
%%% %%%%%%%%%Second half: avalanche distribution%%%% size and duration%%%%%
        data=lattice(:,:,Trans+1:end);
        [avalanche,~]=avalanche(data);
        avsz=avalanche.S;
        nav=avalanche.N;
        xmin=min(avsz);
        xmax=max(avsz);
        dat=avsz(avsz>=xmin);
        n=length(dat);
        expon=1.5;
        refcdf=((xmin:xmax).^(1-expon)-xmin^(1-expon))/(xmax^(1-expon)-xmin^(1-expon)); %reference CDF
        datcdf = cumsum(hist(dat,xmin:xmax)./n); %data CDF
        ndiff=10;
        xlist=round(logspace(log10(xmin*1.1),log10(xmax/1.1),ndiff));
        kappa=1+sum(refcdf(xlist-xmin+1)-datcdf(xlist-xmin+1))/ndiff;
        %plot avalanche size PDF
        binv=unique(round(logspace(log10(xmin),log10(xmax),20*log10(xmax/xmin))));
        nb=length(binv)-1;
        num=histc(avsz,binv);
        Avalanche_data{Cstep}=avsz;
        Kappa_data(Cstep)=kappa;
        if Cstep==1
            loglog(binv(1:nb)+diff(binv),num(1:nb)./diff(binv)/nav,'k')%,'color',kcol(kind,:))
            hold on;
        elseif Cstep==2
            loglog(binv(1:nb)+diff(binv),num(1:nb)./diff(binv)/nav,'r')%,'color',kcol(kind,:))
            hold on;
        elseif Cstep==3
            loglog(binv(1:nb)+diff(binv),num(1:nb)./diff(binv)/nav,'b')%,'color',kcol(kind,:))
            hold on;
        elseif Cstep==4
            loglog(binv(1:nb)+diff(binv),num(1:nb)./diff(binv)/nav,'g')%,'color',kcol(kind,:))
            hold on;
        elseif Cstep==5
            loglog(binv(1:nb)+diff(binv),num(1:nb)./diff(binv)/nav,'c')%,'color',kcol(kind,:))
            hold on;
        end
        clear avalanche
    end
else
    disp('Enter between 1 to 6 only');
end