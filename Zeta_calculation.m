% To estimate zeta_min and k,x0 corresponding to it
% function-> Zeta_calculation
% Input-> h_nume, h_deno: fine scale spike counts
% h_nume_k1, h_deno_k1: course scale spike counts
% k,x0: arrays for transformation scheme parameters
% task: choose option for simple model or others
function [Zeta_min,k_zetamin,x0_zetamin]=Zeta_calculation(h_nume,h_deno,h_nume_k1,h_deno_k1,k,x0,task)
bins=linspace(0,5,6); % n-counts bins
% activation probability at fine scale
for nn=1:length(bins)
    if h_deno(nn)==0
        phi(nn)= 0;
    else
        phi(nn)= h_nume(nn)/h_deno(nn);
    end
end
% activation probability at coarse scale for different transformational
% scheme(different (k,x0) value combinations)
for ii=1:length(k)
    for jj=1:length(x0)
        for nn=1:length(bins)
            if h_deno_k1(nn,ii,jj)==0
                phi_k1(nn,ii,jj)= 0;
            else
                phi_k1(nn,ii,jj)= h_nume_k1(nn,ii,jj)/h_deno_k1(nn,ii,jj);
            end
        end
        % zeta calculations from phi and phi_k1
        for nn=1:length(bins)
            zeta_n(nn) = abs(phi(nn)-phi_k1(nn,ii,jj));
        end
        if task==4
           zeta(ii,jj) = sum(zeta_n);  
        else
           zeta(ii,jj) = sum(zeta_n(2:end-1));
        end 
    end
end
% minima of zeta 
Zeta_min=min(min(zeta));
for ii=1:length(k)
    for jj=1:length(x0)
        if zeta(ii,jj)==Zeta_min
            k_zetamin = k(ii); % k for zeta min
            x0_zetamin = x0(jj); % x0 for zeta min
        end
    end
end