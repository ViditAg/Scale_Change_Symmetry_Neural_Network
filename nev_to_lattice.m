%To convert 1D neuron data to 2D lattice form
%Input-> nev: 1D time series; l,b: dimensions
%Output-> Lattice: 2D lattice time series
function lat = nev_to_lattice(nev,l,b)
for len=1:l
    for bre=1:b
        nev_pos = (len-1)*l + bre;
        lat(len,bre,:) = nev(nev_pos,:);% activity data in form of a lattice                         
    end
end