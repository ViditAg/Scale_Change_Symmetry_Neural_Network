%To calculate the connectivity matrix for realistic model
%Input-> l,b, dimensions of network lattice
%Output-> B: Connectivity Matrix
function B=RealConnect(l,b)
N=l*b; % Number of neurons
% short-range inhibition and long range excitation
xpos=[]; for i=1:l; xpos=[xpos; ones(l,1)*i]; end %x positions
ypos=[]; for i=1:b; ypos=[ypos; (1:b)']; end %y positions
dmat=squareform(pdist([xpos ypos])); %pairwise distances
Ilist = rand(N,1)<0.2;% 20% neurons set to inhibitory
C_E=2;C_I=3; %sigma for (ext./inh.) gaussian weight dist.
B=single(zeros(N));
B(:,~Ilist)=exp(-(dmat(:,~Ilist)/C_E).^2);% exct. weights
B(:,Ilist)=-exp(-(dmat(:,Ilist)/C_I).^2);% inh. weights
for bi=1:N
    B(bi,bi)=0;% no self-activation
end