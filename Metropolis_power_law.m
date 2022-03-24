function [l,p,plm]=Metropolis_power_law(gamma,N,T,Teq,Kc,figure_handle)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  If you use this code, please cite 
%  G. Bianconi  
%  "Statistical physics of exchangeable sparse simple networks, 
%   multiplex networks and simplicial complexes"
%  Physical Review E  105, 034310 (2022). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code that generates exchangeable simple network 
%with power-law degree distribution.
% This code uses: 
% gamma power-law exponent
% N total number of nodes
% T total number of Metropolis-Hastings choices  
% Teq equilibration time after which the degree distribution is sampled
% Kc maximum degree allowed (it should be much smaller 
%     than the structural cutoff)
% figure_handle: if set equal to 1 outputs a figure of the target degree
% distribution (solid line) and the average degree distribution 
% of the model (circles)

% This code has as output:
% l vectors of possible degrees
% p  target degree distribution
% plm average degree distribution of the model
% 
%Example of how to use the code
% [l,p,plm]=Metropolis_power_law(2.5,8000,400000,300000,20,1)
%gamma=2.5
%N=8000;
%T=400000;
%Teq=300000;
%Kc=20;

l=1:Kc;

p=l.^(-gamma);
p=p/sum(p);%%%% Target degree distribution
c=sum(l.*p);%make sure c>1
L=ceil(c*N/2);
a=sparse(N,N);
%%%%%Initial condition %%%%%%%%%%%%%%%%%%%
for i=1:2:N,
    a(i,i+1)=1;
    a(i+1,i)=1;
end

for t1=(N/2):L
    i=ceil(N*rand(1));
    j=ceil(N*rand(1));
    k=sum(a);
    while((a(i,j)>0)||(abs(j-i)==0)||(k(i)>Kc-1)||(k(j)>Kc-1)),
         i=ceil(N*rand(1));
        j=ceil(N*rand(1));
    end
    a(i,j)=1;
    a(j,i)=1;
end
k=sum(a);
 
 l2=(1:1000)/1000;
[I,J,V]=find(triu(a,1));
L=numel(V);
plm=zeros(size(l));
count=0;

%%%%%%%%%%%%Metropolis_Hastings algorithm%%%%%%%%%%%%%
for it=1:T,
k=sum(a);
n=ceil(L*rand(1));
i=I(n);
j=J(n);
in=ceil(N*rand(1));
jn=ceil(N*rand(1));
while((a(in,jn)==1)||(in==jn))
in=ceil(N*rand(1));
jn=ceil(N*rand(1));
end
x=rand(1);
if((k(in)<Kc-1)&&(k(jn)<Kc-1)&&(k(i)>1)&&(k(j)>1)),
Hin=-log(p(k(i))*p(k(j))*p(k(in))*p(k(jn))*factorial(k(i))*factorial(k(j))*factorial(k(in))*factorial(k(jn)));
Hfin=-log(p(k(i)-1)*p(k(j)-1)*p(k(in)+1)*p(k(jn)+1)*factorial(k(i)-1)*factorial(k(j)-1)*factorial(k(in)+1)*factorial(k(jn)+1));
DeltaH=Hfin-Hin;
if (x<exp(-DeltaH))
    a(i,j)=0;
    a(j,i)=0;
    a(in,jn)=1;
    a(jn,in)=1;
    I(n)=in;
    J(n)=jn;
end
end


if (it>Teq)
k2=sum(a);
pl2=hist(k2,l)/N;
plm=plm+pl2;
count=count+1;
end

end
plm=plm/count;
%%%%%%%%%%%%Plot figure%%%%%%%%%%%%%%
if (figure_handle==1)
figure
loglog(l,p,'-',l,plm,'o');
 xlabel('k','FontSize',26,'FontWeight','bold')
ylabel(' p(k)','FontSize',26,'FontWeight','bold')
%ylim([10^(-5) 1]);
box on
set(gca,'FontWeight','bold','FontSize',18);
end
end



