function bgrateM
clc;clear all;close all;
n0=3; % outer medium index
n1=1.1; %first layer refractive index
n2=2.5; %second layer index
n3=4; % sustrate index
h1=0.5e-6;% first layer thickess
h2=0.8e-6;% second layer thickess
lamm=[400:1:1000]*1e-9; % wavelength range

T4 = readtable('N-BK7.txt');
T1 = readtable('Air.txt');
T2 = readtable('SiO2.txt');
%T3 = readtable('Al2O3.txt');
T3 = readtable('TiO2.txt');

function n = Nmatrl(lam, T)
    Tlam=T{:,1};
    Tn=T{:,2};
    n=interp1(Tlam, Tn, lam, 'PCHIP');
    n=interp1(Tlam, Tn, lam, 'PCHIP');
end

str=[0,1,2,1,2,1,2,1,2,1,2,1,3,2,1,2,1,2,1,2,1,2,3];
hb=[0, h1, h2, 0];
nb=[n0, n1, n2, n3];
nbT={T1, T2, T3, T4};
% Medium from which light comes is str first element
% structure order. Numbering shows layer ordering. (0,1,2,...2,3) means 0 layer, 1
% layer, 2 layer ..., 3 other layer it can be a sustrate
% write in str your order


% for ilam=1:length(lamm)
%     lam=lamm(ilam)*1e9
%     nk(ilam)=Nmatrl(lam,nbT{4});
%     lk(ilam)=lam;
% end
% 
% nk2=Nmatrl(1200,T4)
% 
% wlght = T4{:,1};
% Ren = T4{:,2};
% plot(wlght, Ren, lk, nk, 'o', 1200, nk2, 'x');
% ffff
% function res = mfunc(lam)
%     res = lam+2
% end
% d = @mfunc;
% d(10)
% ddddd


for il=1:length(lamm)
    lam=lamm(il);
    il/length(lamm)
    z=0;
    M=[1,0;0,1];
    M2=[1,0;0,1];

for is=2:length(str)
    id = str(is)+1;
    idm = str(is-1)+1;
    Tid = nbT{id};
    Tidm = nbT{idm};
    switch str(is)
        case 0
            ki=2*pi*n0/lam;
            wi=n0;            
        case 1
            ki=2*pi*n1/lam;
            wi=n1;                        
        case 2
            ki=2*pi*n2/lam;
            wi=n2;
        case 3
            ki=2*pi*n3/lam;
            wi=n3;
    end
    switch str(is-1)
        case 0
            kim1=2*pi*n0/lam;
            wim1=n0;
            h=0;            
        case 1
            kim1=2*pi*n1/lam;
            wim1=n1;
            h=h1;            
        case 2
            kim1=2*pi*n2/lam;
            wim1=n2;
            h=h2;
         case 3
            kim1=2*pi*n3/lam;
            wim1=n3;
            h=0;
    end
    kim1p=2*pi/lam*nb(idm);
    hmp=hb(idm);
    wim1p=nb(idm);
    wip=nb(id);
       
    %Tidm{:,1}
    %Tidm{:,2}
    nbim=Nmatrl(lam*1e9,Tidm);
    nbi=Nmatrl(lam*1e9,Tid);
       
    if (nbim||nbi)<0
       nbim
       nbi
       lam
       Tidm
       Tid
       ddddd
    end
    kim1p=2*pi/lam*nbim;
    hmp=hb(idm);
    wim1p=nbim;
    wip=nbi;
      
    z=z+h;  
    
    Ma=0.5*[exp(-1i*(kim1-ki)*z)*(1+wim1/wi),exp(1i*(kim1+ki)*z)*(1-wim1/wi);...
           exp(-1i*(kim1+ki)*z)*(1-wim1/wi),exp(1i*(kim1-ki)*z)*(1+wim1/wi)];
    M=Ma*M;
    
    Ma2=0.5*[exp(-1i*kim1p*hmp)*(1+wim1p/wip),exp(1i*kim1p*hmp)*(1-wim1p/wip);...
           exp(-1i*kim1p*hmp)*(1-wim1p/wip),exp(1i*kim1p*hmp)*(1+wim1p/wip)];
    M2=Ma2*M2;
    
       
end
nouter=Nmatrl(lam*1e9,nbT{1});
nsubst=Nmatrl(lam*1e9,nbT{end});
nbT{1};
nbT{end};

Emin=-1*M(2,1)/M(2,2);
Epout=(M(1,1)-M(1,2)*M(2,1)/M(2,2));
Pmin(il)=Emin*conj(Emin)*n0/n0; % power back
Ppout(il)=Epout*conj(Epout)*n3/n0; % power forward, power that passed

Emin2=-1*M2(2,1)/M2(2,2);
Epout2=(M2(1,1)-M2(1,2)*M2(2,1)/M2(2,2));
Pmin2(il)=Emin2*conj(Emin2)*nouter/nouter; % power back
Ppout2(il)=Epout2*conj(Epout2)*nsubst/nouter; % power forward, power that passed


end
lamm=lamm/1e-9;
figure
plot(lamm, Pmin, lamm, Pmin+Ppout)
figure
plot(lamm, Pmin2, lamm, Pmin2+Ppout2) % plots reflected power

figure
plot(lamm,Pmin2, lamm,Ppout2) % plots passed power

end