function bgrateM
clc;clear all;close all;
% The script works in two cases:
%1 - a layer parameters that doesn't depend from light wavelength. For that case
%comment lines:
%     nbim=Nmatrl(lam*1e9,Tidm);
%     nbi=Nmatrl(lam*1e9,Tid);
%and
%     nouter=Nmatrl(lam*1e9,nbT{1});
%     nsubst=Nmatrl(lam*1e9,nbT{end});

%2 - a layer parameters depend from light wavelength. For that case
%comment lines:
%     nbim=nb(idm);
%     nbi=nb(id);
%and
%     nouter=nb(1); 
%     nsubst=nb(end);
% that code with data files can be downloaded from:
% https://github.com/ArkadiyMastin/Bragg-grating-optical-coating-simulation-

n0=1; % outer medium index
n1=1.1; %first layer refractive index
n2=2.5; %second layer index
n3=1.5; % sustrate index
lamm=[400:1:1000]*1e-9; % wavelength range in nm

T4 = readtable('N-BK7.txt');
T1 = readtable('Air.txt');
T2 = readtable('SiO2.txt');
%T3 = readtable('Al2O3.txt');
T3 = readtable('TiO2.txt');

function n = Nmatrl(lam, T)
    Tlam=T{:,1};
    Tn=T{:,2};
    n=interp1(Tlam, Tn, lam, 'PCHIP');
end

hb = [0, 94.024, 121.971, 29.191, 14.263, 0]*1e-9;
%coating layer thickness. 94.024 - first layer, ..., 14.263 - last closest
%to substrate layer. So here is 4 layer.
nb=[n0, n1, n2, n1, n2, n3];
% coating index that doesn't depends for wavelength. The right
% element is substrate (n3), left element is incident medium (n0).

nbT={T1, T2, T3, T2, T3, T4};
% coating index that depends for wavelength. The right
% element is substrate (T43), left element is incident medium (T1).
% [Air (incident medium) ,SiO2, TiO2, SiO2, TiO2, N-BK7 (substrate) ]


for il=1:length(lamm)
    lam=lamm(il);
    il/length(lamm)
    M2=[1,0;0,1];

for is=2:length(hb)
    Tid = nbT{is};
    Tidm = nbT{is-1};
    
%     nbim=nb(idm); % comment this to activate case 2
%     nbi=nb(id);   % comment this to activate case 2
    nbim=Nmatrl(lam*1e9,Tidm);
    nbi=Nmatrl(lam*1e9,Tid);
    
    kim1p=2*pi/lam*nbim;
    hmp=hb(is-1);
    wim1p=nbim;
    wip=nbi;
    
    Ma2=0.5*[exp(-1i*kim1p*hmp)*(1+wim1p/wip),exp(1i*kim1p*hmp)*(1-wim1p/wip);...
           exp(-1i*kim1p*hmp)*(1-wim1p/wip),exp(1i*kim1p*hmp)*(1+wim1p/wip)];
    M2=Ma2*M2;
    
       
end
%nouter=nb(1); % comment this to activate case 2
%nsubst=nb(end);   % comment this to activate case 2
nouter=Nmatrl(lam*1e9,nbT{1});
nsubst=Nmatrl(lam*1e9,nbT{end});

Emin2=-1*M2(2,1)/M2(2,2);
Epout2=(M2(1,1)-M2(1,2)*M2(2,1)/M2(2,2));
Pmin2(il)=Emin2*conj(Emin2); % power back
Ppout2(il)=Epout2*conj(Epout2)*nsubst/nouter; % power forward, power that passed


end
lamm=lamm/1e-9;
figure
plot(lamm, Pmin2) % plots reflected power
hold on
plot(lamm, Pmin2+Ppout2)

xlabel('Wavelength, nm')
grid on
legend('Reflection back','Sum power forth and back')

figure
plot(lamm,Pmin2) % plots passed power
hold on
plot(lamm,Ppout2)
xlabel('Wavelength, nm')
legend('Reflection back','Transmission')
grid on
end