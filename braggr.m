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

n0=1; % outer medium index
n1=1.1; %first layer refractive index
n2=2.5; %second layer index
n3=1.5; % sustrate index
h1=0.5e-6;% first layer thickess in um
h2=0.8e-6;% second layer thickess in um
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

str=[0,1,2,1,2,1,2,1,2,1,2,1,3,2,1,2,1,2,1,2,1,2,3];
hb=[0, h1, h2, 0];
nb=[n0, n1, n2, n3];

% str is layer numbering.
% layer number in str correspondes to material index in nb array. 
% for example, str(2)=1. That correspondes to nb(1)=n0.
% The same with layer thicknes: str(2)=1, layer thicknes hb(1)=0.

nbT={T1, T2, T3, T4};
% nbT is layer index of real material with index dispertion effect.
% Layer elemnt witn number 3 have str(3)=2, nbT(2)=T2('SiO2.txt'), hb(2)=h1.

for il=1:length(lamm)
    lam=lamm(il);
    il/length(lamm)
    M2=[1,0;0,1];

for is=2:length(str)
    id = str(is)+1;
    idm = str(is-1)+1;
    Tid = nbT{id};
    Tidm = nbT{idm};
    
%     nbim=nb(idm); # comment this to activate case 2
%     nbi=nb(id);   # comment this to activate case 2
    nbim=Nmatrl(lam*1e9,Tidm);
    nbi=Nmatrl(lam*1e9,Tid);
    
    kim1p=2*pi/lam*nbim;
    hmp=hb(idm);
    wim1p=nbim;
    wip=nbi;
    
    Ma2=0.5*[exp(-1i*kim1p*hmp)*(1+wim1p/wip),exp(1i*kim1p*hmp)*(1-wim1p/wip);...
           exp(-1i*kim1p*hmp)*(1-wim1p/wip),exp(1i*kim1p*hmp)*(1+wim1p/wip)];
    M2=Ma2*M2;
    
       
end
% nouter=nb(1); # comment this to activate case 2
% nsubst=nb(end);   # comment this to activate case 2
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