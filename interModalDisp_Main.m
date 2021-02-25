clear;
tic
% close all;
c=3*10^8;
% radius=0.305/2;
radius=0.127/2;
% freq = 2.40E9:0.3125E6:2.5E9;
freq = 59.4E9:1E7:60.4E9;

load('besDerZerMat.mat');
load('besZerMat.mat');

nTE=[];
mTE=[];
fcTE=[];
coWnTE=[];
nTM=[];
mTM=[];
fcTM=[];
coWnTM=[];
for m=1:1001
    for n=1:1000
        tt = besDerZerMat(m,n);
        fc_TE_temp=(c/(2*pi*radius))*besDerZerMat(m,n);
        if fc_TE_temp <= freq(length(freq)) %freq(fi)
            nTE = [nTE; n];
            mTE = [mTE; m-1];
            fcTE = [fcTE; fc_TE_temp];
            coWnTE = [coWnTE; besDerZerMat(m,n)/radius];
        end
        fc_TM_temp=(c/(2*pi*radius))*besZerMat(m,n);
        if fc_TM_temp <= freq(length(freq))
            nTM = [nTM; n];
            mTM = [mTM; m-1];
            fcTM = [fcTM; fc_TM_temp];
            coWnTM = [coWnTM; (besZerMat(m,n))/radius];
        end
    end
end
    
Zo=50;
eta = 377;
% l = 0.035;
l = 0.001;
sigma = 10^6;
mu = 4*pi*10^-7;
R = ((2*pi*freq*mu)./(2*sigma)).^0.5;
k = 2*pi*freq./c;
% WGlenS = 1:1:10;
WGlenS = 4;
channel = [];
pn_dB = [];
att = [];
cdfp = [];
meanCDF = [];
varCDF = [];

[radresTE, radresTM, gammaTE, gammaTM] = radResCyl_multitone(mTE,nTE,mTM,nTM,radius,freq,fcTE,fcTM,c,k,R,eta,l,coWnTE,coWnTM);

radreacTE = imag(hilbert(radresTE));
radreacTM = imag(hilbert(radresTM));

for fi=1:length(freq)
    antresTE(fi) = sum(radresTE(fi,:));
    antresTM(fi) = sum(radresTM(fi,:));
    antreacTE(fi) = sum(radreacTE(fi,:));
    antreacTM(fi) = sum(radreacTM(fi,:));
    antimpTE(fi) = antresTE(fi) +1i*antreacTE(fi);
    antimpTM(fi) = antresTM(fi) +1i*antreacTM(fi);
    for n = 1:length(nTE)
        TEmodeimp(fi,n)=radresTE(fi,n)+1i*radreacTE(fi,n);
    end
    for n = 1:length(nTM)
        TMmodeimp(fi,n)=radresTM(fi,n)+1i*radreacTM(fi,n);
    end    
end

for i = 1:length(WGlenS)
%     [channelT, attT, totPowerS, theta, meanDelayT, rmsDelayT, pn_dBT] = totalDisp(freq,fc_TE,fc_TM,c,m_TE,n_TE,m_TM,n_TM,k,WGlenS(i),radius,Zo,R,eta,l,coWnTE,coWnTM);
    [channel, att, totPowerS, theta, meanDelay, rmsDelay, pn_dB] = interModalDisp(antimpTE,antimpTM,gammaTE,gammaTM,TEmodeimp,TMmodeimp,freq,WGlenS(i),Zo);
%     channel = [channel; channelT];
%     att = [att; attT];
%     meanDelay(i) = meanDelayT;
%     rmsDelay(i) = rmsDelayT;
%     pn_dB = [pn_dB; pn_dBT];
%     attB = attT(attT ~= 0);
%     meanCDF = [meanCDF mean(attB)];
%     varCDF = [varCDF var(attB)];
end
figure
plot(freq,att(1,:),'r');

delayax = 0.001:0.016:3200;
pdf = pn_dB+50;

figure
plot(delayax, pdf);
title('PDF versus Time');

cdfxx = -36:0.5:-31.5;
pd = fitdist(meanCDF','Normal');
figure
cdfplot(meanCDF)
hold on
cdfpT = cdf(pd,cdfxx);
plot(cdfxx,cdfpT);

figure
pdfPT = pdf(pd,cdfxx);
plot(cdfxx,pdfPT,'LineWidth',2)
    
cdfVarx = 1.75:0.25:4;
pdVar = fitdist(varCDF','Normal');
figure
cdfplot(varCDF)
hold on
cdfVarpT = cdf(pdVar,cdfVarx);
plot(cdfVarx,cdfVarpT);

chAtVsDis = acons.*exp(bcons.*WGlenS);
   
figure
plot(WGlenS,meanCDF);
hold on
plot(WGlenS,chAtVsDis);
    
    


% ylim([-45,0]);
% xlim([0.001,2000]);

% figure
% plot(delayax, [pdf(1:150000) -60*ones(1,50000)]);
% title('PDF versus Time');
% ylim([-45,0]);
% xlim([0.001,2000]);

figure
plot(WGlenS, rmsDelay);
title('RMS Delay versus Time');

% figure
% plot(refCoff1n, rmsTemp);
% title('RMS Delay versus Reflection Coefficient');


% figure
% plot(WGlenS, rmstemp);
% title('RMS Delay versus Time');

% cohBW = 1E9./(5*rmsDelay);

cohBWtemp = 1./(5*rmsDelay);

figure
plot(WGlenS, cohBWtemp);
title('Coherence BW versus Time');

%pd = fitdist(B','Normal');
% ecdf(B,'Bounds','on')
% hold off
%pd = makedist('Normal');
%x = -3:.1:0;
toc
