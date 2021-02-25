tic
% close all;
c=3*10^8;
radius=0.305/2;
% radius=0.127/2;
freq = 2.40E9:0.3125E6:2.5E9;
% freq = 59E9:2E6:60E9;
% freq = 59.4E9:2E7:60.4E9;

load('besDerZerMat.mat');
load('besZerMat.mat');

for fi=1:length(freq)
    n_TE=[];
    m_TE=[];
    fc_TE=[];
    coWnTE=[];
    for m=1:1001
        for n=1:1000
            tt = besDerZerMat(m,n);
            fc_TE_temp=(c/(2*pi*radius))*besDerZerMat(m,n);
            if fc_TE_temp <= freq(fi)
                n_TE = [n_TE; n];
                m_TE = [m_TE; m-1];
                fc_TE = [fc_TE; fc_TE_temp];
                coWnTE = [coWnTE; besDerZerMat(m,n)/radius];
            end
        end
    end
end

for fi=1:length(freq)
    n_TM=[];
    m_TM=[];
    fc_TM=[];
    coWnTM=[];
    for m=1:1001
        for n=1:1000
            fc_TM_temp=(c/(2*pi*radius))*besZerMat(m,n);
            if fc_TM_temp <= freq(fi)
                n_TM = [n_TM; n];
                m_TM = [m_TM; m-1];
                fc_TM = [fc_TM; fc_TM_temp];
                coWnTM = [coWnTM; (besZerMat(m,n))/radius];
            end
        end
    end
end

Zo=50;
eta = 377;
l = 0.035;
% l = 0.005;

sigma = 10^6;
mu = 4*pi*10^-7;
R = ((2*pi*freq*mu)./(2*sigma)).^0.5;
k = 2*pi*freq./c;
% WGlenS = 1:1:10;
WGlenS = 4.57;
channel = [];
att = [];

[radresTE, radresTM, gammaTE, gammaTM] = radResCyl_multitone2(m_TE,n_TE,m_TM,n_TM,radius,freq,fc_TE,fc_TM,c,k,R,eta,l,coWnTE,coWnTM);

radreacTE = imag(hilbert(radresTE));
radreacTM = imag(hilbert(radresTM));

for fi=1:length(freq)
    antresTE(fi) = sum(radresTE(fi,:));
    antresTM(fi) = sum(radresTM(fi,:));
    antreacTE(fi) = sum(radreacTE(fi,:));
    antreacTM(fi) = sum(radreacTM(fi,:));
    antimpTE(fi) = antresTE(fi) +1i*antreacTE(fi);
    antimpTM(fi) = antresTM(fi) +1i*antreacTM(fi);
    for n = 1:length(n_TE)
        TEmodeimp(fi,n)=radresTE(fi,n)+1i*radreacTE(fi,n);
    end
    for n = 1:length(n_TM)
        TMmodeimp(fi,n)=radresTM(fi,n)+1i*radreacTM(fi,n);
    end    
end

for i = 1:length(WGlenS)
    [channelT, attT] = interModalDisp2(antimpTE,antimpTM,gammaTE,gammaTM,TEmodeimp,TMmodeimp,freq,WGlenS(i),Zo);
    channel = [channel; channelT];
    att = [att; attT];
    attB = attT(attT ~= 0);
end

figure
plot(freq, att);
title('Coherence BW versus Time');
toc
