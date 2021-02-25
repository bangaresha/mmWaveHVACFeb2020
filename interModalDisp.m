function [channel, att, totPowerS, theta, meanDelay, rmsDelay, pn_dB] = interModalDisp(antimpTE,antimpTM,gammaTE,gammaTM,TEmodeimp,TMmodeimp,freq,WGlen,Zo)

for fi=1:length(freq)
    TsTE = diag(exp(-1*gammaTE(fi,:)*WGlen));
    chTEmode = TEmodeimp(fi,:)*TsTE;
    TsTM = diag(exp(-1*gammaTM(fi,:)*WGlen));
    chTMmode = TMmodeimp(fi,:)*TsTM;
    totPowerS(fi) = 2*(sum(sum(TsTE)) + sum(sum(TsTM)));
    channel(fi) = ((2*Zo)/(abs(antimpTM(fi) + Zo + antimpTE(fi))^2))*...
        (sum(chTEmode)+sum(chTMmode));
    theta(fi) = angle(channel(fi));

    if isnan(channel(fi)) == 1
        channel(fi) = 0;
        att(fi) = 0;
    else
        att(fi) = 10*log10(abs(channel(fi)));
    end
end

channelTime = ifft(channel,200000);

temp = 0;
for i = 1:length(channelTime)
    pn_dB(i) = (10*log10((abs(channelTime(i)))^2)) + 70;
    pn(i) = (abs((channelTime(i)))^2)*1E4;
    temp = temp + (pn(i))*i;
end
maxpn = max(pn_dB);

pndBtemp1=[];
delayTime=[];
for i = 1:length(pn_dB)
    pndBtemp = pn_dB(i) - maxpn;
    if pndBtemp >= -20
        pndBtemp1 = [pndBtemp1 pndBtemp];
        delayTime = [delayTime i];
    end
end
delayTime = delayTime*(3.2E-12)/200000;
meanDelay = sum(delayTime.*pndBtemp1)/sum(pndBtemp1);
meanDelaySq = sum(delayTime.*delayTime.*pndBtemp1)/sum(pndBtemp1);
rmsDelay = 200000*sqrt(meanDelaySq - meanDelay*meanDelay);
end


   


