function [channel, att] = interModalDisp2(antimpTE,antimpTM,gammaTE,gammaTM,TEmodeimp,TMmodeimp,freq,WGlen,Zo)

for fi=1:length(freq)
    TsTE = diag(exp(-1*gammaTE(fi,:)*WGlen));
    chTEmode = TEmodeimp(fi,:)*TsTE;
    TsTM = diag(exp(-1*gammaTM(fi,:)*WGlen));
    chTMmode = TMmodeimp(fi,:)*TsTM;
    channel(fi) = ((2*Zo)/(abs(antimpTM(fi) + Zo + antimpTE(fi))^2))*...
        (sum(chTEmode)+sum(chTMmode));
    if isnan(channel(fi)) == 1
        channel(fi) = 0;
        att(fi) = 0;
    else
        att(fi) = 10*log10(abs(channel(fi)));
    end
end
end


   


