function [radresTE, radresTM, gammaTE, gammaTM] = radResCyl_multitone(m_TE,n_TE,m_TM,n_TM,radius,freq,fc_TE,fc_TM,c,k,R,eta,l,coWnTE,coWnTM)
gammaTE = zeros(length(freq), length(n_TE));
gammaTM = zeros(length(freq), length(n_TM));
radresTE = zeros(length(freq), length(n_TE));
radresTM = zeros(length(freq), length(n_TM));
sumRadResTE = zeros(1,length(freq));
sumRadResTM = zeros(1,length(freq));
 
for fr=1:length(freq)
    for n = 1:length(n_TE)
        gTESq = (coWnTE(n))^2;
        betaTE = (sqrt(k(fr)^2 - gTESq));%*1E-3;
        alphaTempTE = (R(fr)*k(fr))/(eta*radius*betaTE);
        alphaTE=8.85*(alphaTempTE*(((m_TE(n)^2)/...
            (((radius*(coWnTE(n)))^2)- ((m_TE(n))^2)))...
            +(((c*coWnTE(n))/(2*pi*freq(fr)))^2)));
        gammaTE(fr,n)=alphaTE+1i*betaTE;
        constRadTE = eta*(k(fr))*((m_TE(n))^2)/...
            (betaTE*pi*((sin(l*(k(fr))*pi/180))^2)*...
            ((besselj(m_TE(n),radius*coWnTE(n)))^4)*...
            ((((radius*coWnTE(n))^2)-((m_TE(n))^2))^2));
        radTEfunc=@(zeta) ((besselj(m_TE(n),radius.*zeta.*coWnTE(n)).*...
            sin((k(fr).*radius.*((l/radius)-1+zeta))*pi/180))./zeta);
        radTEnum = (integral(radTEfunc,1-(l/radius),1)).^2;
        radresTE(fr,n)= constRadTE*radTEnum;
        if imag(radresTE(fr,n)) ~= 0
            radresTE(fr,n) =0;
            gammaTE(fr,n) = 0;
        end
    end
    for n = 1:length(n_TM)
        gTMSq = (coWnTM(n))^2;
        betaTM = (sqrt(k(fr)^2 - gTMSq));%*1E-3;
        alphaconstTM = (R(fr)*k(fr))/(eta*radius*betaTM);
        alphaTM=8.85*alphaconstTM;
        gammaTM(fr,n)=alphaTM+1i*betaTM;
        constRadTM = eta*betaTM/(pi*((sin((k(fr))*l*pi/180))^2)*...
            (k(fr))*((0.5*((besselj(m_TM(n)-1,radius*coWnTM(n)))-...
            (besselj(m_TM(n)+1,radius*coWnTM(n)))))^2));
        
        radTMfunc=@(zeta) 0.5*((besselj((m_TM(n))-1,zeta.*...
            (coWnTM(n))*radius))-(besselj(m_TM(n)+1,zeta.*...
            (coWnTM(n))*radius))).*...
            sin((k(fr).*radius.*((l./radius)-1+zeta))*pi/180);
        radTMnum =(integral(radTMfunc,1-(l./radius),1)).^2;
        radresTM(fr,n)= constRadTM*radTMnum;
        if imag(radresTM(fr,n)) ~= 0
            radresTM(fr,n) =0;
            gammaTM(fr,n) = 0;
        end
    end
end
pTE_f = 0;
fcTE_f = 0;
mnTE_f = 0;
%nTE_f = 0;
pTE = zeros(length(freq), length(radresTE)); 
pTM = zeros(length(freq), length(radresTM));

for fr=1:length(freq)
    sumRadResTE(fr) = sum(radresTE(fr,:));
    for i = 1:length(radresTE(fr,:))
        pTE(fr,i) = 100*radresTE(fr,i)/sumRadResTE(fr);
        if pTE(fr,i) >= 10
            pTE_f = [pTE_f, pTE(fr,i)];
            fcTE_f = [fcTE_f fc_TE(i)];
            %mTE_f = [mTE_f m_TE(i)];
            mnTE_f = [mnTE_f, strcat('TE',string(m_TE(i)),',' ,string(n_TE(i)))];
            %nTE_f = [nTE_f, strcat('TE',string(n_TE(i)))];
            %nTE_f = [nTE_f n_TE(i)];
        end
    end
end

pTE_f = (pTE_f(2:length(pTE_f)))';
fcTE_f = (fcTE_f(2:length(fcTE_f)))';
mnTE_f = (mnTE_f(2:length(mnTE_f)))';
%nTE_f = (nTE_f(2:length(nTE_f)))';

pTM_f = 0;
fcTM_f = 0;
mnTM_f = 0;
%nTM_f = 0;
for fr=1:length(freq)
    sumRadResTM(fr) = sum(radresTM(fr,:));
    for i = 1:length(radresTM(fr,:))
        pTM(fr,i) = 100*radresTM(fr,i)/sumRadResTM(fr);
        if pTM(fr,i) >= 10
            pTM_f = [pTM_f, pTM(fr,i)];
            fcTM_f = [fcTM_f, fc_TM(i)];
%             mTM_f = [mTM_f, m_TM(i)];
%             nTM_f = [nTM_f, n_TM(i)];
            mnTM_f = [mnTM_f, strcat('TM',string(m_TM(i)), ',', string(n_TM(i)))];
            %nTM_f = [nTM_f, strcat('TM',string(n_TM(i)))];
        end
    end
end

pTM_f = (pTM_f(2:length(pTM_f)))';
fcTM_f = (fcTM_f(2:length(fcTM_f)))';
mnTM_f = (mnTM_f(2:length(mnTM_f)))';

p = [pTE_f', pTM_f'];
fc = [fcTE_f', fcTM_f'];
% sc = stem(fc,p');
% title('Mode Power(%Total Power) Versus Mode Cutoff Frequency');
% modePower = [sum(radRes) counting];
% modePower = [p; counting];
end

