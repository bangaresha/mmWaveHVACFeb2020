distance = [0.81, 1.06, 1.37, 1.72, 2.03, 2.31, 2.61, 2.89, 3.14, 3.47, 3.75, 4.03, 4.29, 4.57, 4.85];
RSSI = [-36, -36, -36, -37, -37, -37, -39, -40, -40, -39, -40, -39, -40, -39, -41];
SNR = [13, 13, 13, 12, 12, 12, 11, 9, 9, 8, 8, 7, 7, 6, 6];
EVM = [-9, -10, -7, -7, -8, -7, -5, -7, -5, -5, -4, -5, -5, -5, -4];
BitRate = [640.7, 640.3, 639.7, 639.8, 638.9, 638.6, 638.9, 638, 638, 637, 637, 635, 634, 634, 632.8];

figure
plot(distance, RSSI);
title('RSSI versus distance')
yticks('auto')
xticks([0.81, 1.06, 1.37, 1.72, 2.03, 2.31, 2.61, 2.89, 3.14, 3.47, 3.75, 4.03, 4.29, 4.57, 4.85])
xlim([0.81 4.85])

figure
plot(distance, SNR);
title('SNR versus distance')
yticks('auto')
xticks([0.81, 1.06, 1.37, 1.72, 2.03, 2.31, 2.61, 2.89, 3.14, 3.47, 3.75, 4.03, 4.29, 4.57, 4.85])
xlim([0.81 4.85])

figure
plot(distance, EVM);
title('EVM versus distance')
yticks('auto')
xticks([0.81, 1.06, 1.37, 1.72, 2.03, 2.31, 2.61, 2.89, 3.14, 3.47, 3.75, 4.03, 4.29, 4.57, 4.85])
xlim([0.81 4.85])

figure
plot(distance, BitRate);
title('Bit Rate versus distance')
yticks('auto')
xticks([0.81, 1.06, 1.37, 1.72, 2.03, 2.31, 2.61, 2.89, 3.14, 3.47, 3.75, 4.03, 4.29, 4.57, 4.85])
xlim([0.81 4.85])

% figure
% subplot(2,2,1)
% plot(distance, RSSI);
% title('RSSI versus distance')
% subplot(2,2,2)
% plot(distance, SNR);
% title('SNR versus distance')
% subplot(2,2,3)
% plot(distance, EVM);
% title('EVM versus distance')
% subplot(2,2,4)
% plot(distance, BitRate);
% title('Bit Rate versus distance')

%Coefficients (with 95% confidence bounds):
        a1 =      0.0823;%  (0.07829, 0.08632)
       b1 =   3.816e-10;%  (fixed at bound)
       c1 =       47.02;%  (46.96, 47.08)
       a2 =     0.02753;%  (0.02321, 0.03185)
       b2 =   5.175e-09;%  (fixed at bound)
       c2 =      -67.59;%  (-67.75, -67.43)
       freq = 59E9:4E6:60E9;
       omega = 2*pi*freq;
       f=[];
       for x = omega
           %f =  [f a1*sin(b1*x+c1)];
           f =  [f a1*sin(b1*x+c1) + a2*sin(b2*x+c2)];
       end
       h = ifft(f);
       figure
       plot(abs(f));
       figure
       plot(abs(h));
       %188.46 -0.0121   1850.82 0.7751 0.0213 0.00919