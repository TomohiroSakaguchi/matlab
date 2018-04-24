%ŒQ’x‰„‚ÆˆÊ‘Š·‚Å’x‰„ŠÔ‚ğZo‚·‚éƒvƒƒOƒ‰ƒ€
clear
close all

[FileName, PathName] = uigetfile('*.wav', '‰¹Œ¹‚ğ‘I‘ğ');
T = strcat(PathName, FileName);
[y,Fs] = audioread(T);

%gd = grpdelay(y,512);
[h,w] = freqz(y,512);
pd = -unwrap(angle(h))./w;

%plot(w/pi,gd,w/pi,pd),grid
plot(w/pi,gd),grid
xlabel 'Normalized Frequency (\times\pi rad/sample)'
ylabel 'Group and Phase delays'

%legend('Group delay','Phase delay')
legend('Phase delay')