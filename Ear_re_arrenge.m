%異なる耳介を用いているKEMARを同じ耳介に組み直すプログラム

clear
close all

[FileName, PathName] = uigetfile('*.wav', '1個目の音源を選択');
T = strcat(PathName, FileName);
[y1,Fs1] = audioread(T);

[FileName, PathName] = uigetfile('*.wav', '2個目の音源を選択');
T = strcat(PathName, FileName);
[y2,Fs2] = audioread(T);

left = y1(:,1);
righ = y2(:,1);

out = zeros(len(y1),2);

