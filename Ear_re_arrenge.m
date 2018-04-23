%異なる耳介を用いているKEMARを同じ耳介に組み直すプログラム

clear
close all

[FileName, PathName] = uigetfile('*.wav', '1個目の音源を選択');
T = strcat(PathName, FileName);
[y1,Fs1] = audioread(T);

[FileName, PathName] = uigetfile('*.wav', '2個目の音源を選択');
T = strcat(PathName, FileName);
[y2,Fs2] = audioread(T);

len_y1 = length(y1);
t_y1 = (len_y1-1)/Fs1;
time1 = 0:1/Fs1:t_y1;

left = y1(:,1);
righ = y2(:,1);

out = zeros(len_y1,2);%組み合わせ用の変数を用意する

out(:,1) = left;
out(:,2) = righ;

figure(1)
subplot(2,1,1)
plot(time1,y1(:,1))

subplot(2,1,2)
plot(time1,y1(:,2))


figure(2)
subplot(2,1,1)
plot(time1,out(:,1))

subplot(2,1,2)
plot(time1,out(:,2))

soundsc(out,Fs1)