%�قȂ鎨���p���Ă���KEMAR�𓯂�����ɑg�ݒ����v���O����

clear
close all

[FileName, PathName] = uigetfile('*.wav', '1�ڂ̉�����I��');
T = strcat(PathName, FileName);
[y1,Fs1] = audioread(T);

[FileName, PathName] = uigetfile('*.wav', '2�ڂ̉�����I��');
T = strcat(PathName, FileName);
[y2,Fs2] = audioread(T);

len_y1 = length(y1);
t_y1 = (len_y1-1)/Fs1;
time1 = 0:1/Fs1:t_y1;

left = y1(:,1);
righ = y2(:,1);

out = zeros(len_y1,2);%�g�ݍ��킹�p�̕ϐ���p�ӂ���

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