%�قȂ鎨���p���Ă���KEMAR�𓯂�����ɑg�ݒ����v���O����

clear
close all

[FileName, PathName] = uigetfile('*.wav', '1�ڂ̉�����I��');
T = strcat(PathName, FileName);
[y1,Fs1] = audioread(T);

[FileName, PathName] = uigetfile('*.wav', '2�ڂ̉�����I��');
T = strcat(PathName, FileName);
[y2,Fs2] = audioread(T);

left = y1(:,1);
righ = y2(:,1);

out = zeros(len(y1),2);

