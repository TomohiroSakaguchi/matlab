clear
%close all
%<include>prewhiten.m</include>

[FileName, PathName] = uigetfile('*.wav', '�c�������������f�[�^��I��');

T = strcat(PathName, FileName);
[signal,Fs] = audioread(T);%�f�[�^�̓ǂݍ���
signal = signal(1:30*Fs,1);%����̕������������o��
len_data = length(signal);
N = 300; %�P��:ms
windowLen = 2^8;         % window length
shiftLen  = windowLen/8; % shifting stepsize
fftLen    = windowLen;   % number of FFT points
rotateFlag    = true; % deciding phase convention of DGT
zeroPhaseFlag = true; % deciding phase convention of window
%% ����c���̗\���̏���
tau = 70;%70~120ms�͈̔͂Œ���(�����c�����ԕ�)
tau_sample = fix((tau/1000)*Fs);
N_sample = fix((N/1000)*Fs);
lp_range = tau_sample:(len_data-1);
%% ���`�\��
%data_w = prewhiten(signal(:,1));%���O���F��
pcc = lpc(signal(lp_range,1),N_sample);%����c���\���t�B���^�[�̌v�Z
LRP = 0.7.*(filter(pcc,1,signal(:,1)));%���f�[�^�ɑ΂��Č���c���t�B���^�[�������s���C����c���𓱏o
%% DGT�̏���
[win,diffWin] = generalizedCosWin(windowLen,'nuttall4termC1'); %���֐��̐���
signal = zeroPaddingForDGT(signal,shiftLen,fftLen); % DGT�����s������O�Ƀ[������
%% DGT(�X�y�N�g���O�����̌v�Z)
spec = DGT(signal,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);%���M���̃X�y�N�g���O�����̌v�Z
LRPspec = DGT(LRP,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);%����c���̃X�y�N�g���O�����̌v�Z
diffSpec = DGT(signal,diffWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);% �u�����g���̌v�Z������ۂɕK�v
figure, imagesc(20*log10(abs(spec))), colorbar, axis xy%
figure, imagesc(20*log10(abs(diffSpec))), colorbar, axis xy%

%% ����c������
X_abso = abs(spec).^2; 
R_abso = abs(LRPspec).^2;
S_hat_abso_square = X_abso - R_abso;
size_dif = size(S_hat_abso_square);
dif_line = size_dif(1);
dif_row = size_dif(2);
for i = 1:dif_line
    for j = 1:dif_row
        if S_hat_abso_square(i,j) < 0
            S_hat_abso_square(i,j) = 0;
        end
    end
end
S_hat_abso = sqrt(S_hat_abso_square);
figure, imagesc(20*log10(abs(S_hat_abso))), colorbar, axis xy

%% �ʑ��C��
IF = calcInstFreq(spec,diffSpec,fftLen,windowLen,rotateFlag);% �u�����g�����v�Z����֐��C���M����DGT�C��������DGT�������̂��g��
iPCspec = instPhaseCorrection(S_hat_abso,IF,shiftLen,fftLen);% �u���ʑ��C���σX�y�N�g���O�������v�Z����֐�

figure, cRange = max(abs(spec(:)))/5;
subplot(1,3,1), imagesc(abs(spec)), axis xy, caxis(cRange*[-1 1])
title 'Amplitude'
subplot(1,3,2), imagesc(real(spec)), axis xy, caxis(cRange*[-1 1])
title 'Real part'
subplot(1,3,3), imagesc(real(iPCspec)), axis xy, caxis(cRange*[-1 1])
title 'Real part (iPC)'
spec = invInstPhaseCorrection(iPCspec,IF,shiftLen,fftLen); % �ʑ��t��]
figure, imagesc(20*log10(abs(spec))), colorbar, axis xy
%% ���ډ��č\��
dualWin = calcCanonicalDualWindow(win,shiftLen);% 
reconst = invDGT(spec,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);% �ʑ��C���t�ϊ��������̂��tDGT��
