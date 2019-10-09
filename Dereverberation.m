clear
%close all
%<include>prewhiten.m</include>

[FileName, PathName] = uigetfile('*.wav', '�c�������������f�[�^��I��');

T = strcat(PathName, FileName);
[signal,Fs] = audioread(T);%�f�[�^�̓ǂݍ���

signal = signal(1:33*Fs,:);%����̕������������o��
len_data = length(signal);
signal_L = signal(:,1);
signal_R = signal(:,2);
%%
windowLen = 2^8;         % window length
shiftLen  = windowLen/8; % shifting stepsize
fftLen    = windowLen;   % number of FFT points
rotateFlag    = true; % deciding phase convention of DGT
zeroPhaseFlag = true; % deciding phase convention of window
%% ����c���̗\���̏���
tau = 70;%70~120ms�͈̔͂Œ���(�����c�����ԕ�)
N = tau/2;
tau_sample = fix((tau/1000)*Fs);
N_sample = fix((N/1000)*Fs);
lp_range = tau_sample:(len_data-1);
%% ���`�\��(�������獶�E��������)
for i=1:2
signal_w = prewhiten(signal(:,i));%���O���F��
pcc = lpc(signal_w(lp_range,1),N_sample);%����c���\���t�B���^�[�̌v�Z
LRP = 0.7.*(filter([0 -pcc(2:end)],1,signal(:,i)));%���f�[�^�ɑ΂��Č���c���t�B���^�[�������s���C����c���𓱏o
%% DGT�̏���
[win,diffWin] = generalizedCosWin(windowLen,'nuttall4termC1'); %���֐��̐���
signal_pad = zeroPaddingForDGT(signal(:,i),shiftLen,fftLen); % DGT�����s������O�Ƀ[�����߁C������ւ񂩂烂�m�����ŏ���
LRP = zeroPaddingForDGT(LRP,shiftLen,fftLen);
%% DGT(�X�y�N�g���O�����̌v�Z)
spec = DGT(signal_pad,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);%���M���̃X�y�N�g���O�����̌v�Z
LRPspec = DGT(LRP,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);%����c���̃X�y�N�g���O�����̌v�Z
diffSpec = DGT(signal_pad,diffWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);% �u�����g���̌v�Z������ۂɕK�v
%% �}�̕\���̎d���͂�����ƍl����
% figure, imagesc(20*log10(abs(spec))), colorbar, axis xy% figure,
% imagesc(20*log10(abs(LRPspec))), colorbar, axis xy%

%% ����c������
X_abso = abs(spec).^2; 
R_abso = abs(LRPspec).^2;
S_hat_abso_square = X_abso - R_abso;
size_dif = size(S_hat_abso_square);
dif_line = size_dif(1);
dif_row = size_dif(2);
for m = 1:dif_line
    for n = 1:dif_row
        if S_hat_abso_square(m,n) < 0
            S_hat_abso_square(m,n) = 0;
        end
    end
end
S_hat_abso = sqrt(S_hat_abso_square);
figure, imagesc(20*log10(abs(S_hat_abso))), colorbar, axis xy
%% �ʑ��̌v�Z
Fai = spec./S_hat_abso;%���f�ʑ��̎Z�o
S_hat = S_hat_abso.*Fai; %S_hat = S_hat_abso * exp(i*Fai);

for m = 1:dif_line
    for n = 1:dif_row
        if isnan(S_hat(m,n))
            S_hat(m,n) = 0;
        end
    end
end
%% �ʑ��C��
IF = calcInstFreq(S_hat,diffSpec,fftLen,windowLen,rotateFlag);% �u�����g�����v�Z����֐��C���M����DGT�C��������DGT�������̂��g��
iPCspec = instPhaseCorrection(S_hat,IF,shiftLen,fftLen);% �u���ʑ��C���σX�y�N�g���O�������v�Z����֐�

figure, cRange = max(abs(spec(:)))/5;
subplot(1,3,1), imagesc(abs(spec(1:40,fix(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Amplitude'
subplot(1,3,2), imagesc(real(spec(1:40,fix(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Real part'
subplot(1,3,3), imagesc(real(iPCspec(1:40,fix(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Real part (iPC)'
recon_spec = invInstPhaseCorrection(iPCspec,IF,shiftLen,fftLen); % �ʑ��t��]
figure, imagesc(20*log10(abs(recon_spec))), colorbar, axis xy
%% ���ډ��č\��
dualWin = calcCanonicalDualWindow(win,shiftLen);% 
reconst = invDGT(recon_spec,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);% �ʑ��C���t�ϊ��������̂��tDGT��

%% ���E����
if i == 1
    len_reconst = length(reconst);
    len_rev = length(LRP);
    direct = zeros(len_reconst,2);
    rev = zeros(len_rev,2);    
end
direct(:,i) = reconst;
rev(:,i) = LRP;
end
%% �I�[�f�B�I�����o��
audiowrite('direct.wav',direct,Fs);
audiowrite('diffuse.wav',rev,Fs);
audiowrite('original.wav',signal,Fs);
