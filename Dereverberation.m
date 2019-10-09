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
%% ���`�\��
signal_w_L = prewhiten(signal_L);%���O���F��
pcc = lpc(signal_w_L(lp_range,1),N_sample);%����c���\���t�B���^�[�̌v�Z
LRP_L = 0.7.*(filter([0 -pcc(2:end)],1,signal_L));%���f�[�^�ɑ΂��Č���c���t�B���^�[�������s���C����c���𓱏o
signal_w_R = prewhiten(signal_R);%���O���F��
pcc = lpc(signal_w_R(lp_range,1),N_sample);%����c���\���t�B���^�[�̌v�Z
LRP_R = 0.7.*(filter([0 -pcc(2:end)],1,signal_R));%���f�[�^�ɑ΂��Č���c���t�B���^�[�������s���C����c���𓱏o
%% DGT�̏���(�������獶�E��������)�܂�����ĂȂ�
[win,diffWin] = generalizedCosWin(windowLen,'nuttall4termC1'); %���֐��̐���
signal_L = zeroPaddingForDGT(signal_L,shiftLen,fftLen); % DGT�����s������O�Ƀ[�����߁C������ւ񂩂烂�m�����ŏ���
signal_R = zeroPaddingForDGT(signal_R,shiftLen,fftLen); % DGT�����s������O�Ƀ[�����߁C������ւ񂩂烂�m�����ŏ���
LRP_L = zeroPaddingForDGT(LRP_L,shiftLen,fftLen);
LRP_R = zeroPaddingForDGT(LRP_R,shiftLen,fftLen);
%% DGT(�X�y�N�g���O�����̌v�Z)
specL = DGT(signal_L,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);%���M���̃X�y�N�g���O�����̌v�Z
specR = DGT(signal_R,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);%���M���̃X�y�N�g���O�����̌v�Z
LRPspecL = DGT(LRP_L,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);%����c���̃X�y�N�g���O�����̌v�Z
LRPspecR = DGT(LRP_R,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);%����c���̃X�y�N�g���O�����̌v�Z
diffSpecL = DGT(signal_L,diffWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);% �u�����g���̌v�Z������ۂɕK�v
diffSpecR = DGT(signal_R,diffWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);% �u�����g���̌v�Z������ۂɕK�v

%% �}�̕\���̎d���͂�����ƍl����
% figure, imagesc(20*log10(abs(spec))), colorbar, axis xy% figure,
% imagesc(20*log10(abs(LRPspec))), colorbar, axis xy%

%% ����c������
% Lch
X_L_abso = abs(specL).^2; 
R_L_abso = abs(LRPspecL).^2;
S_hat_L_abso_square = X_L_abso - R_L_abso;
size_dif = size(S_hat_L_abso_square);
dif_line = size_dif(1);
dif_row = size_dif(2);
for i = 1:dif_line
    for j = 1:dif_row
        if S_hat_L_abso_square(i,j) < 0
            S_hat_L_abso_square(i,j) = 0;
        end
    end
end
S_hat_L_abso = sqrt(S_hat_L_abso_square);
figure, imagesc(20*log10(abs(S_hat_L_abso))), colorbar, axis xy

% ��������Rch
X_R_abso = abs(specR).^2; 
R_R_abso = abs(LRPspecR).^2;
S_hat_R_abso_square = X_R_abso - R_R_abso;
size_dif = size(S_hat_R_abso_square);
dif_line = size_dif(1);
dif_row = size_dif(2);
for i = 1:dif_line
    for j = 1:dif_row
        if S_hat_R_abso_square(i,j) < 0
            S_hat_R_abso_square(i,j) = 0;
        end
    end
end
S_hat_R_abso = sqrt(S_hat_R_abso_square);
figure, imagesc(20*log10(abs(S_hat_R_abso))), colorbar, axis xy
%% �ʑ��̌v�Z
% Lch
Fai_L = specL./S_hat_L_abso;%���f�ʑ��̎Z�o
S_hat_L = S_hat_L_abso.*Fai_L; %S_hat = S_hat_abso * exp(i*Fai);

for i = 1:dif_line
    for j = 1:dif_row
        if isnan(S_hat_L(i,j))
            S_hat_L(i,j) = 0;
        end
    end
end

% Rch
Fai_R = specR./S_hat_R_abso;%���f�ʑ��̎Z�o
S_hat_R = S_hat_R_abso.*Fai_R; %S_hat = S_hat_abso * exp(i*Fai);

for i = 1:dif_line
    for j = 1:dif_row
        if isnan(S_hat_R(i,j))
            S_hat_R(i,j) = 0;
        end
    end
end
%% �ʑ��C��
% Lch
IF_L = calcInstFreq(S_hat_L,diffSpecL,fftLen,windowLen,rotateFlag);% �u�����g�����v�Z����֐��C���M����DGT�C��������DGT�������̂��g��
iPCspec_L = instPhaseCorrection(S_hat_L,IF_L,shiftLen,fftLen);% �u���ʑ��C���σX�y�N�g���O�������v�Z����֐�
IFrev_L = calcInstFreq(LRPspecL,diffSpecL,fftLen,windowLen,rotateFlag);
iPCrev_L = instPhaseCorrection(LRPspecL,IFrev_L,shiftLen,fftLen);

figure, cRange = max(abs(specL(:)))/5;
subplot(1,3,1), imagesc(abs(specL(1:40,(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Amplitude'
subplot(1,3,2), imagesc(real(specL(1:40,(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Real part'
subplot(1,3,3), imagesc(real(iPCspec_L(1:40,(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Real part (iPC)'
recon_spec_L = invInstPhaseCorrection(iPCspec_L,IF_L,shiftLen,fftLen); % �ʑ��t��]
recon_rev_L = invInstPhaseCorrection(iPCrev_L,IF_L,shiftLen,fftLen);
figure, imagesc(20*log10(abs(recon_spec_L))), colorbar, axis xy

% Rch
IF_R = calcInstFreq(S_hat_R,diffSpecR,fftLen,windowLen,rotateFlag);% �u�����g�����v�Z����֐��C���M����DGT�C��������DGT�������̂��g��
iPCspec_R = instPhaseCorrection(S_hat_R,IF_R,shiftLen,fftLen);% �u���ʑ��C���σX�y�N�g���O�������v�Z����֐�
IFrev_R = calcInstFreq(LRPspecR,diffSpecR,fftLen,windowLen,rotateFlag);
iPCrev_R = instPhaseCorrection(LRPspecR,IFrev_R,shiftLen,fftLen);

figure, cRange = max(abs(specR(:)))/5;
subplot(1,3,1), imagesc(abs(specR(1:40,fix(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Amplitude'
subplot(1,3,2), imagesc(real(specR(1:40,fix(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Real part'
subplot(1,3,3), imagesc(real(iPCspec_R(1:40,fix(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Real part (iPC)'
recon_spec_R = invInstPhaseCorrection(iPCspec_R,IF_R,shiftLen,fftLen); % �ʑ��t��]
recon_rev_R = invInstPhaseCorrection(iPCrev_R,IF_R,shiftLen,fftLen);
figure, imagesc(20*log10(abs(recon_spec_R))), colorbar, axis xy
%% ���ډ��č\��
dualWin = calcCanonicalDualWindow(win,shiftLen);% 
reconst_L = invDGT(recon_spec_L,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);% �ʑ��C���t�ϊ��������̂��tDGT��
reconstrev_L = invDGT(recon_rev_L,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);
reconst_R = invDGT(recon_spec_R,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);% �ʑ��C���t�ϊ��������̂��tDGT��
reconstrev_R = invDGT(recon_rev_R,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);

%% ���E����
direct = [reconst_L,reconst_R];
rev = [reconstrev_L,reconstrev_R];

%% �I�[�f�B�I�����o��
audiowrite('direct.wav',direct,Fs);
audiowrite('diffuse.wav',rev,Fs);
audiowrite('original.wav',signal,Fs);
