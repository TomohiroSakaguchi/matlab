clear
%close all
%<include>prewhiten.m</include>

[FileName, PathName] = uigetfile('*.wav', '残響除去したいデータを選択');

T = strcat(PathName, FileName);
[signal,Fs] = audioread(T);%データの読み込み

signal = signal(1:33*Fs,:);%特定の部分だけ抜き出し
len_data = length(signal);
signal_L = signal(:,1);
signal_R = signal(:,2);
%%
windowLen = 2^8;         % window length
shiftLen  = windowLen/8; % shifting stepsize
fftLen    = windowLen;   % number of FFT points
rotateFlag    = true; % deciding phase convention of DGT
zeroPhaseFlag = true; % deciding phase convention of window
%% 後期残響の予測の準備
tau = 70;%70~120msの範囲で調整(初期残響時間分)
N = tau/2;
tau_sample = fix((tau/1000)*Fs);
N_sample = fix((N/1000)*Fs);
lp_range = tau_sample:(len_data-1);
%% 線形予測
signal_w_L = prewhiten(signal_L);%事前白色化
pcc = lpc(signal_w_L(lp_range,1),N_sample);%後期残響予測フィルターの計算
LRP_L = 0.7.*(filter([0 -pcc(2:end)],1,signal_L));%元データに対して後期残響フィルターかけを行い，後期残響を導出
signal_w_R = prewhiten(signal_R);%事前白色化
pcc = lpc(signal_w_R(lp_range,1),N_sample);%後期残響予測フィルターの計算
LRP_R = 0.7.*(filter([0 -pcc(2:end)],1,signal_R));%元データに対して後期残響フィルターかけを行い，後期残響を導出
%% DGTの準備(ここから左右分割処理)まだやってない
[win,diffWin] = generalizedCosWin(windowLen,'nuttall4termC1'); %窓関数の生成
signal_L = zeroPaddingForDGT(signal_L,shiftLen,fftLen); % DGTを実行させる前にゼロ埋め，ここらへんからモノラルで処理
signal_R = zeroPaddingForDGT(signal_R,shiftLen,fftLen); % DGTを実行させる前にゼロ埋め，ここらへんからモノラルで処理
LRP_L = zeroPaddingForDGT(LRP_L,shiftLen,fftLen);
LRP_R = zeroPaddingForDGT(LRP_R,shiftLen,fftLen);
%% DGT(スペクトログラムの計算)
specL = DGT(signal_L,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);%元信号のスペクトログラムの計算
specR = DGT(signal_R,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);%元信号のスペクトログラムの計算
LRPspecL = DGT(LRP_L,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);%後期残響のスペクトログラムの計算
LRPspecR = DGT(LRP_R,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);%後期残響のスペクトログラムの計算
diffSpecL = DGT(signal_L,diffWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);% 瞬時周波数の計算をする際に必要
diffSpecR = DGT(signal_R,diffWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);% 瞬時周波数の計算をする際に必要

%% 図の表示の仕方はちょっと考える
% figure, imagesc(20*log10(abs(spec))), colorbar, axis xy% figure,
% imagesc(20*log10(abs(LRPspec))), colorbar, axis xy%

%% 後期残響除去
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

% ここからRch
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
%% 位相の計算
% Lch
Fai_L = specL./S_hat_L_abso;%複素位相の算出
S_hat_L = S_hat_L_abso.*Fai_L; %S_hat = S_hat_abso * exp(i*Fai);

for i = 1:dif_line
    for j = 1:dif_row
        if isnan(S_hat_L(i,j))
            S_hat_L(i,j) = 0;
        end
    end
end

% Rch
Fai_R = specR./S_hat_R_abso;%複素位相の算出
S_hat_R = S_hat_R_abso.*Fai_R; %S_hat = S_hat_abso * exp(i*Fai);

for i = 1:dif_line
    for j = 1:dif_row
        if isnan(S_hat_R(i,j))
            S_hat_R(i,j) = 0;
        end
    end
end
%% 位相修正
% Lch
IF_L = calcInstFreq(S_hat_L,diffSpecL,fftLen,windowLen,rotateFlag);% 瞬時周波数を計算する関数，源信号のDGT，微分窓をDGTしたものを使う
iPCspec_L = instPhaseCorrection(S_hat_L,IF_L,shiftLen,fftLen);% 瞬時位相修正済スペクトログラムを計算する関数
IFrev_L = calcInstFreq(LRPspecL,diffSpecL,fftLen,windowLen,rotateFlag);
iPCrev_L = instPhaseCorrection(LRPspecL,IFrev_L,shiftLen,fftLen);

figure, cRange = max(abs(specL(:)))/5;
subplot(1,3,1), imagesc(abs(specL(1:40,(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Amplitude'
subplot(1,3,2), imagesc(real(specL(1:40,(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Real part'
subplot(1,3,3), imagesc(real(iPCspec_L(1:40,(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Real part (iPC)'
recon_spec_L = invInstPhaseCorrection(iPCspec_L,IF_L,shiftLen,fftLen); % 位相逆回転
recon_rev_L = invInstPhaseCorrection(iPCrev_L,IF_L,shiftLen,fftLen);
figure, imagesc(20*log10(abs(recon_spec_L))), colorbar, axis xy

% Rch
IF_R = calcInstFreq(S_hat_R,diffSpecR,fftLen,windowLen,rotateFlag);% 瞬時周波数を計算する関数，源信号のDGT，微分窓をDGTしたものを使う
iPCspec_R = instPhaseCorrection(S_hat_R,IF_R,shiftLen,fftLen);% 瞬時位相修正済スペクトログラムを計算する関数
IFrev_R = calcInstFreq(LRPspecR,diffSpecR,fftLen,windowLen,rotateFlag);
iPCrev_R = instPhaseCorrection(LRPspecR,IFrev_R,shiftLen,fftLen);

figure, cRange = max(abs(specR(:)))/5;
subplot(1,3,1), imagesc(abs(specR(1:40,fix(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Amplitude'
subplot(1,3,2), imagesc(real(specR(1:40,fix(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Real part'
subplot(1,3,3), imagesc(real(iPCspec_R(1:40,fix(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Real part (iPC)'
recon_spec_R = invInstPhaseCorrection(iPCspec_R,IF_R,shiftLen,fftLen); % 位相逆回転
recon_rev_R = invInstPhaseCorrection(iPCrev_R,IF_R,shiftLen,fftLen);
figure, imagesc(20*log10(abs(recon_spec_R))), colorbar, axis xy
%% 直接音再構成
dualWin = calcCanonicalDualWindow(win,shiftLen);% 
reconst_L = invDGT(recon_spec_L,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);% 位相修正逆変換したものを逆DGTし
reconstrev_L = invDGT(recon_rev_L,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);
reconst_R = invDGT(recon_spec_R,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);% 位相修正逆変換したものを逆DGTし
reconstrev_R = invDGT(recon_rev_R,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);

%% 左右結合
direct = [reconst_L,reconst_R];
rev = [reconstrev_L,reconstrev_R];

%% オーディオ書き出し
audiowrite('direct.wav',direct,Fs);
audiowrite('diffuse.wav',rev,Fs);
audiowrite('original.wav',signal,Fs);
