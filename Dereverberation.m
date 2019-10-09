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
%% 線形予測(ここから左右分割処理)
for i=1:2
signal_w = prewhiten(signal(:,i));%事前白色化
pcc = lpc(signal_w(lp_range,1),N_sample);%後期残響予測フィルターの計算
LRP = 0.7.*(filter([0 -pcc(2:end)],1,signal(:,i)));%元データに対して後期残響フィルターかけを行い，後期残響を導出
%% DGTの準備
[win,diffWin] = generalizedCosWin(windowLen,'nuttall4termC1'); %窓関数の生成
signal_pad = zeroPaddingForDGT(signal(:,i),shiftLen,fftLen); % DGTを実行させる前にゼロ埋め，ここらへんからモノラルで処理
LRP = zeroPaddingForDGT(LRP,shiftLen,fftLen);
%% DGT(スペクトログラムの計算)
spec = DGT(signal_pad,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);%元信号のスペクトログラムの計算
LRPspec = DGT(LRP,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);%後期残響のスペクトログラムの計算
diffSpec = DGT(signal_pad,diffWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);% 瞬時周波数の計算をする際に必要
%% 図の表示の仕方はちょっと考える
% figure, imagesc(20*log10(abs(spec))), colorbar, axis xy% figure,
% imagesc(20*log10(abs(LRPspec))), colorbar, axis xy%

%% 後期残響除去
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
%% 位相の計算
Fai = spec./S_hat_abso;%複素位相の算出
S_hat = S_hat_abso.*Fai; %S_hat = S_hat_abso * exp(i*Fai);

for m = 1:dif_line
    for n = 1:dif_row
        if isnan(S_hat(m,n))
            S_hat(m,n) = 0;
        end
    end
end
%% 位相修正
IF = calcInstFreq(S_hat,diffSpec,fftLen,windowLen,rotateFlag);% 瞬時周波数を計算する関数，源信号のDGT，微分窓をDGTしたものを使う
iPCspec = instPhaseCorrection(S_hat,IF,shiftLen,fftLen);% 瞬時位相修正済スペクトログラムを計算する関数

figure, cRange = max(abs(spec(:)))/5;
subplot(1,3,1), imagesc(abs(spec(1:40,fix(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Amplitude'
subplot(1,3,2), imagesc(real(spec(1:40,fix(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Real part'
subplot(1,3,3), imagesc(real(iPCspec(1:40,fix(end-end/128):end))), axis xy, caxis(cRange*[-1 1])
title 'Real part (iPC)'
recon_spec = invInstPhaseCorrection(iPCspec,IF,shiftLen,fftLen); % 位相逆回転
figure, imagesc(20*log10(abs(recon_spec))), colorbar, axis xy
%% 直接音再構成
dualWin = calcCanonicalDualWindow(win,shiftLen);% 
reconst = invDGT(recon_spec,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);% 位相修正逆変換したものを逆DGTし

%% 左右結合
if i == 1
    len_reconst = length(reconst);
    len_rev = length(LRP);
    direct = zeros(len_reconst,2);
    rev = zeros(len_rev,2);    
end
direct(:,i) = reconst;
rev(:,i) = LRP;
end
%% オーディオ書き出し
audiowrite('direct.wav',direct,Fs);
audiowrite('diffuse.wav',rev,Fs);
audiowrite('original.wav',signal,Fs);
