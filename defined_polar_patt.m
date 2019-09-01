%% レーダーチャートを表示するプログラム
%インパルス応答を読み込んで指定の周波数に対して1周分のチャートを描画するプログラム
%インパルス応答のファイルは　[種類がわかるような文字]_[角度の値(半角)]_imp.wav　で統一　一箇所のフォルダに格納しておく
%実行はインパルス応答を保存してあるフォルダを選択
%20170824 Sakaguchi修正

%%
clear;
close all;
clc;
%% --音のファイルを選択する処理

m = 360;
num = m-10;
f = 100;%取り出す周波数を設定
%for h = from:1000:to
for v = 0:10:num%0から350まで回す(10°ごとのデータを取得)
    
    tmp = v/10;%10°づつデータを取ったのでデータ数は10で割った数になる
    s = num2str(v);%vを文字配列型に変換(fnameで使うため)
    fname =['100_cut4006_0_' s '_imp.wav'];
    PathName = './';
    T = strcat(PathName, fname);
    [data,Fs] = audioread(T);%データの読み込み
    disp(v)%今どの角度を実行しているのかをコマンドに表示
    data = data(:,1);%モノラルデータがほしいのでLだけデータを取り出す

    %%FFT

    len_y = length(data);
    %長さの調整
    if (len_y < 2*Fs)
        add1 = 2*Fs - len_y;
        z_add1 = zeros(add1,1);
        yn = vertcat(data,z_add1);
    else
        yn = data;
        
    end
    len_yn = length(yn);
    Y = fft(yn);
    %グラフ描画のための準備

    %周波数領域の準備
    len_Y = length(Y);%周波数分解能の計算
    len_Y = len_Y - rem(len_Y,2);
    time_Y = (1:len_Y)/Fs;
    Y = Y(1:len_Y);
    f_vin = Fs/len_Y;
    fre = 0:f_vin:Fs/2; %周波数分解能の計算ここまで

    %指定の周波数を取り出す
    deltaF = ( 0:len_Y-1).*Fs./ len_Y; %横軸
    d = deltaF(3)-deltaF(2);%deltaFの差を求める(deltaFは等差数列)
    n = fix(( f + d ) / d);%n番目に取り出す周波数があるのかを求めている,nは整数なのでfixで丸め込み
    
    
    %deltaF = deltaF(n);
    avr = abs(Y);%複素数で値が出るのでその絶対値を求める
    
    %同じ周波数のデータを結びつける
        if(v == 1)
            avr_h = zeros(1,tmp);%0配列を準備して初期化
            avr_h(1,tmp+1) = avr(n);%tmpは0から始まるが、配列は1から始まるのでtmp+1からスタート
        else
            avr_h(1,tmp+1) = avr(n);
        end
        
    avr_av = 20.*log10(avr_h./(max(avr_h)));
    
end
disp(f)%どの周波数を計算したのか忘れないようにコマンドに表示
%dffr = dffr + 1;
%end
%% 初めてチャート表示するときに実行
flg =1;%フラグ指定
%% チャート表示
%図に表示する
if(flg == 1)
    figure(1)
    polarplot(avr_av,'LineWidth',1.0);
    %f_st =num2str(f);
    %fre =[f_st 'Hz'];
    %legend({fre})
    title('Polar Pattern Chart')
    hold on
    ax = gca;
    pax = gca;
    rticks(pax,[-20 -15 -10 -5 0])
    ax.ThetaTickLabel = {'0°','30°','60°','90°','120°','150°','180°','210°','240°','270°','300°','330°'};
    ax.RTickLabel = {'-20dB','','-10dB','','0dB'};
    ax.ThetaZeroLocation = 'top';
    ax.Color = 'w';
    ax.ThetaDir = 'clockwise';
    %ax.RDir = 'reverse';
    
    %印刷時に見えやすいようにする設定
    ax.GridAlpha = 0.7;
    ax.FontSize = 20;
    
    %表示範囲の設定
    ax.RLim = [-25 0];
    ax.ALim = [0 3];
    flg = 2;%フラグを2に変更
else
    avr_av1 = avr_av;
    polarplot(avr_av1,'-.','LineWidth',1.0);
    %f_st =num2str(f);
    %fre =[f_st 'Hz'];
    %legend({fre})
    hold on
end

%%
lgd = legend({'5kHz','6kHz','7kHz','8kHz','9kHz','10kHz'},'FontSize',12);
%title(lgd);