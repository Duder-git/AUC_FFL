%% 计算AUC_Flat
% 2022/4/27 Duder
% 2023/05/24 重新修改 创建主程序调用

% 说明
%   本程序用于仿真二维图像信号
%   使用贾广老师AUC-FFL方法

% 说明
%   pulsed Gradient Coils   脉冲梯度线圈
%       产生脉冲梯度信号，对FFL进行弛豫编码
%   receivc Coils           接收线圈
%       接收粒子信号，并根据粒子信号计算方波下面积(AUC）
%   AUCflat                 接收信号

% 使用
% step1 计算系统矩阵
% step2 计算接收信号
% step3 解系统矩阵求粒子分布

clc,clear
close all
% 导入函数
AUC_FFL = AUC_FFL_Func;
%% 参数设定
fprintf('参数设定\n');
fs = 2*10^6;            % 采样率  1.25MHz
fx = 2.5e3;             % 激励频率 2.5kHz
Ts = 1/fs;              % 时间间隔
Ns = round(fs/fx);      % 单位周期采样点数
T = 1/fx;               % 单激励周期点数
t = (0:Ts:T);               % 时间序列

EncodStep = 32;             % 离散点数
FFLnum = EncodStep;          
FOVSize = [0.02,0.02];      % 成像范围
x = (-1:2/(EncodStep-1):1)*1/2*FOVSize(2);
y = (-1:2/(FFLnum-1):1)*1/2*FOVSize(1);
x_move = (-1:2/(EncodStep-1):1)*1/2*FOVSize(2);
[X,Y] = meshgrid(x,y);      % 设定x和y方向位置信息



%% 脉冲梯度线圈 pulsed Gradient Coils
Bxmax = 25;                     % 弱梯度场(mT)
GridX = 2*Bxmax/FOVSize(2);     % 脉冲梯度(mT/m)
B_FFL_X = X*GridX;              % 脉冲梯度场(mT)

GridY = 2*Bxmax/FOVSize(1);     % 弛豫编码梯度(mT/m)
B_FFL_Y = Y*GridY;              % 移动梯度场(mT)

B_FFL = B_FFL_X + B_FFL_Y;

% 方波脉冲信号
B_FFL_wave = AUC_FFL.Ladder1(t,B_FFL);
save B_FFL_wave.mat B_FFL_wave;
% figure
% plot(squeeze(B_FFL_wave(EncodStep/2,1,:)))
% hold on
% plot(squeeze(B_FFL_wave(EncodStep/2,8,:)))
% hold on
% plot(squeeze(B_FFL_wave(EncodStep/2,25,:)))
% legend('-25mT','-12.5mT','12.5mT')
% xlabel('time(s)')
% ylabel('弛豫编码场强(mT)')

%% 系统矩阵
fprintf('计算系统矩阵\n');
AUC_rx = AUC_FFL.CalCulate_AUC_RX(t,B_FFL_wave,fs);
save('AUC_rx.mat',"AUC_rx" );



%% 图像设置
Map = imread("5points.png");
Map = im2gray(Map);
Map = imresize(Map,[FFLnum,EncodStep]);
Map = double(Map);

figure
imagesc(Map)
colorbar
title('粒子分布')
axis equal
axis off

%% 计算接收信号
fprintf('计算二维接收信号\n');

% 计算信号
AUC_flat = AUC_FFL.CalCulate_AUC_Flat_2D(t,B_FFL_wave,fs,Map);
save('AUC_flat.mat', "AUC_flat")


%% 解系统矩阵
clc,clear

fprintf('开始求解系统矩阵\n')
load('AUC_rx.mat');
load('AUC_flat.mat');

figure
mesh(AUC_rx);
axis equal

FFLnum = 32;%  FFL数量
[m,n] = size(AUC_rx);
C_img = zeros(FFLnum,n);
for i = 1:FFLnum
    %[C,iter] = ART(AUC_rx,AUC_flat(:,i),zeros(n,1),1e-3);
    C = kaczmarzReg( AUC_rx ,AUC_flat(:,i),100 ,1*10^-6 , 0,1,1);
    C_img(i,:) = C';
    fprintf('梯度编码FFL 第%d条\n',i);
end
figure
imagesc(C_img)
colorbar
title('重建图像')
axis equal
axis off