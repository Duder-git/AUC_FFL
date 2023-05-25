function AUC_FFL = AUC_FFL_Func()
% AUC_FFL 复现所需函数
    AUC_FFL.CalCulate_AUC_RX = @CalCulate_AUC_RX;
    AUC_FFL.CalCulate_AUC_Flat = @CalCulate_AUC_Flat;
    AUC_FFL.CalCulate_AUC_Flat_2D = @CalCulate_AUC_Flat_2D;
    AUC_FFL.Ladder1 = @Ladder1;
    AUC_FFL.receiveRX = @receiveRX;
    AUC_FFL.receiveFFL = @receiveFFL;
%     AUC_FFL.DebyeRelaxation = @DebyeRelaxation;
%     AUC_FFL.DebyeRelaxation2D = @DebyeRelaxation2D;
%     AUC_FFL.CalAUC = @CalAUC;
%     AUC_FFL.CalAUC2D = @CalAUC2D;

end


%% 函数
function AUC_rx = CalCulate_AUC_RX(t,B_FFL_wave,fs)
%% 计算系统矩阵
U = receiveRX(B_FFL_wave,fs);
St = DebyeRelaxation2D(U,t);
AUC_rx = CalAUC2D(St,U);

end


function AUC_flat = CalCulate_AUC_Flat_2D(t,B_FFL_wave,fs,Map)
%% 计算二维接收信号
[m,n,p] = size(B_FFL_wave);
U = zeros(n,p,m);
for i = 1:m
    Mapi = repmat(Map(i,:),m,1);
    U(:,:,i) = receiveFFL(B_FFL_wave,fs,Mapi);
end
U = permute(U,[3,1,2]);
St = DebyeRelaxation2D(U,t);
AUC_flat = CalAUC2D(St,U);
AUC_flat = AUC_flat';

end


function wave1 = Ladder1(t,BY)
%% 用于计算相同压摆率下不同选频强度信号
fs = 1;
slewRate = 25*fs/50;       % 压摆率(mT/s)
[m,n] = size(BY);
p = length(t);      % 末尾
wave1 = zeros(m,n,p);
beginPoint = 200;
k = slewRate*1/fs;
%y = k*(x-126);
for i = 1:n
    for j = 1:m
        wave = zeros(1,p);
        A = BY(i,j);
        a = floor(-A/k+beginPoint); % 第一个点
        b = floor(A/k+beginPoint);  %
        if a<b
            flag = 1;
        else %a>b
            flag = -1;
            tmp = a;
            a = b;
            b = tmp;
        end
        wave(1:a) = -A;
        wave((a+1):b) = flag*k*(((a+1):b)-beginPoint);
        wave((b+1):((p-1)/2+1)) = A;
        tmp = fliplr(wave);
        wave(((p-1)/2+1):p) = tmp(((p-1)/2+1):p);
        wave1(i,j,:) = wave;
    end
end
end

function U = receiveRX(B,fs)
%% 计算系统矩阵
u0 = 4*pi*1e-7;     % 真空磁导率(N/A^2) (T*m/A)
S0 = 9*1e-3;        % 线圈灵敏度(T/A)
U = zeros(size(B));
Mx =MHcurve(B*1e-3);
U(:,:,2:end) = u0*S0*diff(Mx*1e3,1,3)*fs;
end


function St = DebyeRelaxation2D(S,t)
%% 德拜弛豫
tao = 10e-6;   %s          
r = 1e-6/tao*exp(-t/tao);
r = r/sum(r);
r = [zeros(1,length(t)),r];

[m,n,p] = size(S);
for i = 1:m
    for j = 1:n
        St(i,j,:) = conv(squeeze(S(i,j,:)),r,'same');
    end
end
end


function AUC = CalAUC2D(Sigt,Sig)
%% 用于计算AUC
[m,n,p] = size(Sigt);
index = 200:580;
AUC = zeros(m,n);
SigtR = Sigt(:,:,index);
SigR = Sig(:,:,index);
for i = 1:n %FFL上点个数
    for j = 1:m
        SigR(i,j,:) = SigR(i,j,:)./max(SigR(i,j,:),2);
        index = abs(SigR(i,j,:))<=0.1;
        AUC(j,i) = sum(SigtR(j,i,index));
    end
end
end


function U = receiveFFL(B_FFL_wave,fs,Map)
%% 计算单行FFL信号
u0 = 4*pi*1e-7;     % 真空磁导率(N/A^2) (T*m/A)
S0 = 9*1e-3;        % 线圈灵敏度(T/A)

% 循环计算每次梯度编码
[m,n,k] = size(B_FFL_wave);

U = zeros(n,k);
Mx = zeros(m,k);

for i = 1:k
    Bt = B_FFL_wave(:,:,i);
    M = MHcurve(Map.*Bt*1e-3);
    Mx(:,i) = sum(M,2);
end
U(:,2:end) = u0*S0*diff(Mx*1e3,1,2)*fs;
end


function St = DebyeRelaxation(S,t)
%% 德拜弛豫
tao = 20e-6;   %s
r = 1e-6/tao*exp(-t/tao);
r = r/sum(r);
r = [zeros(1,length(t)),r];

[m,p] = size(S);
for i = 1:m
    St(i,:) = conv(squeeze(S(i,:)),r,'same');
end
end


function AUC = CalAUC(Sigt,Sig)
%% 用于计算AUC
index = 200:580;
[m,p] = size(Sigt);
AUC = zeros(m,1);
% 获取信号
SigtR = Sigt(:,index);
SigR = Sig(:,index);


for i = 1:m %FFL行数
    SigR(i,:) = SigR(i,:)./max(SigR(i,:),2);
    index = abs(SigR(i,:))<=0.1;
    AUC(i) = sum(SigtR(i,index));
end
end