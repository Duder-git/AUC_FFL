% test
load B_FFL_wave.mat

a = B_FFL_wave;

Ts = 1/2e6;                 % 时间间隔
T = 1/2.5e3;                % 单激励周期点数
t = (0:Ts:T);               % 时间序列

Nt = 1:32;
Y = zeros(1,801);
figure
for i =1:32
plot3(t,Y+50/32*i,squeeze(a(1,i,:)))
hold on
end
xlabel("时间(s)")
ylabel("梯度编码")
zlabel("磁场强度(mT)")
%%
figure
plot3(t,Y+1,squeeze(a(1,1,:)))
hold on
plot3(t,Y+2,squeeze(a(1,2,:)))