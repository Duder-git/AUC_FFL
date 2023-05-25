function [ x ] = kaczmarzReg( A ,b ,iterations ,lambd ,shuff ,enforceReal ,enforcePositive)
%c_normReguArt = kaczmarzReg(S_truncated(:,:),u_mean_truncated(:),100       , 1*10^-6 , 0    ,1          ,1               )
% regularized Kaczmarz
% As published here: http://iopscience.iop.org/article/10.1088/0031-9155/55/6/003/meta on page 1582.
% Other references : Saunders M 1995: Solution of sparse rectangular
% systems using LSQR and CRAIG  使用LSOQR和CRAIG系统百分比
% or
% From Dax A 1993: On row relaxation methods for large constrained 关于大型约束的行松弛方法
% least squares problems 最小二乘问题

% initialization of the variable
% 变量初始化
[N, M] = size(A);   %[1936,1528]
x = complex(zeros(N,1)); 
residual = complex(zeros(M,1));
rowIndexCycle = 1:M;

% calculate the energy of each frequency component
% 计算每个频率分量的能量
energy = sqrt(sum(abs(A.*A),1));	


% may use a randomized Kaczmarz
% 可能使用随机化的Kaczmarz
if shuff
    rowIndexCycle = randperm(M);
end

% estimate regularization parameter
% 估计正则化参数
lambdZero = sum(energy.^2)/N;
lambdIter = lambd*lambdZero;

for l = 1:iterations        % 迭代次数
    for m = 1:M             % 频率分量
        k = rowIndexCycle(m);   % 随机索引
        
        if energy(k) > 0        % 对应频率分量能量大于0
            tmp = A(:,k).'*x;
            beta = (b(k) - tmp - sqrt(lambdIter)*residual(k)) / (energy(k)^2 + lambdIter);
            x = x + beta*conj(A(:,k)); %共轭复数
            residual(k) = residual(k) + beta*sqrt(lambdIter);
        end

    end
    
    if enforceReal && ~isreal(x)
        x = complex(real(x),0);
    end
    
    if enforcePositive
        x(real(x) < 0) = 0;
    end
end

end