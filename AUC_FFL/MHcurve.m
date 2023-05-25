function M = MHcurve(H,varargin)
% 用于计算特定磁场强度，粒子的磁化强度M
% H 磁场强度(T)
% D 粒子直径(m)     20nm
% c 溶液浓度(kg/L)  1mg/ml
% M 磁化强度(KA/m)

    ip = inputParser;
    addParameter(ip,'D',20e-9);     % 20 nm
    addParameter(ip,'c',1e-3);      % 1 mg/ml
    parse(ip,varargin{:});
    D = ip.Results.D;
    c_Fe = ip.Results.c;

    u0 = 4*pi*1e-7; %真空磁导率
    MS = 446000;    %磁粒子饱和磁化强度(A/m)
    kB = 1.38e-23;  %玻尔兹曼常数
    Tp = 300;       %温度

    V_core = 1/6*pi*D^3;        %磁粒子体积(m^3)
    c = Particles(D,c_Fe);  

    m = V_core*MS;          %粒子磁矩
    k = m*u0/(kB*Tp);       %计算出β作为系数

    M = m*c*langevin(k*H/u0);   %A/m磁化强度M(A/m * m3/L = KA/m)

end

function out = langevin(in)   
    out = zeros(size(in));
    ind = find(in ~= 0);
    out(ind) = coth(in(ind))-1./in(ind);
    out(in==0) = 0;
    
end

function c = Particles(D,cFe)
    V=1/6*pi*D^3;           %计算磁粒子体积(m^3)
    p_Fe3O4 = 5200;         %四氧化三铁密度(Kg/m^3)
    m_core = p_Fe3O4*V;     %每个粒子核质量(Kg)
    mol_Fe3O4 = 0.232;      %四氧化三铁摩尔质量(Kg/mol)
    mol_Fe = 0.056;         %铁摩尔质量(Kg/mol)
    m_Fe = m_core*3*mol_Fe/mol_Fe3O4;%每个核中铁离子质量(Kg)
    c = cFe/m_Fe;           %溶液中四氧化三铁浓度(/L)
end
