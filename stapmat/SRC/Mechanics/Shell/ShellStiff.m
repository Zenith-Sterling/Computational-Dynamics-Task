%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Compute the stiffness matrix of quad                        *
%*                                                                 *
%* - Call procedures:                                              *
%*     ShellStiff.m - InitQuad()                                   *
%*     ./ReadShell.m - ReadQuad()                                  *
%*     SRC/Mechanics/Addres.m - Addres()                           *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Mechanics/GetStiff.m                                    *
%*                                                                 *
%* - Programmed by:                                                *
%*     Zhang  Chen                                                 *
%* *****************************************************************



function ShellStiff()

% Init variables of the element
InitShell();

% Read Material and Elements
ReadShell();

fprintf('Solution phase ...\n\n');

% calculate addresses of diagonal elements
Addres();

% Data check Or Solve
global cdata;
if (cdata.MODEX == 0) 
    cdata.TIM(3,:) = clock;
    cdata.TIM(4,:) = clock;
    cdata.TIM(5,:) = clock;
    return; 
end

% Assemble structure stiffness matrix
Assemble();

AssembleMass();

end


% ----------------------- Functions -----------------------------------

% Init parameters of quad element

function InitShell()
global sdata;
sdata.NNODE = 4;
sdata.NDOF = 3;   %自由度：w，θx，θy

end

% Assemble structure stiffness matrix
function Assemble()
global sdata;
global cdata;
K = zeros(12, 12, 'double');
sdata.STIFF = zeros(sdata.NWK, 1, 'double');

NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
E = sdata.E; nu = sdata.nu; LM = sdata.LM; thick=sdata.thick;%厚度
for N = 1:NUME
    MTYPE = MATP(N);
    
%   compute the matrix D
    D_plane = E(MTYPE)/(1-nu(MTYPE)^2)*[1 nu(MTYPE) 0;nu(MTYPE) 1 0;0 0 (1-nu(MTYPE))/2]; %平面应力矩阵Dplane
    k=5/6; %剪切校正因子，通常为5/6
    D_shear=k*E(MTYPE)/(2*(1+nu(MTYPE)))*[1 0;0 1]; %横向剪切矩阵Dshear


%   compute the matrix B & K 
    K = zeros(12,12);
    P = [XYZ(1, N) XYZ(2, N);XYZ(4, N) XYZ(5, N);
        XYZ(7, N) XYZ(8, N);XYZ(10, N) XYZ(11, N)];   %四节点坐标
    s = [-sqrt(3)/3 sqrt(3)/3];
    t = [-sqrt(3)/3 sqrt(3)/3];   %高斯积分点
    for i = 1:2
        for j = 1:2        %循环遍历所有高斯积分点
                N_st = 0.25*[1+t(j) -1-t(j) -1+t(j) 1-t(j);
                        s(i)+1 1-s(i) -1+s(i) -1-s(i)];         %形函数的导数
            J = N_st*P;              %计算雅可比矩阵J
            detJ = det(J);           % 雅可比行列式
            N_xy = J\N_st;           %计算局部到全局坐标的映射
            B_plane = [0  0         N_xy(1,1) 0 0          N_xy(1,2) 0 0          N_xy(1,3) 0 0          N_xy(1,4);
                       0 -N_xy(2,1) 0         0 -N_xy(2,2) 0         0 -N_xy(2,3) 0         0 -N_xy(2,4) 0;
                       0 -N_xy(1,1) N_xy(2,1) 0 -N_xy(1,2) N_xy(2,2) 0 -N_xy(1,3) N_xy(2,3) 0 -N_xy(1,4) N_xy(2,4)];
            N1 = 0.25*(1+s(i))*(1+t(j));
            N2 = 0.25*(1-s(i))*(1+t(j));
            N3 = 0.25*(1-s(i))*(1-t(j));
            N4 = 0.25*(1+s(i))*(1-t(j));
            B_shear = [N_xy(2,1) -N1 0  N_xy(2,2) -N2 0  N_xy(2,3) -N3 0  N_xy(2,4) -N4 0;
                       N_xy(1,1) 0   N1 N_xy(1,2) 0   N2 N_xy(1,3) 0   N3 N_xy(1,4) 0   N4];
            % 计算 B_plane 的刚度贡献
            K_plane = B_plane' * D_plane * B_plane*2*1/3*(thick^3/8);
            % 计算 B_shear 的刚度贡献
            K_shear = B_shear' * D_shear * B_shear*thick;  % 假设这是一个 12x12 矩阵
            % 组合两个刚度矩阵
            K_combined = K_plane + K_shear;
            K = K + K_combined*detJ;
        end
    end
    
%   SRC/Mechanics/ADDBAN.m
    ADDBAN(K, LM(:, N));
    
end

% The third time stamp
cdata.TIM(3, :) = clock;

end


function AssembleMass()
global sdata;
global cdata;
M1 = zeros(12, 12, 'double');
M2 = zeros(12, 12, 'double');
sdata.Mass1 = zeros(sdata.NWK, 1, 'double');
sdata.Mass2 = zeros(sdata.NWK, 1, 'double');

NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
thick = sdata.thick; rho=sdata.rho;LM = sdata.LM;
    
%   SRC/Mechanics/ADDBAN.m
    %AddM(M1,M2, LM(:, N));
    
end

% % The third time stamp
% cdata.TIM(3, :) = clock;

%end



