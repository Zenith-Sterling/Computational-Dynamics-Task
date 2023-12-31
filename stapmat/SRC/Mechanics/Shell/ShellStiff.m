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
global sdata;
if (cdata.MODEX == 0) 
    cdata.TIM(3,:) = clock;
    cdata.TIM(4,:) = clock;
    cdata.TIM(5,:) = clock;
    return; 
end

% Assemble structure stiffness matrix
Assemble();

AssembleMass();

if sdata.f(1,1) ~= 0
    AssembleR();
end


end


% ----------------------- Functions -----------------------------------

% Init parameters of quad element

function InitShell()
global sdata;
sdata.NNODE = 4;
sdata.NDOF = 6;   %自由度：u,v,w,θx,θy,θz

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
    D = E(MTYPE)/(1-nu(MTYPE)^2)*[1 nu(MTYPE) 0;nu(MTYPE) 1 0;0 0 (1-nu(MTYPE))/2];


%   compute the matrix B & K 
    K = zeros(12,12);
    K_uv = zeros(8,8);

    % 局部坐标系建立
    r1 = [XYZ(4,N)-XYZ(1,N);XYZ(5,N)-XYZ(2,N);XYZ(6,N)-XYZ(3,N)]; %节点1到节点2
    r2 = [XYZ(10,N)-XYZ(1,N);XYZ(11,N)-XYZ(2,N);XYZ(12,N)-XYZ(3,N)]; %节点1到节点4
    r23 = [XYZ(7,N)-XYZ(4,N);XYZ(8,N)-XYZ(5,N);XYZ(9,N)-XYZ(6,N)]; %节点2到节点3
    r3 = cross(r1,r2);
    e1 = r1/norm(r1); %局部坐标系x轴
    e3 = r3/norm(r3); %局部坐标系z轴
    e2 = cross(e3,e1); %局部坐标系y轴
    x2 = norm(r1); %局部坐标系下2点的x坐标
    x3 = x2 + r23'*e1;
    y3 = r23'*e2;
    x4 = r2'*e1;
    y4 = r2'*e2;

    P = [0 0;x2 0;
        x3 y3;x4 y4];   %四节点坐标
    s = [-sqrt(3)/3 sqrt(3)/3];
    t = [-sqrt(3)/3 sqrt(3)/3];   %高斯积分点
    R=25;
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
            B_uv = [N_xy(1,1) 0 N_xy(1,2) 0 N_xy(1,3) 0 N_xy(1,4) 0;
                0 N_xy(2,1) 0 N_xy(2,2) 0 N_xy(2,3) 0 N_xy(2,4);
                N_xy(2,1) N_xy(1,1) N_xy(2,2) N_xy(1,2) N_xy(2,3) N_xy(1,3) N_xy(2,4) N_xy(1,4)];
            
%             K_uv = K_uv+B_uv'*D*B_uv*detJ*(thick+thick^3/(12*R^2));
            K_uv = K_uv+B_uv'*D*B_uv*detJ*(thick);
            % 计算 B_plane 的刚度贡献
%             K_plane = B_plane' * D_plane * B_plane*(1/12*thick^3+1/240*thick^5/(R^2));
            K_plane = B_plane' * D_plane * B_plane*(1/12*thick^3);
            % 组合两个刚度矩阵
            K = K + K_plane*detJ;
        end
    end

    % N_st = 0.25*[1 -1 -1 1;
    %     1 1 -1 -1];         %形函数的导数
    % J = N_st*P;              %计算雅可比矩阵J
    % detJ = det(J);           % 雅可比行列式
    % N_xy = J\N_st;           %计算局部到全局坐标的映射
    % B_plane = [0  0         N_xy(1,1) 0 0          N_xy(1,2) 0 0          N_xy(1,3) 0 0          N_xy(1,4);
    %     0 -N_xy(2,1) 0         0 -N_xy(2,2) 0         0 -N_xy(2,3) 0         0 -N_xy(2,4) 0;
    %     0 -N_xy(1,1) N_xy(2,1) 0 -N_xy(1,2) N_xy(2,2) 0 -N_xy(1,3) N_xy(2,3) 0 -N_xy(1,4) N_xy(2,4)];
    % B_uv = [N_xy(1,1) 0 N_xy(1,2) 0 N_xy(1,3) 0 N_xy(1,4) 0;
    %     0 N_xy(2,1) 0 N_xy(2,2) 0 N_xy(2,3) 0 N_xy(2,4);
    %     N_xy(2,1) N_xy(1,1) N_xy(2,2) N_xy(1,2) N_xy(2,3) N_xy(1,3) N_xy(2,4) N_xy(1,4)];
    % K_uv = K_uv+4*B_uv'*D*B_uv*detJ*thick;
    % % 计算 B_plane 的刚度贡献
    % K_plane = 4*B_plane' * D_plane * B_plane*2*1/3*(thick^3/8);
    % % 组合两个刚度矩阵
    % K = K + K_plane*detJ;
    % 
    % 单点高斯积分
    N_st = 0.25*[1 -1 -1 1;
        1 1 -1 -1];         %形函数的导数
    J = N_st*P;              %计算雅可比矩阵J
    detJ = det(J);           % 雅可比行列式
    N_xy = J\N_st;           %计算局部到全局坐标的映射
    B_shear = [N_xy(2,1) -0.25 0  N_xy(2,2) -0.25 0  N_xy(2,3) -0.25 0  N_xy(2,4) -0.25 0;
                       N_xy(1,1) 0 0.25 N_xy(1,2) 0 0.25 N_xy(1,3) 0 0.25 N_xy(1,4) 0 0.25];
    % 计算 B_shear 的刚度贡献
%     K_shear = 4 * B_shear' * D_shear * B_shear*(thick+thick^3/(12*R^2));
    K_shear = 4 * B_shear' * D_shear * B_shear*(thick);
    K = K + K_shear*detJ;
    

    % 扩展维数
    Q1 = [0 0 1 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 1 0];
    Q2 = [1 0 0 0 0 0;
        0 1 0 0 0 0];
    Q1_total = blkdiag(Q1,Q1,Q1,Q1);
    Q2_total = blkdiag(Q2,Q2,Q2,Q2);

    % 转换矩阵
    T = [e1'*[1;0;0] e1'*[0;1;0] e1'*[0;0;1];
        e2'*[1;0;0] e2'*[0;1;0] e2'*[0;0;1];
        e3'*[1;0;0] e3'*[0;1;0] e3'*[0;0;1]];
    T_total = blkdiag(T,T,T,T,T,T,T,T);


    %板+膜的单元刚度阵
    K = Q1_total'*K*Q1_total;
    K = Q2_total'*K_uv*Q2_total + K;

    % 罚函数法
    alpha = 1000000;
    B = [0 0 0 0 0 1];
    B = [B,B,B,B];
    K = K + alpha*(B'*B);

    %全局坐标系的单元刚度阵
    K = T_total'*K*T_total; % 转置应该是逆但不太确定
%     K = T_total\K*T_total;

%   SRC/Mechanics/ADDBAN.m
    ADDBAN(K, LM(:, N));
    
end

% The third time stamp
cdata.TIM(3, :) = clock;

end


function AssembleMass()
global sdata;
global cdata;
M1 = zeros(24, 24, 'double');
M2 = zeros(24, 24, 'double');
sdata.Mass1 = zeros(sdata.NWK, 1, 'double');
sdata.Mass2 = zeros(sdata.NWK, 1, 'double');

NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
E = sdata.E; nu = sdata.nu; rho=sdata.rho; LM = sdata.LM; thick=sdata.thick;%厚度
for N = 1:NUME
    MTYPE = MATP(N);
   
    % 局部坐标系建立
    r1 = [XYZ(4,N)-XYZ(1,N);XYZ(5,N)-XYZ(2,N);XYZ(6,N)-XYZ(3,N)]; %节点1到节点2
    r2 = [XYZ(10,N)-XYZ(1,N);XYZ(11,N)-XYZ(2,N);XYZ(12,N)-XYZ(3,N)]; %节点1到节点4
    r23 = [XYZ(7,N)-XYZ(4,N);XYZ(8,N)-XYZ(5,N);XYZ(9,N)-XYZ(6,N)]; %节点2到节点3
    r3 = cross(r1,r2);
    e1 = r1/norm(r1); %局部坐标系x轴
    e3 = r3/norm(r3); %局部坐标系z轴
    e2 = cross(e3,e1); %局部坐标系y轴
    x2 = norm(r1); %局部坐标系下2点的x坐标
    x3 = x2 + r23'*e1;
    y3 = r23'*e2;
    x4 = r2'*e1;
    y4 = r2'*e2;

    P = [0 0;x2 0;
        x3 y3;x4 y4];   %四节点坐标


% 质量阵
    s = [-sqrt(3)/3 sqrt(3)/3];
    t = [-sqrt(3)/3 sqrt(3)/3];   %高斯积分点
    M_wtheta = zeros(12,12);
    M_uv = zeros(8,8);
    for i = 1:2
        for j = 1:2        %循环遍历所有高斯积分点
            N_st = 0.25*[1+t(j) -1-t(j) -1+t(j) 1-t(j);
                        s(i)+1 1-s(i) -1+s(i) -1-s(i)];         %形函数的导数
            J = N_st*P;              %计算雅可比矩阵J
            detJ = det(J);           % 雅可比行列式
            N1 = 0.25*(1+s(i))*(1+t(i));
            N2 = 0.25*(1-s(i))*(1+t(i));
            N3 = 0.25*(1-s(i))*(1-t(i));
            N4 = 0.25*(1+s(i))*(1-t(i));
            N_uv = [N1 0 N2 0 N3 0 N4 0;
                0 N1 0 N2 0 N3 0 N4];
            N_wtheta = [N1 0 0 N2 0 0 N3 0 0 N4 0 0;
                0 N1 0 0 N2 0 0 N3 0 0 N4 0;
                0 0 N1 0 0 N2 0 0 N3 0 0 N4];
            I_uv = [rho(MTYPE)*thick 0;
                0 rho(MTYPE)*thick];
            I_wtheta = [rho(MTYPE)*thick 0 0;
                0 rho(MTYPE)*thick^3/12 0;
                0 0 rho(MTYPE)*thick^3/12];
            M_uv = M_uv + N_uv'*I_uv*N_uv*detJ;
            M_wtheta = M_wtheta + N_wtheta'*I_wtheta*N_wtheta*detJ;
        end
    end


    % 扩展维数
    Q1 = [0 0 1 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 1 0];
    Q2 = [1 0 0 0 0 0;
        0 1 0 0 0 0];
    Q1_total = blkdiag(Q1,Q1,Q1,Q1);
    Q2_total = blkdiag(Q2,Q2,Q2,Q2);


    %板+膜的单元质量阵
    M1 = Q1_total'*M_wtheta*Q1_total;
    M1 = Q2_total'*M_uv*Q2_total + M1;

    % 罚函数法
    alpha = 100000;
    B = [0 0 0 0 0 1];
    B1 = [B,B,B,B];
    M1 = M1 + alpha*(B1'*B1);
    S = sum(M1,2);
    M2 = diag(S);


    
%   SRC/Mechanics/ADDBAN.m
    AddM(M1,M2, LM(:, N));
    
end

% The third time stamp
cdata.TIM(3, :) = clock;

end


function AssembleR()
global sdata;
global cdata;

f = sdata.f;
[n,~] = size(f);
XYZ = sdata.XYZ; 
LM = sdata.LM; thick=sdata.thick;%厚度

for k = 1:n
    for N = f(1):f(3):f(2)
        % 局部坐标系建立
        r1 = [XYZ(4,N)-XYZ(1,N);XYZ(5,N)-XYZ(2,N);XYZ(6,N)-XYZ(3,N)]; %节点1到节点2
        r2 = [XYZ(10,N)-XYZ(1,N);XYZ(11,N)-XYZ(2,N);XYZ(12,N)-XYZ(3,N)]; %节点1到节点4
        r23 = [XYZ(7,N)-XYZ(4,N);XYZ(8,N)-XYZ(5,N);XYZ(9,N)-XYZ(6,N)]; %节点2到节点3
        r3 = cross(r1,r2);
        e1 = r1/norm(r1); %局部坐标系x轴
        e3 = r3/norm(r3); %局部坐标系z轴
        e2 = cross(e3,e1); %局部坐标系y轴
        x2 = norm(r1); %局部坐标系下2点的x坐标
        x3 = x2 + r23'*e1;
        y3 = r23'*e2;
        x4 = r2'*e1;
        y4 = r2'*e2;
        %四节点坐标
        P = [0 0;x2 0;
            x3 y3;x4 y4];
        % 转换矩阵
        T = [e1'*[1;0;0] e1'*[0;1;0] e1'*[0;0;1];
            e2'*[1;0;0] e2'*[0;1;0] e2'*[0;0;1];
            e3'*[1;0;0] e3'*[0;1;0] e3'*[0;0;1]];
        g = blkdiag(T,T)*[f(k,4);f(k,5);f(k,6);f(k,7);f(k,8);f(k,9)]; %单元局部坐标系内体积力
        
        % 单元内载荷力r
        r = zeros(24,1);
        s = [-sqrt(3)/3 sqrt(3)/3];
        t = [-sqrt(3)/3 sqrt(3)/3];   %高斯积分点
        for i = 1:2
            for j = 1:2        %循环遍历所有高斯积分点
                N_st = 0.25*[1+t(j) -1-t(j) -1+t(j) 1-t(j);
                        s(i)+1 1-s(i) -1+s(i) -1-s(i)];         %形函数的导数
                J = N_st*P;              %计算雅可比矩阵J
                detJ = det(J);           % 雅可比行列式
                N1 = 0.25*(1+s(i))*(1+t(j));
                N2 = 0.25*(1-s(i))*(1+t(j));
                N3 = 0.25*(1-s(i))*(1-t(j));
                N4 = 0.25*(1+s(i))*(1-t(j));        
                NT = [diag([N1,N1,N1,N1,N1,N1]),diag([N2,N2,N2,N2,N2,N2]),diag([N3,N3,N3,N3,N3,N3]),diag([N4,N4,N4,N4,N4,N4])];
                r = r + NT'*g*thick*detJ;
            end
        end

        R = blkdiag(T',T',T',T',T',T',T',T')*r;

%   SRC/Mechanics/ADDBAN.m
    AddR(R, LM(:, N));

    end


end
  
    


% The third time stamp
    cdata.TIM(3, :) = clock;

end


