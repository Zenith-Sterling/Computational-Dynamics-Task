%* *****************************************************************
%* - Function of STAPMAT in solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate stresses                                       *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Solver/GetStress.m                                      *
%*                                                                 *
%* - Programmed by:                                                *
%*      Zhang Chen                                                 *
%* *****************************************************************

function QuadStress(NUM, NG)

% Get global data
global cdata;
global sdata;

IOUT = cdata.IOUT;
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ;
E = sdata.E; nu = sdata.nu; LM = sdata.LM;thick=sdata.thick;
U = sdata.DIS(:, NUM);

fprintf(IOUT, ['\n\n  S T R E S S  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '       ELEMENT             SIGMAxx         SIGMAyy         SIGMAxy         TAUxz          TAUyz\n' ...
    '       NUMBER\n'], NG);

for N = 1:NUME
    MTYPE = MATP(N);
   
%   displacement
    for i = 1:12
        if LM(i, N) ~= 0
            u(i) = U(LM(i, N));
        else
            u(i) = 0;
        end
    end
   
    d = [u(1);u(2);u(3);u(4);u(5);u(6);u(7);u(8);u(9);u(10);u(11);u(12)];
    
    B1_plane=0;
    B1_shear=0;
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
            B1_plane=B_plane+B1_plane;
            B1_shear=B_shear+B1_shear;
        end
    end
    z=thick/2; %上表面
    eps_plane = z*B1_plane*d; 
    eps_shear = B1_shear*d;
    sigma_plane = D_plane*eps_plane; %σx、σy、σxy
    sigma_shear = D_shear*eps_shear; %τxz、τyz
    
    fprintf(IOUT, ' %10d     %13.6e     %13.6e    %13.6e     %13.6e     %13.6e\n', N, sigma_plane(1), sigma_plane(2) ,sigma_plane(3), sigma_shear(1), sigma_shear(2));
end

end