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
%*     Shuai Zhang                                                 *
%* *****************************************************************

function QuadStress(NUM, NG)

% Get global data
global cdata;
global sdata;

IOUT = cdata.IOUT;
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ;
E = sdata.E; nu = sdata.nu; LM = sdata.LM;
U = sdata.DIS(:, NUM);

fprintf(IOUT, ['\n\n  S T R E S S  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '       ELEMENT             SIGMAxx         SIGMAyy         SIGMAxy\n' ...
    '       NUMBER\n'], NG);

for N = 1:NUME
    MTYPE = MATP(N);
   
%   displacement
    for i = 1:8
        if LM(i, N) ~= 0
            u(i) = U(LM(i, N));
        else
            u(i) = 0;
        end
    end
    d = [u(1);u(2);u(3);u(4);u(5);u(6);u(7);u(8)];
    B1 = 0;

%   compute the matrix D
    D = E(MTYPE)/(1-nu(MTYPE)^2)*[1 nu(MTYPE) 0;nu(MTYPE) 1 0;0 0 (1-nu(MTYPE))/2];

%   compute the matrix B 
    P = [XYZ(1, N) XYZ(2, N);XYZ(4, N) XYZ(5, N);XYZ(7, N) XYZ(8, N);XYZ(10, N) XYZ(11, N)];
    s = [-sqrt(3)/3 sqrt(3)/3];
    t = [-sqrt(3)/3 sqrt(3)/3];
    for i = 1:2
        for j = 1:2
            N_st = 0.25*[t(j)-1 1-t(j) 1+t(j) -1-t(j);s(i)-1 -1-s(i) 1+s(i) 1-s(i)];
            J = N_st*P;
            N_xy = J\N_st;
            B = [N_xy(1,1) 0 N_xy(1,2) 0 N_xy(1,3) 0 N_xy(1,4) 0;0 N_xy(2,1) 0 N_xy(2,2) 0 N_xy(2,3) 0 N_xy(2,4);N_xy(2,1) N_xy(1,1) N_xy(2,2) N_xy(1,2) N_xy(2,3) N_xy(1,3) N_xy(2,4) N_xy(1,4)];    
            B1 = 0.25*B + B1;
        end
    end
    eps = B1*d;
    sigma = D*eps;
    
    fprintf(IOUT, ' %10d     %13.6e     %13.6e    %13.6e\n', N, sigma(1), sigma(2) ,sigma(3));
end

end