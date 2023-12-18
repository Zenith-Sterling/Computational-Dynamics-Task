%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Compute the stiffness matrix of quad                        *
%*                                                                 *
%* - Call procedures:                                              *
%*     QuadStiff.m - InitQuad()                                    *
%*     ./ReadQuad.m - ReadQuad()                                   *
%*     SRC/Mechanics/Addres.m - Addres()                           *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Mechanics/GetStiff.m                                    *
%*                                                                 *
%* - Programmed by:                                                *
%*     Shuai Zhang                                                 *
%* *****************************************************************



function QuadStiff()

% Init variables of the element
InitQuad();

% Read Material and Elements
ReadQuad();

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

function InitQuad()
global sdata;
sdata.NNODE = 4;
sdata.NDOF = 2;

end

% Assemble structure stiffness matrix
function Assemble()
global sdata;
global cdata;
K = zeros(8, 8, 'double');
sdata.STIFF = zeros(sdata.NWK, 1, 'double');

NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
E = sdata.E; nu = sdata.nu; LM = sdata.LM;
for N = 1:NUME
    MTYPE = MATP(N);
    
%   compute the matrix D
    D = E(MTYPE)/(1-nu(MTYPE)^2)*[1 nu(MTYPE) 0;nu(MTYPE) 1 0;0 0 (1-nu(MTYPE))/2];

%   compute the matrix B & K 
    K = zeros(8,8);
    P = [XYZ(1, N) XYZ(2, N);XYZ(4, N) XYZ(5, N);
        XYZ(7, N) XYZ(8, N);XYZ(10, N) XYZ(11, N)];
    s = [-sqrt(3)/3 sqrt(3)/3];
    t = [-sqrt(3)/3 sqrt(3)/3];
    for i = 1:2
        for j = 1:2
            N_st = 0.25*[t(j)-1 1-t(j) 1+t(j) -1-t(j);
                s(i)-1 -1-s(i) 1+s(i) 1-s(i)];
            J = N_st*P;
            N_xy = J\N_st;
            B = [N_xy(1,1) 0 N_xy(1,2) 0 N_xy(1,3) 0 N_xy(1,4) 0;
                0 N_xy(2,1) 0 N_xy(2,2) 0 N_xy(2,3) 0 N_xy(2,4);
                N_xy(2,1) N_xy(1,1) N_xy(2,2) N_xy(1,2) N_xy(2,3) N_xy(1,3) N_xy(2,4) N_xy(1,4)];
            K = K + B'*D*B*det(J);
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
M1 = zeros(8, 8, 'double');
M2 = zeros(8, 8, 'double');
sdata.Mass1 = zeros(sdata.NWK, 1, 'double');
sdata.Mass2 = zeros(sdata.NWK, 1, 'double');

NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
E = sdata.E; rho=sdata.rho;LM = sdata.LM;
for N = 1:NUME
    MTYPE = MATP(N);
    
%   compute the length of truss element
    DX1 = XYZ(1, N) - XYZ(3, N);
    DY1 = XYZ(2, N) - XYZ(4, N);
    XL2 = DX1*DX1 + DY1*DY1;
    a = sqrt(XL2);
    DX1 = XYZ(3, N) - XYZ(5, N);
    DY1 = XYZ(4, N) - XYZ(6, N);
    XL2 = DX1*DX1 + DY1*DY1;
    b = sqrt(XL2);
    DX1 = XYZ(5, N) - XYZ(7, N);
    DY1 = XYZ(6, N) - XYZ(8, N);
    XL2 = DX1*DX1 + DY1*DY1;
    c = sqrt(XL2);
    DX1 = XYZ(7, N) - XYZ(1, N);
    DY1 = XYZ(8, N) - XYZ(2, N);
    XL2 = DX1*DX1 + DY1*DY1;
    d = sqrt(XL2);
    z=(a+b+c+d)/2;
    S=2*sqrt((z-a)*(z-b)*(z-c)*(z-d) );

    W=rho(MTYPE)*S;
    
    M1=W/4*[4/9 0 2/9 0 1/9 0 2/9 0;
       0 4/9 0 2*9 0 1/9 0 2/9;
       2/9 0 4/9 0 2/9 0 1/9 0;
       0 2/9 0 4/9 0 2/9 0 1/9;
       1/9 0 2/9 0 4/9 0 2/9 0;
       0 1/9 0 2/9 0 4/9 0 2/9;
       2/9 0 1/9 0 2/9 0 4/9 0;
       0 2/9 0 1/9 0 2/9 0 4/9];
%     M2=W/2*[1 0 0 0 0 0;
%        0 1 0 0 0 0;
%        0 0 1 0 0 0;
%        0 0 0 1 0 0;
%        0 0 0 0 1 0;
%        0 0 0 0 0 1];
    
%   SRC/Mechanics/ADDBAN.m
    AddM(M1,M2, LM(:, N));
    
end

% The third time stamp
cdata.TIM(3, :) = clock;

end



