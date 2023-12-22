%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Read the element information of quad                        *
%*                                                                 *
%* - Call procedures:                                              *
%*     ReadQuad.m - ReadMaterial()                                 *
%*     ReadQuad.m - ReadElements()                                 *
%*                                                                 *
%* - Called by :                                                   *
%*     ./QuadStiff.m                                               *
%*                                                                 *
%* - Programmed by:                                                *
%*     Shuai Zhang                                                 *
%* *****************************************************************


function ReadShell()

% Read Material information
ReadMaterial()

% Read Element information
ReadElements()

% the second time stamp
global cdata;
cdata.TIM(2,:) = clock;

end

% ----------------------- Functions -----------------------------------
% Read Material information
function ReadMaterial()

global cdata;
global sdata;
% Get file pointers
IIN = cdata.IIN;
IOUT = cdata.IOUT;

if (cdata.NPAR(3) == 0) cdata.NPAR(3) = 1; end
fprintf(IOUT, '\n M A T E R I A L   D E F I N I T I O N\n');
fprintf(IOUT, '\n NUMBER OF DIFFERENT SETS OF MATERIAL\n');
fprintf(IOUT, ' AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . = %10d\n', ...
    cdata.NPAR(3));
fprintf(IOUT, '  SET       YOUNG''S     CROSS-SECTIONAL      THICKNESS\n');
fprintf(IOUT, ' NUMBER     MODULUS          Poisson                    \n');
fprintf(IOUT, '               E             nu                  \n');


% Read material datas
sdata.NUME = cdata.NPAR(2);
sdata.NUMMAT = cdata.NPAR(3);
NUMMAT = cdata.NPAR(3);
sdata.E = zeros(NUMMAT, 1, 'double');
sdata.nu = zeros(NUMMAT, 1, 'double');
sdata.thick = zeros(NUMMAT, 1, 'double');
sdata.rho=zeros(NUMMAT,1,'double');
for I = 1:cdata.NPAR(3)
    tmp = str2num(fgetl(IIN));
    N = round(tmp(1));
    sdata.E(N) = tmp(2);
    sdata.nu(N) = tmp(3);
    sdata.thick(N) = tmp(5);
    sdata.rho(N)=tmp(4);
    fprintf(IOUT, '%5d    %12.5e  %14.6e  %14.6e\n', N, tmp(2), tmp(3), tmp(4));
end

end

% Read elements information
function ReadElements()

global cdata;
global sdata;

% Get file pointer
IIN = cdata.IIN;
IOUT = cdata.IOUT;

fprintf(IOUT, '\n\n E L E M E N T   I N F O R M A T I O N\n');
fprintf(IOUT, '\n      ELEMENT          NODE          NODE          NODE          NODE       MATERIAL\n');
fprintf(IOUT, '      NUMBER-N           I1             I2             I3             I4       SET NUMBER\n');

% Get Position data
NUME = cdata.NPAR(2);
sdata.XYZ = zeros(12, NUME, 'double');
sdata.MATP = zeros(NUME, 1, 'int64');                 % the type of material
sdata.LM = zeros(12, NUME, 'double');                  % connectivity matrix
sdata.MHT = zeros(sdata.NEQ, 1, 'int64');
X = sdata.X; Y = sdata.Y; Z = sdata.Z; ID = sdata.ID;
XYZ = sdata.XYZ; MATP = sdata.MATP; LM = sdata.LM;

for N = 1:NUME
    tmp = str2num(fgetl(IIN));
    I1 = round(tmp(2));
    I2 = round(tmp(3));
    I3 = round(tmp(4));
    I4 = round(tmp(5));
    MTYPE = round(tmp(6));
    
%   Save element information
    XYZ(1, N) = X(I1);
    XYZ(2, N) = Y(I1);
    XYZ(3, N) = Z(I1);
    XYZ(4, N) = X(I2);
    XYZ(5, N) = Y(I2);
    XYZ(6, N) = Z(I2);
    XYZ(7, N) = X(I3);
    XYZ(8, N) = Y(I3);
    XYZ(9, N) = Z(I3);
    XYZ(10, N) = X(I4);
    XYZ(11, N) = Y(I4);
    XYZ(12, N) = Z(I4);
    MATP(N) = MTYPE;
    
    fprintf(IOUT, '%10d      %10d    %10d    %10d   %10d       %5d\n', N, I1, I2, I3, I4, MTYPE);

%   Compute connectivity matrix
    LM(1, N) = ID(3, I1);
    LM(4, N) = ID(3, I2);
    LM(7, N) = ID(3, I3);
    LM(10, N) = ID(3, I4);

    LM(2, N) = ID(1, I1);
    LM(5, N) = ID(1, I2);
    LM(8, N) = ID(1, I3);
    LM(11, N) = ID(1, I4);

    LM(3, N) = ID(2, I1);
    LM(6, N) = ID(2, I2);
    LM(9, N) = ID(2, I3);
    LM(12, N) = ID(2, I4);

%   Updata column heights and bandwidth
    ColHt(LM(:, N))
end
sdata.XYZ = XYZ; sdata.MATP = MATP; sdata.LM = LM;

% Clear the memory of X, Y, Z
sdata.X = double(0);
sdata.Y = double(0);
sdata.Z = double(0);

end



