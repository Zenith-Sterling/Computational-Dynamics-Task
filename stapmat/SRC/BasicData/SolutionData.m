%* *****************************************************************
%* - Basic data class of STAPMAT                                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Storing variables used in solving process                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.20                *
%*                                                                 *
%* *****************************************************************

classdef SolutionData
    properties (Constant)
        % Gauss coord, 1D to 3D
        GC1 = double(0.0);
        GC2 = double([1/sqrt(3),-1/sqrt(3)]);
        GC3 = double([sqrt(0.6), 0.0, -sqrt(0.6)]);
        % Gauss weight, 1D to 3D
        GW1 = double(2.0);
        GW2 = double([1.0, 1.0]);
        GW3 = double([5.0/9.0, 8.0/9.0, 5.0/9.0]);
    end
    properties
        % Basic data
        ID;       % int, ID(3, NUMNP), Boundary condition codes (0=free, 1=deleted)
        IDOrigin; % int, backups of ID after computing of NEQ
        X;        % double, X(NUMNP), X coordinates
        Y;        % double, Y(NUMNP), Y coordinates
        Z;        % double, Z(NUMNP), Z coordinates
        R;        % double, R(NEQ), Load vector
        NOD;      % int, NOD(NLOAD), Node number to which this load is applied (1~NUMNP)
        f;        % double, f(NLCASE,9)=[开始单元,结束单元,步长,x方向力值,y方向力值,z方向力值,x方向力矩,y方向力矩,z方向力矩], Load vector
        ELE_1;      % int, ELE(NLOAD), Element number to which this load is applied (1~NUMEG)
        ELE_2;
        IDIRN;    % int, IDIRN(NLOAD), Degree of freedom number for this load component
                  %                     1 : X-direction;
                  %                     2 : Y-direction;
                  %                     3 : Z-direction;
        FLOAD;    % double, FLOAD(NLOAD), Magnitude of load

        
        
        % Element data
        NUME;     % int, number of elements
        NNODE;    % int, number of nodes in an element
        NINIP;    % int, number of integration points in an element
        NDOF;     % int, the DOF of displacement
        NSTIFF;   % int, the number of number in element stiffness matrix
        XYZ;      % double, XYZ(3*NNODE, NUME), element position
        
        InitCoord;  % double array, integration coordinates
        InitWeight; % double array, integration weights
        
        % Material data
        NUMMAT;     % int, the number of types of material 
        E;          % double array, Young's Modulus
        nu;         % double array, possion ratio
        AREA;       % double array, cross-sectional constants
        rho;
        thick;
        MATP;       % int, MATP(NUME), types of elements
        
        % Solve data
        NEQ;      % int, Number of equations
        NWK;      % Number of matrix elements
        MK;       % Maximum half bandwidth
        MHT;      % int, MHT(NEQ), Vector of column heights
        LM;       % int, LM(6, NUME), Connectivity matrix
        MAXA;     % int, MAXA(NEQ)
        STIFF;    % double ,STIFF(NWK), store the elements of stiffness matrix
        
        % Result data
        DIS;      % double, DIS(NEQ, NLCASE), Displacement of nodes
        STRAIN;   % double, STRAIN(NEQ, NLCASE), Strain
        STRESS;   % double, STRESS(NEQ, NLCASE), Stress
        Mass1;
        Mass2;
        
    end
end