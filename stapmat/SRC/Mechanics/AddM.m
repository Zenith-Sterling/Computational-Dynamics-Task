function AddM(M1,M2, LM)

% Get global data
global sdata;
MAXA = sdata.MAXA; 
Mass1= sdata.Mass1;
Mass2= sdata.Mass2;
ND = sdata.NDOF * sdata.NNODE;
for J = 1:ND
    JJ = LM(J);
    if (JJ > 0)
        for I = 1:J
            II = LM(I);
            if (II > 0)
                if (JJ > II) KK = MAXA(JJ) + JJ - II;
                else KK = MAXA(II) + II - JJ; end
                Mass1(KK) = Mass1(KK) + M1(I, J);
                Mass2(KK) = Mass2(KK) + M2(I, J);
            end
        end
    end
end

sdata.Mass1 = Mass1;
sdata.Mass2 = Mass2;

end