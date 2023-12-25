function AddR(r, LM)

% Get global data
global sdata;
f = sdata.f;
[n,~] = size(f);
R= sdata.R;
ND = sdata.NDOF * sdata.NNODE;
for I = 1:n
    for K = 1:ND
        II = LM(K);
        if (II > 0) 
            R(II, 1) = R(II, 1) + r(K,1); 
        end
    end
end

sdata.R = R;
end