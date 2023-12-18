function dealoutfile()

global cdata
global sdata
% infile=fopen('STAPMAT.OUT',"r");
% fprintf('Open in file\n');
% %读取数据
% fclose(infile);
% fprintf('Close in file\n');

outfile=fopen('.\DATA\stapmatout.dat',"w");
fprintf('Open dat file\n');
fprintf(outfile,'Title="sample finite-element\n');
%fprintf(outfile,'Variables="X","Y","Z", "ux", "uy", "uz", "theta x", "theta y", "theta z"\n');
fprintf(outfile,'Variables="X","Y","Z"\n');
fprintf(outfile,'ZONE T="DISPLACEMENT", N=%d, E=%d, F=FEPOINT, ET=BRICK\n',cdata.NUMNP,cdata.NUMEG);%N和E瞎给的

% D = zeros(3, 1, 'double');
% ID=sdata.ID;
% DIS=sdata.DIS;
% for II = 1:cdata.NUMNP
%     D(:) = 0;
%     if (ID(1, II) ~= 0) D(1) = DIS(ID(1, II)); end
%     if (ID(2, II) ~= 0) D(2) = DIS(ID(2, II)); end
%     if (ID(3, II) ~= 0) D(3) = DIS(ID(3, II)); end
%     
%     fprintf(outfile, ' %10d        %18.6e%18.6e%18.6e\n', II, D(1), D(2), D(3));
% end

X = sdata.X; Y = sdata.Y; Z = sdata.Z;
for i = 1:cdata.NUMNP
    fprintf(outfile, '%13.3f%13.3f%13.3f\n', ...
        X(i), Y(i), Z(i));
end


fclose(outfile);
fprintf('Close dat file\n');


end