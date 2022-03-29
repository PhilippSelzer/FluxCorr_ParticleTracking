function varargout =read_HGS_binary_rain(dim,no_nodes,filename3)

fd = fopen(filename3,'r','a');

fseek(fd,0 +4 + 11*8,'bof'); 
skipbytes = 0;
data = fread(fd,[dim,no_nodes],'float64',skipbytes,'a')';
data_neu = data;
varargout{1}=data_neu;


fclose(fd);

end