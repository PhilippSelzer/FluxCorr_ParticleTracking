function varargout =read_HGS_binary_exchangefluxes(dim,no_nodes,filename3)

fd = fopen(filename3,'r','a');

fseek(fd,92,'bof');
skipbytes = 0;
data = fread(fd,[dim,no_nodes],'float32',skipbytes,'a')';
data_neu = data;
varargout{1}=data_neu;


fclose(fd);

end