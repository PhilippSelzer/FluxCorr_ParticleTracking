function varargout = read_HGS_binary_porosity(dimP,no_elem,filenameP)

fd = fopen(filenameP,'r','a');

fseek(fd,8*11 + 4,'bof'); 
skipbytes = 0;
data = fread(fd,[dimP,no_elem],'float32',skipbytes,'a')';
data_neu = data;
varargout{1}=data_neu;


fclose(fd);

end