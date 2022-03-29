function varargout = read_HGS_binary_Ksat(xyz,no_elements,filename2)

fd = fopen(filename2,'r','a');

fseek(fd,4 + 8*11,'bof');
skipbytes = 0;
data = fread(fd,[xyz,no_elements],'float64',skipbytes,'a')';
data_neu = data;
varargout{1}=data_neu;


fclose(fd);

end