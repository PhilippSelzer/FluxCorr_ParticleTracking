function varargout = read_HGS_binary_coordinates(xyz,no_nodes,filename2)

fd = fopen(filename2,'r','a');

fseek(fd,7+1*8 +1,'bof'); 
skipbytes = 0;
data = fread(fd,[xyz,no_nodes],'float64',skipbytes,'a')';
data_neu = data;
varargout{1}=data_neu;


fclose(fd);

end