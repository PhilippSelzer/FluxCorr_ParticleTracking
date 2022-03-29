function varargout = read_HGS_binary_krels(no_nodes_per_element,no_elements,filename1)

fd = fopen(filename1,'r','a');

fseek(fd,4,'bof'); 
skipbytes = 0;
data = fread(fd,[no_nodes_per_element,no_elements],'float64',skipbytes,'a')';
varargout{1}=data;


fclose(fd);

end