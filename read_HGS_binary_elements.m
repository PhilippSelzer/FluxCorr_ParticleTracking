function varargout = read_HGS_binary_elements(no_nodes_per_element,no_elements,filename1)

fd = fopen(filename1,'r','l');

fseek(fd,28,'bof'); 
skipbytes = 0;
data = fread(fd,[no_nodes_per_element,no_elements],'int32',skipbytes,'l')';
varargout{1}=data;


fclose(fd);

end