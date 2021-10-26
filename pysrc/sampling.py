filename = "Data/0V_1.0g.csv" 
efilename = "Data/sampled_0V_1.0g.csv" 

file = open(filename,'r')
ofile = open(efilename,'w')

i=0;

for line in file:
	if(i%10 == 0):
		ofile.write(line)
	i+=1

ofile.close()