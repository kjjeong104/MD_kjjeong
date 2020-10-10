import os
import numpy as np
avg_list = []
std_list = []
for file in os.listdir("./"):
	if file.startswith("pot_"):
		print("reading {}".format(str(file)))
		data = np.loadtxt(file)
		print(str(file))
		n_data = len(data)
		avg_list.append(np.mean(data[int(n_data/2):]))
		std_list.append(np.std(data[int(n_data/2):]))
		
print(avg_list)
print(std_list)

