import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt

roi_path = 'egfr.roi.txt'
path_dir = 'out'
file_list = os.listdir(path_dir)

roifile = open(roi_path, "r")

lines = roifile.readlines()

for i, line in enumerate(lines) :
	first = line.find(':')
	second = line.find('-')
	contig = line[:first]
	start = int(line[first+1:second])
	stop = int(line[second+1:])
	coverage = [[] for i in range(stop-start+1)]

	for file_name in file_list :
		outlist = list(map(int, open(path_dir+"/"+file_name, "r").readlines()[i].split(',')))
		for j, out in enumerate(outlist) :
			coverage[j].append(out)
			print('\r', file_name+", roi"+str(i+1), end='')
	
	print('\n'+"roi"+str(i+1)+" saving...")
	fig = plt.figure()
	xticks = np.arange(start, stop+1)
	df = pd.DataFrame(list(map(list, zip(*coverage))))#, columns=xticks)
	boxplot = df.boxplot()
	plt.savefig(roi_path+"_"+"roi"+str(i+1)+"_boxplot"+'.png')
	plt.close(fig)
	print(roi_path+"_"+"roi"+str(i+1)+"_boxplot"+'.png saved!')

	fig2 = plt.figure()
	df2 = pd.DataFrame(coverage, index=xticks, columns=None)
	df2.plot(color='black', alpha=0.1)
	plt.legend().remove()
	plt.savefig(roi_path+"_"+"roi"+str(i+1)+"_lineplot"+'.png')
	plt.close(fig)
	print(roi_path+"_"+"roi"+str(i+1)+"_lineplot"+'.png saved!')

roifile.close()

'''
# box plot test
df = pd.DataFrame(samlist, columns = [columnname])
plt.figure(figsize=(7, 6))
boxplot = df.boxplot(column = [columnname])
plt.yticks(np.arange(min(samlist), max(samlist), step=(max(samlist)-min(samlist))/20))

# line plot test
cvindex = np.arange(start-span, stop+span+1)
df2 = pd.DataFrame(cv, index = cvindex)
lines = df2.plot.line()



plt.show()
'''
