import os

id_list = open('gdcs', 'r')
id_lines = id_list.readlines()
for line in id_lines :
	os.system('/data2/gdc/gdc-client download '+line[:-1]+' -t /data2/gdc/gdc_token.txt')

