import os, sys
import pysam, argparse
import numpy as np
import copy

bam_dir = '/data/allbam'
btp_dir = bam_dir+'/190916_work'

#bam_list = os.listdir(bam_dir)
#bam_list = [file for file in bam_list if file.endswith(".bam")]
btp_list = os.listdir(btp_dir)
btp_list = [file for file in btp_list if file.endswith(".bed")]

output_50x = []
output_100x = []

for btp in btp_list :

	btpname = btp[:btp.find('.')]
	print(btpname+' start!')

	output_50x.append({})
	output_100x.append({})
	btp_path = btp_dir+'/'+btp

	slfile = open(btp_dir+'/sample.'+btpname, 'r')
	sl_ls = slfile.readlines()
	bam_list = []
	
	for sl_l in sl_ls :
		bam_list.append(sl_l[:-1]+'.bwamem.sorted.dedup.realn.recal.dedup.bam')

	for bnum, bam in enumerate(bam_list) :

#		if bnum > 2 :
#			break

		sam_path = bam_dir+'/'+bam
		
		if not os.path.isfile(sam_path) :
			continue

		samfile = pysam.AlignmentFile(sam_path, "rb")
		btpfile = open(btp_path, 'r')
		lines = btpfile.readlines()

		bamname = bam[:bam.find('.')]
		output_50x[-1][bamname] = {}
		output_100x[-1][bamname] = {}

		for lnum, line in enumerate(lines) :
	
			splited = line.split('\t')
			contig = splited[0]
			start = int(splited[1])
			stop = int(splited[2])
			gene = splited[3]
			gene = gene[:-1]

			print('\r', bamname+', '+gene+' : '+str(bnum+1)+' / '+str(len(bam_list))+', '+str(lnum+1)+' / '+str(len(lines)), end='\r')
			sys.stdout.write("\033[K")

			cv_original = np.array(samfile.count_coverage(contig, start=start, stop=stop))
			cv_50x = cv_original.sum(axis=0)
			cv_100x = copy.deepcopy(cv_50x)

			cv_50x = [1 if ce > 50 else 0 for ce in cv_50x]
			cv_100x = [1 if ce > 100 else 0 for ce in cv_100x]

			num_50x = sum(cv_50x)
			num_100x = sum(cv_100x)
			
			if gene in output_50x[-1][bamname] :
				output_50x[-1][bamname][gene]['num'] += num_50x
				output_50x[-1][bamname][gene]['len'] += len(cv_50x)
			else :
				output_50x[-1][bamname][gene] = {}
				output_50x[-1][bamname][gene]['num'] = num_50x
				output_50x[-1][bamname][gene]['len'] = len(cv_50x)

			if gene in output_100x[-1][bamname] :
				output_100x[-1][bamname][gene]['num'] += num_100x
				output_100x[-1][bamname][gene]['len'] += len(cv_100x)
			else :
				output_100x[-1][bamname][gene] = {}
				output_100x[-1][bamname][gene]['num'] = num_100x
				output_100x[-1][bamname][gene]['len'] = len(cv_100x)

		btpfile.close()
		samfile.close()

	outfile_50x = open(btp_dir+'/'+btpname+'.50x.tsv', 'w')
	outfile_100x = open(btp_dir+'/'+btpname+'.100x.tsv', 'w')

	key_50x = output_50x[-1].keys()
	key_100x = output_50x[-1].keys()

	gene_50x = output_50x[-1][list(key_50x)[0]].keys()
	gene_100x = output_100x[-1][list(key_100x)[0]].keys()

	outfile_50x.write('GENE')
	for bamn in key_50x :
		outfile_50x.write('\t'+bamn)

	outfile_50x.write('\n')
	for genn in gene_50x :
		outfile_50x.write(genn)
		for bamn in key_50x :
			cnum = output_50x[-1][bamn][genn]['num'] / output_50x[-1][bamn][genn]['len']
			outfile_50x.write('\t'+str(cnum))
		outfile_50x.write('\n')
	
	outfile_50x.close()


	outfile_100x.write('GENE')
	for bamn in key_100x :
		outfile_100x.write('\t'+bamn)

	outfile_100x.write('\n')
	for genn in gene_100x :
		outfile_100x.write(genn)
		for bamn in key_100x :
			cnum = output_100x[-1][bamn][genn]['num'] / output_100x[-1][bamn][genn]['len']
			outfile_100x.write('\t'+str(cnum))
		outfile_100x.write('\n')

	outfile_100x.close()

