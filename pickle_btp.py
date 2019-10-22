import pickle, sys, os, pysam

sample_list_path = '/data/allbam/190916_work/sample.btp5'
bam_dir = '/data/allbam'

slfile = open(sample_list_path, 'r')
sl_ls = slfile.readlines()
bam_list = []
slfile.close()

for sl_l in sl_ls :
	bam_list.append(sl_l[:-1]+'.bwamem.sorted.dedup.realn.recal.dedup.bam')

for bamn, bam in enumerate(bam_list) :

	print('\r', bam, ':', str(bamn)+'/'+str(len(bam_list)), ': start...', end='')
	sys.stdout.write("\033[K")
	sam_path = bam_dir+'/'+bam	

	if not os.path.isfile(sam_path) :
		continue

	print('\r', bam, ':', str(bamn)+'/'+str(len(bam_list)), ': opening bam file...', end='')
	sf = open(sam_path, 'r')
	print('\r', bam, ':', str(bamn)+'/'+str(len(bam_list)), ': loading bam file...', end='')

	samfile = pysam.AlignmentFile(sf)
	with open('pickle/'+bam+'.p', 'wb') as file:
		pickle.dump(samfile, file)

	samfile.close()
