import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec
from operator import add
import pysam, argparse
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, silhouette_samples
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.metrics import auc


'''
Argument Setting
'''

# 파일 이름과 span을 argument로 불러들인다.
parser = argparse.ArgumentParser()
parser.add_argument('bam_dir_path', help='bam파일을 읽어들일 디렉토리를 정합니다.')
parser.add_argument('refseq_path', help='refseq 파일의 경로를 지정합니다.')
parser.add_argument('nmid_to_draw', help='사용할 NMID를 지정합니다.')
parser.add_argument('exon_num', type=int, help='exon number를 지정합니다.(1부터 시작)')
parser.add_argument('output_prefix', help='output 파일명을 정합니다.')
parser.add_argument('--color_h', default='b', help='영역 high 색을 정합니다.')
parser.add_argument('--color_l', default='g', help='영역 low 색을 정합니다.')
parser.add_argument('--alpha', type=float, default=0.5, help='영역 투명도를 정합니다.(0~1)')

args = parser.parse_args()
bam_dir = args.bam_dir_path
refseq_path = args.refseq_path
nmid_to_draw = args.nmid_to_draw
output_prefix = args.output_prefix

color_h = args.color_h
color_l = args.color_l
alpha = args.alpha
exon_num = args.exon_num-1




'''
Reading Refseq Data
'''

# Refseq 불러오기
print('reading refseq data...')
refseq = pd.read_csv(refseq_path,
	sep='\t',
	names=['bin', 'name', 'chrom', 'strand',
		'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
		'exonCount', 'exonStarts', 'exonEnds', 'score',
		'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
	)

# 사용자가 찾으려 하는 nmid의 refseq 정보만 불러온다.
refseq_nm = refseq[refseq.name.str.contains(nmid_to_draw)]

name_nm = refseq_nm.name.tolist()
chrom_nm = refseq_nm.chrom.tolist()
tx_s_nm = refseq_nm.txStart.tolist()
tx_e_nm = refseq_nm.txEnd.tolist()
exon_s_nm = refseq_nm.exonStarts.tolist()
exon_e_nm = refseq_nm.exonEnds.tolist()

# 그 유전자를 토대로 그릴 범위를 정한다.
name_n = name_nm[0]
contig = chrom_nm[0]
start = tx_s_nm[0]
stop = tx_e_nm[0]
exon_s = list(map(int, exon_s_nm[0][:-1].split(',')))
exon_e = list(map(int, exon_e_nm[0][:-1].split(',')))










'''
Bam Information Analysis
'''

# coverage를 저장할 데이터를 초기화한다.
normal_coverage = np.zeros(stop-start+1)
coverage = [[] for i in range(stop-start+1)]
samfile = None


# Cancer Bam
print('\nanalyzing cancer bam information...')

# sample list가 들어있는 파일을 불러온다.
bamd_list = os.listdir(bam_dir)
bamd_list = [file for file in bamd_list if os.path.isdir(bam_dir+'/'+file)]
bam_list = []
for ddd in bamd_list :
	ffff = os.listdir(bam_dir+'/'+ddd)
	ffff = [ddd+'/'+file for file in ffff if file.endswith(".bam")]
	bam_list.extend(ffff)


for bamn, bam in enumerate(bam_list) :

#	if sdfs :
#		break

#	if bamn >= 20 :
#		break

	sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': start...')
	sys.stderr.write("\033[K")
	# Cancer Bam 경로
	sam_path = bam_dir+'/'+bam	

	# Cancer Bam이 파일인지 확인
	if not os.path.isfile(sam_path) :
		continue
	

	sdfs = True


	# cache path
	bam_fn = os.path.basename(bam)
	cache_path = 'cache/'+bam_fn+'_'+name_n+'_'+contig+'_'+str(start)+'_'+str(stop)
	
	cv = []

	if not os.path.isfile(cache_path+'.npy') :

		# Bam파일을 불러온다.
		sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': loading bam file...')
		sys.stderr.write("\033[K")
		samfile = pysam.AlignmentFile(sam_path, "rb")

		# Coverage 계산
		sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': coverage calculating...')
		sys.stderr.write("\033[K")
		assdfe = samfile.count_coverage(contig, start=start, stop=stop+1)

		sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': converting coverage to array...')
		sys.stderr.write("\033[K")
		cv_original = np.array(assdfe)
		samfile.close()

		# Coverage가  A, T, G, C 4성분으로 다 따로 출력되기 때문에,
		# 이를 합쳐주는 작업을 한다.
		sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': coverage adding...')
		sys.stderr.write("\033[K")
		cv = cv_original.sum(axis=0) 

		# caching
		sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': caching...')
		sys.stderr.write("\033[K")
		cv_np = np.array(cv)
		np.save(cache_path, cv_np)
	else :
		sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': load cache...')
		sys.stderr.write("\033[K")
		cv = list(np.load(cache_path+'.npy'))

	# flag2 속성이 켜져 있으면, coverage의 최댓값으로
	# 나머지 값을 나누어주는 작업을 한다.
	# 전체 coverage의 최댓값
	maxcv = max(cv)	
	# 모든 위치에 대해 최댓값으로 나누어 준다.
	for j, out in enumerate(cv) :
		coverage[j].append(out/maxcv)




	





'''
Draw Lineplot
'''

st = exon_s[exon_num]
stop_n = exon_e[exon_num]+1


# x축의 값을 정해진 시작점과 끝점으로 한다.
xticks = np.arange(st, stop_n)

# test dataframe
# df = 
# coverage의 dataframe을 만들고 plot의 윗 부분을 불러온다.
df = pd.DataFrame(coverage[st-start:stop_n-start], index=xticks, columns=None)


'''
Clustering
'''

df_T = df.T
data_lines = df_T.values
kmeans = KMeans(n_clusters=2).fit(data_lines)

df_T['cluster_id'] = kmeans.labels_
df_c1 = df_T[df_T.cluster_id == 0].drop('cluster_id', 1).T
df_c2 = df_T[df_T.cluster_id == 1].drop('cluster_id', 1).T	

df_c1_mean = np.mean(df_c1.values)
df_c2_mean = np.mean(df_c2.values)

labels = []

if df_c1_mean < df_c2_mean :
	labels = [1 if lb == 0 else 0 for lb in kmeans.labels_]
else :
	labels = kmeans.labels_

print('\n\n\n')
print('exon 2 h : '+str(len([1 for lb in labels if lb == 0])))
print('exon 2 l : '+str(len([1 for lb in labels if lb == 1])))




for an in range(0, len(exon_s)) :

	print('\n\n---------------------------------------------------')
	print("exon "+str(exon_num+1)+" and exon "+str(an+1))
	print('---------------------------------------------------\n\n')

	'''
	Another Exon
	'''


	st2 = exon_s[an]
	stop_n2 = exon_e[an]+1


	# x축의 값을 정해진 시작점과 끝점으로 한다.
	xticks2 = np.arange(st2, stop_n2)

	if len(xticks) > len(xticks2) :
		continue

	df2 = pd.DataFrame(coverage[st2-start:stop_n2-start], index=xticks2, columns=None)

	df2_T = df2.T
	data_lines2 = df2_T.values
	kmeans2 = KMeans(n_clusters=2).fit(data_lines2)

	df2_T['cluster_id'] = kmeans2.labels_
	df2_c1 = df2_T[df2_T.cluster_id == 0].drop('cluster_id', 1).T
	df2_c2 = df2_T[df2_T.cluster_id == 1].drop('cluster_id', 1).T	

	df2_c1_mean = np.mean(df2_c1.values)
	df2_c2_mean = np.mean(df2_c2.values)

	labels2 = []

	if df2_c1_mean < df2_c2_mean :
		labels2 = [1 if lb == 0 else 0 for lb in kmeans2.labels_]
	else :
		labels2 = kmeans2.labels_






	#X_train, X_test, Y_train, Y_test = train_test_split(df_T, labels, test_size=0.3)


	Y_train = labels
	X_train = df_T
	#Y_test = labels
	#X_test = df_T
	Y_test = labels2
	X_test = df2_T.iloc[:, :len(df_T.columns)]


	Y_l = []

	log_clf = LogisticRegression()
	log_clf.fit(X_train,Y_train)
	print(log_clf.score(X_test, Y_test))

	y_predicted = log_clf.predict(X_test)
	y_hat = [xt[1] for xt in log_clf.predict_proba(X_test)]

	plt.ylabel('y_hat')
	plt.xlabel('y(true_value)')
	plt.xlim(-0.5, 1.5)
	plt.ylim(-0.1, 1.1)
	plt.xticks([0, 1])
	plt.yticks(np.arange(0, 1.2, 0.2))
	plt.plot(Y_test, y_hat,'o')
	plt.savefig(output_prefix+'_'+str(an+1)+'_.pdf')
	plt.close()


	Y_l = ['h' if lb == 0 else 'l' for lb in kmeans.labels_]



	fpr, tpr, thresholds = roc_curve(Y_test, y_hat)
	print(auc(fpr, tpr))

	plt.plot(fpr, tpr, 'o-')
	plt.plot([0, 1], [0, 1], 'k--')
	plt.ylabel('sensitivity')
	plt.xlabel('1-specificity')
	plt.xlim(-0.1, 1.1)
	plt.ylim(-0.1, 1.1)
	plt.xticks([0, 1])
	plt.yticks(np.arange(0, 1.2, 0.2))
	plt.savefig(output_prefix+'_'+str(an+1)+'_auroc.pdf')
	plt.close()






