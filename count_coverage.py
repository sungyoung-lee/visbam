import pysam, argparse
import numpy as np

# 파일 이름과 span을 argument로 불러들인다.
parser = argparse.ArgumentParser()
parser.add_argument("--bam", help="bam 파일 경로를 지정합니다.")
parser.add_argument("--roi", help="roi 파일 경로를 지정합니다.")
parser.add_argument("--span", help="span 범위를 지정합니다.", type=int, default=0)
parser.add_argument("--out", help="output 파일 경로를 지정합니다.")
args = parser.parse_args()

sam_path = args.bam
roi_path = args.roi
span = args.span
out_path = args.out

# bam 파일과 roi 파일을 읽어들인다.
samfile = pysam.AlignmentFile(sam_path, "rb")
roifile = open(roi_path, 'r')
outfile = open(out_path, 'w')
lines = roifile.readlines()

for line in lines :	
	# 각 roi에 해당하는 bam파일의 coverage를 구하여 출력한다.
	first = line.find(':')
	second = line.find('-')
	contig = line[:first]
	start = int(line[first+1:second])
	stop = int(line[second+1:])

	cv_original = np.array(samfile.count_coverage(contig, start=start-span, stop=stop+span+1))
	cv = cv_original.sum(axis=0) 

	outfile.write(str(cv[0]))
	for coverage in cv[1:] :
		outfile.write(","+str(coverage))
	outfile.write("\n")


outfile.close()
roifile.close()
samfile.close()
