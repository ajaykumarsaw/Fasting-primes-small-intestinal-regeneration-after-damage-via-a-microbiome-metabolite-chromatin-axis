import re
import subprocess
import sys

input2 = sys.argv[1]
target_reads = sys.argv[2]
input0 = sys.argv[3]
output0 = sys.argv[4]
log = sys.argv[5]
CASES = sys.argv[2:]

with open (input2, "r") as f:
    # fifth line contains the number of mapped reads
    line = f.readlines()[6]
    match_number = re.match(r'(\d.+) \+.+', line)
    total_reads = int(match_number.group(1))

target_reads = int(target_reads) # 15million reads  by default, set up in the config.yaml file
if total_reads > target_reads:
    down_rate = target_reads/total_reads
else:
    down_rate = 1

subprocess.run("sambamba view -f bam -t 5 --subsampling-seed=3 -s {rate} {inbam} | samtools sort -m 2G -@ 5 -T {outbam}.tmp > {outbam} 2> {log}".format(rate = down_rate, inbam = input0, outbam = output0, log = log), shell=True)
subprocess.run("samtools index {outbam}".format(outbam = output0), shell=True)
