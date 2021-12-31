# Input file format: /data/tusers/liying/project/CTCF/GM12878/ENCFF356LIU_AS_reads.bed


import random
import argparse
import copy
from collections import defaultdict

# add command parameters
parser = argparse.ArgumentParser()
parser.add_argument("-AS", help = "please give an AS bed file path, eg. ENCFF758RQJ_AS.bed", required = True)
parser.add_argument("-nonAS", help = "please give an non-AS bed file path, eg. ENCFF758RQJ_nonAS.bed", required = True)
parser.add_argument("-t", help = "simulation times you want", required = True)
args = parser.parse_args()

AS_file = args.AS
non_AS_file = args.nonAS
times = int(args.t)

path = '/'.join(AS_file.split('/')[:-1])
name = AS_file.split('/')[-1].split('_')[0]

distances = []
mafs = []
with open(AS_file) as f0:
    for line in f0:
        line = line.strip()
        pos = int(line.split('\t')[1])
        alt = line.split('\t')[4]
        ref = line.split('\t')[3]
        # summit = int(line.split('\t')[8])
        # tss = int(line.split('\t')[10])
        # dis = abs(tss - summit)

        phased_alt = line.split('\t')[14].split(',')
        phased_ref = line.split('\t')[13]

        if alt not in phased_alt:
            continue

        if ref != phased_ref:
            continue

        dis = int(line.split('\t')[9])
        distances.append(dis)

        maf_index = phased_alt.index(alt)
        # get alt allele's matched MAF；
        maf = float(line.split('\t')[-2].split(';')[1].split('=')[1].split(',')[maf_index])

        if maf > 0.5:
            maf = 1 - maf
        mafs.append(maf)

# sort distances of peak summit and MAFs
distances.sort()
mafs.sort()

dis1 = distances[int(len(distances) * (1/4) - 1)]
dis2 = distances[int(len(distances) * (2/4) - 1)]
dis3 = distances[int(len(distances) * (3/4) - 1)]
dis4 = distances[int(len(distances) * (4/4) - 1)]

maf1 = mafs[int(len(mafs) * (1/4) - 1)]
maf2 = mafs[int(len(mafs) * (2/4) - 1)]
maf3 = mafs[int(len(mafs) * (3/4) - 1)]
maf4 = mafs[int(len(mafs) * (4/4) - 1)]



# make 16 AS bins
bins = [[] for i in range(16)]

# classify AS SNPs to 16 groups
with open(AS_file) as f:
    for line in f:
        line = line.strip()
        pos = int(line.split('\t')[1])

        alt = line.split('\t')[4]
        ref = line.split('\t')[3]
        phased_alt = line.split('\t')[14].split(',')
        phased_ref = line.split('\t')[13]

        if alt not in phased_alt:
            continue

        if ref != phased_ref:
            continue

        maf_index = phased_alt.index(alt)
        # get alt allele's matched MAF；
        maf = float(line.split('\t')[-2].split(';')[1].split('=')[1].split(',')[maf_index])

        if maf > 0.5:
            maf = 1 - maf


        dis = int(line.split('\t')[9])

        if dis <= dis1:
            if maf <= maf1:
                bins[0].append(line)
            elif maf1 < maf <= maf2:
                bins[1].append(line)
            elif maf2 < maf <= maf3:
                bins[2].append(line)
            else:
                bins[3].append(line)
        elif dis1 < dis <= dis2:
            if maf <= maf1:
                bins[4].append(line)
            elif maf1 < maf <= maf2:
                bins[5].append(line)
            elif maf2 < maf <= maf3:
                bins[6].append(line)
            else:
                bins[7].append(line)
        elif dis2 < dis <= dis3:
            if maf <= maf1:
                    bins[8].append(line)
            elif maf1 < maf <= maf2:
                bins[9].append(line)
            elif maf2 < maf <= maf3:
                bins[10].append(line)
            else:
                bins[11].append(line)
        else:
            if maf <= maf1:
                bins[12].append(line)
            elif maf1 < maf <= maf2:
                bins[13].append(line)
            elif maf2 < maf <= maf3:
                bins[14].append(line)
            else:
                bins[15].append(line)


# classify non-AS SNPs in 16 groups
non_bins = [[] for j in range(16)]

# classify AS SNPs to 16 groups
with open(non_AS_file) as f:
    for line in f:
        line = line.strip()
        pos = int(line.split('\t')[1])

        alt = line.split('\t')[4]
        ref = line.split('\t')[3]
        phased_alt = line.split('\t')[14].split(',')
        phased_ref = line.split('\t')[13]


        if alt not in phased_alt:
            continue

        if ref != phased_ref:
            continue

        maf_index = phased_alt.index(alt)

        # get alt allele's matched MAF；
        maf = float(line.split('\t')[-2].split(';')[1].split('=')[1].split(',')[maf_index])

        if maf > 0.5:
            maf = 1 - maf

        dis = int(line.split('\t')[9])

        if dis <= dis1:
            if maf <= maf1:
                non_bins[0].append(line)
            elif maf1 < maf <= maf2:
                non_bins[1].append(line)
            elif maf2 < maf <= maf3:
                non_bins[2].append(line)
            else:
                non_bins[3].append(line)
        elif dis1 < dis <= dis2:
            if maf <= maf1:
                non_bins[4].append(line)
            elif maf1 < maf <= maf2:
                non_bins[5].append(line)
            elif maf2 < maf <= maf3:
                non_bins[6].append(line)
            else:
                non_bins[7].append(line)
        elif dis2 < dis <= dis3:
            if maf <= maf1:
                    non_bins[8].append(line)
            elif maf1 < maf <= maf2:
                non_bins[9].append(line)
            elif maf2 < maf <= maf3:
                non_bins[10].append(line)
            else:
                non_bins[11].append(line)
        else:
            if maf <= maf1:
                non_bins[12].append(line)
            elif maf1 < maf <= maf2:
                non_bins[13].append(line)
            elif maf2 < maf <= maf3:
                non_bins[14].append(line)
            else:
                non_bins[15].append(line)


# with open(path + '/' + name + '_global_binscount.txt', 'w') as count_file:
#     count_file.write('\t'.join([str(i) for i in [len(b) for b in bins]]) + '\n')
#     count_file.write('\t'.join([str(i) for i in [len(b) for b in non_bins]]))

for i in range(times):
    new_random_name = path + '/../controls/' + name + '/' + name + '_random_global_' + str(i) + '.txt'
    # new = open('/data/tusers/liying/project/CTCF/GM12878/iterations/' + name + '_random_1.txt', 'w')
    new = open(new_random_name, 'w')

    # generate controls

    for i in range(16):
        as_reads = []
        nonas_reads = []

        as_sample = defaultdict(list)  # 一个read对应所有是这个read的AS信息
        nonas_sample = defaultdict(list)

        intervals_AS = defaultdict(list) # 格式为{'20-24': [sample1的信息，sample2的信息...]，...}
        intervals_nonAS = defaultdict(list)

        for j in bins[i]:
            as_reads.append(j.split('\t')[-1])
            as_sample[j.split('\t')[-1]].append(j)

        for k in non_bins[i]:
            nonas_reads.append(k.split('\t')[-1])
            nonas_sample[k.split('\t')[-1]].append(k)

        if as_reads != []:
            # as_count = Counter(as_reads) # 统计不同的reads数量
            # nonas_count = Counter(nonas_reads)

            AS_list = ['{0}-{1}'.format(x, x) for x in range(6,20)]
            AS_list = AS_list + ['{0}-{1}'.format(x, x + 4) for x in range(20,50,5)]
            AS_list = AS_list + ['{0}-{1}'.format(x, x + 9) for x in range(50,200,10)]
            AS_list = AS_list + ['{0}-{1}'.format(x, x + 19) for x in range(200, max([int(x) for x in as_reads]) + 1, 20)]

            nonAS_list = ['{0}-{1}'.format(x, x) for x in range(6,20)]
            nonAS_list = nonAS_list + ['{0}-{1}'.format(x, x + 4) for x in range(20,50,5)]
            nonAS_list = nonAS_list + ['{0}-{1}'.format(x, x + 9) for x in range(50,200,10)]
            nonAS_list = nonAS_list + ['{0}-{1}'.format(x, x + 19) for x in range(200, max([int(x) for x in nonas_reads]) + 1, 20)]

            for as_inter in AS_list: # AS每个区间计数
                start = as_inter.split('-')[0]
                end = as_inter.split('-')[1]
                for as_num in as_sample:
                    if int(start) <= int(as_num) <= int(end):
                        for tmp_as in as_sample[as_num]:
                            intervals_AS[as_inter].append(tmp_as)

            for nonas_inter in nonAS_list: # nonAS每个区间计数
                start = nonas_inter.split('-')[0]
                end = nonas_inter.split('-')[1]
                for nonas_num in nonas_sample:
                    if int(start) <= int(nonas_num) <= int(end):
                        for tmp_nonas in nonas_sample[nonas_num]:
                            intervals_nonAS[nonas_inter].append(tmp_nonas)


            for as_inter in intervals_AS:
                second_condition = [] # 第二种情况，本区间有没有匹配上但是往后退可以找到的
                third_condition = [] # 第三种情况，也就是在本区间和往后推都match不上的情况
                as_mess = intervals_AS[as_inter] # as_mess是一个read在某个区间内的所有AS的信息的列表
                as_read = len(as_mess)

                nonas_mess = intervals_nonAS[as_inter]
                nonas_read = len(nonas_mess)

                flag = False

                if as_read <= nonas_read: # 匹配上的
                    # new.write(as_inter + '\t' + as_inter + '\t' + str(as_read) + '\t' + str(nonas_read) + '\n')
                    random_match1 = random.sample(nonas_mess, as_read) # 对于这个read区间的从nonAS中随机挑取与AS同样数量的
                    for match1 in random_match1:
                        new.write(match1 + '\n')
                        intervals_nonAS[as_inter].remove(match1)

                else: # 没有匹配上的

                    if int(as_inter.split('-')[0]) > max([int(x) for x in nonas_reads]): # 判断这个AS的范围是不是比nonAS的都大
                        # 对于这种直接往前取, 作为第三种情况 最后判断
                        # new.write(as_inter + '\t' + 'NA' + '\t' + str(as_read) + '\t0' + '\n')
                        pass

                    else: # 对于没有匹配上的
                        second_condition = second_condition + nonas_mess # 可以在本区间match的信息
                        for tmp_match in nonas_mess:
                            new.write(tmp_match + '\n')
                        index_inter = nonAS_list.index(as_inter) # 取到这个区间的索引

                        # num = 0

                        for z in range(index_inter + 1, len(nonAS_list)): # 往后推
                            if nonAS_list[z] in intervals_nonAS:
                                # num += len(intervals_nonAS[nonAS_list[z]])
                                second_condition = second_condition + intervals_nonAS[nonAS_list[z]] # 将后边这个区间的信息加入

                                # 现在这个区间中能匹配几个 减去
                                if len(second_condition) >= as_read:
                                    previous = len(second_condition) - len(intervals_nonAS[nonAS_list[z]]) # 上一次总共多少个
                                    diff = as_read - previous # 还需要多少个匹配
                                    random_match2 = random.sample(intervals_nonAS[nonAS_list[z]], diff)
                                    for match2 in random_match2:
                                        new.write(match2 + '\n')
                                        intervals_nonAS[nonAS_list[z]].remove(match2)

                                    flag = True # 当满足第二个条件后flag才为true，否则是flase，进去第三种情况
                                    break
                                else:
                                    for tmp_match in intervals_nonAS[nonAS_list[z]]:
                                        new.write(tmp_match + '\n')
                                    intervals_nonAS[nonAS_list[z]] = []

                    if flag == False: # 匹配到最后仍然不够，转到第三种情况
                        third_condition = third_condition + second_condition

                       # 对于往后匹配仍然没有匹配上的往前推

                        if as_inter in nonAS_list:
                            index_third = nonAS_list.index(as_inter)
                        else:
                            index_third = len(nonAS_list)

                        for x in range(index_third - 1, -1, -1): # 全都取到的区间是都取出来的，最后一次是随机的
                            if nonAS_list[x] in intervals_nonAS:
                                third_condition = third_condition + intervals_nonAS[nonAS_list[x]]

                                if len(third_condition) >= as_read:
                                    flag2 = True
                                    previous = len(third_condition) - len(intervals_nonAS[nonAS_list[x]]) # 上一次总共多少个
                                    diff = as_read - previous # 还需要多少个匹配
                                    random_match3 = random.sample(intervals_nonAS[nonAS_list[x]], diff)
                                    for match3 in random_match3:
                                        new.write(match3 + '\n')
                                        intervals_nonAS[nonAS_list[x]].remove(match3)
                                    break
                                else:
                                    for tmp_match in intervals_nonAS[nonAS_list[x]]:
                                        new.write(tmp_match + '\n')
                                    intervals_nonAS[nonAS_list[x]] = [] # 取空


    new.close()