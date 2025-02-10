def GetReverseComplement(str):
    if str == 'NA':
        return 
    dict = {}
    dict['A'] = 'T'
    dict['C'] = 'G'
    dict['G'] = 'C'
    dict['T'] = 'A'
    str2 = ''
    for i in range(0, len(str)):
        str2 = str2 + dict[str[i]]
    return str2[::-1]

with open("WT_UV_bothruns_Muts_allSNVs_sorted.bed") as f:  
    data = f.read()
with open("rad16_UV_bothruns_Muts_allSNVs_sorted.bed") as f16:
    data16 = f16.read()
with open("rad30_SNPs_sorted_tandem.bed") as f30:
    data30 = f30.read()
with open("WT_UVA_Muts_allSNVs_sorted.bed") as fA:
    dataUVA = fA.read()
with open("ogg1_UVA_Muts_allSNVs_sorted.bed") as fOg:
    dataOg = fOg.read()
with open("rad26_UV_bothruns_Muts_allSNVs_sorted.bed") as f26:
    data26 = f26.read()
with open("WT_UVB_Muts_allSNVs_sorted.bed") as fB:
    dataUVB = fB.read()
with open("rad16_UVB_Muts_allSNVs_sorted.bed") as fB16:
    dataB16 = fB16.read()
with open('rad26_UVB_Muts_allSNVs_sorted.bed') as fB26:
    dataB26 = fB26.read()
with open('rad30_UVB_Muts_allSNVs_sorted.bed') as fB30:
    dataB30 = fB30.read()
with open("WT_UV_run1_Muts_allSNVs_sorted.bed") as fr1:
    datar1 = fr1.read()
with open("WT_UV_run2_Muts_allSNVs_sorted.bed") as fr2:
    datar2 = fr2.read()
with open("clusters16.bed") as c:
    datac = c.read()
with open("clusters16_UVB.bed") as c1:
    datac1 = c1.read()
with open("clusters.bed") as w:
    dataw = w.read()
with open("clusters_UVB.bed") as w1:
    dataw1 = w1.read()
with open("WT_combined_sorted_tandem.bed") as a1:
    dataa1 = a1.read()
with open("rad16_combined_sorted_tandem.bed") as b1:
    datab1 = b1.read()
with open("rad26_combined_sorted_tandem.bed") as d1:
    datad1 = d1.read()
with open("rad30_combined_sorted_tandem.bed") as e1:
    datae1 = e1.read()
with open("WT_deletions_sorted_tandem.bed") as a2:
    dataa2 = a2.read()
with open("rad16_deletions_sorted_tandem.bed") as b2:
    datab2 = b2.read()
with open("rad26_deletions_sorted_tandem.bed") as d2:
    datad2 = d2.read()
with open("rad30_deletions_sorted_tandem.bed") as e2:
    datae2 = e2.read()
with open("WT_insertions_sorted_tandem.bed") as a3:
    dataa3 = a3.read()
with open("rad16_insertions_sorted_tandem.bed") as b3:
    datab3 = b3.read()
with open("rad26_insertions_sorted_tandem.bed") as d3:
    datad3 = d3.read()
with open("rad30_insertions_sorted_tandem.bed") as e3:
    datae3 = e3.read()

#with open("GSM2109562_0hr_UV2_A1_dipy_bkgd.bed") as a1:
    #dataa1 = a1.read()
#with open("GSM2109560_UV_0hr_A2_dipy_bkgd.bed") as a2:
    #dataa2 = a2.read()
#with open("GSM2109566_UV90J_A3_dipy_bkgd.bed") as a3:
    #dataa3 = a3.read()
#with open("GSM2109565_UV60J_A2_dipy_bkgd.bed") as j60:
    #dataj60 = j60.read()
#with open("2hr_UV2_WT_A3_1bp_dipy_bk.bed") as h2:
    #datah2 = h2.read()
#with open("GSM3763516_rad16_UV0hr_rep2_A4_combined_1bp_dipy_bk.bed") as a4:
    #dataa4 = a4.read()
#with open("GSM3763517_rad16_UV2hr_rep2_A5_combined_1bp_dipy_bk.bed") as a5:
    #dataa5 = a5.read()
#with open("GSM4339207_rad26-0hr-rep2-A1_dipy_bk.bed") as a_r26:
    #dataa_r26 = a_r26.read()
#with open("GSM4339208_rad26-2hr-rep2-A2_dipy_bk.bed") as a3_r26:
    #dataa3_r26 = a3_r26.read()
    
#extract columns from bed files to make arrays
# with matching indices   
def getChromosome(data):       
    lines = [s.strip().split() for s in data.splitlines()]
    chr_counts = {}
    chromosome = []
    for i in range (0, len(lines)):
        if lines[i][0] != "chrXII" or (int(lines[i][2]) < 451430) or (int(lines[i][2]) > 468920 and int(lines[i][2]) < 489928) or (int(lines[i][2]) > 490455):
            chromosome.append(lines[i][0])
            if lines[i][0] not in chr_counts.keys():
                chr_counts[lines[i][0]] = 1
            else:
                chr_counts[lines[i][0]] = chr_counts[lines[i][0]] + 1
    for key in chr_counts.keys():
        print(key + ": " + str(chr_counts[key]))
    return chromosome
chromosome = getChromosome(data)
chromosome16 = getChromosome(data16)
chromosome30 = getChromosome(data30)
chromosomeA = getChromosome(dataUVA)
chromosomeOg = getChromosome(dataOg)
chromosome26 = getChromosome(data26)
chromosomeB = getChromosome(dataUVB)
chromosomeB16 = getChromosome(dataB16)
chromosomeB26 = getChromosome(dataB26)
chromosomeB30 = getChromosome(dataB30)
chromosomer1 = getChromosome(datar1)
chromosomer2 = getChromosome(datar2)
chromosomec = getChromosome(datac)
chromosomec1 = getChromosome(datac1)
chromosomew = getChromosome(dataw)
chromosomew1 = getChromosome(dataw1)
chromosomea1 = getChromosome(dataa1)
chromosomeb1 = getChromosome(datab1)
chromosomed1 = getChromosome(datad1)
chromosomee1 = getChromosome(datae1)
chromosomea2 = getChromosome(dataa2)
chromosomeb2 = getChromosome(datab2)
chromosomed2 = getChromosome(datad2)
chromosomee2 = getChromosome(datae2)
chromosomea3 = getChromosome(dataa3)
chromosomeb3 = getChromosome(datab3)
chromosomed3 = getChromosome(datad3)
chromosomee3 = getChromosome(datae3)

#chromosomea2 = getChromosome(dataa2)
#chromosomea3 = getChromosome(dataa3)
#chromosomej60 = getChromosome(dataj60)
#chromosomeh2 = getChromosome(datah2)
#chromosomea4 = getChromosome(dataa4)
#chromosomea5 = getChromosome(dataa5)
#chromosomea_r26 = getChromosome(dataa_r26)
#chromosomea3_r26 = getChromosome(dataa3_r26)

def getPosition1(data):
    lines = [s.strip().split() for s in data.splitlines()]   
    position1 = []
    for j in range (0, len(lines)):
        if lines[j][0] != "chrXII" or (int(lines[j][2]) < 451430) or (int(lines[j][2]) > 468920 and int(lines[j][2]) < 489928) or (int(lines[j][2]) > 490455):
            position1.append(float(lines[j][1]))
    return position1
position1 = getPosition1(data)
position1_16 = getPosition1(data16)
position1_30 = getPosition1(data30)
position1_A = getPosition1(dataUVA)
position1_Og = getPosition1(dataOg)
position1_26 = getPosition1(data26)
position1_B = getPosition1(dataUVB)
position1_B16 = getPosition1(dataB16)
position1_B26 = getPosition1(dataB26)
position1_B30 = getPosition1(dataB30)
position1_r1 = getPosition1(datar1)
position1_r2 = getPosition1(datar2)
position1c = getPosition1(datac)
position1c1 = getPosition1(datac1)
position1w = getPosition1(dataw)
position1w1 = getPosition1(dataw1)
position1_a1 = getPosition1(dataa1)
position1_b1 = getPosition1(datab1)
position1_d1 = getPosition1(datad1)
position1_e1 = getPosition1(datae1)
position1_a2 = getPosition1(dataa2)
position1_b2 = getPosition1(datab2)
position1_d2 = getPosition1(datad2)
position1_e2 = getPosition1(datae2)
position1_a3 = getPosition1(dataa3)
position1_b3 = getPosition1(datab3)
position1_d3 = getPosition1(datad3)
position1_e3 = getPosition1(datae3)
#position1_a2 = getPosition1(dataa2)
#position1_a3 = getPosition1(dataa3)
#position1_j60 = getPosition1(dataj60)
#position1_h2 = getPosition1(datah2)
#position1_a4 = getPosition1(dataa4)
#position1_a5 = getPosition1(dataa5)
#position1_ar26 = getPosition1(dataa_r26)
#position1_a3r26 = getPosition1(dataa3_r26)


def getPosition2(data):
    lines = [s.strip().split() for s in data.splitlines()]   
    position2 = []
    for k in range (0, len(lines) ):
        if lines[k][0] != "chrXII" or (int(lines[k][2]) < 451430) or (int(lines[k][2]) > 468920 and int(lines[k][2]) < 489928) or (int(lines[k][2]) > 490455):
            position2.append(float(lines[k][2]))
    return position2
position2 = getPosition2(data)
position2_16 = getPosition2(data16)
position2_30 = getPosition2(data30)
position2_A = getPosition2(dataUVA)
position2_Og = getPosition2(dataOg)
position2_26 = getPosition2(data26)
position2_B = getPosition2(dataUVB)
position2_B16 = getPosition2(dataB16)
position2_B26 = getPosition2(dataB26)
position2_B30 = getPosition2(dataB30)
position2_r1 = getPosition2(datar1)
position2_r2 = getPosition2(datar2)
position2c = getPosition2(datac)
position2c1 = getPosition2(datac1)
position2w = getPosition2(dataw)
position2w1 = getPosition2(dataw1)
position2_a1 = getPosition2(dataa1)
position2_b1 = getPosition2(datab1)
position2_d1 = getPosition2(datad1)
position2_e1 = getPosition2(datae1)
position2_a2 = getPosition2(dataa2)
position2_b2 = getPosition2(datab2)
position2_d2 = getPosition2(datad2)
position2_e2 = getPosition2(datae2)
position2_a3 = getPosition2(dataa3)
position2_b3 = getPosition2(datab3)
position2_d3 = getPosition2(datad3)
position2_e3 = getPosition2(datae3)
#position2_a2 = getPosition2(dataa2)
#position2_a3 = getPosition2(dataa3)
#position2_j60 = getPosition2(dataj60)
#position2_h2 = getPosition2(datah2)
#position2_a4 = getPosition2(dataa4)
#position2_a5 = getPosition2(dataa5)
#position2_ar26 = getPosition2(dataa_r26)
#position2_a3r26 = getPosition2(dataa3_r26)


def getTrinucleotide(data):
    lines = [s.strip().split() for s in data.splitlines()]   
    trinucleotide = []
    for m in range (0, len(lines)):
        if lines[m][0] != "chrXII" or (int(lines[m][2]) < 451430) or (int(lines[m][2]) > 468920 and int(lines[m][2]) < 489928) or (int(lines[m][2]) > 490455):
            trinucleotide.append(lines[m][3].replace('-', ''))
    return trinucleotide
trinucleotide = getTrinucleotide(data)
trinucleotide16 = getTrinucleotide(data16)
trinucleotide30 = getTrinucleotide(data30)
trinucleotideA = getTrinucleotide(dataUVA)
trinucleotideOg = getTrinucleotide(dataOg)
trinucleotide26 = getTrinucleotide(data26)
trinucleotideB = getTrinucleotide(dataUVB)
trinucleotideB16 = getTrinucleotide(dataB16)
trinucleotideB26 = getTrinucleotide(dataB26)
trinucleotideB30 = getTrinucleotide(dataB30)
trinucleotider1 = getTrinucleotide(datar1)
trinucleotider2 = getTrinucleotide(datar2)
trinucleotidec = getTrinucleotide(datac)
trinucleotidec1 = getTrinucleotide(datac1)
trinucleotidew = getTrinucleotide(dataw)
trinucleotidew1 = getTrinucleotide(dataw1)
trinucleotidea1 = getTrinucleotide(dataa1)
trinucleotideb1 = getTrinucleotide(datab1)
trinucleotided1 = getTrinucleotide(datad1)
trinucleotidee1 = getTrinucleotide(datae1)
trinucleotidea2 = getTrinucleotide(dataa2)
trinucleotideb2 = getTrinucleotide(datab2)
trinucleotided2 = getTrinucleotide(datad2)
trinucleotidee2 = getTrinucleotide(datae2)
trinucleotidea3 = getTrinucleotide(dataa3)
trinucleotideb3 = getTrinucleotide(datab3)
trinucleotided3 = getTrinucleotide(datad3)
trinucleotidee3 = getTrinucleotide(datae3)


def getCounts(data):
    lines = [s.strip().split() for s in data.splitlines()]   
    damage_counts = []
    
    for n in range (0, len(lines)):
        #if lines[n][0] != "chrXII" or (lines[n][2] < 451430) or (lines[n][2] > 468920 and lines[n][2] < 489928) or (lines[n][2] > 490455):
            damage_counts.append(lines[n][3])
    return damage_counts
#countsa1 = getCounts(dataa1)
#countsa2 = getCounts(dataa2)
#countsa3 = getCounts(dataa3)
#countsj60 = getCounts(dataj60)
#countsh2 = getCounts(datah2)
#countsa4 = getCounts(dataa4)
#countsa5 = getCounts(dataa5)
#countsa_r26 = getCounts(dataa_r26)
#countsa3_r26 = getCounts(dataa3_r26)


with open("Raghuraman_replication_timing_saccer3.bed") as t:
    data2 = t.read()
lines2 = [u.strip().split() for u in data2.splitlines()]
chr = []
for v in range (0, len(lines2)):
    #if lines2[v][0] != "chrXII" or (int(lines2[v][2]) < 451430) or (int(lines2[v][2]) > 468920 and int(lines2[v][2]) < 489928) or (int(lines2[v][2]) > 490455):
        chr.append(lines2[v][0])
       
    
start = []
for w in range (0, len(lines2)):
    #if lines2[w][0] != "chrXII"or (int(lines2[w][2]) < 451430) or (int(lines2[w][2]) > 468920 and int(lines2[w][2]) < 489928) or (int(lines2[w][2]) > 490455):
        start.append(float(lines2[w][2]))

#scores need to be sorted in order to 
#separate into early, late, and middle replication times
#but unsorted to assign times to mutations
#since indices must match
scores = []
sorted_scores = []
for x in range (0, len(lines2)):
    #if lines2[x][0] != "chrXII"or (int(lines2[x][2]) < 451430) or (int(lines2[x][2]) > 468920 and int(lines2[x][2]) < 489928) or (int(lines2[x][2]) > 490455):
        scores.append(lines2[x][3])
        sorted_scores.append(lines2[x][3])
    
    

for sc in range(0, len(scores)):
    scores[sc] = float(scores[sc])

for ss in range(0, len(sorted_scores)):
    sorted_scores[ss] = float(sorted_scores[ss])

sorted_scores.sort()

#Determines replication time at which a mutation occurs
def calculateTime(time1, time2, pos, pos1, pos2):
    fraction = (pos -pos1)/(pos2 - pos1)
    time = time1 + (time2- time1) * fraction
    return float(time)

for s1 in range(0, len(start)):
    start[s1] = float(start[s1])




#Define 1/3 of genome as early replicating,
#1/3 as middle replicating, and 1/3 as late 
#Separate further into 3 bins each containing 1/9
#of the genome
early = {}
early1 = {}
early2 = {}
early3 = {}
middle = {}
middle1 = {}
middle2 = {}
middle3 = {}
late = {}
late1 =  {}
late2 = {}
late3 = {}

fe = open('early_regions.bed', 'w+')
fm = open('middle_regions.bed', 'w+')
fl = open('late_regions.bed', 'w+')



for y in range(0, len(sorted_scores)):
    if y <= len(sorted_scores)/3:
        early[y] = sorted_scores[y]
        if y <= len(sorted_scores)/9:
            early1[y] = sorted_scores[y]
        elif y <= (2* len(sorted_scores))/ 9:
            early2[y] = sorted_scores[y]
        elif y <= len(sorted_scores)/3:
            early3[y] = sorted_scores[y]
early_1 = 0
early1_1 = early1[int(len(sorted_scores)/9)]
early2_1 = early2[int((2* len(sorted_scores))/ 9)]
early_2 = early[len(early) - 1]
for z in range(0, len(sorted_scores)):
    if z >=  (2 * len(sorted_scores))/3 and z <= len(sorted_scores):
        late[z] = sorted_scores[z]
        if z >= (8 * len(sorted_scores))/9:
            late3[z] = sorted_scores[z]
        elif z >= (7 * len(sorted_scores))/9:
             late2[z] = sorted_scores[z]
        elif z >=  (2 * len(sorted_scores))/3:
            late1[z] = sorted_scores[z] 
       
late_1 = late[int(2 * len(sorted_scores)/3 + 1)]
late2_1 = late2[int((7 * len(sorted_scores))/9)+ 1]
late3_1 = late3[int((8 * len(sorted_scores))/9) + 1]
late_2 = late[int(len(sorted_scores)- 1)]
for q in range(0, len(sorted_scores) ):
    if q > len(sorted_scores)/3 and q < (2 * len(sorted_scores))/3:
        middle[q] = sorted_scores[q]
        if q > len(sorted_scores)/3 and q <= (4 * len(sorted_scores))/9:
            middle1[q] = sorted_scores[q]
        elif q > len(sorted_scores)/3 and q <= (5 * len(sorted_scores))/9:
            middle2[q] = sorted_scores[q]
        elif q > len(sorted_scores)/3 and q < (2 * len(sorted_scores))/3:
            middle3[q] = sorted_scores[q]
middle_1 = early_2 
middle1_1 = middle1[int((4 * len(sorted_scores))/9)]
middle2_1 = middle2[int((5 * len(sorted_scores))/9)]
middle_2 = late_1


#for s in range(0, len(scores)):
    #if scores[s] <= early_2:
        #fe.write(chr[s] + '\t' + str(int(start[s]) - 1) + '\t' + str(start[s]) + '\t' + str(scores[s]))
        #fe.write('\n')
    #if scores[s] <= middle_2:
        #fm.write(chr[s] + '\t' + str(int(start[s]) - 1) + '\t' + str(start[s]) + '\t' + str(scores[s]))
        #fm.write('\n')
    #else:
        #fl.write(chr[s] + '\t' + str(int(start[s]) - 1) + '\t' + str(start[s]) + '\t' + str(scores[s]))
        #fl.write('\n')
#fe.close()
#fm.close()
#fl.close()

#percentage of early, middle, and late replicating regions in chromosome
def ChromosomePercentages(chr_name):
    ch_total = 0
    ch_early = 0
    ch_middle = 0
    ch_late = 0
    for s6 in range(0, len(scores)):
        if chr[s6] == chr_name:
            if scores[s6] <= middle_1:
                ch_total = ch_total + 1
                ch_early = ch_early + 1
            if scores[s6] <= late_1:
                ch_total = ch_total + 1
                ch_middle = ch_middle + 1
            else:
                ch_total = ch_total + 1
                ch_late = ch_late + 1
    early_percentage = (ch_early / ch_total) * 100
    print(early_percentage)
    middle_percentage = (ch_middle / ch_total) * 100
    print(middle_percentage)
    late_percentage = (ch_late / ch_total) * 100
    print(late_percentage)
    return ([early_percentage, middle_percentage, late_percentage])
percentages = ChromosomePercentages('chrIX')
print(percentages)


f1 = open("replication_time_definition", 'w+')
f1.write('Replication Time Definitions')
f1.write('\n')
f1.write('Early=' + str(early_1) +'-' + str(early_2))
f1.write('\n')
f1.write('Bin1=' + str(early_1) +'-' + str(early1_1))
f1.write('\n')
f1.write('Bin2=' + str(early1_1) +'-' + str(early2_1))
f1.write('\n')
f1.write('Bin3=' + str(early2_1) +'-' + str(early_2))
f1.write('\n')
f1.write('Middle=' + str(middle_1) +'-' + str(middle_2))
f1.write('\n')
f1.write('Bin4=' + str(middle_1) +'-' + str(middle1_1))
f1.write('\n')
f1.write('Bin5=' + str(middle1_1) +'-' + str(middle2_1))
f1.write('\n')
f1.write('Bin6=' + str(middle2_1) +'-' + str(middle_2))
f1.write('\n')
f1.write('Late=' + str(late_1) +'-' + str(late_2))
f1.write('\n')
f1.write('Bin7=' + str(late_1) +'-' + str(late2_1))
f1.write('\n')
f1.write('Bin8=' + str(late2_1) +'-' + str(late3_1))
f1.write('\n')
f1.write('Bin9=' + str(late3_1) +'-' + str(late_2))

f1.close()



#Goes through position and start arrays in parallel
#to assign replication times for each mutation based on scores
def AssignTimes(chromosome, position2, start, scores):
    times = []
    po = 0
    st = 0

    while po < len(position2)  and st <= len(start) - 2:
        
        if position2[po] >= start[st] and position2[po] <= start[st + 1] and chromosome[po] == chr[st]:
            time1 = (scores[st])
            time2 = (scores[st + 1])
            pos = (position2[po])
            pos1 = (start[st])
            pos2 = (start[st + 1])
            time = calculateTime(time1, time2, pos, pos1, pos2)
            times.append(time)
            po = po + 1

        elif chromosome[po] < chr[st]:
            po = po + 1
            times.append('NA')
        elif chromosome[po] > chr[st]:
            st = st + 1
    
        elif position2[po] >= start[st] :
            st = st + 1
        
        elif position2[po] <= start[st]:
            po = po + 1
            times.append('NA')
    return times

times = AssignTimes(chromosome, position2, start, scores)
times16 = AssignTimes(chromosome16, position2_16, start, scores)
times30 = AssignTimes(chromosome30 , position2_30, start, scores)
times_A = AssignTimes(chromosomeA, position2_A, start, scores)
times_Og = AssignTimes(chromosomeOg, position2_Og, start, scores)
times26 = AssignTimes(chromosome26, position2_26, start, scores)
times_B = AssignTimes(chromosomeB, position2_B, start, scores)
timesB16 = AssignTimes(chromosomeB16, position2_B16, start, scores)
timesB26 = AssignTimes(chromosomeB26, position2_B26, start, scores)
timesB30 = AssignTimes(chromosomeB30, position2_B30, start, scores)
timesr1 = AssignTimes(chromosomer1, position2_r1, start, scores)
timesr2 = AssignTimes(chromosomer2, position2_r2, start, scores)
timesc = AssignTimes(chromosomec, position2c, start, scores)
timesc1 = AssignTimes(chromosomec1, position2c1, start, scores)
timesw = AssignTimes(chromosomew, position2w, start, scores)
timesw1 = AssignTimes(chromosomew1, position2w1, start, scores)
timesa1 = AssignTimes(chromosomea1, position2_a1, start, scores)
timesb1 = AssignTimes(chromosomeb1, position2_b1, start, scores)
timesd1 = AssignTimes(chromosomed1, position2_d1, start, scores)
timese1 = AssignTimes(chromosomee1, position2_e1, start, scores)
timesa2 = AssignTimes(chromosomea2, position2_a2, start, scores)
timesb2 = AssignTimes(chromosomeb2, position2_b2, start, scores)
timesd2 = AssignTimes(chromosomed2, position2_d2, start, scores)
timese2 = AssignTimes(chromosomee2, position2_e2, start, scores)
timesa3 = AssignTimes(chromosomea3, position2_a3, start, scores)
timesb3 = AssignTimes(chromosomeb3, position2_b3, start, scores)
timesd3 = AssignTimes(chromosomed3, position2_d3, start, scores)
timese3 = AssignTimes(chromosomee3, position2_e3, start, scores)
#timesa2 = AssignTimes(chromosomea2, position2_a2, start, scores)
#timesa3 = AssignTimes(chromosomea3, position2_a3, start, scores)
#timesj60 = AssignTimes(chromosomej60, position2_j60, start, scores)
#timesh2 = AssignTimes(chromosomeh2, position2_h2, start, scores)
#timesa4 = AssignTimes(chromosomea4, position2_a4, start, scores)
#timesa5 = AssignTimes(chromosomea5, position2_a5, start, scores)
#timesa_r26 = AssignTimes(chromosomea_r26, position2_ar26, start, scores)
#timesa3_r26 = AssignTimes(chromosomea3_r26, position2_a3r26, start, scores)


#Numbers and percentages of mutations in early,middle,
#and late replication times. Adds to dictionaries and separates
#into three files
def CountMutations(filename, times, chromosome, chr_name):
    early_mut = {}
    early1_mut = {}
    early2_mut = {}
    early3_mut = {}
    middle_mut = {}
    middle1_mut = {}
    middle2_mut = {}
    middle3_mut = {}
    late_mut = {}
    late1_mut = {}
    late2_mut = {}
    late3_mut = {}
    chr_early = {}
    chr_middle = {}
    chr_late = {}

    for tm in range(0, len(times) ):
        if times[tm] != 'NA' and times[tm ] >= 0 and times[tm] < early_2:
            early_mut[tm] = times[tm]
            if times[tm ] >= 0 and times[tm] < early1_1:
                early1_mut[tm] = times[tm]
            if times[tm] >= early1_1 and times[tm] < early2_1:
                early2_mut[tm] = times[tm]
            if times[tm] >= early2_1 and times[tm] <early_2:
                early3_mut[tm] = times[tm]
            if chromosome[tm] == chr_name:
                chr_early[tm] = times[tm]
        if times[tm] != 'NA' and times[tm ] >= middle_1 and times[tm] < middle_2:
            middle_mut[tm] = times[tm]
            if times[tm] >= middle_1 and times[tm] < middle1_1:
                middle1_mut[tm] = times[tm]
            if times[tm] >= middle1_1 and times[tm] < middle2_1:
                middle2_mut[tm] = times[tm]
            if times[tm] >= middle2_1 and times[tm] < middle_2:
                middle3_mut[tm] = times[tm]
            if chromosome[tm] == chr_name:
                chr_middle[tm] = times[tm]
        if times[tm] != 'NA' and times[tm ] >= late_1 :
            late_mut[tm] = times[tm]
            if times[tm] >= late_1 and times[tm] < late2_1:
                late1_mut[tm] = times[tm]
            if times[tm] >= late2_1 and times[tm] < late3_1:
                late2_mut[tm] = times[tm]
            if times[tm] >= late3_1:
                late3_mut[tm] = times[tm]
            if chromosome[tm] == chr_name:
                chr_late[tm] = times[tm]
    total = len(early_mut)  + len(middle_mut) + len(late_mut)
    total_chr = len(chr_early)  + len(chr_middle) + len(chr_late)
    if total != 0:
        percent_early = (len(early_mut)/ total) * 100
        percent_middle = (len(middle_mut)/ total) * 100
        percent_late = (len(late_mut)/ total) * 100
    else:
        percent_early = 0
        percent_middle = 0
        percent_late = 0
    if total_chr != 0:
        percent_early_chr = (len(chr_early)/ total_chr) * 100
        percent_middle_chr = (len(chr_middle)/ total_chr) * 100
        percent_late_chr = (len(chr_late)/ total_chr) * 100
    else:
        percent_early_chr = 'NA'
        percent_middle_chr = 'NA'
        percent_late_chr = 'NA'
    f1 = open(filename, 'w+')
    f1.write('Early Mutation Number: '+ str(len(early_mut)) + '  Percentage: ' + str(percent_early) + '%')
    f1.write('\n')
    f1.write('Middle Mutation Number: '+ str(len(middle_mut)) + '  Percentage: ' + str(percent_middle) + '%')
    f1.write('\n')
    f1.write('Late Mutation Number: ' + str(len(late_mut)) + '  Percentage: ' + str(percent_late) + '%') 
    f1.write('\n')
    f1.write('Early Mutation Number for '+ chr_name + ': ' + str(len(chr_early)) + '  Percentage: ' + str(percent_early_chr) + '%')
    f1.write('\n')
    f1.write('Middle Mutation Number for '+ chr_name + ': ' + str(len(chr_middle)) + '  Percentage: ' + str(percent_middle_chr) + '%')
    f1.write('\n')
    f1.write('Late Mutation Number for '+ chr_name + ': ' + str(len(chr_late)) + '  Percentage: ' + str(percent_late_chr) + '%')


    f1.close()
    return [early_mut, middle_mut, late_mut, early1_mut, early2_mut, early3_mut, middle1_mut, middle2_mut, middle3_mut, late1_mut, late2_mut, late3_mut]

muts = CountMutations("WT_mutations", times, chromosome, 'chrXVI')
early_mut = muts[0]
middle_mut = muts[1]
late_mut = muts[2]
muts_16= CountMutations("rad16_mutations", times16, chromosome16, 'chrXVI')
early_mut_16 = muts_16[0]
middle_mut_16 = muts_16[1]
late_mut_16 = muts_16[2]
muts_30 = CountMutations("rad30_mutations", times30, chromosome30, 'chrII')
early_mut_30 = muts_30[0]
middle_mut_30 = muts_30[1]
late_mut_30 = muts_30[2]
muts_A = CountMutations("UVA_mutations", times_A, chromosomeA, 'chrII')
early_mut_A = muts_A[0]
middle_mut_A = muts_A[1]
late_mut_A = muts_A[2]
muts_Og = CountMutations("ogg1_mutations", times_Og, chromosomeOg, 'chrII')
early_mut_Og = muts_Og[0]
middle_mut_Og = muts_Og[1]
late_mut_Og = muts_Og[2]
muts_26= CountMutations("rad26_mutations", times26, chromosome26, 'chrII')
early_mut_26 = muts_26[0]
middle_mut_26 = muts_26[1]
late_mut_26 = muts_26[2]
muts_B = CountMutations("UVB_mutations", times_B, chromosomeB, 'chrXVI')
early_mut_B = muts_B[0]
middle_mut_B = muts_B[1]
late_mut_B = muts_B[2]
muts_B16 = CountMutations("UVB_rad16_mutations", timesB16, chromosomeB16, 'chrXVI')
early_mut_B16 = muts_B16[0]
middle_mut_B16 = muts_B16[1]
late_mut_B16 = muts_B16[2]
muts_B26 = CountMutations("UVB_rad26_mutations", timesB26, chromosomeB26, 'chrII')
early_mut_B26 = muts_B26[0]
middle_mut_B26 = muts_B26[1]
late_mut_B26 = muts_B26[2]
muts_B30 = CountMutations("UVB_rad30_mutations", timesB30, chromosomeB30, 'chrII')
early_mut_B30 = muts_B30[0]
middle_mut_B30 = muts_B30[1]
late_mut_B30 = muts_B30[2]
muts_r1 = CountMutations("WT_r1_mutations", timesr1, chromosomer1, 'chrII')
early_mut_r1 = muts_r1[0]
middle_mut_r1 = muts_r1[1]
late_mut_r1 = muts_r1[2]
muts_r2 = CountMutations("WT_r2_mutations", timesr2, chromosomer2, 'chrII')
early_mut_r2 = muts_r2[0]
middle_mut_r2 = muts_r2[1]
late_mut_r2 = muts_r2[2]
mutsc = CountMutations("complex_rad16_mutations", timesc, chromosomec, 'chrVI')
early_mutc = mutsc[0]
middle_mutc = mutsc[1]
late_mutc = mutsc[2]
mutsc1 = CountMutations("complex_rad16_UVB_mutations", timesc1, chromosomec1, 'chrVI')
early_mutc1 = mutsc1[0]
middle_mutc1 = mutsc1[1]
late_mutc1 = mutsc1[2]
mutsw = CountMutations("complex_mutations", timesw, chromosomew, 'chrVI')
early_mutw = mutsw[0]
middle_mutw = mutsw[1]
late_mutw = mutsw[2]
mutsw1 = CountMutations("complex_UVB_mutations", timesw1, chromosomew1, 'chrVI')
early_mutw1 = mutsw1[0]
middle_mutw1 = mutsw1[1]
late_mutw1 = mutsw1[2]
muts_a1 = CountMutations("WT_indels", timesa1, chromosomea1, 'chrVI')
early_mut_a1 = muts_a1[0]
middle_mut_a1 = muts_a1[1]
late_mut_a1 = muts_a1[2]
muts_b1 = CountMutations("rad16_indels", timesb1, chromosomeb1, 'chrVI')
early_mut_b1 = muts_b1[0]
middle_mut_b1 = muts_b1[1]
late_mut_b1 = muts_b1[2]
muts_d1 = CountMutations("rad26_indels", timesd1, chromosomed1, 'chrVI')
early_mut_d1 = muts_d1[0]
middle_mut_d1 = muts_d1[1]
late_mut_d1 = muts_d1[2]
muts_e1 = CountMutations("rad30_indels", timese1, chromosomee1, 'chrVI')
early_mut_e1 = muts_e1[0]
middle_mut_e1 = muts_e1[1]
late_mut_e1 = muts_e1[2]
muts_a2 = CountMutations("WT_deletions", timesa2, chromosomea2, 'chrVI')
early_mut_a2 = muts_a2[0]
middle_mut_a2 = muts_a2[1]
late_mut_a2 = muts_a2[2]
muts_b2 = CountMutations("rad16_deletions", timesb2, chromosomeb2, 'chrVI')
early_mut_b2 = muts_b2[0]
middle_mut_b2 = muts_b2[1]
late_mut_b2 = muts_b2[2]
muts_d2 = CountMutations("rad26_deletions", timesd2, chromosomed2, 'chrVI')
early_mut_d2 = muts_d2[0]
middle_mut_d2 = muts_d2[1]
late_mut_d2 = muts_d2[2]
muts_e2 = CountMutations("rad30_deletions", timese2, chromosomee2, 'chrVI')
early_mut_e2 = muts_e2[0]
middle_mut_e2 = muts_e2[1]
late_mut_e2 = muts_e2[2]
muts_a3 = CountMutations("WT_insertions", timesa3, chromosomea3, 'chrVI')
early_mut_a3 = muts_a3[0]
middle_mut_a3 = muts_a3[1]
late_mut_a3 = muts_a3[2]
muts_b3 = CountMutations("rad16_insertions", timesb3, chromosomeb3, 'chrVI')
early_mut_b3 = muts_b3[0]
middle_mut_b3 = muts_b3[1]
late_mut_b3 = muts_b3[2]
muts_d3 = CountMutations("rad26_insertions", timesd3, chromosomed3, 'chrVI')
early_mut_d3 = muts_d3[0]
middle_mut_d3 = muts_d3[1]
late_mut_d3 = muts_d3[2]
muts_e3 = CountMutations("rad30_insertions", timese3, chromosomee3, 'chrVI')
early_mut_e3 = muts_e3[0]
middle_mut_e3 = muts_e3[1]
late_mut_e3 = muts_e3[2]
#muts_a2 = CountMutations("WT_a2_mutations", timesa2 chromosomea2, 'chrVI')
#early_mut_a2 = muts_a2[0]
#middle_mut_a2 = muts_a2[1]
#late_mut_a2 = muts_a2[2]
#muts_a3 = CountMutations("WT_a3_mutations", timesa3, chromosomea3, 'chrVI')
#early_mut_a3 = muts_a3[0]
#middle_mut_a3 = muts_a3[1]
#late_mut_a3 = muts_a3[2]
#muts_j60 = CountMutations("WT_60J_mutations", timesj60, chromosomej60, 'chrVI')
#early_mut_j60 = muts_j60[0]
#middle_mut_j60 = muts_j60[1]
#late_mut_j60 = muts_j60[2]
#muts_h2 = CountMutations("WT_h2_mutations", timesh2, chromosomeh2, 'chrVI')
#early_mut_h2 = muts_h2[0]
#middle_mut_h2 = muts_h2[1]
#late_mut_h2 = muts_h2[2]
#muts_a4 = CountMutations("rad16_a4_mutations", timesa4, chromosomea4, 'chrVI')
#early_mut_a4 = muts_a4[0]
#middle_mut_a4 = muts_a4[1]
#late_mut_a4 = muts_a4[2]
#muts_a5 = CountMutations("rad16_a5_rep2_mutations", timesa5, chromosomea5, 'chrVI')
#early_mut_a5 = muts_a5[0]
#middle_mut_a5 = muts_a5[1]
#late_mut_a5 = muts_a5[2]
#muts_ar26 = CountMutations("rad26_0hr_mutations", timesa_r26, chromosomea_r26, 'chrVI')
#early_mut_ar26 = muts_ar26[0]
#middle_mut_ar26 = muts_ar26[1]
#late_mut_ar26 = muts_ar26[2]
#muts_a3r26 = CountMutations("rad26_2hr_mutations", timesa3_r26, chromosomea3_r26, 'chrVI')
#early_mut_a3r26 = muts_a3r26[0]
#middle_mut_a3r26 = muts_a3r26[1]
#late_mut_a3r26 = muts_a3r26[2]

#Calculates average times for each of nine bins for use in linear model
def AverageTimes(muts, file):
    early1_mut = muts[3]
    early2_mut = muts[4]
    early3_mut = muts[5]
    middle1_mut = muts[6]
    middle2_mut = muts[7]
    middle3_mut = muts[8]
    late1_mut = muts[9]
    late2_mut = muts[10]
    late3_mut = muts[11]
    sum1 = 0
    for mut1 in early1_mut.keys():
        sum1 = sum1 + early1_mut[mut1]
    average1 = sum1/len(early1_mut)
    sum2 = 0
    for mut2 in early2_mut.keys():
        sum2 = sum2 + early2_mut[mut2]
    average2 = sum2/len(early2_mut)
    sum3 = 0
    for mut3 in early3_mut.keys():
        sum3 = sum3 + early3_mut[mut3]
    average3 = sum3/len(early3_mut)
    sum4 = 0
    for mut4 in middle1_mut.keys():
        sum4 = sum4 + middle1_mut[mut4]
    average4 = sum4/len(middle1_mut)
    sum5 = 0
    for mut5 in middle2_mut.keys():
        sum5 = sum5 + middle2_mut[mut5]
    average5 = sum5/len(middle2_mut)
    sum6 = 0
    for mut6 in middle3_mut.keys():
        sum6 = sum6 + middle3_mut[mut6]
    average6 = sum6/len(middle3_mut)
    sum7 = 0
    for mut7 in late1_mut.keys():
        sum7 = sum7 + late1_mut[mut7]
    average7 = sum7/len(late1_mut)
    sum8 = 0
    for mut8 in late2_mut.keys():
        sum8 = sum8 + late2_mut[mut8]
    average8 = sum8/len(late2_mut)
    sum9 = 0
    for mut9 in late3_mut.keys():
        sum9 = sum9 + late3_mut[mut9]
    average9 = sum9/len(late3_mut)
    f = open(file, 'w+')
    f.write('Bin1: ' + str(average1))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(early1_mut)))
    f.write('\n')
    f.write('Bin2: ' + str(average2))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(early2_mut)))
    f.write('\n')
    f.write('Bin3: ' + str(average3))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(early3_mut)))
    f.write('\n')
    f.write('Bin4: ' + str(average4))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(middle1_mut)))
    f.write('\n')
    f.write('Bin5: ' + str(average5))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(middle2_mut)))
    f.write('\n')
    f.write('Bin6: ' + str(average6))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(middle3_mut)))
    f.write('\n')
    f.write('Bin7: ' + str(average7))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(late1_mut)))
    f.write('\n')
    f.write('Bin8: ' + str(average8))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(late2_mut)))
    f.write('\n')
    f.write('Bin9: ' + str(average9))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(late3_mut)))
    f.close()
AverageTimes(muts, 'average_times.txt')
AverageTimes(muts_16, 'average_times_rad16.txt')
AverageTimes(muts_26, 'average_times_rad26.txt')
AverageTimes(muts_30, 'average_times_rad30.txt')
AverageTimes(muts_B, 'average_times_UVB.txt')
AverageTimes(muts_B16, 'average_times_rad16_UVB.txt')
AverageTimes(muts_B26, 'average_times_rad26_UVB.txt')
AverageTimes(muts_B30, 'average_times_rad30_UVB.txt')
AverageTimes(muts_A, 'average_times_UVA.txt')
AverageTimes(muts_Og, 'average_times_ogg1.txt')
AverageTimes(mutsc, 'average_times_rad16_complex')
AverageTimes(mutsc1, 'average_times_rad16_complex_UVB')
AverageTimes(mutsw, 'average_times_complex')
AverageTimes(mutsw1, 'average_times_complex_UVB')
#AverageTimes(muts_a1, 'average_times_indels')
#AverageTimes(muts_b1, 'average_times_indels_rad16')
#AverageTimes(muts_d1, 'average_times_indels_rad26')
#AverageTimes(muts_e1, 'average_times_indels_rad30')
#AverageTimes(muts_a2, 'average_times_deletions')
#AverageTimes(muts_b2, 'average_times_deletions_rad16')
#AverageTimes(muts_d2, 'average_times_deletions_rad26')
#AverageTimes(muts_e2, 'average_times_deletions_rad30')
#AverageTimes(muts_a3, 'average_times_insertions')
#AverageTimes(muts_b3, 'average_times_insertions_rad16')
#AverageTimes(muts_d3, 'average_times_insertions_rad26')
#AverageTimes(muts_e3, 'average_times_insertions_rad30')

#Uses times array to add a column for replication times to bed files
def AddTimeColumn(filename, file1, file2, file3, data, times, early_mut, middle_mut, late_mut):
    fi = open(filename, 'w+')
    fa = open(file1 , 'w+')
    fb = open(file2, 'w+')
    fc = open(file3, 'w+')   
    e = 0   
    for d in data.splitlines():
        d.strip()
        if e < len(times):
     
            d = d + '\t' + str(times[e])
        
        
            if e in early_mut.keys():
                fa.write(d)
                fa.write('\n')
            if e in middle_mut.keys():
                fb.write(d)
                fb.write('\n')
            if e in late_mut.keys():
                fc.write(d)
                fc.write('\n')
            e = e + 1
       
        
        fi.write(d)
        fi.write('\n')
    
   
    fi.close()
    fa.close()
    fb.close()
    fc.close()
AddTimeColumn("WT_UV_SNPs_sorted.bed", "early_mutations.bed", "middle_mutations.bed", "late_mutations.bed",data, times, early_mut, middle_mut, late_mut )
AddTimeColumn("rad16_UV_bothruns_allSNVs_sorted.bed", "early_mutations_rad16.bed", "middle_mutations_rad16.bed", "late_mutations_rad16.bed",data16, times16, early_mut_16, middle_mut_16, late_mut_16 )
AddTimeColumn("rad30_UV_bothruns_allSNVs_sorted.bed", "early_mutations._rad30bed", "middle_mutations_rad30.bed", "late_mutations_rad30.bed",data30, times30, early_mut_30, middle_mut_30, late_mut_30 )
AddTimeColumn("WT_UVA_allSNVs_sorted.bed", "early_mutations_UVA.bed", "middle_mutations_UVA.bed", "late_mutations_UVA.bed",dataUVA, times_A, early_mut_A, middle_mut_A, late_mut_A)
AddTimeColumn("ogg1_UVA_allSNVs_sorted.bed", "early_mutations_ogg1.bed", "middle_mutations_ogg1.bed", "late_mutations_ogg1.bed",dataOg, times_Og, early_mut_Og, middle_mut_Og, late_mut_Og)
AddTimeColumn("rad26_UV_bothruns_allSNVs_sorted.bed", "early_mutations_rad26.bed", "middle_mutations_rad26.bed", "late_mutations_rad26.bed",data26, times26, early_mut_26, middle_mut_26, late_mut_26 )
AddTimeColumn("WT_UVB_allSNVs_sorted.bed", "early_mutations_UVB.bed", "middle_mutations_UVB.bed", "late_mutations_UVB.bed",dataUVB, times_B, early_mut_B, middle_mut_B, late_mut_B)
AddTimeColumn("Rad16_UVB_allSNVs_sorted.bed", "early_mutations_rad16UVB.bed", "middle_mutations_rad16UVB.bed", "late_mutations_rad16UVB.bed",dataB16, timesB16, early_mut_B16, middle_mut_B16, late_mut_B16)
AddTimeColumn("Rad26_UVB_allSNVs_sorted.bed", "early_mutations_rad26UVB.bed", "middle_mutations_rad26UVB.bed", "late_mutations_rad26UVB.bed",dataB26, timesB26, early_mut_B26, middle_mut_B26, late_mut_B26)
AddTimeColumn("Rad30_UVB_allSNVs_sorted.bed", "early_mutations_rad30UVB.bed", "middle_mutations_rad30UVB.bed", "late_mutations_rad30UVB.bed",dataB30, timesB30, early_mut_B30, middle_mut_B30, late_mut_B30)
AddTimeColumn("WT_UV_run1_allSNVs_sorted.bed", "early_mutations_r1.bed", "middle_mutations_r1.bed", "late_mutations_r1.bed",datar1, timesr1, early_mut_r1, middle_mut_r1, late_mut_r1 )
AddTimeColumn("WT_UV_run2_allSNVs_sorted.bed", "early_mutations_r2.bed", "middle_mutations_r2.bed", "late_mutations_r2.bed",datar2, timesr2, early_mut_r2, middle_mut_r2, late_mut_r2 )
AddTimeColumn("clusters_rad16_allSNVs_sorted.bed", "clusters_rad16_early_mutations.bed", "clusters_rad16_middle_mutations.bed", "clusters_rad16_late_mutations.bed",datac, timesc, early_mutc, middle_mutc, late_mutc)
AddTimeColumn("clusters_rad16_UVB_allSNVs_sorted.bed", "clusters_rad16_early_mutations_UVB.bed", "clusters_rad16_middle_mutations_UVB.bed", "clusters_rad16_late_mutations_UVB.bed",datac1, timesc1, early_mutc1, middle_mutc1, late_mutc1)
AddTimeColumn("clusters_allSNVs_sorted.bed", "clusters_early_mutations.bed", "clusters_middle_mutations.bed", "clusters_late_mutations.bed",dataw, timesw, early_mutw, middle_mutw, late_mutw)
AddTimeColumn("clusters_UVB_allSNVs_sorted.bed", "clusters_early_mutations_UVB.bed", "clusters_middle_mutations_UVB.bed", "clusters_late_mutations_UVB.bed",dataw1, timesw1, early_mutw1, middle_mutw1, late_mutw1)
AddTimeColumn("WT_indels.bed", "early_indels.bed", "middle_indels.bed", "late_indels.bed",dataa1, timesa1, early_mut_a1, middle_mut_a1, late_mut_a1 )
AddTimeColumn("rad16_indels.bed", "early_indels_rad16.bed", "middle_indels_rad16.bed", "late_indels_rad16.bed",datab1, timesb1, early_mut_b1, middle_mut_b1, late_mut_b1 )
AddTimeColumn("rad26_indels.bed", "early_indels_rad26.bed", "middle_indels_rad26.bed", "late_indels_rad26.bed",datad1, timesd1, early_mut_d1, middle_mut_d1, late_mut_d1 )
AddTimeColumn("rad30_indels.bed", "early_indels_rad30.bed", "middle_indels_rad30.bed", "late_indels_rad30.bed",datae1, timese1, early_mut_e1, middle_mut_e1, late_mut_e1 )
AddTimeColumn("WT_deletions.bed", "early_deletions.bed", "middle_deletions.bed", "late_deletions.bed",dataa2, timesa2, early_mut_a2, middle_mut_a2, late_mut_a2 )
AddTimeColumn("rad16_deletions.bed", "early_deletions_rad16.bed", "middle_deletions_rad16.bed", "late_deletions_rad16.bed",datab2, timesb2, early_mut_b2, middle_mut_b2, late_mut_b2)
AddTimeColumn("rad26_deletions.bed", "early_deletions_rad26.bed", "middle_deletions_rad26.bed", "late_deletions_rad26.bed",datad2, timesd2, early_mut_d2, middle_mut_d2, late_mut_d2)
AddTimeColumn("rad30_deletions.bed", "early_deletions_rad30.bed", "middle_deletions_rad30.bed", "late_deletions_rad30.bed",datae2, timese2, early_mut_e2, middle_mut_e2, late_mut_e2)
AddTimeColumn("WT_insertions.bed", "early_insertions.bed", "middle_insertions.bed", "late_insertions.bed",dataa3, timesa3, early_mut_a3, middle_mut_a3, late_mut_a3 )
AddTimeColumn("rad16_insertions.bed", "early_insertions_rad16.bed", "middle_insertions_rad16.bed", "late_insertions_rad16.bed",datab3, timesb3, early_mut_b3, middle_mut_b3, late_mut_b3)
AddTimeColumn("rad26_insertions.bed", "early_insertions_rad26.bed", "middle_insertions_rad26.bed", "late_insertions_rad26.bed",datad3, timesd3, early_mut_d3, middle_mut_d3, late_mut_d3)
AddTimeColumn("rad30_insertions.bed", "early_insertions_rad30.bed", "middle_insertions_rad30.bed", "late_insertions_rad30.bed",datae3, timese3, early_mut_e3, middle_mut_e3, late_mut_e3)
#AddTimeColumn("UV2_A1_dipy_bkgd.bed", "early_mutations_UV2_A1_dipy_bkgd.bed", "middle_mutations_UV2_A1_dipy_bkgd.bed", "late_mutations_UV2_A1_dipy_bkgd.bed",dataa1, timesa1, early_mut_a1, middle_mut_a1, late_mut_a1 )
#AddTimeColumn("UV_0hr_A2_dipy_bkgd.bed", "early_mutations_UV_0hr_A2_dipy_bkgd.bed", "middle_mutations_UV_0hr_A2_dipy_bkgd.bed", "late_mutations_UV_0hr_A2_dipy_bkgd.bed",dataa2, timesa2, early_mut_a2, middle_mut_a2, late_mut_a2 )
#AddTimeColumn("UV90J_A3_dipy_bkgd.bed", "early_mutations_UV90J_A3_dipy_bkgd.bed", "middle_mutations_UV90J_A3_dipy_bkgd.bed", "late_mutations_UV90J_A3_dipy_bkgd.bed",dataa3, timesa3, early_mut_a3, middle_mut_a3, late_mut_a3 )
#AddTimeColumn("UV60J_A2_dipy_bkgd.bed", "early_mutations_UV60J_A2_dipy_bkgd.bed", "middle_mutations_UV60J_A2_dipy_bkgd.bed", "late_mutations_UV60J_A2_dipy_bkgd.bed",dataj60, timesj60, early_mut_j60, middle_mut_j60, late_mut_j60 )
#AddTimeColumn("2hr_UV2_WT_A3_1bp_dipy_bk.bed", "early_mutations_2hr_UV2_WT_A3_1bp_dipy_bk.bed", "middle_mutations_2hr_UV2_WT_A3_1bp_dipy_bk.bed", "late_mutations_2hr_UV2_WT_A3_1bp_dipy_bk.bed",datah2, timesh2, early_mut_h2, middle_mut_h2, late_mut_h2 )
#AddTimeColumn("rad16_0hr_A4_rep2_dipy_bkgd.bed", "early_mutations_rad16_0hr_A4_rep2_dipy_bkgd.bed", "middle_mutations_rad16_0hr_A4_rep2_dipy_bkgd.bed", "late_mutations_rad16_0hr_A4_rep2_dipy_bkgd.bed",dataa4, timesa4, early_mut_a4, middle_mut_a4, late_mut_a4)
#AddTimeColumn("rad16_2hr_rep2_A5_1bp_dipy_bk.bed", "early_mutations_rad16_2hr_A5_rep2_1bp_dipy_bk.bed", "middle_mutations_rad16_2hr_A5_rep2_1bp_dipy_bk.bed", "late_mutations_rad16_2hr_A5_rep2_1bp_dipy_bk.bed",dataa5, timesa5, early_mut_a5, middle_mut_a5, late_mut_a5 )
#AddTimeColumn("rad26_0hr_rep2_dipy_bkgd.bed", "early_mutations_rad26_0hr_rep2_dipy_bkgd.bed", "middle_mutations_rad26_0hr_rep2_dipy_bkgd.bed", "late_mutations_rad26_0hr_rep2_dipy_bkgd.bed",dataa_r26, timesa_r26, early_mut_ar26, middle_mut_ar26, late_mut_ar26)
#AddTimeColumn("rad26_2hr_rep2_dipy_bk.bed", "early_mutations_rad26_2hr_rep2_dipy_bk.bed", "middle_mutations_rad26_2hr_rep2_dipy_bk.bed", "late_mutations_rad26_2hr_rep2_dipy_bk.bed",dataa3_r26, timesa3_r26, early_mut_a3r26, middle_mut_a3r26, late_mut_a3r26 )

#Calculates number of mutations for each codon
def MutationRate(file1, trinucleotides, muts, file2):
    mutationRates = {}
    base_percentages = {}
    mutationRates_early = {}
    mutationRates_middle = {}
    mutationRates_late = {}
    base_percentages_early = {}
    base_percentages_middle = {}
    base_percentages_late = {}
    fm = open(file1)
    datam = fm.read()
    fm2 = open(file2, 'w+')
    for line in datam.splitlines():
        line = datam.strip().split()
    for d in range(0, len(line)):
        mutationRates[line[d]] = 0
        base_percentages[line[d]] = 0
        mutationRates_early[line[d]] = 0
        mutationRates_middle[line[d]] = 0
        mutationRates_late[line[d]] = 0
        base_percentages_early[line[d]] = 0
        base_percentages_middle[line[d]] = 0
        base_percentages_late[line[d]] = 0
    for key in mutationRates.keys():
        for t in range (0, len(trinucleotides)):
            if trinucleotides[t][1] == 'C' or trinucleotides[t][1] == 'T':
                if trinucleotides[t] == key and t in muts[0].keys():
                    mutationRates_early[key] = mutationRates_early[key] + 1
                    base_percentages_early[key] = (mutationRates_early[key] / len(muts[0].keys()))
                    mutationRates[key] = mutationRates[key] + 1
                    base_percentages[key] = (mutationRates[key] / len(trinucleotides))
                if trinucleotides[t] == key and t in muts[1].keys():
                    mutationRates_middle[key] = mutationRates_middle[key] + 1
                    base_percentages_middle[key] = (mutationRates_middle[key] / len(muts[1].keys()))
                    mutationRates[key] = mutationRates[key] + 1
                    base_percentages[key] = (mutationRates[key] / len(trinucleotides))
                if trinucleotides[t] == key and t in muts[2].keys():
                    mutationRates_late[key] = mutationRates_late[key] + 1
                    base_percentages_late[key] = (mutationRates_late[key] / len(muts[2].keys()))
                    mutationRates[key] = mutationRates[key] + 1
                    base_percentages[key] = (mutationRates[key] / len(trinucleotides))
            if trinucleotides[t][1] == 'A' or trinucleotides[t][1] == 'G':
                if GetReverseComplement(trinucleotides[t]) == key:
                    mutationRates[key] = mutationRates[key] + 1
                    base_percentages[key] = (mutationRates[key] / len(trinucleotides))
                    if t in muts[0].keys():
                        mutationRates_early[key] = mutationRates_early[key] + 1
                        base_percentages_early[key] = (mutationRates_early[key] / len(muts[0].keys()))
                        mutationRates[key] = mutationRates[key] + 1
                        base_percentages[key] = (mutationRates[key] / len(trinucleotides))
                    if t in muts[1].keys():
                        mutationRates_middle[key] = mutationRates_middle[key] + 1
                        base_percentages_middle[key] = (mutationRates_middle[key] / len(muts[1].keys()))
                        mutationRates[key] = mutationRates[key] + 1
                        base_percentages[key] = (mutationRates[key] / len(trinucleotides))
                    if t in muts[2].keys():
                        mutationRates_late[key] = mutationRates_late[key] + 1
                        base_percentages_late[key] = (mutationRates_late[key] / len(muts[2].keys()))
                        mutationRates[key] = mutationRates[key] + 1
                        base_percentages[key] = (mutationRates[key] / len(trinucleotides))
    for key in mutationRates.keys():
        fm2.write(key + ": Number = " + str(mutationRates[key])+ ' Percentage =' + str(base_percentages[key] * 100) + "%")
        fm2.write('\n')
    fm2.write('Early Replication: ')
    fm2.write('\n')
    for key in mutationRates_early.keys():
        fm2.write(key + ": Number = " + str(mutationRates_early[key])+ ' Percentage =' + str(base_percentages_early[key] * 100) + "%")
        fm2.write('\n')
    fm2.write('Middle Replication: ')
    fm2.write('\n')
    for key in mutationRates_middle.keys():
        fm2.write(key + ": Number = " + str(mutationRates_middle[key])+ ' Percentage =' + str(base_percentages_middle[key] * 100) + "%")
        fm2.write('\n')
    fm2.write('Late Replication: ')
    fm2.write('\n')
    for key in mutationRates_late.keys():
        fm2.write(key + ": Number = " + str(mutationRates_late[key])+ ' Percentage =' + str(base_percentages_late[key] * 100) + "%")
        fm2.write('\n')
    return(mutationRates)

mutationRates = MutationRate('mutationrates.txt', trinucleotide, muts, 'mutationrate_results')
mutationRates16 = MutationRate('mutationrates.txt', trinucleotide16, muts_16, 'mutationrate16_results')
mutationRatesr1 = MutationRate('mutationrates.txt', trinucleotider1, muts_r1, 'mutationrate_r1_results')
mutationRatesr2 = MutationRate('mutationrates.txt', trinucleotider2, muts_r2, 'mutationrate_r2_results')
mutationRates30 = MutationRate('mutationrates.txt', trinucleotide30, muts_30, 'mutationrate30_results')
mutationRatesUVA = MutationRate('mutationrates.txt', trinucleotideA,muts_A, 'mutationrateUVA_results')  
mutationRatesOg = MutationRate('mutationrates.txt', trinucleotideOg, muts_Og, 'mutationrateOgg1_results')
mutationRates26 = MutationRate('mutationrates.txt', trinucleotide26,muts_26, 'mutationrate26_results') 
mutationRatesUVB = MutationRate('mutationrates.txt', trinucleotideB, muts_B, 'mutationrateUVB_results')
mutationRatesB16 = MutationRate('mutationrates.txt', trinucleotideB16, muts_B16, 'mutationrateUVB16_results')
mutationRatesB26 = MutationRate('mutationrates.txt', trinucleotideB26, muts_B26, 'mutationrateUVB26_results')
mutationRatesB30 = MutationRate('mutationrates.txt', trinucleotideB30, muts_B30, 'mutationrateUVB30_results')
mutationRatesc = MutationRate('mutationrates.txt', trinucleotidec, mutsc, 'mutationrate16_complex_sub_results')
mutationRatesa1 = MutationRate('mutationrates.txt', trinucleotidea1, muts_a1, 'mutationrate_indel_results')
mutationRatesb1 = MutationRate('mutationrates.txt', trinucleotideb1, muts_b1, 'mutationrate16_indel_results')
mutationRatesd1 = MutationRate('mutationrates.txt', trinucleotided1, muts_d1, 'mutationrate26_indel_results')
mutationRatese1 = MutationRate('mutationrates.txt', trinucleotidee1, muts_e1, 'mutationrate30_indel_results')
mutationRatesa2 = MutationRate('mutationrates.txt', trinucleotidea2, muts_a2, 'mutationrate_deletion_results')
mutationRatesb2 = MutationRate('mutationrates.txt', trinucleotideb2, muts_b2, 'mutationrate16_deletion_results')
mutationRatesd2 = MutationRate('mutationrates.txt', trinucleotided2, muts_d2, 'mutationrate26_deletion_results')
mutationRatese2 = MutationRate('mutationrates.txt', trinucleotidee2, muts_e2, 'mutationrate30_deletion_results')
mutationRatesa3 = MutationRate('dinucleotide_rates.txt', trinucleotidea3, muts_a3, 'mutationrate_insertion_results')
mutationRatesb3 = MutationRate('dinucleotide_rates.txt', trinucleotideb3, muts_b3, 'mutationrate16_insertion_results')
mutationRatesd3 = MutationRate('dinucleotide_rates.txt', trinucleotided3, muts_d3, 'mutationrate26_insertion_results')
mutationRatese3 = MutationRate('dinucleotide_rates.txt', trinucleotidee3, muts_e3, 'mutationrate30_insertion_results')


rep_times = []
rep_times.append(early_1)
rep_times.append(early_2)
rep_times.append(middle_1)
rep_times.append(middle_2)
rep_times.append(late_1)
rep_times.append(late_2)


def ChromosomeAverageTimes(chromosome, file):
    value1 = float(0)
    value2 = float(0)
    value3 = float(0)
    value4 = float(0)
    value5 = float(0)
    value6 = float(0)
    value7 = float(0)
    value8 = float(0)
    value9 = float(0)
    value10 = float(0)
    value11 = float(0)
    value12 = float(0)
    value13 = float(0)
    value14 = float(0)
    value15 = float(0)
    value16 = float(0)
    valuea = 0
    valueb = 0
    valuec = 0
    valued = 0
    valuee = 0
    valuef = 0
    valueg = 0
    valueh = 0
    valuei = 0
    valuej = 0
    valuek = 0
    valuel = 0

    valuem = 0
    valuen = 0
    valueo = 0
    valuep = 0
    for t in range(0, len(scores)):
        if chromosome[t] == 'chrI':
            value1 = value1 + float(scores[t])
            valuea = valuea + 1
        if chromosome[t] == 'chrII':
            value2 = value2 + float(scores[t])
            valueb = valueb + 1
        if chromosome[t] == 'chrIII':
            value3 = value3 + float(scores[t])
            valuec = valuec + 1
        if chromosome[t] == 'chrIV':
            value4 = value4 + float(scores[t])
            valued = valued + 1
        if chromosome[t] == 'chrV':
            value5 = value5 + float(scores[t])
            valuee = valuee + 1
        if chromosome[t] == 'chrVI' :
            value6 = value6 + float(scores[t])
            valuef = valuef + 1
        if chromosome[t] == 'chrVII':
            value7 = value7 + float(scores[t])
            valueg = valueg + 1
        if chromosome[t] == 'chrVIII':
            value8 = value8 + float(scores[t])
            valueh = valueh + 1
        if chromosome[t] == 'chrIX':
            value9 = value9 + float(scores[t])
            valuei = valuei + 1
        if chromosome[t] == 'chrX':
            value10 = value10 + float(scores[t])
            valuej = valuej + 1
        if chromosome[t] == 'chrXI':
            value11 = value11 + float(scores[t])
            valuek = valuek + 1
        if chromosome[t] == 'chrXII':
            value12 = value12 + float(scores[t])
            valuel = valuel + 1
        if chromosome[t] == 'chrXIII':
            value13 = value13 + float(scores[t])
            valuem = valuem + 1
        if chromosome[t] == 'chrXIV':
            value14 = value14 + float(scores[t])
            valuen = valuen + 1
        if chromosome[t] == 'chrXV':
            value15 = value15 + float(scores[t])
            valueo = valueo + 1
        if chromosome[t] == 'chrXVI':
            value16 = value16 + float(scores[t])
            valuep = valuep + 1
 
    ft = open(file, 'w+')
    ft.write('Chromosome 1: ' + str(value1/valuea))
    ft.write('\n')
    ft.write('Chromosome 2: ' + str(value2/valueb))
    ft.write('\n')
    ft.write('Chromosome 3: ' + str(value3/valuec))
    ft.write('\n')
    ft.write('Chromosome 4: ' + str(value4/valued))
    ft.write('\n')
    ft.write('Chromosome 5: ' + str(value5/valuee))
    ft.write('\n')
    ft.write('Chromosome 6: ' + str(value6/valuef))
    ft.write('\n')
    ft.write('Chromosome 7: ' + str(value7/valueg))
    ft.write('\n')
    ft.write('Chromosome 8: ' + str(value8/valueh))
    ft.write('\n')
    ft.write('Chromosome 9: ' + str(value9/valuei))
    ft.write('\n')
    ft.write('Chromosome 10: ' + str(value10/valuej))
    ft.write('\n')
    ft.write('Chromosome 11: ' + str(value11/valuek))
    ft.write('\n')
    ft.write('Chromosome 12: ' + str(value12/valuel))
    ft.write('\n')
    ft.write('Chromosome 13: ' + str(value13/valuem))
    ft.write('\n')
    ft.write('Chromosome 14: ' + str(value14/valuen))
    ft.write('\n')
    ft.write('Chromosome 15: ' + str(value15/valueo))
    ft.write('\n')
    ft.write('Chromosome 16: ' + str(value16/valuep))
    ft.write('\n')

ChromosomeAverageTimes(chr, 'chromosome_average_times.txt')
        

early_TrinucCounts = {}
middle_TrinucCounts = {}
late_TrinucCounts = {}
early_DinucCounts = {}
middle_DinucCounts = {}
late_DinucCounts = {}

def GetDipyrimidineCounts(di_counts, seq, value):
    if seq == 'CC' or seq == 'GG':
        di_counts[value]['CC'] = di_counts[value]['CC'] + 1
    if seq == 'CT' or seq == 'AG':
        di_counts[value]['CT'] = di_counts[value]['CT'] + 1
    if seq == 'TC' or seq == 'GA':
        di_counts[value]['TC'] = di_counts[value]['TC'] + 1
    if seq == 'TT' or seq == 'AA':
        di_counts[value]['TT'] = di_counts[value]['TT'] + 1
    return di_counts

def PrintTotalDipyrimidineCounts(di_counts, file):
    file.write('Early Trinucleotide Counts: ')
    file.write('\n')
    early_sum = 0
    for key in di_counts[0]:
        early_sum = early_sum + di_counts[0][key]
        file.write(key + ": " + str(di_counts[0][key]))
        file.write('\n')
    file.write('Total: ' + str(early_sum))
    file.write('\n')
    file.write('Middle Trinucleotide Counts: ')
    file.write('\n')
    middle_sum = 0
    for key in di_counts[1]:
        middle_sum = middle_sum + di_counts[1][key]
        file.write(key + ": " + str(di_counts[1][key]))
        file.write('\n')
    file.write('Total: ' + str(middle_sum))
    file.write('\n')
    file.write('Late Trinucleotide Counts: ')
    file.write('\n')
    late_sum = 0
    for key in di_counts[2]:
        late_sum = late_sum + di_counts[2][key]
        file.write(key + ": " + str(di_counts[2][key]))
        file.write('\n')
    file.write('Total: ' + str(late_sum))
    file.write('\n')
    file.close()





dicts = [early_TrinucCounts, middle_TrinucCounts, late_TrinucCounts]
dinuc_dicts = [early_DinucCounts, middle_DinucCounts, late_DinucCounts]
di_dict = {0:{'CC':0, 'CT':0, 'TC':0, 'TT':0}, 1:{'CC':0, 'CT':0, 'TC':0, 'TT':0}, 2:{'CC':0, 'CT':0, 'TC':0, 'TT':0}}
#Counts numbers of trinucleotides in early, middle, and late replicating regions (only regions included in the timing map)
#based on input from yeast genome assembly R64
def TrinucleotideContext(file, name, start, scores, early, middle, dicts, num, counts):
    early_TrinucCounts = dicts[0]
    middle_TrinucCounts = dicts[1]
    late_TrinucCounts = dicts[2]
    
    
    
    chr_scores = []
    ft = open(file)
    sequence = ft.read()
    sequence = str(sequence.strip().upper())
    sequence = ''.join(sequence.split('\n'))
    di_counts = counts
       
            
    for i in range(1, len(scores)):
        if chr[i] == name and chr[i - 1] == name:
            
            for q in range(int(start[i -1]) - 1, int(start[i]) - 1):
                if chr[i] != "chrXII" or (q < 451430) or (q  > 468920 and q  < 489928) or (q  > 490455):
                    if sequence[(q - 1):(q+num - 1)][1] == 'A' or sequence[(q-1):(q+num -1)][1] == 'G':
                        seq1 =GetReverseComplement (sequence[(q-1):(q + num - 1)])
                    else:
                        seq1 = sequence[(q-1):(q + num - 1)]

                    time = calculateTime(scores[i - 1], scores[i],q, start[i - 1] - 1, start[i] - 1)
                    if time < early:
                        if seq1 not in early_TrinucCounts.keys():
                            early_TrinucCounts[seq1] = 1
                        else:
                            early_TrinucCounts[seq1] = early_TrinucCounts[seq1] + 1
                        di_counts = GetDipyrimidineCounts(di_counts, sequence[q-1:q+1], 0)
                    if time >= early and time < middle:
                        if seq1 not in middle_TrinucCounts.keys():
                            middle_TrinucCounts[seq1] = 1
                        else:
                            middle_TrinucCounts[seq1] = middle_TrinucCounts[seq1] + 1
                        di_counts = GetDipyrimidineCounts(di_counts, sequence[q-1:q+1], 1)

                    if time >= middle:
                        if seq1 not in late_TrinucCounts.keys():
                            late_TrinucCounts[seq1] = 1
                        else:
                            late_TrinucCounts[seq1] = late_TrinucCounts[seq1] + 1
                        di_counts = GetDipyrimidineCounts(di_counts, sequence[q-1:q+1], 2)
                    


                                         
    ft.close()
    return [early_TrinucCounts, middle_TrinucCounts, late_TrinucCounts], di_counts
di_dict1 = TrinucleotideContext('chr1.txt', 'chrI', start, scores, early_2, middle_2, dinuc_dicts, 2, di_dict)
di_dict2 = TrinucleotideContext('chr2.txt', 'chrII', start, scores, early_2, middle_2, di_dict1[0], 2, di_dict1[1])
di_dict3 = TrinucleotideContext('chr3.txt', 'chrIII', start, scores, early_2, middle_2, di_dict2[0], 2, di_dict2[1])
di_dict4 = TrinucleotideContext('chr4.txt', 'chrIV', start, scores, early_2, middle_2, di_dict3[0], 2, di_dict3[1])
di_dict5 = TrinucleotideContext('chr5.txt', 'chrV', start, scores, early_2, middle_2, di_dict4[0], 2, di_dict4[1])
di_dict6 = TrinucleotideContext('chr6.txt', 'chrVI', start, scores, early_2, middle_2, di_dict5[0], 2, di_dict5[1])
di_dict7 = TrinucleotideContext('chr7.txt', 'chrVII', start, scores, early_2, middle_2, di_dict6[0], 2, di_dict6[1])
di_dict8 = TrinucleotideContext('chr8.txt', 'chrVIII', start, scores, early_2, middle_2, di_dict7[0], 2, di_dict7[1])
di_dict9 = TrinucleotideContext('chr9.txt', 'chrIX', start, scores, early_2, middle_2, di_dict8[0], 2, di_dict8[1])
di_dict10 = TrinucleotideContext('chr10.txt', 'chrX', start, scores, early_2, middle_2, di_dict9[0], 2, di_dict9[1])
di_dict11 = TrinucleotideContext('chr11.txt', 'chrXI', start, scores, early_2, middle_2, di_dict10[0], 2, di_dict10[1])
di_dict12 = TrinucleotideContext('chr12.txt', 'chrXII', start, scores, early_2, middle_2, di_dict11[0], 2, di_dict11[1])
di_dict13 = TrinucleotideContext('chr13.txt', 'chrXIII', start, scores, early_2, middle_2, di_dict12[0], 2, di_dict12[1])
di_dict14 = TrinucleotideContext('chr14.txt', 'chrXIV', start, scores, early_2, middle_2, di_dict13[0], 2, di_dict13[1])
di_dict15 = TrinucleotideContext('chr15.txt', 'chrXV', start, scores, early_2, middle_2, di_dict14[0], 2, di_dict14[1])
dinuc_counts = TrinucleotideContext('chr16.txt', 'chrXVI', start, scores, early_2, middle_2, di_dict15[0], 2, di_dict15[1])


dict1 = TrinucleotideContext('chr1.txt', 'chrI', start, scores, early_2, middle_2, dicts, 3, di_dict)
dict2 = TrinucleotideContext('chr2.txt', 'chrII', start, scores, early_2, middle_2, dict1[0], 3, dict1[1])
dict3 = TrinucleotideContext('chr3.txt', 'chrIII', start, scores, early_2, middle_2, dict2[0], 3, dict2[1])
dict4 = TrinucleotideContext('chr4.txt', 'chrIV', start, scores, early_2, middle_2, dict3[0], 3, dict3[1])
dict5 = TrinucleotideContext('chr5.txt', 'chrV', start, scores, early_2, middle_2, dict4[0], 3, dict4[1])
dict6 = TrinucleotideContext('chr6.txt', 'chrVI', start, scores, early_2, middle_2, dict5[0], 3, dict5[1])
dict7 = TrinucleotideContext('chr7.txt', 'chrVII', start, scores, early_2, middle_2, dict6[0], 3, dict6[1])
dict8 = TrinucleotideContext('chr8.txt', 'chrVIII', start, scores, early_2, middle_2, dict7[0], 3, dict7[1])
dict9 = TrinucleotideContext('chr9.txt', 'chrIX', start, scores, early_2, middle_2, dict8[0], 3, dict8[1])
dict10 = TrinucleotideContext('chr10.txt', 'chrX', start, scores, early_2, middle_2, dict9[0], 3, dict9[1])
dict11 = TrinucleotideContext('chr11.txt', 'chrXI', start, scores, early_2, middle_2, dict10[0], 3, dict10[1])
dict12 = TrinucleotideContext('chr12.txt', 'chrXII', start, scores, early_2, middle_2, dict11[0], 3, dict11[1])
dict13 = TrinucleotideContext('chr13.txt', 'chrXIII', start, scores, early_2, middle_2, dict12[0], 3, dict12[1])
dict14 = TrinucleotideContext('chr14.txt', 'chrXIV', start, scores, early_2, middle_2, dict13[0], 3, dict13[1])
dict15 = TrinucleotideContext('chr15.txt', 'chrXV', start, scores, early_2, middle_2, dict14[0], 3, dict14[1])
counts = TrinucleotideContext('chr16.txt', 'chrXVI', start, scores, early_2, middle_2, dict15[0], 3, dict15[1])




total_dicounts = 0
early_dicounts = 0
middle_dicounts = 0
late_dicounts = 0
print('Early: ')

for key in counts[1][0].keys():
    print(key + ": " + str(counts[1][0][key]))
    total_dicounts = total_dicounts + counts[1][0][key]
    early_dicounts = early_dicounts + counts[1][0][key]
    
print('Middle: ')

for key in counts[1][1].keys():
    print(key + ": " +str(counts[1][1][key]))
    total_dicounts = total_dicounts + counts[1][1][key]
    middle_dicounts = middle_dicounts + counts[1][1][key]
    
print('Late: ')

for key in counts[1][2].keys():
    print(key + ": " + str(counts[1][2][key]))
    total_dicounts = total_dicounts + counts[1][2][key]
    late_dicounts = late_dicounts + counts[1][2][key]


print(early_dicounts/total_dicounts)
print(middle_dicounts/total_dicounts)
print(late_dicounts/total_dicounts)

early_Trinuc = {}
middle_Trinuc = {}
late_Trinuc = {}

chr_dicts = [early_Trinuc, middle_Trinuc, late_Trinuc]
chr_di_dict = {0:{'CC':0, 'CT':0, 'TC':0, 'TT':0}, 1:{'CC':0, 'CT':0, 'TC':0, 'TT':0}, 2:{'CC':0, 'CT':0, 'TC':0, 'TT':0}}
chr16_counts =   TrinucleotideContext('chr16.txt', 'chrXVI', start, scores, early_2, middle_2, chr_dicts, 3, chr_di_dict)  
#file = open('chrXVI_dipyrimidines', 'w+')
#PrintTotalDipyrimidineCounts(chr16_counts[1], file)

def GetNewMutationRates(mutationRates, count):
    early_totals = count[0] 
    middle_totals = count[1]
    late_totals = count[2] 

    for key in mutationRates.keys():
        
    #mutated over total for each codon
        if key in late_totals.keys() and key in middle_totals.keys() and key in early_totals.keys():
            mutationRates[key] = mutationRates[key]/(count[0][key] + count[1][key] + count[2][key])
        elif key not in late_totals.keys():
            mutationRates[key] = mutationRates[key]/(count[0][key] + count[1][key]) 
        elif key not in middle_totals.keys():
            mutationRates[key] = mutationRates[key]/(count[0][key] + count[2][key])
        elif key not in early_totals.keys():
            mutationRates[key] = mutationRates[key]/(count[1][key] + count[2][key])
    return mutationRates

new_mutationRates = GetNewMutationRates(mutationRates, counts[0])
new_mutationRates16 = GetNewMutationRates(mutationRates16, counts[0])
new_mutationRates26 = GetNewMutationRates(mutationRates26, counts[0])
new_mutationRates30 = GetNewMutationRates(mutationRates30, counts[0])
new_mutationRatesUVB = GetNewMutationRates(mutationRatesUVB, counts[0])
new_mutationRatesB16 = GetNewMutationRates(mutationRatesB16, counts[0])
new_mutationRatesB26 = GetNewMutationRates(mutationRatesB26, counts[0])
new_mutationRatesB30 = GetNewMutationRates(mutationRatesB30, counts[0])
new_mutationRatesUVA = GetNewMutationRates(mutationRatesUVA, counts[0])
new_mutationRatesOg =  GetNewMutationRates(mutationRatesOg, counts[0])
new_mutationRatesr1 = GetNewMutationRates(mutationRatesr1, counts[0])
new_mutationRatesr2 = GetNewMutationRates(mutationRatesr2, counts[0])
#new_mutationRatesa1 = GetNewMutationRates(mutationRatesa1, counts[0])
#new_mutationRatesb1 = GetNewMutationRates(mutationRatesb1, counts[0])
#new_mutationRatesd1 = GetNewMutationRates(mutationRatesd1, counts[0])
#new_mutationRatese1 = GetNewMutationRates(mutationRatese1, counts[0])
#new_mutationRatesa2 = GetNewMutationRates(mutationRatesa2, counts[0])
#new_mutationRatesb2 = GetNewMutationRates(mutationRatesb2, counts[0])
#new_mutationRatesd2 = GetNewMutationRates(mutationRatesd2, counts[0])
#new_mutationRatese2 = GetNewMutationRates(mutationRatese2, counts[0])
#new_mutationRatesa3 = GetNewMutationRates(mutationRatesa3, dinuc_counts[0])
#new_mutationRatesb3 = GetNewMutationRates(mutationRatesb3, dinuc_counts[0])
#new_mutationRatesd3 = GetNewMutationRates(mutationRatesd3, dinuc_counts[0])
#new_mutationRatese3 = GetNewMutationRates(mutationRatese3, dinuc_counts[0])


def FindExpected(mutationRates, count, file):
    early_totals = count[0] 
    middle_totals = count[1] 
    late_totals = count[2] 
    expected_early = 0
    expected_middle = 0
    expected_late = 0
    f = open(file, 'w+')
    for key in mutationRates.keys():
        if key in early_totals.keys():    
            expected_early = expected_early + (early_totals[key] * mutationRates[key])
        if key in middle_totals.keys():
            expected_middle = expected_middle + (middle_totals[key] * mutationRates[key])
        if key in late_totals.keys():
            expected_late = expected_late + (late_totals[key] * mutationRates[key])
    total = expected_early + expected_middle + expected_late
    f.write('Early')
    f.write('\n')
    for key in early_totals.keys():
       
        f.write(key + ': ' + str((early_totals[key] )))
        f.write('\n')
       
    f.write('Middle')
    f.write('\n')
    for key in middle_totals.keys():
        f.write(key + ': ' + str((middle_totals[key])))
        f.write('\n')
       
    f.write('Late')
    f.write('\n')
    for key in late_totals.keys():
        
        f.write(key + ': ' + str((late_totals[key])))
        f.write('\n')
        
    print(total) 
    f.close()     

    percentage_early = expected_early/total
    percentage_middle = expected_middle/total
    percentage_late = expected_late/total
    print(percentage_early)
    print(percentage_middle)
    print(percentage_late)
    print(expected_early)
    print(expected_middle)
    print(expected_late)
    return[percentage_early, percentage_middle, percentage_late, expected_early, expected_middle, expected_late]

FindExpected(new_mutationRates, counts[0], 'Rates.txt')
FindExpected(new_mutationRates16, counts[0], 'Rates16.txt')
FindExpected(new_mutationRatesr1, counts[0], 'Rates_r1.txt')
FindExpected(new_mutationRatesr2, counts[0], 'Rates_r2.txt')
FindExpected(new_mutationRates30, counts[0], 'Rates30.txt')
FindExpected(new_mutationRates26, counts[0], 'Rates26.txt')
FindExpected(new_mutationRatesOg, counts[0], 'RatesOg.txt')
FindExpected(new_mutationRatesUVA, counts[0], 'RatesA.txt')
FindExpected(new_mutationRatesUVB, counts[0], 'Rates_UVB.txt')
FindExpected(new_mutationRatesB16, counts[0], 'Rates_UVB16.txt')
FindExpected(new_mutationRatesB26, counts[0], 'Rates_UVB26.txt')
FindExpected(new_mutationRatesB30,  counts[0], 'Rates_UVB30.txt')
#WT_deletions = FindExpected(new_mutationRatesa2, counts[0], 'Rates_deletions.txt')
#rad16_deletions = FindExpected(new_mutationRatesb1, counts[0], 'Rates16_deletions.txt')
#rad26_deletions =FindExpected(new_mutationRatesd2, counts[0], 'Rates26_deletions.txt')
#rad30_deletions = FindExpected(new_mutationRatese2, counts[0], 'Rates30_deletions.txt')
#WT_insertions = FindExpected(new_mutationRatesa3, dinuc_counts[0], 'Rates_insertions.txt')
#rad26_insertions =FindExpected(new_mutationRatesd3, dinuc_counts[0], 'Rates26_insertions.txt')
#rad16_insertions = FindExpected(new_mutationRatesb3, dinuc_counts[0], 'Rates16_insertions.txt')
#rad30_insertions = FindExpected(new_mutationRatese3, dinuc_counts[0], 'Rates30_insertions.txt')


#print(WT_deletions[3] + WT_insertions[3])
#print(WT_deletions[4] + WT_insertions[4])
#print(WT_deletions[5] + WT_insertions[5])
#print(rad16_deletions[3] + rad16_insertions[3])
#print(rad16_deletions[4] + rad16_insertions[4])
#print(rad16_deletions[5] + rad16_insertions[5])
#print(rad26_deletions[3] + rad26_insertions[3])
#print(rad26_deletions[4] + rad26_insertions[4])
#print(rad26_deletions[5] + rad26_insertions[5])
#print(rad30_deletions[3] + rad30_insertions[3])
#print(rad30_deletions[4] + rad30_insertions[4])
#print(rad30_deletions[5] + rad30_insertions[5])
#FindExpected(mutationRatesc,  counts[0], 'Rates16_complexsubs.txt')



expected_percentages = FindExpected(new_mutationRates, chr16_counts[0], 'Rates_chr16.txt')
expected_percentages16 = FindExpected(new_mutationRates16, chr16_counts[0], 'Rates16_chr16.txt')
expected_percentagesUVB = FindExpected(new_mutationRatesUVB, chr16_counts[0], 'Rates_UVB_chr16.txt')
#expected_percentages26 = FindExpected(new_mutationRates26, chr1_counts[0], 'Rates26_chr1.txt')
#expected_percentagesB26 = FindExpected(new_mutationRatesB26, chr1_counts[0], 'Rates26_UVB_chr1.txt')
expected_percentagesB16 = FindExpected(new_mutationRatesB16, chr16_counts[0], 'Rates_UVB16_chr16.txt')
