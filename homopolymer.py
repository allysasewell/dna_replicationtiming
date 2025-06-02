def GetReverseComplement(str):
    if str == 'NA':
        return 
    if str == '-':
        return '-'
    dict = {}
    dict['A'] = 'T'
    dict['C'] = 'G'
    dict['G'] = 'C'
    dict['T'] = 'A'
    dict['-'] = '-'
    dict['*'] = '*'
    str2 = ''
    for i in range(0, len(str)):
        str2 = str2 + dict[str[i]]
    return str2[::-1]




def GetSequence(file):
    f1 = open(file)
    sequence1 = f1.read()
    sequence1 = str(sequence1.strip().upper())
    sequence1 = ''.join(sequence1.split('\n'))
    return sequence1


sequences = []
sequence1 = GetSequence('chr1.txt')
sequences.append(sequence1)
sequence2 = GetSequence('chr2.txt')
sequences.append(sequence2)
sequence3 = GetSequence('chr3.txt')
sequences.append(sequence3)
sequence4 = GetSequence('chr4.txt')
sequences.append(sequence4)
sequence5 = GetSequence('chr5.txt')
sequences.append(sequence5)
sequence6 = GetSequence('chr6.txt')
sequences.append(sequence6)
sequence7 = GetSequence('chr7.txt')
sequences.append(sequence7)
sequence8 = GetSequence('chr8.txt')
sequences.append(sequence8)
sequence9 = GetSequence('chr9.txt')
sequences.append(sequence9)
sequence10 = GetSequence('chr10.txt')
sequences.append(sequence10)
sequence11 = GetSequence('chr11.txt')
sequences.append(sequence11)
sequence12 = GetSequence('chr12.txt')
sequences.append(sequence12)
sequence13 = GetSequence('chr13.txt')
sequences.append(sequence13)
sequence14 = GetSequence('chr14.txt')
sequences.append(sequence14)
sequence15 = GetSequence('chr15.txt')
sequences.append(sequence15)
sequence16 = GetSequence('chr16.txt')
sequences.append(sequence16)



sequence = "GGGGTTTTT"

def GetHomopolymerLength(sequence, mut_pos1, mut_pos, mut_position):
    #mut_position = int(len(sequence)/2) + 1
    mut_pos1 = int(mut_pos1)
    mut_pos = int(mut_pos) 
    difference = mut_pos - mut_pos1 - 1
    mut_position2 = mut_position - difference
    end_position = -1
    lengths = {}
    for i in range(0, len(sequence)):
        flag = True
        length = 1
        n = 0
        while flag == True and (i + n + 1) < len(sequence):
            if sequence[i + n] == sequence[i + n + 1]:
                end_position = i + n + 2
                length = length + 1 
                n = n + 1
            else:
                flag = False
                if mut_position <= end_position  and mut_position >= end_position - (length - 1 ) or mut_position2 <= end_position  and mut_position2 >= end_position - (length - 1 ):
                     lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length
    
        if (mut_position ) <= end_position   and (mut_position ) > end_position - (length - 1 ) or mut_position <= end_position  and mut_position >= end_position - (length - 1 ):
            lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position) )] = length
           
    if len(lengths) == 0:
        return 1, (mut_pos, mut_pos)
    for key in lengths:
        if lengths[key] == max(lengths.values()):
            return lengths[key], key
        


with open("UV_yeast_Merged_SNP_DIPs_sorted_anz2_tandem.txt") as f:
    data = f.read()
lines = [s.strip().split() for s in data.splitlines()]
genotypes = {'WT': 0, 'rad16':0, 'rad26':0, 'rad30': 0}
WT_chromosomes = {}
rad16_chromosomes = {}
rad26_chromosomes = {}
rad30_chromosomes = {}
WT_indels = 0
rad16_indels = 0
rad26_indels = 0
rad30_indels = 0
genotype = []
background = []
chromosome = []
position1 = []
position2 = []
mutation_type = []
length = []
allele = []
mutation = []
frequency = []
left_reference = []
right_reference = []
left_consensus = []
right_consensus = []
trinucleotide = []
strand = []
isolate = []
isolate_list = {}
lines.pop(0)
for line in lines:
   
    if 'Mito' not in line[6] and 'ySR128' not in line[3] and line[0] != 'dRad16_25':#and line[10] == "DIP" or line[10] == "Complex DIP":
        isolate.append(line[0])
    
        genotype.append(line[4])
    
        mutation_type.append(line[10])
        if 'WT' in line[4]:
            genotypes['WT'] = genotypes['WT'] + 1
        if 'rad16' in line[4]:
            genotypes['rad16'] = genotypes['rad16'] + 1
        if 'rad26' in line[4]:
            genotypes['rad26'] = genotypes['rad26'] + 1
        if 'rad30' in line[4]:
            genotypes['rad30'] = genotypes['rad30'] + 1
        chromosome.append('chr' + str(line[6]))
        if 'WT' in line[4]:
            if str(line[6]) in WT_chromosomes.keys():
                WT_chromosomes[str(line[6])] = WT_chromosomes[str(line[6])] + 1

            else:
                WT_chromosomes[str(line[6])] = 1
        if 'rad16' in line[4]:
            if str(line[6]) in rad16_chromosomes.keys():
                rad16_chromosomes[str(line[6])] = rad16_chromosomes[str(line[6])] + 1
            else:
                rad16_chromosomes[str(line[6])] = 1
        if 'rad26' in line[4]:
            if str(line[6]) in rad26_chromosomes.keys():
                rad26_chromosomes[str(line[6])] = rad26_chromosomes[str(line[6])] + 1
            else:
                rad26_chromosomes[str(line[6])] = 1
        if 'rad30' in line[4]:
            if str(line[6]) in rad30_chromosomes.keys():
                rad30_chromosomes[str(line[6])] = rad30_chromosomes[str(line[6])] + 1
            else:
                rad30_chromosomes[str(line[6])] = 1
        position1.append(int(line[8]) - 1)
       
        position2.append(int(line[9]))
        length.append(line[11])
        allele.append(line[12])
        if len(line) >= 22: 
            if '-' not in  line[18] and line[18].isalpha() == False:
                    if float(line[18]) >= 90:
                        frequency.append('Homozygous')
                    else:
                        frequency.append('Heterozygous')
            elif float(line[19]) >= 90:
                    frequency.append('Homozygous')
            else:
                    frequency.append('Heterozygous')
            left_reference.append(line[23])
            right_reference.append(line[24])
            left_consensus.append(line[25])
            right_consensus.append(line[26])
            if line[27].isalpha() == False and len(line[24]) <= 10:
                    if line[24][1] == "G" or line[24][1] == "A":
                        #trinucleotide.append(GetReverseComplement(line[24]))
                        strand.append('-')
                        if '-' not in line[13]:
                            mutation.append(GetReverseComplement(line[13]))
                        else:
                            mutation.append(line[13])
                    else:
                        #trinucleotide.append(line[24])
                        mutation.append(line[13])
                        strand.append('+')
            elif len(line[27]) <= 10:
                    if line[27][1] == "A" or line[27][1] == "G":
                        #trinucleotide.append(GetReverseComplement(line[27]))
                        strand.append('-')
                        if '-' not in line[13]:
                            mutation.append(GetReverseComplement(line[13]))
                        else:
                            mutation.append(line[13])
                    else:
                        #trinucleotide.append(line[27])
                        mutation.append(line[13])
                        strand.append('+')
            else:
                    if line[28][1] == "A" or line[28][1] == "G":
                        #trinucleotide.append(GetReverseComplement(line[28]))
                        strand.append('-')
                        if '-' not in line[13]:
                            mutation.append(GetReverseComplement(line[13]))
                        else:
                            mutation.append(line[13])
                    else:
                        #trinucleotide.append(line[28])
                        mutation.append(line[13])
                        strand.append('+')
        


    
        else:
            frequency.append('NA')
            left_reference.append('NA')
            right_reference.append('NA')
            left_consensus.append('NA')
            right_consensus.append('NA')
            #trinucleotide.append('NA')
            mutation.append(line[13])
            if allele[len(allele) - 1] == 'A' or allele[len(allele) - 1] == 'G':
                strand.append('-')
            else:
                strand.append('+')


for pos in range(0, len(position2)):
    #if mutation_type_2[pos] == 'Insertion' or mutation_type_2[pos] == 'Deletion':# or mutation_type_2[pos] == 'Deletion' and len(allele2[pos]) == 1:
        if chromosome[pos] == 'chrI':
            sequence = GetSequence('chr1.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrII':
            sequence = GetSequence('chr2.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrIII':
            sequence = GetSequence('chr3.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrIV':
            sequence = GetSequence('chr4.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrV':
            sequence = GetSequence('chr5.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrVI':
            sequence = GetSequence('chr6.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrVII':
            sequence = GetSequence('chr7.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrVIII':
            sequence = GetSequence('chr8.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrIX':
            sequence = GetSequence('chr9.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrX':
            sequence = GetSequence('chr10.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXI':
            sequence = GetSequence('chr11.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXII':
            sequence = GetSequence('chr12.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXIII':
            sequence = GetSequence('chr13.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXIV':
            sequence = GetSequence('chr14.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXV':
            sequence = GetSequence('chr15.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXVI':
            sequence = GetSequence('chr16.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        else:
            trinucleotide.append('NA')
            strand.append('+')
for m in range(0, len(mutation)):
    if strand[m] == '-':
        mutation[m] = GetReverseComplement(mutation[m])


for i1 in range(0, len(isolate)):

    if (isolate[i1], genotype[i1]) not in isolate_list.keys():
            isolate_list[(isolate[i1], genotype[i1])] = [(position2[i1], 'chr' + str(chromosome[i1]), trinucleotide[i1], mutation[i1], strand[i1], mutation_type[i1], isolate[i1])]
    else:
        isolate_list[(isolate[i1], genotype[i1])] .append((position2[i1], 'chr' + str(chromosome[i1]), trinucleotide[i1], mutation[i1], strand[i1], mutation_type[i1], isolate[i1]))


f1a = open("WT_substitutions_sorted_tandem.bed", 'w+')
f1b = open("rad16_substitutions_sorted_tandem.bed", 'w+')
f1c = open("rad26_substitutions_sorted_tandem.bed", 'w+')
f1d = open("rad30_substitutions_sorted_tandem.bed", 'w+')


def WriteChr(chr_name, file1a, file1b, file1c, file1d, chromosome, position1, position2,
             trinucleotide, allele, mutation, genotype, mutation_type ,isolate, strand):
    file_data = []
    chromosome1_dict = {}
    chromosome1 =[]
    for x in range(0, len(mutation_type)):
        if mutation_type[x] == 'SNV' or mutation_type[x] == 'SNP':
                if trinucleotide[x].isalpha():# len(trinucleotide[x]) == 3 and trinucleotide[x].isalpha():
                    if chromosome[x] == chr_name: 
                        chromosome1_dict[int(position1[x]), isolate[x]] = x
                


    for key1 in chromosome1_dict.keys():
        chromosome1.append((int(key1[0]), key1[1]))
    sorted_chromosome1= sorted(chromosome1)
    #sorted_chromosome1.sort()
    for y in range(0, len(sorted_chromosome1)):
        file_data.append(chromosome[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + str(int(position1[chromosome1_dict[sorted_chromosome1[y]]])) + '\t' + str(position2[chromosome1_dict[sorted_chromosome1[y]]]) + '\t' + trinucleotide[chromosome1_dict[sorted_chromosome1[y]]]  + '\t' + mutation[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + strand[chromosome1_dict[sorted_chromosome1[y]]] + '\t' )
        #print(chromosome[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + str(int(position1[chromosome1_dict[sorted_chromosome1[y]]]) - 1) + '\t' + position2[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + trinucleotide[chromosome1_dict[sorted_chromosome1[y]]]  + '\t' + mutation[chromosome1_dict[sorted_chromosome1[y]]])
    for fd in range(0, len(file_data)):
        
            if genotype[chromosome1_dict[sorted_chromosome1[fd]]] == 'WT':
                file1a.write(file_data[fd])
               
                file1a.write(isolate[chromosome1_dict[sorted_chromosome1[fd]]])
                file1a.write('\n')
                
            if 'rad16' in genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad16' in  genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
                file1b.write(file_data[fd])
                
                file1b.write(isolate[chromosome1_dict[sorted_chromosome1[fd]]])
                file1b.write('\n')
              
            if 'rad26' in genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad26' in genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
    
                file1c.write(file_data[fd])
               
                file1c.write(isolate[chromosome1_dict[sorted_chromosome1[fd]]])
                file1c.write('\n')
              
            if 'rad30' in genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad30' in  genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
                file1d.write(file_data[fd])
              
                file1d.write(isolate[chromosome1_dict[sorted_chromosome1[fd]]])
                file1d.write('\n')
               

       

WriteChr('chrI',f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrII', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrIII',  f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrIV', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrIX', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrV', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrVI', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrVII', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrVIII', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrX', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXI', f1a, f1b, f1c, f1d,  chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXII', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXIII', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXIV', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand) 
WriteChr('chrXV', f1a, f1b, f1c, f1d,  chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXVI', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)

f1a.close()
f1b.close()
f1c.close()
f1d.close()


def MakePentaNucFile(file_1, file_2):
    new_data = file_1.read()
    new_chromosome = []
    new_position = []
    new_position1 = []
    new_mutation = []
    new_allele = []
    new_trinucleotide = []
    new_strand = []
    new_isolate = []
    lengths = []
    new_lines = [s.strip().split() for s in new_data.splitlines()]
    for line in new_lines:
        new_chromosome.append(line[0])
        new_position1.append(int(line[1]))
        new_position.append(int(line[2]))
        new_mutation.append(line[4])
        new_allele.append(line[3][1])
        new_trinucleotide.append(line[3])
        new_strand.append(line[5])
        new_isolate.append(line[6])

    #Add two more nucleotides after the 3' base
    pentanucleotide = []
    for pos in range(0, len(new_position)):
        if new_chromosome[pos] == 'chrI':
            sequence = GetSequence('chr1.txt')
            #if strand[pos] == '+':
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrII':
            sequence = GetSequence('chr2.txt')
            #if strand[pos] == '+':
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrIII':
            sequence = GetSequence('chr3.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
           
            


        elif new_chromosome[pos] == 'chrIV':
            sequence = GetSequence('chr4.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrV':
            sequence = GetSequence('chr5.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrVI':
            sequence = GetSequence('chr6.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrVII':
            sequence = GetSequence('chr7.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrVIII':
            sequence = GetSequence('chr8.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrIX':
            sequence = GetSequence('chr9.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrX':
            sequence = GetSequence('chr10.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrXI':
            sequence = GetSequence('chr11.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrXII':
            sequence = GetSequence('chr12.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrXIII':
            sequence = GetSequence('chr13.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrXIV':
            sequence = GetSequence('chr14.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrXV':
            sequence = GetSequence('chr15.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrXVI':
            sequence = GetSequence('chr16.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0] )
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            
            
    
    for i in range(0, len(new_position)):
        if new_strand[i] == '-':
            file_2.write(new_chromosome[i] + '\t' + str(int(new_position[i]) - 1) + '\t'  +str(new_position[i]) + '\t' + new_trinucleotide[i] + '\t' + new_mutation[i] + '\t' + new_strand[i]  + '\t'+ GetReverseComplement(pentanucleotide[i]) + '\t' +str(GetHomopolymerLength(GetReverseComplement(pentanucleotide[i]), new_position1[i], new_position[i], int(len(pentanucleotide[i])/2) + 1)[0]) + '\t' + new_isolate[i])
            file_2.write('\n')
        if new_strand[i] == '+':
            file_2.write(new_chromosome[i] + '\t' + str(int(new_position[i]) - 1)  + '\t' + str(new_position[i]) + '\t' + new_trinucleotide[i] + '\t' +  new_mutation[i] + '\t' + new_strand[i] + '\t' + pentanucleotide[i] + '\t' + str(GetHomopolymerLength(pentanucleotide[i], new_position1[i], new_position[i], int(len(pentanucleotide[i])/2) + 1)[0]) + '\t' + new_isolate[i])
            file_2.write('\n')

    return pentanucleotide, new_strand, new_chromosome, new_position, new_isolate, new_trinucleotide, new_mutation, new_allele, lengths, new_position1

def GetPolymerCounts(pentanucleotide,  new_chromosome, new_position, new_isolate, new_trinucleotide, new_mutation, new_allele, lengths):
    homopolymers = {}
    homopolymer_list = {}
    polymer_counts = []
    for j in range(0, 4):
        polymer_counts.append({})
    for i in range (0, len(lengths)):
        if len(new_mutation[i]) == 1:
            print(pentanucleotide[i])
            print(new_trinucleotide[i])
            key = lengths[i]
            
            #elif strand[i] == '-':
                #key = GetHomopolymerLength(GetReverseComplement(sequence), )[1]
            
            if key >= 4:
                homopolymers[(new_isolate[i], new_chromosome[i], new_position[i])] = key
                #homopolymer_list[key[1]] = [(new_chromosome[i], key[0][0], key[0][1])]
                print(new_trinucleotide[i])
                print(pentanucleotide[i])
                print(key)

            if new_allele[i] == 'C' or new_allele[i] == 'G':

                if key in polymer_counts[0].keys():
                    polymer_counts[0][key] = polymer_counts[0][key] + 1
                else:
                    polymer_counts[0][key] = 1
            if new_allele[i] == 'A' or new_allele[i] == 'T':

                if key in polymer_counts[1].keys():
                    polymer_counts[1][key] = polymer_counts[1][key] + 1
                else:
                    polymer_counts[1][key] = 1
            if new_mutation[i] == 'C' or new_mutation[i] == 'G':

                if key in polymer_counts[2].keys():
                    polymer_counts[2][key] = polymer_counts[2][key] + 1
                else:
                    polymer_counts[2][key] = 1
            if new_mutation[i] == 'A' or new_mutation[i] == 'T':

                if key in polymer_counts[3].keys():
                    polymer_counts[3][key] = polymer_counts[3][key] + 1
                else:
                    polymer_counts[3][key] = 1
    return polymer_counts, homopolymers, homopolymer_list

f2 = open("WT_substitutions_sorted_tandem.bed")
f2a = open("WT_substitutions_sorted_multinucleotide.bed", 'w+')
f3 = open("rad16_substitutions_sorted_tandem.bed")
f3a = open("rad16_substitutions_multinucleotide.bed", 'w+')
#f4 = open("rad26_UV_bothruns_Muts_allSNVs_sorted.bed")
#f4a = open("rad26_substitutions_multinucleotide.bed")

WT_data = MakePentaNucFile(f2, f2a)
WT_isolate = WT_data[4]
WT_chromosome = WT_data[2]
WT_position1 = WT_data[9]
WT_position = WT_data[3]
WT_trinucleotide = WT_data[5]
WT_mutation = WT_data[6]
WT_strand = WT_data[1]
WT_multinucleotide = WT_data[0] 

rad16_data = MakePentaNucFile(f3, f3a)
WT_isolate_list = {}
for i1 in range(0, len(WT_isolate)):

    if WT_isolate[i1] not in WT_isolate.keys():
            WT_isolate[WT_isolate[i1]] = [(WT_position[i1], 'chr' + str(WT_chromosome[i1]), WT_trinucleotide[i1], WT_mutation[i1], WT_strand[i1], WT_multinucleotide[i1], WT_isolate[i1], WT_position1[i1])]
    else:
        WT_isolate_list[isolate[i1]] .append((WT_position[i1], 'chr' + str(WT_chromosome[i1]), WT_trinucleotide[i1], WT_mutation[i1], WT_strand[i1], WT_multinucleotide[i1], WT_isolate[i1], WT_position1[i1]))
WT_polymer_counts = GetPolymerCounts(WT_data[0],  WT_data[2], WT_data[3], WT_data[4], WT_data[5], WT_data[6], WT_data[7], WT_data[8])[0]
rad16_polymer_counts = GetPolymerCounts(rad16_data[0],  rad16_data[2], rad16_data[3], rad16_data[4], rad16_data[5], rad16_data[6], rad16_data[7], rad16_data[8])[0]
f2.close()
f2a.close()
f3.close()
f3a.close()



def NormalizeHomopolymer(sequences, lengths, file):
    homopolymer_dict = {}
    homopolymer_counts = {}
    homopolymer_dict['11-15'] = 0
    homopolymer_dict['16-20'] = 0
    homopolymer_dict['21+'] = 0

    for l in range(0, len(lengths)):
        if lengths[l] not in homopolymer_counts.keys():
            homopolymer_counts[lengths[l]] = 1
            homopolymer_dict[lengths[l]] = 0
        else:
            homopolymer_counts[lengths[l]] = homopolymer_counts[lengths[l]] + 1

    homopolymer_dict[1] = 0
    for seq in sequences:
        sequence = seq
        if sequence[len(sequence) - 1] != sequence[len(sequence) - 2]:
            homopolymer_dict[1] = homopolymer_dict[1] + 1
        b = 0
        while b < len(sequence) - 1:
            polymer_range = {}
            length = 1
            if sequence[b] != 'N':# and b not in coordinates:
                i = 1
                polymer_range[0] = b
                while sequence[b] == sequence[b + i]:
                    length = length + 1
                    polymer_range[1] = b + i
                    if (b + i) < len(sequence) - 1:
                        i = i + 1
                    else:
                        break
     
                polymer_count = length
          
                if polymer_count in homopolymer_dict.keys():
                            homopolymer_dict[polymer_count] = homopolymer_dict[polymer_count] + length
                            
                           
                                  

    

                else:
                            homopolymer_dict[polymer_count] = length
                
                if polymer_count >= 11:
                    if polymer_count >= 21:
                        homopolymer_dict['21+'] = homopolymer_dict['21+'] + length
                    elif polymer_count >= 16:
                        homopolymer_dict['16-20'] = homopolymer_dict['16-20'] + length
                    else:
                        homopolymer_dict['11-15'] = homopolymer_dict['11-15'] + length
                            
                b = b + length 
            else:
                b = b + 1
                    
    for key in homopolymer_dict.keys():
        #homopolymer_dict[key] = homopolymer_counts[key]/homopolymer_dict[key]
        file.write(str(key) + ':' + str(homopolymer_dict[key]))
        file.write('\n')
    return homopolymer_dict




def PrintFrequencies(homopolymer_dict, lengths, mutation_type, allele, mutation, file, file2):
    homopolymer_counts = {}
    homopolymer_counts['21+'] = 0
    homopolymer_counts['16-20'] = 0
    homopolymer_counts['11-15'] = 0
    for l in range(0, len(lengths)):
        if allele[l] != 'NA' and lengths[l] != 'NA':
            #if (len(allele[l]) == 1 and mutation_type[l] == "deletion") or (len(mutation[l]) == 1 and mutation_type[l] == "insertion") :
                if lengths[l] not in homopolymer_counts.keys():
                    homopolymer_counts[lengths[l]] = 1
                    #homopolymer_dict[lengths[l]] = 0
                else:
                    homopolymer_counts[lengths[l]] = homopolymer_counts[lengths[l]] + 1
                if int(lengths[l]) >= 11:
                    if int(lengths[l]) >= 21:
                            homopolymer_counts['21+'] = homopolymer_counts['21+'] + 1
                    elif int(lengths[l]) >= 16:
                            homopolymer_counts['16-20'] = homopolymer_counts['16-20'] + 1
                    else:
                            homopolymer_counts['11-15'] = homopolymer_counts['11-15'] + 1
                
    for key in homopolymer_counts.keys():
        if homopolymer_dict[key] != 0:
            homopolymer_counts[key] = homopolymer_counts[key]/homopolymer_dict[key]
            file.write(str(key) + ':' + str(homopolymer_counts[key]))
            file.write('\n')
   

def index(tuple):
    return tuple[1]   

def FindClusters( isolate_list):
    lengths = {}
    clusters = {}
    nonclusters = {}
    total_clusters = 0
    for key in isolate_list.keys():
            sorted_isolate_list = sorted(isolate_list[key])
            sorted_list = sorted(sorted_isolate_list, key = index)
            if sorted_list[0][1] == sorted_list[1][1] and (int(sorted_list[1][0]) <= int((sorted_list[ 0][0])) + 10):
                if sorted_list[0] not in clusters.keys():
                    clusters[sorted_list[0]] = [sorted_list[0]]
                   
                if sorted_list[1] not in clusters.keys():
                    clusters[sorted_list[1]] = [sorted_list[1]]
                    
            else:#if strand[0] == '+':
                end_position = GetHomopolymerLength(sorted_list[0][5], int(sorted_list[0][7]), int(sorted_list[0][0]), int(len(sorted_list[0][5])/2) + 1)[1][1]
                if int(sorted_list[1][0]) <= (end_position + 10) and int(sorted_list[0][1]) == sorted_list[1][1]:
                    if sorted_list[0] not in clusters.keys():
                        clusters[sorted_list[0]] = [sorted_list[0]]
                    if sorted_list[1] not in clusters.keys():
                        clusters[sorted_list[1]] = [sorted_list[1]]
                start_position = GetHomopolymerLength(sorted_list[1][5], int(sorted_list[1][7]), int(sorted_list[1][0]), int(len(sorted_list[1][5])/2) + 1)[1][0]
                if start_position <= (sorted_list[0][0] + 10) and int(sorted_list[0][1]) == sorted_list[1][1]:
                    if sorted_list[0] not in clusters.keys():
                        clusters[sorted_list[0]] = [sorted_list[0]]
                    if sorted_list[1] not in clusters.keys():
                        clusters[sorted_list[1]] = [sorted_list[1]]

            if sorted_list[len (sorted_list) - 2][1] == sorted_list[len (sorted_list) - 1][1] and (int(sorted_list[len (sorted_list) - 1][0]) <= int((sorted_list[ len (sorted_list) - 2][0])) + 10):
                if sorted_list[len(sorted_list) - 1] not in clusters.keys():
                    clusters[sorted_list[len(sorted_list) - 1]] = [sorted_list[len(sorted_list) - 1]]
                  
                if sorted_list[len(sorted_list) - 2] not in clusters.keys() and (sorted_list[len(sorted_list) - 2][7] == "DIP" or sorted_list[len(sorted_list) - 2][7] == "Complex DIP" or sorted_list[len(sorted_list) - 2][7] == "Deletion" or sorted_list[len(sorted_list) - 2][7] == "Insertion"):
                    clusters[sorted_list[len(sorted_list) - 2]] = [sorted_list[len(sorted_list) - 2]]
                    
            else:#if strand[0] == '+':
                end_position = GetHomopolymerLength(sorted_list[len(sorted_list) - 2][5], int(sorted_list[len(sorted_list) - 2][7]), int(sorted_list[len(sorted_list) - 2][0]), int(len(sorted_list[len(sorted_list) - 2][5])/2) + 1)[1][1]
                if int(sorted_list[len(sorted_list) - 1][0]) <= (end_position + 10) and int(sorted_list[len(sorted_list) - 1][1]) == sorted_list[len(sorted_list) - 2][1]:
                    if sorted_list[len(sorted_list) - 1] not in clusters.keys():
                        clusters[sorted_list[len(sorted_list) - 1]] = [sorted_list[len(sorted_list) - 1]]
                    if sorted_list[len(sorted_list) - 2] not in clusters.keys():
                        clusters[sorted_list[len(sorted_list) - 2]] = [sorted_list[len(sorted_list) - 2]]
                start_position = GetHomopolymerLength(sorted_list[len(sorted_list) - 1][5], int(sorted_list[len(sorted_list) - 1][7]), int(sorted_list[len(sorted_list) - 1][0]), int(len(sorted_list[len(sorted_list) - 1][5])/2) + 1)[1][0]
                if start_position <= (sorted_list[len(sorted_list) - 2][0] + 10) and int(sorted_list[len(sorted_list) - 2][1]) == sorted_list[len(sorted_list) - 1][1]:
                    if sorted_list[len(sorted_list) - 1] not in clusters.keys():
                        clusters[sorted_list[len(sorted_list) - 1]] = [sorted_list[len(sorted_list) - 1]]
                    if sorted_list[1] not in clusters.keys():
                        clusters[sorted_list[len(sorted_list) - 2]] = [sorted_list[len(sorted_list) - 2]]
                   
            for j in range (1,len (sorted_list) - 1):
                    if (sorted_list[j][1] == sorted_list[j - 1][1] and (int(sorted_list[j][0]) <= int((sorted_list[j - 1][0])) + 10)) or (sorted_list[j + 1][1] == sorted_list[j][1] and (int(sorted_list[j + 1][0]) <= int((sorted_list[j][0])) + 10)): #or ((isolate_list[key][j + n][6], int(isolate_list[key][j + n][0]), isolate_list[key][j + n][1]) in homopolymers.keys() and (int(isolate_list[key][j + n][0]) <= int((isolate_list[key][j + n - 1][0])) + 10 + homopolymers[key])) : 
                            if sorted_list[j] not in clusters.keys():
                                #if j == 0 or isolate_list[key][j - n] not in clusters.keys():
                                clusters[(sorted_list[j][6], sorted_list[j][1], sorted_list[j][0])] = [sorted_list[j]]
                    else:
                        end_position = GetHomopolymerLength(sorted_list[j][5], int(sorted_list[j][7]), int(sorted_list[j][0]), int(len(sorted_list[j][5])/2) + 1)[1][1]
                        if int(sorted_list[j + 1][0]) <= (end_position + 10) and int(sorted_list[j][1]) == sorted_list[j + 1][1]:
                            if sorted_list[j] not in clusters.keys():
                                clusters[sorted_list[j]] = [sorted_list[j]]
                            if sorted_list[j + 1] not in clusters.keys():
                                clusters[sorted_list[j + 1]] = (sorted_list[j + 1])
                        start_position = GetHomopolymerLength(sorted_list[j][5], int(sorted_list[j][7]), int(sorted_list[j][0]), int(len(sorted_list[j][5])/2) + 1)[1][0]
                        if start_position <= (sorted_list[j - 1][0] + 10) and int(sorted_list[j - 1][1]) == sorted_list[j][1]:
                            if sorted_list[j - 1] not in clusters.keys():
                                clusters[sorted_list[j - 1]] = [sorted_list[j - 1]]
                            if sorted_list[j] not in clusters.keys():
                                clusters[sorted_list[j]] = [sorted_list[j]]           
                    #elif (sorted_list[j][1] == sorted_list[j + 1][1]) and(sorted_list[j][6], int(sorted_list[j][0]), sorted_list[j][1]) in homopolymers.keys() and (int(sorted_list[j + 1][0]) <= int((sorted_list[j][0])) + 10 + (homopolymers[(sorted_list[j][6], int(sorted_list[j][0]), sorted_list[j][1])][0][1] - sorted_list[j][0])):
                        #if sorted_list[j] not in clusters.keys() and (sorted_list[j][7] == "DIP" or sorted_list[j][7] == "Complex DIP" or sorted_list[key][j][7] == "Deletion" or sorted_list[j][7] == "Insertion"):
                            #clusters[(sorted_list[j][6], sorted_list[j][1], sorted_list[j][0])] = [sorted_list[j]]
                        #if sorted_list[j + 1] not in clusters.keys() and (sorted_list[j + 1][7] == "DIP" or sorted_list[j + 1][7] == "Complex DIP" or sorted_list[j + 1][7] == "Deletion" or sorted_list[j + 1][7] == "Insertion"):
                            #clusters[(sorted_list[j + 1][6], sorted_list[j + 1][1], sorted_list[j + 1][0])] = [sorted_list[j + 1]]
                        
                    #elif (sorted_list[j][1] == sorted_list[j - 1][1]) and (sorted_list[j][6], sorted_list[j][0], sorted_list[j][1]) in homopolymers.keys() and (int(sorted_list[j][0]) <= int((sorted_list[j-1][0])) + 10 + (sorted_list[j][0] - homopolymers[sorted_list[j][6], int(sorted_list[j][0]), sorted_list[j][1]][0][0] )):
                        #if sorted_list[j] not in clusters.keys() and (sorted_list[j][7] == "DIP" or sorted_list[j][7] == "Complex DIP" or sorted_list[j][7] == "Deletion" or sorted_list[j][7] == "Insertion"):
                            #clusters[(sorted_list[j][6], sorted_list[j][1], sorted_list[j][0])] = [sorted_list[j]]
                        #if sorted_list[j - 1] not in clusters.keys() and (sorted_list[j - 1][7] == "DIP" or sorted_list[j - 1][7] == "Complex DIP" or sorted_list[j - 1][7] == "Deletion" or sorted_list[j - 1][7] == "Insertion"):
                            #clusters[(sorted_list[j - 1][6], sorted_list[j - 1][1], sorted_list[j - 1][0])] = [sorted_list[j - 1]]
        
   

    complex_count = []
    for key in homopolymers.keys():
        if key in clusters.keys():
            complex_count.append(key)
    
    return clusters