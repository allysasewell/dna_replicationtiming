


di_counts = {}

def DinucleotideContext(file, di_counts):
   
    ft = open(file)
    sequence = ft.read()
    sequence = str(sequence.strip().upper())
    sequence = ''.join(sequence.split('\n'))
    
    for i in range(1, len(sequence) - 1):
        key = sequence[i-1:i+1]
        if key in di_counts.keys():
            di_counts[key] = di_counts[key] + 1
        else:
            di_counts[key] = 1
       
        


                                         
    ft.close()
    return di_counts
di_dict1 = DinucleotideContext('chr1.txt',  di_counts)
di_dict2 = DinucleotideContext('chr2.txt',  di_dict1)
di_dict3 = DinucleotideContext('chr3.txt', di_dict2)
di_dict4 = DinucleotideContext('chr4.txt',  di_dict3)
di_dict5 = DinucleotideContext('chr5.txt', di_dict4)
di_dict6 = DinucleotideContext('chr6.txt',  di_dict5)
di_dict7 = DinucleotideContext('chr7.txt', di_dict6)
di_dict8 = DinucleotideContext('chr8.txt',  di_dict7)
di_dict9 = DinucleotideContext('chr9.txt',  di_dict8)
di_dict10 = DinucleotideContext('chr10.txt',  di_dict9)
di_dict11 = DinucleotideContext('chr11.txt', di_dict10)
di_dict12 = DinucleotideContext('chr12.txt',  di_dict11)
di_dict13 = DinucleotideContext('chr13.txt',  di_dict12)
di_dict14 = DinucleotideContext('chr14.txt',  di_dict13)
di_dict15 = DinucleotideContext('chr15.txt',  di_dict14)
dinuc_counts = DinucleotideContext('chr16.txt',  di_dict15)

