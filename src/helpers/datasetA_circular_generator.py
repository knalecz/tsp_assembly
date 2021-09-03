import os
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


read_len_before_mut = 3000
overlap_lengths = [2700, 2500, 2300]
errors_percentage = [0, 0.5, 1, 1.5]
result_path = "/home/knalecz/Pulpit/QuantumComputing/quantum_gsp/data/datasetAprim"
reads_num = 50
prefix = f'random_circular'
cov=15

   
if __name__=='__main__':
    initial_seq_len = int((read_len_before_mut*reads_num)/cov)
    seq = ''.join(random.choices('ACGT', k=initial_seq_len))    
    print(initial_seq_len)
    seq += seq
        
    for error in errors_percentage:
        for overlap in overlap_lengths:
            reads = []
            for read_num in range(reads_num):
                pos = read_num*(read_len_before_mut-overlap)
                
                if pos > initial_seq_len:
                    pos = pos % initial_seq_len
                                    
                read = seq[pos:pos+read_len_before_mut]
                
                if len(read) < 10:
                    print(f'read_num = {read_num}')
                    print(f'pos = {pos}; pos+read_len_before_mut = {pos+read_len_before_mut}')
                    print(f'read = {read}')
                
                output = os.popen(f'echo {read} | msbar -filter -auto -count {int((error*read_len_before_mut)/100)} -point 1').read()
                read = ''.join(output.split('\n')[1:-1])
                read_len_after_mut = len(read)

                reads.append(SeqRecord(Seq(read),
                                       id=f'{prefix}_{read_num+1:09d}_L{read_len_after_mut:09d}:{pos:09d}-{pos+read_len_after_mut-1:09d}:F',
                                       description=''))          
            result_filename = prefix + f'.L{initial_seq_len}.er{error:.1f}.ovl{overlap}.cov{cov}.rl{read_len_before_mut}.fasta'
            SeqIO.write(reads, result_path + "/" + result_filename, "fasta")


