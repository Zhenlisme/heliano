#!_INTERPRETERPYTHON_PATH_

import os, re, subprocess, sys, argparse, shutil, random
from Bio import SeqIO
from multiprocessing.pool import ThreadPool
from collections import defaultdict
import pybedtools as BT

class Consensus_making:
    def __init__(self, genome, wkdir, represent_bed, process_num):
        self.genome = genome
        self.genome_dict = SeqIO.parse(genome, 'fasta')
        self.genome_dict = {k.id: k.seq.upper() for k in self.genome_dict}
        self.process_num=int(process_num)
        self.wkdir = wkdir
        self.repsenbed = represent_bed
        if not os.path.exists(self.wkdir):
            os.mkdir(self.wkdir)
            os.chdir(self.wkdir)
        else:
            os.chdir(self.wkdir)

        self.genome_size = 'Genome.size'
        genome_size = [[i, len(self.genome_dict[i])] for i in self.genome_dict]
        genome_size = sorted(genome_size, key=lambda x: x[0])
        with open(self.genome_size, 'w') as F:
            F.writelines([''.join([i[0], '\t', str(i[1]), '\n']) for i in genome_size])

    def cdhitest_clust(self, input_fa):
        cons_name = '.'.join([input_fa, 'reduce.temp'])
        cluster_file = '.'.join([cons_name, 'clstr'])
        run_cluster = subprocess.Popen(
            ['cd-hit-est', '-i', input_fa, '-o', cons_name, '-d', '0', '-aS', '0.8', '-aL', '0.8', '-c', '0.8', '-G', '1', '-g',
             '1', '-b', '500', '-T', str(self.process_num), '-M', '0'], stdout=subprocess.DEVNULL)
        run_cluster.wait()

        cluster_dict = {}
        with open(cluster_file, 'r') as F:
            for line in F:
                if line.startswith('>'):
                    cluster_name = line.strip('>\n').replace(' ', '_')
                else:
                    insertion_name = line.split('...')[0].split(', >')[1]
                    cluster_dict[insertion_name] = cluster_name

        if os.path.exists(cons_name):
            os.remove(cons_name)
        os.remove(cluster_file)
        reversed_cluster_dict = defaultdict(list)
        [reversed_cluster_dict[cluster_dict[key]].append(key) for key in cluster_dict]
        return reversed_cluster_dict

    def consencus(self, mfa):
        basename = os.path.basename(mfa).split('.')[0]
        ## need to cluster before making consensus
        cluster_dict = self.cdhitest_clust(mfa)   #{pairname:[insertion1, insertion2, ...]}
        mfa_dict = SeqIO.parse(mfa, 'fasta')
        mfa_dict = {k.id: k.seq.upper() for k in mfa_dict}

        for pairname in cluster_dict:
            insertion_name_list = cluster_dict[pairname]
            if len(insertion_name_list) < 2:
                conseq = str(mfa_dict[insertion_name_list[0]])
                self.consensus_dict[basename].append(conseq)
                continue

            submfa = ''.join(['Consensus/', basename, '.mfa'])
            subaln = ''.join([submfa, '.fa'])
            with open(submfa, 'w') as F:
                ## only first five insertions used for consensus making
                for insertion in insertion_name_list[:5]:
                    F.write(''.join(['>', insertion, '\n', str(mfa_dict[insertion]), '\n']))

            mul_aln_run = subprocess.Popen(['dialign2-2', '-n', '-fa', '-mask', submfa],
                                              stdout=subprocess.DEVNULL)
            mul_aln_run.wait()

            ## To remove non utf-8 codes
            with open(subaln, 'rb') as F:
                text = F.read().decode('utf-8', 'ignore')
            with open(subaln, 'w') as F:
                F.write(text)

            consensus_file = ''.join([submfa, '.con.fa'])
            consencus_task = subprocess.Popen(["cons", "-sequence", subaln, '-outseq', consensus_file],
                                              stdout = subprocess.DEVNULL)
            consencus_task.wait()

            consencus_seq = {}
            with open(consensus_file, 'r') as F:
                for line in F:
                    if line.startswith('>'):
                        key = line.strip('>\n')
                        consencus_seq[key] = ''
                    else:
                        consencus_seq[key] += line.rstrip()
            #os.remove(consensus_file)
            #os.remove(subaln)
            conseq = list(consencus_seq.values())[0]
            conseq = re.findall('[ATCG]{5,}.*[ATCG]{5,}', conseq)[0].replace('n', '').replace('N', '').replace('*', '').replace('x', '')  ## delete gap region
            self.consensus_dict[basename].append(conseq)
    
    def extract_seq(self, pairlist):
        subwkdir = 'Consensus'
        pairname = pairlist[0][3]
        pairfa = ''.join([subwkdir, '/', pairname, '.fa'])
        ## To extend both ends for 50 bp. 
        distance = 50
        with open(pairfa, 'w') as PF:
            for line in pairlist:
                chrmid, start, stop, name, score, strand = line[:6]
                seq_start = int(start) - distance if int(start) - distance > 0 else 0
                seq_stop = int(stop) + distance
                seq = self.genome_dict[chrmid][seq_start-1:seq_stop]
                if strand == '+':
                    PF.write(''.join(['>', chrmid, '-', str(seq_start), '-', str(seq_stop), '\n']))
                    PF.write(str(seq))
                    PF.write('\n')
                else:
                    PF.write(''.join(['>', chrmid, '-', str(seq_start), '-', str(seq_stop), '\n']))
                    PF.write(str(seq.reverse_complement()))
                    PF.write('\n')
        self.mfalist.append(pairfa)

    def main(self):
        sys.stdout.write('Begin to make consensus sequences.\n')
        subwkdir = 'Consensus/'
        if not os.path.exists(subwkdir):
            os.mkdir(subwkdir)
        else:
            shutil.rmtree(subwkdir)
            os.mkdir(subwkdir)
        representative_dict = defaultdict(list)
        with open(self.repsenbed, 'r') as F:
            for line in F:
                splitlines = line.rstrip().split('\t')
                representative_dict[splitlines[3]].append(splitlines)
        representative_list = [representative_dict[pairname] for pairname in representative_dict]

        ## To output the multiple sequences
        self.mfalist = []
        planpool = ThreadPool(int(self.process_num))
        for pairlist in representative_list:
            planpool.apply_async(self.extract_seq, args=(pairlist, ))
        planpool.close()
        planpool.join()

        self.consensus_dict = defaultdict(list)
        ## To make consensus sequences for each subfamily.
        planpool = ThreadPool(int(self.process_num))
        for mfa in self.mfalist:
            opfile=''
            planpool.apply_async(self.consencus, args=(mfa, ))
        planpool.close()
        planpool.join()

        with open('RC.representative.cons.fa', 'w') as F:
            for name in self.consensus_dict:
                init = 1
                for seq in self.consensus_dict[name]:
                    seqname = ''.join(['>', name, '.', str(init)])
                    F.write(''.join([seqname, '\n', seq, '\n']))
                    init += 1
        sys.stdout.write('Consensus sequences got constructed!.\n')

                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Making consensus for Helitron-like sequences. Please visit https://github.com/Zhenlisme/heliano/ for more information. Email us: zhen.li3@universite-paris-saclay.fr")
    parser.add_argument("-g", "--genome", type=str, required=True, help="The genome file in fasta format.")
    parser.add_argument("-r", "--repsenbed", type=str, required=True, help="The representative bed file.")
    parser.add_argument("-o", "--opdir", type=str, required=True, help="The output directory.")
    parser.add_argument("-n", "--process", type=int, default=2, required=False, help="Maximum of threads to be used.")
    parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0.2')
    Args = parser.parse_args()
    makeconsenus = Consensus_making(os.path.abspath(Args.genome), os.path.abspath(Args.opdir), os.path.abspath(Args.repsenbed), Args.process)
    makeconsenus.main()

