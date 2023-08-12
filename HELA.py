#!_INTERPRETERPYTHON_PATH_

import os, re, subprocess, sys, argparse, shutil, random, gc
from Bio import SeqIO
from multiprocessing.pool import ThreadPool
from collections import defaultdict
import pybedtools as BT

"""
This program is supposed to detect and classify different variants of Helitron-like elements: Helitron, Helentron and Helitron2. Please follow the homepage of this software for more information: https://github.com/Zhenlisme/HELA.   
"""

# define Structure_search class for motif identification, including stem loop and terminal inverted repeats.
class Structure_search:
    def __init__(self, genome, START=0):
        self.START = int(START)
        self.genome = genome
        self.genome_dict = SeqIO.parse(genome, 'fasta')
        self.genome_dict = {k.id: k.seq.upper() for k in self.genome_dict}
        self.maximum_length = sorted([len(self.genome_dict[i]) for i in self.genome_dict])[-1]

    def stem_loop(self, stem_loop_description, minus_tailone=1):
        ## Function to find hairpin structures. Start coord is 1 not 0
        rnabobopt = ''.join([os.path.basename(self.genome), '.stemloop.txt'])
        with open(rnabobopt, 'w') as rnabf:
            rnabob_program = subprocess.Popen(["rnabob", "-c", "-q", "-F", "-s", stem_loop_description, self.genome],
                                              stderr=subprocess.DEVNULL, stdout=rnabf)
            rnabob_program.wait()
        if not os.path.exists(rnabobopt):
            return []
        stem_loop_loc = []

        complement_dict = {'A': "T", "T": "A", "G": "C", "C": "G",
                           "K": "M", "M": "K", "Y": "R", "R": "Y", "S": "S", "W": "W",
                           "B": "V", "V": "B", "H": "D", "D": "H", "N": "N", "X": "X"}

        with open(rnabobopt, 'r') as F:
            for line in F:
                line = line.strip()
                if re.match('\d', line):
                    splitline = re.split('\s+', line)[:3]
                    chrid = splitline[2]
                    ## positive strand
                    if int(splitline[0]) < int(splitline[1]):
                        strand = '+'
                        length = int(splitline[1]) - int(splitline[0]) + 1
                        start = int(splitline[0]) + self.START
                        end = start + length - 1
                        end = end - 1 if minus_tailone else end  ## To remove the T nucleotide
                    ## negative strand
                    else:
                        strand = '-'
                        length = int(splitline[0]) - int(splitline[1]) + 1
                        start = int(splitline[0]) + self.START - length + 1
                        end = start + length - 1
                        start = start + 1 if minus_tailone else start  ## To remove the T nucleotide
                else:
                    seq = line.strip('|').split('|')
                    helix_seq1, loop_seq, helix_seq2, tail_seq = seq
                    stem_len = len(helix_seq1)
                    loop_len = len(loop_seq)
                    ## To revise the rnabob output. rnabob sometimes does not return as long as possible of helix. need to revise it.
                    midpoint = int(len(loop_seq) / 2)
                    for i in range(midpoint):
                        ## if the left nucleotide is reverse-complementary to the right nucleotide
                        if loop_seq[i] == complement_dict[loop_seq[-i - 1]]:
                            stem_len += 1
                            loop_len -= 2
                    # The loop should exist.
                    if loop_len >= 1:
                        stem_loop_loc.append([chrid, str(start), str(end), str(stem_len), str(loop_len), strand])
        os.remove(rnabobopt)
        stem_loop_loc = sorted(stem_loop_loc, key=lambda x: [x[0], int(x[1])])
        return stem_loop_loc

    def regularexpression_match(self, pattern, strand='+'):
        ## Use helitronscanner lcv file to detect terminal region of helitron
        coord_record = []
        for chrm in self.genome_dict:  ##  start coord is 1 not 0
            if strand == '+':
                genom_seq = str(self.genome_dict[chrm]).upper()
                pCT_start_list = re.finditer(pattern, genom_seq)
                for p_coord in pCT_start_list:
                    start = str(p_coord.start() + self.START + 1)
                    end = str(p_coord.end() + self.START)
                    coord_record.append([chrm, start, end])
            else:
                genom_seq = str(self.genome_dict[chrm].reverse_complement()).upper()
                sequence_length = len(genom_seq)
                pCT_start_list = re.finditer(pattern, genom_seq)
                for p_coord in pCT_start_list:
                    end = str(sequence_length - p_coord.start() + self.START)
                    start = str(sequence_length - p_coord.end() + self.START + 1)
                    coord_record.append([chrm, start, end])
        coord_record = sorted(coord_record, key=lambda x: [x[0], int(x[1])])
        return coord_record

    def inverted_detection(self, sequencefile, minitirlen, maxtirlen, mintirdist, maxtirdist, seed):
        ## start coord is 1 not 0
        dbname = ''.join([os.path.basename(sequencefile), '.invdb'])
        invttirfile = ''.join([os.path.basename(sequencefile), '.inv.txt'])
        ## build database
        mkinvdb = subprocess.Popen(
            ['gt', 'suffixerator', '-db', sequencefile, '-indexname', dbname, '-mirrored', '-dna', '-suf', '-lcp',
             '-bck'], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        mkinvdb.wait()
        ## run tirvish
        with open(invttirfile, 'w') as invf:
            runinvsearch = subprocess.Popen(
                ['gt', 'tirvish', '-index', dbname, '-mintirlen', str(minitirlen), '-maxtirlen', str(maxtirlen),
                 '-similar', str(Args.simtir), '-mintirdist', str(mintirdist), '-maxtirdist', str(maxtirdist), '-mintsd', '0',
                 '-seed', str(seed), '-vic', '1', '-overlaps', 'all', '-xdrop', '0'], stderr=subprocess.DEVNULL, stdout=invf)
            runinvsearch.wait()

        invt_list = []
        ## The default output is in gff format, extract coord information and do filteration.
        with open(invttirfile, 'r') as F:
            for line in F:
                if line.startswith('#'):
                    continue
                splitlines = line.rstrip().split('\t')
                if splitlines[2] == 'repeat_region':
                    chrmid = splitlines[0]
                    id = splitlines[8].replace('ID=', '')
                    t = 1
                elif splitlines[2] == 'terminal_inverted_repeat_element':
                    sim = re.findall('tir_similarity=(\d+\.\d+)', splitlines[8])[0]
                elif splitlines[2] == 'terminal_inverted_repeat':
                    if t == 1:
                        left_start = str(int(splitlines[3]) + self.START)
                        left_end = str(int(splitlines[4]) + self.START)
                        left_expand = '-'.join([left_start, left_end])
                        invt_length_left = int(splitlines[4]) - int(splitlines[3]) + 1
                        t += 1
                    else:
                        right_start = str(int(splitlines[3]) + self.START)
                        right_end = str(int(splitlines[4]) + self.START)
                        right_expand = '-'.join([right_start, right_end])
                        invt_length_right = int(splitlines[4]) - int(splitlines[3]) + 1
                        ## length of inverted sequences should be greater than 11 and shorter than 16
                        if invt_length_left >= 12 and invt_length_right >= 12 and invt_length_left <= 15 and invt_length_right <= 15:
                            invt_list.append([chrmid, str(left_start), str(right_end), left_expand, right_expand,
                                              (invt_length_right + invt_length_left) / 2, sim])

        invt_list = sorted(invt_list, key=lambda x: int(x[1]))
        os.remove(invttirfile)
        os.system('rm %s*' % dbname)
        return invt_list

# define Homologous_search class to find Helitron-like transposase domain and theri auto/non-auto relatives
class Homologous_search:
    def __init__(self, rep_hel_hmm, genome, wkdir, headerfile, window, distance_domain, pvalue, process_num):
        self.rep_hel_hmm = rep_hel_hmm
        self.genome = genome
        self.genome_dict = SeqIO.parse(genome, 'fasta')
        self.genome_dict = {k.id: k.seq.upper() for k in self.genome_dict}
        self.process_num = int(process_num)
        self.wkdir = wkdir
        self.headerpatternfile = headerfile
        self.window = window
        self.distance_domain = distance_domain
        self.pvalue = float(pvalue)
        # To transform the header file to regular expression.
        with open(headerfile, 'r') as F:
            headerpattern_list = F.read().rstrip().split('\n')
            self.headerpattern = '|'.join(headerpattern_list) if not Args.IS1 else '|'.join([''.join(['A', i]) for i in headerpattern_list])
        if not os.path.exists(self.wkdir):
            os.mkdir(self.wkdir)
            os.chdir(self.wkdir)
        else:
            os.chdir(self.wkdir)
        self.bedtoolstmp = os.path.abspath('BedtoolsTMP')
        if not os.path.exists(self.bedtoolstmp):
            os.mkdir(self.bedtoolstmp)
        BT.set_tempdir(self.bedtoolstmp)

        CWD = os.getcwd()
        self.genome_size = '%s/Genome.size' % CWD
        self.chrm_size = {i:len(self.genome_dict[i]) for i in self.genome_dict}
        genome_size = list(self.chrm_size.items())
        genome_size = sorted(genome_size, key=lambda x: x[0])
        ## To determine the evalue for short-sequence blastn, set the bit-score cutoff as 30, the evalue cutoff should follow the formula: m*n/(2**30)
        sum_genomesize = sum([i[1] for i in genome_size])
        self.evalue_blastn = sum_genomesize * 30/(2**30)
        with open(self.genome_size, 'w') as F:
            F.writelines([''.join([i[0], '\t', str(i[1]), '\n']) for i in genome_size])
        
        # To define stem_loop structure ending with CTRR motif of Helitron 
        self.CTRR_stem_loop_description = '%s/CTRR_stem_loop.descr' % CWD
        CTRR_description = """r1 s1 r1' s2\nr1 1:1 NNNNN[10]:[10]NNNNN TGCA\ns1 0 N[7]\ns2 0 N[15]CTRR%s\n"""
        # Add 'T' in the end if user limitted the 'A-T' insertion site for Helitron.
        CTRR_description = CTRR_description % 'T' if Args.IS1 else CTRR_description % ''
        with open(self.CTRR_stem_loop_description, 'w') as F:
            F.write(CTRR_description)

        # To define stem_loop structure of  Helentron/Helitron2
        self.subtir_description = '%s/subtir_stem_loop.descr' % CWD
        # Add 'T' in the end if user limitted the 'T-T' insertion site for Helentron/Helitron2.
        if Args.IS2:
            subtir_description = """r1 s1 r1' s2\nr1 1:1 NNNNN[10]:[10]NNNNN TGCA\ns1 0 N[15]\ns2 0 NNNNN[10]T\n"""
        else:
            subtir_description = """r1 s1 r1' s2\nr1 1:1 NNNNN[10]:[10]NNNNN TGCA\ns1 0 N[15]\ns2 0 NNNNNN[2]\n"""
        with open(self.subtir_description, 'w') as F:
            F.write(subtir_description)

        dbdir = 'GenomeDB/'
        if not os.path.exists(dbdir):
            os.mkdir(dbdir)
        self.genomedb = ''.join([CWD, '/', dbdir, os.path.basename(self.genome), '.blastndb'])
        makeblastndb = subprocess.Popen(['makeblastdb', '-dbtype', 'nucl', '-in', self.genome, '-out', self.genomedb],
                                        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        makeblastndb.wait()

    def hmmsearch(self, subgenome):
        # Run hmmersearch program to search for Helitron-like transposase
        orf_file = ''.join([subgenome, '.orf'])
        hmm_opt = ''.join([subgenome, '.hmmsearch.out'])

        #The index of getorf output starts from 1, not 0
        # Use getorf to predicte open reading frames for a given genome
        get_orf = subprocess.Popen(['getorf', '-sequence', subgenome, '-outseq', orf_file, '-minsize', '100'],
                                   stderr=subprocess.DEVNULL)
        get_orf.wait()

        Rep_dict, Hel_dict = defaultdict(list), defaultdict(list)
        Rep_opline, Hel_opline = [], []
        if os.path.getsize(orf_file):
            run_hmmsearch = subprocess.Popen(
                ['hmmsearch', '--domtblout', hmm_opt, '--noali', '-E', '1e-3', self.rep_hel_hmm, orf_file],
                stdout=subprocess.DEVNULL)
            run_hmmsearch.wait()
            os.remove(orf_file)
        else:
            return Rep_opline, Hel_opline
        if not os.path.exists(hmm_opt):
            return Rep_opline, Hel_opline

        # To parser hmmsearch output
        with open(hmm_opt, 'r') as F:
            for line in F:
                if line.startswith('#'):
                    continue
                splitlines = re.split('\s+', line.rstrip())
                domain, sub_class = splitlines[3].split('_')
                subchrname = "_".join(splitlines[0].split('_')[:-1])
                chrm_name, START = subchrname.split('startat')
                start, end = re.findall('\[(\d+)\s+-\s+(\d+)\]', line)[0]
                start = str(int(start) + int(START))
                end = str(int(end) + int(START))

                aa_start, aa_end = splitlines[19:21]
                score = splitlines[7]
                c_evalue, i_evalue = splitlines[11:13]
                if float(c_evalue) > 1e-5 or float(i_evalue) > 1e-5:
                    continue
                ## To transform amino acide coord to nucleotide coord
                if int(end) > int(start):
                    strand = '+'
                    nuc_start = int(aa_start) * 3 - 3 + int(start)
                    nuc_end = int(aa_end) * 3 + int(start) - 1
                    orf_loc = '-'.join([start, end])
                else:
                    nuc_end = int(start) - 3 * int(aa_start) + 3
                    nuc_start = int(start) - 3 * int(aa_end) + 1
                    strand = '-'
                    orf_loc = '-'.join([end, start])
                if domain.startswith('Hel'):
                    Hel_dict[splitlines[0]].append(
                        [chrm_name, str(nuc_start), str(nuc_end), sub_class, score, strand, orf_loc])
                else:
                    Rep_dict[splitlines[0]].append(
                        [chrm_name, str(nuc_start), str(nuc_end), sub_class, score, strand, orf_loc])
        for key in Hel_dict:
            hel_candidate = sorted(Hel_dict[key], key=lambda x: float(x[4]))[-1]  ## Select the case with highest score.
            Hel_opline.append(hel_candidate)
        for key in Rep_dict:
            rep_candidate = sorted(Rep_dict[key], key=lambda x: float(x[4]))[-1]  ## Select the case with highest score.
            Rep_opline.append(rep_candidate)
        try:
            Hel_bed = BT.BedTool([BT.create_interval_from_list(line) for line in Hel_opline]).sort()
            Rep_bed = BT.BedTool([BT.create_interval_from_list(line) for line in Rep_opline]).sort()
            return Rep_bed, Hel_bed
        except:
            return [], []

    def intersect(self, location1, location2, slip=0, lportion=0.0, rportion=0.0, bool_and=1):
        # Define intersect function to check if two intervals are intersected or not, similar to bedtools intersect
        location1 = sorted([int(i) for i in location1])
        location2 = sorted([int(i) for i in location2])
        if location1[0] - slip > location2[1] or location1[1] < location2[0] - slip:
            return False
        else:
            total_list = sorted([location1[0], location1[1], location2[0], location2[1]])
            portion1 = (total_list[2] - total_list[1] + 1) / (
                    location1[1] - location1[0] + 1)  ## how much proportion the intersected sequence occupiedonseq1
            portion2 = (total_list[2] - total_list[1] + 1) / (
                    location2[1] - location2[0] + 1)  ## how much proportion the intersected sequence occupiedonseq2
            if bool_and:
                if portion1 >= lportion and portion2 >= rportion:
                    return True
                else:
                    return False
            else:
                if portion1 >= lportion or portion2 >= rportion:
                    return True
                else:
                    return False

    def merge_bedfile(self, BedInput, orf=True, window=1500, splitid=False, S=True):
        # Define function to merge two distance-close genomic features
        if S:
            cluster_bed = BedInput.cluster(d=window, s=True)
        else:
            cluster_bed = BedInput.cluster(d=window)

        merge_dict = defaultdict(list)
        merge_list = []
        for line in cluster_bed:
            line = list(line)
            cluster = line[-1]
            merge_dict[cluster].append(line[:-1])

        for cluster in merge_dict:
            ## To record the domain location.
            coordlist = merge_dict[cluster]
            coord_set = [int(i[1]) for i in coordlist]
            coord_set.extend([int(i[2]) for i in coordlist])
            coord_set = sorted(coord_set)
            start = coord_set[0]
            stop = coord_set[-1]
            strand = coordlist[0][5]
            if orf:
                ## To determain sub_class, select the case with highest score
                sub_class = sorted(coordlist, key=lambda x: float(x[4]))[-1][3]
                ## To merge the ORF location
                orf_list = sorted([int(b) for i in coordlist for b in i[-1].split('-')])
                orf_coord = '-'.join([str(orf_list[0]), str(orf_list[-1])])
            else:
                sub_class = '.'
                bitscore = [float(i[4]) for i in coordlist]
                orf_coord = str(sum(bitscore) / len(bitscore))
            if splitid:
                chrm_id, sub_class = coordlist[0][0].split('#')
            else:
                chrm_id = coordlist[0][0]
            merge_list.append([chrm_id, str(start), str(stop), sub_class, orf_coord, strand])
        # merge_list = sorted(merge_list, key=lambda x: [x[0], int(x[1])])
        if merge_list:
            merge_bed = BT.BedTool([BT.create_interval_from_list(line) for line in merge_list]).sort()
            return merge_bed
        else:
            return 0

    def parser_hmmsearch(self, Rep_bed, Hel_bed, subgenome):
        # To find Rep-Hel structure which might imply a possible Helitron-like transposases.
        if not Rep_bed or not Hel_bed:
            return []
        
        ## To merge helicase or rep domain splicing sites (helitron-like transposase contain introns)
        merge_hel = self.merge_bedfile(Hel_bed, orf=True, window=1500)
        merge_rep = self.merge_bedfile(Rep_bed, orf=True, window=1500)
        if not merge_hel or not merge_rep:  ## Either hel or rep data is null
            return []

        ## To find rep and helicase gene pairs that rep is less than self.distance_domain bp upstream of hel
        joint_rephel = merge_rep.window(merge_hel, l=0, r=int(self.distance_domain), sm=True, sw=True)
        bedlist = []
        for line in joint_rephel:
            splitlines = list(line)
            strand = splitlines[5]
            chrm = splitlines[0]
            rep_start, rep_end = splitlines[1:3]
            hel_start, hel_end = splitlines[7:9]

            # If the helicase and rep domain have a intersection, skip
            if self.intersect([int(rep_start), int(rep_end)], [int(hel_start), int(hel_end)], lportion=0.2, rportion=0.2):
                continue
            rep_orf = splitlines[4].split('-')
            hel_orf = splitlines[10].split('-')
            loc = sorted([rep_start, rep_end, hel_start, hel_end], key=lambda x: int(x))
            start, end = loc[0], loc[-1]  ## They are REP and Helicase domain region
            bedlist.append((chrm, int(start), int(end), '-'.join([rep_start, rep_end]), '-'.join([hel_start, hel_end]),
                            strand, splitlines[3], splitlines[9], 'NA'))
        bedlist = list(set(bedlist))  ## To avoid duplicates
        bedlist = sorted(bedlist, key=lambda x: [x[0], x[1]])
        return bedlist

    def heltentron_terminal(self, helentron_bed):
        # Define function to try to recover Helentron terminal region (stem-loop structure) which is behind the right part of TIRs
        opbed_list = []
        extend_seq = []
        extend_file = 'helentron.extend.fa'
        extend_dict = {}
        extend_dict = defaultdict(list)

        ## To output all extend rigions into one single file
        #init_num = 1
        for feature in helentron_bed:
            feature = list(feature)
            chrmid, start, stop, name, score, strand, pvalue, Bscore = feature
            start, stop = int(start), int(stop)
            if float(Bscore) == 0:  ## without terminal signals
                opbed_list.append([chrmid, start, stop, name, score, strand, pvalue, Bscore])
                continue

            if strand == '+':
                extend_id = '-'.join(['seq', str(start), str(stop), 'p'])
                detect_seq = str(self.genome_dict[chrmid][stop: stop + 80])
                extend_dict[extend_id].append([chrmid, start, stop, name, score, strand, pvalue, Bscore,
                                          stop])  ## The last element is the initial start for terminal detection region
                extend_seq.append(''.join(['>', extend_id, '\n', detect_seq, '\n']))
            else:
                extend_id = '-'.join(['seq', str(start), str(stop), 'n'])
                terminal_start = 0 if start - 80 < 0 else start - 80
                extend_dict[extend_id].append([chrmid, start, stop, name, score, strand, pvalue, Bscore,
                                          terminal_start])  ## The last element is the initial start for terminal detection region
                detect_seq = str(self.genome_dict[chrmid][terminal_start: start])
                extend_seq.append(''.join(['>', extend_id, '\n', detect_seq, '\n']))
        if not extend_seq:  ## means empty
            return helentron_bed
        with open(extend_file, 'w') as F:
            F.writelines(extend_seq)
        ## stem loop detection
        stem_loop_list = Structure_search(extend_file, START=0).stem_loop(self.subtir_description, minus_tailone=int(Args.IS2))
        stem_loop_dict = defaultdict(list)
        [stem_loop_dict[i[0]].append(i) for i in stem_loop_list]

        ## To select the nearest candidate
        for extend_id in stem_loop_dict:
            for sublist in extend_dict[extend_id]:
                chrmid, start, stop, name, score, strand, pvalue, Bscore, e_start = sublist
                if extend_id.endswith('p'):
                    stem_loop = [i for i in stem_loop_dict[extend_id] if i[5] == '+']
                    ## order by start position (closer), stem length (longer), loop length (shorter),  total length (shorter)
                    stem_loop = sorted(stem_loop, key = lambda x: [x[0], int(x[1]), -int(x[3]), int(x[4]), int(x[2]) - int(x[1])])
                    if not stem_loop:
                        opbed_list.append([chrmid, start, stop, name, score, strand, pvalue, Bscore])
                        continue
                    stem_start, stem_stop = stem_loop[0][1:3]
                    stem_stop = int(stem_stop) + int(e_start)
                    opbed_list.append([chrmid, start, stem_stop, name, score, strand, pvalue, Bscore])
                else:
                    stem_loop = [i for i in stem_loop_dict[extend_id] if i[5] == '-']
                    ## order by end position (longer), stem length (longer), loop length (shorter), total length (shorter)
                    stem_loop = sorted(stem_loop, key=lambda x: [x[0], -int(x[2]), -int(x[3]), int(x[4]), int(x[2]) - int(x[1])])
                    if not stem_loop:
                        opbed_list.append([chrmid, start, stop, name, score, strand, pvalue, Bscore])
                        continue
                    stem_start, stem_stop = stem_loop[0][1:3]
                    stem_start = int(stem_start) + int(e_start)
                    opbed_list.append([chrmid, stem_start, stop, name, score, strand, pvalue, Bscore])

        ### To complement the candidates whose stem loop signal doesn't exist.
        remained_cases = set(extend_dict.keys()) - set(stem_loop_dict.keys())
        for key in remained_cases:
            for sublist in extend_dict[key]:
                opbed_list.append(sublist[:8])
        ### To creat bedtools objective
        opbed = BT.BedTool([BT.create_interval_from_list(line) for line in opbed_list])
        os.remove(extend_file)
        return opbed

    def intergrated_program(self, subgenome):
        # This function is used to recover terminal signals of Helitron-like elements (TIRs for Helentron/Helitron2; TC... motif and ...CTRR motif for Helitron)
        rep_hmmsearch_opt, hel_hmmsearch_opt = self.hmmsearch(subgenome)
        ORF_list = self.parser_hmmsearch(rep_hmmsearch_opt, hel_hmmsearch_opt, subgenome)
        sys.stdout.write('Find %s rep-hel blocks in %s.\n' % (str(len(ORF_list)), os.path.basename(subgenome).replace('.fa', '')))
        RC_total_candidate = []
        for Helitron_candidate in ORF_list:
            ORF_chrmid = Helitron_candidate[0]
            ORF_start = int(Helitron_candidate[1])
            ORF_stop = int(Helitron_candidate[2])
            rep_loc = Helitron_candidate[3]
            hel_loc = Helitron_candidate[4]
            strand = Helitron_candidate[5]
            chrm_limit = len(self.genome_dict[ORF_chrmid])
            rep_name, hel_name = Helitron_candidate[6:8]

            ## To produce orf id which will be used as sole identifier of terminal signals.
            ORFID = '-'.join([ORF_chrmid, str(ORF_start), str(ORF_stop)])
            temp_name_for_helitron = '-'.join([ORF_chrmid, str(ORF_start), str(ORF_stop)])
            expansion_all_seqname = ''.join([temp_name_for_helitron, '.expansion.fa'])
            left_seqname = ''.join([temp_name_for_helitron, '.left.fa'])
            right_seqname = ''.join([temp_name_for_helitron, '.right.fa'])

            # Define regions to search for terminal signals
            left_boundary = ORF_start - self.window
            right_boundary = ORF_stop + self.window

            ## left and right boundary should be within chromosome ranges
            if left_boundary <= 0:
                left_boundary = 1
            if right_boundary >= chrm_limit:
                right_boundary = chrm_limit

            ## To avoid get big tandem heltron-like elements
            Helitron_candidate_index = ORF_list.index(Helitron_candidate)
            # Left boundary should not touch last Rep-Hel region
            if Helitron_candidate_index >= 1:  ## not the fist one
                last_candidate = ORF_list[Helitron_candidate_index - 1]
                last_stop = last_candidate[2]
                left_boundary = last_stop if left_boundary < last_stop else left_boundary
            # Right boundary should not touch the next Rep-Hel region
            if Helitron_candidate_index < len(ORF_list) - 1:  ## not the last one
                next_candidate = ORF_list[Helitron_candidate_index + 1]
                next_start = next_candidate[1]
                right_boundary = next_start if right_boundary > next_start else right_boundary
           
            # To output extended sequences to fasta files
            left_seq = str(self.genome_dict[ORF_chrmid][left_boundary - 1: ORF_start])
            right_seq = str(self.genome_dict[ORF_chrmid][ORF_stop: right_boundary])
            expansion_all_seq = str(self.genome_dict[ORF_chrmid][left_boundary - 1: right_boundary])
            with open(expansion_all_seqname, 'w') as F:
                F.write(''.join(['>', ORF_chrmid, '\n', expansion_all_seq, '\n']))
            with open(left_seqname, 'w') as F:
                F.write(''.join(['>', ORF_chrmid, '\n', left_seq, '\n']))
            with open(right_seqname, 'w') as F:
                F.write(''.join(['>', ORF_chrmid, '\n', right_seq, '\n']))
            
            # Candiadte is Helitron, try to search for left terminal signals in left extension and right terminal signals in right extension.
            if hel_name == 'Helitron' and rep_name == 'Helitron':
                if strand == '+':
                    # left terminal signals are TC... like
                    TC_list = Structure_search(genome=left_seqname, START=left_boundary - 1).regularexpression_match(self.headerpattern, '+')
                    # right terminal signals are stem-loop structures ending with CTRR motif.
                    Stem_loop_list = Structure_search(genome=right_seqname, START=ORF_stop).stem_loop(
                        self.CTRR_stem_loop_description, minus_tailone=int(Args.IS1))
                    Stem_loop_list = [i for i in Stem_loop_list if i[-1] == '+']  ## To select the positive strand motif
                    # To reduce one if user set 'AT' insertion because the A was added at the begaining of header regular expression
                    if Args.IS1:
                        TC_list = [[line[0], str(int(line[1])+1), line[2]] for line in TC_list]
                else:
                    TC_list = Structure_search(genome=right_seqname, START=ORF_stop).regularexpression_match(
                        self.headerpattern, '-')
                    Stem_loop_list = Structure_search(genome=left_seqname, START=left_boundary - 1).stem_loop(
                        self.CTRR_stem_loop_description, minus_tailone=int(Args.IS1))
                    Stem_loop_list = [i for i in Stem_loop_list if i[-1] == '-']  ## To select the negative strand motif
                    if Args.IS1:
                        TC_list = [[line[0], line[1], str(int(line[2]) - 1)] for line in TC_list]
                RC_total_candidate.append((
                                          ORF_chrmid, str(ORF_start), str(ORF_stop), strand, 'Helitron', tuple(TC_list),
                                          tuple(Stem_loop_list), ORFID))

            # Candidate is Helentron or Helitron2
            else:
                if hel_name in ['Helentron', 'Helitron2'] and rep_name in ['Helentron', 'Helitron2']:
                    class_name = rep_name if rep_name == hel_name else '_or_'.join([rep_name, hel_name])
                else:
                    # unable to distinguish, will take it as Helentron-like
                    class_name = '_or_'.join([rep_name, hel_name])

                # Size of terminal inverted sequences should be greater than the size of predicted transposase
                mini_dist_tir = ORF_stop - ORF_start
                # Size of terminal inverted sequences should be less than the size of extension
                max_dist_tir = 2 * int(self.window) + ORF_stop - ORF_start

                ## The expected length of helentron/helitron2 should be from 12 to 15, set the maximum length to 20 and exclude the long tir (>15 nt) after.
                invt_list = Structure_search(genome=expansion_all_seqname, START=left_boundary - 1).inverted_detection(
                    expansion_all_seqname, 9, 20, mini_dist_tir, max_dist_tir, 8)
                ## To keep the case that fully covered with ORF region
                invt_list = [i for i in invt_list if int(i[3].split('-')[1]) - ORF_start <= 0 and int(i[4].split('-')[0]) - ORF_stop >= 0]
                left_list, right_list = [], []
                for inv in invt_list:
                    if strand == '+':
                        left_loc = inv[3].split('-')
                        right_loc = inv[4].split('-')
                        if Args.IS2:
                            is_start = int(left_loc[0]) - 2   ## To get the insertion site index. The real index starts from 1, need to transform to python index.
                            is_seq = str(self.genome_dict[ORF_chrmid][is_start:is_start+1]).upper()
                            if is_seq == 'T':
                                left_list.append((ORF_chrmid, int(left_loc[0]), int(left_loc[1])))
                                right_list.append((ORF_chrmid, int(right_loc[0]), int(right_loc[1])))
                        else:
                            left_list.append((ORF_chrmid, int(left_loc[0]), int(left_loc[1])))
                            right_list.append((ORF_chrmid, int(right_loc[0]), int(right_loc[1])))
                    else:
                        left_loc = inv[4].split('-')
                        right_loc = inv[3].split('-')
                        if Args.IS2:
                            is_start = int(left_loc[1])
                            is_seq = str(self.genome_dict[ORF_chrmid][is_start:is_start + 1]).upper()
                            if is_seq == 'A':   ## The complementary of T is A
                                left_list.append((ORF_chrmid, int(left_loc[0]), int(left_loc[1])))
                                right_list.append((ORF_chrmid, int(right_loc[0]), int(right_loc[1])))
                        else:
                            left_list.append((ORF_chrmid, int(left_loc[0]), int(left_loc[1])))
                            right_list.append((ORF_chrmid, int(right_loc[0]), int(right_loc[1])))
                RC_total_candidate.append((ORF_chrmid, str(ORF_start), str(ORF_stop), strand, class_name,
                                           tuple(left_list), tuple(right_list), ORFID))
            os.remove(expansion_all_seqname)
            os.remove(left_seqname)
            os.remove(right_seqname)
        return RC_total_candidate

    def cdhitest_clust(self, input_fa, opfile, helitron_type, id=0.8):
        # This function is for clustering of highly identity sequences
        cons_name = '.'.join([helitron_type, 'reduce.temp'])
        cluster_file = '.'.join([cons_name, 'clstr'])
        run_cluster = subprocess.Popen(
            ['cd-hit-est', '-i', input_fa, '-o', cons_name, '-d', '0', '-aS', '0.8', '-c', str(id), '-G', '1', '-g',
             '1', '-b', '500', '-T', str(self.process_num), '-M', '0'], stdout=subprocess.DEVNULL)
        run_cluster.wait()
        
        # To get classification information
        cluster_dict = {}
        with open(cluster_file, 'r') as F:
            for line in F:
                if line.startswith('>'):
                    cluster_name = '_'.join([helitron_type, line.strip('>\n').split(' ')[1]])
                else:
                    insertion_name = line.split('...')[0].split(', >')[1]
                    cluster_dict[insertion_name] = cluster_name

        opseq = ''
        with open(cons_name, 'r') as F:
            for line in F:
                if line.startswith(">"):
                    insertion_name = line.strip('>\n')
                    cluster_name = cluster_dict[insertion_name]
                    opseq += ''.join(['>', cluster_name, '\n'])
                else:
                    opseq += line
        with open(opfile, 'w') as F:
            F.write(opseq)
        if os.path.exists(cons_name):
            os.remove(cons_name)
        os.remove(cluster_file)
        cluster_file = '.'.join([helitron_type, 'clust.info2'])
        with open(cluster_file, 'w') as F:
            F.writelines(["".join([k, "\t", cluster_dict[k], '\n']) for k in cluster_dict])
        return cluster_dict

    def split_list(self, numbers, num_groups):
        # Calculate target sum for each group
        total_sum = sum([i[1] for i in numbers])
        target_sum = total_sum / num_groups

        # Sort numbers in descending order
        numbers = sorted(numbers, key=lambda x: -x[1])
        # Split numbers into groups with similar sums
        groups = [[] for i in range(num_groups)]
        group_sums = [0] * num_groups
        for number in numbers:
            # Find the group with the smallest current sum and add the number to it
            min_sum_index = group_sums.index(min(group_sums))
            groups[min_sum_index].append(number)
            group_sums[min_sum_index] += number[1]
        return groups

    def split_genome(self, chunk_size=200000000, flanking_size=50000, num_groups=2):
        # To split big genomes into small chunks
        if not os.path.exists('genomes'):
            os.mkdir('genomes')
        subgenome_list = []
        for chrm in self.genome_dict:
            seq_len = len(self.genome_dict[chrm])
            if seq_len < 1000:  ##Skip chrms whose length is shorter than 1000 bp
                sys.stdout.write(
                    "Chrm %s will not be used to detect autonomous Helitron/Helentron as its length is shorter than 1000 bp\n" % chrm)
                continue
            subgenome_list.append((chrm, seq_len))

        num_groups = num_groups if num_groups <= len(subgenome_list) else len(subgenome_list)
        ###  To split the genomes into several files.
        # Calculate target sum for each group
        total_sum = sum([i[1] for i in subgenome_list])
        target_sum = total_sum / num_groups
        # Sort numbers in descending order
        numbers = sorted(subgenome_list, key=lambda x: -x[1])
        # Split numbers into groups with similar sums
        groups = [[] for i in range(num_groups)]
        group_sums = [0] * num_groups
        for number in numbers:
            # Find the group with the smallest current sum and add the number to it
            min_sum_index = group_sums.index(min(group_sums))
            groups[min_sum_index].append(number)
            group_sums[min_sum_index] += number[1]

        ## To split big chrms into smaller chunks.
        subgenome_list = []
        init_num = 1
        for subgroup in groups:
            subgenome = ''.join(['genomes/subgenome', str(init_num), '.fa'])
            init_num += 1
            with open(subgenome, 'w') as F:
                for chrminfo in subgroup:
                    chrid = chrminfo[0]
                    seq = self.genome_dict[chrid]
                    seq_len = len(seq)
                    num = seq_len // chunk_size
                    for i in range(num + 1):
                        start, stop = i * chunk_size, (i + 1) * chunk_size + flanking_size
                        if start >= seq_len:
                            continue
                        if stop > seq_len:
                            stop = seq_len
                        subchrm = 'startat'.join([chrid, str(start)])
                        chunk_seq = str(self.genome_dict[chrid][start:stop])
                        F.write(''.join(['>', subchrm, '\n']))
                        F.write(chunk_seq)
                        F.write('\n')
            subgenome_list.append(subgenome)
        return subgenome_list

    def autonomous_detect(self):
        # main program to search for transposae and terminal signals.
        subgenome_list = self.split_genome(chunk_size=200000000, flanking_size=20000, num_groups=200)
        if len(subgenome_list) < self.process_num:
            processnum = len(subgenome_list)
        else:
            processnum = self.process_num
        
        # Use python multiple threading
        planpool = ThreadPool(processnum)
        run_result = []
        for subgenome in subgenome_list:
            run_result.append(planpool.apply_async(self.intergrated_program, args=(subgenome,)))
        planpool.close()
        planpool.join()
        Helitron_list = []
        for result in run_result:
            result_get = result.get()
            if result_get:
                Helitron_list.extend(result_get)
        Helitron_list = sorted(Helitron_list, key=lambda x: [x[0], int(x[1]), int(x[2])])
        return Helitron_list

    def blastn(self, query_file, merge=True):
        # To search for homologous of given sequences
        blastn_opt = ''.join([query_file, '.tbl'])
        blastn_bed = ''.join([query_file, '.bed6'])
        mergeblastn_bed = ''.join([query_file, '.merge.bed6'])
        blastn_pro = subprocess.Popen(
            ['blastn', '-db', self.genomedb, '-query', query_file, '-num_threads', str(self.process_num),
             '-max_target_seqs', '999999999', '-evalue', str(self.evalue_blastn), '-task', 'blastn-short',
             '-outfmt', '6 qseqid sseqid pident qstart qend sstart send evalue qlen bitscore', '-out', blastn_opt])
        blastn_pro.wait()
        oplines = []
        with open(blastn_opt, 'r') as F:
            for line in F:
                splitlines = line.rstrip().split('\t')
                query_name = splitlines[0]
                chrm, identity, qstart, qend, sstart, send, evalue, qlen, bitscore = splitlines[1:10]
                coverage = round((abs(int(qstart) - int(qend)) + 1) / abs(int(qlen)), 2)
                if int(sstart) < int(send):
                    START = sstart
                    END = send
                    strand = '+'
                else:
                    START = send
                    END = sstart
                    strand = '-'
                if merge:
                    oplines.append(['#'.join([chrm, query_name]), START, END, query_name, bitscore, strand])
                else:
                    oplines.append([chrm, START, END, query_name, bitscore, strand])
        blastn_bed = BT.BedTool([BT.create_interval_from_list(line) for line in oplines]).sort()
        if merge:
            mergeblastn_bed = self.merge_bedfile(blastn_bed, orf=False, window=100, splitid=True, S=True)  ## To merge the tandem repeats
            return mergeblastn_bed
        else:
            return blastn_bed

    def merge_overlaped_intervals(self, bedinput, classname, mobile_type):
        # To merge overlaped candidates
        bedinput = sorted([list(line) for line in bedinput],
                          key=lambda x: [x[0], x[5], int(x[1])])  ## sort by chrm, strand, and coord
        if not bedinput:
            return []
        compared_line = (bedinput[0][0], bedinput[0][1], bedinput[0][2], bedinput[0][5])

        blockname_init = 1
        blockname_dict = defaultdict(list)
        alternative_bedlines = []
        for line in bedinput:
            chrmid, start, stop, pairname, count, strand, pvalue, Bscore = line
            newline = [chrmid, start, stop, pairname, count, strand, pvalue, Bscore, classname, mobile_type]
            if compared_line[0] == chrmid and self.intersect(compared_line[1:3], line[1:3], lportion=0.8, rportion=0.8, bool_and=0) and compared_line[3] == strand:
                merged_coord = sorted([compared_line[1], compared_line[2], start, stop], key=lambda x: int(x))
                compared_line_s, compared_line_e = merged_coord[0], merged_coord[-1]
                #compared_line = [chrmid, compared_line_s, compared_line_e, strand, mobile_type]
                compared_line = [chrmid, start, compared_line_e, strand, mobile_type]
                alternative_bedlines.append(newline)
            else:
                block_name = '_'.join(['insertion', classname, mobile_type, str(blockname_init)])
                blockname_dict[block_name] = alternative_bedlines
                blockname_init += 1
                alternative_bedlines = [newline]
                compared_line = [chrmid, start, stop, strand, mobile_type]

        ## in case the last overlaped series not saved.
        block_name = '_'.join(['insertion', classname, mobile_type, str(blockname_init)])
        blockname_dict[block_name] = alternative_bedlines

        alternative_list = []
        for blockname in blockname_dict:
            ## To select the significant candidates: sort by pvalue, length, count
            candidate_list = blockname_dict[blockname]
            filtered_candidate_list = [line for line in candidate_list if int(line[4]) >=2]
            candidate_list = filtered_candidate_list if filtered_candidate_list else candidate_list
            [i.append(blockname) for i in candidate_list]
            alternative_list.extend(candidate_list)
        # chrmid, start, stop, pairname, count, strand, pvalue, classname, mobile_type, blockname
        return alternative_list

    def split_joint(self, left_bed, right_bed, combind_id_list, sub_bed_dir, half_distance):
        # To split big bedfiles 
        leftflank_dir = ''.join([sub_bed_dir, '/', 'leftflank'])
        rightflank_dir = ''.join([sub_bed_dir, '/', 'rightflank'])
        jointflank_file = ''.join([sub_bed_dir, '/', 'left_right.path.join'])

        if os.path.exists(sub_bed_dir):
            shutil.rmtree(sub_bed_dir)
        os.mkdir(sub_bed_dir)
        os.mkdir(leftflank_dir)
        os.mkdir(rightflank_dir)

        left_dict = defaultdict(list)
        right_dict = defaultdict(list)

        for line in left_bed:
            left_dict[line[3]].append(line)
        for line in right_bed:
            right_dict[line[3]].append(line)

        ## To find uniq left and right terminal markers
        left_list, right_list = [], []
        for combinied in combind_id_list:
            left_list.append(combinied[0])
            right_list.append(combinied[1])
        left_list = list(set(left_list))
        right_list = list(set(right_list))
        ## To output sub-bed files
        for left_name in left_list:
            leftflank_file = ''.join([leftflank_dir, '/', left_name, '.bed'])
            subbed = BT.BedTool(left_dict[left_name]).sort()
            subbed.flank(g=self.genome_size, l=0, r=half_distance, s=True).sort().saveas(leftflank_file)
        for right_name in right_list:
            rightflank_file = ''.join([rightflank_dir, '/', right_name, '.bed'])
            subbed = BT.BedTool(right_dict[right_name]).sort()
            subbed.flank(g=self.genome_size, l=half_distance, r=0, s=True).sort().saveas(rightflank_file)
        ## To output join file
        file_list = ["".join([leftflank_dir, '/', combin[0], '.bed', '\t', rightflank_dir, '/', combin[1], '.bed', '\n']) for
                     combin in combind_id_list]
        with open(jointflank_file, 'w') as F:
            F.writelines(file_list)
        return jointflank_file

    def prepare_terminal_seq(self, Helitron_list, pair=False, classname='Helitron'):
        if not os.path.exists(classname):
            os.mkdir(classname)
        os.chdir(classname)
        
        left_exist, right_exist = 0, 0

        ORF_bedfile = ''.join([classname, '_orf.bed'])
        left_ter_file = ''.join([classname, '_left.fa'])  ## the file name will be split by '.' in vsearch function.
        left_ter_reduce_file = ''.join([classname, '.reduce.left.fa'])

        right_ter_file = ''.join([classname, '_right.fa'])
        right_ter_reduce_file = ''.join([classname, '.reduce.right.fa'])

        ORF_length_list = []

        ## To build container for both left and right pair
        EXTEND=30
        Helitron_pair_list = []
        with open(ORF_bedfile, 'w') as orfF, open(left_ter_file, 'w') as leftF, open(right_ter_file, 'w') as rightF:
            for line in Helitron_list:
                left_ter_name, right_ter_name = [], []
                strand = line[3]
                orfF.write(''.join([line[0], '\t', line[1], '\t', line[2], '\t', line[7], '\t', line[4], '\t', strand, '\n']))
                ORF_length_list.append(int(line[2]) - int(line[1]) + 1)
                left_list = line[5]
                left_seq = []
                init = 0
                for left in left_list:
                    left_exist += 1
                    id = '.'.join([line[7], 'left', str(init)])
                    init += 1
                    left_ter_name.append(id)
                    if strand == '+':
                        loc_s = int(left[1]) - 1
                        loc_s = loc_s if loc_s >= 0 else 0
                        seq = str(self.genome_dict[line[0]][loc_s:int(left[1]) + EXTEND-1].upper())
                    else:
                        seq = self.genome_dict[line[0]][int(left[2]) - EXTEND:int(left[2])]
                        seq = str(seq.reverse_complement().upper())
                    left_seq.append(''.join(['>', id, '\n', seq, '\n']))

                right_list = line[6]
                right_seq = []
                init = 0
                for right in right_list:
                    right_exist += 1
                    id = '.'.join([line[7], 'right', str(init)])
                    init += 1
                    right_ter_name.append(id)
                    if strand == '-':
                        loc_s = int(right[1]) - 1
                        loc_s = loc_s if loc_s >= 0 else 0
                        seq = self.genome_dict[line[0]][loc_s:int(right[1]) + EXTEND-1]
                        seq = str(seq.reverse_complement().upper())
                    else:
                        seq = str(self.genome_dict[line[0]][int(right[2]) - EXTEND:int(right[2])].upper())
                    right_seq.append(''.join(['>', id, '\n', seq, '\n']))
                    
                if classname == 'Helitron' and pair:  ## pair the left and right terminals that ever appears on the same Helitron region.
                    [Helitron_pair_list.append((left, tuple(right_ter_name))) for left in left_ter_name]  ## all possibilities.
                elif classname != 'Helitron':
                    [Helitron_pair_list.append((left_ter_name[i], right_ter_name[i])) for i in range(len(left_ter_name))]
                leftF.writelines(left_seq)
                rightF.writelines(right_seq)

        if not left_exist or not right_exist:  ## means that either left or right terminal signals does not exist, so skip blastn
            os.chdir('../')
            return [], []

        # To redunce redundency via cd-hit-est
        left_cluster_dict = self.cdhitest_clust(left_ter_file, left_ter_reduce_file, helitron_type=''.join([classname, '_left']), id=0.9)
        right_cluster_dict = self.cdhitest_clust(right_ter_file, right_ter_reduce_file, helitron_type=''.join([classname, '_right']), id=0.9)

        ### To get the collapsed cluster pairs
        if pair:
            collapsed_pair_dict = defaultdict(list)
            terminal_pair_dict = dict(Helitron_pair_list)
            if classname == 'Helitron':  ## pair the left and right terminals that ever appears on the same Helitron region.
                for left_name in terminal_pair_dict:
                    right_name_list = terminal_pair_dict[left_name]
                    for right_name in right_name_list:
                        collapsed_left_clust = left_cluster_dict[left_name]
                        collapsed_right_clust = right_cluster_dict[right_name]
                        collapsed_pair_dict[collapsed_left_clust].append(collapsed_right_clust)
                for key in collapsed_pair_dict:
                    collapsed_pair_dict[key]=list(set(collapsed_pair_dict[key]))
            else:  ## pair inverted repeats
                for left_name in terminal_pair_dict:
                    right_name = terminal_pair_dict[left_name]
                    collapsed_left_clust = left_cluster_dict[left_name]
                    collapsed_right_clust = right_cluster_dict[right_name]
                    collapsed_pair_dict[collapsed_left_clust].append(collapsed_right_clust)
        ## To creat bedtool boject
        ORF_bed = BT.BedTool(ORF_bedfile)
        left_bedfile_name, right_bedfile_name = '%s.left.bed' % classname, '%s.right.bed' % classname

        left_bed = self.blastn(left_ter_reduce_file).saveas(left_bedfile_name)
        right_bed = self.blastn(right_ter_reduce_file).saveas(right_bedfile_name)

        ### to pair left and right terminal signals
        max_orf_length = sorted(ORF_length_list)[-1]
        distance_na = int(self.window) * 2 + max_orf_length
        sys.stdout.write('The length of %s is expected to be shorter than %s.\n' % (classname, str(distance_na+100)))
        joint_terminal_bed = left_bed.window(right_bed, l=0, r=distance_na, sm=True, sw=True).saveas('%s.joint' % classname)

        joint_terminal_dict = defaultdict(list)

        ## Do filteration, should filter out the case that occur only once.
        for line in joint_terminal_bed:
            ## To remove the case whose left and right signal are seriously overlaped
            if self.intersect(line[1:3], line[7:9]):
                continue
            left_name, right_name = line[3], line[9]
            if pair:
                if right_name not in collapsed_pair_dict[left_name]:  ## This will be considered as fake joint as they don't form terminal inverted repeats.
                    continue
            combiname = (left_name, right_name)
            loc = sorted([int(line[1]), int(line[2]), int(line[7]), int(line[8])])
            pairname = '-'.join([line[3], line[9]])
            score = str(float(line[4]) + float(line[10]))
            new_line = BT.create_interval_from_list([line[0], loc[0], loc[-1], pairname, 1, line[5], 1, score])
            joint_terminal_dict[combiname].append(new_line)

        ## To filter out the case whose occurence time is less than two
        valid_combind_list = []
        joint_terminal_list, monomer_terminal_list = [], []
        for combinid in joint_terminal_dict:
            sub_list = joint_terminal_dict[combinid]
            if len(sub_list) > 1:
                joint_terminal_list.extend(sub_list)
                valid_combind_list.append(combinid)
            else:
                monomer_terminal_list.extend(sub_list)

        if not joint_terminal_list and not monomer_terminal_list:
            os.chdir('../')
            return [], []

        del joint_terminal_dict
        gc.collect()

        ## To output the filtered window joint bed file
        joint_file = '%s.joint.filtered.bed' % classname
        BT.BedTool([line for line in joint_terminal_list]).saveas(joint_file)
        monomer_terminal_file = '%s.monomer.filtered.bed' % classname
        monomer_joint_bed = BT.BedTool([line for line in monomer_terminal_list]).sort().saveas(monomer_terminal_file)

        #### To find autonomous helitrons ###########
        significant_joint_bed = BT.BedTool([])
        ORF_sig_list = []
        if os.path.getsize(joint_file):
            opt_pvalue = '%s.joint.pvalue.bed' % classname
            sys.stdout.write("Begin to run fisher's exact test for %s!\n" % classname)
            subed_dir = '_'.join([classname, 'SubBed'])
            half_distance = int(round(int(distance_na) / 2, 0))

            joint_filepath_df = self.split_joint(left_bed, right_bed, valid_combind_list, subed_dir, half_distance)
            
            # To run fisher's exact text to select the co-occured left and right terminal signals
            fisher_program = subprocess.Popen(
                ['Rscript', FISHER_PRO, BEDTOOLS_PATH,
                self.genome_size, joint_filepath_df, opt_pvalue, str(self.process_num), str(self.pvalue), joint_file],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            fisher_program.wait()
            sys.stdout.write("Fisher's exact test finished for %s!\n" % classname)

            if os.path.exists(opt_pvalue):
                significant_joint_bed = BT.BedTool(opt_pvalue)
                ## candidates whose terminal signals could be found.
                ## ORF intersect with significant signals.
                ORF_sig_intersection = ORF_bed.intersect(significant_joint_bed, f=1, wo=True, s=True)
                ORF_sig_list.extend(list({BT.create_interval_from_list(line[6:14]) for line in ORF_sig_intersection}))

        BT.BedTool(ORF_sig_list).sort().saveas('%s_RC.auto.bed' % classname)

        ## candidates whose terminal signals occurred once.
        ORF_without_sigter_bed = ORF_bed.intersect(significant_joint_bed, v=True, s=True, wa=True)
        monomer_bed_intersection = ORF_without_sigter_bed.intersect(monomer_joint_bed, f=1, wo=True, s=True)
        monomer_bed_dict = defaultdict(list)
        [monomer_bed_dict[i.name].append(list(i[6:14])) for i in monomer_bed_intersection]

        for orfid in monomer_bed_dict:
            interval_list = monomer_bed_dict[orfid]
            candidate = sorted(interval_list, key=lambda x: (int(x[2]) - int(x[1])))[0]  ## Try to select the shortest one.
            ORF_sig_list.append(BT.create_interval_from_list(candidate))
            
        ## candidates whose terminal signals could not be found.
        monomer_bed = BT.BedTool(ORF_sig_list).sort()
        ORF_without_ter_bed=ORF_without_sigter_bed.intersect(monomer_bed, v=True, s=True, wa=True)
        [ORF_sig_list.append(BT.create_interval_from_list([line[0], line[1], line[2], line[3], 1, line[5], 1, '0'])) for line in ORF_without_ter_bed]

        RC_with_orf = BT.BedTool(ORF_sig_list).sort().saveas('%s_RC.auto.bed' % classname)
        non_autonomous = significant_joint_bed.intersect(RC_with_orf, f=0.8, wa=True, v=True).saveas('%s_RC.nonauto.bed' % classname)

        if pair and classname != 'Helitron':  ## To find stem loop structure of Helentron
            RC_with_orf = self.heltentron_terminal(RC_with_orf).saveas('%s_RC.auto.bed' % classname)
            non_autonomous = self.heltentron_terminal(non_autonomous).saveas('%s_RC.nonauto.bed' % classname)

        auto_alternative_list = self.merge_overlaped_intervals(RC_with_orf, classname, 'auto')
        nonauto_alternative_list = self.merge_overlaped_intervals(non_autonomous, classname, 'nonauto')

        nonauto_alternative_list.extend(auto_alternative_list)
        os.chdir('../')
        return nonauto_alternative_list

    def malign(self, fafile):
        if os.path.getsize(fafile):
            aln_file = ''.join([fafile, '.aln'])
            with open(aln_file, 'w') as mf:
                mul_aln = subprocess.Popen(["mafft", "--auto", "--quiet", fafile], stdout=mf)
                mul_aln.wait()

    def flanking_seq(self, pairlist):
        subwkdir = 'boundary_align'
        pairname = pairlist[0][3]
        left_extend_file = ''.join([subwkdir, '/', pairname, '.left.fa'])
        right_extend_file = ''.join([subwkdir, '/', pairname, '.right.fa'])

        candidate_bed = BT.BedTool(
            [BT.create_interval_from_list([line[0], line[1], line[2], line[3], line[7], line[5]]) for line in pairlist])
        candidate_bed = self.merge_bedfile(candidate_bed, orf=False, window=100)
        candidate_list = sorted([list(line) for line in candidate_bed], key=lambda x: -float(x[4]))
        if len(candidate_list) < 2:
            return 0

        self.filepath_list.append(left_extend_file)
        self.filepath_list.append(right_extend_file)

        selected_list = candidate_list[:20]
        with open(left_extend_file, 'w') as LF, open(right_extend_file, 'w') as RF:
            for line in selected_list:
                chrmid, start, stop, name, score, strand = line
                if strand == '+':
                    if int(start) - 50 >= 0:
                        seq_stop = int(start)
                        seq_start = int(start) - 50
                        seq = self.genome_dict[chrmid][seq_start:seq_stop]
                        LF.write(''.join(['>', chrmid, '-', str(seq_start), '-', str(seq_stop), '\n']))
                        LF.write(str(seq))
                        LF.write('\n')

                    if int(stop) + 50 <= self.chrm_size[chrmid]:
                        seq_start = int(stop)
                        seq_stop = int(stop) + 50
                        seq = self.genome_dict[chrmid][seq_start:seq_stop]
                        RF.write(''.join(['>', chrmid, '-', str(seq_start), '-', str(seq_stop), '\n']))
                        RF.write(str(seq))
                        RF.write('\n')
                else:
                    if int(stop) + 50 <= self.chrm_size[chrmid]:
                        seq_start = int(stop)
                        seq_stop = int(stop) + 50
                        seq = self.genome_dict[chrmid][seq_start:seq_stop].reverse_complement()
                        LF.write(''.join(['>', chrmid, '-', str(seq_start), '-', str(seq_stop), '\n']))
                        LF.write(str(seq))
                        LF.write('\n')
                    if int(start) - 50 >= 0:
                        seq_stop = int(start)
                        seq_start = int(start) - 50
                        seq = self.genome_dict[chrmid][seq_start:seq_stop].reverse_complement()
                        RF.write(''.join(['>', chrmid, '-', str(seq_start), '-', str(seq_stop), '\n']))
                        RF.write(str(seq))
                        RF.write('\n')

    def MakeSelection(self, alternative_list):
        # This function is to filter out the candidates who might insert into other superfamily of transposons. 
        sys.stdout.write('Begin to run boundary check program.\n')
        subwkdir = 'boundary_align'
        if not os.path.exists(subwkdir):
            os.mkdir(subwkdir)
        else:
            shutil.rmtree(subwkdir)
            os.mkdir(subwkdir)
        alternative_dict = defaultdict(list)
        for line in alternative_list:
            alternative_dict[line[3]].append(line)

        new_alternative_list = [alternative_dict[pairname] for pairname in alternative_dict]

        ## To output the flanking sequences
        self.filepath_list = []
        planpool = ThreadPool(int(self.process_num))
        for pairlist in new_alternative_list:
            planpool.apply_async(self.flanking_seq, args=(pairlist,))
        planpool.close()
        planpool.join()

        filepath_list = [i for i in self.filepath_list if os.path.getsize(i)]
        planpool = ThreadPool(int(self.process_num))
        for fafile in filepath_list:
            planpool.apply_async(self.malign, args=(fafile,))
        planpool.close()
        planpool.join()

        boundary_identity_tbl = 'Boundary.identity.tbl'
        Boundary_check_pro = subprocess.Popen(
            ['Rscript', BOUNDARY_PRO, os.path.abspath(subwkdir), os.path.abspath(boundary_identity_tbl)],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        Boundary_check_pro.wait()
        sys.stdout.write("Boundary check was finished!\n")

        name_iden_dict = defaultdict(lambda :defaultdict(lambda :0))
        if os.path.exists(boundary_identity_tbl):
            with open(boundary_identity_tbl, 'r') as F:
                F.readline() ## Skip header
                for line in F:
                    name, direction, iden = line.rstrip().split('\t')
                    name_iden_dict[name][direction] = float(iden)

        insertion_name_dict = defaultdict(list)
        [insertion_name_dict[line[10]].append(line) for line in alternative_list]
        RC_replist = []

        non_auto_pairname_list = []
        for insertion in insertion_name_dict:  ## To do selection for autonomous firstly, because the nonautonomous counterparts will be dertermined by autonomous ones.
            candidate_list = []
            for line in insertion_name_dict[insertion]:
                pairname = line[3]
                classname= line[8]
                if name_iden_dict[pairname]['left'] < 0.7:
                    if classname == 'Helitron':
                        if name_iden_dict[pairname]['right'] < 0.7:
                            candidate_list.append(line)
                    else:  ## For helentron just detect the left side.
                        candidate_list.append(line)

            if not candidate_list and line[9] == 'auto':
                ## Might represent a truncated helitron insertion. Just keep the autonomous ones and remove the nonautonomous ones.
                ## Because we are sure that there should be one insertion in autonomous region.
                candidate_list = insertion_name_dict[insertion]
            candidate_list = sorted(candidate_list, key=lambda x:[float(x[6]), int(x[1])-int(x[2]), -float(x[7])])   ## pvalue, length, then bitscore
            if candidate_list:
                RC_replist.append(candidate_list[0])
        return RC_replist

    def main(self):
        RC_total_candidate = self.autonomous_detect()
        RC_total_candidate_dict = defaultdict(list)
        [RC_total_candidate_dict[line[4]].append(line) for line in RC_total_candidate]
        total_alternative_list = []

        if os.path.exists('Pairname.siginfo.tbl'):
            os.remove('Pairname.siginfo.tbl')
            
        for class_name in RC_total_candidate_dict:
            RC_list = RC_total_candidate_dict[class_name]
            if re.findall('Helentron|Helitron2', class_name):  ## need to pair the terminal repeats
                alter_list = self.prepare_terminal_seq(RC_list, pair=True, classname=class_name)
            else:
                if not Args.pair_helitron:
                    alter_list = self.prepare_terminal_seq(RC_list, pair=False, classname=class_name)
                else:
                    alter_list = self.prepare_terminal_seq(RC_list, pair=True, classname=class_name)
            total_alternative_list.extend(alter_list)

        total_alternative_list = [i for i in total_alternative_list if i]
        total_alternative_list = sorted(total_alternative_list, key=lambda x: [x[0], int(x[1])])
        BT.BedTool([BT.create_interval_from_list(line) for line in total_alternative_list]).saveas('RC.alternative.bed')
        RC_rep_list = self.MakeSelection(total_alternative_list)
        RC_rep_list = sorted(RC_rep_list, key=lambda x: [x[0], int(x[1])])
        RC_repbed = BT.BedTool([BT.create_interval_from_list(line) for line in RC_rep_list]).saveas('RC.representative.bed')
        for dir in os.listdir('./'):
            if os.path.isdir(dir) and dir.startswith('Hel'):
                orffile = ''.join([dir, '/', dir, '_orf.bed'])
                if os.path.exists(orffile):
                    orfbed = BT.BedTool(orffile).sort()
                    remained_orf_bed = orfbed.intersect(RC_repbed, v=True, f=1, wa=True)
                    i = 0
                    for line in remained_orf_bed:
                        RC_rep_list.append([line[0], str(line[1]), str(line[2]), ''.join([dir, '_orf', str(i)]), '0', line[5], '1', '0', dir, 'orf_only', ''.join(['orf_only', str(i)])])
                        i+=1
        RC_rep_list = sorted(RC_rep_list, key=lambda x: [x[0], int(x[1])])
        RC_repbed = BT.BedTool([BT.create_interval_from_list(line) for line in RC_rep_list]).saveas(
            'RC.representative.bed')
        os.remove(self.CTRR_stem_loop_description)
        #os.remove(self.subtir_description)
        shutil.rmtree(self.bedtoolstmp)
        if os.path.exists('GenomeDB'):
            shutil.rmtree('GenomeDB')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="HELA can detect and classify different variants of Helitron-like elements: Helitron, Helentron and Helitron2. Please visit https://github.com/Zhenlisme/HELA/ for more information. Email us: zhen.li3@universite-paris-saclay.fr")
    parser.add_argument("-g", "--genome", type=str, required=True, help="The genome file in fasta format.")
    parser.add_argument("-w", "--window", type=int, default=10000, required=False,
                        help="To check terminal signals within a given window bp upstream and downstream of ORF ends, default is 10 kb.")
    parser.add_argument("-dm", "--distance_domain", type=int, default=2500, required=False,
                        help="The distance between HUH and Helicase domain, default is 2500.")
    parser.add_argument("-pt", "--pair_helitron", type=int, default=0, required=False, choices=[0, 1],
                        help="For Helitron, its 5' and 3' terminal signal pairs should come from the same autonomous helitorn or not. 0: no, 1: yes. default no.")
    parser.add_argument("-is1", "--IS1", type=int, default=0, required=False, choices=[0, 1],
                        help="Set the insertion site of autonomous Helitron as A and T. 0: no, 1: yes. default no.")
    parser.add_argument("-is2", "--IS2", type=int, default=0, required=False, choices=[0, 1],
                        help="Set the insertion site of autonomous Helentron/Helitron2 as T and T. 0: no, 1: yes. default no.")
    parser.add_argument("-sim_tir", "--simtir", type=int, default=100, required=False, choices=[100, 90, 80],
                        help="Set the simarity between short inverted repeats(TIRs) of Helitron2/Helentron. Default 100.")
    parser.add_argument("-p", "--pvalue", type=float, required=False, default=1e-3, help="The p-value for fisher's exact test.")
    parser.add_argument("-o", "--opdir", type=str, required=True, help="The output directory.")
    parser.add_argument("-n", "--process", type=int, default=2, required=False, help="Maximum of threads to be used.")
    parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0.0')
    Args = parser.parse_args()

    ## To set and check dependency file path
    HMMFILE = '_HMM_'
    HEADERFILE = '_HEADER_'
    FISHER_PRO = '_FISHER_'
    BOUNDARY_PRO = '_BOUNDARY_'
    try:
        BEDTOOLS_PATH = subprocess.check_output("which bedtools", shell=True).decode().rstrip()
        BEDTOOLS_PATH = '/'.join(BEDTOOLS_PATH.split('/')[:-1])
    except:
        print("Could not find bedtools path.")
        exit(0)

    try:
        subprocess.check_output("which rnabob", shell=True)
    except:
        print("Could not find rnabob path.")
        exit(0)

    try:
        subprocess.check_output("which cd-hit-est", shell=True)
    except:
        print("Could not find cd-hit-est path.")
        exit(0)
    
    try:
        subprocess.check_output("which mafft", shell=True)
    except:
        print("Could not find mafft path.")
        exit(0)
    try:
        subprocess.check_output("which hmmsearch", shell=True)
    except:
        print("Could not find hmmsearch path.")
        exit(0)

    try:
        subprocess.check_output("which getorf", shell=True)
    except:
        print("Could not find getorf path.")
        exit(0)

    try:
        subprocess.check_output("which gt", shell=True)
    except:
        print('Could not find genometools path.')
        exit(0)

    try:
        subprocess.check_output("which dialign2-2", shell=True)
    except:
        print('Could not find dialign2 path.')
        exit(0)

    try:
        subprocess.check_output("which blastn", shell=True)
    except:
        print('Could not find genometools path.')
        exit(0)

    if not os.path.exists(HMMFILE):
        print("Hmmer model file not found!")
        exit(0)
    if not os.path.exists(HEADERFILE):
        print('header lcv file not found!')
        exit(0)
    if not os.path.exists(BOUNDARY_PRO):
        print('Boundary check program not found!')
        exit(0)
    if not os.path.exists(FISHER_PRO):
        print("Fisher's exact test program not found!")
        exit(0)

    HomoSearch = Homologous_search(HMMFILE, os.path.abspath(Args.genome), os.path.abspath(Args.opdir),
                                   HEADERFILE, Args.window, Args.distance_domain, Args.pvalue, Args.process)
    HomoSearch.main()

