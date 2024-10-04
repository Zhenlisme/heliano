#!_INTERPRETERPYTHON_PATH_

import os, re, subprocess, sys, argparse, shutil, random, gc
from Bio import SeqIO
from multiprocessing.pool import ThreadPool
from collections import defaultdict
import pybedtools as BT

"""
This program is supposed to detect and classify different variants of Helitron-like elements: HLE1 and HLE2. Please follow the homepage of this software for more information: https://github.com/Zhenlisme/heliano.   
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
                    ## To avoid rnabob bugs
                    if int(splitline[1]) < 0:
                        print(splitline)
                        continue
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
                        ## length of inverted sequences should be greater than 11 and shorter than 18
                        if invt_length_left >= 12 and invt_length_right >= 12 and invt_length_left <= 17 and invt_length_right <= 17:
                            invt_list.append([chrmid, str(left_start), str(right_end), left_expand, right_expand,
                                              (invt_length_right + invt_length_left) / 2, sim])

        invt_list = sorted(invt_list, key=lambda x: int(x[1]))
        os.remove(invttirfile)
        os.system('rm %s*' % dbname)
        return invt_list

# define Homologous_search class to find Helitron-like transposase domain and theri auto/non-auto relatives
class Homologous_search:
    def __init__(self, rep_hel_hmm, genome, wkdir, headerfile, window, distance_domain, distance_na, pvalue, process_num):
        self.rep_hel_hmm = rep_hel_hmm
        self.genome = genome
        self.genome_dict = SeqIO.parse(genome, 'fasta')
        self.genome_dict = {k.id: k.seq.upper() for k in self.genome_dict}
        self.process_num = int(process_num)
        self.wkdir = wkdir
        self.headerpatternfile = headerfile
        self.window = window
        self.distance_domain = distance_domain
        self.distance_na = defaultdict(lambda :int(distance_na))
        self.pvalue = float(pvalue)
        self.cutoff_flank = float(Args.flank_sim)
        if Args.terminal_sequence:
            self.pairfile = os.path.abspath(Args.terminal_sequence)
            sys.stdout.write('You added the pairfile %s.\n' % self.pairfile)
            if not os.path.isfile(self.pairfile):
                sys.stderr.write("Error: The pair list file doesn't exist, please check!\n")
                exit(0)
        # To transform the header file to regular expression.
        with open(headerfile, 'r') as F:
            headerpattern_list = F.read().rstrip().split('\n')
            self.headerpattern = '|'.join(headerpattern_list) if not Args.IS1 else '|'.join([''.join(['A', i]) for i in headerpattern_list])
            #self.headerpattern = ''.join(['(?=(', self.headerpattern, '))'])
        if not os.path.exists(self.wkdir):
            os.mkdir(self.wkdir)
            os.chdir(self.wkdir)
        else:
            sys.stderr.write('Error: Directory %s exists. Please change!\n' % self.wkdir)
            exit(0)

        ## To check the pair list file
        self.terminalfile_dict = defaultdict(lambda: defaultdict(dict))
        self.prepair_dict = defaultdict(list)
        if Args.terminal_sequence:
            pairdict = defaultdict(list)
            with open(self.pairfile, 'r') as F:
                for line in F:
                    splitline = line.rstrip().split('\t')
                    pairdict[splitline[0]].append(splitline[1:])
            pairdir = 'Pre_pair/'
            if not os.path.exists(pairdir):
                os.mkdir(pairdir)
            for classname in pairdict:
                leftfile = os.path.abspath(''.join([pairdir, classname, '.left.pre.fa']))
                rightfile = os.path.abspath(''.join([pairdir, classname, '.right.pre.fa']))
                self.terminalfile_dict[classname]['left'] = leftfile
                self.terminalfile_dict[classname]['right'] = rightfile
                recorder_dict = {}
                with open(leftfile, 'w') as left_w, open(rightfile, 'w') as right_w:
                    for line in pairdict[classname]:
                        leftname, leftseq, rightname, rightseq = line
                        leftname, rightname = ''.join([leftname, 'pre']), ''.join([rightname, 'pre'])
                        ## To avoid repeatly writing
                        if leftname not in recorder_dict:
                            left_w.write(''.join(['>', leftname, '\n', leftseq, '\n']))
                        if rightname not in recorder_dict:
                            right_w.write(''.join(['>', rightname, '\n', rightseq, '\n']))
                        self.prepair_dict[leftname].append(rightname)
                        recorder_dict[leftname] = 1
                        recorder_dict[rightname] = 1

        if Args.dis_denovo:
            if not self.terminalfile_dict:
                sys.stderr.write('Error: The pair list file is either not specified or empty. See parameter "-ts".\n')
                exit(0)
            else:
                sys.stdout.write(
                    'You will not search for the terminal structures of HLE in a de-novo way, but by using the pair file: %s.\n' % self.pairfile)
        self.bedtoolstmp = os.path.abspath('BedtoolsTMP')
        if not os.path.exists(self.bedtoolstmp):
            os.mkdir(self.bedtoolstmp)
        BT.set_tempdir(self.bedtoolstmp)

        CWD = os.getcwd()
        self.genome_size = '%s/Genome.size' % CWD
        self.chrm_size = {i:len(self.genome_dict[i]) for i in self.genome_dict}
        genome_size = list(self.chrm_size.items())
        genome_size = sorted(genome_size, key=lambda x: x[0])
        with open(self.genome_size, 'w') as F:
            F.writelines([''.join([i[0], '\t', str(i[1]), '\n']) for i in genome_size])

        ## To determine the evalue for short-sequence blastn, set the bit-score cutoff as 30, the evalue cutoff should follow the formula: m*n/(2**30)
        sum_genomesize = sum([i[1] for i in genome_size])
        self.evalue_blastn = sum_genomesize * 30 / (2 ** int(Args.score))

        # To define stem_loop structure ending with CTRR motif of Helitron
        self.CTRR_stem_loop_description = '%s/CTRR_stem_loop.descr' % CWD
        CTRR_description = """r1 s1 r1' s2\nr1 1:1 NNNNN[10]:[10]NNNNN TGCA\ns1 0 N[7]\ns2 0 N[15]CTRR%s\n"""
        # Add 'T' in the end if user limitted the 'A-T' insertion site for Helitron.
        CTRR_description = CTRR_description % 'T' if Args.IS1 else CTRR_description % ''
        with open(self.CTRR_stem_loop_description, 'w') as F:
            F.write(CTRR_description)

        # To define stem_loop structure of  HLE2
        self.subtir_description = '%s/subtir_stem_loop.descr' % CWD
        # Add 'T' in the end if user limitted the 'T-T' insertion site for HLE2.
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

    def merge_bedfile(self, BedInput, window=1500):
        # Define function to merge two distance-close genomic features
        cluster_bed = BedInput.cluster(d=window, s=True)
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
            ## To determain sub_class, select the case with highest score
            sub_class = sorted(coordlist, key=lambda x: float(x[4]))[-1][3]
            ## To merge the ORF location
            orf_list = sorted([int(b) for i in coordlist for b in i[-1].split('-')])
            orf_coord = '-'.join([str(orf_list[0]), str(orf_list[-1])])
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
        merge_hel = self.merge_bedfile(Hel_bed, window=1500)
        merge_rep = self.merge_bedfile(Rep_bed, window=1500)
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
        extend_seq = []
        opbed_list = []
        extend_file = ''.join([helentron_bed, '.fa'])
        extend_dict = {}
        extend_dict = defaultdict(list)
        ## To output all extend rigions into one single file
        #init_num = 1
        with open(helentron_bed, 'r') as F:
            for line in F:
                feature = line.rstrip().split('\t')
                chrmid, start, stop, name, score, strand, pvalue, Bscore, classname, mobile_type, insertion_name = feature
                start, stop = int(start), int(stop)
                if float(Bscore) == 0:  ## without terminal signals
                    opbed_list.append([chrmid, str(start), str(stop), name, score, strand, pvalue, Bscore, classname, mobile_type, insertion_name])
                    continue

                if strand == '+':
                    extend_id = '-'.join([chrmid, str(start), str(stop), 'p'])
                    if '.2' not in classname:
                        detect_seq = str(self.genome_dict[chrmid][stop: stop + 80])
                        StemStart = stop
                    else:  ## Helentron.2, the stem loop is at 5'end
                        terminal_start = 0 if start - 80 < 0 else start - 80
                        detect_seq = str(self.genome_dict[chrmid][terminal_start: start])
                        StemStart = terminal_start
                    ## To remove short seequence
                    if len(detect_seq) < 17:
                        continue
                    ## To remove N rich seq
                    N_count = detect_seq.count('N')
                    if N_count >= 5:
                        continue
                    ## The last element is the initial start for terminal detection region
                    extend_dict[extend_id].append([chrmid, str(start), str(stop), name, score, strand, pvalue, Bscore,
                                              classname, mobile_type, insertion_name, StemStart])
                    extend_seq.append(''.join(['>', extend_id, '\n', detect_seq, '\n']))
                else:
                    extend_id = '-'.join([chrmid, str(start), str(stop), 'n'])
                    if '.2' not in classname:
                        terminal_start = 0 if start - 80 < 0 else start - 80
                        detect_seq = str(self.genome_dict[chrmid][terminal_start: start])
                        StemStart = terminal_start
                    else:  ## Helentron.2, the stem loop is at 3'end
                        detect_seq = str(self.genome_dict[chrmid][stop: stop + 80])
                        StemStart = stop
                    ## To remove short seequence
                    if len(detect_seq) < 17:
                        continue
                    ## To remove N rich seq
                    N_count = detect_seq.count('N')
                    if N_count >= 5:
                        continue
                    ## The last element is the initial start for terminal detection region
                    extend_dict[extend_id].append([chrmid, str(start), str(stop), name, score, strand, pvalue, Bscore,
                                              classname, mobile_type, insertion_name, StemStart])
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
                chrmid, start, stop, name, score, strand, pvalue, Bscore, classname, mobile_type, insertion_name, e_start = sublist
                if extend_id.endswith('p'):
                    stem_loop = [i for i in stem_loop_dict[extend_id] if i[5] == '+']
                    ## order by start position (closer), stem length (longer), loop length (shorter),  total length (shorter)
                    if '.2' not in classname:
                        stem_loop = sorted(stem_loop, key = lambda x: [x[0], int(x[1]), -int(x[3]), int(x[4]), int(x[2]) - int(x[1])])
                        if stem_loop:
                            stem_start, stem_stop = stem_loop[0][1:3]
                            stem_stop = int(stem_stop) + int(e_start)
                            opbed_list.append([chrmid, start, str(stem_stop), name, score, strand,
                                               pvalue, Bscore, classname, mobile_type, insertion_name])
                    else:
                        stem_loop = sorted(stem_loop, key=lambda x: [x[0], -int(x[2]), -int(x[3]), int(x[4]), int(x[2]) - int(x[1])])
                        if stem_loop:
                            stem_start, stem_stop = stem_loop[0][1:3]
                            stem_start = int(stem_start) + int(e_start)
                            opbed_list.append([chrmid, str(stem_start), stop, name, score, strand,
                                               pvalue, Bscore, classname, mobile_type, insertion_name])
                    if not stem_loop:
                        opbed_list.append([chrmid, start, stop, name, score, strand, pvalue, Bscore, classname, mobile_type, insertion_name])
                        continue
                else:
                    stem_loop = [i for i in stem_loop_dict[extend_id] if i[5] == '-']
                    ## order by end position (longer), stem length (longer), loop length (shorter), total length (shorter)
                    if '.2' not in classname:
                        stem_loop = sorted(stem_loop, key=lambda x: [x[0], -int(x[2]), -int(x[3]), int(x[4]), int(x[2]) - int(x[1])])
                        if stem_loop:
                            stem_start, stem_stop = stem_loop[0][1:3]
                            stem_start = int(stem_start) + int(e_start)
                            opbed_list.append([chrmid, str(stem_start), stop, name, score, strand,
                                               pvalue, Bscore, classname, mobile_type, insertion_name])
                    else:
                        stem_loop = sorted(stem_loop, key=lambda x: [x[0], int(x[1]), -int(x[3]), int(x[4]),
                                                                     int(x[2]) - int(x[1])])
                        if stem_loop:
                            stem_start, stem_stop = stem_loop[0][1:3]
                            stem_stop = int(stem_stop) + int(e_start)
                            opbed_list.append([chrmid, start, str(stem_stop), name, score, strand,
                                               pvalue, Bscore, classname, mobile_type, insertion_name])
                    if not stem_loop:
                        opbed_list.append([chrmid, start, stop, name, score, strand, pvalue, Bscore, classname, mobile_type, insertion_name])
                        continue
        ### To complement the candidates whose stem loop signal doesn't exist.
        remained_cases = set(extend_dict.keys()) - set(stem_loop_dict.keys())
        for key in remained_cases:
            for sublist in extend_dict[key]:
                opbed_list.append(sublist[:11])
        ### To creat bedtools objective
        #opbed = BT.BedTool([BT.create_interval_from_list(line) for line in opbed_list])
        with open(helentron_bed, 'w') as F:
            F.write('\n'.join(['\t'.join(line) for line in opbed_list]))
            F.write('\n')
        os.remove(extend_file)

    def intergrated_program(self, subgenome):
        # This function is used to recover terminal signals of Helitron-like elements (TIRs for HLE2; TC... motif and ...CTRR motif for Helitron)
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

            ## To decide class name
            if hel_name == 'HLE1' and rep_name == 'HLE1':
                classname = 'HLE1'
            elif hel_name == 'HLE2' and rep_name == 'HLE2':
                classname = 'HLE2'
            else:
                classname = '_or_'.join([rep_name, hel_name])
            ## To produce orf id which will be used as sole identifier of terminal signals.
            ORFID = '-'.join([ORF_chrmid, str(ORF_start), str(ORF_stop)])

            ## To decide to run denovo structural search or not. if not, just save the ORFs
            if Args.dis_denovo:
                RC_total_candidate.append(
                    (ORF_chrmid, str(ORF_start), str(ORF_stop), strand, classname, (), (), ORFID))
                continue

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
            if classname == 'HLE1':
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
                                          ORF_chrmid, str(ORF_start), str(ORF_stop), strand, classname, tuple(TC_list),
                                          tuple(Stem_loop_list), ORFID))

            # Candidate is HLE2
            else:
                # Size of terminal inverted sequences should be greater than the size of predicted transposase
                mini_dist_tir = ORF_stop - ORF_start
                # Size of terminal inverted sequences should be less than the size of extension
                max_dist_tir = 2 * int(self.window) + ORF_stop - ORF_start

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
                RC_total_candidate.append((ORF_chrmid, str(ORF_start), str(ORF_stop), strand, classname,
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
                    "Chrm %s will not be used to detect autonomous HLEs as its length is shorter than 1000 bp\n" % chrm)
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
        sys.stdout.write('Start to search for HLE rep-hel blocks...\n')
        # Use python multiple threading
        planpool = ThreadPool(processnum)
        #processnum = processnum if self.cpu_count > processnum else self.cpu_count
        #planpool = Pool(processnum)
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

    def blastn(self, query_file, optdir):
        # To search for homologous of given sequences
        blastn_opt = ''.join([query_file, '.tbl'])
        blastn_pro = subprocess.Popen(
            ['blastn', '-db', self.genomedb, '-query', query_file, '-num_threads', str(self.process_num),
             '-max_target_seqs', '999999999', '-evalue', str(self.evalue_blastn), '-task', 'blastn-short',
             '-outfmt', '6 qseqid sseqid pident qstart qend sstart send evalue qlen bitscore', '-out', blastn_opt])
        blastn_pro.wait()
        ## if the blastn_opt is empty, return
        if not os.path.getsize(blastn_opt):
            return 0
        ## The blastn output is sorted by queryname. Split blastn outfmt-6 table into bed files.
        with open(blastn_opt, 'r') as F:
            ## init reading 
            splitlines = F.readline().rstrip().split('\t')
            init_qname = splitlines[0]
            chrm, identity, qstart, qend, sstart, send, evalue, qlen, bitscore = splitlines[1:10]
            if int(sstart) < int(send):
                START = sstart
                END = send
                strand = '+'
            else:
                START = send
                END = sstart
                strand = '-'
            deposit_list = [''.join([chrm, '\t', START, '\t', END, '\t', init_qname, '\t', bitscore, '\t', strand, '\n'])]
            for line in F:
                splitlines = line.rstrip().split('\t')
                chrm, identity, qstart, qend, sstart, send, evalue, qlen, bitscore = splitlines[1:10]
                if int(sstart) < int(send):
                    START = sstart
                    END = send
                    strand = '+'
                else:
                    START = send
                    END = sstart
                    strand = '-'
                query_name = splitlines[0]
                if query_name == init_qname:
                    deposit_list.append(''.join([chrm, '\t', START, '\t', END, '\t', init_qname, '\t', bitscore, '\t', strand, '\n']))
                else:
                    subfile = ''.join([optdir, '/', init_qname, '.bed'])
                    with open(subfile, 'w') as wF:
                        wF.writelines(deposit_list)
                    ## reinit
                    init_qname = query_name
                    deposit_list = [''.join([chrm, '\t', START, '\t', END, '\t', init_qname, '\t', bitscore, '\t', strand, '\n'])]
        ## The last deposit need to be saved manually
        subfile = ''.join([optdir, '/', init_qname, '.bed'])
        with open(subfile, 'w') as wF:
            wF.writelines(deposit_list)
        os.remove(blastn_opt)
        return 1

    def merge_overlaped_intervals(self, bedinput, represent_bed, alt_optbed, rep_type = True):
        # To merge overlaped nonautonomous candidates
        True_pair_list = []
        with open(bedinput, 'r') as RF, open(represent_bed, 'a') as WFr, open(alt_optbed, 'a') as WFalt:
            ### To output alternative insertions and extract first line
            EMPTY_DEDUCE=0
            for line in RF:
                EMPTY_DEDUCE=1
                splitlines = line.rstrip().split('\t')
                chrmid, start, stop, pairname, count, strand, pvalue, Bscore, classname, mobile_type = splitlines
                mobile_type, altype = mobile_type.split('-')
                if altype == 'alt' and rep_type:
                    WFalt.write(line)
                    continue
                else:
                    break
            if not EMPTY_DEDUCE:
                return {}
            alternative_bedlines = [[chrmid, start, stop, pairname, count, strand, pvalue, Bscore, classname, mobile_type]]
            compared_line = (chrmid, start, stop, strand)
            blockname_init = 1

            output_recorder = {}
            for line in RF:
                chrmid, start, stop, pairname, count, strand, pvalue, Bscore, classname, mobile_type = line.rstrip().split('\t')
                mobile_type, altype = mobile_type.split('-')
                ### To output alternative insertions
                if altype == 'alt' and rep_type:
                    WFalt.write(line)
                    continue

                newline = [chrmid, start, stop, pairname, count, strand, pvalue, Bscore, classname, mobile_type]
                Intersection_deduce = self.intersect(compared_line[1:3], [start, stop], lportion=0.8, rportion=0.8, bool_and=0)
                if compared_line[0] == chrmid and Intersection_deduce:
                    merged_coord = sorted([compared_line[1], compared_line[2], start, stop], key=lambda x: int(x))
                    compared_line_s, compared_line_e = merged_coord[0], merged_coord[-1]
                    compared_line = [chrmid, start, compared_line_e, strand]
                    alternative_bedlines.append(newline)
                else:
                    block_name = '_'.join(['insertion', classname, mobile_type, str(blockname_init)])
                    output_recorder[block_name]=1
                    ## Begin to output the former to alternative file
                    [i.append(block_name) for i in alternative_bedlines]
                    if mobile_type=='auto':
                        print('blocks', alternative_bedlines)
                    ## Try to select a representative from these alternatives.
                    RC_list = self.filter(alternative_bedlines, classname=classname, mobile_type=mobile_type)
                    if RC_list:
                        for line in RC_list:
                            WFr.write('\t'.join(line))
                            WFr.write('\n')
                            True_pair_list.append(line[3])

                    blockname_init += 1
                    alternative_bedlines = [newline]
                    compared_line = [chrmid, start, stop, strand]

            ## in case the last overlaped series not saved.
            block_name = '_'.join(['insertion', classname, mobile_type, str(blockname_init)])
            if block_name not in output_recorder:
                [i.append(block_name) for i in alternative_bedlines]
                RC_list = self.filter(alternative_bedlines, classname=classname, mobile_type=mobile_type)
                if RC_list:
                    for line in RC_list:
                        WFr.write('\t'.join(line))
                        WFr.write('\n')
                        True_pair_list.append(line[3])
        return set(True_pair_list)

    def merge_overlaped_autos(self, bedinput, represent_bed, alt_optbed, rep_type = True):
        True_pair_list = []

        with open(bedinput, 'r') as RF, open(represent_bed, 'a') as WF, open(alt_optbed, 'a') as WFalt:
            ### To output alternative insertions and extract first line
            EMPTY_DEDUCE=0
            for line in RF:
                EMPTY_DEDUCE=1
                splitlines = line.rstrip().split('\t')
                chrmid, start, stop, pairname, count, strand, pvalue, Bscore, classname, mobile_type = splitlines[6:16]
                mobile_type, altype = mobile_type.split('-')
                splitlines[15] = mobile_type
                if altype == 'alt' and rep_type:
                    WFalt.write(line)
                    continue
                else:
                    break
            if not EMPTY_DEDUCE:
                return {}
        
            last_orfid = splitlines[3]
            alternative_bedlines = [splitlines[6:16]]
            blockname_init = 1
            output_recorder = {}
            for line in RF:
                splitlines = line.rstrip().split('\t')
                orfid = splitlines[3]
                mobile_type=splitlines[15]
                mobile_type, altype = mobile_type.split('-')
                splitlines[15] = mobile_type
                ### To output alternative insertions
                if altype == 'alt' and rep_type:
                    WFalt.write(line)
                    continue
                if orfid == last_orfid:
                    alternative_bedlines.append(splitlines[6:16])
                else:
                    block_name = '_'.join(['insertion', classname, mobile_type, str(blockname_init)])
                    output_recorder[block_name] = 1
                    ## Begin to output the former to alternative file
                    [i.append(block_name) for i in alternative_bedlines]
                    ## Try to select a representative from these alternatives.
                    RC_list = self.filter(alternative_bedlines, classname=classname, mobile_type=mobile_type)
                    if RC_list:
                        for line in RC_list:
                            WF.write('\t'.join(line))
                            WF.write('\n')
                            True_pair_list.append(line[3])
                    blockname_init += 1
                    last_orfid = orfid
                    alternative_bedlines = [splitlines[6:16]]
            ## in case the last overlaped series not saved.
            block_name = '_'.join(['insertion', classname, mobile_type, str(blockname_init)])
            if block_name not in output_recorder:
                [i.append(block_name) for i in alternative_bedlines]
                RC_list = self.filter(alternative_bedlines, classname=classname, mobile_type=mobile_type)
                if RC_list:
                    for line in RC_list:
                        WF.write('\t'.join(line))
                        WF.write('\n')
                        True_pair_list.append(line[3])
        return set(True_pair_list)

    def add_unique_pre_ts(self, pre_lts_fa, pre_rts_fa, new_lts_fa, new_rts_fa):
        pre_lts_dict = SeqIO.parse(pre_lts_fa, 'fasta')
        pre_rts_dict = SeqIO.parse(pre_rts_fa, 'fasta')
        pre_lts_dict = {k.id: k.seq.upper() for k in pre_lts_dict}
        pre_rts_dict = {k.id: k.seq.upper() for k in pre_rts_dict}
        if os.path.isfile(new_lts_fa) and os.path.isfile(new_rts_fa):
            sum_lts_length = sum([len(pre_lts_dict[i]) for i in pre_lts_dict])
            lts_evalue_blastn = sum_lts_length * 30 / (2 ** 30)
            lts_blastn_opt = 'lts.blastn.tbl'
            lts_blastn_pro = subprocess.Popen(
                ['blastn', '-subject', new_lts_fa, '-query', pre_lts_fa,
                '-max_target_seqs', '5', '-evalue', str(lts_evalue_blastn), '-task', 'blastn-short',
                '-outfmt', '6 qseqid sseqid pident qstart qend sstart send evalue qlen bitscore', '-out', lts_blastn_opt])
            lts_blastn_pro.wait()

            sum_rts_length = sum([len(pre_rts_dict[i]) for i in pre_rts_dict])
            rts_evalue_blastn = sum_rts_length * 30 / (2 ** 30)
            rts_blastn_opt = 'rts.blastn.tbl'
            rts_blastn_pro = subprocess.Popen(
                ['blastn', '-subject', new_rts_fa, '-query', pre_rts_fa,
                '-max_target_seqs', '5', '-evalue', str(rts_evalue_blastn), '-task', 'blastn-short',
                '-outfmt', '6 qseqid sseqid pident qstart qend sstart send evalue qlen bitscore', '-out', rts_blastn_opt])
            rts_blastn_pro.wait()

            with open(lts_blastn_opt, 'r') as F:
                overlap_lts = {line.split('\t')[0] for line in F}
            with open(rts_blastn_opt, 'r') as F:
                overlap_rts = {line.split('\t')[0] for line in F}
        else: ## the de-novo structural files are empty
            overlap_lts = set()
            overlap_rts = set()

        unique_lts = list(set(pre_lts_dict.keys()) - overlap_lts)
        unique_rts = list(set(pre_rts_dict.keys()) - overlap_rts)

        ## To find unique pairs
        unique_pair = []
        for lts in unique_lts:
            for rts in self.prepair_dict[lts]:
                if rts in unique_rts:
                    unique_pair.append((lts, rts))
        ## To pend
        with open(new_lts_fa, 'a') as lts_w, open(new_rts_fa, 'a') as rts_w:
            for pair in unique_pair:
                left, right = pair
                lts_w.write(''.join(['>', left, '\n', str(pre_lts_dict[left]), '\n']))
                rts_w.write(''.join(['>', right, '\n', str(pre_rts_dict[right]), '\n']))
        return unique_pair

    def transform_orfbedfile(self, orffile, classname):
        orf_alternative_file = ''.join([classname, '.orfonly.alternative.bed'])
        with open(orffile, 'r') as F, open(orf_alternative_file, 'w') as wf:
            init = 1
            for line in F:
                splitlines = line.rstrip().split('\t')
                name = ''.join(['insertion_', classname, '_orfonly_', str(init)])
                newline = '\t'.join(
                    [splitlines[0], splitlines[1], splitlines[2], splitlines[3], '1', splitlines[5], '1', '0', classname, 'orfonly', name])
                wf.write(newline)
                wf.write('\n')
                init += 1

    def prepare_terminal_seq(self, Helitron_list, pair=False, classname='HLE1'):
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
                    
                if classname == 'HLE1' and pair:  ## pair the left and right terminals that ever appears on the same Helitron region.
                    [Helitron_pair_list.append((left, tuple(right_ter_name))) for left in left_ter_name]  ## all possibilities.
                elif classname != 'HLE1':
                    [Helitron_pair_list.append((left_ter_name[i], right_ter_name[i])) for i in range(len(left_ter_name))]
                leftF.writelines(left_seq)
                rightF.writelines(right_seq)

        ## if run denovo search for motifs:
        if not Args.dis_denovo:
            if not left_exist or not right_exist: ## means that neither left nor right terminal signals exist, so skip blastn
                self.transform_orfbedfile(ORF_bedfile, classname=classname)
                os.chdir('../')
                return []
            else: ### To redunce redundency via cd-hit-est
                left_cluster_dict = self.cdhitest_clust(left_ter_file, left_ter_reduce_file,
                                                        helitron_type=''.join([classname, '_left']), id=0.9)
                right_cluster_dict = self.cdhitest_clust(right_ter_file, right_ter_reduce_file,
                                                         helitron_type=''.join([classname, '_right']), id=0.9)
        else: ## No need to run de-novo terminal structural detection
            left_cluster_dict = {}
            right_cluster_dict = {}

        pre_left_file = self.terminalfile_dict[classname]['left']
        pre_right_file = self.terminalfile_dict[classname]['right']
        unique_pairlist = []
        if pre_left_file: ### To pend unique-pre-ts signals
            unique_pairlist = self.add_unique_pre_ts(pre_left_file, pre_right_file, left_ter_reduce_file, right_ter_reduce_file)
        else:
            if Args.dis_denovo:
                ## This means no terminal structures are avaliable. Just output the orf information.
                self.transform_orfbedfile(ORF_bedfile, classname=classname)
                os.chdir('../')
                return []
        ## To make left-right pairs
        combinid_file = '%s.combinid.txt' % classname
        if pair:
            ### To get the collapsed cluster pairs
            collapsed_pair_dict = defaultdict(list)
            terminal_pair_dict = dict(Helitron_pair_list)
            if classname == 'HLE1':  ## pair the left and right terminals that ever appears on the same Helitron region.
                for left_name in terminal_pair_dict:
                    right_name_list = terminal_pair_dict[left_name]
                    for right_name in right_name_list:
                        collapsed_left_clust = left_cluster_dict[left_name]
                        collapsed_right_clust = right_cluster_dict[right_name]
                        collapsed_pair_dict[collapsed_left_clust].append(collapsed_right_clust)
            else:  ## pair inverted repeats
                for left_name in terminal_pair_dict:
                    right_name = terminal_pair_dict[left_name]
                    collapsed_left_clust = left_cluster_dict[left_name]
                    collapsed_right_clust = right_cluster_dict[right_name]
                    collapsed_pair_dict[collapsed_left_clust].append(collapsed_right_clust)
            # To reduce redundency
            for key in collapsed_pair_dict:
                collapsed_pair_dict[key] = list(set(collapsed_pair_dict[key]))
            with open(combinid_file, 'w') as F:
                for left in collapsed_pair_dict:
                    for right in collapsed_pair_dict[left]:
                        F.write(''.join([left, '\t', right, '\n']))
                ## To output unique-pre-ts pairs
                for pair in unique_pairlist:
                    F.write(''.join([pair[0], '\t', pair[1], '\n']))
        else:
            ## Just make all possibilities
            left_reduced_list = set(left_cluster_dict.values())
            right_cluster_dict = set(right_cluster_dict.values())
            with open(combinid_file, 'w') as F:
                for left in left_reduced_list:
                    for right in right_cluster_dict:
                        F.write(''.join([left, '\t', right, '\n']))
                ## To output unique-pre-ts pairs, not mixing.
                for pair in unique_pairlist:
                    F.write(''.join([pair[0], '\t', pair[1], '\n']))

        ## To evaluate the size range
        max_orf_length = sorted(ORF_length_list)[-1]
        distance_na = int(self.window) * 2 + max_orf_length
        half_distance = int(round(int(distance_na) / 2, 0))
        self.distance_na[classname] = distance_na if not self.distance_na[classname] else self.distance_na[classname]  ## If specified, will use the specified one.
        sys.stdout.write('The length of autonomous %s is expected to be shorter than %s.\n' % (classname, str(distance_na + 100)))
        sys.stdout.write('The length of nonautonomous %s is expected to be shorter than %s.\n' % (classname, str(self.distance_na[classname] + 100)))
        ## To run blastn to get homologies
        left_beddir, right_beddir = 'SubBlastnBed/%s_left' % classname, 'SubBlastnBed/%s_right' % classname
        if not os.path.exists(left_beddir):
            os.makedirs(left_beddir)
        if not os.path.exists(right_beddir):
            os.makedirs(right_beddir)
        lblastn_status = self.blastn(left_ter_reduce_file, left_beddir)
        rblastn_status = self.blastn(right_ter_reduce_file, right_beddir)
        ## To report if blastn runs well
        if lblastn_status and rblastn_status:
            sys.stdout.write('blastn runs well!\n')
        else:
            sys.stdout.write('No similar hits found for LTS/RTS for %s!\n' % classname)
        ## Prepare for fisher's exact test (avoid overloading, to split bedfiles )
        sys.stdout.write("Prepare for windowing for %s ...\n" % classname)
        subed_dir = '_'.join([classname, 'Windowing'])
        try:
            split_joint_program = subprocess.Popen(
                ['Rscript', SPLIT_JOINT_PRO, left_beddir, right_beddir, combinid_file, subed_dir,
                self.genome_size, str(half_distance), BEDTOOLS_PATH, str(self.process_num)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            split_joint_program.wait()
            joint_filepath_blastn = ''.join([subed_dir, '/', 'left_right.path.join'])
        except:
            sys.stderr.write("Error: Windowing program failed...\n")
            exit(0)

        ## To run fisher's exact text to select the co-occured left and right terminal signals
        sys.stdout.write("Begin to run fisher's exact test for %s!\n" % classname)

        bed_optdir = '%s_BedFisher' % classname
        if os.path.exists(bed_optdir):
            shutil.rmtree(bed_optdir)
        os.mkdir(bed_optdir)

        Strategy = '1' if Args.nearest else '0'
        try:
            fisher_program = subprocess.Popen(
                ['Rscript', FISHER_PRO, BEDTOOLS_PATH, self.genome_size, joint_filepath_blastn, bed_optdir,
                str(self.process_num), str(self.pvalue), ORF_bedfile, Strategy], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            fisher_program.wait()
            sys.stdout.write("Fisher's exact test finished for %s!\n" % classname)
        except:
            sys.stderr.write("Error: Fisher's exact test failed...\n")
            exit(0)

        FisherBed_files = ['/'.join([bed_optdir, file]) for file in os.listdir(bed_optdir)]

        ## To merge and sort all Fisher bed files
        fisher_pvalue_file = '%s.joint.pvalue.bed' % classname
        fisher_filelist = [''.join([bed_optdir, '/', file]) for file in os.listdir(bed_optdir) if file.endswith('.bed')]
        with open(fisher_pvalue_file, 'w') as WF:
            for subfisher in fisher_filelist:
                with open(subfisher, 'r') as F:
                    content = F.read()
                    WF.write(content)
        try:
            sort_pro = subprocess.Popen(
                ['bash', SORT_PRO, fisher_pvalue_file, str(self.process_num)],
                stdout=subprocess.DEVNULL)
            sort_pro.wait()
            sys.stdout.write("Merge and sort fisher result files finished for %s!\n" % classname)
        except:
            sys.stdout.write("Merge and sort fisher result files failed for %s!\n" % classname)
            exit(0)

        #############To annotate auto or non-autonomous ############
        ORF_bed = BT.BedTool(ORF_bedfile)
        ## To get autonomous candidates.
        fisher_pvalue_bed = BT.BedTool(fisher_pvalue_file)
        ## candidates whose significant terminal signals are able to be found.
        orf_fisher_bed = fisher_pvalue_bed.intersect(ORF_bed, nonamecheck=True, F=1, wa=True, s=True, c=True).saveas('Fisher_with_ORF.bed')
        ## To select the intervals that contain only one ORF
        auto_alternative_file = ''.join([classname, '.auto.alternative.bed'])
        
        with open(auto_alternative_file, 'w') as WF:
            with open('Fisher_with_ORF.bed', 'r') as f1:
                for line in f1:
                    splitlines = line.rstrip().split('\t')
                    if int(splitlines[-1])==1:
                        splitlines = line.rstrip().split('\t')
                        newline = '\t'.join(
                            [splitlines[0], splitlines[1], splitlines[2], splitlines[3], splitlines[4], splitlines[5],
                             splitlines[6], splitlines[7], classname, 'auto-%s'%splitlines[8]])
                        WF.write(newline)
                        WF.write('\n')

        ## To get non-autonomous candidates
        non_autonomous_bed = fisher_pvalue_bed.intersect(ORF_bed, nonamecheck=True, F=0.5, wa=True, v=True).saveas('Nonauto.bed')
        ## candidates whose significant terminal signals are unable to be found.
        auto_alternative_bed = BT.BedTool(auto_alternative_file)
        ORF_without_ter_bed = ORF_bed.intersect(auto_alternative_bed, nonamecheck=True, v=True, s=True, wa=True).saveas('OnlyORF.bed')
        self.transform_orfbedfile('OnlyORF.bed', classname=classname)
        os.chdir('../')

    def malign(self, fafile):
        if os.path.getsize(fafile):
            aln_file = ''.join([fafile, '.aln'])
            with open(aln_file, 'w') as mf:
                mul_aln = subprocess.Popen(["mafft", "--auto", "--quiet", fafile], stdout=mf)
                mul_aln.wait()

    def flanking_seq(self, fisherfile):
        subwkdir = 'boundary_align'
        pairname = os.path.basename(fisherfile).replace('.bed', '')
        left_extend_file = ''.join([subwkdir, '/', pairname, '.left.fa'])
        right_extend_file = ''.join([subwkdir, '/', pairname, '.right.fa'])
        terminal_file = ''.join([subwkdir, '/', pairname, '.terminal'])

        candidate_bed = BT.BedTool(fisherfile).sort()
        candidate_bed = candidate_bed.merge(d=100, c='4,5,6', o='first,mean,first')
        candidate_list = sorted([list(line) for line in candidate_bed], key=lambda x: -float(x[4]))
        if len(candidate_list) < 2:
            return 0

        self.filepath_list.append(left_extend_file)
        self.filepath_list.append(right_extend_file)
        self.filepath_list.append(terminal_file)

        selected_list = candidate_list[:20]
        with open(left_extend_file, 'w') as LF, open(right_extend_file, 'w') as RF, open(terminal_file, 'w') as TF:
            for line in selected_list:
                chrmid, start, stop, name, score, strand = line
                TF.write(''.join(['>', chrmid, '-', str(start), '-', str(stop), strand.replace('-', '-n').replace('+', '-p'), '\n']))
                if int(stop) - int(start) > 200:
                    left_terminal_seq = str(self.genome_dict[chrmid][int(start):int(start)+100])
                    right_terminal_seq = str(self.genome_dict[chrmid][int(stop) -100:int(stop)])
                    TF.write(''.join([left_terminal_seq, 'N'*10, right_terminal_seq, '\n']))
                else:
                    TF.write(''.join([str(self.genome_dict[chrmid][int(start):int(stop)]), '\n']))

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

    def TIR_denection(self, sequencefile):
        sequence_length_dict = SeqIO.parse(sequencefile, 'fasta')
        sequence_length_dict = {k.id: len(k.seq) for k in sequence_length_dict}
        Terminal_dict = defaultdict(lambda :0)
        pairname = os.path.basename(sequencefile).replace('.terminal', '')
        ## start coord is 1 not 0
        dbname = ''.join([os.path.basename(sequencefile), '.invdb'])
        invttirfile = ''.join([os.path.basename(sequencefile), '.inv.txt'])
        ## build database
        mkinvdb = subprocess.Popen(
            ['gt', 'suffixerator', '-db', sequencefile, '-indexname', dbname, '-mirrored', '-dna', '-suf', '-lcp', '-bck'],
            stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        mkinvdb.wait()
        ## run tirvish
        with open(invttirfile, 'w') as invf:
            runinvsearch = subprocess.Popen(['gt', 'tirvish', '-index', dbname, '-mintirlen', '20', '-maxtirlen', '150',
                                             '-similar', '85', '-mintirdist', '2', '-maxtirdist', '200', '-mintsd', '0',
                                             '-seed', '12', '-vic', '1', '-overlaps', 'all', '-xdrop', '0'],
                                            stderr=subprocess.DEVNULL, stdout=invf)
            runinvsearch.wait()

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
                elif splitlines[2] == 'terminal_inverted_repeat':
                    if t == 1:
                        left_start = str(int(splitlines[3]))
                        t += 1
                    else:
                        right_end = str(int(splitlines[4]))
                        if self.intersect([1, int(sequence_length_dict[chrmid])], [left_start, right_end], lportion=0.7):
                            ## means long terminal inverted repeats detected
                            Terminal_dict[chrmid]+=1

        os.remove(invttirfile)
        os.system('rm %s*' % dbname)
        Inv_count = len([i for i in Terminal_dict if Terminal_dict[i]>0])
        Total_num = len(sequence_length_dict.keys())
        self.Terminal_dict[pairname]=Inv_count/Total_num

    def MakeSelection(self, FisherFile_pathlist):
        # This function is to filter out the candidates who might insert into other superfamily of transposons. 
        sys.stdout.write('Begin to run boundary check program.\n')
        subwkdir = 'boundary_align'
        if not os.path.exists(subwkdir):
            os.mkdir(subwkdir)
        else:
            shutil.rmtree(subwkdir)
            os.mkdir(subwkdir)
        ## To output the flanking sequences
        self.filepath_list = []

        planpool = ThreadPool(self.process_num)
        for fisherfile in FisherFile_pathlist:
            planpool.apply_async(self.flanking_seq, args=(fisherfile,))
        planpool.close()
        planpool.join()
        filepath_list = [i for i in self.filepath_list if os.path.getsize(i) and i.endswith('.fa')]
        planpool = ThreadPool(self.process_num)

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

        ## insertions that only contain orf regions will automatically pass this filter

        if os.path.exists(boundary_identity_tbl):
            with open(boundary_identity_tbl, 'r') as F:
                F.readline() ## Skip header
                for line in F:
                    name, direction, iden = line.rstrip().split('\t')
                    self.Boundary_iden_dict[name][direction] = float(iden)

        ## To filter out the candidates who might contain long tirs.
        terminal_file_list = [i for i in self.filepath_list if os.path.getsize(i) and i.endswith('.terminal')]
        planpool = ThreadPool(self.process_num)
        for terminal_file in terminal_file_list:
            planpool.apply_async(self.TIR_denection, args=(terminal_file,))
        planpool.close()
        planpool.join()
        sys.stdout.write("Inverted repeats detection was finished!\n")

    def filter(self, alternative_list, classname, mobile_type):
        ## To do selection for autonomous firstly, because the nonautonomous counterparts will be dertermined by autonomous ones.
        candidate_list = []
        RC_replist = []
        if classname == 'HLE1':
            for line in alternative_list:
                pairname = line[3]
                ## insertions that only contain orf regions will automatically pass this filter, because default value is 0.
                if self.Boundary_iden_dict[pairname]['left'] < self.cutoff_flank and self.Boundary_iden_dict[pairname]['right'] < self.cutoff_flank:
                    candidate_list.append(line)
        else: ## means HLE2
            for line in alternative_list:
                pairname = line[3]
                left_iden, right_iden = self.Boundary_iden_dict[pairname]['left'], self.Boundary_iden_dict[pairname]['right']
                if left_iden < self.cutoff_flank:
                    ## means regular HLE2, try to find stem-loop terminal markers at 5'end
                    candidate_list.append(line)
                elif right_iden < self.cutoff_flank and left_iden >= self.cutoff_flank:
                    ## unregular HLE2, try to find stem-loop terminal markers at 3'end
                    line[-3] = ''.join([line[-3], '.2'])
                    candidate_list.append(line)
        if not candidate_list:
            if mobile_type != 'nonauto':
                ## Might represent a truncated helitron insertion. Just keep the autonomous ones and remove the nonautonomous ones.
                ## Because we are sure that there should be one insertion in autonomous region.
                candidate_list = sorted(alternative_list, key=lambda x: [int(x[2]) - int(x[1]), float(x[6]), -float(x[7])])  ## length (shorter), ## pvalue, then bitscore
        else:
            candidate_list = sorted(candidate_list, key=lambda x: [float(x[6]), int(x[1]) - int(x[2]), -float(x[7])])  ## pvalue, length, then bitscore
        ## To further filter out TIRs
        candidate_list2 = [i for i in candidate_list if self.Terminal_dict[i[3]] <= 0.2]
        if candidate_list2:
            candidate_list = candidate_list2
        else:
            if mobile_type == 'nonauto':
                return []
        if candidate_list:
            RC_replist.append(candidate_list[0])
        return RC_replist

    def OutputSequence(self, bedinput, faoutput):
        with open(bedinput, 'r') as RF, open(faoutput, 'w') as WF:
            for line in RF:
                splitlines = line.rstrip().split('\t')
                chrmid, start, stop, pairname, score, strand = splitlines[:6]
                insertion_name = splitlines[-1]
                start = int(start) - 1
                start = 0 if start <0 else start
                seq = self.genome_dict[chrmid][start:int(stop)]
                if strand == '+':
                    WF.write(''.join(['>', insertion_name, '\n']))
                    WF.write(str(seq))
                    WF.write('\n')
                else:
                    WF.write(''.join(['>', insertion_name, '\n']))
                    WF.write(str(seq.reverse_complement()))
                    WF.write('\n')

    def Convienced_LTS_RTS(self, RCfile, opfile):
        ## To output convinced pair list which can be used as query in another close species genome.
        with open(RCfile, 'r') as F:
            pairlist = {line.rstrip().split("\t")[3] for line in F}
        convienced_dict = defaultdict(list)
        for pair in pairlist:
            tirvalue = self.Terminal_dict[pair]
            if float(tirvalue) <= 0.2:
                if pair in self.Boundary_iden_dict:
                    left_value = float(self.Boundary_iden_dict[pair]['left'])
                    right_value = float(self.Boundary_iden_dict[pair]['right'])
                    classname = pair.split('_left')[0]
                    if 'HLE2' in pair:
                        if left_value < self.cutoff_flank or right_value < self.cutoff_flank:
                            convienced_dict[classname].append(pair)
                    else:  ## means Helitron
                        if left_value < self.cutoff_flank and right_value < self.cutoff_flank:
                            convienced_dict[classname].append(pair)

        oplist = []
        for classname in convienced_dict:
            ## To read the left terminal sequence
            left_terminal_file = ''.join([classname, '/', classname, '.reduce.left.fa'])
            left_terminal_dict = SeqIO.parse(left_terminal_file, 'fasta')
            left_terminal_dict = {k.id: k.seq.upper() for k in left_terminal_dict}
            ## To read the right terminal sequence
            right_terminal_file = ''.join([classname, '/', classname, '.reduce.right.fa'])
            right_terminal_dict = SeqIO.parse(right_terminal_file, 'fasta')
            right_terminal_dict = {k.id: k.seq.upper() for k in right_terminal_dict}

            for pair in convienced_dict[classname]:
                left, right = pair.split('-')
                left_seq = str(left_terminal_dict[left])
                right_seq = str(right_terminal_dict[right])
                oplist.append(''.join([classname, '\t', left, '\t', left_seq, '\t', right,'\t', right_seq, '\n']))
        with open(opfile, 'w') as F:
            F.writelines(oplist)

    def Pre_ts_solo(self, classname, pre_lts_file, pre_rts_file):
        if not os.path.exists(classname):
            os.mkdir(classname)
        os.chdir(classname)
        ORF_bedfile = ''.join([classname, '_orf.bed'])
        ## To evaluate the size range
        distance_na = int(self.window) * 2
        half_distance = int(round(int(distance_na) / 2, 0))
        self.distance_na[classname] = distance_na if not self.distance_na[classname] else self.distance_na[classname]  ## If specified, will use the specified one.
        sys.stdout.write('The length of autonomous %s is expected to be shorter than %s.\n' % (classname, str(distance_na + 100)))
        sys.stdout.write('The length of nonautonomous %s is expected to be shorter than %s.\n' % (classname, str(self.distance_na[classname] + 100)))
        ## To run blastn to get homologies
        left_beddir, right_beddir = 'SubBlastnBed/%s_left' % classname, 'SubBlastnBed/%s_right' % classname
        if not os.path.exists(left_beddir):
            os.makedirs(left_beddir)
        if not os.path.exists(right_beddir):
            os.makedirs(right_beddir)
        self.blastn(pre_lts_file, left_beddir)
        self.blastn(pre_rts_file, right_beddir)
        ## To output combine id
        combinid_file = '%s.combinid.txt' % classname
        with open(combinid_file, 'w') as F:
            F.writelines([''.join([left, '\t', right, '\n']) for left in self.prepair_dict for right in self.prepair_dict[left] if re.match('%s_left'%classname, left)])

        ## Prepare for fisher's exact test (avoid overloading, to split bedfiles )
        sys.stdout.write("Prepare for windowing for %s ...\n" % classname)
        subed_dir = '_'.join([classname, 'Windowing'])
        try:
            split_joint_program = subprocess.Popen(
                ['Rscript', SPLIT_JOINT_PRO, left_beddir, right_beddir, combinid_file, subed_dir,
                 self.genome_size, str(half_distance), BEDTOOLS_PATH, str(self.process_num)], stdout=subprocess.DEVNULL)
            split_joint_program.wait()
            joint_filepath_blastn = ''.join([subed_dir, '/', 'left_right.path.join'])
        except:
            sys.stderr.write("Error: Windowing program failed...\n")
            exit(0)

        ## To run fisher's exact text to select the co-occured left and right terminal signals
        sys.stdout.write("Begin to run fisher's exact test for %s!\n" % classname)

        bed_optdir = '%s_BedFisher' % classname
        if os.path.exists(bed_optdir):
            shutil.rmtree(bed_optdir)
        os.mkdir(bed_optdir)
        Strategy = '1' if Args.nearest else '0'
        try:
            fisher_program = subprocess.Popen(
                ['Rscript', FISHER_PRO, BEDTOOLS_PATH, self.genome_size, joint_filepath_blastn, bed_optdir,
                 str(self.process_num), str(self.pvalue), ORF_bedfile, Strategy], stdout=subprocess.DEVNULL)
            fisher_program.wait()
            sys.stdout.write("Fisher's exact test finished for %s!\n" % classname)
        except:
            sys.stderr.write("Error: Fisher's exact test failed...\n")
            exit(0)

        FisherBed_files = ['/'.join([bed_optdir, file]) for file in os.listdir(bed_optdir)]

        ## To merge and sort all Fisher bed files
        fisher_pvalue_file = '%s.joint.pvalue.bed' % classname
        fisher_filelist = [''.join([bed_optdir, '/', file]) for file in os.listdir(bed_optdir) if file.endswith('.bed')]
        with open(fisher_pvalue_file, 'w') as WF:
            for subfisher in fisher_filelist:
                with open(subfisher, 'r') as F:
                    content = F.read()
                    WF.write(content)
        try:
            sort_pro = subprocess.Popen(
                ['bash', SORT_PRO, fisher_pvalue_file, str(self.process_num)],
                stdout=subprocess.DEVNULL)
            sort_pro.wait()
            sys.stdout.write("Merge and sort fisher result files finished for %s!\n" % classname)
        except:
            sys.stdout.write("Merge and sort fisher result files failed for %s!\n" % classname)
            exit(0)
        fisher_pvalue_bed = BT.BedTool(fisher_pvalue_file).saveas('Nonauto.bed')
        os.chdir('../')
        
    def main(self):
        self.Boundary_iden_dict = defaultdict(lambda: defaultdict(lambda: 1))
        self.Terminal_dict = defaultdict(lambda: 1)

        RC_total_candidate_dict = defaultdict(list)
        RC_total_candidate = self.autonomous_detect()
        [RC_total_candidate_dict[line[4]].append(line) for line in RC_total_candidate]
        orf_classname_list = list(RC_total_candidate_dict.keys())
        unique_pre_ts_classname = []
        classname_list = []
        classname_list.extend(orf_classname_list)
        if self.prepair_dict:
            ## means add pre-ts
            unique_pre_ts_classname = list(set(self.terminalfile_dict.keys()) - set(orf_classname_list))
            classname_list.extend(unique_pre_ts_classname)

        Repsentative_file_list = []
        Alternative_file_list = []
        ## Loop for each variants.
        for class_name in classname_list:
            op_representative = '%s.representative.bed' % class_name
            opauto_nest_alt = '%s.auto.nest.alt.bed' % class_name
            opnonauto_nest_alt = '%s.nonauto.nest.alt.bed' % class_name
            op_nest_rep = '%s.nest.bed' % class_name

            if os.path.exists(op_representative):
                os.remove(op_representative)
            Repsentative_file_list.append(op_representative)
            if os.path.exists(op_nest_rep):
                os.remove(op_nest_rep)
            Alternative_file_list.append(op_nest_rep)

            FisherFile_list = []
            Auto_file_list = []
            ORFonly_file_list = []
            RC_list = RC_total_candidate_dict[class_name]
            if RC_list:
                if re.findall('HLE2', class_name):  ## need to pair the terminal repeats
                    self.prepare_terminal_seq(RC_list, pair=True, classname=class_name)
                else:
                    if not Args.pair_helitron:
                        self.prepare_terminal_seq(RC_list, pair=False, classname=class_name)
                    else:
                        self.prepare_terminal_seq(RC_list, pair=True, classname=class_name)
            else:
                self.Pre_ts_solo(classname=class_name, pre_lts_file=self.terminalfile_dict[class_name]['left'],
                                 pre_rts_file=self.terminalfile_dict[class_name]['right'])

            fisher_beddir = ''.join([class_name, '/', class_name, '_BedFisher'])
            if os.path.exists(fisher_beddir):
                FisherFile_list.extend([''.join([fisher_beddir, '/', file]) for file in os.listdir(fisher_beddir)])
                ## To check that if the boundary of each family can be well aligned.
                ## To check that if the families contain long terminal inverted repeats.
                self.MakeSelection(FisherFile_list)

            ## To merge autonomous insertions and select representatives
            pairnamelist = set()
            if RC_list:
                ORF_bed = BT.BedTool(''.join([class_name, '/', class_name, '_orf.bed']))
                Auto_fisher_file = ''.join([class_name, '/', '%s.auto.alternative.bed' % class_name])

                if os.path.isfile(Auto_fisher_file):
                    Auto_fisher_bed = BT.BedTool(Auto_fisher_file)
                else:
                    Auto_fisher_bed = BT.BedTool([])
                Auto_ORF_intersect_file = ''.join([class_name, '/', class_name, 'ORFinterAUTO.bed'])
                ORF_bed.intersect(Auto_fisher_bed, nonamecheck=True, f=1, wo=True, s=True).saveas(Auto_ORF_intersect_file)

                if os.path.exists(Auto_fisher_file):
                    pairnamelist = self.merge_overlaped_autos(Auto_ORF_intersect_file, op_representative, opauto_nest_alt, rep_type=True)
                    ## To output alternative nest insertions
                    if os.path.exists(opauto_nest_alt):
                        self.merge_overlaped_autos(opauto_nest_alt, op_nest_rep, 'tmp.txt', rep_type=False)

            ## To get non-autonomous candidates
            nonauto_bedfile = '%s/Nonauto.bed' % class_name
            nonauto_alternative_file = ''.join([class_name, '/', class_name, '.nonauto.alternative.bed'])
            if os.path.exists(nonauto_bedfile):
                with open(nonauto_alternative_file, 'w') as WF, open(nonauto_bedfile, 'r') as RF:
                    for line in RF:
                        splitlines = line.rstrip().split('\t')
                        ## To filter out ultra-large nonautonomous
                        if int(splitlines[2]) - int(splitlines[1]) + 1 > self.distance_na[class_name] + 100:
                            continue
                        ## To limit outputing nonautonomous candidates who shares the same autonomous boundaries.
                        if Args.multi_ts:
                            newline = '\t'.join(
                                [splitlines[0], splitlines[1], splitlines[2], splitlines[3], splitlines[4],
                                 splitlines[5], splitlines[6], splitlines[7], class_name, 'nonauto-%s'%splitlines[8]])
                            WF.write(newline)
                            WF.write('\n')
                        else:
                            if splitlines[3] in pairnamelist or 'pre' in splitlines[3]:
                                newline = '\t'.join([splitlines[0], splitlines[1], splitlines[2], splitlines[3], splitlines[4],
                                                    splitlines[5], splitlines[6], splitlines[7], class_name, 'nonauto-%s'%splitlines[8]])
                                WF.write(newline)
                                WF.write('\n')

            ## To merge nonautonomous intervals and do selection.
            if os.path.exists(nonauto_alternative_file):
                self.merge_overlaped_intervals(nonauto_alternative_file, op_representative, opnonauto_nest_alt, rep_type=True)
                ## To output alternative nest insertions
                if os.path.exists(opnonauto_nest_alt):
                    self.merge_overlaped_intervals(opnonauto_nest_alt, op_nest_rep, 'tmp.txt', rep_type=False)
            ## To integrate orf file into representative file
            # To obtain orf only insertions
            Orfonly_file = ''.join([class_name, '/', '%s.orfonly.alternative.bed' % class_name])
            if os.path.exists(Orfonly_file):
                with open(op_representative, 'a') as WF, open(Orfonly_file, 'r') as RF:
                    WF.write(RF.read())
            ## Remove intermediate files
            if os.path.exists(opauto_nest_alt):
                os.remove(opauto_nest_alt)
            if os.path.exists(opnonauto_nest_alt):
                os.remove(opnonauto_nest_alt)

        ## To output boundary check file and tir check file
        with open('Boundary.tbl', 'w') as BF:
            for id in self.Boundary_iden_dict:
                for direction in self.Boundary_iden_dict[id]:
                    BF.write(''.join([id, '\t', direction, '\t', str(self.Boundary_iden_dict[id][direction]), '\n']))
        with open('TIR_count.tbl', 'w') as TF:
            for id in self.Terminal_dict:
                TF.write(''.join([id, '\t', str(self.Terminal_dict[id]), '\n']))

        ## To add terminal markers for HLE2
        for helen_file in Repsentative_file_list:
            if re.findall('HLE2', helen_file):
                self.heltentron_terminal(helentron_bed=helen_file)

        if os.path.exists('tmp.txt'):
            os.remove('tmp.txt')

        op_representative = 'RC.representative.bed'
        with open(op_representative, 'w') as wf:
            for file in Repsentative_file_list:
                if os.path.exists(file):
                    with open(file, 'r') as rf:
                        [wf.write(line) for line in rf]
                    os.remove(file)

        op_alternative = 'RC.alternative.bed'
        with open(op_alternative, 'w') as wf:
            for file in Alternative_file_list:
                if os.path.exists(file):
                    with open(file, 'r') as rf:
                        [wf.write(line) for line in rf]
                    os.remove(file)

        ## To sort big file
        if not os.path.exists(op_representative):
            with open(op_representative, 'w') as F:
                F.write('')
        if not os.path.exists(op_alternative):
            with open(op_alternative, 'w') as F:
                F.write('')
        try:
            sort_pro = subprocess.Popen(['bash', SORT_PRO, op_representative, str(self.process_num)],
                                        stdout=subprocess.DEVNULL)
            sort_pro.wait()
            sys.stdout.write("Sort representative files finished ...\n")

            sort_pro = subprocess.Popen(['bash', SORT_PRO, op_alternative, str(self.process_num)],
                                        stdout=subprocess.DEVNULL)
            sort_pro.wait()
            sys.stdout.write("Sort alternative files finished ...\n")

        except:
            sys.stdout.write("Sort representative files failed ...\n")
            exit(0)
        sequence_fa = 'RC.representative.fa'
        ## Delete alternative file if empty
        if not os.path.getsize(op_alternative):
            os.remove(op_alternative)
        ## To output fasta format file.
        self.OutputSequence(op_representative, sequence_fa)
        ## To output convinced pair list.
        Convinced_pairlist = 'pairlist.tbl'
        self.Convienced_LTS_RTS(op_representative, Convinced_pairlist)

        ## To remove tempary files
        os.remove(self.CTRR_stem_loop_description)
        os.remove(self.subtir_description)
        os.remove(self.genome_size)
        shutil.rmtree(self.bedtoolstmp)
        if os.path.exists('GenomeDB'):
            shutil.rmtree('GenomeDB')
        if os.path.exists('genomes'):
            shutil.rmtree('genomes')
        #for class_name in RC_total_candidate_dict:
        #    if os.path.exists(class_name):
        #       shutil.rmtree(class_name)
        if os.path.exists("boundary_align"):
            shutil.rmtree("boundary_align")
        if os.path.exists('Boundary.identity.tbl'):
            os.remove('Boundary.identity.tbl')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="heliano can detect and classify different variants of Helitron-like elements: HLE1 and HLE2. Please visit https://github.com/Zhenlisme/heliano/ for more information. Email us: zhen.li3@universite-paris-saclay.fr")
    parser.add_argument("-g", "--genome", type=str, required=True, help="The genome file in fasta format.")
    parser.add_argument("-w", "--window", type=int, default=10000, required=False,
                        help="To check terminal signals within a given window bp upstream and downstream of ORF ends, default is 10 kb.")
    parser.add_argument("-dm", "--distance_domain", type=int, default=2500, required=False,
                        help="The distance between HUH and Helicase domain, default is 2500.")
    parser.add_argument("-dn", "--distance_ts", type=int, default=0, required=False,
                        help="The maximum distance between LTS and RTS. If not specified, HELIANO will set it as two times window size plus the maximum ORF length.")
    parser.add_argument("-pt", "--pair_helitron", type=int, default=1, required=False, choices=[0, 1],
                        help="For HLE1, its 5' and 3' terminal signal pairs should come from the same autonomous helitorn or not. 0: no, 1: yes. default yes.")
    parser.add_argument("-is1", "--IS1", type=int, default=0, required=False, choices=[0, 1],
                        help="Set the insertion site of autonomous HLE1 as A and T. 0: no, 1: yes. default yes.")
    parser.add_argument("-is2", "--IS2", type=int, default=0, required=False, choices=[0, 1],
                        help="Set the insertion site of autonomous HLE2 as T and T. 0: no, 1: yes. default yes.")
    parser.add_argument("-sim_tir", "--simtir", type=int, default=100, required=False, choices=[100, 90, 80],
                        help="Set the simarity between short inverted repeats(TIRs) of HLE2. Default 100.")
    parser.add_argument("-flank_sim", "--flank_sim", type=float, default=0.5, required=False, choices=[0.4, 0.5, 0.6, 0.7],
                        help="The cut-off to define false positive LTS/RTS. The lower the value, the more strigent. Default 0.5.")
    parser.add_argument("-p", "--pvalue", type=float, required=False, default=1e-5, help="The p-value for fisher's exact test. default is 1e-5.")
    parser.add_argument("-s", "--score", type=int, required=False, default=32,
                        help="The minimum bitscore of blastn for searching for homologous sequences of terminal signals. From 30 to 55, default is 32.")
    parser.add_argument("--nearest", action='store_true', required=False,
                        help="If you use this parameter, you will use the reciprocal-nearest LTS-RTS pairs as final candidates. By default, HELIANO will try to use the reciprocal-farthest pairs.")
    parser.add_argument("-ts", "--terminal_sequence", type=str, required=False, default='', help="The terminal sequence file. You can find it in the output of previous run (named as pairlist.tbl).")
    parser.add_argument("--dis_denovo", action='store_true', required=False,
                        help="If you use this parameter, you refuse to search for LTS/RTS de novo, instead you will only use the LTS/RTS information described in the terminal sequence file.")
    parser.add_argument("--multi_ts", action='store_true', required=False,
                        help="To allow an auto HLE to have multiple terminal sequences. If you enable this, you might find nonauto HLEs coming from the same auto HLE have different terminal sequences.")
    parser.add_argument("-o", "--opdir", type=str, required=True, help="The output directory.")
    parser.add_argument("-n", "--process", type=int, default=2, required=False, help="Maximum number of threads to be used.")
    parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.2.1')
    Args = parser.parse_args()
    if int(Args.score) < 30 or int(Args.score) >= 55:
        sys.stderr.write("Error: The bitscore value should be greater than 30 and less than 55.\n")
        exit(0)
    if int(Args.distance_ts) < 0 or int(Args.window) < 0 or int(Args.distance_domain) < 0:
        sys.stderr.write("Error: Parameter value should not be negative.\n")
        exit(0)
    ## To set and check dependency file path
    HMMFILE = '_HMM_'
    HEADERFILE = '_HEADER_'
    FISHER_PRO = '_FISHER_'
    BOUNDARY_PRO = '_BOUNDARY_'
    SPLIT_JOINT_PRO = '_SPLIT_JOINT_'
    SORT_PRO = '_SORTPRO_'
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
    if not  os.path.exists(SPLIT_JOINT_PRO):
        print("Split bed file program not found!")
        exit(0)

    HomoSearch = Homologous_search(HMMFILE, os.path.abspath(Args.genome), os.path.abspath(Args.opdir),
                                   HEADERFILE, Args.window, Args.distance_domain, Args.distance_ts, Args.pvalue,
                                   Args.process)
    HomoSearch.main()

