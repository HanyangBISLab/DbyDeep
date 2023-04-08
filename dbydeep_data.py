import time
import pandas as pd
import json
import re
import argparse
from ahocorapy.keywordtree import KeywordTree
import warnings
warnings.filterwarnings(action='ignore')


class preprocessor():
    def __init__(self):
        pass

    def get_protein(self, data_path: str):
        '''
        protein database to protein sequence dictionary
        Input : fasta file
        Output : { protein : sequence } dict
        '''
        prot_seq = dict()  # init
        seq = ''
        NAME_IDX = 0
        FIRST_LINE = True
        with open(data_path) as f:
            while f:
                l = f.readline().replace('\n', '')
                if l == '':
                    prot_seq[prot_name] = seq
                    break

                # get protein name and sequence
                if '>' in l:
                    if FIRST_LINE:
                        prot_name = l.replace('>','').split(' ')[NAME_IDX]
                        FIRST_LINE = False
                        continue
                    prot_seq[prot_name] = seq
                    # new init for next protein
                    prot_name = l.replace('>','').split(' ')[NAME_IDX]
                    seq = ''
                else:
                    seq += l
        self.prot_seq = prot_seq

    def get_precursor(self, data_path: str, tool_name:str ='mgf'):
        '''
        database search result to identified proteins and precursors dictionary
        Input : database search result
        Output : prot_precursor { protein : [ (peptide, charge)s ] } dict
        '''
        if tool_name == 'comet':
            pattern = r'\[[^]]*\]'
            prot_precursor = dict()
            with open(data_path) as f:
                for idx, line in enumerate(f):
                    if idx == 0:
                        header = line.replace('\n', '').split('\t')
                    else:
                        line = line.replace('\n', '').split('\t')
                        PSMId = line[0]
                        score = line[1]
                        q_value = line[2]
                        posterior_error_prob = line[3]
                        peptide = line[4]
                        proteinIds = line[5:]  # use all protein ids

                        if float(q_value) <= 0.01:
                            peptide = re.sub(pattern=pattern, repl='', string=peptide).split('.')[1]  # except PTM
                            peptide = ''.join([_ for _ in peptide if _.isupper()])
                            precursor = (peptide, 2)  # charge 2 fix
                            for pep_prot in proteinIds:
                                if 'DECOY_' not in pep_prot:  # exclude decoy protein
                                    if pep_prot not in prot_precursor:
                                        prot_precursor[pep_prot] = set()
                                    prot_precursor[pep_prot].add(precursor)
            self.prot_precursor = prot_precursor
        
        elif tool_name == 'msgfplus':
            pattern = r'\([^)]*\)'
            prot_precursor = dict()
            df=pd.read_csv(data_path, sep='\t')
            peptides = df.Peptide.values
            proteins = df.Protein.values
            for peptide, proteinIds in zip(peptides, proteins):
                proteinIds=proteinIds.split(';')
                peptide = ''.join([_ for _ in peptide if ord(_) in range(65, 91)])  # except PTM
                peptide = ''.join([_ for _ in peptide if _.isupper()])
                precursor = (peptide, 2)  # charge 2 fix
                for pep_prot in proteinIds:
                    if 'XXX' not in pep_prot:  # exclude decoy protein
                        pep_prot = re.sub(pattern=pattern, repl='', string=pep_prot)
                        if pep_prot not in prot_precursor:
                            prot_precursor[pep_prot] = set()
                        prot_precursor[pep_prot].add(precursor)
            self.prot_precursor = prot_precursor

    def merge_precursor_cnt(self):
        '''
        Merge Human Spectral Library(HSL) Candidate to HSL for spectral count per precursor
        Input : protein2precursor, precursor2count
        Output : protein2precursor_cnt { protein : [ (peptide, charge, spectral_cnt)s ] }
        '''
        prot_precursor_cnt = dict()
        for prot_name, precursors in self.prot_precursor.items():
            precursors_cnt = [(*precursor, 2) for precursor in precursors]  # all peptide has spectral count of 2 (for deploy version)
            prot_precursor_cnt[prot_name] = precursors_cnt
        self.prot_precursor = prot_precursor_cnt

    def merge_protein_precursor(self):
        '''
        Merge Uniprot and Human Spectral Library(HSL) for identified peptide counting
        Input : prot_precursor_cnt, prot_seq
        output : prot_seq_precursor { protein : {sequence: ARNDCEQG, 
                                                precursors: [ (peptide, charge, spectral_cnt)s ]
                                                } 
                                    }
        '''
        identified_prot = set(self.prot_precursor.keys())
        uniprot_prot = set(self.prot_seq.keys())
        inter_prot_names = identified_prot.intersection(uniprot_prot)
        prot_seq_precursor = dict()
        for prot_name in inter_prot_names:
            prot_seq_precursor[prot_name] = {'sequence':self.prot_seq[prot_name],
                                            'precursors':self.prot_precursor[prot_name]}
        self.prot_seq_precursor = prot_seq_precursor

    def get_peptide_cnt(self):
        '''
        Counting of striped peptide's spectral count (sum)
        Input : prot_seq_precursor
        Output : prot_seq_peptide { protein : {sequence: ARNDCEQG, 
                                            peptides: { peptide : count }
                                            } 
                                    }
        '''
        prot_seq_peptide = dict()
        for prot_name, prot_info in self.prot_seq_precursor.items():
            prot_seq = prot_info['sequence']
            prot_precursors = prot_info['precursors']

            prot_pep_cnt = dict()
            for precursor in prot_precursors:
                pre_ptm_seq, pre_charge, pre_cnt = precursor
                pre_strip_seq = ''.join(filter(lambda x: ord(x) in range(65, 91), pre_ptm_seq))
                if pre_strip_seq not in prot_pep_cnt:
                    prot_pep_cnt[pre_strip_seq] = 0
                prot_pep_cnt[pre_strip_seq] += pre_cnt
            if prot_name not in prot_seq_peptide:
                prot_seq_peptide[prot_name] = dict()
            prot_seq_peptide[prot_name]['sequence'] = prot_seq
            prot_seq_peptide[prot_name]['peptides'] = prot_pep_cnt
        
        self.prot_seq_peptide = prot_seq_peptide

    def get_peptide_tree(self):
        '''
        Contructing peptide tree
        Input : prot_seq_peptide
        Output : prot_seq_peptide { protein : {sequence: ARNDCEQG, 
                                            peptides: { peptide : count }
                                            peptide_tree: tree_instance
                                            } 
                                    }
        '''
        for prot_name, prot_info in self.prot_seq_peptide.items():
            prot_peptideas = prot_info['peptides']
            tree = KeywordTree()
            for peptide in prot_peptideas:
                tree.add(peptide)
            tree.finalize()
            self.prot_seq_peptide[prot_name]['peptide_tree'] = tree

    def get_peptide_from_protein(self, DIGEST_MERS:int = 7):
        '''
        Getting peptides from proteins
        Input : prot_seq_peptide
        Output : prot_seq_peptide { protein : {sequence: ARNDCEQG, 
                                            cleavage_cnt: [1, 0, ..., 4],
                                            miss_cleavage_cnt: [0, 1, ..., 0],
                                            peptides: { peptide : count }
                                            peptide_tree: tree_instance,
                                            coverage: 0.5,
                                            labelled_peptides : { peptide : {n: AA,
                                                                                c: AA,
                                                                                m1: AA,
                                                                                m2: AA
                                                                    }
                                            }
                                    }
        '''
        
        TRYPTIC_SITE = 'KR'
        for prot_name, prot_info in self.prot_seq_peptide.items():
            # for labelling
            pep_tree = prot_info['peptide_tree']
            pep_cnt = prot_info['peptides']

            # slicing
            prot_seq = prot_info['sequence']
            self.prot_seq_peptide[prot_name]['labelled_peptides'] = dict()  # init
            cleavage_sites = [AA_IDX for AA_IDX, aa in enumerate(prot_seq) if aa in TRYPTIC_SITE]
            MISS_CLEAVAGES = [0, 1, 2]
            for MISS_CLEAVAGE_CNT in MISS_CLEAVAGES:

                # consider miss cleavage sites and c-terminal site
                cleavage_range = range(len(cleavage_sites) - MISS_CLEAVAGE_CNT - 1)  
                for CS_IDX in cleavage_range:
                    

                    # N terminal of peptide
                    condi_prot_nterm_pep_nterm = cleavage_sites[CS_IDX] < DIGEST_MERS
                    condi_prot_cterm_pep_nterm = cleavage_sites[CS_IDX] > len(prot_seq) - 1 - (DIGEST_MERS)
                    nterm_pad_pep_nterm = 'Z' * (DIGEST_MERS - cleavage_sites[CS_IDX]) 
                    cterm_pad_pep_nterm = 'Z' * (DIGEST_MERS - (len(prot_seq) - 1 - cleavage_sites[CS_IDX]))
                    pep_nterm_start_idx = cleavage_sites[CS_IDX] - DIGEST_MERS
                    pep_nterm_end_idx = cleavage_sites[CS_IDX] + DIGEST_MERS + 1
                    if condi_prot_nterm_pep_nterm and condi_prot_cterm_pep_nterm:  # ex) --MNQKLLK-- : both n and c term of protein are insufficient
                        # pep_seq = prot_seq[:pep_nterm_end_idx]  # test ...
                        pep_seq = prot_seq
                        pep_nterm = nterm_pad_pep_nterm + pep_seq + cterm_pad_pep_nterm
                    elif condi_prot_nterm_pep_nterm:
                        pep_seq = prot_seq[:pep_nterm_end_idx]
                        pep_nterm = nterm_pad_pep_nterm + pep_seq
                    elif condi_prot_cterm_pep_nterm:
                        pep_seq = prot_seq[pep_nterm_start_idx:]
                        pep_nterm = pep_seq + cterm_pad_pep_nterm
                    else:
                        pep_nterm = prot_seq[pep_nterm_start_idx : pep_nterm_end_idx]
                    
                    # C terminal of peptide
                    condi_prot_nterm_pep_cterm = cleavage_sites[CS_IDX + MISS_CLEAVAGE_CNT + 1] < DIGEST_MERS
                    condi_prot_cterm_pep_cterm = cleavage_sites[CS_IDX + MISS_CLEAVAGE_CNT + 1] > len(prot_seq) - 1 - (DIGEST_MERS)
                    nterm_pad_pep_cterm = 'Z' * (DIGEST_MERS - cleavage_sites[CS_IDX + MISS_CLEAVAGE_CNT + 1])
                    cterm_pad_pep_cterm = 'Z' * (DIGEST_MERS - (len(prot_seq) - 1 - cleavage_sites[CS_IDX + MISS_CLEAVAGE_CNT + 1]))
                    pep_cterm_start_idx = cleavage_sites[CS_IDX + MISS_CLEAVAGE_CNT + 1] - DIGEST_MERS
                    pep_cterm_end_idx = cleavage_sites[CS_IDX + MISS_CLEAVAGE_CNT + 1] + (DIGEST_MERS + 1)
                    if condi_prot_nterm_pep_cterm and condi_prot_cterm_pep_cterm:
                        # pep_seq = prot_seq[:pep_cterm_end_idx]  # test ...
                        pep_seq = prot_seq
                        pep_cterm = nterm_pad_pep_cterm + pep_seq + cterm_pad_pep_cterm
                    elif condi_prot_nterm_pep_cterm:
                        pep_seq = prot_seq[:pep_cterm_end_idx]
                        pep_cterm = nterm_pad_pep_cterm + pep_seq
                    elif condi_prot_cterm_pep_cterm:
                        pep_seq = prot_seq[pep_cterm_start_idx:]
                        pep_cterm = pep_seq + cterm_pad_pep_cterm
                    else:
                        pep_cterm = prot_seq[pep_cterm_start_idx : pep_cterm_end_idx]

                    # missed cleavage site of peptide
                    pep_miss = {'pep_miss1': 'Z' * DIGEST_MERS + 'Z' + 'Z' * DIGEST_MERS,  # init
                                'pep_miss2': 'Z' * DIGEST_MERS + 'Z' + 'Z' * DIGEST_MERS}
                    for mcc in range(1, MISS_CLEAVAGE_CNT + 1):
                        condi_prot_nterm_pep_miss = cleavage_sites[CS_IDX + mcc] < DIGEST_MERS
                        condi_prot_cterm_pep_miss = cleavage_sites[CS_IDX + mcc] > len(prot_seq) - 1 - (DIGEST_MERS)
                        nterm_pad_pep_miss = 'Z' * (DIGEST_MERS - cleavage_sites[CS_IDX + mcc])
                        cterm_pad_pep_miss = 'Z' * (DIGEST_MERS - (len(prot_seq) - 1 - cleavage_sites[CS_IDX + mcc]))
                        pep_miss_start_idx = cleavage_sites[CS_IDX + mcc] - DIGEST_MERS
                        pep_miss_end_idx = cleavage_sites[CS_IDX + mcc] + DIGEST_MERS + 1
                        if condi_prot_nterm_pep_miss and condi_prot_cterm_pep_miss:
                            pep_seq = prot_seq
                            pep_miss['pep_miss' + str(mcc)] = nterm_pad_pep_miss + pep_seq + cterm_pad_pep_miss
                        if condi_prot_nterm_pep_miss:
                            pep_seq = prot_seq[:pep_miss_end_idx]
                            pep_miss['pep_miss' + str(mcc)] = nterm_pad_pep_miss + pep_seq
                        elif condi_prot_cterm_pep_miss:
                            pep_seq = prot_seq[pep_miss_start_idx:]
                            pep_miss['pep_miss' + str(mcc)] = pep_seq + cterm_pad_pep_miss
                        else:
                            pep_miss['pep_miss' + str(mcc)] = prot_seq[pep_miss_start_idx : pep_miss_end_idx]
                    pep_miss1, pep_miss2 = pep_miss['pep_miss1'], pep_miss['pep_miss2']

                    # peptide
                    pep_body_start_idx = cleavage_sites[CS_IDX] + 1
                    pep_body_end_idx = cleavage_sites[CS_IDX + MISS_CLEAVAGE_CNT + 1] + 1
                    pep_body = prot_seq[pep_body_start_idx : pep_body_end_idx]
                    if not 6 <= len(pep_body) <= 40:  # following identified peptide (massIVE-KB) length distribution
                        continue  # not save this peptide

                    # filtration of U, X (Amino Acid)
                    ux_filter = lambda seq: True if 'U' in seq or 'X' in seq else False
                    if ux_filter(pep_body) or ux_filter(pep_nterm) or ux_filter(pep_cterm) or ux_filter(pep_miss1) or ux_filter(pep_miss2):
                        continue  # not save this peptide

                    # labelling
                    PEP_IDX = 0
                    PEP_START_IDX = 1
                    identified_peps = [p[PEP_IDX] for p in pep_tree.search_all(pep_body)]
                    if pep_body in identified_peps:
                        if pep_cnt[pep_body] > 1:  # positive : peptide has spectral count over 1 (at least 2)
                            label = True
                        else:
                            continue  # not save this peptide (neither positive nor negative peptide)
                    else:  # negative : not identified peptide
                        label = False
                    
                    # save
                    pep_set = ' '.join([pep_body, pep_nterm, pep_cterm, pep_miss1, pep_miss2])
                    self.prot_seq_peptide[prot_name]['labelled_peptides'][pep_set] = label

    def divide_digestability_detectability(self):
        '''
        dividing peptides for compare with AP3
        Input : prot_seq_peptide
        Output : prot_digest, prot_detect (same with input, just divide)
        '''
        
        prot_name_idx = {pn:idx for idx, pn in enumerate(self.prot_seq_peptide.keys())}
        prot_name_idx_rev = {idx:pn for pn, idx in prot_name_idx.items()}
        digest_prot_name = {prot_name_idx_rev[digest_idx] for digest_idx in range(len(prot_name_idx_rev)//2)}
        detect_prot_name = {prot_name_idx_rev[detect_idx] for detect_idx in range(len(prot_name_idx_rev)//2, len(prot_name_idx_rev))}
        prot_digest = {pn:dict(self.prot_seq_peptide[pn]) for pn in digest_prot_name}
        prot_detect = {pn:dict(self.prot_seq_peptide[pn]) for pn in detect_prot_name}
        self.prot_digest = prot_digest
        self.prot_detect = prot_detect
           
    def merge_digestability_detectability(self):
        '''
        Merging peptide to dataframe
        Input : prot_seq_peptide
        Output : pandas dataframe (peptide, nterm, cterm, miss1, miss2, label columns)
        '''

        for_detect_peptide_dupli = dict()
        for prot_name, prot_info in self.prot_detect.items():
            for pep_set, label in prot_info['labelled_peptides'].items():
                p, n, c, m1, m2 = pep_set.split()
                pep_set = (p, n, c, m1, m2)
                if pep_set not in for_detect_peptide_dupli:
                    for_detect_peptide_dupli[pep_set] = set()
                for_detect_peptide_dupli[pep_set].add(label)
        for_detect_peptide = set()
        for pep_set, labels in for_detect_peptide_dupli.items():
            label = max(labels)
            p, n, c, m1, m2 = pep_set
            for_detect_peptide.add((p, n, c, m1, m2, label))

        for_digest_peptide_dupli = dict()
        for prot_name, prot_info in self.prot_digest.items():
            for pep_set, label in prot_info['labelled_peptides'].items():
                p, n, c, m1, m2 = pep_set.split()
                pep_set = (p, n, c, m1, m2)
                if pep_set not in for_digest_peptide_dupli:
                    for_digest_peptide_dupli[pep_set] = set()
                for_digest_peptide_dupli[pep_set].add(label)
        for_digest_peptide = set()
        for pep_set, labels in for_digest_peptide_dupli.items():
            label = max(labels)
            p, n, c, m1, m2 = pep_set
            for_digest_peptide.add((p, n, c, m1, m2, label))

        df_detect = pd.DataFrame(for_detect_peptide,
                                columns = ['peptide', 'nterm', 'cterm', 'miss1', 'miss2', 'label'])
        df_digest = pd.DataFrame(for_digest_peptide,
                                columns = ['peptide', 'nterm', 'cterm', 'miss1', 'miss2', 'label'])
        
        # remove duplicated peptide context
        df_detect.drop_duplicates(inplace=True, ignore_index=True)
        df_digest.drop_duplicates(inplace=True, ignore_index=True)
        df_data = pd.concat([df_detect, df_digest], axis=0).reset_index(drop=True)
        df_data_key = df_data.drop(['label'], axis=1)
        detect_idx = set(df_data.iloc[:len(df_detect)].index)
        digest_idx = set(df_data.iloc[len(df_detect):].index)
        df_data_key.drop_duplicates(inplace=True, ignore_index=False)
        remain_idx = set(df_data_key.index)
        detect_idx = remain_idx.intersection(detect_idx)
        digest_idx = remain_idx.intersection(digest_idx)
        df_detect = df_data.loc[detect_idx]
        df_digest = df_data.loc[digest_idx]

        # divide train test dataset
        df_test = df_detect.sample(frac = 0.4, random_state = 2022)
        df_train_detect = df_detect.drop(df_test.index, axis=0)
        df_train = pd.concat([df_train_detect, df_digest], axis=0).reset_index(drop=True)

        # remove shared peptide in train dataset
        df_key = pd.DataFrame([[p, True] for p in df_test['peptide'].unique()],
                            columns = ['peptide', 'drop'])
        df_train = df_train.merge(df_key, on='peptide', how='left')
        df_train = df_train.loc[df_train['drop'].isnull()].drop('drop', axis=1).reset_index(drop=True)

        # sampling of negative data
        positive_tmp = df_test.loc[df_test.label==True]
        negative_tmp = df_test.loc[df_test.label==False].sample(len(positive_tmp), random_state=2022)
        df_test_sampled = pd.concat([positive_tmp, negative_tmp] ,axis=0).reset_index(drop=True)

        positive_tmp = df_train.loc[df_train.label==True]
        negative_tmp = df_train.loc[df_train.label==False].sample(len(positive_tmp), random_state=2022)
        df_train_sampled = pd.concat([positive_tmp, negative_tmp] ,axis=0).reset_index(drop=True)
        return df_train_sampled, df_test_sampled, df_train, df_test


def main(save_path, protein_path, peptide_path, tool_name):

    p = preprocessor()
    start_time = time.time()

    # preprocessing of Uniprot
    print('### get protein ... \t\t\t', end='\r')
    p.get_protein(protein_path)
    elapsed_time = int((time.time()-start_time)/60)
    print(f'### finish getting of proteins ... time elapsed {elapsed_time} min \t\t\t')

    # preprocessing of peptide (massIVE-KB)
    print('### get spectral count ...\t\t\t', end='\r')
    p.get_precursor(peptide_path, tool_name)
    p.merge_precursor_cnt()
    elapsed_time = int((time.time()-start_time)/60)
    print(f'### finish getting of precursors ... time elapsed {elapsed_time} min \t\t\t')

    # preprocessing of massIVE-KB and Uniprot
    print('### merge files ... \t\t\t', end='\r')
    p.merge_protein_precursor()
    elapsed_time = int((time.time()-start_time)/60)
    print(f'### finish merging of proteins and precursors ... time elapsed {elapsed_time} min \t\t\t')

    # preprocessing of dataset for train
    print('### matching proteins and peptides ... \t\t\t', end='\r')
    p.get_peptide_cnt()
    p.get_peptide_tree()        
    elapsed_time = int((time.time()-start_time)/60)
    print(f'### finish matching proteins and peptides ... time elapsed {elapsed_time} min \t\t\t')

    # Counting spectral count
    print('### Counting of peptide\'s spectral experiments ... \t\t\t', end='\r')
    elapsed_time = int((time.time()-start_time)/60)
    print(f'### finish counting and Labelling ... time elapsed {elapsed_time} min \t\t\t')

    # Getting peptides from protein
    print('### Getting peptides from proteins ... \t\t\t', end='\r')
    p.get_peptide_from_protein()
    elapsed_time = int((time.time()-start_time)/60)
    print(f'### finish getting and labelling of peptides ... time elapsed {elapsed_time} min \t\t\t')

    print('### Dividing and Merging for train, val, test ... \t\t\t', end='\r')
    p.divide_digestability_detectability()
    df_train, df_test, df_train_whole, df_test_whole = p.merge_digestability_detectability()
    elapsed_time = int((time.time()-start_time)/60)
    print(f'### finish merging peptide to train, val, and test dataset ... time elapsed {elapsed_time} min \t\t\t')
    
    pd.concat([df_train_whole, df_test_whole]).to_csv(save_path+'data.csv', index=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--save-path', default='./data/', type=str, help='save folder path')
    # protein DB
    parser.add_argument('--protein-fasta', type=str, help='fasta. human protein (uniprot) file path')
    # peptide search result for test
    parser.add_argument('--peptide-tsv', type=str, help='tsv. peptide search result file path')
    parser.add_argument('--tool-name', type=str, help='tool name')
    opt = parser.parse_args()
    print(opt)
    
    main(
        opt.save_path, 
        opt.protein_fasta, 
        opt.peptide_tsv, 
        opt.tool_name
        )