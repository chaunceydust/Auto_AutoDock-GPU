#!/usr/bin/env python3

import os
import subprocess
import multiprocessing
import time

import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.SimDivFilters import rdSimDivPickers

import numpy as np
import pandas as pd


class Misc:

    def __init__ (self, args):

        self.args = args

    def parallelize (self, inp_df, func):

        num_partitions = self.args.np

        df_split = np.array_split(inp_df, num_partitions)
        pool = multiprocessing.Pool(num_partitions)
        out_df = pd.concat(pool.map(func, df_split))
        pool.close()
        pool.join()

        return out_df

    def listsplit (self):

        list_file = os.path.abspath(self.args.listpath)
        split_num = self.args.splitnum

        ls_file_mod = list_file.replace('.txt', '')

        with open (list_file, 'r') as file:
            list_lines = file.readlines()
            protein_path = list_lines[0]
            ligand_path = list_lines[1:]

            num = round(len(ligand_path) / split_num)

            if num%2==1:
                num = num + 1 # odd to even
            else:
                pass # even

            for i in range(split_num):
                j = i + 1
                with open (f'{ls_file_mod}_{j}.txt', 'a') as split:
                    split.writelines(protein_path)
                    split.writelines(ligand_path[i * num : (i + 1) * num])


class Parallel:

    def __init__ (self, args):

        self.args = args

    def obabel_func (self, df):

        if self.args.qtpath == None:
            qt_path = self.args.resultpath + '_pdbqt'
        else:
            qt_path = os.path.abspath(self.args.qtpath)

        df2_ls = []
        pdbqt_list = df['description'].to_list()

        for i, j in enumerate(pdbqt_list):
            cmd = ['obabel', '-i', 'pdbqt', f'{qt_path}/{j}', '-o', 'smi']
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            smi = stdout.decode('UTF-8').split()[0]

            df2 = pd.DataFrame([[j, smi]], columns=['description', 'SMILES'])
            df2_ls.append(df2)

        new_df = pd.concat(df2_ls, ignore_index=True)

        return new_df


    def qedcalc_func (self, df):

        if self.args.qtpath == None:
            qt_path = self.args.resultpath + '_pdbqt'
        else:
            qt_path = os.path.abspath(self.args.qtpath)


        df2_ls = []

        col = 'SMILES'

        if col in df.columns:
            smi_list = df[col].to_list()
            pdbqt_list = df['description'].to_list()
        else:
            print ('The csv file should have the "SMILES" and "description" column')
            quit ()

        for i, j in zip (smi_list, pdbqt_list):
            try:
                canon_smi = Chem.CanonSmiles(i)
                mol = Chem.MolFromSmiles(canon_smi)
                qmo = rdkit.Chem.QED.qed(mol)
            except:
                qmo = int(0)

            df2 = pd.DataFrame([[j, qmo]], columns=['description', 'QED_mean'])
            df2_ls.append(df2)

        new_df = pd.concat(df2_ls, ignore_index=True)

        return new_df


    def clustering_prep (self, df):

        if 'SMILES' in df.columns and 'description' in df.columns:
            smi_list = df['SMILES'].to_list()
            pdbqt_list = df['description'].to_list()
        else:
            print ('The csv file should have the "SMILES" and "description" column')
            quit()

        qt_smi_dict = {qt:smi for qt, smi in zip(pdbqt_list, smi_list)}

        qt_smi_dict_pass = {}

        for qt, smi in qt_smi_dict.items():
            try:
                canon_smi = Chem.CanonSmiles(smi)
                qt_smi_dict_pass[qt] = canon_smi
            except:
                pass
            # if smi == None:
            #     pass
            # else:
            #     qt_smi_dict_pass[qt] = smi

        ms = {qt:Chem.MolFromSmiles(smi) for qt, smi in qt_smi_dict_pass.items()}
        fps = {qt:rdMolDescriptors.GetMorganFingerprintAsBitVect(mol,2,2048) for qt, mol in ms.items()}
        
        return fps


    def clustering_func (self, dict):

        # This function was made by referring to the rdkit website.
        # https://greglandrum.github.io/rdkit-blog/similarity/tutorial/2020/11/18/sphere-exclusion-clustering.html

        csv = os.path.abspath(self.args.csv)
        df = pd.read_csv(csv, index_col = 0)

        start = time.time()
        fps = dict

        print ('-----------------------------')
        print ('Subset generating ...')

        lp = rdSimDivPickers.LeaderPicker()
        threshold = 0.65 # clustering threshold
        picks = lp.LazyBitVectorPick(list(fps.values()), len(fps), threshold)

        print ('Num. subsets: ',len(picks))
        mid_time_1 = time.time()
        print ('Elapsed time of subset generation: ', mid_time_1 - start, 'sec')

        def assignPointsToClusters(picks, fps):
            from collections import defaultdict

            clusters = defaultdict(list)
            
            for i,idx in enumerate(picks):
                clusters[i].append(idx)
            
            sims = np.zeros((len(picks),len(fps)))
            for i in range(len(picks)):
                pick = picks[i]
                sims[i,:] = DataStructs.BulkTanimotoSimilarity(fps[pick],fps)
                sims[i,i] = 0 

            best = np.argmax(sims,axis=0)

            for i , idx in enumerate(best):
                if i not in picks:
                    clusters[idx].append(i)
                    
            return clusters

        print ('-----------------------------')
        print ('Mol. clustering ...')

        fps_mol = list(fps.values())

        clusters = assignPointsToClusters(picks, fps_mol)

        df2_ls = []

        for clust, ind in clusters.items():
            for i in ind:
                _description = list(fps.keys())[i]

                df2 = pd.DataFrame([[_description, clust]], columns=['description', 'clusters'])
                df2_ls.append(df2)

        new_df = pd.concat(df2_ls, ignore_index=True)

        name = str(self.args.fn)
        new_dst = os.path.splitext(csv)[0] + f'_{name}.csv'

        df_merge = pd.merge(df, new_df, on='description')
        df_merge.to_csv(csv)

        mid_time_2 = time.time()
        print ('Elapsed time of Mol. clustering: ', mid_time_2 - mid_time_1, 'sec')



    def parallelize (self, func):

        csv = os.path.abspath(self.args.csv)
        df = pd.read_csv(csv, index_col = 0)

        func_name = func.__name__

        num_partitions = self.args.np

        df_split = np.array_split(df, num_partitions)
        pool = multiprocessing.Pool(processes=num_partitions)

        if func_name != 'clustering_prep':
            out_df = pd.concat(pool.map(func, df_split))

            pool.close()
            pool.join()

            df_merge = pd.merge(df, out_df, on='description')

            df_merge.to_csv(csv)

        elif func_name == 'clustering_prep':
            out_dicts = pool.map(func, df_split)
            pool.close()
            pool.join()

            new_dict = {}
            for i in out_dicts:
                new_dict.update(i)
        
            return new_dict



class ParallelRun:

    def __init__ (self, args):

        self.args = args
        self.parallel = Parallel(args)
    
    def __print (self):
        if self.args.fn != '':
            print ('-----------------------------')
            print (f'{self.args.fn} running ...')
        else:
            pass

    def obabel_pararun (self):

        self.__print()
        self.parallel.parallelize (self.parallel.obabel_func)

    def qedcalc_pararun (self):

        self.__print()
        self.parallel.parallelize (self.parallel.qedcalc_func)

    def clustering_pararun (self):

        start = time.time()

        self.__print()
        print ('-----------------------------')
        print ('Mol. fingerprint extrating ...')
        new_dict = self.parallel.parallelize(self.parallel.clustering_prep)
        print ('Elapsed time of fingerprine extraction: ', time.time() - start, 'sec')
        self.parallel.clustering_func(new_dict)