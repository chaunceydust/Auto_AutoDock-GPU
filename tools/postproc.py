#!/usr/bin/env python3

import os
import time
import subprocess
import multiprocessing
import shutil

import numpy as np
import pandas as pd

# for ZINC parsing
from xml.etree.ElementTree import parse
import requests
from bs4 import BeautifulSoup

# for etc.
import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import QED
from rdkit.SimDivFilters import rdSimDivPickers



class ResultArrange:

    def __init__ (self, args):

        self.args = args

    def dlg2qt (self):

        print ('-----------------------------')
        print ('dlg2qt running ...')

        resultpath = os.path.abspath(self.args.resultpath)
        qt_path = resultpath + '_pdbqt'

        try:
            os.mkdir(qt_path)
        except FileExistsError:
            pass

        dlg_path = os.listdir(resultpath)
        dlg_list = [file for file in dlg_path if file.endswith('.dlg')]

        for i in dlg_list:
            dst = os.path.splitext(i)[0]
            os.system(f'grep "^DOCKED" {resultpath}/{i} | cut -c9- > {qt_path}/{dst}.pdbqt')

    def result2df (self):

        print ('-----------------------------')
        print ('result2df running ...')

        resultpath = os.path.abspath(self.args.resultpath)

        a = 1
        results = os.listdir(resultpath)
        results_xmls = [file for file in results if file.endswith('.xml')]

        leng = len(results_xmls) + 1

        df = pd.DataFrame (columns=['description', 'Lowest_binding_energy', 'Mean_binding_energy', 'ZINC'])
        df2_ls = []

        for i in results_xmls:
            tree = parse(f'{resultpath}/{i}')

            try:
                cluhis = tree.findall('result/clustering_histogram')
                score = [x.find('cluster').attrib for x in cluhis][0]
            except IndexError:
                cluhis = tree.findall('clustering_histogram')
                score = [x.find('cluster').attrib for x in cluhis][0]

            lowest_E = float(score['lowest_binding_energy'])
            mean_E = float(score['mean_binding_energy'])

            name_pdbqt = i.replace('xml', 'pdbqt')

            name_dlg = i.replace ('xml', 'dlg')

            with open (f'{resultpath}/{name_dlg}', 'r') as data:

                dlg_lines = data.readlines()

                lines1 = dlg_lines[78] # len 52
                lines2 = dlg_lines[81]

                if len(lines1) == 52:
                    lines = lines1
                elif len(lines2) == 52:
                    lines = lines2

                ZINC_code = lines[35:].replace('\n', '')

            df2 = pd.DataFrame([[name_pdbqt, lowest_E, mean_E, ZINC_code]], columns=['description', 'Lowest_binding_energy', 'Mean_binding_energy', 'ZINC'])
            df2_ls.append(df2)

            a += 1
            ratio = a/leng * 100
            print(f'{ratio:.2f} %')

        df = pd.concat(df2_ls, ignore_index=True)
        df = df[['ZINC', 'description', 'Lowest_binding_energy', 'Mean_binding_energy']]
    
        csv_file = os.path.abspath(self.args.csv)
        df.reset_index(drop=True).to_csv(csv_file)


class Parallel:

    def __init__ (self, args):

        self.args = args

    def obabel_func (self, df):

        if 'description' in df.columns:
            pdbqt_list = df['description'].to_list()
        else:
            print ('The csv file should have the "description" column that contains pdbqt file names')
            quit()

        print ('-----------------------------')
        print ('Openbabel running ...')           

        if self.args.qtpath == None:
            qt_path = self.args.resultpath + '_pdbqt'
        else:
            qt_path = os.path.abspath(self.args.qtpath)

        df2_ls = []

        for i in pdbqt_list:
            cmd = ['obabel', '-i', 'pdbqt', f'{qt_path}/{i}', '-o', 'smi']
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            smi = stdout.decode('UTF-8').split()[0]

            df2 = pd.DataFrame([[i, smi]], columns=['description', 'SMILES'])
            df2_ls.append(df2)

        new_df = pd.concat(df2_ls, ignore_index=True)

        return new_df

    def znparsing_func(self, df):

        if 'ZINC' in df.columns and 'description' in df.columns:
            zinc = df['ZINC'].to_list()
            pdbqt = df['description'].to_list()
        else:
            print ('The csv file should have the "ZINC" and "description" column')
            quit()

        row_leng = len(zinc)
        i = 1

        df2_ls = []
        for _zinc_ls, _pdbqt in zip(zinc, pdbqt):

            try:
                response = requests.get(f'http://zinc15.docking.org/substances/{str(_zinc_ls)}/').text
                soup = BeautifulSoup(response, 'html.parser')

                # data
                smiles = soup.select_one('#substance-smiles-field')['value']
                Mwt = soup.select_one('body > div > div > div > div:nth-child(2) > div.col-sm-9 > div:nth-child(3) > table:nth-child(1) > tbody > tr > td:nth-child(4)').text
                logP = soup.select_one('body > div > div > div > div:nth-child(2) > div.col-sm-9 > div:nth-child(3) > table:nth-child(1) > tbody > tr > td:nth-child(5)').text
                Rbonds = soup.select_one('body > div > div > div > div.protomers.panel.panel-default.row.panel- > table > tbody > tr > td:nth-child(6)').text
                Hdonors = soup.select_one('body > div > div > div > div.protomers.panel.panel-default.row.panel- > table > tbody > tr > td:nth-child(3)').text
                Hacceptors = soup.select_one('body > div > div > div > div.protomers.panel.panel-default.row.panel- > table > tbody > tr > td:nth-child(4)').text
                Psa = soup.select_one('body > div > div > div > div.protomers.panel.panel-default.row.panel- > table > tbody > tr > td:nth-child(5)').text
                NetC = soup.select_one('body > div > div > div > div.protomers.panel.panel-default.row.panel- > table > tbody > tr > td:nth-child(2)').text

                mol = Chem.MolFromSmiles(smiles)
                chir = Chem.FindMolChiralCenters(mol,force=True,includeUnassigned=True,useLegacyImplementation=True)

                num_chirals = len(chir)
                chirals = str(chir)

                df2 = pd.DataFrame([[_pdbqt, _zinc_ls, Mwt, logP, Rbonds, Hdonors, Hacceptors, Psa, NetC, num_chirals, chirals, smiles]], columns=['description', 'ZINC', 'Mwt', 'LogP', 'Rotatable bonds', 'H-donors', 'H-acceptors', 'PSA', 'Net charge', 'Num chirals', 'chirals', 'SMILES'])

                df2_ls.append(df2)
                print (f'[{i}/{row_leng}]: {_zinc_ls} --- Done !')
                i += 1

            except (AttributeError, TypeError, ValueError):
                print (f'[{i}/{row_leng}]: {_zinc_ls} --- Error !')
                i += 1


        new_df = pd.concat(df2_ls, ignore_index=True)

        return new_df

    def qedcalc_func (self, df):

        df2_ls = []

        if 'SMILES' in df.columns and 'description' in df.columns:
            smi_list = df['SMILES'].to_list()
            pdbqt_list = df['description'].to_list()
        else:
            print ('The csv file should have the "SMILES" and "description" column')
            quit ()

        print ('-----------------------------')
        print ('QED calculating ...')

        for i, j in zip (smi_list, pdbqt_list):
            try:
                canon_smi = Chem.CanonSmiles(i)
                mol = Chem.MolFromSmiles(canon_smi)
                qmo = QED.qed(mol)
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
        threshold = 0.85 # clustering threshold
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

        df_merge = pd.merge(df, new_df, on='description')

        # if self.args.fn == '' or self.args.fn == 'postproc':
        #     dst = csv.replace('.csv', '') + '_postproc.csv'

        # else:
        #     dst = csv.replace('.csv', '') + '_clustering.csv'

        df_merge.to_csv(csv)

        mid_time_2 = time.time()
        print ('Elapsed time of Mol. clustering: ', mid_time_2 - mid_time_1, 'sec')

    def parallelize (self, func):

        csv = os.path.abspath(self.args.csv)
        func_name = func.__name__

        if self.args.fn == '' or self.args.fn == 'postproc':
            dst = csv.replace('.csv', '') + '_postproc.csv'
        else:
            dst = csv.replace('.csv', '') + '_' + func_name.split('_')[0].csv

        if os.path.isfile(dst) == True:
            pass
        else:
            shutil.copyfile (csv, dst)
            
        df = pd.read_csv(dst, index_col = 0)

        num_partitions = self.args.np

        df_split = np.array_split(df, num_partitions)
        pool = multiprocessing.Pool(processes=num_partitions)

        # if self.args.fn == '' or self.args.fn == 'postproc':
        #     dst = csv.replace('.csv', '') + '_postproc.csv'
        # else:
        #     dst = csv.replace('.csv', '') + '_' + func_name.split('_')[0].csv

        if func_name != 'clustering_prep':
            out_df = pd.concat(pool.map(func, df_split))

            pool.close()
            pool.join()

            df_merge = pd.merge(df, out_df, on='description')

            df_merge.to_csv(dst)

        elif func_name == 'clustering_prep':
            out_dicts = pool.map(func, df_split)
            pool.close()
            pool.join()

            new_dict = {}
            for i in out_dicts:
                new_dict.update(i)
        
            return new_dict, dst



class ParallelRun:

    def __init__ (self, args):

        self.args = args
        self.parallel = Parallel(args)
    
    def __print (self):
        if self.args.fn != '' and self.args.fn != 'postproc':
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

    def znparsing_pararun (self):

        self.__print()
        self.parallel.parallelize (self.parallel.znparsing_func)

    def clustering_pararun (self):

        start = time.time()

        self.__print()
        print ('-----------------------------')
        print ('Mol. fingerprint extrating ...')
        new_dict, dst = self.parallel.parallelize(self.parallel.clustering_prep)
        print ('Elapsed time of fingerprine extraction: ', time.time() - start, 'sec')
        self.args.csv = dst
        self.parallel.clustering_func(new_dict)