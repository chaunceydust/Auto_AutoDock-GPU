#!/usr/bin/env python3

import os
import requests
import numpy as np
import pandas as pd


from xml.etree.ElementTree import parse
from bs4 import BeautifulSoup

import multiprocessing

from rdkit import Chem
from qed import qed


class PostProc:

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
            root = tree.getroot()
            cluhis = root.findall ('clustering_histogram')

            score = [x.find('cluster').attrib for x in cluhis][0]
            lowest_E = float(score['lowest_binding_energy'])
            mean_E = float(score['mean_binding_energy'])

            name_pdbqt = i.replace('xml', 'pdbqt')

            name_dlg = i.replace ('xml', 'dlg')

            with open (f'{resultpath}/{name_dlg}', 'r') as data:
                lines = data.readlines()[78]
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


class ParallelizeParsing:
    
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

    def znparsing(self, df): # return pandas dataframe

        zinc_ls = df.iloc[:,1] # ZINC_code
        pdbqt = df.iloc[:,2] # pdbqt
        lowe = df.iloc[:,3] # lowest binding energy
        mean = df.iloc[:,4] # mean binding energy

        row_leng = len(zinc_ls)
        i = 1

        df2_ls = []
        for _zinc_ls, _pdbqt, _lowe, _mean in zip(zinc_ls, pdbqt, lowe, mean):

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

                qmo = str(qed.weights_mean(mol))

                num_chirals = len(chir)
                chirals = str(chir)

                df2 = pd.DataFrame([[_zinc_ls, _pdbqt, qmo, _lowe, _mean, Mwt, logP, Rbonds, Hdonors, Hacceptors, Psa, NetC, num_chirals, chirals, smiles]], columns=['ZINC', 'pdbqt', 'QED_mean', 'Lowest_E', 'Mean_E', 'Mwt', 'LogP', 'Rotatable bonds', 'H-donors', 'H-acceptors', 'PSA', 'Net charge', 'Num chirals', 'chirals', 'SMILES'])

                df2_ls.append(df2)
                print (f'[{i}/{row_leng}]: {_zinc_ls} --- Done !')
                i += 1

            except (AttributeError, TypeError, ValueError):
                print (f'[{i}/{row_leng}]: {_zinc_ls} --- Error !')
                i += 1

        df = pd.concat(df2_ls, ignore_index=True)

        return df

    
class ZINCParsing:

    def __init__ (self, args):

        self.args = args

    def fxns_znparsing(self):

        multi = ParallelizeParsing(self.args)

        csv_file = os.path.abspath(self.args.csv)
        processed_csv_file = csv_file.replace('.csv', '') + '_postproc.csv'

        data_original = pd.read_csv(csv_file)
        dt = multi.parallelize (data_original, multi.znparsing)
        dt.reset_index(drop=True).to_csv(processed_csv_file)


