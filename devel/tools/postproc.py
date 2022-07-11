#!/usr/bin/env python3

import os
import time
import pandas as pd
from xml.etree.ElementTree import parse


class PostProc:

    def __init__ (self, args):

        self.args = args

    def result2df (self):

        a = 1
        results = os.listdir(self.args.resultpath)
        results_xmls = [file for file in results if file.endswith('.xml')]

        leng = len(results_xmls) + 1

        df = pd.DataFrame (columns=['description', 'Lowest_binding_energy', 'Mean_binding_energy', 'ZINC'])
        df2_ls = []

        for i in results_xmls:
            tree = parse(f'{self.args.resultpath}/{i}')
            root = tree.getroot()
            cluhis = root.findall ('clustering_histogram')

            score = [x.find('cluster').attrib for x in cluhis][0]
            lowest_E = float(score['lowest_binding_energy'])
            mean_E = float(score['mean_binding_energy'])

            name_pdbqt = i.replace('xml', 'pdbqt')

            name_dlg = i.replace ('xml', 'dlg')

            with open (f'{self.args.resultpath}/{name_dlg}', 'r') as data:
                lines = data.readlines()[78]
                ZINC_code = lines[35:].replace('\n', '')

            df2 = pd.DataFrame([[name_pdbqt, lowest_E, mean_E, ZINC_code]], columns=['description', 'Lowest_binding_energy', 'Mean_binding_energy', 'ZINC'])
            df2_ls.append(df2)

            a += 1
            ratio = a/leng * 100
            print(f'{ratio:.2f} %')

        df = pd.concat(df2_ls, ignore_index=True)

        df = df[['ZINC', 'description', 'Lowest_binding_energy', 'Mean_binding_energy']]
        
        # df sorting 추가하기

        df.reset_index(drop=True).to_csv(f'{os.getcwd()}/results.csv')
