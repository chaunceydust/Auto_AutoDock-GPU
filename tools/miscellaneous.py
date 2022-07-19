#!/usr/bin/env python3

import os
import shutil
import multiprocessing

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

    def value_cutoff (self):

        csv = os.path.abspath(self.args.csv)

        if self.args.fn == '' or self.args.fn == 'postproc':
            inp_csv = os.path.splitext(csv)[0] + '_postproc.csv'
            csv = os.path.abspath(inp_csv)
        else:
            pass

        csv_file_name = os.path.basename(csv)
        csv_file_dir = os.path.dirname(csv)
        csv_cutoff_file_name = csv_file_name.replace('.csv','') + '_original.csv'
        dst = os.path.join(csv_file_dir, csv_cutoff_file_name)

        shutil.copyfile(csv, dst)

        df = pd.read_csv(csv, index_col = 0)

        condition = (df.Lowest_binding_energy < (-7.0)) & (df.QED_mean > 0.7)

        df = df[condition].sort_values(by=['Lowest_binding_energy'], ascending=True).reset_index(drop=True)

        df.to_csv(csv)


    def scatterplot(self):
        
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm 


        filename = 'scatter.png'

        csv = os.path.abspath(self.args.csv)
        dirname = os.path.dirname(csv)
        df = pd.read_csv(csv)

        fig, ax = plt.subplots()

        X = df['Lowest_binding_energy']
        Y = df['QED_mean']

        # color = cm.rainbow(np.linspace(0, 1, len(df['clusters'])))

        # plotting
        ax.scatter(X, Y, 80, color='black', linewidths=1.0)
        ax.scatter(X, Y, 80, color='white', linewidths=0)
        ax.scatter(X, Y, 80, color='powderblue', linewidths=0, alpha=0.8)


        # setting
        ax.set_title('Lowest_E - QED')
        ax.set_xlabel('Lowest binding energy', size=12)
        ax.set_ylabel('QED', size=12)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

        dst = os.path.join(dirname, filename)

        fig.savefig(f'{dst}', bbox_inches='tight', dpi=300, transparent=True)


    def copypdbqt(self):

        csv = os.path.abspath(self.args.csv)
        qt_path = os.path.abspath(self.args.qtpath)

        if os.path.isdir(qt_path) == True:
            pass
        else:
            print (f'Please check the {qt_path}')
            quit()

        qt_path_new = qt_path + '_copy'

        df = pd.read_csv(csv)

        if 'description' in df.columns:
            pdbqt_list = df['description'].to_list()
        else:
            print ('The csv file should have the "description" column that contains pdbqt file names')
            quit()

        try:
            os.mkdir(qt_path_new)
        except FileExistsError:
            print (f'Please replace the directory [{qt_path_new}]')
            quit()

        for i in pdbqt_list:
            shutil.copyfile (f'{qt_path}/{i}', f'{qt_path_new}/{i}')
