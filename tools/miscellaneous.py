#!/usr/bin/env python3

import os
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