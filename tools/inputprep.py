#!/usr/bin/env python3

import multiprocessing
import os
import shutil
import subprocess

from tqdm import tqdm
import numpy as np


class InputPrep:

    def __init__ (self, args):

        self.args = args        


    def smi2pdbqt_list (self):

        if os.path.isdir (self.args.ligandpath) == True:
            ligandpath = os.path.abspath(self.args.ligandpath)
        else:
            print ('Please check the ligandpath')
            quit()

        smi_listdir = os.listdir(ligandpath)
        smi_file_list = [file for file in smi_listdir if file.endswith('.smi')]

        newpath = ligandpath + '_pdbqt'
        new_smi_path = ligandpath + '_smi_original'

        try:
            os.mkdir (newpath)
        except FileExistsError:
            print (f'Please move the directory {newpath} to other path')
            quit()

        try:
            os.mkdir (new_smi_path)
        except FileExistsError:
            print (f'Please move the directory {new_smi_path} to other path')
            quit()

        for i in tqdm(smi_file_list):
            with open (f'{ligandpath}/{i}', 'r') as a_smi_file:
                ls = a_smi_file.readlines()[1:]
                j = i.replace('.smi', '')
                with open (f'{ligandpath}/{j}_mod.smi', 'a') as file:
                    file.writelines(ls)
                shutil.move(f'{ligandpath}/{i}', f'{new_smi_path}/{i}')
                file.close()
            a_smi_file.close()

        mod_smi_listdir = os.listdir(ligandpath)
        mod_smi_file_list = [file for file in mod_smi_listdir if file.endswith('.smi')]

        return mod_smi_file_list

    def smi2pdbqt_paraset(self, mod_smi_file_list):

        ligandpath = os.path.abspath(self.args.ligandpath)
        newpath = ligandpath + '_pdbqt'

        for i in tqdm(mod_smi_file_list):
            j = i.replace('.smi', '')
            cmd = ['obabel', '-i', 'smi', f'{ligandpath}/{i}', '-o', 'pdbqt', '--gen3d', '--conformer', '--nconf', '20', '--score', 'energy', '--AddPolarH', '--partialcharge', 'gasteiger', '-O', f'{newpath}/{j}.pdbqt']
            subprocess.run(cmd, capture_output=True)

    def smi2pdbqt (self):

        num_partitions = self.args.np
        mod_smi_file_list = self.smi2pdbqt_list()

        ls_split = np.array_split(mod_smi_file_list, num_partitions)
        pool = multiprocessing.Pool(processes=num_partitions)
        pool.map(self.smi2pdbqt_paraset, ls_split)
        pool.close()
        pool.join()


    def split_ligs (self):

        print ('-----------------------------')
        print ('splitligs running ...')

        ligandpath = os.path.abspath(self.args.ligandpath)
        ligs = os.listdir(ligandpath)
        ligs_list = [file for file in ligs if file.endswith('.pdbqt')]

        vinapath = os.path.abspath(self.args.vinapath)

        if os.path.basename(vinapath) == 'vina_split':
            pass
        else:
            print ('Please check the vinapath (.../vina_split (bin_file))')
            quit()

        newpath = f'{ligandpath}_original'

        try:
            os.mkdir(newpath)
        except FileExistsError:
            pass

        for i in ligs_list:
            os.system(f'{vinapath} --input {ligandpath}/{i}')
            os.rename(f'{ligandpath}/{i}', f'{newpath}/{i}')


    def listgen (self):

        print ('-----------------------------')
        print ('listgen running ...')

        proteinpath = os.path.abspath(self.args.proteinpath)
        ligandpath = os.path.abspath(self.args.ligandpath)
        listpath = os.path.abspath(self.args.listpath)

        ligs = os.listdir(ligandpath)
        ligs_list = [file for file in ligs if file.endswith('.pdbqt')]

        if len (ligs_list) == 0:
            print (f'There is no ligand files in the {ligandpath} ! ')
            quit ()

        try:
            file = open (f'{listpath}', 'a')
            file.write(proteinpath + '\n')
        except FileExistsError:
            print (f'Please remove the {listpath} file')

        for i in ligs_list:
            j = i.replace('.pdbqt', '')
            file.write(ligandpath + '/' + i + '\n')
            file.write(j + '\n')

        file.close()
