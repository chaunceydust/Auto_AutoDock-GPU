#!/usr/bin/env python3

import os

class InputPrep:

    def __init__ (self, args):

        self.args = args        

    def listgen (self):

        def smi2pdbqt (self):

            if os.path.isdir (self.args.ligandpath) == True:
                ligandpath = os.path.abspath(self.args.ligandpath)
            else:
                print ('Please check the ligandpath')
                quit()

            ligandpath = os.path.abspath(self.args.ligandpath)
            ligs = os.listdir(ligandpath)
            ligs_list = [file for file in ligs if file.endswith('smi')]

            newpath = ligandpath + '_'.join(['pdbqt'])

            for i in ligs_list:
                j = i.replace('.smi', '')

                try:
                    os.mkdir (f'{ligandpath}_pdbqt')
                except FileExistsError:
                    pass

                os.system(f'obabel -i smi {ligandpath}/{i} -o {newpath}/{j}.pdbqt --AddPolarH --partialcharge mmff94')

            return newpath


        print ('-----------------------------')
        print ('listgen running ...')

        proteinpath = os.path.abspath(self.args.proteinpath)
        ligandpath = os.path.abspath(self.args.ligandpath)
        listpath = os.path.abspath(self.args.listpath)

        if self.args.smi_dest == True:
            newpath = smi2pdbqt (self.args.ligandpath)
            ligandpath = newpath
        else:
            pass

        ligs = os.listdir(ligandpath)
        ligs_list = [file for file in ligs if file.endswith('pdbqt')]

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


    def split_ligs (self):

        print ('-----------------------------')
        print ('splitligs running ...')


        ligandpath = os.path.abspath(self.args.ligandpath)
        ligs = os.listdir(ligandpath)
        ligs_list = [file for file in ligs if file.endswith('pdbqt')]

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