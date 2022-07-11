#!/usr/bin/env python3

import os
import time


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

        try:
            file = open (f'{listpath}/list.txt', 'a')
            file.write(proteinpath + '\n')
        except FileExistsError:
            print (f'Please remove the list.txt file in "{listpath}"')

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

        if os.path.isdir(vinapath) == True:
            pass
        else:
            print ('Please check the vinapath')
            quit()

        newpath = f'{ligandpath}_original'

        try:
            os.mkdir(newpath)
        except FileExistsError:
            pass

        for i in ligs_list:
            j = i.replace('.pdbqt', '')
            os.system(f'{vinapath}/vina_split --input {ligandpath}/{i}')
            os.rename(f'{ligandpath}/{i}', f'{newpath}/{i}')

        # return newpath



# class PathCheck:

#     def __init__ (self, args):
        
#         self.args = args

#     def path_check (self):

#         path_arg_dict = {'proteinpath':self.args.proteinpath, 'ligandpath':self.args.ligandpath, 'vinapath':self.args.vinapath, 'listpath':self.args.listpath}

#         new_dict = {}

#         none_list = []
#         for key, values in path_arg_dict.items():
#             if values == None:
#                 none_list.append(key)
#             else:
#                 new_dict[key] = values

#         non_path = ', '.join(none_list)
#         print (f'Unspecified path: {non_path}')

#         fn_arg_dict = {'listgen':['proteinpath', 'ligandpath', 'listpath'], 'splitligs':['ligandpath', 'vinapath'], 'run_docking':['listpath', 'resultpath']}

#         if self.args.fn == '': # check path for autorun
#             for key, values in new_dict.items():
#                 if key == 'proteinpath':
#                     if os.path.isfile(str(self.args.proteinpath)) == True:
#                         pass
#                     else:
#                         print ('Please check the proteinpath')
#                         print ('example: --proteinpath /path/to/protein.maps.fld (=.maps.fld file)')
#                         quit()

#                 else:
#                     if os.path.isdir(str(values)) == True:
#                         pass
#                     else:
#                         _dst = key.replace('path','')
#                         print (f'Please check the {key}')
#                         print (f'example: --{key} /path/to/{_dst} (=directory)')
#                         quit()

#         else:
#             for fn, arg_list in fn_arg_dict.items():
#                 if self.args.fn == fn:
#                     for i in arg_list:
#                         if i == 'proteinpath':
#                             if os.path.isfile(str(self.args.proteinpath)) == True:
#                                 pass
#                             else:
#                                 print ('Please check the proteinpath')
#                                 print ('example: --proteinpath /path/to/protein.maps.fld (=.maps.fld file)')
#                                 quit()
#                         else:
#                             if os.path.isdir(str(values)) == True:
#                                 pass
#                             else:
#                                 _dst = key.replace('path','')
#                                 print (f'Please check the {key}')
#                                 print (f'example: --{key} /path/to/{_dst} (=directory)')
#                                 quit()

#                 else:
#                     pass

#         print ('Input arguments are correct !')
#         print ('Running will be continued in 5 seconds ...')
#         time.sleep (5)