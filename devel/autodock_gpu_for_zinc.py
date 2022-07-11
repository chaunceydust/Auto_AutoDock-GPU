#!/usr/bin/env python3


import os
import time
import argparse
# import multiprocessing


from tools import inputprep
from tools import rundock
from tools import postproc

parser = argparse.ArgumentParser (description='Tools for high-throughput docking screening using autodock-gpu')

# path
parser.add_argument('-p', '--proteinpath', required=False, type=str, default=None, help='Set the path to maps.fld file for receptor')
parser.add_argument('-l', '--ligandpath', required=False, type=str, default=None, help='Set the path to directory containing ligands')
parser.add_argument('-v', '--vinapath', required=False, type=str, default=None, help='Path to autodock_vina')
parser.add_argument('-ls', '--listpath', required=False, type=str, default=None, help='Path to list.txt')
parser.add_argument('-d', '--resultpath', required=False, type=str, default=None, help='Set the path to result directory')


# setting
parser.add_argument('-lf', '--ligandformat', required=False, type=str, default='pdbqt', help='Format of the ligand files')
parser.add_argument('-smi', '--smi', dest='smi_dest', action='store_true')
parser.add_argument('-sl', '--splitligs', dest='splitligs_dest', action='store_true')


# etc
parser.add_argument('--np', required=False, default='8', help='Number of cores for execute')
parser.add_argument('--fn', required=False, type=str, default='', help='A function name to use')

args = parser.parse_args()


###################################################################################

class RunScript:

    def __init__ (self, args):

        self.args = args

    def fxn_select (self):

        start = time.time ()

        inp = inputprep.InputPrep(args)
        run_docking = rundock.RunDock(args)

        fn_dict = {'listgen':inp.listgen, 'splitligs':inp.split_ligs, 'run_docking':run_docking.autodock_gpu_run}

        if args.fn in fn_dict:
            fn_ = fn_dict[args.fn]
            fn_()
            print ('Elapsed time: ', time.time() - start, 'sec')
            quit()
        else:
            fxns = '''
            listgen (proteinpath, ligandpath, listpath)
            splitligs (ligandpath, vinapath) 
            run_docking (listpath)
            '''

            print ('Please select one of the functions below and enter it.')
            print (fxns)
            quit()

    def auto_run (self):

        start = time.time ()

        inp = inputprep.InputPrep(args)
        run_docking = rundock.RunDock(args)
        post_proc = postproc.PostProc(args)

        if self.args.splitligs_dest == True:
            new_ligpath = inp.split_ligs()
            # self.args.ligandpath = new_ligpath
        else:
            pass
            
        inp.listgen()
        run_docking.autodock_gpu_run()
        post_proc.result2df()

        print ('Elapsed time: ', time.time() - start, 'sec')
        quit()




class PathCheck:

    def __init__ (self, args):
        
        self.args = args

    def path_check (self):

        path_arg_dict = {'proteinpath':self.args.proteinpath, 'ligandpath':self.args.ligandpath, 'vinapath':self.args.vinapath, 'listpath':self.args.listpath}

        new_dict = {}

        none_list = []
        for key, values in path_arg_dict.items():
            if values == None:
                none_list.append(key)
            else:
                new_dict[key] = values

        non_path = ', '.join(none_list)
        print (f'Unspecified path: {non_path}')

        fn_arg_dict = {'listgen':['proteinpath', 'ligandpath', 'listpath'], 'splitligs':['ligandpath', 'vinapath'], 'run_docking':['listpath', 'resultpath']}

        if self.args.fn == '': # check path for autorun
            for key, values in new_dict.items():
                if key == 'proteinpath':
                    if os.path.isfile(str(self.args.proteinpath)) == True:
                        pass
                    else:
                        print ('Please check the proteinpath')
                        print ('example: --proteinpath /path/to/protein.maps.fld (=.maps.fld file)')
                        quit()

                else:
                    if os.path.isdir(str(values)) == True:
                        pass
                    else:
                        _dst = key.replace('path','')
                        print (f'Please check the {key}')
                        print (f'example: --{key} /path/to/{_dst} (=directory)')
                        quit()

        else:
            for fn, arg_list in fn_arg_dict.items():
                if self.args.fn == fn:
                    for i in arg_list:
                        if i == 'proteinpath':
                            if os.path.isfile(str(self.args.proteinpath)) == True:
                                pass
                            else:
                                print ('Please check the proteinpath')
                                print ('example: --proteinpath /path/to/protein.maps.fld (=.maps.fld file)')
                                quit()
                        else:
                            if os.path.isdir(str(values)) == True:
                                pass
                            else:
                                _dst = key.replace('path','')
                                print (f'Please check the {key}')
                                print (f'example: --{key} /path/to/{_dst} (=directory)')
                                quit()

                else:
                    pass

        print ('Input arguments are correct !')
        print ('Running will be continued in 5 seconds ...')
        time.sleep (5)

###################################################################################


if __name__ == '__main__':

    run_script = RunScript(args)
    check = PathCheck(args)

    check.path_check()

    if args.fn == '':
        pass
    else:
        run_script.fxn_select()

    run_script.auto_run()