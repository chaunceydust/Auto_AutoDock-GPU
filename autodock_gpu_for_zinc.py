#!/usr/bin/env python3

import time
import argparse
# import multiprocessing

from tools import inputprep
from tools import rundock
from tools import postproc
from tools import check
from tools import miscellaneous

parser = argparse.ArgumentParser (description='Tools for high-throughput docking screening using autodock-gpu')

# path
parser.add_argument('-p', '--proteinpath', required=False, type=str, default=None, help='Set the path to maps.fld file for receptor (file)')
parser.add_argument('-l', '--ligandpath', required=False, type=str, default=None, help='Set the path to directory containing ligands (directory)')
parser.add_argument('-v', '--vinapath', required=False, type=str, default=None, help='Path to autodock_vina (directory)')
parser.add_argument('-ls', '--listpath', required=False, type=str, default=None, help='Path to list.txt (file)')
parser.add_argument('-d', '--resultpath', required=False, type=str, default=None, help='Set the path to result directory (directory)')

# setting
parser.add_argument('-lf', '--ligandformat', required=False, type=str, default='pdbqt', help='Format of the ligand files')
parser.add_argument('-bin', '--autodockbin', required=False, type=str, default='autodock_gpu_128wi', help='Binary file of AutoDock-GPU (bin_file)')
parser.add_argument('--gpu', required=False, type=int, default=1, help='Which GPU do you use ? (starts at 1)')

# switch
parser.add_argument('-smi', '--smi', dest='smi_dest', action='store_true')
parser.add_argument('-sl', '--splitligs', dest='splitligs_dest', action='store_true')
parser.add_argument('-c', '--continue', dest='continue_dest', action='store_true')

# etc
parser.add_argument('--np', required=False, type=int, default=8, help='Number of cores for execute')
parser.add_argument('--fn', required=False, type=str, default='', help='A function name to use')
parser.add_argument('--csv', required=False, type=str, default='results.csv', help='csv file for data parsing from ZINC database')
parser.add_argument('--splitnum', required=False, type=int, default=4, help='How many divided list files you want ?')

args = parser.parse_args()


class RunScript:

    def __init__ (self, args):

        self.args = args

    def fxn_select (self):

        start = time.time()

        inp = inputprep.InputPrep(args)
        run_docking = rundock.RunDock(args)
        post_proc = postproc.PostProc(args)
        path_check = check.PathCheck(args)
        zn_parsing = postproc.ZINCParsing(args)
        misc = miscellaneous.Misc(args)

        fn_dict = {'listgen':inp.listgen, 'splitligs':inp.split_ligs, 'run_docking':run_docking.autodock_gpu_run, 'dlg2qt':post_proc.dlg2qt, 'result2df':post_proc.result2df, 'znparsing':zn_parsing.fxns_znparsing, 'listsplit':misc.listsplit}

        if args.fn in fn_dict:
            fn_ = fn_dict[args.fn]
            new_path_dict = path_check.pathcheck()
            args.vinapath = new_path_dict['vinapath']
            args.resultpath = new_path_dict['resultpath']
            args.listpath = new_path_dict['listpath']
            fn_()
            print ('Elapsed time: ', time.time() - start, 'sec')
            quit()

        else:
            fxns = '''
            - Essential fxns
            1) splitligs (ligandpath, vinapath) 
            2) listgen (proteinpath, ligandpath, listpath)
            3) run_docking (listpath, autodockbin, resultpath, gpu)
            4) dlg2qt (resultpath)
            5) result2df (resultpath)
            6) znparsing (np, csv)

            - Miscellaneous fxns
            listsplit (listpath, splitnum)
            '''

            print ('Please select/enter the one of functions below and enter it.')
            print (fxns)
            quit()


    def auto_run (self):

        start = time.time ()

        inp = inputprep.InputPrep(self.args)
        run_docking = rundock.RunDock(self.args)
        post_proc = postproc.PostProc(self.args)
        zn_parsing = postproc.ZINCParsing(self.args)

        if self.args.splitligs_dest == True:
            inp.split_ligs()
        else:
            pass
            
        inp.listgen()
        run_docking.autodock_gpu_run()
        post_proc.dlg2qt()
        post_proc.result2df()
        zn_parsing.fxns_znparsing()

        print ('Elapsed time: ', time.time() - start, 'sec')
        quit()



if __name__ == '__main__':

    run_script = RunScript(args)
    path_check = check.PathCheck(args)

    if args.fn == '':
        pass
    else:
        run_script.fxn_select()

    new_path_dict = path_check.pathcheck()
    args.vinapath = new_path_dict['vinapath']
    args.resultpath = new_path_dict['resultpath']
    args.listpath = new_path_dict['listpath']


    args_dict = vars(args)
    for key, value in args_dict.items():
        print (f'[{key}] {value}')
    print ('-----------------------------')
    cont_ = input('Please check your inputs before continuing !! (y / n): ')
    
    if cont_ == 'n':
        quit()
    else:
        pass

    run_script.auto_run()