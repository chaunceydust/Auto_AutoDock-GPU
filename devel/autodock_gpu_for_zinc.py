#!/usr/bin/env python3

import time
import argparse
# import multiprocessing

from tools import inputprep
from tools import rundock
from tools import postproc
from tools import check

parser = argparse.ArgumentParser (description='Tools for high-throughput docking screening using autodock-gpu')

# path
parser.add_argument('-p', '--proteinpath', required=False, type=str, default=None, help='Set the path to maps.fld file for receptor')
parser.add_argument('-l', '--ligandpath', required=False, type=str, default=None, help='Set the path to directory containing ligands')
parser.add_argument('-v', '--vinapath', required=False, type=str, default='/opt/vina/', help='Path to autodock_vina')
parser.add_argument('-ls', '--listpath', required=False, type=str, default=None, help='Path to list.txt')
parser.add_argument('-d', '--resultpath', required=False, type=str, default=None, help='Set the path to result directory')

# setting
parser.add_argument('-lf', '--ligandformat', required=False, type=str, default='pdbqt', help='Format of the ligand files')

# switch
parser.add_argument('-smi', '--smi', dest='smi_dest', action='store_true')
parser.add_argument('-sl', '--splitligs', dest='splitligs_dest', action='store_true')

# etc
parser.add_argument('--np', required=False, type=int, default=8, help='Number of cores for execute')
parser.add_argument('--fn', required=False, type=str, default='', help='A function name to use')
parser.add_argument('--csv', required=False, type=str, default='results.csv', help='csv file for data parsing from ZINC database')

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

        fn_dict = {'listgen':inp.listgen, 'splitligs':inp.split_ligs, 'run_docking':run_docking.autodock_gpu_run, 'dlg2qt':post_proc.dlg2qt, 'result2df':post_proc.result2df, 'znparsing':zn_parsing.fxns_znparsing}

        if args.fn in fn_dict:
            fn_ = fn_dict[args.fn]
            vina_path, resultpath = path_check.pathcheck()
            args.vinapath = vina_path
            fn_()
            print ('Elapsed time: ', time.time() - start, 'sec')
            quit()

        else:
            fxns = '''
            listgen (proteinpath, ligandpath, listpath)
            splitligs (ligandpath, vinapath) 
            run_docking (listpath)
            dlg2qt (resultpath)
            result2df (resultpath)
            znparsing (np, csv)
            '''

            print ('Please select one of the functions below and enter it.')
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

    vina_path, resultpath = path_check.pathcheck()
    args.vinapath = vina_path
    args.resultpath = resultpath

    run_script.auto_run()