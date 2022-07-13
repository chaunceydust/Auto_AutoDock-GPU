#!/usr/bin/env python3

import os
import time


class RunDock:

    def __init__ (self, args):

        self.args = args

    def autodock_gpu_run (self):

        if os.path.isfile (self.args.listpath) == True:
            listpath = os.path.abspath(self.args.listpath)
        else:
            print ('[rundock] Please check the listpath')
            quit()

        listpath = os.path.abspath(self.args.listpath)


        if os.path.isfile(f'{listpath}'):
            print ('-----------------------------')
            print ('AutoDock-GPU running ...')

            try:
                os.mkdir (self.args.resultpath)
            except FileExistsError:
                pass

            cur_path = os.getcwd()
            os.chdir(self.args.resultpath)
            os.system (f'{self.args.autodockbin} --nrun 1 -filelist {listpath} --devnum {self.args.gpu}')
            os.chdir(cur_path)

        else:
            print ('[rundock] list.txt file for running AutoDock-GPU is not existed.')
            exit()

