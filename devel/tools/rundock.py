#!/usr/bin/env python3

import os
import time


class RunDock:

    def __init__ (self, args):

        self.args = args

    def autodock_gpu_run (self):

        if os.path.isdir (self.args.listpath) == True:
            listpath = os.path.abspath(self.args.listpath)
        else:
            print ('Please check the listpath')
            quit()

        listpath = os.path.abspath(self.args.listpath)


        if os.path.isfile(f'{listpath}/list.txt'):
            print ('-----------------------------')
            print ('AutoDock-GPU will be run in 5 seconds ...')

            try:
                os.mkdir (self.args.resultpath)
            except FileExistsError:
                pass

            time.sleep (5)
            cur_path = os.getcwd()
            os.chdir(self.args.resultpath)
            os.system (f'autodock_gpu_128wi -filelist {listpath}/list.txt')
            os.chdir(cur_path)

        else:
            print ('list.txt file for running AutoDock-GPU is not existed.')
            exit()

