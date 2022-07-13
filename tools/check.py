#!/usr/bin/env python3

import os
import time
from types import NoneType

class PathCheck:

    def __init__ (self, args):
        
        self.args = args


    def pathcheck (self):

        args_dict = vars(self.args)

        predesignated_paths = {'vinapath':'/opt/vina', 'resultpath':'./result', 'listpath':'./list.txt'}

        if self.args.fn == '': # autorun
            
            none_list = [] # NoneType list

            print (' ')
            for key, values in args_dict.items():
                if key == 'proteinpath':
                    if os.path.isfile(values) == True:
                        pass
                    else:
                        print ('Please check the proteinpath')
                        print ('example: --proteinpath /path/to/protein.maps.fld (=.maps.fld file)')
                        quit()

                if key == 'vinapath' or key == 'resultpath' or key == 'listpath':
                    if isinstance (values, NoneType) == True:
                        print (f'* {key} automatically set to "{predesignated_paths[key]}"')
                    else:
                        if key == 'vinapath':
                            predesignated_paths[key] = self.args.vinapath
                            print (f'* {key} set to "{predesignated_paths[key]}" by user')
                        elif key == 'resultpath':
                            predesignated_paths[key] = self.args.resultpath
                            print (f'* {key} set to "{predesignated_paths[key]}" by user')                           
                        elif key == 'listpath':
                            predesignated_paths[key] = self.args.listpath
                            print (f'* {key} set to "{predesignated_paths[key]}" by user')                               
                        

                else:
                    if values == None:
                        none_list.append(key)

                if len(none_list) == 0:
                    pass
                else:
                    non_list_to_str = ', '.join(none_list)
                    print (f'Unspecified path: {non_list_to_str}')
                    print (f'Please check input arguments')
                    quit()

            print (' ')
            print ('Input arguments are correct !')
            print ('-----------------------------')


        else:

            fn_arg_dict = {'listgen':['proteinpath', 'ligandpath', 'listpath'], 'splitligs':['ligandpath', 'vinapath'], 'run_docking':['listpath', 'autodockbin', 'resultpath'], 'dlg2qt':['resultpath'], 'result2df':['resultpath'], 'znparsing':['np', 'csv'], 'listsplit':['listpath']}

            for fxn, arg_list in fn_arg_dict.items():
                if self.args.fn == fxn:
                    for i in arg_list:
                        if i == 'proteinpath':
                            if os.path.isfile(str(self.args.proteinpath)) == True:
                                pass
                            else:
                                print ('Please check the proteinpath')
                                print ('example: --proteinpath /path/to/protein.maps.fld (=.maps.fld file)')
                                quit()

                        elif i == 'listpath': # return path
                            if fxn == 'listgen':
                                if args_dict[i] == None:
                                    pass
                                else:
                                    if os.path.splitext(self.args.listpath)[1] == '.txt':
                                        predesignated_paths[i] = self.args.listpath
                                    else:
                                        print (f'Please check the listpath // input: {args_dict[i]}')
                                        print ('example: --listpath /path/to/list.txt')
                                        quit()       
                            else:
                                if os.path.isfile(str(self.args.listpath)) == True:
                                    predesignated_paths[i] = self.args.listpath
                                else:
                                    print (f'Please check the listpath // input: {args_dict[i]}')
                                    print ('example: --listpath /path/to/list.txt')
                                    quit()

                        elif i == 'vinapath': # return path
                            if args_dict[i] == None:
                                pass
                            elif os.path.isdir(str(self.args.vinapath)) == True:
                                predesignated_paths[i] = self.args.vinapath
                            else:
                                print (f'Please check the vinapath // input: {args_dict[i]}')
                                print ('example: --vinapath /path/to/vina_bin (=directory)')
                                quit()
 
                        elif i == 'resultpath': # return path
                            if fxn == 'run_docking':
                                if args_dict[i] == None:
                                    pass
                                else:
                                    if os.path.splitext(self.args.resultpath)[1] == '':
                                        predesignated_paths[i] = self.args.resultpath
                                    else:
                                        print (f'Please check the resultpath // input: {args_dict[i]}')
                                        print ('example: --resultpath /path/to/result_files (=directory)')
                                        quit() 
                            else:
                                if os.path.isdir(str(self.args.resultpath)) == True:
                                    predesignated_paths[i] = self.args.resultpath
                                else:
                                    print (f'Please check the resultpath // input: {args_dict[i]}')
                                    print ('example: --resultpath /path/to/result_files (=directory)')
                                    quit()

                        elif i == 'autodockbin':
                            if os.path.isfile(str(self.args.autodockbin)) == True:
                                pass
                            else:
                                print (f'Please check the AutoDock-GPU binary // input: {args_dict[i]}')
                                print ('example: --autodockbin /path/to/autodock_gpu_binary (=binary file)')
                                quit()
                                
                        elif i == 'np':
                            pass

                        elif i == 'csv':
                            if os.path.isfile(str(self.args.csv)) == True:
                                pass
                            else:
                                print (f'Please check the csv file // input: {args_dict[i]}')
                                print ('example: --csv /path/to/csv_file.csv (=csv file)')
                                quit()
                        else:
                            if os.path.isdir(str(args_dict[i])) == True:
                                pass
                            else:
                                _dst = i.replace('path','')
                                print (f'Please check the {i} // input: {args_dict[i]}')
                                print (f'example: --{i} /path/to/{_dst} (=directory)')
                                quit()
            
            for arg, path in predesignated_paths.items():
                print (f'* {arg} set to {path}')

        return predesignated_paths