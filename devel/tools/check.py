#!/usr/bin/env python3

import os
import time

class PathCheck:

    def __init__ (self, args):
        
        self.args = args

    def pathcheck (self):

        path_arg_dict = {'proteinpath':self.args.proteinpath, 'ligandpath':self.args.ligandpath, 'vinapath':self.args.vinapath, 'listpath':self.args.listpath, 'resultpath':self.args.resultpath}

        new_dict = {} # input arguments dict

        none_list = [] # NoneType list

        for key, values in path_arg_dict.items():
            if not key == 'resultpath':
                if values == None:
                    none_list.append(key)
                else:
                    new_dict[key] = values
            elif key == 'resultpath':
                if values == None:
                    none_list.append(key)
                else:
                    pass

        if len(none_list) == 0:
            pass

        else:
            if 'vinapath' in none_list:
                print ('Vinapath automatically set to /opt/vina/')
                self.args.vinapath = '/opt/vina/' # return
                none_list.remove('vinapath')

            elif 'resultpath' in none_list:
                print ('Resultpath automatically set to ./result')
                self.args.resultpath = './result' # return
                none_list.remove('resultpath')

        if len(none_list) != 0:
            non_path = ', '.join(none_list)
            print (f'Unspecified path: {non_path}')
            print (f'Please check input arguments')
            quit()


        fn_arg_dict = {'listgen':['proteinpath', 'ligandpath', 'listpath'], 'splitligs':['ligandpath', 'vinapath'], 'run_docking':['listpath', 'resultpath'], 'dlg2qt':['resultpath'], 'result2df':['resultpath']} # 'znparsing':['np', 'csv']}

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
                        print (f'Please check the {key} // input: {path_arg_dict[key]}')
                        print (f'example: --{key} /path/to/{_dst} (=directory)')
                        quit()

        else:
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
                        else:
                            try:
                                if os.path.isdir(str(path_arg_dict[i])) == True:
                                    pass
                                else:
                                    _dst = i.replace('path','')
                                    print (f'Please check the {i} // input: {path_arg_dict[i]}')
                                    print (f'example: --{i} /path/to/{_dst} (=directory)')
                                    quit()
                            except NameError:
                                pass

                else:
                    pass

        print ('Input arguments are correct !')
        print ('Running will be continued in 5 seconds ...')

        time.sleep (5)

        return self.args.vinapath, self.args.resultpath


    def pathcheck_test (self):
        if self.args.fn == '':
            pass
        else:
            pass