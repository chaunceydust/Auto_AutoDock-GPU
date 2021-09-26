# sys
import os
import time
import argparse

# data proc
import pandas as pd
import numpy as np

# parsing
import requests
from bs4 import BeautifulSoup
from xml.etree.ElementTree import parse

# multiproc
import multiprocessing
import time

# rdkit _ chirality
from rdkit import Chem

parser = argparse.ArgumentParser (description="Tools for high-throughput docking screening using autodock-gpu")

# General arguments
parser.add_argument("--proteinpath", required=False, default=None, help="Set the path to maps.fld file for receptor")
parser.add_argument("--ligandpath", required=False, default=None, help="Set the path to directory containing ligands")
parser.add_argument("--dlgpath", required=False, default=None, help="Path to dlg files")
parser.add_argument("--np", required=False, default='4', help="Number of cores for execute")
parser.add_argument("--dst", required=False, default='./ZN_parsed.csv', help="Path to result")

# Default setting
parser.add_argument("--ligandfmt", required=False, default="pdbqt", help="Format of the ligand files")
parser.add_argument("--vinapath", required=False, default="/opt/vina/", help="Path to autodock_vina")
parser.add_argument("--dfpath", required=False, default="./", help="Path to the directory for result.df")
parser.add_argument("--listpath", required=False, default="./", help="Path to list.txt")
parser.add_argument("--setnum", required=False, default="1", help="Path to list.txt")

# Select 
parser.add_argument("--listgen", required=False, default="n", help="Generate list for docking (y/n)")
parser.add_argument("--dlg2qt", required=False, default="n", help="Dlg-to-pdbqt convert only (y/n)")
parser.add_argument("--splitligs", required=False, default="n", help="Ligand split (y/n)")
parser.add_argument("--result2df", required=False, default="n", help="Merge data (y/n)")
parser.add_argument("--znparsing", required=False, default="n", help="Parse ZINC db (y/n)")
parser.add_argument("--rearr", required=False, default="n", help="Rearrange result (y/n)")
args = parser.parse_args()


# initial setup
global num_cores
num_cores = int(args.np)

# Fns
def listgen (proteinpath=args.proteinpath, ligandpath=args.ligandpath, ligandfmt=args.ligandfmt, listpath=args.listpath):

    path = ligandpath
    file_list = os.listdir(path)
    file_list_ligands = [file for file in file_list if file.endswith(ligandfmt)]

    with open (f'{listpath}/list.txt', 'a') as file:
        file.write(proteinpath + "\n")
        for i in file_list_ligands:
            j = i.replace("." + ligandfmt, "")
            file.write(ligandpath + '/' + i + "\n")
            file.write(j + "\n")


    print ("* Listgen - Done !")


def dlg2qt(dlgpath = args.dlgpath):

    path = dlgpath
    file_list3 = os.listdir(path)
    file_list3_outputs = [file for file in file_list3 if file.endswith(".dlg")]
    
    for ligs3 in file_list3_outputs:
        dst = os.path.splitext(ligs3)[0]
        os.system(f"grep '^DOCKED' {ligs3} | cut -c9- > {dst}.pdbqt")


    print ("* Dlg-to-pdbqt - Done !")


def splitligs (vinapath = args.vinapath, ligandpath=args.ligandpath):

    path = ligandpath
    file_list = os.listdir(path)
    file_list_ligands = [file for file in file_list if file.endswith(".pdbqt")]

    os.mkdir (f"{ligandpath}/ligand_original/")

    for lig in file_list_ligands:
        os.system(f"{vinapath}/vina_split --input {ligandpath}/{lig}")
        os.system(f"mv {ligandpath}/{lig} {ligandpath}/ligand_original/")


    print ("* Ligand split - done !")


def result2df (dfpath = args.dfpath, dlgpath = args.dlgpath, set = args.setnum, rearr = args.rearr):
    i = 1
    path = dlgpath
    file_list = os.listdir(path)
    file_list_outputs = [file for file in file_list if file.endswith(".xml")]
    leng = len(file_list_outputs) + 1

    if dfpath != "./":
        os.mkdir(f"{dfpath}")
    
    else:
        pass

    df = pd.DataFrame (columns=['file name', 'Lowest_binding_energy', 'Mean_binding_energy', 'ZINC', 'set'])
    df2_ls = []

    for ligs in file_list_outputs:
        tree = parse(f'{ligs}')
        root = tree.getroot()
        cluhis = root.findall ('clustering_histogram')


        score = [x.find('cluster').attrib for x in cluhis][0]
        lb = float(score['lowest_binding_energy'])
        mb = float(score['mean_binding_energy'])
        l = ligs.replace("xml", "pdbqt")


        k = ligs.replace("xml", "dlg")
        with open (f"./{k}", "r") as data2:
            
            lines = data2.readlines()[78]
            lines2 = lines[35:].replace('\n', '')

        # df = df.append(pd.DataFrame([[l, lb, mb, lines2]], columns=['file name', 'Lowest_binding_energy', 'Mean_binding_energy', 'ZINC']), ignore_index=True)
        # l: pdbqt
        # lb: lowest binding energy
        # mb: mean binding energy
        # lines2: ZINC code
        df2 = pd.DataFrame([[l, lb, mb, lines2, set]], columns=['file name', 'Lowest_binding_energy', 'Mean_binding_energy', 'ZINC', 'set'])
        df2_ls.append(df2)

        i += 1
        ratio = i/leng * 100
        print (ratio, '%')

    df = pd.concat(df2_ls, ignore_index=True)

    if rearr == 'n':
        pass
        
    else:
        # data rearrange
        df = pd.concat([df[df['Lowest_binding_energy'] < 0].sort_values(by=['Lowest_binding_energy']), df[df['Lowest_binding_energy'] >= 0].sort_values(by=['Lowest_binding_energy'], ascending=True)])
        
    df.reset_index(drop=True).to_csv(f'./{dfpath}/result_merged.csv')

    print ("* Result merge - Done !")


def zincparsing (data):

    pdbqt = data.iloc[:,1] # pdbqt
    lowe = data.iloc[:,2] # lowest binding energy
    mean = data.iloc[:,3] # mean binding energy
    ls = data.iloc[:,4] # ZINC code
    set = data.iloc[:,5] # set num

    df = pd.DataFrame (columns=['ZINC', 'pdbqt', 'Lowest_E', 'Mean_E', 'Set', 'Mwt', 'LogP', 'Rotatable bonds', 'H-donors', 'H-acceptors', 'PSA', 'Net charge', 'Num chirals', 'chirals', 'SMILES'])

    leng = len(ls)
    i = 1

    df2_ls = []
    for zinc, l, m, n, o in zip(ls, lowe, mean, set, pdbqt):
        
        try:
            response = requests.get(f'http://zinc.docking.org/substances/{zinc}/').text
            soup = BeautifulSoup(response, 'html.parser')

            # data
            smiles = soup.select_one('#substance-smiles-field')['value']
            Mwt = soup.select_one('body > div > div > div > div:nth-child(2) > div.col-sm-9 > div:nth-child(3) > table:nth-child(1) > tbody > tr > td:nth-child(4)').text
            logP = soup.select_one('body > div > div > div > div:nth-child(2) > div.col-sm-9 > div:nth-child(3) > table:nth-child(1) > tbody > tr > td:nth-child(5)').text
            Rbonds = soup.select_one('body > div > div > div > div.protomers.panel.panel-default.row.panel- > table > tbody > tr > td:nth-child(6)').text
            Hdonors = soup.select_one('body > div > div > div > div.protomers.panel.panel-default.row.panel- > table > tbody > tr > td:nth-child(3)').text
            Hacceptors = soup.select_one('body > div > div > div > div.protomers.panel.panel-default.row.panel- > table > tbody > tr > td:nth-child(4)').text
            Psa = soup.select_one('body > div > div > div > div.protomers.panel.panel-default.row.panel- > table > tbody > tr > td:nth-child(5)').text
            NetC = soup.select_one('body > div > div > div > div.protomers.panel.panel-default.row.panel- > table > tbody > tr > td:nth-child(2)').text
            
            mol = Chem.MolFromSmiles(smiles)
            chir = Chem.FindMolChiralCenters(mol,force=True,includeUnassigned=True,useLegacyImplementation=True)

            num_chirals = len(chir)
            chirals = str(chir)

            df2 = pd.DataFrame([[zinc, o, l, m, n, Mwt, logP, Rbonds, Hdonors, Hacceptors, Psa, NetC, num_chirals, chirals, smiles]], columns=['ZINC', 'pdbqt', 'Lowest_E', 'Mean_E', 'Set', 'Mwt', 'LogP', 'Rotatable bonds', 'H-donors', 'H-acceptors', 'PSA', 'Net charge', 'Num chirals', 'chirals', 'SMILES'])

            df2_ls.append(df2)
            
            print (f'[{i}/{leng}]: {zinc} --- Done !')
            # log.write(f'[{i}/{leng}]: {zinc} --- Done !' + '\n')
            i += 1
            
        
        except (AttributeError, TypeError, ValueError):
            print (f'[{i}/{leng}]: {zinc} --- Error !')
            # log.write(f'[{i}/{leng}]: {zinc} --- Error !' + '\n')
            i += 1

    df = pd.concat(df2_ls, ignore_index=True)
        
    return df

# paralleization
def parallelize_dataframe(dt, func):
    num_partitions = num_cores #number of partitions to split dataframe
    df_split = np.array_split(dt, num_partitions)
    pool = multiprocessing.Pool(num_cores)
    dt = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return dt


############################################################

if __name__ == '__main__':

    start = time.time()

    if args.listgen == 'y' and args.dlg2qt == 'n' and args.splitligs == 'n' and args.result2txt == 'n' and args.znparsing == 'n':
        print ('*** Requirements: proteinpath / ligandpath / ligandfmt')
        listgen (
            proteinpath = args.proteinpath, 
            ligandpath=args.ligandpath,
            ligandfmt=args.ligandfmt
        )

    elif args.listgen == 'n' and args.dlg2qt == 'y' and args.splitligs == 'n' and args.result2df == 'n' and args.znparsing == 'n':
        print ('*** Requirement: dlgpath')
        dlg2qt (
            dlgpath = args.dlgpath
        )

    elif args.listgen == 'n' and args.dlg2qt == 'n' and args.splitligs == 'y' and args.result2df == 'n' and args.znparsing == 'n':
        print ('*** Requirement: ligandpath')
        splitligs (
            vinapath = args.vinapath, 
            ligandpath=args.ligandpath
        )

    elif args.listgen == 'n' and args.dlg2qt == 'n' and args.splitligs == 'n' and args.result2df == 'y' and args.znparsing == 'n':
        print ('*** Requirements: dfpath / dlgpath')
        result2df (
            dfpath = args.dfpath, 
            dlgpath = args.dlgpath
        )

    elif args.listgen == 'n' and args.dlg2qt == 'n' and args.splitligs == 'n' and args.result2df == 'n' and args.znparsing == 'y':
        print ('*** Requirements: dfpath / np / dst')

        data_original = pd.read_csv(args.dfpath)

        dt = parallelize_dataframe(data_original, zincparsing)

        if args.rearr == 'n':
            pass

        else:
            dt = pd.concat([dt[dt['Lowest_E'] < 0].sort_values(by=['Lowest_E']), dt[dt['Lowest_E'] >= 0].sort_values(by=['Lowest_E'], ascending=True)])

        dt.reset_index(drop=True).to_csv(f'{args.dst}')

        print ('* Parsing - Done !')


    else:
        print ("*** Please choose one out of [splitligs / listgen / result2df / dlg2qt]")

    
    print ('Elapsed time: ', time.time() - start, 'sec')