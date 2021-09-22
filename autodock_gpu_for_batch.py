import os
import argparse
from xml.etree.ElementTree import parse
import pandas as pd
import time


parser = argparse.ArgumentParser (description="Tools for high-throughput docking screening using autodock-gpu")

# General arguments
parser.add_argument("--proteinpath", required=False, default=None, help="Set the path to maps.fld file for receptor")
parser.add_argument("--ligandpath", required=False, default=None, help="Set the path to directory containing ligands")
parser.add_argument("--dlgpath", required=False, default=None, help="Path to dlg files")

# Default setting
parser.add_argument("--ligandfmt", required=False, default="pdbqt", help="Format of the ligand files")
parser.add_argument("--vinapath", required=False, default="/opt/vina/", help="Path to autodock_vina")
parser.add_argument("--dfpath", required=False, default="./", help="Path to the directory for result.df")
parser.add_argument("--listpath", required=False, default="./", help="Path to list.txt.txt")

# Select 
parser.add_argument("--listgen", required=False, default="n", help="Generate list for docking (y/n)")
parser.add_argument("--dlg2qt", required=False, default="n", help="Dlg-to-pdbqt convert only (y/n)")
parser.add_argument("--splitligs", required=False, default="n", help="Ligand split (y/n)")
parser.add_argument("--result2df", required=False, default="n", help="Merge data (y/n)")

args = parser.parse_args()




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




def result2df (dfpath = args.dfpath, dlgpath = args.dlgpath):
    i = 1
    path = dlgpath
    file_list = os.listdir(path)
    file_list_outputs = [file for file in file_list if file.endswith(".xml")]
    leng = len(file_list_outputs)

    if dfpath != "./":
        os.mkdir(f"{dfpath}")
    
    else:
        pass

    df = pd.DataFrame (columns=['file name', 'Lowest_binding_energy', 'Mean_binding_energy', 'ZINC'])
    df2_ls = []

    for ligs in file_list_outputs:
        tree = parse(f'{ligs}')
        root = tree.getroot()
        cluhis = root.findall ('clustering_histogram')


        score = [x.find('cluster').attrib for x in cluhis][0]
        lb = score['lowest_binding_energy']
        mb = score['mean_binding_energy']
        l = ligs.replace("xml", "pdbqt")


        k = ligs.replace("xml", "dlg")
        with open (f"./{k}", "r") as data2:
            
            lines = data2.readlines()[78]
            lines2 = lines[35:]

        # df = df.append(pd.DataFrame([[l, lb, mb, lines2]], columns=['file name', 'Lowest_binding_energy', 'Mean_binding_energy', 'ZINC']), ignore_index=True)

        df2 = pd.DataFrame([[l, lb, mb, lines2]], columns=['file name', 'Lowest_binding_energy', 'Mean_binding_energy', 'ZINC'])
        df2_ls.append(df2)

        i += 1
        ratio = i/leng * 100
        print (ratio, '%')

    df = pd.concat(df2_ls, ignore_index=True)

    df.to_csv(f'./{dfpath}/result_merged.csv')

    print ("* Result merge - Done !")




############################################################

if __name__ == '__main__':

    start = time.time()

    if args.listgen == 'y' and args.dlg2qt == 'n' and args.splitligs == 'n' and args.result2txt == 'n':
        print ('*** Requirements: proteinpath / ligandpath / ligandfmt')
        listgen (
            proteinpath = args.proteinpath, 
            ligandpath=args.ligandpath,
            ligandfmt=args.ligandfmt
        )

    elif args.listgen == 'n' and args.dlg2qt == 'y' and args.splitligs == 'n' and args.result2df == 'n':
        print ('*** Requirement: dlgpath')
        dlg2qt (
            dlgpath = args.dlgpath
        )

    elif args.listgen == 'n' and args.dlg2qt == 'n' and args.splitligs == 'y' and args.result2df == 'n':
        print ('*** Requirement: ligandpath')
        splitligs (
            vinapath = args.vinapath, 
            ligandpath=args.ligandpath
        )

    elif args.listgen == 'n' and args.dlg2qt == 'n' and args.splitligs == 'n' and args.result2df == 'y':
        print ('*** Requirements: dfpath / dlgpath')
        result2df (
            dfpath = args.dfpath, 
            dlgpath = args.dlgpath
        )


    else:
        print ("*** Please choose one out of [splitligs / listgen / result2df / dlg2qt]")

    
    print ('Elapsed time: ', time.time() - start)