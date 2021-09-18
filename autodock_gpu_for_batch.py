import os
import argparse



parser = argparse.ArgumentParser (description="Tools for high-throughput docking screening using autodock-gpu")

# General arguments
parser.add_argument("--proteinpath", required=False, default=None, help="Set the path to maps.fld file for receptor")
parser.add_argument("--ligandpath", required=False, default=None, help="Set the path to directory containing ligands")
parser.add_argument("--dlgpath", required=False, default=None, help="Path to dlg files")

# Default setting
parser.add_argument("--ligandfmt", required=False, default="pdbqt", help="Format of the ligand files")
parser.add_argument("--vinapath", required=False, default="/opt/vina/", help="Path to autodock_vina")
parser.add_argument("--txtpath", required=False, default="./", help="Path to the directory for result.txt")
parser.add_argument("--listpath", required=False, default="./", help="Path to list.txt.txt")

# Select 
parser.add_argument("--listgen", required=False, default="n", help="Generate list for docking (y/n)")
parser.add_argument("--dlg2qt", required=False, default="n", help="Dlg-to-pdbqt convert only (y/n)")
parser.add_argument("--splitligs", required=False, default="n", help="Ligand split (y/n)")
parser.add_argument("--result2txt", required=False, default="n", help="Merge data (y/n)")

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




def result2txt (txtpath = args.txtpath, dlgpath = args.dlgpath):

    path = dlgpath
    file_list = os.listdir(path)
    file_list_outputs = [file for file in file_list if file.endswith(".xml")]

    if txtpath != "./":
        os.mkdir(f"{txtpath}")
    
    else:
        pass

    with open (f"{txtpath}/result_merged.txt", "wt") as file:

         for ligs in file_list_outputs:
            with open (f"./{ligs}", "r") as data:
                l = ligs.replace("xml", "pdbqt")

                lines = data.readlines()[3]
                lines1 = lines.rstrip("\n")
                
                file.write(l)
                file.write("-----")
                file.writelines(lines1[28:57])
                file.write("-----")
                file.writelines(lines1[66:-19])

                k = ligs.replace("xml", "dlg")
            with open (f"./{k}", "r") as data2:
                
                lines2 = data2.readlines()[78]

                file.write("-----")
                file.writelines(lines2[35:])


    print ("* Result merge - Done !")



############################################################

if __name__ == '__main__':
    if args.listgen == 'y' and args.dlg2qt == 'n' and args.splitligs == 'n' and args.result2txt == 'n':
        print ('*** Requirements: proteinpath / ligandpath / ligandfmt')
        listgen (
            proteinpath = args.proteinpath, 
            ligandpath=args.ligandpath,
            ligandfmt=args.ligandfmt
        )

    elif args.listgen == 'n' and args.dlg2qt == 'y' and args.splitligs == 'n' and args.result2txt == 'n':
        print ('*** Requirement: dlgpath')
        dlg2qt (
            dlgpath = args.dlgpath
        )

    elif args.listgen == 'n' and args.dlg2qt == 'n' and args.splitligs == 'y' and args.result2txt == 'n':
        print ('*** Requirement: ligandpath')
        splitligs (
            vinapath = args.vinapath, 
            ligandpath=args.ligandpath
        )

    elif args.listgen == 'n' and args.dlg2qt == 'n' and args.splitligs == 'n' and args.result2txt == 'y':
        print ('*** Requiremenst: txtpath / dlgpath')
        result2txt (
            txtpath = args.txtpath, 
            dlgpath = args.dlgpath
        )

    else:
        print ("*** Please choose one out of [splitligs / listgen / result2txt / dlg2qt]")