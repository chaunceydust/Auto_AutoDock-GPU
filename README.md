# ***AutoDock-GPU_for_ZINC***

This script helps preparing input files for protein-ligand(from ZINC db) docking using AutoDock-GPU in batch mode, as well as merging results and parsing detailed information about ligands.

*[AutoDock-GPU GitHub](https://github.com/ccsb-scripps/AutoDock-GPU)*

<br/>

### Requirements
```
AutoDock-GPU
Conda
autodock-vina
```

or you can use the preset docker image for this repository that is using the nvidia-docker.

```
docker pull jongseopark/autodock_gpu_js:latest
```

### Set conda env. for docking
```
conda env create -f requirements.yml
conda activate autodock_gpu
```

<br>
<br>

## How to use

### 1. Download ligand files from ZINC database

You can download the pdbqt formatted files from the ZINC database *(not described in this repository)*.

*[ZINC database](https://zinc.docking.org)*

(Tranches > 3D > select you want > download as pdbqt)

<br>

### 2. Generate a protein.maps.fld file that containing the information about docking box


In the case of protein files, you can prepare the maps.fld formatted file of protein structure using MGLtools and Autogrid softwares.

Please refer to the explanation below.

<br/>

First, create a environment of conda and install mgltools and autodock as follows.

    conda create -n proteinprep

    conda activate proteinprep

    conda install -c bioconda mgltools
    conda install -c hcc autodock

Now, when you type **pmv** in the terminal, the GUI program will be opened.

![pmv](https://github.com/jongseo-park/AutoDock-GPU_for_ZINC/blob/master/images/a.png)


You can create a pdbqt formatted file of protein and gpf formatted file containing several parameters for docking grid. Before you begin, I recommend deleting all waters in the pdb file.

    Edit > Hydrogen > Add > OK (All hydrogens / noBondOrder / yes)

    Click Autodock icon 

    Grid > Macromolecule > Choose > Select Molecule with a selection > ignore an error (click 'OK')

Then, the new window for saving of pdbqt file will be opened. Please save it to the working directory.

![Save_PDBQT](https://github.com/jongseo-park/AutoDock-GPU_for_ZINC/blob/master/images/c.gif)

<br/>

In that programv (pmv), you are also able to obtain the gpf formatted file as follows.

    Grid > Grid box > set grid you want
    File > Close saving current
    Grid > Output > Save GPF

Please refer the gif image below. Then, you can save the gpf file for the grid you want.

![Save_GPF](https://github.com/jongseo-park/AutoDock-GPU_for_ZINC/blob/master/images/b.gif)

<br/>

The pdbqt and gpf files should be located in the ***same directory***.

<br/>

In the terminal, enter the directory containing pdbqt and gpf files, and use this command.

    autogrid4 -p yourprotein.gpf -l yourprotein.glg

Now, finally, the input files for receptor will be created in the same directory. 

<br/>

### 3. Directory setting

When I run this script, I prepare the directories named `protein` and `ligands`.

In the `protein` directory, the prepared **fld** formatted file with several **map** files are located,

and in the `ligands` directory, the ligand files formatted in pdbqt are located.


    Working_dir/

        autodock_gpu_for_zinc.py          # The script in this repository
        tools /
            ...py
            ...py

        protein/                             
            protein.maps.fld
            protein.A.map
            protein.C.map
            protein.Br.map
            ...

        ligands/
            ligand_1.pdbqt
            ligand_2.pdbqt
            ...
            ligand_n.pdbqt


<br/>

### 4. Run docking

*flags*

```
usage: autodock_gpu_for_zinc.py [-h] [-p PROTEINPATH] [-l LIGANDPATH] [-v VINAPATH] [-ls LISTPATH] [-d RESULTPATH]
                                [-lf LIGANDFORMAT] [-bin AUTODOCKBIN] [--gpu GPU] [-smi] [-sl] [-c] [--np NP]
                                [--fn FN] [--csv CSV] [--splitnum SPLITNUM]




Tools for high-throughput docking screening using autodock-gpu

options:
  -h, --help            show this help message and exit
  -p PROTEINPATH, --proteinpath PROTEINPATH
                        Set the path to maps.fld file for receptor (file)
  -l LIGANDPATH, --ligandpath LIGANDPATH
                        Set the path to directory containing ligands (directory)
  -v VINAPATH, --vinapath VINAPATH
                        Path to autodock_vina (directory)
  -ls LISTPATH, --listpath LISTPATH
                        Path to list.txt (file)
  -d RESULTPATH, --resultpath RESULTPATH
                        Set the path to result directory (directory)
  -lf LIGANDFORMAT, --ligandformat LIGANDFORMAT
                        Format of the ligand files
  -bin AUTODOCKBIN, --autodockbin AUTODOCKBIN
                        Binary file of AutoDock-GPU (bin_file)
  --gpu GPU             Which GPU do you use ? (starts at 1)
  -smi, --smi
  -sl, --splitligs
  -c, --continue
  --np NP               Number of cores for execute
  --fn FN               A function name to use
  --csv CSV             csv file for data parsing from ZINC database
  --splitnum SPLITNUM   How many divided list files you want ?
```

<br>

*Functions*

When you type the command as belows, you can find the list of functions.

```
python3 autodock_gpu_for_zinc.py --fn x
```

<br>

```
Please select/enter the one of functions below and enter it.

            - Essential fxns
            1) splitligs (ligandpath, vinapath) 
            2) listgen (proteinpath, ligandpath, listpath)
            3) run_docking (listpath, autodockbin, resultpath, gpu)
            4) dlg2qt (resultpath)
            5) result2df (resultpath)
            6) znparsing (np, csv)

            - Miscellaneous fxns
            listsplit (listpath, splitnum)
```

This script has six tools, `splitligs` / `listgen` / `run_docking` / `dlg2qt` / `result2df` / `znparsing`, for docking and organizing results.

<br/>

`splitligs`: Split merged ligand files downloaded from the ZINC database into each file.

`listgen`: Generate the list.txt file that is essential for the running of AutoDock-GPU for batch mode.

`run_docking`: Run the AutoDock-GPU

`dlg2qt`: Convert the dlg formatted files (result files) to the pdbqt files that enable to open in molecular visualization programs.

`result2df`: Organize the highest score among each result into a pandas dataframe and save it to csv file.

`znparsing`: Collect the detailed information about ligands from ZINC database using web crawling tool, and organize into the csv file.

You can use these functions one by one when you enter the specific function name using argument `--fn`

If you want to perform protein-ligand docking from ligand split to data arrangement, 

then do not enter the `--fn` argument.

<br>

For example ...

If you want to use only one function such as `splitligs`, not entire process,

```
python3 autodock_gpu_for_zinc.py \
    -fn splitligs
    -l ./ligands
    -v /opt/programs/autodock_vina/bin/
```

or it you want to run the entire process, then use as follows.

```
python3 autodock_gpu_for_zinc.py \
    -p ./protein/protein.maps.fld
    -l ./ligands
    -v /opt/programs/autodock_vina/bin/
    -bin /opt/programs/AutoDock-GPU/bin/autodock_gpu_64wi
    -sl             # when your ligand file is merged version, then use -sl switch.
```

When you enter the command for entire process, you can see the check list in your terminal window.

```
 
* vinapath set to "/opt/programs/autodock_vina/bin/" by user
* listpath automatically set to "./list.txt"
* resultpath automatically set to "./result"
 
Input arguments are correct !
-----------------------------
[proteinpath] ./protein/protein.maps.fld
[ligandpath] ./ligands
[vinapath] /opt/programs/autodock_vina/bin/
[listpath] ./list.txt
[resultpath] ./result
[ligandformat] pdbqt
[autodockbin] /opt/programs/AutoDock-GPU/bin/autodock_gpu_64wi
[gpu] 1
[smi_dest] False
[splitligs_dest] True
[continue_dest] False
[np] 8
[fn] 
[csv] results.csv
[splitnum] 4
-----------------------------
Please check your inputs before continuing !! (y / n): 

```

Before run the script, you can lastly check all of the settings.

If they are correct, enter `y` in the terminal window, then docking will begin.

<br/>

When the job is completed, then you can find the csv file in the working directory

that contains several detailed information such as

`QED` / `Lowest_E` / `LogP` / `Mwt` / `Rotatable bonds` / `H-donors` / `H-acceptors` / `PSA` / `Net charge` / `Chiral centers & nums` / `Smile string` .

<br>

### 5. etc

#### *Manually parallelize*
Since the running of AutoDock-GPU is not parallelized, if you want to use multiple gpus,

first, you have to split your list.txt file through the `listsplit` function (--fn listsplit).

At that time, you can specify a split number using the argument `--splitnum`, then the list.txt file will be divided by the given number.

<br>

For multiple GPUs, you can specify the GPU number to use through the argument `--gpu` when run the AutoDock-GPU using the function `--run_docking`.

Keep in mind that the input number for `--gpu` starts at 1, not 0.

<br>

Collectively ...

1) Through `listsplit` function, you can divide the list.txt file to several list files such as list_1.txt, list_2.txt, etc.

2) After that, by setting the GPU number using `--gpu` argument, you can manually parallelize the docking.

<br>

#### *ZINC parsing*
When you use the function `znparsing`, you can set the number of processes to parsing information through the argument `--np`.

This argument does not mean the number of cores, but it means how many jobs will be running to parsing data from ZINC database.

You can set this value as a high enough number such as 50, 100, ...

- - -