# ***AutoDock-GPU_for_ZINC***

This script helps preparing input files for protein-ligand(from ZINC db) docking using AutoDock-GPU in batch mode, as well as merging results and parsing detailed information about ligands.

*[AutoDock-GPU GitHub](https://github.com/ccsb-scripps/AutoDock-GPU)*

<br/>

Before start, the Anaconda *(Miniconda)* and Nvidia-docker must be installed in your system.

In addition, you have to prepared input files, `protein.maps.fld` and `ligand.pdbqt` files.


In the case of protein files, you can prepare the maps.fld formatted file of protein structure using MGLtools and Autogrid softwares.

Please refer to the explanation below.

<br/>

In the case of ligand files, you can download the pdbqt formatted files from the ZINC database *(not described in this repository)*.

*[ZINC database](https://zinc.docking.org)*

(Tranches > 3D > select you want > download as pdbqt)

- - -

This script has five tools, `splitligs` / `listgen` / `result2df` / `dlg2qt` / `znparsing`, for docking and organizing results.

<br/>

`splitligs`: Split the merged ligand files downloaded from the ZINC database into each file.

`listgen`: Generate the list.txt file which is essential for the running of AutoDock-GPU for batch mode.

`result2df`: Organize the highest score among each result into a pandas dataframe csv file.

`dlg2qt`: Convert the dlg formatted files (result files) to the pdbqt files that enable to open in pymol.

`znparsing`: Collect the detailed information about ligands from ZINC database using web crawling tool, and organize into the csv file.


- - -

<br/>

## ***How to use***
<br/>

### **1. Pull the docker image**

    docker pull jongseopark/autodock_gpu_js:latest

<br/>

### **2. Prepare several files**

First, create a environment of Anaconda and install **mgltools** and **autodock** as follows.

    conda create -n autodock

    conda activate autodock

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

### **3. Directory setting**

When I run this script, I prepare the directories named ***'protein'*** and ***'ligand'***.

In the ***'protein'*** directory, the prepared **'fld'** formatted file with several **'map'** files are located,

and in the ***'ligands'*** directory, the ligand files formatted in pdbqt are located.


    Working_dir/

        result/
            autodock_gpu_for_zinc.py          # The script in this repository

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

        merged_ligs/
            merged_ligands.pdbqt
            ...

<br/>

### **4. Run the docker container from the image with mounting of the working directory**

    docker run -it --gpus all -v "/path/to/working/directory/:/home/run/" jongseopark/autodock_gpu_js:latest

<br/>

### **5. Run the script**

Run the script with various arguments.


You have to choose one argument out of `splitligs` / `listgen` / `result2df` / `dlg2qt` / `znparsing`.

The workflow is ...

1. Run the `splitligs` to split all of the ligand files into each file.
2. Generate the `list.txt` file for running AutoDock-GPU in batch mode using `listgen`.
3. Run the AutoDock-GPU in batch mode (without this python script).
4. Organize the highest score among each result into a text file using `result2df`.
If you need, you can convert the .dlg files to .pdbqt files using `dlg2qt`.
5. Parsing detailed information about ligands using `znparsing`.

<br/>

you can pull the script as follows in the docker container.

    wget https://raw.githubusercontent.com/jongseo-park/AutoDock-GPU_for_batch/master/autodock_gpu_for_zinc.py

<br/>

### **step 1**

When you donwload the pdbqt formatted ligand files from ZINC database,

you can see that all of the files are composed of multiple ligand bundles.

Then, you can use `splitligs` as follows, and the merged ligand files will be splitted into each file.

I have confirmed that `splitligs` works well with the merged pdbqt files,

but in the case of other file formats such as sdf, I cannot guarantee the operation.


    python3 autodock_gpu_for_zinc.py
    --splitligs y
    --ligandpath ./merged_ligs/                # Path to the directory 
                                               # containing merged ligands

<br/>

### **step 2**

The list.txt file is essential file for running the AutoDock-GPU in batch mode.

For creating the list.txt file, use `listgen` as follows.



    python3 autodock_gpu_for_zinc.py
    --listgen y
    --proteinpath ./protein/protein.maps.fld   # Path to the maps.fld file
    --ligandpath ./ligands/                    # Path to the directory containing ligands


The list file is comprised of the path to protein, path to ligand, and the job name.

Generated list.txt file comprises of a line of the path to protein,

and several lines of the path to ligand and job names.

    ./protein/protein.maps.fld                 # Path to protein.maps.fld
    ./ligands/ligand_1.pdbqt                   # Path to ligand
    ligand_1                                   # Job name
    ./ligands/ligand_2.pdbqt                   # Path to ligand
    ligand_2                                   # Job name
    ./ligands/ligand_3.pdbqt                   # ...
    ligand_3                                   # ...
    ...
    ...
    ./ligands/ligand_n.pdbqt
    ligand_n
    
You can also specify the path to list.txt files using `--listpath` argument.

Default path is `./`.

<br/>

### **step 3**

For running the AutoDock-GPU in batch mode, use this command in the result directory.

    autodock_gpu_128wi -filelist /path/to/list
    
    # example
    autodock_gpu_128wi -filelist ../list.txt

Then, the result files (xml , dlg) will be generated in the result directory.


<br/>

### **step 4**

After finish the doking, you can use `result2df` to organize results.

    python3 autodock_gpu_for_zinc.py
    --result2df y
    --dfpath ./result/                         # Set the path you want
    --dlgpath ./result/                        # Path to the directory containing dlg files
    --rearr y                                  # Sort Lowest_binding_energy in ascending order (non-essential)
    --setnum 1                                 # The value that distinguish sets (non-essential)

After some time, you can find a pandas dataframe csv file (result_merged.csv) in the `result` directory declared by `--dfpath` argument.

If you set the `--rearr` argument to `y`, then the *Lowest_binding_energy* will be sorted in ascending order (default = n).

You can use the `--setnum` argument to distinguish several sets. 

The value you set using `--setnum` is added to the last column of the csv file.

If you don't need to distinguish `--setnum`, you don't need to enter it (default = 1).

<br/>

You can also use `dlg2qt` tools as follows, and the pdbqt files (result) will be generated in the `dlgpath`.

    python3 autodock_gpu_for_zinc.py
    --dlg2qt y
    --dlgpath ./result/                        # Path to the dlg files (result files)


### **step 5**

You can collect the detailed information about all of the ligands using `znparsing` tool.

    python3 autodock_gpu_for_zinc.py
    --znparsing y
    --dfpath ./result/result_merged.csv
    --dst ./result/dst/result_parsed.csv
    --rearr y                                  # Sort Lowest_binding_energy in ascending order.
    --np 60

Then, you can obtain one csv file containing the detailed information such as 

`Lowest_E` / `LogP` / `Mwt` / `Rotatable bonds` / `H-donors` / `H-acceptors` / `PSA` / `Net charge` / `Chiral centers & nums`.

You have to specify the path of `result_merged.csv` file generated by the `result2df` tool.

In addition, you can set the path to result file containing parsed data from ZINC database using `--dst` argument.

If you set the `--rearr` argument to `y`, then the *Lowest_E* will be sorted in ascending order (default = n).

The number of cores for running the `znparsing` can be declared using `--np` argument.

I recommend the finding of the optimal np value using small dataset. 


<br/>

You can check all of the arguments with the explanation using this command.

    python3 autodock_gpu_for_zinc.py -h


<br/>

- - -

*Thanks.*
