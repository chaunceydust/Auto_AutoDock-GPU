# ***Auto AutoDock-GPU***

This script helps preparing input files for protein-ligand docking using AutoDock-GPU in batch mode, as well as merging results.

*[AutoDock-GPU GitHub](https://github.com/ccsb-scripps/AutoDock-GPU)*

Before start, the Anaconda (Miniconda) and Nvidia-docker must be installed in your system.

In addition, maps.fld formatted file of protein structure have to be prepared.

The maps.fld file is able to be prepared using MGLtools and Autogrid softwares.

Please refer to the explanation below.

- - -

This script has four tools, `listgen` / `dlg2qt` / `splitligs` / `result2txt`, for docking.

<br/>

`listgen`: Generate the list.txt file which is essential for the running of AutoDock-GPU for batch mode.

`dlg2qt`: Convert the dlg formatted files (result files) to the pdbqt files that enable to open in pymol.

`splitligs`: Split the merged ligand files downloaded from the ZINC database into each file.

`result2txt`: Organize the highest score among each result into a text file.

- - -

<br/>

## ***How to use***
<br/>

### **1. Pull the docker image**

    docker pull jongseopark/autodock_gpu_js

<br/>

### **2. Prepare several files**

First, create a environment of Anaconda and install **mgltools** and **autodock** as follows.

    conda create -n autodock

    conda activate autodock

    conda install -c bioconda mgltools
    conda install -c hcc autodock


Now, when you type **pmv** in the terminal, the GUI program will be opened.


![pmv](https://github.com/jongseo-park/Auto_AutoDock-GPU/blob/master/images/a.png)


You can create a pdbqt formatted file of protein and gpf formatted file containing several parameters for docking grid. Before you begin, I recommend deleting all waters in the pdb file.

    Edit > Hydrogen > Add > OK (All hydrogens / noBondOrder / yes)

    Click Autodock icon 

    Grid > Macromolecule > Choose > Select Molecule with a selection > ignore an error (click 'OK')

Then, the new window for saving of pdbqt file will be opened. Please save it to the working directory.

![Save_PDBQT](https://github.com/jongseo-park/Auto_AutoDock-GPU/blob/master/images/c.gif)

<br/>

In that programv (pmv), you are also able to obtain the gpf formatted file as follows.

    Grid > Grid box > set grid you want
    File > Close saving current
    Grid > Output > Save GPF

Please refer the gif image below. Then, you can save the gpf file for the grid you want.

![Save_GPF](https://github.com/jongseo-park/Auto_AutoDock-GPU/blob/master/images/b.gif)

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

and in the ***'ligand'*** directory, the ligand files formatted in pdbqt are located.


    Working_dir/

        autodock_gpu_for_batch.py

        result/
            # empty_directory         

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

    docker run -it --gpus all -v "/path/to/working/directory/:/home/run/" jongseopark/autodock_gpu_js:v01

<br/>

### **5. Run the script**

Run the script with various arguments.


You have to choose one argument out of `splitligs` / `listgen` / `result2txt` / `dlg2qt`.

The workflow is ...

1. If you need, run the `splitligs`.
2. Generate the list.txt file for running AutoDock-GPU in batch mode using `listgen`.
3. Run the AutoDock-GPU in batch mode.
4. Organize the highest score among each result into a text file using `result2txt`.
If you need, you can convert the .dlg files to .pdbqt files using `dlg2qt`.



<br/>

You can use `splitligs` tools as follows.

    python3 autodock_gpu_for_batch.py
    --splitligs y
    --ligandpath ./merged_ligs/                # Path to the directory 
                                               # containing merged ligands


You can use `listgen` tools as follows.

    python3 autodock_gpu_for_batch.py
    --listgen y
    --proteinpath ./protein/protein.maps.fld   # Path to the maps.fld file
    --ligandpath ./ligands/                    # Path to the directory containing ligands
    --ligandfmt pdbqt                          # Format of ligand files

For running the AutoDock-GPU in batch mode, use this command in the result directory.

    autodock_gpu_128wi -filelist ../list.txt

After finish the doking, you can use `result2txt` tools to organize results.

    python3 autodock_gpu_for_batch.py
    --result2txt y
    --txtpath ./result/                        # Set the path you want
    --dlgpath ./result/                        # Path to the directory containing dlg files


You can also use `dlg2qt` tools as follows.

    python3 autodock_gpu_for_batch.py
    --dlg2qt y
    --dlgpath ./result/                        # Path to the dlg files (result files)




<br/>

You can check all of the arguments with the explanation using this command.

    python3 autodock_gpu_for_batch.py -h


<br/>


- - -
## ***And ...***

Since I am a structural biologist, but not a computational biologist, this script which I wrote might not be fancy or cause error. So, feel free to contact me, and tell or give advice about the coding.

<br/>
<br/>

*Thanks.*
