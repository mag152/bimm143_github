\#OUTPUT NEEDS TO BE “gitbub\_document”

\#\#Setup for Docking

We will first prepare our HIV-Pr system for drug docking by making a
protein only PDB format file (i.e. we will remove water, existing ligand
etc.)

1.1 Obtaining and inspecting our input structure.

``` r
#Load the Bio3D package and use the get.pdb() function to download the "1hsg" PDB entry into your RStudio Project directory

library(bio3d)
file <- get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

``` r
#Next use the read.pdb()function to read this PDB file into R so we can prepare it for further analysis.
hiv <- read.pdb(file)

#You can get a quick summary of the contents of this pdb structure object by typing the name ofthe object we just created.
hiv
```

    ## 
    ##  Call:  read.pdb(file = file)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
#Q1: What is the name of the two non protein resid values in this structure? What does resid correspond to and how would you get a listing of all reside values in this structure? 
#Ans: HOH, MK1
```

1.2 Prepare initial protein and ligand input files

Use the atom.select()/trim.pdb() function to make protein-only and
ligand only objects called prot and lig that you can then write out to
new PDB format files in your RStudio project directory.

``` r
prot <- atom.select(hiv, "protein", value = TRUE)
write.pdb(prot, file = "1hsg_protein.pdb")
prot
```

    ## 
    ##  Call:  trim.pdb(pdb = pdb, sele)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1514,  XYZs#: 4542  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 0  (residues: 0)
    ##      Non-protein/nucleic resid values: [ none ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
lig <- atom.select(hiv, "ligand", value=TRUE)
write.pdb(lig, file = "1hsg_ligand.pdb")
lig
```

    ## 
    ##  Call:  trim.pdb(pdb = pdb, sele)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 45,  XYZs#: 135  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 45  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

1.3 Using AutoDockTools to setup protein docking input 1.4 Prepare the
ligand (This is uses AutoDock Tools via MGLTools)

1.5 Prepare a docking configuration file Use to create an input file
first to define protein, ligand, and searc parameters -&gt; this is
stored in config.txt

Section 2: Docking ligands into HIV-1 protease This section uses
Autodock Vina for docking -&gt; creates files all.pdbqt and log.txt

2.3 Inspecting your docking results In order to visualize the docks and
compare to the crystal conformation of the ligand we will process the
all.pdbqt to a PDB format file that can be loaded into VMD. To do this
we will use R and the Bio3D package.

``` r
#Q5: Qualitatively, how good are the docks? Is the crystal binding mode reproduced? Is it the best conformation according to AutoDock Vina?
#The docks look accurate
```

To assess the results quantitatively we will calculate the RMSD (root
mean square distance) between each of the docking results and the known
crystal structure using the bio3d package. Back in RStudio read the
original ligand with added hydrogens that you produced earlier and use
the rmsd() function to compare to your docking results.

``` r
library(bio3d)
res <- read.pdb("all.pdbqt", multi=TRUE)
ori <- read.pdb("ligand.pdbqt")
rmsd(ori,res)
```

    ##  [1] 15.563 16.143 15.935 14.106 13.986 15.875 16.047 16.055 14.707 15.344
    ## [11] 16.261 15.093 16.255 16.469 13.596 13.384 16.531 18.590 14.312 16.738

``` r
write.pdb(res,"results.pdb")
```

Section 4 \#\# Normal Mode analysis for flexibility prediction

``` r
pdb <- read.pdb("1hel")
```

    ##   Note: Accessing on-line PDB file

``` r
modes <- nma( pdb )
```

    ##  Building Hessian...     Done in 0.019 seconds.
    ##  Diagonalizing Hessian...    Done in 0.105 seconds.

``` r
m7 <- mktrj(modes, mode=7, file="mode_7.pdb")
```
