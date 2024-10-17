#!/bin/bash

# Warning about the first step
echo "Center the trajectory and remove PBC artifacts, select 1 for protein and 0 for the system."

# Step 1: Center the trajectory and remove PBC artifacts
gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -n index.ndx -center -pbc mol -ur compact

echo "RMSD analysis, select 4 for backbone and 13 for ligand, depending on the analysis being performed."

# Step 2: RMSD analysis
gmx rms -s md.tpr -f md_center.xtc -o rmsd_pro.xvg -tu ns
export LC_NUMERIC="en_US.UTF-8"

# Step 3: Plot RMSD using xmgrace
xmgrace rmsd_pro.xvg

echo "RMSF analysis, select c-alpha (3), depending on the desired analysis."

# Step 4: RMSF analysis
gmx rmsf -f md_center.xtc -s md.tpr -o rmsf.xvg -res

# Step 5: Plot RMSF using xmgrace
xmgrace rmsf.xvg

echo "SASA (Solvent Accessible Surface Area) analysis, select protein (1) or depending on the desired analysis."

# Step 6: SASA analysis
gmx sasa -s md.tpr -f md_center.xtc -o sasa.xvg -tu ns

# Step 7: Plot SASA using xmgrace
xmgrace sasa.xvg

echo "Radius of gyration analysis, select main chain (5), this is flexible too."

# Step 8: Radius of gyration analysis
gmx gyrate -s md.tpr -f md_center.xtc -o gyration.xvg

# Step 9: Modify the first column of gyration.xvg (divide by 1000)
if [ -f "gyration.xvg" ]; then
    awk '
    /^#/ || /^@/ {print $0; next}
    {
        $1 = $1 / 1000
        printf "%.6f ", $1
        for (i = 2; i <= NF; i++) {
            printf "%s ", $i
        }
        print ""
    }
    ' "gyration.xvg" > modified_gyration.xvg
    
    # Replace "ps" with "ns" in the modified file
    sed -i 's/ps/ns/g' modified_gyration.xvg

    echo "The first column of gyration.xvg has been modified and 'ps' has been replaced by 'ns' in modified_gyration.xvg."
else
    echo "gyration.xvg file not found!"
fi

# Step 10: Plot the modified gyration using xmgrace
xmgrace modified_gyration.xvg

echo "Hydrogen bond analysis, select protein (1) and ligand (13)."

# Step 11: Hydrogen bond analysis
gmx hbond -s md.tpr -f md_center.xtc -num hbond.xvg -tu ns 

# Step 12: Plot hydrogen bond data using xmgrace
xmgrace hbond.xvg

# End message and script credits
echo "Analysis completed successfully."
echo "Script written by Purnawan Pontana Putra, Email: purnawanpp@phar.unand.ac.id"
echo "Please cite this code as: Putra PP, Indradi RB, Yuniarta TA, Hanifa D, Syaban MFR, Suharti N. Computational investigation of Pluchea indica mechanism targeting peroxisome proliferator-activated receptor gamma. J Herbmed Pharmacol. 2024;13(4):630-639. doi: 10.34172/jhp.2024.52544."

