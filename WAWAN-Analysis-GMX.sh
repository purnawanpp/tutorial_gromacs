#! /bin/bash

echo "Center the trajectory and remove PBC artifacts choose 1 protein and 0 sistem."

# Step 1: Center the trajectory and remove PBC artifacts choose 1 protein dan 0 sistem
gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -n index.ndx -center -pbc mol -ur compact
echo "RMSD analysis, choose backbone 4 and ligand 13 analysis, this depends on what will be analyzed later."

# Step 2: RMSD analysis, Select backbone analysis 4 and ligand 13, this depends on what will be analyzed later
gmx rms -s md.tpr -f md_center.xtc -o rmsd_pro.xvg -tu ns

# Step 3: Plot RMSD using xmgrace 
xmgrace rmsd_pro.xvg
echo "RMSF analysis, choose c-alfa 3, this depends on what will be analyzed later."

# Step 4: RMSF analysis, choose c-alfa 3, this is flexible according to what will be analyzed
gmx rmsf -f md_center.xtc -s md.tpr -o rmsf.xvg -res

# Step 5: Plot RMSF using xmgrace
xmgrace rmsf.xvg
echo "SASA (Solvent Accessible Surface Area) analysis, choose protein 1 or depending on what you want it can be backbone, or c alpha."

# Step 6: SASA (Solvent Accessible Surface Area) analysis, Choose protein 1 or flexible as desired, bacbone or c alpha
gmx sasa -s md.tpr -f md_center.xtc -o sasa.xvg -tu ns

# Step 7: Plot SASA using xmgrace
xmgrace sasa.xvg
echo "Radius of gyration analysis, choose main chain 5, this is flexible too."

# Step 8: Radius of gyration analysis, Choose main chain 5, it's flexible too
gmx gyrate -s md.tpr -f md_center.xtc -o gyration.xvg


# Step 9: Plot gyration using xmgrace
xmgrace gyration.xvg
echo "Hydrogen bond analysis, hydrogen bond for protein 1 and ligand 13 analysis."


# Step 10: Hydrogen bond analysis, hydrogen bond for protein 1 and ligand 13 analysis
gmx hbond -s md.tpr -f md_center.xtc -num hbond.xvg -tu ns 


# Step 11: Plot hydrogen bond data using xmgrace
xmgrace hbond.xvg

# End message and script credits
echo "Analysis completed successfully."
echo "Script written by Purnawan Pontana Putra, Email: purnawanpp@phar.unand.ac.id"

