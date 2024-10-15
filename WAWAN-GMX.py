import os
import subprocess

def run_command(command):
    """Run a shell command and wait for it to finish."""
    process = subprocess.Popen(command, shell=True)
    process.wait()

# Loop over all directories
for dir in os.listdir():
    if os.path.isdir(dir):
        # Change to directory
        os.chdir(dir)
        
        # Process ligand files
        ligand_pdb = "lig*.pdb"
        run_command(f"cp {ligand_pdb} LIG.pdb")
        run_command("acpype -i LIG.pdb -c gas")
        
        # Process protein files
        protein_pdb = "protein*.pdb"
        run_command(f"gmx pdb2gmx -f {protein_pdb} -o protein.pdb -ignh")
        
        # Create complex.pdb by combining protein and ligand
        run_command("grep -h ATOM protein.pdb LIG.acpype/LIG_NEW.pdb >| complex.pdb")
        run_command("cp LIG.acpype/LIG_GMX.itp LIG.itp")
        
        # Edit topol.top
        with open("topol.top", "r") as file:
            topol_content = file.read()
        
        topol_content = topol_content.replace('forcefield.itp"', '''forcefield.itp"\n#include "LIG.itp"\n\n; Ligand position restraints\n#ifdef POSRES\n#include "LIG-posre.itp"\n#endif\n''')
        
        with open("topol2.top", "w") as file:
            file.write(topol_content)
        
        run_command("mv topol2.top topol.top")
        run_command("echo 'LIG   1' >> topol.top")
        
        # Set up the box
        run_command("gmx editconf -f complex.pdb -o box.pdb -bt dodecahedron -d 1.0 -c")
        run_command("gmx solvate -cp box.pdb -cs spc216.gro -p topol.top -o solv.pdb")
        
        # Ion addition
        run_command("gmx grompp -f ../ions.mdp -c solv.pdb -p topol.top -o ions.tpr -maxwarn 1")
        run_command("gmx genion -s ions.tpr -o ions.pdb -p topol.top -pname NA -nname CL -neutral")
        
        # Energy minimization
        run_command("gmx grompp -f ../em.mdp -c ions.pdb -p topol.top -o em.tpr -maxwarn 1")
        run_command("gmx mdrun -v -deffnm em")
        run_command("gmx energy -f em.edr -o potential.xvg")
        
        # Create index for ligand
        run_command("gmx make_ndx -f LIG.acpype/LIG_NEW.pdb -o LIG-index.ndx")
        
        # Generate ligand position restraints
        run_command("gmx genrestr -f LIG.acpype/LIG_NEW.pdb -n LIG-index.ndx -o LIG-posre.itp -fc 1000 1000 1000")
        
        # Create general index
        run_command("gmx make_ndx -f em.gro -o index.ndx")
        
        # NVT equilibration
        run_command("gmx grompp -f ../nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 1")
        run_command("gmx mdrun -v -s nvt.tpr -deffnm nvt")
        
        # NPT equilibration
        run_command("gmx grompp -f ../npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 1")
        run_command("gmx mdrun -v -s npt.tpr -deffnm npt")
        
        # Molecular dynamics run
        run_command("gmx grompp -f ../md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr -maxwarn 1")
        run_command("gmx mdrun -v -s md.tpr -deffnm md")
        
        # Move back to the parent directory
        os.chdir("..")

print("Successfully completed all steps.")
print("Code compiled by Purnawan Pontana Putra, Email: purnawanpp@phar.unand.ac.id")
print ("Please cite this code as: Putra PP, Indradi RB, Yuniarta TA, Hanifa D, Syaban MFR, Suharti N. Computational investigation of Pluchea indica mechanism targeting peroxisome proliferator-activated receptor gamma. J Herbmed Pharmacol. 2024;13(4):630-639. doi: 10.34172/jhp.2024.52544.")

