# Command to running GROMACS (Input File From CHARMM-GUI)
1. gmx grompp -f step4.0_minimization.mdp -o step4.0_minimization.tpr -c step3_input.gro -r step3_input.gro -p topol.top -n index.ndx -maxwarn -1
2. gmx mdrun -v -deffnm step4.0_minimization
3. gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_input.gro -p topol.top -n index.ndx -maxwarn -1
4. gmx mdrun -v -deffnm step4.1_equilibration
5. export GMX_MAXCONSTRWARN=-1
6. gmx grompp -f step5_production.mdp -o step5_1.tpr -c step4.1_equilibration.gro -p topol.top -n index.ndx -maxwarn -1
7. gmx mdrun -v -deffnm step5_1


# If an error occurs within a certain time, it can be continued without repeating it from the beginning
1. gmx mdrun -s step5_1.tpr -cpi step5_1.cpt -append -deffnm step5_1 -nb GPU

# Remove PBC with no Jump
1. gmx trjconv -s md.tpr -f md.xtc -pbc nojump -o analisis.xtc


# MMGBSA Calculation
![Screenshot_1](https://github.com/purnawanpp/tutorial_gromacs/assets/77323253/94249ebe-ca27-4064-b746-cdb02b73fd57)

# MMGBSA Calculation
1. gmx_MMPBSA -O -i mmpbsa.in -cs md.tpr -ci index.ndx -cg 1 13 -ct md_noPBC.xtc -cp topol.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

# PLOT PCA
1. gmx covar -s step5_1.tpr -f analisis.xtc -o eigenvalues.xvg -v eigenvectors.trr
2. Select protein 1
3. Select ligand 15 (check code your ligand)
4. gmx anaeig -s step5_1.tpr -f analisis.xtc -v eigenvectors.trr -first 1 -last 3 -tu ns -proj projection.xvg
5. gmx anaeig -v eigenvectors.trr -f analisis.xtc -s step5_1.tpr -n new_index.ndx -comp comp.xvg -rmsf eigrmsf.xvg -2d 2d.xvg -b 1 -tu ns -first 1 -last 3
