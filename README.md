# Command to running GROMACS (Input FIle From CHARMM-GUI)
1. gmx grompp -f step4.0_minimization.mdp -o step4.0_minimization.tpr -c step3_input.gro -r step3_input.gro -p topol.top -n index.ndx -maxwarn -1
2. gmx mdrun -v -deffnm step4.0_minimization
3. gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_input.gro -p topol.top -n index.ndx -maxwarn -1
4. gmx mdrun -v -deffnm step4.1_equilibration
5. export GMX_MAXCONSTRWARN=-1
6. gmx grompp -f step5_production.mdp -o step5_1.tpr -c step4.1_equilibration.gro -p topol.top -n index.ndx -maxwarn -1
7. gmx mdrun -v -deffnm step5_1
