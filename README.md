# tutorial_gromacs
Tutorial Gromacs Protein-Ligand
Semua input file dapat didownload, setelah didownload dapat diextract.

# Solusi Eror Gromacs memakai CHARMM-GUI
1. Tambahkan perintah -maxwarn -1 diterminal tiap equilibrasi dan produksi, contoh: **gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_input.gro -p topol.top -n index.ndx -maxwarn -1** dan **gmx grompp -f step5_production.mdp -o step5_1.tpr -c step4.1_equilibration.gro -p topol.top -n index.ndx -maxwarn -1**
2. Jika terjadi eror berupa pesan Step 36273, time 145.092 (ps)  LINCS WARNING
relative constraint deviation after LINCS tambahakn perintah berikut sebelum menjalankan gromacs diterminal: **export GMX_MAXCONSTRWARN=-1**
