#! /bin/tcsh -q
#$ -cwd -j y -N bio32 -t 1:10
echo $HOSTNAME
echo $SGE_TASK_ID
sleep $SGE_TASK_ID
python zm_sphere2_refine.py target=blq.dat qmax=0.3 rmax=80 pdb=pdb_model.pdb splat=8 prefix=snbr1.80.$SGE_TASK_ID nmax=16 n_trial=5 > log_nbr1.80.$SGE_TASK_ID
