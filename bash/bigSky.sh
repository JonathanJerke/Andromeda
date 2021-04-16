#!/bin/bash
echo "#!/bin/bash"
echo "#----------------------------------------------------"
echo "# Example SLURM job script to run pure OpenMP applications"
echo "#----------------------------------------------------"
echo "#SBATCH -J omp_job        # Job name"
echo "#SBATCH -o omp_job.o%j    # Name of stdout output file"
echo "#SBATCH -e omp_job.o%j    # Name of stderr output file"
echo "#SBATCH -p normal  # Queue name"

if [ -e hosts ] 
then
	let hosts=`cat hosts`
else
	echo "set hosts for speedup"
	let hosts=1
fi
let MPI=48*hosts
echo "#SBATCH -N $hosts              # Total number of nodes requested"

echo "#SBATCH -n $MPI	          # Total number of mpi tasks requested"
echo "#SBATCH -t $1:00:00       # Run time (hh:mm:ss) - 1.5 hours"
echo "# The next line is required if the user has more than one project"
echo "#SBATCH -A ScalIT  # Project/allocation number"
echo "# This example will run on 1 node with 24 OpenMP threads"
echo "# Please do set the number of threads by yourself!"
echo "export OMP_NUM_THREADS=48"
echo "module load intel launcher hdf5 ooops"
echo "# Launch the pure OpenMP application directly"
echo "$LAUNCH/../bash/launcher.sh $2 $3"
