##########fibro  scopen/mytemp
#!/usr/bin/zsh
  
source ~/.zshrc
conda activate  scopen # r-4.0.3

start=$(date +'%s')

time scopen --input "../Fibroblast/filtered_peak_bc_matrix/" --input_format 10X --output_dir "./" --output_prefix "scOpen" --output_format dense --verbose 0 --nc 90 --estimate_rank --no_impute --min_n_components 2 --max_n_components 12 --step_n_components 2

echo "It took $(($(date +'%s') - $start)) seconds"


zsh run_fibroblast.zsh







#############step2 data dimension reduction

import os
import subprocess

input_dir = "../../../ArchR/filtered_peak_bc_matrix/"
output_dir = "../"

job_name = "scOpen"
output_prefix = "scOpen"
subprocess.run(["sbatch", "-J", job_name,
                "-o", f"./cluster_out/{job_name}.txt",
                "-e", f"./cluster_err/{job_name}.txt",
                "--time", "120:00:00",
                "--mem", "180G",
                "-A", "rwth0233",
                "-c", "12",
                "run.zsh", input_dir, output_dir, output_prefix])






#########run.zsh

#!/usr/local_rwth/bin/zsh

source ~/.zshrc
conda activate r-4.0.3

start=$(date +'%s')

time scopen --input $1 --input_format 10X --output_dir $2 --output_prefix $3 --output_format dense --verbose 0 --nc 1 --estimate_rank --no_impute --min_n_components 10 --max_n_components 40 --step_n_components 2

echo "It took $(($(date +'%s') - $start)) seconds"

