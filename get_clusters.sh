#!/bin/bash


# If running on LC, execute with: ./get_clusters.sh LLNL-LC /path/to/training_data.xyzf
# Otherwise, if on UM-ARC, run with: ./get_clusters.sh /path/to/training_data.xyzf
#
# Dependencies:
# - helpers.py
# - extract_clusters.cpp
# - training_data.xyzf

################################################################3

# Total frames in the training set:

nf=1116
traj=$2 "/p/lustre2/rlindsey/cluUQ/remaining_frames//training_data.xyzf"


# For LC:

if [ "$1" == "LLNL-LC" ]] ; then 
    module load python/3.10.8
fi

mpiicc -O3 -o extract_clusters extract_clusters.cpp

    

if [ 2 -eq 1 ] ; then

    cp $traj training_data.xyzf
    #cp ~/Codes/al_driver-myLLfork/src/helpers.py . -- manually copied it here
    python -c "import helpers; helpers.break_apart_xyz($nf,\"training_data.xyzf\")"
    rm -f *FORCES*
fi 

# Takes about 10s/frame

curr=1

echo "Building the submit script and submitting..."

jobs_per_block=10 # How many jobs per block
task=0              # Current job
block=-1            # Current block
prev_block=-2       # Previous block


if [ 1 -eq 1 ] ; then
    time for i in `ls training_data_*xyzf`
    do

        tag=${i%*.xyzf}
        tag=${tag#*#}
    
    
        # make sure there weren't any wierd formatting issues
        
        if [ -e ${tag}.2b_clu-s.txt ] && [ -e ${tag}.3b_clu-s.txt ] && [ -e ${tag}.4b_clu-s.txt ]; then
        
            problem_2b=`awk '{print NF}' ${tag}.2b_clu-s.txt | sort | uniq -c | wc -l` 
            problem_3b=`awk '{print NF}' ${tag}.3b_clu-s.txt | sort | uniq -c | wc -l` 
            problem_4b=`awk '{print NF}' ${tag}.4b_clu-s.txt | sort | uniq -c | wc -l` 
            
            problematic=`echo "$problem_2b + $problem_3b + $problem_4b" | bc`
            
            if [ $problematic -gt 3 ] ; then
            
                echo "DETECTED PROBLEMATIC CLUSTER: ${tag}.*b_clu-s.txt " 
                echo "2B clusters:"; awk '{print NF}' ${tag}.2b_clu-s.txt | sort | uniq -c 
                echo "3B clusters:"; awk '{print NF}' ${tag}.3b_clu-s.txt | sort | uniq -c 
                echo "4B clusters:"; awk '{print NF}' ${tag}.4b_clu-s.txt | sort | uniq -c 
                rm -f ${tag}.*b_clu-s.txt
                echo " ====== "
            else
                echo "Skipping: $tag"
                continue
            fi
        fi
    
        echo "Working on tag: $tag"

        if [ `echo "${task} % ${jobs_per_block}" | bc`  == 0 ] ; then
    
            let block=block+1
            let prev_block=prev_block+1
        
            
            if [ $prev_block -ge 0 ]; then      
                echo "  ... Finished writing block ${prev_block} cmd file ... submitting!"
                sbatch run-partition-${prev_block}.cmd
            fi
            
            # Create the .cmd script
            
            rm -f  run-partition-${block}.cmd         
            echo "#!/bin/bash                                 "  >> run-partition-${block}.cmd
            echo "                                            "  >> run-partition-${block}.cmd
            echo "#SBATCH -J block-${block}                   "  >> run-partition-${block}.cmd
            echo "#SBATCH --nodes 1                           "  >> run-partition-${block}.cmd         
            echo "#SBATCH --ntasks-per-node 36                "  >> run-partition-${block}.cmd         
            echo "#SBATCH -t 0:30:00                          "  >> run-partition-${block}.cmd         



            if [ "$1" == "LLNL-LC" ]] ; then 
                echo "#SBATCH -p pdebug                           "  >> run-partition-${block}.cmd       
                echo "#SBATCH -A pbronze                          "  >> run-partition-${block}.cmd        
            else
                echo "#SBATCH -p standard                         "  >> run-partition-${block}.cmd        
                echo "#SBATCH -A rklinds1                     "  >> run-partition-${block}.cmd           
            fi

            echo "#SBATCH -V                                  "  >> run-partition-${block}.cmd
            echo "#SBATCH -o stdoutmsg                        "  >> run-partition-${block}.cmd          

            if [ "$1" != "LLNL-LC" ]] ; then 
                echo "source /home/rklinds/codes/chimes_lsq-myLLfork/modfiles/UM-ARC.mod "  >> run-partition-${block}.cmd    
            fi
    	
            echo "rm -f run-${block}.log                                             "  >> run-partition-${block}.cmd    
        fi
     
        echo "time srun -n 36 ./extract_clusters $i $block >> run-${block}.log" >> run-partition-${block}.cmd 

        echo "rm -f tmp-block-${block}; for m in {0..35}; do cat 2b_clu-r-rank-\${m}-block-${block}.txt >> tmp-block-${block}; done; mv tmp-block-${block} ${tag}.2b_clu-r.txt" >> run-partition-${block}.cmd
        echo "rm -f tmp-block-${block}; for m in {0..35}; do cat 3b_clu-r-rank-\${m}-block-${block}.txt >> tmp-block-${block}; done; mv tmp-block-${block} ${tag}.3b_clu-r.txt" >> run-partition-${block}.cmd
        echo "rm -f tmp-block-${block}; for m in {0..35}; do cat 4b_clu-r-rank-\${m}-block-${block}.txt >> tmp-block-${block}; done; mv tmp-block-${block} ${tag}.4b_clu-r.txt" >> run-partition-${block}.cmd
    
        echo "rm -f tmp-block-${block}; for m in {0..35}; do cat 2b_clu-s-rank-\${m}-block-${block}.txt >> tmp-block-${block}; done; mv tmp-block-${block} ${tag}.2b_clu-s.txt" >> run-partition-${block}.cmd
        echo "rm -f tmp-block-${block}; for m in {0..35}; do cat 3b_clu-s-rank-\${m}-block-${block}.txt >> tmp-block-${block}; done; mv tmp-block-${block} ${tag}.3b_clu-s.txt" >> run-partition-${block}.cmd
        echo "rm -f tmp-block-${block}; for m in {0..35}; do cat 4b_clu-s-rank-\${m}-block-${block}.txt >> tmp-block-${block}; done; mv tmp-block-${block} ${tag}.4b_clu-s.txt" >> run-partition-${block}.cmd
    
        echo " rm -f *b_clu-*-rank-*-block-${block}.txt tmp-block-${block}" >> run-partition-${block}.cmd
    
        echo "" >> run-partition-${block}.cmd
    
        let curr=curr+1  
        let task=task+1

    done
 
    echo "  .... Finished writing the last block ${block} cmd file ... submitting!"
    sbatch run-partition-${block}.cmd

fi
