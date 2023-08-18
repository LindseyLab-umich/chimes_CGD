#!/bin/bash

tag_list=""

# Figure out what frames we have to analyze

for i in `ls *.4b_clu-s.hist`
do
    tag=${i%%-*}
    taglist="${taglist} $tag"
    
done

tagarr=($taglist)

echo "Compiling the code..."
g++ -O3 -o CluUQ_similarity CluUQ_similarity.cpp


echo "Running the tasks"

ntags=${#tagarr[@]}

rm -f dist_matrix.dat

itag=""    
jtag=""

prefix=""

time for (( i=0; i<$ntags; i++ ))
do

    itag=`printf "%05d" "$i"`
    prefix="./"

    echo "Running loop over i= $i"
    for (( j=i+1; j<ntags; j++ ))
    do

        jtag=`printf "%05d" "$j"`
        
        # Check to make sure all the histograms were generated - this could be an issue if the gen_histograms.sh's *.cmd's timed out before all the jobs finished

        
        if [ ! -e ${prefix}${itag}-${itag}.2b_clu-s.hist ] ; then echo "WARNING: File ${prefix}${itag}-${itag}.2b_clu-s.hist does not exist!"; fi
        if [ ! -e ${prefix}${itag}-${itag}.3b_clu-s.hist ] ; then echo "WARNING: File ${prefix}${itag}-${itag}.3b_clu-s.hist does not exist!"; fi
        if [ ! -e ${prefix}${itag}-${itag}.4b_clu-s.hist ] ; then echo "WARNING: File ${prefix}${itag}-${itag}.4b_clu-s.hist does not exist!"; fi
        
        if [ ! -e ${jtag}-${jtag}.2b_clu-s.hist ] ; then echo "WARNING: File ${jtag}-${jtag}.2b_clu-s.hist does not exist!"; fi
        if [ ! -e ${jtag}-${jtag}.3b_clu-s.hist ] ; then echo "WARNING: File ${jtag}-${jtag}.3b_clu-s.hist does not exist!"; fi
        if [ ! -e ${jtag}-${jtag}.4b_clu-s.hist ] ; then echo "WARNING: File ${jtag}-${jtag}.4b_clu-s.hist does not exist!"; fi
        
        # Calculate the dissimilarities
            

        del_2b=`./CluUQ_similarity ${prefix}${itag}-${itag}.2b_clu-s.hist ${jtag}-${jtag}.2b_clu-s.hist`
        del_3b=`./CluUQ_similarity ${prefix}${itag}-${itag}.3b_clu-s.hist ${jtag}-${jtag}.3b_clu-s.hist`
        del_4b=`./CluUQ_similarity ${prefix}${itag}-${itag}.4b_clu-s.hist ${jtag}-${jtag}.4b_clu-s.hist`
        
        del_all=`echo "$del_2b + $del_3b + $del_4b" | bc -l`
        
        echo "$i $j $del_2b $del_3b $del_4b $del_all" >> dist_matrix.dat
    done
    
    # Needed to make the file gnuplot heatmap compatible
    
    echo "" >> dist_matrix.dat
done
