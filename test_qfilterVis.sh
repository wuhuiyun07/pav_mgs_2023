#!/bin/bash
module load r/4.3.2/gcc-9.3.0
        # cd report/fastp

        echo -e "\nGenerating quality filtering results file qfilter.stats: ... "
        for folder in report/fastp/*/;do
            for file in $folder*.json;do
                ID=$(echo $file|sed 's|/.*$||g')
                readsBF=$(head -n 25 $file|grep total_reads|cut -d ':' -f2|sed 's/,//g'|head -n 1)
                readsAF=$(head -n 25 $file|grep total_reads|cut -d ':' -f2|sed 's/,//g'|tail -n 1)
                basesBF=$(head -n 25 $file|grep total_bases|cut -d ':' -f2|sed 's/,//g'|head -n 1)
                basesAF=$(head -n 25 $file|grep total_bases|cut -d ':' -f2|sed 's/,//g'|tail -n 1)
                q20BF=$(head -n 25 $file|grep q20_rate|cut -d ':' -f2|sed 's/,//g'|head -n 1)
                q20AF=$(head -n 25 $file|grep q20_rate|cut -d ':' -f2|sed 's/,//g'|tail -n 1)
                q30BF=$(head -n 25 $file|grep q30_rate|cut -d ':' -f2|sed 's/,//g'|head -n 1)
                q30AF=$(head -n 25 $file|grep q30_rate|cut -d ':' -f2|sed 's/,//g'|tail -n 1)
                percent=$(awk -v RBF="$readsBF" -v RAF="$readsAF" 'BEGIN{{print RAF/RBF}}' )
                echo "$ID $readsBF $readsAF $basesBF $basesAF $percent $q20BF $q20AF $q30BF $q30AF" >> qfilter.stats
                echo "Sample $ID retained $percent * 100 % of reads ... "
            done
        done

        echo "Done summarizing quality filtering results ... \nMoving to /stats/ folder and running plotting script ... "
        # mv qfilter.stats {config[path][root]}/{config[folder][stats]}
        # cd {config[path][root]}/{config[folder][stats]}

        Rscript workflow/scripts/qfilterVis.R
        echo "Done. "
        rm Rplots.pdf
