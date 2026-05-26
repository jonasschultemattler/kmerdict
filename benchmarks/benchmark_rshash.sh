#!/bin/bash

PROGRAM="../build/source/rshash"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
OUT="progout.txt"
CSV="rshash-results-$today.csv"

run()
{
    f="$1"
    BASENAME=$(echo "$f" | sed 's/\.fa.gz$//')
    k=31    
        if [[ "$BASENAME" == *"bacterial"* ]]; then
          params=( 18 0 0 128 0 0 1  19 0 0 128 0 0 1 )
        elif [[ "$BASENAME" == *"human"* ]]; then
          params=( 17 21 25 64 64 128 3  17 0 0 128 0 0 1 )
        elif [[ "$BASENAME" == *"cod"* ]]; then
          params=( 16 20 24 64 64 128 3  16 0 0 128 0 0 1 )
        elif [[ "$BASENAME" == *"kestrel"* ]]; then
          params=( 17 20 0 64 128 0 2  18 0 0 128 0 0 1 )
        fi
        # if [[ "$BASENAME" == *"bacterial"* ]]; then
        #   params=( 18 0 0 128 0 0 1 )
        # elif [[ "$BASENAME" == *"human"* ]]; then
        #   params=( 17 21 25 64 64 128 3 )
        # elif [[ "$BASENAME" == *"cod"* ]]; then
        #   params=( 16 20 24 64 64 128 3 )
        # elif [[ "$BASENAME" == *"kestrel"* ]]; then
        #   params=( 17 20 0 64 128 0 2 )
        # fi
        # if [[ "$BASENAME" == *"bacterial"* ]]; then
        #   params=( 18 0 0 128 0 0 1 )
        # elif [[ "$BASENAME" == *"human"* ]]; then
        #   params=( 17 0 0 128 0 0 1 )
        # elif [[ "$BASENAME" == *"cod"* ]]; then
        #   params=( 16 0 0 128 0 0 1 )
        # elif [[ "$BASENAME" == *"kestrel"* ]]; then
        #   params=( 18 0 0 128 0 0 1 )
        # fi
        # if [[ "$BASENAME" == *"bacterial"* ]]; then
        #   params=( 18 0 0 128 0 0 1  19 0 0 128 0 0 1 )
        # elif [[ "$BASENAME" == *"human"* ]]; then
        #   params=( 17 19 21 32 32 32 3  17 19 21 64 64 64 3  17 19 21 128 128 128 3  17 19 0 32 32 0 2  17 19 0 64 64 0 2  17 19 0 128 128 0 2  17 0 0 32 0 0 1  17 0 0 64 0 0 1  17 0 0 128 0 0 1 )
        # elif [[ "$BASENAME" == *"cod"* ]]; then
        #   params=( 16 20 24 64 64 128 3  16 0 0 128 0 0 1 )
        # elif [[ "$BASENAME" == *"kestrel"* ]]; then
        #   params=( 17 20 0 64 128 0 2  18 0 0 128 0 0 1 )
        # fi

        for ((i=0; i<${#params[@]}; i+=7)); do
          m1=${params[i]}
          m2=${params[i+1]}
          m3=${params[i+2]}
          t1=${params[i+3]}
          t2=${params[i+4]}
          t3=${params[i+5]}
          l=${params[i+6]}

          echo $BASENAME 
          echo $m1 $t1 $m2 $t2 $m3 $t3
          echo $f >> $LOG
          echo $m1 $t1 $m2 $t2 $m3 $t3 >> $LOG

          # /usr/bin/time -v -o time.txt $PROGRAM build -i "$f" -d "${BASENAME}_.dict" -k $k --level $l --m1 $m1 --m2 $m2 --m3 $m3 --t1 $t1 --t2 $t2 --t3 $t3 > OUT 2>&1
          /usr/bin/time -v -o time.txt $PROGRAM build -i "$f" -d "${BASENAME}_.dict" -k $k --level $l --m1 $m1 --m2 $m2 --m3 $m3 --t1 $t1 --t2 $t2 --t3 $t3 --loc -t 255 > OUT 2>&1

          cat OUT >> $LOG

          file_size=$(stat -c%s "${BASENAME}_.dict")
          buildtime=$(grep "User time" time.txt | awk -F': ' '{print $2}')
          buildmem=$(grep "Maximum resident set size" time.txt | awk -F': ' '{print $2}')
          textlength=$(grep "text length: " OUT | sed -E 's/.*text length: ([0-9]+).*/\1/')
          density_r1=$(awk '/density r1/ {print $3}' OUT | tr -d '%')
          density_r2=$(awk '/density r2/ {print $3}' OUT | tr -d '%')
          density_r3=$(awk '/density r3/ {print $3}' OUT | tr -d '%')
          density_s1=$(awk '/density s1/ {print $3}' OUT | tr -d '%')
          density_s2=$(awk '/density s2/ {print $3}' OUT | tr -d '%')
          density_s3=$(awk '/density s3/ {print $3}' OUT | tr -d '%')
          density_ht=$(grep "no kmers HT:" OUT | awk '{print $(NF)}' | sed 's/%//')
          no_kmers=$(awk '/textkmers/ {print $2}' OUT)
          no_minimiser1=$(awk '/no minimiser1/ {print $3}' OUT)
          no_minimiser2=$(awk '/no minimiser2/ {print $3}' OUT)
          no_minimiser3=$(awk '/no minimiser3/ {print $3}' OUT)
          no_distinct_minimiser1=$(awk '/no distinct minimiser1/ {print $4}' OUT)
          no_distinct_minimiser2=$(awk '/no distinct minimiser2/ {print $4}' OUT)
          no_distinct_minimiser3=$(awk '/no distinct minimiser3/ {print $4}' OUT)
          spaceoffsets1=$(grep '^offsets1:' OUT | cut -d':' -f2 | xargs)
          spaceoffsets2=$(grep '^offsets2:' OUT | cut -d':' -f2 | xargs)
          spaceoffsets3=$(grep '^offsets3:' OUT | cut -d':' -f2 | xargs)
          spaceht=$(grep '^Hashtable:' OUT | cut -d':' -f2 | xargs)
          spacer1=$(grep '^R_1:' OUT | cut -d':' -f2 | xargs)
          spacer2=$(grep '^R_2:' OUT | cut -d':' -f2 | xargs)
          spacer3=$(grep '^R_3:' OUT | cut -d':' -f2 | xargs)
          spaces1=$(grep '^S_1:' OUT | cut -d':' -f2 | xargs)
          spaces2=$(grep '^S_2:' OUT | cut -d':' -f2 | xargs)
          spaces3=$(grep '^S_3:' OUT | cut -d':' -f2 | xargs)
          spacetotal=$(grep '^total:' OUT | cut -d':' -f2 | xargs)
          memtotal=$(grep '^total:' OUT | cut -d':' -f2 | xargs)
          spacetotal=$(echo "scale=2; $file_size / $no_kmers * 8" | bc)

          parent_dir=$(dirname "$f")
          file_name=$(basename "$f")

          # low hitrates
          if [[ "$BASENAME" == *"bacterial"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR5833294.fastq.gz")
          elif [[ "$BASENAME" == *"cod"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR11449743_1.fastq.gz")
          elif [[ "$BASENAME" == *"kestrel"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR12858649.fastq.gz")
          elif [[ "$BASENAME" == *"human"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR5901135_1.fastq.gz")
          else
            queries=()
          fi

          for query in "${queries[@]}"; do
              echo $f >> $LOG
              echo $query >> $LOG
              echo $query
              # /usr/bin/time -v -o time.txt $PROGRAM lookup -d "${BASENAME}_.dict" -q $query > OUT 2>&1
              /usr/bin/time -v -o time.txt $PROGRAM locate -d "${BASENAME}_.dict" -q $query > OUT 2>&1

              cat OUT >> $LOG
              
              querytime=$(grep "User time" time.txt | awk -F': ' '{print $2}')
              querymem=$(grep "Maximum resident set size" time.txt | awk -F': ' '{print $2}')
              k_mers=$(grep "num_kmers" OUT | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
              found=$(grep "num_positive_kmers" OUT | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
              positions=$(grep "num_positions" OUT | sed -E 's/.*num_positions = ([0-9]+).*/\1/')
              # querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
              querytimekmer=$(grep 'time_per_kmer' OUT | awk -F'=' '{print $2}' | awk '{print $1}' | sed 's/ns//')
              extensions=$(grep 'extensions' OUT | awk -F'=' '{print $2}' | awk '{print $1}')

              /usr/bin/time -v -o time.txt $PROGRAM bench -d "${BASENAME}_.dict" > OUT 2>&1
              cat OUT >> $LOG
              lookup_pos=$(grep '^pos_time_per_kmer' OUT | cut -d= -f2 | xargs)
              lookup_neg=$(grep '^neg_time_per_kmer' OUT | cut -d= -f2 | xargs)

              echo "lookup_neg: $lookup_neg"
              echo "lookup_pos: $lookup_pos"

              echo "$f,$query,$k,$m1,$m2,$m3,$t1,$t2,$t3,$buildtime,$buildmem",$file_size,$spaceoffsets1,$spaceoffsets2,$spaceoffsets3,$spacer1,$spacer2,$spacer3,$spaces1,$spaces2,$spaces3,$spaceht,$density_r1,$density_r2,$density_r3,$density_s1,$density_s2,$density_s3,$density_ht,$no_minimiser1,$no_minimiser2,$no_minimiser3,$no_distinct_minimiser1,$no_distinct_minimiser2,$no_distinct_minimiser3,$spacetotal,$memtotal,$querytimekmer,$querymem,$extensions,$k_mers",$found,$positions,$lookup_neg" >> "$CSV"

              # rm -f time.txt
              # rm -f OUT
          done

          # high hitrates
          if [[ "$BASENAME" == *"bacterial"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR5901135_1.fastq.gz")
          elif [[ "$BASENAME" == *"cod"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR12858649.fastq.gz")
          elif [[ "$BASENAME" == *"kestrel"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR11449743_1.fastq.gz")
          elif [[ "$BASENAME" == *"human"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR5833294.fastq.gz")
          else
            queries=()
          fi

          for query in "${queries[@]}"; do
              echo $f >> $LOG
              echo $query >> $LOG
              echo $query
              # /usr/bin/time -v -o time.txt $PROGRAM lookup -d "${BASENAME}_.dict" -q $query > OUT 2>&1
              /usr/bin/time -v -o time.txt $PROGRAM locate -d "${BASENAME}_.dict" -q $query > OUT 2>&1

              cat OUT >> $LOG
              
              querytime=$(grep "User time" time.txt | awk -F': ' '{print $2}')
              querymem=$(grep "Maximum resident set size" time.txt | awk -F': ' '{print $2}')
              k_mers=$(grep "num_kmers" OUT | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
              found=$(grep "num_positive_kmers" OUT | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
              positions=$(grep "num_positions" OUT | sed -E 's/.*num_positions = ([0-9]+).*/\1/')
              # querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
              querytimekmer=$(grep 'time_per_kmer' OUT | awk -F'=' '{print $2}' | awk '{print $1}' | sed 's/ns//')
              extensions=$(grep 'extensions' OUT | awk -F'=' '{print $2}' | awk '{print $1}')
              
              echo "$f,$query,$k,$m1,$m2,$m3,$t1,$t2,$t3,$buildtime,$buildmem",$file_size,$spaceoffsets1,$spaceoffsets2,$spaceoffsets3,$spacer1,$spacer2,$spacer3,$spaces1,$spaces2,$spaces3,$spaceht,$density_r1,$density_r2,$density_r3,$density_s1,$density_s2,$density_s3,$density_ht,$no_minimiser1,$no_minimiser2,$no_minimiser3,$no_distinct_minimiser1,$no_distinct_minimiser2,$no_distinct_minimiser3,$spacetotal,$memtotal,$querytimekmer,$querymem,$extensions,$k_mers",$found,$positions,$lookup_pos" >> "$CSV"

              # rm -f time.txt
              # rm -f OUT
          done

        
        done

}

# FILES=( "../DNA_datasets/Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs.fa.ust.fa.gz" "../DNA_datasets/Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs.fa.ust.fa.gz" "../DNA_datasets/Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs.fa.ust.fa.gz" "../DNA_datasets/bacterial.genome.fixed.fa.unitigs.fa.ust.fa.gz" )
# FILES=( "../../datasets/cod.k31.unitigs.fa.ust.fa.gz" "../../datasets/kestrel.k31.unitigs.fa.ust.fa.gz" "../../datasets/human.k31.unitigs.fa.ust.fa.gz")
# FILES=( "../../datasets/human.k31.unitigs.fa.ust.fa.gz" )
# FILES=( "../../datasets/cod.genome.fixed.fa.gz" "../../datasets/kestrel.genome.fixed.fa.gz" "../../datasets/human.genome.fixed.fa.gz" "../../datasets/bacterial.genome.fixed.fa.gz")
# FILES=( "../../datasets/cod.genome.fixed.fa" "../../datasets/kestrel.genome.fixed.fa.gz" "../../datasets/human.genome.fixed.fa.gz" )
FILES=( "../../datasets/cod.genome.preproc.fa" "../../datasets/kestrel.genome.preproc.fa.gz" "../../datasets/human.genome.preproc.fa.gz" )
echo "textfile,queryfile,k,m1,m2,m3,t1,t2,t3,buildtime [s],buildmem [B],indexsize [B],spaceoffsets1 [bits/kmer],spaceoffsets2 [bits/kmer],spaceoffsets3 [bits/kmer],spaceR1 [bits/kmer],spaceR2 [bits/kmer],spaceR3 [bits/kmer],spaceS1 [bits/kmer],spaceS2 [bits/kmer],spaceS3 [bits/kmer],HT [bits/kmer],density_r1 [%],density_r2 [%],density_r3 [%],density_s1 [%],density_s2 [%],density_s3 [%],kmers HT [%],no minimizer1,no minimizer2,no minimizer3, no distinct minimizer1,no distinct minimizer2,no distinct minimizer3,space theo [bits/kmer],mem theo [bits/kmer],querytime [ns/kmer],querymem [B],extensions,kmers,found,positions,lookup" > "$CSV"
for f in "${FILES[@]}"; do
  run $f
done
