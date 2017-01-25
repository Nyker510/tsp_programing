#!/bin/sh

g++ -std=gnu++11 -O3 $1 -o tsp_mine
g++ -std=gnu++11 -O3 original_tsp.cpp -o tsp_original

root_tsp="./new_tsp_instances"
RUNORIGIN=0
output_filepath=$2
output_log=log.txt

rm $output_filepath
echo "p,n,time,optimal_value,branches,leaves,stop_forced_branch,stop_delete_branch,stop_bounding,sop_nonedge_can_select,method" >> $output_filepath

make_csv_string(){
  if [ $# -eq 1 ]; then
    echo $1
  else
    str=$1
    for s in $*; do
      if [ $s != $str ]; then
        str="$str,$s"
      fi
    done
    echo $str
  fi
}

for p in 0.3; do
  for n in 6 7 8 9 10 12 14 16 18 20; do
    fp="$root_tsp/p_$p/n_$n"
    for f_name in $fp/*; do
      echo "$f_name"

      if [ $RUNORIGIN -eq 1 ]; then
        result=`./tsp_original < $f_name`
        make_csv_string $p $n $result "original" >> $output_filepath
        printf "original\t$result\n"
      fi

      result=`./tsp_mine < $f_name`
      make_csv_string $p $n $result "mine" >> $output_filepath
      printf "mine\t$result\n"
    done
  done
done
