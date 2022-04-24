#!/usr/bin/env bash

for i in 1
do
  for mode in "top"
  do
    for gname in "fs"
#    for gname in "lj"
    #for gname in "email" "dblp" "youtube" "orkut" "lj" "fs"
#    for gname in "email" "dblp"
    do
        for func in "size_sum" "size_avg"
#        for func in "size_sum"
#        for func in "sum"
        do
          if [ $func = "size_avg" ];
          then
            for alg in "casual" "climb"
#            for alg in "climb" "casual"
            do
              for k in 4 6 8 10
              do
                for r in 5 10 15 20
                do
                  for s in 5 10 15 20
                  do
                    ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r -s=$s 
                  done
                done
              done
            done
          elif [ $func = "size_sum" ];
          then
            for alg in "casual" "climb"
            do
              for k in 4 6 8 10
              do
                for r in 5 10 15 20
                do
                  for s in 5 10 15 20
                  do
                    ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r -s=$s
                  done
                done
              done
            done
          fi
        done
    done
  done
done
