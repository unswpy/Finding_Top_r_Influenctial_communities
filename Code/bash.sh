#!/usr/bin/env bash

for i in 1
do
  for mode in "top"
  do
    #for gname in "fs"
#    for gname in "lj"
#    for gname in "email" "dblp" "youtube" "orkut" "lj"
    for gname in "email" "dblp"
    do
      if [ $gname = "email" ];
      then
        for func in "size_sum" "size_avg"
#        for func in "size_sum"
#        for func in "sum"
        do
          if [ $func = "sum" ];
          then
            for alg in "naive" "improve" "approx"
#            for alg in "improve"
            do
              for k in 4 6 8 10
              do
                for r in 5 10 15 20
                do
                  if [ $alg = "approx" ];
                  then
                    for eps in 0.01 0.05 0.1 0.2 0.3 0.5 0.8
                    do
                      ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r -eps=$eps
                    done
                  else
                    ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r
                  fi
                done
              done
            done
          elif [ $func = "size_avg" ];
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
      elif [ $gname = "dblp" ];
      then
        for func in "size_sum" "size_avg"
        do
          if [ $func = "sum" ];
          then
            for alg in "naive" "improve" "approx"
#            for alg in "improve"
            do
              for k in 4 6 8 10
              do
                for r in 5 10 15 20
                do
                  if [ $alg = "approx" ];
                  then
                    for eps in 0.01 0.05 0.1 0.2 0.3 0.5 0.8
                    do
                      ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r -eps=$eps
                    done
                  else
                    ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r
                  fi
                done
              done
            done
          elif [ $func = "size_avg" ];
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
      elif [ $gname = "youtube" ];
      then
        for func in "size_sum" "size_avg"
        do
          if [ $func = "sum" ];
          then
            for alg in "naive" "improve" "approx"
#            for alg in "improve"
            do
              for k in 4 6 8 10
              do
                for r in 5 10 15 20
                do
                  if [ $alg = "approx" ];
                  then
                    for eps in 0.01 0.05 0.1 0.2 0.3 0.5 0.8
                    do
                      ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r -eps=$eps
                    done
                  else
                    ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r
                  fi
                done
              done
            done
          elif [ $func = "size_avg" ];
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
      elif [ $gname = "orkut" ];
      then
        for func in "size_sum" "size_avg"
        do
          if [ $func = "sum" ];
          then
            for alg in "naive" "improve" "approx"
#            for alg in "improve"
            do
              for k in 40 50 100 150 200
              do
                for r in 5 10 15 20
                do
                  if [ $alg = "approx" ];
                  then
                    for eps in 0.01 0.05 0.1 0.2 0.3 0.5 0.8
                    do
                      ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r -eps=$eps
                    done
                  else
                    ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r
                  fi
                done
              done
            done
          elif [ $func = "size_avg" ];
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
      elif [ $gname = "lj" ];
      then
        for func in "size_sum" "size_avg"
        do
          if [ $func = "sum" ];
          then
            for alg in "naive" "improve" "approx"
#            for alg in "improve"
            do
              for k in 40 50 100 150 200
              do
                for r in 5 10 15 20
                do
                  if [ $alg = "approx" ];
                  then
                    for eps in 0.01 0.05 0.1 0.2 0.3 0.5 0.8
                    do
                      ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r -eps=$eps
                    done
                  else
                    ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r
                  fi
                done
              done
            done
          elif [ $func = "size_avg" ];
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
      elif [ $gname = "fs" ];
      then
        for func in "size_avg" "size_sum"
        do
          if [ $func = "sum" ];
          then
#            for alg in "improve" "approx" "naive"
            for alg in "naive"
            do
              for k in 200 150 100 50 40
#              for k in 100 125 150 175 200
              do
                if [ $alg = "approx" ];
                then
                  for r in 5 10 15 20
                  do
                    for eps in 0.01 0.05 0.1 0.2 0.3 0.5 0.8
                    do
                      timeout 86400 ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r -eps=$eps
                    done
                  done
                elif [ $alg = "naive" ]
                then
                  for r in 5
                  do
                    ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r
                  done
                else
                  for r in 5 10 15 20
                  do
                    timeout 86400 ./Weight_Community -gname=$gname -mode=$mode -alg=$alg -func=$func -k=$k -r=$r
                  done
                fi
              done
            done
          elif [ $func = "size_avg" ];
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
      fi
    done
  done
done
