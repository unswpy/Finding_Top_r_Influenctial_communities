echo "sh avg_greedy_test.sh dataname size"
./Weight_Community -gname=email -mode=top -alg=casual -func=size_avg -k=4 -r=5 -s=$2
./Weight_Community -gname=email -mode=top -alg=climb -func=size_avg -k=4 -r=5 -s=$2
