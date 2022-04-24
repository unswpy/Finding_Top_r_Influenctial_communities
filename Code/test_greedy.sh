echo "sh test_greedy.sh dataname"
./Weight_Community -gname=$1 -mode=top -alg=casual -func=size_sum -k=4 -r=5 -s=20
./Weight_Community -gname=$1 -mode=top -alg=climb -func=size_sum -k=4 -r=5 -s=20
