#!/usr/bin/python
# -*- coding: UTF-8 -*-

import sys
print "python pre_process.py *.txt"
i = 0;
out = open("result.txt", "w");
for name in sys.argv:
    if i == 0:
        i+=1;
        continue;
    i+=1;
    f = open(name);
    lines = f.readlines();
    out.write(name + "\t" + lines[1].strip('\n') + "\t" + lines[3].strip('\n') + "\n")
    f.close();
out.close();

