#!/usr/bin/env python
# -*- coding:utf-8 -*-

###################################
# concatenate consecutive domains #
# input: bed file                 #
# output: bed file                #
# e.g                             #
# input:                          #
#     chr3  6001  9000  red       #
#     chr3  9001  12000 red       #
#     chr3  12001 15000 red       #
# output:                         #
#     chr3  6001 15000 red        #
###################################

# test run
# ./domainIt.py test_data/test_domainIt.bed


import sys

with open(sys.argv[1]) as f:
  preLine = f.readline()
  preItem = preLine.rstrip().split('\t')
  preChr = preItem[0]
  preStart = preItem[1]
  preTag = preItem[3]
  for line in f:
    item = line.rstrip().split('\t')
    if item[3] != preTag or item[0] != preChr:
      sys.stdout.write(preLine)
      preLine = line
      preChr = item[0]
      preStart = item[1]
      preEnd = item[2]
      preTag = item[3]
    elif item[3] == preTag and item[0] == preChr:
      preLine = "%s\t%d\t%d\t%s\n" %(
          item[0], int(preStart), 
          (int(item[2])- (int(preStart)-1)) + (int(preStart)-1), 
          item[3] )
sys.stdout.write(preLine)

