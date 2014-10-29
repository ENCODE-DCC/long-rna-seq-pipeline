#!/bin/bash
# AWK script provided by Alex Dobin to compare 2 bedGraph files 

awk 'BEGIN {print "Lines present in 2nd file but not in the 1st:"} \
     {if (ARGIND==1) { \
         s[$1 "_" $2 "_" $3]=$4 \
     } else { \
         a=$1 "_" $2 "_" $3; \
         if (a in s) { \
             d1=s[a]-$4; d1=d1>0 ? d1 : -d1; d=d>d1?d:d1; s[a]=-1 \
         } else { \
             print \
         } \
     }} END { \
         print "Maximum difference in the values=" d "\nLines present in 1st file but not in the 2nd:"; \
         for (ii in s) { \
             if (s[ii]!=-1) print ii,s[ii] \
         }; print "Completed!"}' $1 $2

