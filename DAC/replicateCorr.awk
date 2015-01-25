BEGIN {
}

{
  #a[ARGIND][$6 " " $1]=$1;
  if (FNR==1) next;  # to skip RSEM header line for each file
  b[ARGIND,$1]=log($6+1e-3); # gene TP log-scaled. RSEM min non-zero TPM=0.01
  G[$1]=0;
}

END {
  for (gg in G) {
    for (ii=1;ii<=2;ii++) {
       s[ii]+=b[ii,gg];
       s2[ii]+=b[ii,gg]^2;
    };
    xy+=b[1,gg]*b[2,gg];
  };

  n=length(G);
  P=(n*xy-s[1]*s[2])/sqrt((n*s2[1]-s[1]^2)*(n*s2[2]-s[2]^2));
  
  print P; #Pearson

}
