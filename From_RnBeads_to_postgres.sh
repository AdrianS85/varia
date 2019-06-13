#Appearently postgres doesnt work well with dots in header? Also "end" seem to fuck it up. So, lets change all the dots in headers into underscores using this:

# Prepare blueprint file for parallel, e.g. like this, doesnt matter
ls *.csv > p1; sed s/.csv/.ass/g p1 >p2; paste p1 p2 >pairs;

parallel --colsep '\t' "sed -e '1 s/[\.]/_/g' -e '1 s/End/End_/g' {1} >{2}" :::: pairs
parallel --colsep  '\t'  "mv {2} {1}" :::: pairs
