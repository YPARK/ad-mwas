#!/usr/bin/awk -f

BEGIN{

    if(length(ROWS) > 0) {
	split(ROWS, rows, ",")
    } else if(length(ROWSF) > 0) {
	getline < ROWSF
	split($0, rows, ",")
    }

    mtot = length(rows)

    for(j=1; j<=mtot; ++j) {
	row2pos[j] = rows[j]
	pos2row[rows[j]] = j
    }

    FS = "\t"
    CUTOFF = 0.01  # p-value cutoff
}
(NR in pos2row){
    row = pos2row[NR]

    { # first one
	j = 38
	v = ($(j+4) < CUTOFF && length($j) > 0)? $j : "NA"
	row2data[row] = v
    }

    for(j=47; j<=NF; j+=9 ){
	v = ($(j+4) < CUTOFF && length($j) > 0)? $j : "NA"
	row2data[row] = row2data[row] FS v
    }
}
END{
    for(j = 1; j<= mtot; ++j)
	print row2data[j]
}
