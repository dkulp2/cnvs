# NBINS = e.g. 10. Creates an input file specifying coordinates of NBINS width.
# MAXLEN = e.g. 2000. Skip segments that are too long, indicating a gap not to span.
{
    if (NR > 1) {
        bin = NR-1;
        binIndex = bin - (NBINS-1)*int(bin/(NBINS-1));
        if (binIndex == 0) { binIndex = NBINS-1 };
        if (bin >= NBINS && ($4-STARTS[binIndex]) < MAXLEN ) {
            id = "SEG_" $2 "_" STARTBIN[binIndex] "_" STARTS[binIndex] "_" $4;
            print id, $2, STARTS[binIndex], $4;
        }
        STARTS[binIndex] = $3
	STARTBIN[binIndex] = substr($1,2)
    }
}
