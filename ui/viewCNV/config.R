# These are OS and file system-dependent settings that must be manually configured

shell <- "C:\\cygwin64\\bin\\bash.exe -lc"

extent.offset <- 12*100/2 # FIXME: outputs should set map positions centered in their windows
extent.offfset <- 0

profile.fn <- "d:/mccarroll/gpc_wave2_batch1/profile_seq_20_100.dat.gz"
profile.fn <- '/home/dkulp/data/gpc_wave2_batch1/profile_seq_20_100.dat.gz'
#profile.fn <- '/home/dkulp/mccarroll/gpc_wave2_batch1/unmasked/profile_seq_20_100.dat.gz'

data.dir <- "d:/mccarroll/cnv_seg.12.500"
data.dir <- "C:/cygwin64/home/dkulp/data/out/cnv_seg.B12.L500.Q13.4"
#data.dir <- '/home/dkulp/data/out/cnv_seg.B12.L500.Q13.4'

tmp.dir <- "C:\\cygwin64\\tmp"

# source(input.dump.fn); input <- input.dump
input.dump.fn <- 'C:\\cygwin64\\home\\dkulp\\data\\tmp\\inputdump.R'

db.conn.str <- 'dkulp:localhost:5432:seq'

