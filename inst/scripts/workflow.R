

# to use: set DBSNP_VCF_PATH to location of gzipped VCF and .gz.tbi
# set DBSNP_PARQ_PATH to folder where you want converted parquet to live


## PHASE 1 -- PARSE SELECTED VCF to PARQUET

source("dbsnp_build.R", echo=TRUE)
path_to_vcf = Sys.getenv("DBSNP_VCF_PATH")
path_to_parq = Sys.getenv("DBSNP_PARQ_PATH")

stopifnot(file.exists(path_to_vcf))
stopifnot(dir.exists(path_to_parq))

dbsnp_chrom_to_parquet(path_to_vcf, chrom="17", outdir=path_to_parq, release="b157", build="GRCh38")

# assuming
#-rw-rw-r-- 1 exouser exouser     3120044 Nov 20 18:45 GCF_000001405.40.gz.tbi
#-rw-rw-r-- 1 exouser exouser 29552227779 Nov 20 18:21 GCF_000001405.40.gz
# are available at DBSNP_VCF_PATH

# this call
# dbsnp_chrom_to_parquet(path_to_vcf, chrom="17", outdir=path_to_parq, release="b157", build="GRCh38")
# leads to:
#
#wrote 35376823 rows per table to newparq/b157/GRCh38
#>
#[4]+  Stopped                 R
#exouser@anno324:~/DBSNP/dbSNP4Bioc/inst/scripts$ tree newparq
#newparq  # 415 Mb
#└── b157
#    └── GRCh38
#        ├── by_pos
#        │   └── 17.parquet
#        ├── by_rsid
#        │   └── 17.parquet
#        └── manifest.json
#
#5 directories, 3 files

## PHASE 2 -- cache the parquet for reuse

source("dbsnp_client.R")          # copy dbsnp_client.R to that machine too
library(BiocFileCache)

bfc   <- BiocFileCache(ask = FALSE)          # default cache; pass a path to use another
rel   <- "b157"
build <- "GRCh38"

base <- file.path(Sys.getenv("DBSNP_PARQ_PATH"), rel, build)

for (f in c("manifest.json", "by_pos/17.parquet", "by_rsid/17.parquet"))
  bfcadd(bfc, rname = sprintf("dbsnp/%s/%s/%s", rel, build, f),
         fpath = file.path(base, f), rtype = "local", action = "copy")

dbsnp_cached(rel, build, bfc)                 # should list chrom 17, by_pos and by_rsid
# db <- dbsnp_open(rel, build, source = "cache", cache = bfc)   # offline from now on

## PHASE 3 -- query the cached resource

library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)

eg  <- mapIds(org.Hs.eg.db, "ORMDL3", "ENTREZID", "SYMBOL")
gr  <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, filter = list(gene_id = eg))
win <- resize(gr, width(gr) + 2 * 50000, fix = "center")   # gene body +/- 50 kb
as.character(seqnames(win)); start(win); end(win)         # sanity check: chr17

db  <- dbsnp_open("b157", "GRCh38", chroms = "17", source = "cache")  # or "remote"
snp <- snps_in_range(db, "17", start(win), end(win))
nrow(snp); head(snp)

snp2 <- snps_in_range(db, "17", start(win), end(win), biallelic_only = TRUE)   # n_alt == 1 only
dbsnp_close(db)
