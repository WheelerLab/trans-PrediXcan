#!/usr/bin/python3
import gzip
#gene summary file
mydir = "/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/NCBI_Gene_Summaries/"

summfile = mydir + "hgFixed.refSeqSummary.gz" #available at http://hgdownload.soe.ucsc.edu/goldenPath/hgFixed/database/refSeqSummary.txt.gz
gen2refseqfile = mydir + "gencode.v19.metadata.RefSeq.gz"
gencodefile = "/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/hg19_GENCODE19/gencode.v19.annotation.gtf.gz"
outfile = mydir + "hgFixed.refSeqSummary_gencode.v19.annotations.txt.gz"

#dict of keys=refseq IDs, values=ENST IDs
refdict = {}
for line in gzip.open(gen2refseqfile,'rb'):
    arr = line.strip().split()
    refseq = arr[1].split(b'.')[0] #only keep id up to period to match summfile
    refdict[refseq] = arr[0]

#dict of transcripts from gencode annotation file
#keys=ENST ID, values=ENSG id <TAB> all info
gendict = {}
for line in gzip.open(gencodefile,'rb'):
    if(line.startswith(b'##')): #need b b/c opening .gz file in bytes mode
            continue
    arr = line.strip().split(b'\t')
    entry = arr[2]
    if entry == b"transcript":
        col9 = arr[8]
        col9list = col9.split(b";")
        geneidlist = col9list[0].split(b'"')
        geneid = geneidlist[1]
        txidlist = col9list[1].split(b'"')
        txid = txidlist[1]
        val = geneid + b"\t" + line
        gendict[txid] = val

       
#output file for joined annotations
out = gzip.open(outfile,'wb')
#join on refseq ID first, then on ENST ID
for line in gzip.open(summfile,'rb'):
    arr = line.strip().split(b'\t')
    #only keep RefSeq IDs with a gene summary
    if(len(arr) == 3):
        refseq = arr[0]
        summary = arr[2]
        if refseq in refdict:
            enst = refdict[refseq]
            if enst in gendict:
                geninfo = gendict[enst]
                ensg = geninfo.split(b'\t')[0]
                out.write(summary + b'\t' + refseq + b'\t' + geninfo) #newline already there

out.close()
           
