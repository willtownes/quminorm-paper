#useful functions for data loading
library(data.table)
fp<-file.path

parse_sra_csv<-function(bp,kallisto_tech,species_map,kallisto_threads=8,gz=TRUE){
  #bp is the base path for the dataset such as real/macosko_2015
  #kallisto_tech is the droplet protocol (run kallisto bus -l for available options)
  #species_map: a list where the names are unique values of the SRA run info ScientificName and the values are the file names of kallisto transcriptome indices.
  #set gz=FALSE if planning to extract FASTQ from SRA using fasterq-dump instead of fastq-dump
  #saves files under extdata that are useful for automated data processing
  srp <- read.csv(fp(bp,"extdata/SraRunInfo.csv"),stringsAsFactors=FALSE)
  coldata <- srp[,c("SampleName","Run","avgLength","Experiment","Sample","BioSample","ScientificName","download_path")]
  colnames(coldata)[1] <- "geo_accession"
  rownames(coldata) <- coldata$Run
  srr<-coldata$Run
  write(srr,file=fp(bp,"extdata/sraFiles.txt"))
  write(coldata$download_path,file=fp(bp,"extdata/sraFilesPath.txt"))
  fwrite(coldata,file=fp(bp,"extdata/coldata.txt"),sep=" ")
  
  #generate file with list of kallisto commands
  cmd_base<-paste("kallisto bus -t",kallisto_threads,"-x",kallisto_tech)
  idx<-fp("../../resources",unlist(species_map[coldata$ScientificName]))
  kof<-fp("data/original/kallisto",srr)
  ending<-".fastq"
  if(gz){ ending<-paste0(ending,".gz") }
  fq1<-fp("data/original/fastq",paste0(srr,"_1",ending))
  fq2<-fp("data/original/fastq",paste0(srr,"_2",ending))
  cmd<-paste(cmd_base,"-i",idx,"-o",kof,fq1,fq2)
  write(cmd,file=fp(bp,"extdata/kallisto_cmd.txt"))
}

ec2counts_inner<-function(g,fc){
  #function to apply to a single barcode x UMI combo
  #the length of g,fc is the number of equivalence classes
  #g is list of genes (a list where each element is a char vector)
  #fc is vector of fragment counts
  gene<-Reduce(intersect,g)
  #gene is intersection of genes across all ec's
  n<-length(gene)
  if(n>0){
    umi_ct<-1/n
    read_ct<-sum(fc)*umi_ct
  } else {
    umi_ct<-read_ct<-numeric()
  }
  mget(c("gene","umi_ct","read_ct"))
}

ec2counts<-function(pth,sp=c("Mus musculus","Homo sapiens"),ensembl_version=95,bc_wl=NULL,verbose=FALSE){
  #bc_wl a character vector of valid barcodes, all others will be excluded
  sp<-gsub(" ","_",match.arg(sp),fixed=TRUE)
  gxname<-paste0("tr2g_ensembl",ensembl_version,"_",sp,".txt")
  tr2g<-fread(fp("./resources",gxname))
  ecg<-as.data.table(BUSpaRse::EC2gene(tr2g,pth,verbose=verbose))
  colnames(ecg)<-c("ec","transcripts","genes")
  if(verbose){ print("reading bus file") }
  bus<-fread(fp(pth,"output_sort.txt"))
  colnames(bus)<-c("bc","umi","ec","frag_ct")
  if(!is.null(bc_wl)){
    if(verbose){ print("barcode whitelist detected, will subset barcodes") }
    bus<-bus[bc %in% bc_wl,]
  } else {
    if(verbose){ print("barcode whitelist not detected, will use all barcodes") }
  }
  res<-merge(bus,ecg[,.(ec,genes)],all.x=TRUE,all.y=FALSE,sort=FALSE)
  #data.table with columns ec,bc,umi,frag_ct and genes.
  #the 'genes' column is a "list column"
  #collapse ec, expand genes so each row is a bc/umi/gene combination
  if(verbose){ print("collapsing over equivalence classes, expanding over genes") }
  res<-res[,ec2counts_inner(genes,frag_ct),.(bc,umi)]
  #collapse over UMIs
  if(verbose){ print("collapsing over UMIs") }
  res[,.(umi_ct=sum(umi_ct),read_ct=sum(read_ct)),.(bc,gene)]
}

bus2genecount_all<-function(bp="./",ensembl_version=95,bc_filter=FALSE){
  coldata<-fread(fp(bp,"extdata/coldata.txt"))
  #check for barcode whitelist
  bcname<-fp(bp,"extdata/barcodes.txt")
  if(bc_filter && file.exists(bcname)){ 
    bc_wl<-scan(bcname,character()) 
  } else { 
    bc_wl<-NULL 
  }
  odir<-fp(bp,"data/original/genecounts")
  if(!dir.exists(odir)){ dir.create(odir,recursive=TRUE) }
  for(i in 1:nrow(coldata)){
    srr<-coldata$Run[i]
    print(srr)
    pth<-fp(bp,"data/original/kallisto",srr)
    ofile<-fp(odir,paste0(srr,".txt"))
    if(!file.exists(ofile)){
      res<-ec2counts(pth,coldata$ScientificName[i],ensembl_version=ensembl_version,bc_wl=bc_wl)
      fwrite(res,file=ofile,sep=" ")
    }
  }
}
