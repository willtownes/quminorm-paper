#Miscellaneous functions used by other scripts

##### Data Loading functions #####

get_10x_readcounts<-function(counts_folder,mol_info_h5){
  #provide two file paths:
  #counts_folder, containing "barcodes.tsv","genes.tsv","matrix.mtx"
  #mol_info_h5, an hdf5 file with the 10x molecule information file
  #output: a SingleCellExperiment with two sparse Matrix assays:
  #"counts" (UMI counts), and "read_counts"
  sce0<-DropletUtils::read10xCounts(counts_folder)
  colnames(sce0)<-colData(sce0)$Barcode
  m<-assay(sce0,"counts")
  gg<-Matrix::rowSums(m)>0 #remove genes that are all zero
  sce<-sce0[gg,]
  m<-m[gg,]
  cm<-colData(sce0)
  bc<-substr(cm$Barcode,1,nchar(cm$Barcode)-2) #assumes barcode ends in "-1"
  cm$bc_enc<-DropletUtils::encodeSequences(bc)
  mi0<-as.data.frame(rhdf5::h5dump(mol_info_h5))
  if(!all(cm$bc_enc %in% mi0$barcode)){
    stop("Some count matrix barcodes not found in molecule information file!")
  }
  gidx<-read.table(file.path(counts_folder,"genes.tsv"),stringsAsFactors=FALSE)[,1]
  mi<-subset(mi0, barcode %in% cm$bc_enc & gene<length(gidx))
  mi$gene_symbol<-gidx[mi$gene+1]
  cm2<-cm[,c("bc_enc","Barcode")]
  colnames(cm2)<-c("barcode","barcode_str")
  mi<-merge(mi,cm2,by="barcode")
  rc<-DropletUtils::makeCountMatrix(mi$gene_symbol,mi$barcode_str,value=mi$reads)
  rc<-rc[rownames(m),colnames(m)]
  #umi<-DropletUtils::makeCountMatrix(mi$gene_symbol,mi$barcode_str)
  #umi<-umi[rownames(m),colnames(m)]
  #umi counts from molecule info file consistent with sce
  #all(m==umi)
  if(!all((rc>0) == (m>0) )){
    stop("zero pattern inconsistent between counts matrix and molecule info")
  }
  if(!all(rc>=m)){
    stop("Read counts not all >= UMI counts")
  }
  assay(sce,"read_counts")<-rc
  sce
}

##### ID conversion functions #####

ensembl2symbol<-function(genes,sp=c("hsapiens","mmusculus")){
  #genes is a character vector of ensembl IDs
  #sp=the species
  #returns a data frame with the ensembl ID mapped to symbols
  #if the gene didn't have a corresponding symbol, it is removed from the data frame.
  sp<-match.arg(sp)
  symb<-list(hsapiens="hgnc",mmusculus="mgi")
  bm<-biomaRt::useMart("ensembl")
  mart<-biomaRt::useDataset(paste0(sp,"_gene_ensembl"),bm)
  bml<-biomaRt::getBM(filters="ensembl_gene_id",attributes=c("ensembl_gene_id",paste0(symb[[sp]],"_symbol")),values=genes,mart=mart)
  #make sure order matches the input genes order
  #if no matching symbol was found, use NA
  genes<-as.data.frame(genes,stringsAsFactors=FALSE)
  colnames(genes)<-"ensembl_gene_id"
  plyr::join(genes,bml,by="ensembl_gene_id",type="left",match="first")
}

##### Downsampling functions #####

Down_Sample_Matrix<-function(expr_mat,min_lib_size=NULL){
  #adapted from https://hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html#normalisations
  min_sz<-min(colSums(expr_mat))
  if(is.null(min_lib_size)){
    min_lib_size<-min_sz
  } else {
    stopifnot(min_lib_size<=min_sz)
  }
  down_sample<-function(x){
    prob <- min_lib_size/sum(x)
    unlist(lapply(x,function(y){rbinom(1, y, prob)}))
  }
  apply(expr_mat, 2, down_sample)
}
