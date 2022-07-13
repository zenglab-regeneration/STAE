library('reticulate')
packages_path<-function(){
  packages_dir<-system.file(package = "STAE")
  #packages_dir <- c('E:/work/sunhang/code/package_0708')
  return(packages_dir)
}
#'@title env_python_set
#'@description
#'Set the path of the conda environment you use
#'Many python programs are called in our program,
#'so you need to set the path of the python interpreter,
#' and you can set it directly to the conda environment here.
#'@param py_path path of the conda environment
#'@import reticulate
#'@export
env_python_set<-function(py_path){
  # library(reticulate)
  use_condaenv(py_path)
  #use_python(py_path)
}
#'@title env_test
#'Detect environment dependencies of python
#'@description
#'You can use this function to detect if a package is missing from a dependent python environment.
#'@return bool TRUE or FALSE
#'@import reticulate
#'@export
env_test<-function(){
  library(reticulate)
  package_flag<-TRUE
  package_detect<-c('time',
                    'plotly',
                    'pandas',
                    'numpy',
                    'networkx',
                    'scanpy',
                    'sklearn',
                    'anndata',
                    'torch',
                    'scipy'
                    )
  for(package_name in package_detect){
    if(!py_module_available(package_name)){
      py_install(package_name, pip = T)
    }
  }
  for(package_name in package_detect){
    if(!py_module_available(package_name)){
      error_message<-paste(package_name,"is not ready",sep=" ")
      print(error_message)
    }
  }
  return(package_flag)
}


#'@import reticulate
call_python_program<-function(pyname){
  #package_path<-c('E:/work/sunhang/code/package_0708')
  package_path_dir <- packages_path()
  os<-import('os')
  os$chdir(package_path_dir)
  print(pyname)
  py_dir=paste(package_path_dir,'/code/',pyname,'.py',sep="")
  source_python(py_dir)
}

#'@title data_deal
#'@param bimr before_iterative_mapping_result
#'@param aimr after_iterative_mapping_result
#'@param bsd before_sc_data
#'@param asd after_sc_data
#'@param mg marker_gene
#'@param pse pseudotime
#'@export
data_deal<-function(bimr,aimr,bsd,asd,mg,pse){
  before_iterative_mapping_result <- bimr
  after_iterative_mapping_result <- aimr
  before_sc_data <- bsd
  after_sc_data <- asd
  marker_gene <- mg
  pseudotime <- pse
  data_path <- paste(packages_path(),'/data',sep = '')
  if(!dir.exists(data_path)){
    dir.create(data_path)
  }
  write.csv(before_iterative_mapping_result,file = paste(data_path,'/before_iterative_mapping_result.csv',sep = ''))
  write.csv(after_iterative_mapping_result,file = paste(data_path,'/after_iterative_mapping_result.csv',sep = ''))
  write.csv(before_sc_data,file = paste(data_path,'/before_sc_data.csv',sep = ''))
  write.csv(after_sc_data,file = paste(data_path,'/after_sc_data.csv',sep = ''))
  write.csv(marker_gene,file = paste(data_path,'/marker_gene.csv',sep = ''))
  write.csv(pseudotime,file = paste(data_path,'/pseudotime.csv',sep = ''))
}
#'@title stae_main 
#'@description main programe
#'@export
stae_main <- function(){
  call_python_program('move_center')
  call_python_program('TL_pic_distance_new')
  call_python_program('TL_gene_similar_comp')
  call_python_program('AE')
  call_python_program('hvg_gene_distance')
  call_python_program('prepare_comp')
  call_python_program('TL_get')
}


