library('reticulate')
packages_path<-function(){
  packages_dir<-system.file(package = "STAE")
  #packages_dir <- c('E:/work/sunhang/code/package_0708')
  return(packages_dir)
}
#'@title PythonEnvSet
#'@description
#'Set the path of the conda environment you use
#'Many python programs are called in our program,
#'so you need to set the path of the python interpreter,
#' and you can set it directly to the conda environment here.
#'@param py_path path of the conda environment
#'@import reticulate
#'@export
pythonEnvSet <- function(py_path){
  # library(reticulate)
  use_condaenv(py_path)
  #use_python(py_path)
}
#'@title TestEnv
#'Detect environment dependencies of python
#'@description
#'You can use this function to detect if a package is missing from a dependent python environment.
#'@return bool TRUE or FALSE
#'@import reticulate
#'@export
testEnv<-function(){
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
                    'scipy',
                    'rich',
                    'kaleido'
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

#'@title SetParameters
#'Hyperparameter setting
#'@description
#'Set Resolution and Cell Columns,resolution is the multiple of program amplification,
#'Cell Columns is the specified column name
#'@param pdr position_distance_ratio
#'@param dam Differentiation and migration of cell type
#' @param pseflag Whether there is pseudo-temporal data(True/False) : True
#'@export
setParameters <- function(pdr,pseflag){
  parameter_settings_csv<-c(pdr,pseflag)
  dir=packages_path()
  parameter_settings_path = paste(dir,"/data/parameter_settings.csv",sep="")
  write.table(parameter_settings_csv,file = parameter_settings_path,row.names = FALSE,col.names = 'parameter')
}

#'@import reticulate
callPythonProgram<-function(pyname){
  #package_path<-c('E:/work/sunhang/code/package_0708')
  package_path_dir <- packages_path()
  os<-import('os')
  os$chdir(package_path_dir)
  print(pyname)
  py_dir=paste(package_path_dir,'/code/',pyname,'.py',sep="")
  source_python(py_dir)
}

#'@title dataProcessing
#'@param bimr before_iterative_mapping_result
#'@param aimr after_iterative_mapping_result
#'@param bsd before_sc_data
#'@param asd after_sc_data
#'@param pse pseudotime
#'@import data.table
#'@export
dataProcessing<-function(bimr,aimr,bsd,asd,pse = 1){
  library(data.table)
  # before_iterative_mapping_result <- bimr
  # after_iterative_mapping_result <- aimr
  # before_sc_data <- bsd
  # after_sc_data <- asd
  # marker_gene <- mg
  # pseudotime <- pse
  data_path <- paste(packages_path(),'/data',sep = '')
  if(!dir.exists(data_path)){
    dir.create(data_path)
  }
  fwrite(as.data.frame(bimr),file = paste(data_path,'/before_iterative_mapping_result.csv',sep = ''))
  fwrite(as.data.frame(aimr),file = paste(data_path,'/after_iterative_mapping_result.csv',sep = ''))
  fwrite(as.data.frame(bsd),file = paste(data_path,'/before_sc_data.csv',sep = ''))
  fwrite(as.data.frame(asd),file = paste(data_path,'/after_sc_data.csv',sep = ''))
  fwrite(as.data.frame(pse),file = paste(data_path,'/pseudotime.csv',sep = ''))
}
#'@title stae
#'@description main programe
#'@export
stae <- function(pdr,pseflag = FALSE){
  setParameters(pdr = pdr,pseflag = pseflag)
  callPythonProgram('move_center')
  callPythonProgram('TL_pic_distance_new')
  callPythonProgram('TL_sample_get_adata')
  callPythonProgram('AE')
  callPythonProgram('TL_get')
  reault_path <- paste(packages_path(),'/data/result/tow_time_TL_edges.csv',sep = "")
  stae_result <- fread(input = reault_path)
  return(stae_result)
}
#'@title staePlot
#' @description draw result
#' @export
staePlot <- function(dam){
  Parameter_settings(pdr = dam,pseflag = FALSE)
  callPythonProgram('pic')
}


