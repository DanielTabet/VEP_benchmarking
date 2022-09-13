loadConfig = function(configFilePath) {
  cat("Loading configuration file\n")
  library(yaml)
  config = read_yaml(configFilePath)
  
  return(config)
}