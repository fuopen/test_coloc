#################################
# this script is used for parsing the configuration json file
# require library rjson

library(rjson)
parseJson<-function(jsonfile){
	jp<-fromJSON(file=jsonfile)
	jp.pipeline=jp$PipelineSetting
	jp.pstructure=jp$PipelineStructure
	jp.input=jp$Inputsetting
#########################################
# check if the path in configuration file exists
	for(check.dir in c('InputPath','OutPath')){
		if(!dir.exists(jp.input[[check.dir]])){
			stop(paste0("path of ",check.dir,":",jp.pipeline[[check.dir]]," doesn't exist"))
		}
	}
	for(check.file in c('InputDataFile','InputLDMatFile')){
		if(!file.exists(paste0(jp.input$InputPath,'/',jp.input[[check.file]]))){
			stop(paste0("path of ",check.file,':',jp.pipeline[[check.file]]," doesn't exist"))
		}
	}
	#if(!file.exists(paste0(jp.input$OutPath,'/Coloc_pipeline.config'))){
	#	write.table(jp.pipeline.config,paste0(jp.pipeline$OutDir,'/Coloc_pipeline.config'),row.names=F,col.names=F,quote=F)
	#}
###########################################
####generate exposure and outcome input, update the header for downstream harmonisation
}

args<-commandArgs(T)
coloc.json.file<-args[1]
#coloc.json.file<-"coloc_configure.json"
print(coloc.json.file)
parseJson(coloc.json.file)
q(save='no')
