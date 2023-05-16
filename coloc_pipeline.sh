#!/bin/bash

Help () {

	echo "you need to provide the following arguments to run the coloc pipeline"
	echo "-c, --configure:"
	echo "coloc configuration json file,this option is mandatory!"

	echo "-h, --help"
	echo "help information"

	echo "-v, --verbose"
	echo "verbose setting, enabled by default"
}

if test $# -eq 0;then
	echo "No arguments provided, quit!"
	Help
	exit 1
fi

sft=1
coloc_CONFIG=
vb=
while ! test $# -eq 0;do
	case $1 in
	-h | --help)
		Help
		exit 0
	;;
	-c | --configure)
		coloc_CONFIG=$2
		if ! test -f $coloc_CONFIG;then
			echo "the json configuration file doesn't exist! quit"
			exit 1
		fi
		sft=2
	;;
	-v | --verbose)
		vb=1
		sft=1
	;;
	*)
		echo "unrecongnised arguments,quit"
		Help
		exit 1
	;;
	esac
	shift $sft
done

###############################################
####step1 parse json configuration file
if ! test -z $vb;then
	echo "step1 start, parsing the json configuration!"
fi

if ! test -f $coloc_CONFIG;then
    echo "Coloc config json file doesn't exist!"
    exit 1
else
    echo $coloc_CONFIG
fi
colocCodeDir=$(jq -r '.PipelineSetting.CodePath' $coloc_CONFIG)

if ! test -z $vb;then
	echo "Find the parseJson script at $colocCodeDir"
fi

colocJsonScript="$colocCodeDir/coloc_Json_parse.R"

if ! test -f $colocJsonScript;then
	echo "coloc_Json_parse.R doesn't exist at $colocCodeDir !"
	exit 1
fi

Rscript $colocJsonScript $coloc_CONFIG
colocInputDir=$(jq -r '.Inputsetting.InputPath' $coloc_CONFIG)
colocOutDir=$(jq -r '.Inputsetting.OutPath' $coloc_CONFIG)
Author=$(jq -r '.Inputsetting.user' $coloc_CONFIG)
studyName=$(jq -r '.Inputsetting.study_name' $coloc_CONFIG)
colocplotLDFile=$(jq -r '.Inputsetting.PlotLDPath' $coloc_CONFIG)
colocInput=$(jq -r '.Inputsetting.InputDataFile' $coloc_CONFIG)
colocInputFile=${colocInputDir}/${colocInput}
colocLDMatInput=$(jq -r '.Inputsetting.InputLDMatFile' $coloc_CONFIG)
colocLDMatFile=${colocInputDir}/${colocLDMatInput}
leadSNPID=$(jq -r '.Inputsetting.leadSNP' $coloc_CONFIG)
GeneRegFile=$(jq -r '.Inputsetting.GeneRegFile' $coloc_CONFIG)
GeneticMapFile=$(jq -r '.Inputsetting.GeneticMap' $coloc_CONFIG)
###################################################################
####Run step2: coloc(Suise) analysis

if ! test -z $vb;then
	echo "Run the coloc analysis and generate report!"
fi

if ! test -z $vb;then
	echo "Rscript ${colocCodeDir}/coloc_report.R ${colocOutDir} $Author $studyName ${colocCodeDir}/coloc_report.Rmd $colocInputFile $colocLDMatFile $colocplotLDFile $GeneRegFile $GeneticMapFile $leadSNPID"
fi

Rscript ${colocCodeDir}/coloc_report.R ${colocOutDir} $Author $studyName "${colocCodeDir}/coloc_report.Rmd" $colocInputFile $colocLDMatFile $colocplotLDFile $GeneRegFile $GeneticMapFile $leadSNPID

echo "All coloc analysis finished!"
