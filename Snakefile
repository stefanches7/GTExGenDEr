import pathlib
import wbuild
config['wBuildPath'] =  str(pathlib.Path(wbuild.__file__).parent)

configfile: "wbuild.yaml"
include: config['wBuildPath'] + "/wBuild.snakefile"

htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"

rule all:
	input: rules.Index.output, htmlOutputPath + "/readme.html"
	output: touch("Output/all.done")

rule generateTissues:
	input: "data/generalTissues.txt", "data/detailedTissues.txt"
	output: touch("Output/dirStructure.done")
	shell: "while read line; do mkdir -p \'Scripts/Tissues/General/$line\'; done < data/generalTissues.txt && while read line; do mkdir -p \'Scripts/Tissues/Detailed/$line\'; done < data/detailedTissues.txt"
