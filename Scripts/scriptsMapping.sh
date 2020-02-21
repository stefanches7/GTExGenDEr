# maps template script to single tissue directories
# these are parametrized using wBuild wildcards
template=Scripts/_Template/singleTissueAnalysis.R
linkFileName=singleTissueAnalysis.ln.R

for d in Scripts/GeneralTissues/*/; do # only directories in the wildcard
	cp $template "$d/$linkFileName"
done

for d in Scripts/DetailedTissues/*/; do # only directories in the wildcard
	cp $template "$d/$linkFileName"
done


