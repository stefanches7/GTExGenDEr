while read line; do
        undersc=${line//[ -]/_}	# - is not allowed in snakemake rules names
	undersc=${undersc//(/8_}
	undersc=${undersc//)/_9}
	mkdir -p "Scripts/GeneralTissues/$undersc"
 done < data/generalTissues.txt && 
while read line; do 
        undersc=${line//[ -]/_}	
	undersc=${undersc//(/8_}
	undersc=${undersc//)/_9}
	mkdir -p "Scripts/DetailedTissues/$undersc"
 done < data/detailedTissues.txt
