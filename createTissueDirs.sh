while read line; do
        undersc=${line//[ ]/_}	
	mkdir -p "Scripts/GeneralTissues/$undersc"
 done < data/generalTissues.txt && 
while read line; do 
        undersc=${line//[ ]/_}	
	mkdir -p "Scripts/DetailedTissues/$undersc"
 done < data/detailedTissues.txt
