################################
### If gene finished, copy your whole gene_id folder back to the
### central storage

FOLDER=`pwd`
BASENAME=$(basename $FOLDER)
var=$(echo $BASENAME | awk -F"." '{print $1,$2,$3}')
set -- $var
TYPE=$2

mkdir -p $EBOLA_MOUNT/holy_folder/"$TYPE"_out/$BASENAME
OUT=$EBOLA_MOUNT/holy_folder/"$TYPE"_out/$BASENAME/

cp -r $FOLDER/snapshots $OUT
cp *.tex $OUT

#echo "Your gene folder with id $BASENAME was succesfully copied to $EBOLA_MOUNT/holy_folder/"$TYPE"_out/"
chmod 775 $OUT
for i in `find $OUT/*`; do if [ -d "$i" ]; then chmod 775 $i; else chmod 664 $i; fi; done
