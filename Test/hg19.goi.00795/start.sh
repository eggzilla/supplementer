###############EXAMPLE PATHS, USE YOUR ONE ENVIRONMENT VARIABLES!
#EBOLA_WORK="/home/hoelzer/00HOLY_GENES_EBOLA/"
#EBOLA_IGV="/home/fall/Offlinework/Ebola_Jena/00HOLY_EBOLA_GENES/IGV_2.3.39/igv.jar"
#EBOLA_MOUNT="/mnt/mahlzeitlocal/fight_against_ebola"
#EBOLA_HACKER="Martin"

################# IF CHECKER = 'TRUE' load folder anyway and write line to CENTRAL_FILE
CHECKER=$1

if [ "$CHECKER" == "TRUE" ]; then
############################### COPY GENE ANYWAY AND WRITE IT TO CENTRAL FILE

    ruby $EBOLA_MOUNT"/scripts/write_gene_central.rb" "hg19.goi.00795" $EBOLA_HACKER $EBOLA_MOUNT

    FOLDER="hg19.goi.00795"
    OLD_FOLDER=$EBOLA_MOUNT"/holy_folder/goi/"$FOLDER

    EBOLA_WORK=$EBOLA_WORK"/"
    EBOLA_MOUNT=$EBOLA_MOUNT"/"

    echo "$(basename $OLD_FOLDER) copied and paths adjusted."
    mkdir -p $EBOLA_WORK/$FOLDER/

    if [ -e $OLD_FOLDER/*.readcounts.csv ]; then
	cp $OLD_FOLDER/*.readcounts.csv $EBOLA_WORK/$FOLDER
    fi

    for batch in $OLD_FOLDER/*.batch; do
        cp $batch $EBOLA_WORK/$FOLDER
    done
    for sh in $OLD_FOLDER/*.sh; do
        cp $sh $EBOLA_WORK/$FOLDER
    done

    cp $OLD_FOLDER/*.tex $EBOLA_WORK/$FOLDER

        for batch in $OLD_FOLDER/*.batch; do
            basename=$(basename $batch)
            awk -v x=$EBOLA_MOUNT -v y=$EBOLA_WORK/$FOLDER/snapshots/ '{gsub(/\/mnt\/mahlzeitlocal\/fight_against_ebola\/holy_folder\/goi\/.*\/snapshots\//,y,$0); gsub(/\/mnt\/mahlzeitlocal\/fight_against_ebola/,x,$0); print $0}'  $batch > $EBOLA_WORK/$FOLDER/$basename
        done

        for sh in $OLD_FOLDER/*.sh; do
            basename=$(basename $sh)
            awk -v x=$EBOLA_WORK -v y=$EBOLA_IGV '{gsub(/\/mnt\/mahlzeitlocal\/fight_against_ebola\/software\/IGV_2.3.39\/igv.jar/,y,$0); gsub(/\/mnt\/mahlzeitlocal\/fight_against_ebola\/holy_folder\/goi\//,x,$0);print $0}'  $sh > $EBOLA_WORK/$FOLDER/$basename
        done

else

    ################################ CHECK IF GENE IS ALREADY PRESENT
    ruby $EBOLA_MOUNT"/scripts/ruby_mine_project/central_gene_checker.rb" "hg19.goi.00795" $EBOLA_HACKER $EBOLA_MOUNT 2> $EBOLA_WORK"/tmp"

    tmp=`grep ATTENTION $EBOLA_WORK"/tmp"`
    if ! [[ -z "$tmp" ]]; then
        cat $EBOLA_WORK"/tmp"
    else
        cat $EBOLA_WORK"/tmp"

        FOLDER="hg19.goi.00795"
        OLD_FOLDER=$EBOLA_MOUNT"/holy_folder/goi/"$FOLDER

        EBOLA_WORK=$EBOLA_WORK"/"
        EBOLA_MOUNT=$EBOLA_MOUNT"/"

        echo "$(basename $OLD_FOLDER) copied and paths adjusted."
        mkdir -p $EBOLA_WORK/$FOLDER/

        if [ -e $OLD_FOLDER/*.readcounts.csv ]; then
            cp $OLD_FOLDER/*.readcounts.csv $EBOLA_WORK/$FOLDER
        fi

        for batch in $OLD_FOLDER/*.batch; do
            cp $batch $EBOLA_WORK/$FOLDER
        done
        for sh in $OLD_FOLDER/*.sh; do
             cp $sh $EBOLA_WORK/$FOLDER
        done

        cp $OLD_FOLDER/*.tex $EBOLA_WORK/$FOLDER

        for batch in $OLD_FOLDER/*.batch; do
            basename=$(basename $batch)
            awk -v x=$EBOLA_MOUNT -v y=$EBOLA_WORK/$FOLDER/snapshots/ '{gsub(/\/mnt\/mahlzeitlocal\/fight_against_ebola\/holy_folder\/goi\/.*\/snapshots\//,y,$0); gsub(/\/mnt\/mahlzeitlocal\/fight_against_ebola/,x,$0); print $0}'  $batch > $EBOLA_WORK/$FOLDER/$basename
        done

        for sh in $OLD_FOLDER/*.sh; do
            basename=$(basename $sh)
            awk -v x=$EBOLA_WORK -v y=$EBOLA_IGV '{gsub(/\/mnt\/mahlzeitlocal\/fight_against_ebola\/software\/IGV_2.3.39\/igv.jar/,y,$0); gsub(/\/mnt\/mahlzeitlocal\/fight_against_ebola\/holy_folder\/goi\//,x,$0);print $0}'  $sh > $EBOLA_WORK/$FOLDER/$basename
        done
    fi
fi
