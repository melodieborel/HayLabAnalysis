while [ 1 ];
do
#    count=`curl -s "https://c3n-cn.fr/2025/01/31/classements-dadmissibilite-au-concours-cnrs-2025/" | grep -c "admissibles"`
    count=`curl -s "https://mts-ncomms.nature.com/cgi-bin/main.plex?form_type=status_details&j_id=18&ms_id=600463&ms_rev_no=0&ms_id_key=ftd83pQaHjSF2skLkRY18gYg"  | grep -c "under consideration"`
    echo "$count"
    if [ "$count" != 0 ]
    then
       echo "Manuscript updated!"
       afplay HavaNagila.mp3
       exit 0
    fi
    sleep 10
done