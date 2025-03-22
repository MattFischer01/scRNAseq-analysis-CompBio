while read -r id; do
    nohup fasterq-dump "$id" && gzip "$id"*.fastq &
done < srr_ids.txt
