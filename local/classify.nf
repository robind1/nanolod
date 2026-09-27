workflow CLASSIFY {

    take:
    ch_reads     
    ch_kraken2_db
    ch_diamond_db

    main:

    emit:
    versions = Channel.empty()
}
