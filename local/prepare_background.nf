
workflow PREPARE_BACKGROUND {

    take:
    ch_reads 
    ch_host  

    main:

    emit:
    versions = Channel.empty()
}
