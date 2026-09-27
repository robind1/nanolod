
workflow SIMULATE_AND_SPIKE {

    take:
    ch_target     
    ch_background 
    ch_abundance  

    main:

    emit:
    versions = Channel.empty()
}
