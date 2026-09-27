
workflow BUILD_ABLATED_DB {

    take:
    ch_target    
    ch_reference 

    main:

    emit:
    versions = Channel.empty()
}
