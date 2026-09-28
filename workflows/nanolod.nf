
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { workflowVersionToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nanolod_pipeline'

workflow NANOLOD {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    Channel.topic('versions')
        .map { process, tool, version ->
            [ process.tokenize(':')[-1], tool, version.toString() ]
        }
        .unique()
        .map { process, tool, version ->
            new org.yaml.snakeyaml.Yaml().dumpAsMap([ (process): [ (tool): version ] ]).trim()
        }
        .mix(Channel.of(workflowVersionToYAML()))
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nanolod_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    versions       = ch_collated_versions        // channel: [ path(versions.yml) ]

}
