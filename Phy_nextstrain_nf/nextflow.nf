#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

/*
 * pipeline input parameters
 */

Channel
  .fromPath(params.sequences, checkIfExists: true)
  .set {sequences}
Channel
  .fromPath(params.metadata, checkIfExists: true)
  .set {metadata}

process align {
    publishDir "${params.outdir}/01.alignment", mode: 'copy', pattern: "*"
    input:
      path sequences

    output:
      path "alignment.fasta",emit:alignseq

    script:
    """
      augur align \
          --sequences ${sequences} \
          --reference-sequence ${params.reference} \
          --output alignment.fasta \
          --nthreads 2 \
          --remove-reference | tee logs

    """
}

process tree {
    publishDir "${params.outdir}/02.tree", mode: 'copy', pattern: "*"
    input:
      path alignseq

    output:
      path "tree_raw.nwk",emit:rawtree

    script:
      """
      FastTreeMP \
          -boot 1000 \
          -gtr \
          -nt ${alignseq} > tree_raw.nwk | tee logs

      """
}
process assign_genotype {
    publishDir "${params.outdir}/02.tree/", mode: 'copy', pattern: "*"
    input:
      path sequence_file
      path tree_file

    output:
      path "assign_geno.tsv",emit:genotypeTab

    script:
      """
      python3 ${baseDir}/scripts/assign_genotype.py \
          --sequence_file ${sequence_file} \
          --tree_file ${tree_file} \
          --dataframe_file ${params.metadata} \
          --output_file assign_geno.tsv
      """
}

// Rule refine
process refine {
    publishDir "${params.outdir}/02.tree/", mode: 'copy', pattern: "*"
    input:
      path rawtree
      path alignseq
      path genotypeTab

    output:
      path "tree.nwk",emit:refinetree
      path "branch_lengths.json",emit:node_length_data

    script:
// you need set your root-seq by using --root params;
//for example: --root hCoV-19/Wuhan/Hu-1/2019 hCoV-19/Wuhan/WH01/2019
      """
      augur refine \
          --tree ${rawtree} \
          --alignment ${alignseq} \
          --metadata ${genotypeTab} \
          --output-node-data branch_lengths.json \
          --output-tree tree.nwk \
          --date-inference joint \
          --timetree \
          --divergence-units ${params.units} \
          --coalescent COALESCENT\
          --branch-length-inference joint\
          --root ${params.root}  | tee logs
      """
}

process ancestral {
    publishDir "${params.outdir}/02.tree/", mode: 'copy', pattern: "*"
    input:
      path refinetree
      path alignseq
    output:
      path "nt_muts.json",emit:node_nt_data

    script:
      """
      augur ancestral \
          --tree $refinetree \
          --alignment $alignseq \
          --output-node-data nt_muts.json \
          --inference joint \
          --keep-ambiguous | tee logs

      """
}

process translate {
    publishDir "${params.outdir}/02.tree/", mode: 'copy', pattern: "*"
    input:
      path refinetree
      path node_nt_data

    output:
      path "aa_muts.json",emit:node_aa_data

    script:
      """
      augur translate \
          --tree ${refinetree} \
          --ancestral-sequences ${node_nt_data} \
          --reference-sequence ${params.referencegb} \
          --output-node-data aa_muts.json


      """
}

process traits {
    echo true
    publishDir "${params.outdir}/03.traits/", mode: 'copy', pattern: "*"

    input:
      path refinetree
      path genotypeTab
    output:
      path "traits.json",emit:node_traits_data

    script:
      """
      augur traits \
          --tree ${refinetree} \
          --metadata ${genotypeTab} \
          --output traits.json \
          --columns ${params.columns} \
          --confidence | tee logs

      """
}

process colors {
    echo true
    publishDir "${params.outdir}/03.traits/", mode: 'copy', pattern: "*"
    input:
      path genotypeTab
    output:
      path "colors.tsv",emit:colors

    script:
      """
      python3 ${baseDir}/scripts/assign-colors.py \
          --ordering ${params.ordering} \
          --color-schemes ${params.color_schemes} \
          --output colors.tsv \
          --metadata ${genotypeTab}
      """
}


process export {
    publishDir "${params.outdir}/04.auspice/", mode: 'copy', pattern: "*"
    input:
      path refinetree
      path node_data
      path colors
      path genotypeTab

    output:
      path "result_accessions.json", emit:auspice_json

    script:
      """
      augur export v2 \
          --tree ${refinetree} \
          --metadata ${genotypeTab} \
          --node-data ${node_data} \
          --auspice-config ${params.auspice_config} \
          --colors $colors \
          --lat-longs ${params.lat_longs} \
          --output result_accessions.json
      """
}

// Rule finalize
process finalize {
    publishDir "${params.outdir}/04.auspice/", mode: 'copy', pattern: "*"
    input:
      path auspice_json

    output:
      path "final.json",emit:seqtree

    script:
      """
      python3 ${baseDir}/scripts/fix-colorings.py \
          --input ${auspice_json} \
          --output final.json
      """
}


workflow {
        align=align(sequences)
        tree=tree(align.alignseq)
        assign_genotype=assign_genotype(align.alignseq,tree.rawtree)
        refine=refine(tree.rawtree,align.alignseq,assign_genotype.genotypeTab)
        ancestral=ancestral(refine.refinetree,align.alignseq)
        translate=translate(refine.refinetree,ancestral.node_nt_data)
        traits=traits(refine.refinetree, assign_genotype.genotypeTab)
        colors=colors(assign_genotype.genotypeTab)
        export = export(refine.refinetree,refine.node_length_data
          .mix(ancestral.node_nt_data)
          .mix(translate.node_aa_data)
          .mix(traits.node_traits_data).collect()
          ,colors.colors, assign_genotype.genotypeTab)
        final_file = finalize(export.auspice_json)
}
