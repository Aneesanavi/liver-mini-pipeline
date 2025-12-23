nextflow.enable.dsl = 2

params.counts_path = "data/GSE126848/GSE126848_Gene_counts_raw.txt"
params.series_path = "data/GSE126848_series_matrix.txt"
params.outdir      = "results"
params.run_ml      = true
params.run_shap    = true

workflow {

  // PREPARE (emits a tuple via 'prepared')
  PREPARE_DATA(
    file(params.counts_path),
    file(params.series_path)
  )

  // channel: tuple(raw_counts, metadata)
  prepared_ch = PREPARE_DATA.out.prepared

  // QC
  QC_PCA(prepared_ch)

  // DESEQ2 (consumes prepared tuple)
  DESEQ2(prepared_ch)

  // Use process outputs directly
  MAKE_HEATMAP(
    DESEQ2.out.res,
    DESEQ2.out.norm,
    DESEQ2.out.meta
  )

  ANNOT_ENRICH(DESEQ2.out.res)

  if (params.run_ml) {
    ML(
      MAKE_HEATMAP.out.hm_matrix,
      MAKE_HEATMAP.out.meta_used,
      params.run_shap
    )
  }
}


process PREPARE_DATA {
  tag "prepare_data"
  publishDir ".", mode: 'copy'

  input:
    path counts_path
    path series_path

  output:
    tuple path("data/raw_counts.tsv"), path("data/metadata.tsv"), emit: prepared

  script:
  """
  mkdir -p data/GSE126848
  Rscript ${projectDir}/pipeline/scripts/02_prepare_data.R \
    ${counts_path} \
    ${series_path} \
    data/raw_counts.tsv \
    data/metadata.tsv \
    data/GSE126848/sample_map.tsv
  """
}

process QC_PCA {
  tag "qc_pca"
  publishDir params.outdir, mode: 'copy'

  input:
    tuple path(raw_counts), path(metadata)

  output:
    path("qc_summary.tsv")
    path("qc_library_size.png")
    path("qc_detected_genes.png")
    path("pca_logcpm.png")
    path("pca_logcpm_control_vs_disease.png")
    path("sample_correlation.png")

  script:
  """
  Rscript ${projectDir}/pipeline/scripts/03_qc_pca.R ${raw_counts} ${metadata} .
  """
}

process DESEQ2 {
  tag "deseq2"
  publishDir params.outdir, mode: 'copy'

  input:
    tuple path(raw_counts), path(metadata)

  output:
    path("deseq2_results.tsv"), emit: res
    path("deseq2_normalized_counts.tsv"), emit: norm
    path("metadata_used.tsv"), emit: meta
    path("deseq2_significant_padj_lt_0.05.tsv")
    path("deseq2_MAplot.png")
    path("deseq2_volcano.png")
    path("deseq2_summary.txt")

  script:
  """
  Rscript ${projectDir}/pipeline/scripts/04_deseq2.R ${raw_counts} ${metadata} .
  """
}

process MAKE_HEATMAP {
  tag "heatmap"
  publishDir params.outdir, mode: 'copy'

  input:
    path deseq_res
    path norm_counts
    path meta_used

  output:
    path("top50_heatmap_matrix_scaled.tsv"), emit: hm_matrix
    path("top50_genes_for_heatmap.tsv")
    path("top50_heatmap.png")
    path("meta_for_ml.tsv"), emit: meta_used

  script:
  """
  cp ${meta_used} meta_for_ml.tsv
  Rscript ${projectDir}/pipeline/scripts/06_heatmap_top50.R \
    ${deseq_res} \
    ${norm_counts} \
    meta_for_ml.tsv \
    .
  """
}

process ANNOT_ENRICH {
  tag "annot_enrich"
  publishDir params.outdir, mode: 'copy'

  input:
    path deseq_res

  output:
    path("deseq2_results_annotated.tsv")
    path("top10_genes_by_padj.tsv")
    path("top50_genes_by_padj.tsv")
    path("enrichGO_BP.tsv")
    path("enrichGO_BP_dotplot.png", optional: true)
    path("enrichKEGG.tsv", optional: true)
    path("enrichKEGG_dotplot.png", optional: true)

  script:
  """
  Rscript ${projectDir}/pipeline/scripts/05_annot_enrich.R ${deseq_res} .
  """
}

process ML {
  tag "ml"
  publishDir params.outdir, mode: 'copy'

  input:
    path hm_matrix
    path meta_for_ml
    val  run_shap

  output:
    path("ml_features_used.tsv")
    path("ml_holdout_logreg_summary.txt")
    path("ml_cv_logreg_summary.txt")
    path("logreg_coefficients.tsv")
    path("shap_importance.tsv", optional: true)
    path("shap_values_train.tsv", optional: true)
    path("shap_importance_bar.png", optional: true)

  script:
  """
  Rscript ${projectDir}/pipeline/scripts/07_ml_pca_shap.R \
    ${meta_for_ml} \
    ${hm_matrix} \
    . \
    ${run_shap}
  """
}
