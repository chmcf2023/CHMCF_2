docker run -it --rm \
    -v ~/projects/chmcf/scenic/stb:/data \
	-v ~/projects/share/SCENIC/:/share \
    aertslab/pyscenic:0.10.0 pyscenic grn \
        --num_workers 32 \
		--transpose \
        -o /data/expr_mat.adjacencies.tsv \
        /data/input/stb.counts.csv \
        /share/hs_hgnc_tfs.txt

docker run -it --rm \
    -v ~/projects/chmcf/scenic/stb:/data \
	-v ~/projects/share/SCENIC/:/share \
    aertslab/pyscenic:0.10.0 pyscenic ctx \
        /data/expr_mat.adjacencies.tsv \
        /share/hg19-tss-centered-10kb-7species.mc9nr.feather \
        --annotations_fname /share/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
        --expression_mtx_fname /data/input/stb.counts.csv \
        --mode "dask_multiprocessing" \
        --output /data/regulons.csv \
		--transpose \
        --num_workers 32

docker run -it --rm \
    -v ~/projects/chmcf/scenic/stb:/data \
    -v ~/projects/share/SCENIC/:/share \
    aertslab/pyscenic:0.10.0 pyscenic aucell \
        /data/input/stb.counts.csv \
        /data/regulons.csv \
        -o /data/auc_mtx.csv \
		--transpose \
        --num_workers 32
