// Define variables
params.IMAGE_DIR="/PATH/TO/IMAGES/"
params.OUTPUT_DIR="/PATH/TO/preprocessing/qc"
params.PIPELINE = "/PATH/TO/preprocessing/2024-07-20_QC.cpproj"
params.IMAGE_RANGES = "/PATH/TO/preprocessing/image_ranges_qc.txt"

/*
 * Run cell profiler
 */
process run_cellprofiler {
    maxForks 20  // Adjust based on your memory test
    cache 'true'
	conda "/DKEC/users/mihat/.conda/envs/cp4"

    input:
        tuple val(first_image), val(last_image)

    output:
        //path("$params.OUTPUT_DIR")

    script:
    """
    OMP_NUM_THREADS=2 cellprofiler -c -r -p $params.PIPELINE \
        -i $params.IMAGE_DIR -o $params.OUTPUT_DIR \
        -f $first_image -l $last_image
    """
}

workflow {
    image_ranges = Channel.fromPath(params.IMAGE_RANGES)
                      .splitCsv(sep: ' ', strip: true)
                      .map { row -> tuple(row[0].toInteger(), row[1].toInteger()) }

    run_cellprofiler(image_ranges)
}
