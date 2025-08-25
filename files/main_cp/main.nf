// Define variables
params.IMAGE_DIR = "/PATH/TO/IMAGES/"
params.OUTPUT_DIR = "/PATH/TO/main/results"
params.PIPELINE = "/PATH/TO/main/2024-07-16_main_pipeline_nocut.cpproj"
params.IMAGE_RANGES = "/PATH/TO/main/image_ranges_main.txt"

/*
 * Run cell profiler
 */
process run_cellprofiler {
    maxForks 20  // Adjust based on your memory test
    cache 'true'
    conda "/DKEC/users/mihat/.conda/envs/cp4"

    input:
        tuple val(first_image), val(last_image), val(output_subdir)

    output:
        //path("${output_subdir}")

    script:
    """
    mkdir -p ${output_subdir}
    OMP_NUM_THREADS=2 cellprofiler -c -r -p $params.PIPELINE \
        -i $params.IMAGE_DIR -o ${output_subdir} \
        -f $first_image -l $last_image
    """
}

workflow {
    image_ranges = Channel.fromPath(params.IMAGE_RANGES)
                      .splitCsv(sep: ' ', strip: true)
                      .map { row -> tuple(row[0].toInteger(), row[1].toInteger()) }
                      .map { first_image, last_image -> tuple(first_image, last_image, "${params.OUTPUT_DIR}/instance_${first_image}_${last_image}") }

    run_cellprofiler(image_ranges)
}
