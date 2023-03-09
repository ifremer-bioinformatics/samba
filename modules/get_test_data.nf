process get_test_data {

    label 'internet_access'

    output:
        val('data_is_ready'), emit: data4test_ready
        path('*_sample_file.xls'), emit: xls

    script:
    def datatype = params.data_type
    """
    00_get_test_data.sh ${baseDir} data_is_ready ${datatype} ${datatype}_sample_file.xls &> 00_get_test_data.log 2>&1
    """

}
