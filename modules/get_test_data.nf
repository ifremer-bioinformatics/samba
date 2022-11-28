process get_test_data {

    label 'internet_access'

    output:
        val('data_is_ready'), emit: data4test_ready
        path('excel_sample_file.xls'), emit: xls

    script:
    def datatype = params.longreads ? "longreads" : "shortreads"
    """
    00_get_test_data.sh ${baseDir} data_is_ready ${datatype} excel_sample_file.xls &> 00_get_test_data.log 2>&1
    """

}
