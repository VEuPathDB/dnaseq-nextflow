#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process runTest {
  container = 'veupathdb/dnaseqanalysis'

  input:
    path test

  output:
    path 'check.txt'

  script:
    """
    prove -l $test > check.txt
    """

  stub:
    """
    touch check.txt
    """
}


workflow tests {

  take:

    tests_qch
    
  main:
    
    runTest(tests_qch) | collectFile(storeDir: params.outputDir, name: 'tests.txt')

}
