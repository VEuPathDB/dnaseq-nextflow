#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process runTest {

  input:
    path test

  output:
    path 'check.txt'

  script:
    """
    prove -l $test > check.txt
    """
}


workflow runTests {

  take:

    tests_qch
    
  main:
    
    runTest(tests_qch) | collectFile(storeDir: params.outputDir, name: 'tests.txt')

}