#!groovy

@Library('pipelib@github-creds')
import org.veupathdb.lib.Builder

node('centos8') {
  def builder = new Builder(this)

  checkout scm
  builder.buildContainers([
    [ name: 'dnaseqAnalysis' ]
  ])
}
