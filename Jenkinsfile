pipeline {
    agent any
    options {
        timeout(time: 30, unit: 'MINUTES')
    }
    stages{
        stage('Build') {
            steps {
                dir('.') {
                    sh 'mkdir build'
                }
                dir('build'){
                    deleteDir()
                    sh 'cmake ..'
                    sh 'make'
                }
            }
        }
    }


}