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
                }
            }
        }
        stage('Stage2') {
            steps {
                dir('.'){
                    sh 'mkdir plop'
                }
            }
        }
    }


}