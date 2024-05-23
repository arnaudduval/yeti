pipeline {
    agent any
    options {
        timeout(time: 30, unit: 'MINUTES')
    }
    stages{
        stage('Build') {
            steps {
                echo 'Building...'
                sh 'which python'
                sh 'which cmake'
            }
        }
    }


}