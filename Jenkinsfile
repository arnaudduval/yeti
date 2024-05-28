pipeline {
    agent any
    options {
        timeout(time: 30, unit: 'MINUTES')
    }
    stages{
        stage('prepare') {
            steps {
                dir('.') {
                    sh 'mkdir build'
                }
                dir('build') {
                    deleteDir()
                    sh 'python3 -m venv .venv'
                    sh '. .venv/bin/activate'
                    sh 'pip install numpy scipy matplotlib nlopt'
                }
            }
        }
        stage('configure') {
            steps {
                dir('build') {
                    sh 'cmake ..'
                }
            }
        }
        stage('build') {
            steps {
                dir('build') {
                    sh 'make -j4'
                }
            }
        }
        stage('test') {
            steps {
                dir('build') {
                    sh 'ctest'
                }
            }
        }
        stage('doc') {
            steps {
                dir('build') {
                    sh 'pip install sphinx sphinx-rtd-theme sphinxcontrib-bibtex'
                }
            }
        }

    }


}