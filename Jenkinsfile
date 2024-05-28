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
                    sh '. .venv/bin/activate && pip install numpy scipy matplotlib nlopt'
                }
            }
        }
        stage('configure') {
            steps {
                dir('build') {
                    sh '. .venv/bin/activate && cmake ..'
                }
            }
        }
        stage('build') {
            steps {
                dir('build') {
                    sh '. .venv/bin/activate && make -j4'
                }
            }
        }
        stage('test') {
            steps {
                dir('build') {
                    sh '. .venv/bin/activate && ctest'
                }
            }
        }
        stage('doc') {
            steps {
                dir('build') {
                    sh '. .venv/bin/activate && pip install sphinx sphinx-rtd-theme sphinxcontrib-bibtex'
                    sh '. .venv/bin/activate && make yeti-user-manual'
                }
            }
        }

    }


}