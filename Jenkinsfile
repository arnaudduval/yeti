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
                    sh '. .venv/bin/activate && pip install numpy\\<=1.26.4 scipy'
                    sh '. .venv/bin/activate && pip install sphinx sphinx-rtd-theme sphinxcontrib-bibtex tomli'
                    sh '. .venv/bin/activate && pip install cibuildwheel'
                    sh 'lsb_release -a'
                    sh 'docker version'
                }
            }
        }
        stage('configure legacy') {
            steps {
                dir('build') {
                    sh '. .venv/bin/activate && cmake ..'
                }
            }
        }
        stage('build legacy') {
            steps {
                dir('build') {
                    sh '. .venv/bin/activate && make -j4'
                }
            }
        }
        stage('test legacy') {
            steps {
                dir('build') {
                    sh '. .venv/bin/activate && ctest'
                }
            }
        }
        stage('build wheel') {
            steps {
                dir('.') {
                    sh '. build/.venv/bin/activate && cibuildwheel --platform linux --output-dir wheelhouse'
                }
            }
        }
        stage('doc') {
            steps {
                dir('build') {
                    sh '. .venv/bin/activate && make yeti-user-manual'
                }
            }
        }

    }


}