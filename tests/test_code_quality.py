"""
Test global code quality with pylint
"""

import os
import subprocess
import sys

def get_python_files(directory):
    """
    Recursively get python files from a given directory
    """
    python_files = []
    for root, dirs, files in os.walk(directory):
        # Filter to avoid hidden directories
        dirs[:] = [d for d in dirs if not d.startswith('.')]
        for file in files:
            if file.endswith(".py"):
                python_files.append(os.path.join(root, file))
    return python_files

def run_pylint(file_path):
    """
    Run pylint on a file and return score
    """
    try:
        # Exécute pylint en capturant la sortie
        result = subprocess.run(
            ["pylint", file_path, "--score=y", "--output-format=text"],
            capture_output=True,
            text=True,
        )
        # Extraire le score de la sortie
        for line in result.stdout.splitlines():
            if line.strip().startswith("Your code has been rated at"):
                score = float(line.split(" ")[6].split("/")[0])
                return score
    except Exception as e:
        print(f"Erreur avec {file_path}: {e}")
        return None
    return None


def main(directory, threshold, percentage):
    """
    Get files, run pylint and compute percentage of files under a given
    threshold score
    """
    python_files = get_python_files(directory)
    if not python_files:
        print("No Python file found.")
        sys.exit(1)

    total_files = len(python_files)
    below_threshold = 0

    print(f"Analyze {total_files} Python files in '{directory}'...")
    for file_path in python_files:
        score = run_pylint(file_path)
        if score is not None:
            print(f"{file_path}: {score:.2f}")
            if score < threshold:
                below_threshold += 1
        else:
            print(f"Cannot get a score for {file_path}.")

    # Calculer le pourcentage de fichiers sous le seuil
    below_percentage = (below_threshold / total_files) * 100

    print(f"\nFiles under thershold ({threshold}): {below_threshold}/{total_files}")
    print(f"% of files under thershold: {below_percentage:.2f}%")

    # Vérifier si le pourcentage dépasse la limite
    if below_percentage > percentage:
        print(f"Fail : {below_percentage:.2f}% of files are under threshold.")
        sys.exit(1)
    else:
        print("Success")
        sys.exit(0)

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.realpath(__file__))
    directory_to_check = f'{script_dir}/..'
    score_threshold = 7.0
    max_percentage_below = 20.0

    main(directory_to_check, score_threshold, max_percentage_below)