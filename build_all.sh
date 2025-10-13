#!/bin/sh

for v in 3.9.23 3.10.18 3.11.13 3.12.11 3.13.8 3.14.0; do
    if [ -f build_${v}.log ]; then
      rm build_${v}.log
    fi
    py="$HOME/.pyenv/versions/$v/bin/python"

    echo creating venv for Python $v
    venv_dir=".venv_$v"
    $py -m venv $venv_dir
    . $venv_dir/bin/activate

    unset PYTHONPATH
    unset PYTHONHOME
    unset PYTHON_PATH

    python --version
    echo installing requirements
    python -m pip install --upgrade pip >> build_${v}.log
    python -m pip install -r requirements.txt >> build_${v}.log
    python setup.py build >> build_${v}.log 2>&1

    deactivate

    rm -rf $venv_dir

    echo "completed build for Python $v"
    echo
done
