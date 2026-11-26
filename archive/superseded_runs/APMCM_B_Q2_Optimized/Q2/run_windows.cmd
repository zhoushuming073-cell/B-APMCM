@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  py -3.12 -m venv .venv
  if errorlevel 1 goto fail
)
".venv\Scripts\python.exe" -m pip install -r requirements.txt
if errorlevel 1 goto fail
set OPENBLAS_NUM_THREADS=1
set OMP_NUM_THREADS=1
".venv\Scripts\python.exe" code\run.py
if errorlevel 1 goto fail
echo Completed. Open the tables and figures folders.
pause
exit /b 0
:fail
echo Run failed. Check that Python 3.12 is installed and pip can reach PyPI.
pause
exit /b 1
