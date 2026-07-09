@echo off
rem Launch the LineDetector pipeline GUI with the ArcGIS base Python.
rem (Swap python.exe -> pythonw.exe to hide this console window once you trust it.)
set PY="C:\Program Files\ArcGIS\Pro\bin\Python\envs\arcgispro-py3\python.exe"
%PY% "%~dp0pipeline_gui.py"
if errorlevel 1 pause
