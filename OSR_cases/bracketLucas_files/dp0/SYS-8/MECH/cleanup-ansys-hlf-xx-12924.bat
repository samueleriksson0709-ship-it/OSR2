@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="hlf-xx" (taskkill /f /pid 6160)
if /i "%LOCALHOST%"=="hlf-xx" (taskkill /f /pid 3140)
if /i "%LOCALHOST%"=="hlf-xx" (taskkill /f /pid 9672)
if /i "%LOCALHOST%"=="hlf-xx" (taskkill /f /pid 12924)

del /F cleanup-ansys-hlf-xx-12924.bat
