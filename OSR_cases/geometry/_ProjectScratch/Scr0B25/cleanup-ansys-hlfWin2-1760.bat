@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 5700)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 16796)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 17204)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 1760)

del /F cleanup-ansys-hlfWin2-1760.bat
