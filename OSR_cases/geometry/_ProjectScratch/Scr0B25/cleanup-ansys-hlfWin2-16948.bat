@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 5740)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 6572)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 17496)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 16948)

del /F cleanup-ansys-hlfWin2-16948.bat
