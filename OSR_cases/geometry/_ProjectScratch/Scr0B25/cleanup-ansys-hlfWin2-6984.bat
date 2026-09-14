@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 14088)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 13904)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 3036)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 6984)

del /F cleanup-ansys-hlfWin2-6984.bat
