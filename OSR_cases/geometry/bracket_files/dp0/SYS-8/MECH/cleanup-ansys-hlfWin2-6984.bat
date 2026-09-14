@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 14276)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 13072)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 17092)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 6984)

del /F cleanup-ansys-hlfWin2-6984.bat
