@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 18864)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 16652)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 15844)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 16676)

del /F cleanup-ansys-hlfWin2-16676.bat
