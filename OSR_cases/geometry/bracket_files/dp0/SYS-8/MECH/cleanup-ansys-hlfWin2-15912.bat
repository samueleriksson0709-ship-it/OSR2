@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 16264)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 15308)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 1928)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 15912)

del /F cleanup-ansys-hlfWin2-15912.bat
