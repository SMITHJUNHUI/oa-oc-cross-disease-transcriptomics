@echo off
setlocal
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_project.ps1" %*
if errorlevel 1 (
  echo.
  echo Project failed. See results\logs for details.
  pause
  exit /b 1
)
echo.
echo Project completed.
pause

