@echo off
chcp 65001 >nul
cd /d "%~dp0"
where julia >nul 2>nul
if errorlevel 1 (
  echo No encuentro julia en el PATH.
  pause
  exit /b 1
)
echo ETAPA C: DE con objetivo suavizado y sin buscador local (des),
echo en las 4 conjeturas refutadas donde la suavizacion cambia algo (29, 30, 166, 189).
echo 30 corridas de hasta 300 s. Completa la comparacion de la etapa A.
echo.
julia -t 8 scripts\experiment.jl --ids 29,30,166,189 --algs des --runs 30 --time 300
julia scripts\summary.jl
julia scripts\wilcoxon.jl
pause
