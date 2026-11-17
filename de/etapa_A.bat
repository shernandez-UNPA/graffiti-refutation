@echo off
chcp 65001 >nul
cd /d "%~dp0"
where julia >nul 2>nul
if errorlevel 1 (
  echo No encuentro julia en el PATH. Instalalo con juliaup o agrega julia al PATH.
  pause
  exit /b 1
)
echo Verificando el codigo...
julia test\runtests.jl
if errorlevel 1 (
  echo Fallaron las pruebas. No se lanza el experimento.
  pause
  exit /b 1
)
echo.
echo ETAPA A: 13 conjeturas ya refutadas, 4 algoritmos, 30 corridas, 300 s maximo por corrida.
echo Se puede cortar con Ctrl+C y volver a lanzar: retoma donde quedo.
echo No dejes que la laptop entre en suspension mientras corre.
echo.
julia -t 8 scripts\experiment.jl --only refuted --runs 30 --time 300
julia scripts\summary.jl
echo.
echo Listo. Resultados en results\runs.csv, results\summary.csv y results\summary.tex
pause
