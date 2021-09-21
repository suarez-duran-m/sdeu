#!/bin/bash


nl=$(wc -l $HOME/2021/sdeu/fullUubStationsListVert.txt | awk '{print $1}')

cp sdeuMonitoring_base.tex tmp.tex

for st in $(seq 1 ${nl});
do
  stId=$(head -${st} $HOME/2021/sdeu/fullUubStationsListVert.txt | tail -1)
    
  echo "\begin{frame} ">>tmp.tex
  echo "  \frametitle{Station ${stId}}">>tmp.tex
  echo "  \begin{center}">>tmp.tex
  echo "    \begin{columns}">>tmp.tex
  echo "      \begin{column}{0.33\textwidth}">>tmp.tex
  echo "        \includegraphics[width=1.15\textwidth]{../plots/qpkRmsCdasStSt${stId}Pmt1.pdf}">>tmp.tex
  echo "      \end{column}">>tmp.tex
  echo "      \begin{column}{0.33\textwidth}">>tmp.tex
  echo "        \includegraphics[width=1.15\textwidth]{../plots/qpkRmsCdasStSt${stId}Pmt2.pdf}">>tmp.tex
  echo "      \end{column}">>tmp.tex
  echo "      \begin{column}{0.33\textwidth}">>tmp.tex
  echo "        \includegraphics[width=1.15\textwidth]{../plots/qpkRmsCdasStSt${stId}Pmt3.pdf}">>tmp.tex
  echo "      \end{column}">>tmp.tex
  echo "    \end{columns}">>tmp.tex
  echo "  \end{center}">>tmp.tex
  echo "">>tmp.tex
  echo "  \begin{center}">>tmp.tex
  echo "    \begin{columns}">>tmp.tex
  echo "      \begin{column}{0.33\textwidth}">>tmp.tex
  echo "        \includegraphics[width=1.15\textwidth]{../plots/qpkRmsOffStSt${stId}Pmt1.pdf}">>tmp.tex
  echo "      \end{column}">>tmp.tex
  echo "      \begin{column}{0.33\textwidth}">>tmp.tex
  echo "        \includegraphics[width=1.15\textwidth]{../plots/qpkRmsOffStSt${stId}Pmt2.pdf}">>tmp.tex
  echo "      \end{column}">>tmp.tex
  echo "      \begin{column}{0.33\textwidth}">>tmp.tex
  echo "        \includegraphics[width=1.15\textwidth]{../plots/qpkRmsOffStSt${stId}Pmt3.pdf}">>tmp.tex
  echo "      \end{column}">>tmp.tex
  echo "    \end{columns}">>tmp.tex
  echo "  \end{center}">>tmp.tex
  echo "\end{frame}">>tmp.tex
  echo "" >> tmp.tex
  echo "" >> tmp.tex
done

echo "\end{document}" >> tmp.tex

mv tmp.tex sdeuMonitoring.tex
