#!/bin/bash
cwd=$(pwd)

echo "Removing files in directory $cwd"
find . -name '*.tex' -or -name '*.ps' -or -name '*.pdf' -or -name '*.fig' -or -name '*.log' | xargs rm -f
FILES=$(find -L . -maxdepth 1 -type f -iname '*.gpl')
for f in $FILES
do
    echo "Processing file $f ..."
    gnuplot $f
done

FILES=$(find -L . -maxdepth 1 -type f -iname '*.fig')

for i in $FILES;
do
    fig2dev -L pdftex -p dummy $i ${i%???}pdf;
    fig2dev -L pdftex_t -F -p ${i%????} $i ${i%???}tex;
    #cp ${i%???}pdf $cwd/../fig/
    #cp ${i%???}tex $cwd/../fig/

done

FILES=$(find -L . -maxdepth 1 -type f -iname 'tex_*.tex')
for i in $FILES;
do
    tmp=${i%.*};
    pdflatex $i;
    echo "$tmp";
    pdfcrop ${tmp}.pdf;
    #mv ${tmp}-crop.pdf ${tmp}.pdf;
    #cp $tmp.pdf ../pics/${tmp#tex_}.pdf

done

mv *.pdf ./doc/
mv *.tex ./doc/
mv *.fig ./doc/
