# Code for Week 3 and 4: using git

`week3.ipynb` is the jupyter notebook for week 3, which capable of solving all questions.

`week3-question4.py` is customized for question 4 to enable it to be run on the command line.

To run it on all local fna data through command line, use `week3-question4.py` file, then run `find . -type f -name "*.fna" -exec python week3-question4.py {} ;` (For Windows please use `gci -r *.fna | % { python week3-question4.py $_.FullName }` instead).

