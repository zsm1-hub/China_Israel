FILE_NAME="Parcels_Iceland_wave"

# transfer to py and carry out
jupyter nbconvert --to python "$FILE_NAME.ipynb"
/home/simingzhang/anaconda3/bin/mpirun -np 1 python "$FILE_NAME.py"

