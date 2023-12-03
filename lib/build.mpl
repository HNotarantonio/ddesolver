##############
# PATH refers to the path to the file were ddesolver.mla is located.
# It can be found by using the terminal, going to the location of interest,
#    and by typing pwd
##############
system("rm PATH_TO/ddesolver/lib/ddesolver.mla"):
read "PATH_TO/ddesolver/lib/ddesolver.mpl":
LibraryTools[Create]("PATH_TO/ddesolver/lib/ddesolver.mla"):
savelib('ddesolver',"PATH_TO/ddesolver/lib/ddesolver.mla");
quit;