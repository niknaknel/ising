subdir=$1

make
mkdir -p ${subdir}

# T=1.5, var init, time_seed=true
./ising 0 > "${subdir}/out01.csv"
./ising 1 > "${subdir}/out02.csv"
./ising 2 > "${subdir}/out03.csv"

Rscript ising_plots.R ${subdir}
