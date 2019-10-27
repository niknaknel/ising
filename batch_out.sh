subdir=$1

make
mkdir -p ${subdir}

# T=0, var init, time_seed=true
./ising 1.5 0 1 > "${subdir}/out_001_0.csv"
./ising 1.5 1 1 > "${subdir}/out_011_0.csv"
./ising 1.5 2 1 > "${subdir}/out_021_0.csv"

# T=1, var init, time_seed=true
./ising 1.8 0 1 > "${subdir}/out_101_0.csv"
./ising 1.8 1 1 > "${subdir}/out_111_0.csv"
./ising 1.8 2 1 > "${subdir}/out_121_0.csv"

# T=2, var init, time_seed=true
./ising 2 0 1 > "${subdir}/out_201_0.csv"
./ising 2 1 1 > "${subdir}/out_211_0.csv"
./ising 2 2 1 > "${subdir}/out_221_0.csv"

# T=3, var init, time_seed=true
./ising 2.3 0 1 > "${subdir}/out_301_0.csv"
./ising 2.3 1 1 > "${subdir}/out_311_0.csv"
./ising 2.3 2 1 > "${subdir}/out_321_0.csv"

# T=4, var init, time_seed=true
./ising 2.8 0 1 > "${subdir}/out_401_0.csv"
./ising 2.8 1 1 > "${subdir}/out_411_0.csv"
./ising 2.8 2 1 > "${subdir}/out_421_0.csv"

Rscript ising_plots.R ${subdir}
