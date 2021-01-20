import pandas_plink as pp
from dask.diagnostics import ProgressBar
for i in range(1, 24):
        G = pp.read_plink1_bin("PATH_TO_CHROMOSOME_DATA_" + str(i) + '.bed')
        G = G.astype('int8')
        G = G.to_dataset()
        with ProgressBar():
            G.compute()
        with ProgressBar():
            G.to_zarr('PATH_TO_ZARR_'+str(i))