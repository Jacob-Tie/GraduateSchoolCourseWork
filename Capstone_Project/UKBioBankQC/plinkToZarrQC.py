import pandas_plink as pp
from dask.diagnostics import ProgressBar
for i in range(1, 23):
        G = pp.read_plink1_bin("PATH_TO_BED_TRAIN" + str(i) + '.bed')
        G = G.astype('int8')
        G = G.to_dataset()
        with ProgressBar():
            G.compute()
        with ProgressBar():
            G.to_zarr('PATH_TO_TRAINING_DATA'+str(i))