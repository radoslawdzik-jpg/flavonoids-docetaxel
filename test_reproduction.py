from pathlib import Path
import numpy as np,pandas as pd
from src.extended import bh
from src.pipeline import calculate_bliss

def test_bh():
    assert np.allclose(bh([.01,.04,.03]),[.03,.04,.04])
def test_bio_hierarchy():
    d=pd.read_csv('data/processed/du145_biological_means.csv')
    assert set(d.bio_experiment)=={1,2,3}
    assert d.groupby(['compound','time_h','assay','flavonoid_conc_uM','docetaxel_nM']).size().eq(3).all()
def test_bliss():
    d=pd.read_csv('data/processed/du145_biological_means.csv')
    x=calculate_bliss(d)
    assert len(x)==144
